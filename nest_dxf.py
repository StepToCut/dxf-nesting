"""nest_dxf.py - true-shape nesting of DXF parts on rectangular sheets.

One part per DXF file. Each part's closed outline becomes a polygon, the polygon a raster, and an FFT convolution
finds every free position on the sheet at once. Among the free positions the program takes the one that extends
the used sheet length least. Output: one DXF with the sheets and the parts (original geometry, rotated and moved),
plus a CSV report.

    python nest_dxf.py PARTS_DIR --sheet 3000x1500 --border 10 --gap 10 --out nest.dxf

Quantity: a file named  bracket__q4.dxf  means 4 pieces (default 1).
Tested with Python 3.14, ezdxf 1.4.4, shapely 2.1.2, scipy 1.17.1, numpy 2.4.6.
"""
import argparse, csv, math, sys, time
from pathlib import Path

import numpy as np
import ezdxf
from ezdxf import path as ezpath
from ezdxf.math import Matrix44
from scipy.signal import fftconvolve
from scipy.spatial import cKDTree
from shapely import affinity, contains_xy
from shapely.geometry import LineString, Polygon
from shapely.ops import polygonize, unary_union

FLATTEN = 0.1          # mm, chord error when arcs and splines become polylines


def read_part(file):
    """Closed outline of the part (holes removed) and the original modelspace entities."""
    doc = ezdxf.readfile(file)
    msp = doc.modelspace()
    lines, ents = [], []
    for e in msp:
        if e.dxftype() in ("TEXT", "MTEXT", "DIMENSION", "HATCH", "POINT", "INSERT"):
            continue  # text and marks are copied but do not take part in the collision
        try:
            pth = ezpath.make_path(e)
            pts = [(v.x, v.y) for v in pth.flattening(FLATTEN)]
        except Exception:
            continue
        closed = pth.is_closed or getattr(e, "closed", False) or getattr(e, "is_closed", False)
        if len(pts) >= 2 and (closed or math.dist(pts[0], pts[-1]) < 0.5) and pts[0] != pts[-1]:
            pts.append(pts[0])   # close the ring: a 0.1 mm gap is enough to lose the outline
        if len(pts) >= 2:
            lines.append(LineString(pts))
    for e in msp:
        ents.append(e)
    blocks = sum(1 for e in msp if e.dxftype() == "INSERT")
    if blocks:   # geometry inside a block is not read: it would be lost without a word
        raise ValueError(f"{blocks} block(s) in the file: explode them before export")
    # bridge gaps under 0.5 mm between the open ends of different entities (lines drawn one by one)
    ends = [pt for ln in lines if not ln.is_closed for pt in (ln.coords[0], ln.coords[-1])]
    if ends:
        for a, b in cKDTree(ends).query_pairs(0.5):
            if ends[a] != ends[b]:
                lines.append(LineString([ends[a], ends[b]]))
    polys = list(polygonize(unary_union(lines)))
    if not polys:
        raise ValueError("no closed outline: open contour or gap over 0.5 mm")
    # the outline is the face with the largest OUTER area: a frame's hole can have more area than the frame itself
    outer = max(polys, key=lambda p: Polygon(p.exterior).area)
    outline = Polygon(outer.exterior)
    if not outline.buffer(0.5).contains(unary_union(lines)):
        raise ValueError("geometry outside the outline: open contour, gap over 0.5 mm, or two parts in one file")
    holes = [p for p in polys if p is not outer and outline.contains(p.representative_point())]
    shape = outline
    for h in holes:
        shape = shape.difference(Polygon(h.exterior))
    return outline, shape, doc, ents


def raster(poly, res, gap):
    """Boolean mask of the part grown by half the gap, at res mm per cell, with its origin in mm.

    A cell is marked when its centre falls inside the outline grown by gap/2 + half a cell diagonal: every cell the
    grown outline touches is then marked. So two masks that do not overlap guarantee a real distance >= gap.
    (Our first version dilated the raster by gap/2 with binary_dilation's default, cross-shaped structure: on diagonal
    edges parts came out 7.07 mm apart for a 10 mm gap.)
    """
    grown = poly.buffer(gap / 2 + res * math.sqrt(2) / 2, quad_segs=8)
    minx, miny, maxx, maxy = grown.bounds
    ox, oy = minx - res, miny - res
    w = int(math.ceil((maxx - ox) / res)) + 1
    h = int(math.ceil((maxy - oy) / res)) + 1
    xs = ox + (np.arange(w) + 0.5) * res
    ys = oy + (np.arange(h) + 0.5) * res
    gx, gy = np.meshgrid(xs, ys)
    mask = contains_xy(grown, gx.ravel(), gy.ravel()).reshape(h, w)
    return mask, ox, oy


def orientations(outline):
    """0, 90, 180, 270 degrees plus the angle of the minimum rotated rectangle."""
    angles = [0.0, 90.0, 180.0, 270.0]
    rect = outline.minimum_rotated_rectangle
    c = list(rect.exterior.coords)
    a = math.degrees(math.atan2(c[1][1] - c[0][1], c[1][0] - c[0][0])) % 90.0
    if 0.5 < a < 89.5:
        angles += [-a, 90.0 - a]
    return angles


def nest(parts, sheet_w, sheet_h, border, gap, res):
    # The marked cells of a part always contain the part grown by gap/2 (see raster), so the grid may reach gap/2
    # beyond the cuttable area and every outline still ends up inside it. Without this, parts kept ~23 mm from a 10 mm border.
    g0 = border - gap / 2
    H = int((sheet_h - 2 * g0) // res)   # rows = sheet height
    W = int((sheet_w - 2 * g0) // res)   # cols = sheet length
    sheets = []                              # each: dict(occ, frontier, placed)
    for p in sorted(parts, key=lambda q: -q["outline"].area):
        cands = []
        for ang in orientations(p["outline"]):
            rot = affinity.rotate(p["outline"], ang, origin=(0, 0))
            m, ox, oy = raster(rot, res, gap)
            if m.shape[0] > H or m.shape[1] > W:
                continue
            cands.append((ang, m, ox, oy))
        if not cands:
            raise ValueError(f"{p['name']} does not fit on the sheet")
        placed = False
        for s in sheets + [None]:
            if s is None:
                s = {"occ": np.zeros((H, W), bool), "frontier": 0, "placed": []}
                sheets.append(s)
            best = None
            for ang, m, ox, oy in cands:
                mh, mw = m.shape
                lim = max(mw, min(W, s["frontier"] + mw + 2))
                conv = fftconvolve(s["occ"][:, :lim].astype(np.float64), m[::-1, ::-1].astype(np.float64), mode="valid")
                free = np.argwhere(conv < 0.5)    # float64: in float32 the rounding gave a 107% yield
                if free.size == 0:
                    continue
                order = np.lexsort((free[:, 0], free[:, 1] + mw))  # least extension of the sheet, then lowest
                j, i = free[order[0]]
                cost = (int(i) + mw, int(j))
                if best is None or cost < best[0]:
                    best = (cost, ang, m, ox, oy, int(i), int(j))
            if best is None:
                continue
            _, ang, m, ox, oy, i, j = best
            s["occ"][j:j + m.shape[0], i:i + m.shape[1]] |= m
            s["frontier"] = max(s["frontier"], i + m.shape[1])
            dx = g0 + i * res - ox
            dy = g0 + j * res - oy
            s["placed"].append({**p, "angle": ang, "dx": dx, "dy": dy})
            placed = True
            break
        if not placed:
            raise ValueError(f"{p['name']} could not be placed")
    return sheets


def used_length(s, border):
    """Where the sheet can be cut off: the rightmost point of the real parts plus the border."""
    return max(affinity.translate(affinity.rotate(q["outline"], q["angle"], origin=(0, 0)), q["dx"], q["dy"]).bounds[2]
               for q in s["placed"]) + border


def check(sheets, sheet_w, sheet_h, border, gap):
    """Independent of the raster: real polygons, shapely distances. Returns a list of problems."""
    problems = []
    for k, s in enumerate(sheets):
        usable = Polygon([(border, border), (sheet_w - border, border), (sheet_w - border, sheet_h - border), (border, sheet_h - border)])
        shapes = []
        for q in s["placed"]:
            g = affinity.translate(affinity.rotate(q["outline"], q["angle"], origin=(0, 0)), q["dx"], q["dy"])
            if not usable.buffer(1e-6).contains(g):
                problems.append(f"sheet {k + 1}: {q['name']} outside the usable area")
            shapes.append((q["name"], g))
        for a in range(len(shapes)):
            for b in range(a + 1, len(shapes)):
                d = shapes[a][1].distance(shapes[b][1])
                if d < gap - 1e-6:
                    problems.append(f"sheet {k + 1}: {shapes[a][0]} - {shapes[b][0]} {d:.2f} mm apart (gap {gap})")
    return problems


def write_dxf(sheets, out, sheet_w, sheet_h, spacing):
    doc = ezdxf.new("R2010")
    msp = doc.modelspace()
    for k, s in enumerate(sheets):
        oy = -k * (sheet_h + spacing)
        msp.add_lwpolyline([(0, oy), (sheet_w, oy), (sheet_w, oy + sheet_h), (0, oy + sheet_h)], close=True, dxfattribs={"layer": "SHEET"})
        msp.add_text(f"SHEET {k + 1}", dxfattribs={"layer": "SHEET", "height": 40}).set_placement((0, oy + sheet_h + 20))
        for q in s["placed"]:
            m = Matrix44.chain(Matrix44.z_rotate(math.radians(q["angle"])), Matrix44.translate(q["dx"], q["dy"] + oy, 0))
            for e in q["ents"]:
                try:
                    c = e.copy()
                    c.transform(m)        # transform() also rotates text: do not rotate it again
                    msp.add_foreign_entity(c)
                except Exception as ex:
                    print(f"  ! {q['name']}: {e.dxftype()} not copied ({ex})")
    doc.saveas(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("parts_dir")
    ap.add_argument("--sheet", default="3000x1500", help="length x height in mm")
    ap.add_argument("--border", type=float, default=10.0, help="non-cuttable border, mm")
    ap.add_argument("--gap", type=float, default=10.0, help="minimum distance between parts, mm")
    ap.add_argument("--raster", type=float, default=5.0, help="raster cell, mm (smaller = tighter and slower)")
    ap.add_argument("--out", default="nest.dxf")
    ap.add_argument("--bbox", action="store_true", help="nest bounding rectangles instead of true shapes (to compare)")
    a = ap.parse_args()
    sheet_w, sheet_h = (float(v) for v in a.sheet.lower().split("x"))

    parts = []
    for f in sorted(Path(a.parts_dir).glob("*.dxf")):
        name, _, q = f.stem.partition("__q")
        if q and not q.isdigit():
            sys.exit(f"{f.name}: quantity after __q must be a number")
        try:
            outline, shape, doc, ents = read_part(f)
        except Exception as ex:
            sys.exit(f"{f.name}: {ex}")
        for n in range(int(q) if q else 1):
            parts.append({"name": f"{name}#{n + 1}", "outline": outline.envelope if a.bbox else outline, "area": shape.area, "doc": doc, "ents": ents})
    if not parts:
        sys.exit("no DXF files found")

    t0 = time.time()
    try:
        sheets = nest(parts, sheet_w, sheet_h, a.border, a.gap, a.raster)
    except ValueError as ex:
        sys.exit(f"stopped: {ex}")
    secs = time.time() - t0
    problems = check(sheets, sheet_w, sheet_h, a.border, a.gap)
    write_dxf(sheets, a.out, sheet_w, sheet_h, spacing=200)

    rep = Path(a.out).with_suffix(".csv")
    with open(rep, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["sheet", "parts", "used_length_mm", "parts_area_m2", "yield_on_used_length_%"])
        tot_len = tot_area = 0.0
        for k, s in enumerate(sheets):
            used = used_length(s, a.border)
            area = sum(q["area"] for q in s["placed"])
            tot_len += used
            tot_area += area
            w.writerow([k + 1, len(s["placed"]), round(used), round(area / 1e6, 3), round(100 * area / (used * sheet_h), 1)])
    print(f"{len(parts)} parts on {len(sheets)} sheet(s), used length {tot_len / 1000:.2f} m, "
          f"yield on used length {100 * tot_area / (tot_len * sheet_h):.1f}%, {secs:.1f} s")
    print(f"check: {'OK, 0 overlaps, 0 parts outside' if not problems else str(len(problems)) + ' PROBLEMS'}")
    for p in problems[:20]:
        print("  " + p)
    print(f"written {a.out} and {rep}")
    sys.exit(1 if problems else 0)


if __name__ == "__main__":
    main()
