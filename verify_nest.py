"""verify_nest.py - check a nest DXF without trusting the program that made it.

Reads the output file only: sheet rectangles on layer SHEET, every other curve rebuilt into outlines from scratch.
Reports the smallest distance between two parts and from a part to the sheet edge, and any overlap.

    python verify_nest.py nest.dxf --gap 10 --border 10 --parts 33
"""
import argparse, sys
import ezdxf
from ezdxf import path as ezpath
from shapely.geometry import LineString, Polygon, box
from shapely.ops import polygonize, unary_union
from scipy.spatial import cKDTree


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dxf")
    ap.add_argument("--gap", type=float, default=10.0)
    ap.add_argument("--border", type=float, default=10.0)
    ap.add_argument("--parts", type=int, required=True, help="number of parts the nesting printed: a stacked part shows only in the count")
    a = ap.parse_args()
    msp = ezdxf.readfile(a.dxf).modelspace()
    sheets, lines = [], []
    for e in msp:
        if e.dxftype() in ("TEXT", "MTEXT", "DIMENSION", "HATCH", "POINT", "INSERT"):
            continue
        try:
            pts = [(v.x, v.y) for v in ezpath.make_path(e).flattening(0.1)]
        except Exception:
            continue
        if len(pts) < 2:
            continue
        if e.dxf.layer == "SHEET":
            sheets.append(Polygon(pts).envelope)
        else:
            if pts[0] != pts[-1] and (getattr(e, "closed", False) or getattr(e, "is_closed", False)):
                pts.append(pts[0])
            lines.append(LineString(pts))
    ends = [pt for ln in lines if not ln.is_closed for pt in (ln.coords[0], ln.coords[-1])]
    if ends:                              # same 0.5 mm bridging as the nesting script
        for i, j in cKDTree(ends).query_pairs(0.5):
            if ends[i] != ends[j]:
                lines.append(LineString([ends[i], ends[j]]))
    faces = list(polygonize(unary_union(lines)))
    # outlines = faces that are not inside another face (holes are inside their part)
    outlines = [Polygon(f.exterior) for f in faces]
    outlines = [o for o in outlines if not any(p is not o and p.contains(o.representative_point()) and p.area > o.area for p in outlines)]
    min_gap, min_edge, bad, outside = float("inf"), float("inf"), 0, 0
    orphans = [o for o in outlines if not any(s.contains(o.representative_point()) for s in sheets)]
    outside += len(orphans)
    for s in sheets:
        usable = box(*s.bounds).buffer(-a.border, join_style=2)
        mine = [o for o in outlines if s.contains(o.representative_point())]
        for o in mine:
            if usable.contains(o):
                min_edge = min(min_edge, usable.exterior.distance(o))
            else:
                outside += 1
        for i in range(len(mine)):
            for j in range(i + 1, len(mine)):
                d = mine[i].distance(mine[j])
                min_gap = min(min_gap, d)
                if mine[i].intersects(mine[j]):
                    bad += 1
    print(f"{len(sheets)} sheet(s), {len(outlines)} part outline(s) rebuilt from the file")
    print(f"smallest distance between parts: {min_gap:.2f} mm (required {a.gap})")
    print(f"smallest distance from a part to the cuttable edge: {min_edge:.2f} mm")
    print(f"parts outside the cuttable area: {outside}, overlapping pairs: {bad}")
    ok = bad == 0 and outside == 0 and min_gap >= a.gap - 1e-6
    # a part stacked on or inside another merges into one outline: only the count shows it
    print(f"part count: {len(outlines)} found, {a.parts} expected")
    ok = ok and len(outlines) == a.parts
    print("RESULT: " + ("OK" if ok else "PROBLEMS: see the lines above"))
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
