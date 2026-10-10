"""Invented test parts (none from a real job): plates with holes, angles, flanges, gussets, rounded covers."""
import ezdxf, math
from pathlib import Path
D = Path(__file__).parent / "test_parts"; D.mkdir(exist_ok=True)
for f in D.glob("*.dxf"): f.unlink()
def save(name, build, q):
    doc = ezdxf.new("R2010"); msp = doc.modelspace(); build(msp); doc.saveas(D / f"{name}__q{q}.dxf")
def plate(w, h, holes=()):
    def b(msp):
        msp.add_lwpolyline([(0,0),(w,0),(w,h),(0,h)], close=True)
        for x, y, r in holes: msp.add_circle((x, y), r)
    return b
def ell(a, b_, t):
    def b(msp): msp.add_lwpolyline([(0,0),(a,0),(a,t),(t,t),(t,b_),(0,b_)], close=True)
    return b
def gusset(a, b_, c):
    def b(msp): msp.add_lwpolyline([(0,0),(a,0),(a,c),(c,b_),(0,b_)], close=True)
    return b
def flange(r, rh, n, rb):
    def b(msp):
        msp.add_circle((0,0), r); msp.add_circle((0,0), rh)
        for k in range(n):
            a = 2*math.pi*k/n; msp.add_circle((0.72*r*math.cos(a), 0.72*r*math.sin(a)), rb)
    return b
def rounded(w, h, r):
    def b(msp):  # rectangle with rounded corners: bulge on the corner segments
        bl = math.tan(math.radians(90)/4)
        msp.add_lwpolyline([(r,0,0,0,0),(w-r,0,0,0,bl),(w,r,0,0,0),(w,h-r,0,0,bl),(w-r,h,0,0,0),(r,h,0,0,bl),(0,h-r,0,0,0),(0,r,0,0,bl)], format="xyseb", close=True)
        msp.add_text("PL-07", dxfattribs={"height": 20}).set_placement((w/2-30, h/2))
    return b
save("base_plate", plate(600, 400, [(50,50,11),(550,50,11),(50,350,11),(550,350,11)]), 4)
save("long_strip", plate(1400, 120, [(60,60,9),(1340,60,9)]), 3)
save("angle_L", ell(500, 350, 90), 6)
save("gusset", gusset(260, 260, 40), 12)
save("flange", flange(180, 80, 8, 9), 5)
save("cover", rounded(700, 300, 40), 3)
print(sorted(p.name for p in D.glob("*.dxf")))
