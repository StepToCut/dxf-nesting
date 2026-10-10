"""Two awkward test parts: a frame whose hole is bigger than its material, and a plate drawn with separate LINEs that stop 0.1 mm short of the corners, with holes."""
import ezdxf, shutil
from pathlib import Path
D = Path(__file__).parent / "edge"
shutil.rmtree(D, ignore_errors=True); D.mkdir()
def new(): d = ezdxf.new("R2010"); return d, d.modelspace()
d, m = new(); m.add_lwpolyline([(0,0),(1000,0),(1000,1000),(0,1000)], close=True); m.add_lwpolyline([(50,50),(950,50),(950,950),(50,950)], close=True); d.saveas(D/"frame__q3.dxf")
d, m = new()
pts=[(0,0),(600,0),(600,400),(0,400)]
for k in range(4):
    a=pts[k]; b=pts[(k+1)%4]
    # each LINE stops 0.1 mm short of the next corner
    dx=(b[0]-a[0]); dy=(b[1]-a[1]); L=(dx*dx+dy*dy)**0.5
    m.add_line(a,(b[0]-dx/L*0.1,b[1]-dy/L*0.1))
for x,y in [(50,50),(550,50),(50,350),(550,350)]: m.add_circle((x,y),11)
d.saveas(D/"gappy_plate__q3.dxf")
