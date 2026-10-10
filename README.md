# DXF nesting

True-shape nesting of DXF parts on rectangular sheets, in Python, plus a separate script that checks the result.

- `nest_dxf.py` reads a folder of DXF files (one part per file), places the parts on sheets and writes `nest.dxf` and `nest.csv`.
- `verify_nest.py` re-reads `nest.dxf` without trusting the nesting script and measures the real distance between parts, parts outside the cuttable area, and the part count.

How it works, setup step by step, measured results and limits: [How to nest cut parts on steel sheets with Python, from DXF profiles](https://steptocut.com/blog/nesting-sheet-metal-python-dxf/).

## Quick start

```
pip install ezdxf shapely scipy numpy
python nest_dxf.py parts --sheet 3000x1500 --border 10 --gap 10 --out nest.dxf
python verify_nest.py nest.dxf --gap 10 --border 10 --parts 33
```

One part per DXF file, in mm, with a closed outer contour. `bracket__q4.dxf` means 4 pieces. Explode blocks before export. `--parts` is the number of parts the nesting printed. Look at the nest in your CAD before cutting.

`make_test_parts.py` writes the 33 test parts used on the page; `make_edge_cases.py` writes the awkward ones (a frame whose hole is bigger than its material, a plate drawn with separate lines and 0.1 mm gaps).

## Tested

Windows 11 (Python 3.14, ezdxf 1.4.4, shapely 2.1, scipy 1.17, numpy 2.4) and Ubuntu (Python 3.14, ezdxf 1.4.4, shapely 2.2, scipy 1.18, numpy 2.5). On the 33 test parts, sheet 3000 x 1500, gap 10: 3.34 m of sheet with true shapes against 3.97 m with bounding rectangles, smallest distance between parts 14.1 mm, 0 overlaps. The output opens in AutoCAD 2025 with 0 errors in AUDIT. Not tested on macOS.

## Limits

Greedy order (largest part first), no parts inside holes, no common-line cutting. Dimensions are not copied to the nest. The gap is the only clearance: set it for your kerf and process.

Written with Claude (AI), tested by [StepToCut](https://steptocut.com). No warranty: check every nest before cutting.

## License

MIT
