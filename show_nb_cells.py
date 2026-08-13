# -*- coding: utf-8 -*-
"""打印 notebook 的 cell 地图：序号 / 类型 / 行数 / 首行 / 产出的图。

用法：python show_nb_cells.py PoD_esti_v30.ipynb
"""
import json
import re
import sys

sys.stdout.reconfigure(encoding="utf-8")
path = sys.argv[1] if len(sys.argv) > 1 else "PoD_esti_v30.ipynb"
nb = json.load(open(path, encoding="utf-8"))
for i, c in enumerate(nb["cells"]):
    src = "".join(c.get("source", []))
    first = next((l.strip() for l in src.splitlines() if l.strip()), "")
    kind = "MD  " if c["cell_type"] == "markdown" else "CODE"
    figs = re.findall(r'savefig\("([^"]+)"', src)
    tag = ("   figs=" + ",".join(figs)) if figs else ""
    print(f"{i:3d} {kind} {len(src.splitlines()):4d}L  {first[:72]}{tag}")
print(f"\n合计 {len(nb['cells'])} cells")
