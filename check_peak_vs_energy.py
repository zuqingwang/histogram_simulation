# -*- coding: utf-8 -*-
"""无头执行 `peak_vs_energy_v01.ipynb` 的全部 code cell，并把每张图存成 PNG。

用来在不开 Jupyter 的情况下确认 cell 能跑通、图能画出来、字体没缺字。
    python -u check_peak_vs_energy.py
"""
from __future__ import annotations

import json
import os
import sys
import traceback

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

NB = "peak_vs_energy_v01.ipynb"
PREFIX = "pve_"

_saved: list[str] = []
_cell_no = [0]


def _show(*a, **k):
    for num in plt.get_fignums():
        f = plt.figure(num)
        path = f"{PREFIX}m{_cell_no[0]:02d}_{len(_saved)}.png"
        f.savefig(path, dpi=105, bbox_inches="tight")
        _saved.append(path)
    plt.close("all")


plt.show = _show

BUILDER = "build_peak_vs_energy.py"
if os.path.exists(BUILDER) and os.path.getmtime(BUILDER) > os.path.getmtime(NB):
    # 构建脚本一旦语法出错，.ipynb 不会被重写，校验就会悄悄跑在旧 notebook 上，
    # 图看起来「没改动」，实际是根本没重建。这里直接拦下来。
    raise SystemExit(
        f"[中止] {BUILDER} 比 {NB} 新，说明上一次构建没成功。\n"
        f"        先单独跑 `python {BUILDER}` 看报错，构建成功后再校验。")

nb = json.load(open(NB, encoding="utf-8"))
cells = [(i, "".join(c["source"])) for i, c in enumerate(nb["cells"])
         if c["cell_type"] == "code"]

g: dict = {"__name__": "__main__"}
fail = 0
for i, src in cells:
    _cell_no[0] = i
    print("\n" + "=" * 96)
    print(f"cell {i}")
    print("=" * 96)
    try:
        exec(compile(src, f"<cell {i}>", "exec"), g)
    except Exception:                                            # noqa: BLE001
        fail += 1
        traceback.print_exc()
        print(f"[cell {i} 失败]")

# 缺字形自查：matplotlib 找不到字形时会发 UserWarning，这里再显式核一遍
from matplotlib import font_manager as _fm

_font = _fm.findfont(_fm.FontProperties(family=plt.rcParams["font.sans-serif"]))
print("\n" + "=" * 96)
print(f"实际使用字体：{_font}")
print(f"生成图片 {len(_saved)} 张：")
for p in _saved:
    print(f"  {p}")
print(f"\n{'[全部 cell 通过]' if fail == 0 else f'[{fail} 个 cell 失败]'}")
sys.exit(1 if fail else 0)
