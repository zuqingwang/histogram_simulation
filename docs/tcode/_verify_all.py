# -*- coding: utf-8 -*-
"""核验 docs/tcode 下所有码表：约束 cost + 两种 ratio 的残留（细扫 2m）。"""
import os
import sys
from pathlib import Path

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import solve_tcode as S
from _compress_zero import rebuild, residual, GAP, SEP

rebuild()
rows = []
for p in sorted(Path(HERE).glob("tcode_table*.py")):
    ns = {}
    ns_src = p.read_text(encoding="utf-8")
    exec(compile(ns_src, str(p), "exec"), ns)
    tbl = ns["TCODE_TABLE"]
    if any(v not in tbl for v in S.VARS):
        print(f"{p.name:<38} 变量不匹配，跳过")
        continue
    u = np.array([tbl[v] for v in S.VARS], dtype=np.int32)
    cost = S.total_cost(u)
    out = {"name": p.name, "max": int(u.max()), "cost": cost}
    for r in (1.5, 2.5):
        ga, gb, kill, tb, clean = residual(u, r, 2.0)
        out[r] = (ga, gb, kill, clean)
    rows.append(out)

print(f"\nSEP={SEP} gap={GAP}，细扫 5~600m 步长 2m")
print(f"{'码表':<38} {'max':>4} {'cost':>5} "
      f"{'r1.5 残留':>12} {'r2.5 残留':>12}")
print("-" * 78)
for o in sorted(rows, key=lambda x: x["max"]):
    def fmt(r):
        ga, gb, kill, clean = o[r]
        return ("零残留" if clean and kill == 0 else f"{ga}/{gb}")
    print(f"{o['name']:<38} {o['max']:>4} {o['cost']:>5} "
          f"{fmt(1.5):>12} {fmt(2.5):>12}")
