# -*- coding: utf-8 -*-
import os, sys
sys.path.insert(0, os.path.dirname(__file__))
import numpy as np
import solve_tcode as S
import gen_tcode_figures as G

path = os.path.join(os.path.dirname(__file__), "tcode_table_zero_r1.5_32ns.py")
ns = {}
exec(compile(open(path, encoding="utf-8").read(), path, "exec"), ns)
tbl = ns["TCODE_TABLE"]
S.CROSSTALK_MAX_GAP = 2
S.set_sep(12)
S.CLASSES = S.build_classes(gap=2)
S.AVOID_PAIRS = S.build_avoid_true(gap=2)
u = np.array([tbl[v] for v in S.VARS], dtype=np.int32)
print("max", int(u.max()), "cost", S.total_cost(u))
fn = S.make_code_fn(u)
fr = G.build_firings(fn)
for ratio in (1.5, 2.5):
    ga = gb = kill = tb = 0
    for D in np.arange(5.0, 601.0, 2.0):
        r = G.simulate(D, fr, max_gap=2, ratio=ratio)
        ga += r["ga"]; gb += r["gb"]; kill += r["kill"]; tb += r["tb"]
    tag = "零残留" if ga == 0 and kill == 0 else "有残留"
    print(f"ratio={ratio}: 鬼{gb}->{ga} 误杀{kill}/{tb}  {tag}")
