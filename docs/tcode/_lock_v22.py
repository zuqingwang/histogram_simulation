# -*- coding: utf-8 -*-
"""锁定 v22 码表：优先 36ns（已探到可行），找不到则退到 40/48。"""
import sys
import numpy as np
sys.path.insert(0, ".")
import solve_tcode as S

SEP = 12
S.set_sep(SEP)


def collect(B, n_seed=30, steps=40000):
    sols = []
    for seed in range(n_seed):
        u = S.solve(B, np.random.default_rng(seed), max_steps=steps)
        if u is not None:
            sols.append(u)
            print(f"  B={B} seed={seed} OK max={int(max(u))}", flush=True)
            if len(sols) >= 5:
                break
    return sols


chosen_B, cands = None, []
for B in (36, 40, 48):
    print(f"搜索预算 {B} ns ...", flush=True)
    cands = collect(B)
    if cands:
        chosen_B = B
        break
if not cands:
    raise SystemExit("36/40/48 均无解")

print(f"采用预算 {chosen_B}，评估 {len(cands)} 组解 ...", flush=True)
scored = []
for u in cands:
    t = S.evaluate(S.make_code_fn(u), 1.6, 20.0)
    scored.append((t["ga"], t["kill"], int(max(u)), u))
scored.sort()
best = scored[0][3]
fn = S.make_code_fn(best)
ok, fails = S.check_avoid_true(fn, SEP)
print(f"最优 max={max(best)} 避真={ok} 违规={len(fails)}", flush=True)
for ratio in (1.6, 2.5):
    t = S.evaluate(fn, ratio, step=2.0)
    print(f"  ratio={ratio}: {t['gb']}->{t['ga']} "
          f"({t['ga']/max(t['gb'],1):.3%}) kill {t['kill']}/{t['tb']}", flush=True)
S.dump_table(best, "tcode_table_v22.py", SEP, chosen_B,
             note=f"v22: 散开+避真 SEP={SEP}ns；探测最小可行约36ns，本表预算{chosen_B}ns")
for l in S.LASER_IDS:
    print(f"  L{l:<3d} " + " ".join(
        f"K{k}={int(best[S.VIDX[(l,k)]])}" for k in S.KICKS_OF[l]))
