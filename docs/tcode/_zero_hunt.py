# -*- coding: utf-8 -*-
"""
_zero_hunt.py —— 用差分抽样求解器批量产可行解，筛出「零残留」的最短预算。

判据顺序：约束可行 → 1~600m 细扫零残留（鬼影 0 且误杀 0）→ 再压预算。

用法：python docs/tcode/_zero_hunt.py 1.5 40 --seconds 240
"""
import argparse
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import solve_tcode as S
import fast_search_v22 as F
from _compress_zero import rebuild, residual, GAP, SEP


def hunt(B, ratio, seconds, seed0=0):
    rebuild(sep=SEP)
    F.S.set_sep(SEP)
    groups = F.laser_groups()
    quads = F.valid_quads(B, SEP)
    if not quads:
        return None, "低于硬下界"
    rng = np.random.default_rng(seed0)
    t0 = time.time()
    n_try = n_feas = 0
    best_ga = None
    while time.time() - t0 < seconds:
        n_try += 1
        codes, ok = {}, True
        for ks, ls in groups:
            sol = None
            for _ in range(80):
                sol = F.try_group(B, SEP, ks, ls, quads, rng)
                if sol is not None:
                    break
            if sol is None:
                ok = False
                break
            codes.update(sol)
        if not ok:
            continue
        u = np.zeros(len(S.VARS), dtype=np.int32)
        for (l, k), vi in S.VIDX.items():
            u[vi] = codes[(l, k)]
        if S.total_cost(u) != 0:
            continue
        n_feas += 1
        ga, gb, kill, tb, clean = residual(u, ratio, 10.0)
        if clean:
            ga, gb, kill, tb, clean = residual(u, ratio, 2.0)
            if clean:
                return u, (f"可行{n_feas}/{n_try}，细扫零残留"
                           f"（鬼{gb} 真{tb}）{time.time()-t0:.0f}s")
        if best_ga is None or ga < best_ga:
            best_ga = ga
    return None, (f"可行{n_feas}/{n_try}，最好残留={best_ga}，"
                  f"{time.time()-t0:.0f}s")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ratio", type=float)
    ap.add_argument("budgets", type=int, nargs="+")
    ap.add_argument("--seconds", type=float, default=180.0)
    a = ap.parse_args()

    print(f"零残留优先 · ratio={a.ratio} SEP={SEP} gap={GAP}")
    for B in a.budgets:
        u, msg = hunt(B, a.ratio, a.seconds, seed0=7000 + B * 13)
        if u is None:
            print(f"  B={B:>3d}ns 失败  {msg}", flush=True)
            continue
        print(f"  B={B:>3d}ns ★零残留  {msg}", flush=True)
        name = f"tcode_table_zero_r{a.ratio}_{int(u.max())}ns.py"
        S.dump_table(u, os.path.join(HERE, name), SEP, int(u.max()),
                     note=f"零残留优先再压预算；ratio={a.ratio} gap={GAP} SEP={SEP}")
        for l in S.LASER_IDS:
            print("    L%-2d " % l + " ".join(
                f"K{k}={int(u[S.VIDX[(l,k)]]):>2d}" for k in S.KICKS_OF[l]))


if __name__ == "__main__":
    main()
