# -*- coding: utf-8 -*-
"""
_zero_batch.py —— 纯 fast_search 狂产可行解，粗筛零残留（对话里秒级出解那套）

用法：
  python docs/tcode/_zero_batch.py --ratio 1.5 --budget 36 --seconds 120 --want 3
"""
import argparse
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import solve_tcode as S
import gen_tcode_figures as G
import fast_search_v22 as F

SEP, GAP = 12, 2


def rebuild():
    S.CROSSTALK_MAX_GAP = GAP
    S.set_sep(SEP)
    S.CLASSES = S.build_classes(gap=GAP)
    S.AVOID_PAIRS = S.build_avoid_true(gap=GAP)
    S.VAR2CLS = [[] for _ in S.VARS]
    for ci, (_, prs) in enumerate(S.CLASSES):
        for ia, ib in prs:
            S.VAR2CLS[ia].append(ci)
            S.VAR2CLS[ib].append(ci)
    S.VAR2CLS = [sorted(set(c)) for c in S.VAR2CLS]
    S.VAR2AVOID = [[] for _ in S.VARS]
    for pi, (ia, ib) in enumerate(S.AVOID_PAIRS):
        S.VAR2AVOID[ia].append(pi)
        S.VAR2AVOID[ib].append(pi)
    F.GAP = GAP
    F.LASER_IDS = list(S.LASER_IDS)
    F.KICKS_OF = {l: tuple(S.KICKS_OF[l]) for l in S.LASER_IDS}
    F.VIDX = S.VIDX


def residual(u, ratio, step):
    fn = S.make_code_fn(u)
    fr = G.build_firings(fn)
    ga = gb = kill = tb = 0
    for D in np.arange(5.0, 601.0, step):
        r = G.simulate(D, fr, max_gap=GAP, ratio=ratio)
        ga += r["ga"]; gb += r["gb"]; kill += r["kill"]; tb += r["tb"]
        if ga or kill:
            return ga, gb, kill, tb, False
    return ga, gb, kill, tb, True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ratio", type=float, default=1.5)
    ap.add_argument("--budget", type=int, required=True)
    ap.add_argument("--seconds", type=float, default=120.0)
    ap.add_argument("--want", type=int, default=1)
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()

    rebuild()
    B = a.budget
    groups = F.laser_groups()
    quads = F.valid_quads(B, SEP)
    print(f"B={B}ns ratio={a.ratio} 合法4元组={len(quads)} 组={len(groups)}",
          flush=True)
    if not quads:
        print("低于硬下界"); return

    rng = np.random.default_rng(a.seed + B * 31)
    t0 = time.time()
    n_try = n_feas = 0
    best_ga = None
    zeros = []

    while time.time() - t0 < a.seconds and len(zeros) < a.want:
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
        ga, gb, kill, tb, clean = residual(u, a.ratio, 10.0)
        if not clean:
            if best_ga is None or ga < best_ga:
                best_ga = ga
            if n_feas % 50 == 0:
                print(f"  [{time.time()-t0:5.1f}s] 可行{n_feas}/{n_try} "
                      f"最好残留={best_ga}", flush=True)
            continue
        ga, gb, kill, tb, clean = residual(u, a.ratio, 2.0)
        if clean:
            zeros.append(u.copy())
            print(f"  [{time.time()-t0:5.1f}s] ★零残留#{len(zeros)} "
                  f"可行{n_feas}/{n_try} 鬼{gb} 真{tb}", flush=True)
            name = f"tcode_table_zero_r{a.ratio}_{int(u.max())}ns.py"
            S.dump_table(u, os.path.join(HERE, name), SEP, int(u.max()),
                         note=f"零残留；ratio={a.ratio} gap={GAP} SEP={SEP} "
                              f"（_zero_batch / fast_search）")

    print(f"\n结束：尝试{n_try} 可行{n_feas} 零残留{len(zeros)} "
          f"最好残留={best_ga} {time.time()-t0:.0f}s", flush=True)
    if zeros:
        u = zeros[0]
        for l in S.LASER_IDS:
            print(f"  L{l:<2d} " + " ".join(
                f"K{k}={int(u[S.VIDX[(l,k)]]):>2d}" for k in S.KICKS_OF[l]))


if __name__ == "__main__":
    main()
