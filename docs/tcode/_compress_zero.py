# -*- coding: utf-8 -*-
"""
_compress_zero.py —— 在「鬼影零残留」前提下压缩 tcode 预算

判据顺序（不可颠倒）：
  1) 约束可行：散开 + 避真（cost=0）
  2) 验收零残留：1~600m 扫描，XM 后鬼影残留 = 0 且真峰误杀 = 0
  3) 再压预算

用法：
    python docs/tcode/_compress_zero.py 1.5
    python docs/tcode/_compress_zero.py 2.5
"""
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import solve_tcode as S
import gen_tcode_figures as G

SEP = 12
GAP = 2


def rebuild(gap=GAP, sep=SEP):
    S.CROSSTALK_MAX_GAP = gap
    S.set_sep(sep)
    S.CLASSES = S.build_classes(gap=gap)
    S.AVOID_PAIRS = S.build_avoid_true(gap=gap)
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


def hunt(B, ratio, restarts, seed0, screen=10.0, verify=2.0):
    """在预算 B 内找零残留解；返回 (u, 说明) 或 (None, 统计)。"""
    n_feas = 0
    for i in range(restarts):
        rng = np.random.default_rng(seed0 + i * 977)
        u = S.solve(B, rng, max_steps=80000, plateau=2500)
        if u is None:
            continue
        n_feas += 1
        ga, gb, kill, tb, clean = residual(u, ratio, screen)
        if not clean:
            continue
        ga, gb, kill, tb, clean = residual(u, ratio, verify)
        if clean:
            return u, f"可行{n_feas}/{i+1}次，细扫零残留（鬼{gb} 真{tb}）"
    return None, f"可行{n_feas}/{restarts}次，无零残留解"


def main():
    ratio = float(sys.argv[1]) if len(sys.argv) > 1 else 1.5
    rebuild()
    print(f"XM_RATIO={ratio}  SEP={SEP}ns  gap={GAP}")
    print(f"散开类={len(S.CLASSES)} 避真对={len(S.AVOID_PAIRS)}")
    print("目标：先零残留，再压预算\n")

    best = None
    # 从大到小压，直到某档找不到零残留解
    for B in (80, 72, 64, 56, 52, 48, 44, 40, 36, 32, 28, 24):
        t0 = time.time()
        u, msg = hunt(B, ratio, restarts=24, seed0=1000 + B * 31)
        tag = "★零残留" if u is not None else "  失败  "
        print(f"  B={B:>3d}ns {tag}  {msg}  [{time.time()-t0:.1f}s]", flush=True)
        if u is not None:
            best = (B, u.copy())
        else:
            break

    if best is None:
        print("\n未找到零残留解")
        return
    B, u = best
    name = f"tcode_table_zero_r{ratio}_{int(u.max())}ns.py"
    S.dump_table(u, os.path.join(HERE, name), SEP, int(u.max()),
                 note=f"零残留优先再压预算；ratio={ratio} gap={GAP} SEP={SEP}")
    print(f"\n最短零残留预算 = {B}ns（实际 max={int(u.max())}ns）")
    for l in S.LASER_IDS:
        print(f"  L{l:<2d} " + " ".join(
            f"K{k}={int(u[S.VIDX[(l,k)]]):>2d}" for k in S.KICKS_OF[l]))


if __name__ == "__main__":
    main()
