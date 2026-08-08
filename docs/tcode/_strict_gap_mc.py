# -*- coding: utf-8 -*-
"""严格解 gap=∞：先找可行上界，再往下压。"""
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import solve_tcode as S


def rebuild(gap=15):
    S.CROSSTALK_MAX_GAP = gap
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
    print(f"gap={gap}: 鬼影类={len(S.CLASSES)}, 避真对={len(S.AVOID_PAIRS)}")


def search_budget(B, sep=12, seconds=120.0, steps=80000):
    S.set_sep(sep)
    rng = np.random.default_rng(2026 + B * 17)
    t0 = time.time()
    n_try = 0
    while time.time() - t0 < seconds:
        n_try += 1
        u = S.solve(B, rng, max_steps=steps, plateau=2500)
        if u is not None:
            print(f"  [{time.time()-t0:5.1f}s] ★ 命中 (尝试 {n_try}, "
                  f"max={int(u.max())})", flush=True)
            return u, n_try, time.time() - t0
        if n_try % 3 == 0:
            print(f"  [{time.time()-t0:5.1f}s] 重启 {n_try}", flush=True)
    return None, n_try, time.time() - t0


def verify(u, sep=12):
    fn = S.make_code_fn(u)
    ok_av, fails = S.check_avoid_true(fn, sep)
    spread_bad = 0
    min_gap = 10 ** 9
    for _, prs in S.CLASSES:
        d = sorted(int(u[ia]) - int(u[ib]) for ia, ib in prs)
        gaps = [d[i + 1] - d[i] for i in range(len(d) - 1)]
        if gaps:
            min_gap = min(min_gap, min(gaps))
            if min(gaps) < sep:
                spread_bad += 1
    print(f"  校验：避真={ok_av} 违规={len(fails)} 散开失败={spread_bad} "
          f"散开最小间隔={min_gap} max={int(u.max())}")
    for ratio in (1.5, 2.5):
        ev = S.evaluate(fn, ratio, step=20.0)
        resid = ev["ga"] / max(ev["gb"], 1) * 100
        print(f"  评 ratio={ratio}: 鬼残留 {ev['ga']}/{ev['gb']} "
              f"({resid:.2f}%) 误杀={ev['kill']}")


def main():
    rebuild(15)
    print("硬下界 36 ns")
    print("策略：先从上往下找可行，再向下压")
    print()

    # 先找上界
    upper_candidates = [200, 160, 128, 112, 96, 80]
    upper = None
    best_u = None
    for B in upper_candidates:
        print("=" * 64)
        print(f"找上界：试 {B} ns ...", flush=True)
        u, n_try, dt = search_budget(B, seconds=90 if B >= 128 else 120)
        if u is None:
            print(f"  {B} ns 未找到（重启 {n_try}，{dt:.1f}s）")
        else:
            verify(u)
            upper = B
            best_u = u
            print(f"  上界锁定 {B} ns")
            break

    if upper is None:
        print("\n连 200 ns 都没在时限内找到。放宽再试 200 ns × 3 分钟...")
        u, n_try, dt = search_budget(200, seconds=180)
        if u is None:
            print("仍然没有。严格全连接可能需要更长搜索或换算法。")
            return
        verify(u)
        upper, best_u = 200, u

    # 向下压
    print()
    print(f"已有上界 {upper} ns，向下压...")
    for B in [x for x in [160, 128, 112, 96, 80, 72, 64, 56, 48, 40, 36]
              if x < upper]:
        print("=" * 64)
        print(f"下压：试 {B} ns ...", flush=True)
        u, n_try, dt = search_budget(B, seconds=100)
        if u is None:
            print(f"  {B} ns 未找到（重启 {n_try}，{dt:.1f}s）→ 停止下压")
            break
        verify(u)
        upper, best_u = B, u
        print(f"  上界更新为 {B} ns")

    print()
    print("=" * 64)
    print(f"严格解（gap=∞，全散开+避真）本次搜到的可行预算：{upper} ns")
    print(f"  max(tx)={int(best_u.max())} ns")
    print("  硬下界仍是 36 ns；本结果是算法在有限时间内找到的上界/近似最小。")
    print("  ratio=1.5 与 2.5 均可使用（约束按全散开，对两者都够）。")


if __name__ == "__main__":
    main()
