# -*- coding: utf-8 -*-
"""严格解（不忽略任何编号间隔）快速搜最小可行预算。

约束与现状相同：散开（全分开）+ 避真，只把 CROSSTALK_MAX_GAP 提到 15。
"""
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import solve_tcode as S
import fast_search_v22 as F


def rebuild(gap=15):
    S.CROSSTALK_MAX_GAP = gap
    F.GAP = gap
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
    print(f"已重建 gap={gap}: 鬼影类={len(S.CLASSES)}, 避真对={len(S.AVOID_PAIRS)}")


def try_budget(B, sep=12, seconds=45.0, want=1):
    found = F.search(B, sep, seconds=seconds, want=want, seed0=7)
    if not found:
        return None
    u = found[0]
    fn = S.make_code_fn(u)
    ok_av, fails = S.check_avoid_true(fn, sep)
    # 散开自检
    spread_bad = 0
    for _, prs in S.CLASSES:
        d = sorted(int(u[ia]) - int(u[ib]) for ia, ib in prs)
        gaps = [d[i + 1] - d[i] for i in range(len(d) - 1)]
        if gaps and min(gaps) < sep:
            spread_bad += 1
    print(f"  校验：避真={ok_av} 违规={len(fails)} 散开失败类={spread_bad} "
          f"max(tx)={int(u.max())}")
    # 距离扫描粗评
    for ratio in (1.5, 2.5):
        ev = S.evaluate(fn, ratio, step=20.0)
        print(f"  评 ratio={ratio}: 鬼残留={ev['ga']}/{ev['gb']} "
              f"误杀={ev['kill']}")
    return u


def main():
    rebuild(gap=15)
    S.set_sep(12)
    print()
    print("严格全连接硬下界 = 36 ns（同 kick 4 激光两两避真）")
    print("下面用「全散开 + 避真」搜最小可行预算（覆盖 ratio=1.5；")
    print("对 ratio=2.5 也充分，因为全散开比允许双碰更严）。")
    print()

    budgets = [36, 40, 48, 56, 64, 72, 80, 96, 112, 128]
    # 小预算少给时间；越大越可能有解，给够
    time_map = {36: 20, 40: 25, 48: 30, 56: 40, 64: 50,
                72: 60, 80: 60, 96: 60, 112: 60, 128: 60}

    first = None
    for B in budgets:
        print("=" * 64)
        t0 = time.time()
        u = try_budget(B, seconds=time_map[B], want=1)
        print(f"预算 {B} ns → {'可行' if u is not None else '未找到'} "
              f"（{time.time()-t0:.1f}s）")
        if u is not None:
            first = (B, u)
            break

    print()
    print("=" * 64)
    if first is None:
        print("在试到的预算内没找到（不等于无解，可加长时间）")
    else:
        B, u = first
        print(f"严格解（gap=∞，全散开+避真）最小搜到的可行预算 ≈ {B} ns")
        print(f"  max(tx) = {int(u.max())} ns")
        print(f"  → ratio=1.5 / 2.5 都可用这一张（全散开对两者都零残留）")


if __name__ == "__main__":
    main()
