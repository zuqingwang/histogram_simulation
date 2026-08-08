# -*- coding: utf-8 -*-
"""
子问题精确判定 —— 回答「32ns 到底有没有解」
==============================================
取一条最紧的链：L5 - L7 - L9 - L11
  （编号间隔均为 2 ⟹ 三对都构成串扰；且四者共享 kick K0/K1/K2/K6）

只保留两类【必须满足】的约束：
  避真：同一 kick 上，相邻对 |c_a[k] - c_b[k]| >= SEP
  散开：每一对的 4 个码差两两相差 >= SEP

这是完整问题的【松弛】（约束更少）。因此：
  子问题无解  ==> 完整问题必定无解   （严格证明）
  子问题有解  ==> 完整问题未必有解

算法：
  记 d1=c5-c7, d2=c7-c9, d3=c9-c11。三个坐标各自要求 4 值两两相隔 >= SEP，
  彼此只通过「存在 c7 使四个码都落在 [0,B]」这一条耦合。于是：
    ① 枚举 d1 的升序 4 元组（kick 之间对称，升序即可去重）
    ② 枚举 d2 到 4 个 kick 的赋值
    ③ 对每个 kick 算出 d3 的可行区间，再解一个 4 变量小 CSP
  规模从 C(n,4) 降到几万个小 CSP。
"""
SEP = 12


def valid_quads(B, sep):
    """所有满足 |d|>=sep 且两两相隔 >=sep 的升序 4 元组。"""
    vals = [d for d in range(-B, B + 1) if abs(d) >= sep]
    out = []

    def rec(start, cur):
        if len(cur) == 4:
            out.append(tuple(cur))
            return
        for i in range(start, len(vals)):
            v = vals[i]
            if cur and v - cur[-1] < sep:
                continue
            rec(i + 1, cur + [v])

    rec(0, [])
    return out


def perms(seq):
    from itertools import permutations
    return list(permutations(seq))


def d3_window(B, d1, d2):
    """给定 d1,d2，返回 d3 的可行闭区间；不可行返回 None。
       约束：存在 t=c7，使 t, t+d1, t-d2, t-d2-d3 全在 [0,B]。"""
    lo = max(0, -d1, d2)
    hi = min(B, B - d1, B + d2)
    if lo > hi:
        return None
    return (lo - d2 - B, hi - d2)


def solve_d3(B, sep, windows):
    """每个 kick 给一个 d3 允许区间，选 4 个值两两相隔 >=sep。"""
    cands = []
    for (w_lo, w_hi) in windows:
        s = [d for d in range(max(-B, w_lo), min(B, w_hi) + 1) if abs(d) >= sep]
        if not s:
            return False
        cands.append(s)

    chosen = []

    def bt(i):
        if i == len(cands):
            return True
        for v in cands[i]:
            if all(abs(v - p) >= sep for p in chosen):
                chosen.append(v)
                if bt(i + 1):
                    return True
                chosen.pop()
        return False

    return bt(0)


def feasible(B, sep=SEP, report=False):
    quads = valid_quads(B, sep)
    if len(quads) == 0:
        return False, 0, "单个坐标就放不下 4 个值"

    n_inner = 0
    for q1 in quads:                       # d1 升序（kick 对称性已去重）
        for q2 in quads:
            for a2 in perms(q2):           # d2 到 4 个 kick 的赋值
                windows = []
                ok = True
                for i in range(4):
                    w = d3_window(B, q1[i], a2[i])
                    if w is None:
                        ok = False
                        break
                    windows.append(w)
                if not ok:
                    continue
                n_inner += 1
                if solve_d3(B, sep, windows):
                    return True, n_inner, "找到可行解"
    return False, n_inner, "穷举完毕，确实无解"


if __name__ == "__main__":
    print(f"链 L5-L7-L9-L11（共享 4 个 kick），SEP={SEP} ns")
    print("这是完整问题的松弛 —— 判为无解即【证明】完整问题无解\n")
    print(f"  {'预算':>7} {'合法4元组':>10} {'结论':>8}   说明")
    print("  " + "-" * 60)
    for B in (20, 22, 23, 24, 26, 28, 30, 32, 36):
        quads = valid_quads(B, SEP)
        ok, n_inner, why = feasible(B)
        tag = "有解" if ok else "【无解】"
        print(f"  {B:>5d}ns {len(quads):>10d} {tag:>8}   {why}", flush=True)
