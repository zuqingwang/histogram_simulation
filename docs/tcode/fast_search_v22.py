# -*- coding: utf-8 -*-
"""
fast_search_v22.py —— 码差抽样 + 差分求解（针对 24~36ns）
==========================================================
对每个发光组（4 激光 × 4 kick）：
  1. 给每条避真对随机抽一组合法码差 4 元组并打乱到 4 个 kick
  2. 每个 kick 上用差分约束求绝对码（带 [0,B] 界）
  3. 检查自身混叠 / 跨 kick 散开
四组合龙后验全局 cost，为 0 即成功。

用法：
    python docs/tcode/fast_search_v22.py --budget 32
    python docs/tcode/fast_search_v22.py --budget 28 --seconds 60 --want 5
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

LASER_IDS = G.LASER_IDS
KICKS_OF = {l: tuple(ks) for l, ks in G.KICKS_OF.items()}
VIDX = S.VIDX
GAP = S.CROSSTALK_MAX_GAP


def laser_groups():
    buckets = {}
    for l, ks in KICKS_OF.items():
        buckets.setdefault(ks, []).append(l)
    return [(ks, sorted(ls)) for ks, ls in buckets.items()]


def valid_quads(B, sep):
    vals = [d for d in range(-B, B + 1) if abs(d) >= sep]
    out = []

    def rec(s, cur):
        if len(cur) == 4:
            out.append(tuple(cur))
            return
        for i in range(s, len(vals)):
            v = vals[i]
            if cur and v - cur[-1] < sep:
                continue
            rec(i + 1, cur + [v])

    rec(0, [])
    return out


def solve_kick_path(B, lasers, pair_diffs):
    """lasers 按编号排序；pair_diffs[(i,j)] = d 表示 c[lasers[i]]-c[lasers[j]]=d
       实际上我们存 pair_diffs[(a,b)] with a<b as c[a]-c[b]=d.
       返回 {laser: code} 或 None。"""
    # Union-Find with offset: c[x] = c[root] + off[x]
    parent = {l: l for l in lasers}
    off = {l: 0 for l in lasers}

    def find(x):
        if parent[x] != x:
            r, o = find(parent[x])
            parent[x] = r
            off[x] += o
            return parent[x], off[x]
        return x, 0

    def unite(a, b, w):
        """c[a] - c[b] = w"""
        ra, oa = find(a)
        rb, ob = find(b)
        # c[a]=c[ra]+oa, c[b]=c[rb]+ob
        # c[ra]+oa - c[rb]-ob = w => c[ra]-c[rb] = w - oa + ob
        if ra == rb:
            return (oa - ob) == w
        # attach ra under rb: c[ra] = c[rb] + delta, delta = w - oa + ob
        parent[ra] = rb
        off[ra] = w - oa + ob
        return True

    for (a, b), d in pair_diffs.items():
        if not unite(a, b, d):
            return None

    comps = {}
    for l in lasers:
        r, o = find(l)
        comps.setdefault(r, []).append((l, o))

    codes = {}
    for r, mem in comps.items():
        offs = [o for _, o in mem]
        lo, hi = -min(offs), B - max(offs)
        if lo > hi:
            return None
        rv = (lo + hi) // 2
        for l, o in mem:
            codes[l] = rv + o
    return codes


def try_group(B, sep, ks, lasers, quads, rng):
    """尝试一次：返回 {(l,k):c} 或 None。"""
    if not quads:
        return None
    pairs = [(a, b) for i, a in enumerate(lasers)
             for b in lasers[i + 1:] if abs(a - b) <= GAP]
    nK = len(ks)
    # 每对抽码差
    diff_map = {}  # (a,b) -> list d per kick index
    for a, b in pairs:
        q = quads[int(rng.integers(0, len(quads)))]
        perm = rng.permutation(4)
        diff_map[(a, b)] = [q[int(p)] for p in perm[:nK]]

    # 逐 kick 解绝对码
    codes = {}
    for ki, k in enumerate(ks):
        pd = {(a, b): diff_map[(a, b)][ki] for a, b in pairs}
        sol = solve_kick_path(B, lasers, pd)
        if sol is None:
            return None
        for l, c in sol.items():
            codes[(l, k)] = c

    # 自身混叠：同一 l，c[k-1]-c[k]
    k_index = {k: i for i, k in enumerate(ks)}
    self_links = [(k_index[k - 1], ki) for ki, k in enumerate(ks) if k - 1 in k_index]
    if len(self_links) >= 2:
        for l in lasers:
            ds = sorted(codes[(l, ks[kp])] - codes[(l, ks[ki])] for kp, ki in self_links)
            if any(ds[t + 1] - ds[t] < sep for t in range(len(ds) - 1)):
                return None

    # 跨 kick 串扰
    if len(self_links) >= 2:
        for a, b in pairs:
            ds = sorted(codes[(a, ks[kp])] - codes[(b, ks[ki])] for kp, ki in self_links)
            if any(ds[t + 1] - ds[t] < sep for t in range(len(ds) - 1)):
                return None
            ds = sorted(codes[(b, ks[kp])] - codes[(a, ks[ki])] for kp, ki in self_links)
            if any(ds[t + 1] - ds[t] < sep for t in range(len(ds) - 1)):
                return None

    return codes


def search(B, sep, seconds=30.0, want=3, seed0=0):
    S.set_sep(sep)
    groups = laser_groups()
    quads = valid_quads(B, sep)
    print(f"预算 {B} ns | SEP={sep} | 合法4元组 {len(quads)} | "
          f"{len(groups)} 组 | 时限 {seconds:.0f}s", flush=True)
    if not quads:
        print("低于硬下界，无解")
        return []

    rng = np.random.default_rng(seed0)
    t0 = time.time()
    found, n_try = [], 0

    while time.time() - t0 < seconds and len(found) < want:
        n_try += 1
        codes = {}
        ok = True
        for ks, ls in groups:
            # 每组多抽几次
            sol = None
            for _ in range(80):
                sol = try_group(B, sep, ks, ls, quads, rng)
                if sol is not None:
                    break
            if sol is None:
                ok = False
                break
            codes.update(sol)
        if not ok:
            if n_try % 100 == 0:
                print(f"  [{time.time()-t0:5.1f}s] 尝试 {n_try}，命中 {len(found)}",
                      flush=True)
            continue

        u = np.zeros(len(S.VARS), dtype=np.int32)
        for (l, k), vi in VIDX.items():
            u[vi] = codes[(l, k)]
        cost = S.total_cost(u)
        if cost == 0:
            found.append(u.copy())
            print(f"  [{time.time()-t0:5.1f}s] ★ 第 {len(found)} 组 "
                  f"(尝试 {n_try}, max={int(u.max())})", flush=True)
        elif n_try % 100 == 0:
            print(f"  [{time.time()-t0:5.1f}s] 尝试 {n_try} cost={cost}", flush=True)

    print(f"结束：尝试 {n_try}，命中 {len(found)}，用时 {time.time()-t0:.1f}s",
          flush=True)
    return found


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--budget", type=int, required=True)
    ap.add_argument("--sep", type=int, default=12)
    ap.add_argument("--seconds", type=float, default=60.0)
    ap.add_argument("--want", type=int, default=3)
    ap.add_argument("--out", type=str, default=None)
    a = ap.parse_args()

    found = search(a.budget, a.sep, a.seconds, a.want)
    if not found:
        print("未找到。")
        sys.exit(1)

    print(f"\n评估 {len(found)} 组 ...", flush=True)
    scored = []
    for u in found:
        t = S.evaluate(S.make_code_fn(u), 1.6, 20.0)
        scored.append((t["ga"], t["kill"], int(max(u)), u))
    scored.sort()
    best = scored[0][3]
    fn = S.make_code_fn(best)
    ok, fails = S.check_avoid_true(fn, a.sep)
    print(f"最优 max={int(max(best))} 避真={ok} 违规={len(fails)}")
    for ratio in (1.6, 2.5):
        t = S.evaluate(fn, ratio, step=2.0)
        print(f"  ratio={ratio}: {t['gb']}->{t['ga']} "
              f"({t['ga']/max(t['gb'],1):.3%}) kill {t['kill']}/{t['tb']}")
    out = a.out or f"tcode_table_v22_{a.budget}ns.py"
    S.dump_table(best, os.path.join(HERE, out), a.sep, a.budget,
                 note="v22 码差抽样 fast_search_v22.py")
    for l in LASER_IDS:
        print(f"  L{l:<3d} " + " ".join(
            f"K{k}={int(best[VIDX[(l,k)]])}" for k in KICKS_OF[l]))


if __name__ == "__main__":
    main()
