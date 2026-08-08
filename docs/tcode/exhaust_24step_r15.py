# -*- coding: utf-8 -*-
"""
exhaust_24step_r15.py —— {0,24,48,72,96} + ratio=1.5 系统枚举

结构：
  1. 4 个激光组约束互不跨组（56 类 / 40 避真对均组内；global cost = 四组局部之和）
  2. enum_group_complete.py 用 kick 回溯 + spread 剪枝收集组内局部 cost=0 赋值
  3. 4 组笛卡尔积（或阻塞式 min-conflicts 作备用）→ 1~600m 粗/细扫验收零残留

用法：
  python docs/tcode/exhaust_24step_r15.py
  python docs/tcode/exhaust_24step_r15.py --skip-group-enum --group-cache exhaust_groups.json
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import fast_search_v22 as F
import gen_tcode_figures as G
import solve_tcode as S

ALPH = np.array([0, 24, 48, 72, 96], dtype=np.int32)
SEP, GAP, BUDGET = 12, 2, 96


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


def group_var_indices(ks, lasers):
    return [S.VIDX[(l, k)] for l in lasers for k in ks]


def matches_blocked(u, blocked_list, indices):
    ui = u[indices]
    for b in blocked_list:
        if np.array_equal(ui, b):
            return True
    return False


def solve_group(ks, lasers, blocked_list, indices, rng, max_steps=80000, plateau=3000):
    """min-conflicts，解不等于任一 blocked 组赋值。"""
    u_full = ALPH[rng.integers(0, len(ALPH), len(S.VARS))].astype(np.int32)
    cost = S.total_cost(u_full)
    best_cost, stale = cost, 0

    for _ in range(max_steps):
        if cost == 0 and not matches_blocked(u_full, blocked_list, indices):
            return { (l, k): int(u_full[S.VIDX[(l, k)]])
                     for l in lasers for k in ks }
        bad_c = [ci for ci in range(len(S.CLASSES)) if S.class_cost(u_full, ci) > 0]
        bad_a = [pi for pi in range(len(S.AVOID_PAIRS)) if S.avoid_cost_one(u_full, pi)]
        if not bad_c and not bad_a:
            if not matches_blocked(u_full, blocked_list, indices):
                return { (l, k): int(u_full[S.VIDX[(l, k)]])
                         for l in lasers for k in ks }
            # 撞 blocked：随机翻一个组内变量
            vi = int(rng.choice(indices))
            u_full[vi] = ALPH[rng.integers(0, len(ALPH))]
            cost = S.total_cost(u_full)
            stale = 0
            continue
        if bad_a and (not bad_c or rng.random() < 0.45):
            pi = int(bad_a[rng.integers(0, len(bad_a))])
            ia, ib = S.AVOID_PAIRS[pi]
            vi = int(ia if rng.random() < 0.5 else ib)
        else:
            prs = S.CLASSES[int(rng.choice(bad_c))][1]
            vi = int(prs[int(rng.integers(0, len(prs)))][int(rng.integers(0, 2))])
        touched_c = S.VAR2CLS[vi]
        touched_a = S.VAR2AVOID[vi]
        base = (cost
                - sum(S.class_cost(u_full, ci) for ci in touched_c)
                - sum(S.avoid_cost_one(u_full, pi) for pi in touched_a))
        best, bestc = [], 1 << 30
        old = int(u_full[vi])
        for val in ALPH:
            av = sum(int(abs(int(val) - int(u_full[ib if ia == vi else ia])) < SEP)
                     for pi in touched_a for ia, ib in [S.AVOID_PAIRS[pi]])
            c = (base
                 + sum(S._class_cost_with(u_full, ci, vi, int(val)) for ci in touched_c)
                 + av)
            if c < bestc:
                best, bestc = [int(val)], c
            elif c == bestc:
                best.append(int(val))
        pick = old if (bestc == cost and rng.random() < 0.2) else int(best[rng.integers(0, len(best))])
        u_full[vi] = pick
        cost = (base
                + sum(S.class_cost(u_full, ci) for ci in touched_c)
                + sum(S.avoid_cost_one(u_full, pi) for pi in touched_a))
        if cost < best_cost:
            best_cost, stale = cost, 0
        else:
            stale += 1
            if stale >= plateau:
                return None
    return None


def enumerate_group_all(ks, lasers, seed0=0):
    indices = group_var_indices(ks, lasers)
    blocked, solutions = [], []
    rng = np.random.default_rng(seed0)
    t0 = time.time()
    tries = 0
    while True:
        tries += 1
        sol = solve_group(ks, lasers, blocked, indices, rng)
        if sol is None:
            print(f"    组{lasers} 阻塞枚举结束  唯一可行={len(solutions)}  "
                  f"尝试={tries}  [{time.time()-t0:.1f}s]", flush=True)
            break
        key = tuple(sorted(sol.items()))
        solutions.append(sol)
        blocked.append(np.array([sol[(l, k)] for l in lasers for k in ks], dtype=np.int32))
        if len(solutions) % 500 == 0:
            print(f"    组{lasers} 已收 {len(solutions)} … [{time.time()-t0:.0f}s]", flush=True)
    return solutions


def load_group_caches(cache_path):
    """加载合并缓存或四份 per-group 缓存。"""
    if os.path.isfile(cache_path):
        with open(cache_path, encoding="utf-8") as f:
            cached = json.load(f)
        return [[dict((tuple(k), v) for k, v in s.items()) for s in g]
                for g in cached["groups"]]
    merged = []
    for gi in range(4):
        p = os.path.join(HERE, f"exhaust_groups_r15_g{gi}.json")
        if not os.path.isfile(p):
            return None
        with open(p, encoding="utf-8") as f:
            merged.append([dict((tuple(k), v) for k, v in s)
                           for s in json.load(f)["solutions"]])
    return merged


def merge_to_u(codes_list):
    u = np.zeros(len(S.VARS), dtype=np.int32)
    for codes in codes_list:
        for (l, k), v in codes.items():
            u[S.VIDX[(l, k)]] = v
    return u


def check_one(u, ratio):
    ga, kill, clean = residual(u, ratio, 10.0)
    if not clean:
        return ga, kill, False
    return residual(u, ratio, 2.0)


def iterate_merge(all_group, ratio, mode, sample_n, rng):
    """mode: 'full' | 'sample'"""
    counts = [len(g) for g in all_group]
    prod = int(np.prod(counts))
    n_checked = best_ga = 0
    best_u = None
    t1 = time.time()

    if mode == "full":
        for c0 in all_group[0]:
            for c1 in all_group[1]:
                for c2 in all_group[2]:
                    for c3 in all_group[3]:
                        n_checked += 1
                        u = merge_to_u((c0, c1, c2, c3))
                        ga, kill, clean = check_one(u, ratio)
                        if best_ga is None or ga < best_ga:
                            best_ga, best_u = ga, u.copy()
                        if clean:
                            return n_checked, best_ga, best_u, u, True
                        if n_checked % max(1, prod // 20) == 0:
                            print(f"  全局 {n_checked:,}/{prod:,}  最好残留={best_ga}  "
                                  f"[{time.time()-t1:.0f}s]", flush=True)
        return n_checked, best_ga, best_u, None, False

    seen = set()
    while n_checked < sample_n:
        picks = tuple(int(rng.integers(0, c)) for c in counts)
        if picks in seen:
            continue
        seen.add(picks)
        n_checked += 1
        u = merge_to_u([all_group[i][picks[i]] for i in range(4)])
        ga, kill, clean = check_one(u, ratio)
        if best_ga is None or ga < best_ga:
            best_ga, best_u = ga, u.copy()
        if clean:
            return n_checked, best_ga, best_u, u, True
        if n_checked % max(1, sample_n // 20) == 0:
            print(f"  抽样 {n_checked:,}/{sample_n:,}  最好残留={best_ga}  "
                  f"[{time.time()-t1:.0f}s]", flush=True)
    return n_checked, best_ga, best_u, None, False


def residual(u, ratio, step):
    fn = S.make_code_fn(u)
    fr = G.build_firings(fn)
    ga = kill = 0
    for d in np.arange(5.0, 601.0, step):
        r = G.simulate(d, fr, max_gap=GAP, ratio=ratio)
        ga += r["ga"]
        kill += r["kill"]
        if ga or kill:
            return ga, kill, False
    return 0, 0, True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ratio", type=float, default=1.5)
    ap.add_argument("--out", default=os.path.join(HERE, "tcode_table_v40_r1.5_L5_24step_96ns.py"))
    ap.add_argument("--group-cache", default=os.path.join(HERE, "exhaust_groups_r15_24step.json"))
    ap.add_argument("--skip-group-enum", action="store_true")
    ap.add_argument("--max-full-merge", type=int, default=10_000_000,
                    help="笛卡尔积超过此值则改用随机抽样合并")
    ap.add_argument("--sample-merge", type=int, default=500_000,
                    help="抽样合并时检查的组合数")
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()

    rebuild()
    groups = [(list(ks), list(ls)) for ks, ls in F.laser_groups()]

    print("=" * 64)
    print(f"全局合并验收  alphabet={list(map(int, ALPH))}  ratio={a.ratio}")
    print("=" * 64, flush=True)

    t0 = time.time()
    if a.skip_group_enum:
        all_group = load_group_caches(a.group_cache)
        if all_group is None:
            print("未找到合并缓存或 exhaust_groups_r15_g{0..3}.json")
            sys.exit(1)
        print(f"已加载组缓存", flush=True)
    else:
        all_group = []
        for ks, ls in groups:
            print(f"\n阻塞枚举组 lasers={ls} kicks={tuple(ks)}", flush=True)
            all_group.append(enumerate_group_all(ks, ls))
        with open(a.group_cache, "w", encoding="utf-8") as f:
            json.dump({"alphabet": list(map(int, ALPH)),
                       "groups": [[list(map(list, s.items())) for s in g]
                                  for g in all_group]}, f, ensure_ascii=False)

    counts = [len(g) for g in all_group]
    prod = int(np.prod(counts))
    mode = "full" if prod <= a.max_full_merge else "sample"
    print(f"\n各组可行: {counts}  笛卡尔积={prod:,}  合并模式={mode}", flush=True)

    if any(c == 0 for c in counts):
        print("某组无可行解。")
        sys.exit(2)

    if mode == "sample":
        print(f"笛卡尔积过大，随机抽样 {a.sample_merge:,} 组做 1~600m 验收", flush=True)

    n_checked, best_ga, best_u, zero_u, found = iterate_merge(
        all_group, a.ratio, mode, a.sample_merge, np.random.default_rng(a.seed))

    if found:
        print(f"\n★ 零残留  已检={n_checked:,}", flush=True)
        S.dump_table(zero_u, a.out, SEP, BUDGET,
                     note=f"24*(0..4) ratio={a.ratio} exhaust")
        with open(a.out, "a", encoding="utf-8") as f:
            f.write("\nTCODE_ALPHABET=[0, 24, 48, 72, 96]\nTCODE_N_LEVELS=5\n")
        print(f"  写入 {a.out}")
        sys.exit(0)

    label = "全枚举" if mode == "full" else f"抽样({n_checked:,})"
    print(f"\n{label}完成：零残留=0  最好残留={best_ga}")
    print(f"总用时 {time.time()-t0:.1f}s")
    sys.exit(2)


if __name__ == "__main__":
    main()
