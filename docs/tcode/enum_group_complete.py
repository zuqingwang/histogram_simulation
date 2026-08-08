# -*- coding: utf-8 -*-
"""
enum_group_complete.py —— 单组完整 kick 回溯（字母表 5 档）

每组 4 kick × 400 合法赋值/kick，spread 剪枝后叶节点验 **组内局部 cost=0**
（56 类 + 40 避真对均不跨组，global cost = 四组局部 cost 之和；不可用 total_cost，
否则未赋值变量为 0 会误计他组约束）。

结果供 exhaust_24step_r15.py --skip-group-enum 笛卡尔积合并。

用法：
  python docs/tcode/enum_group_complete.py              # 4 组顺序跑
  python docs/tcode/enum_group_complete.py --group 0  # 只跑第 0 组
  python docs/tcode/enum_group_complete.py --group 0 1 2 3  # 多组
"""
from __future__ import annotations

import argparse
import itertools
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import fast_search_v22 as F
import solve_tcode as S

ALPH = (0, 24, 48, 72, 96)
ALPH_I = {v: i for i, v in enumerate(ALPH)}
SEP, GAP = 12, 2
CACHE = os.path.join(HERE, "exhaust_groups_r15_24step.json")


def rebuild():
    S.CROSSTALK_MAX_GAP = GAP
    S.set_sep(SEP)
    S.CLASSES = S.build_classes(gap=GAP)
    S.AVOID_PAIRS = S.build_avoid_true(gap=GAP)


def group_constraint_indices(groups):
    """按 laser_groups 划分各类/避真对索引（已证均不跨组）。"""
    gof = {}
    for gi, (ks, ls) in enumerate(groups):
        for l in ls:
            for k in ks:
                gof[S.VIDX[(l, k)]] = gi
    cls_g = [[] for _ in groups]
    for ci in range(len(S.CLASSES)):
        gs = {gof[ia] for ia, ib in S.CLASSES[ci][1]}
        cls_g[list(gs)[0]].append(ci)
    av_g = [[] for _ in groups]
    for pi, (ia, ib) in enumerate(S.AVOID_PAIRS):
        gs = {gof[ia], gof[ib]}
        av_g[list(gs)[0]].append(pi)
    return cls_g, av_g


def group_local_cost(u, gi, cls_g, av_g):
    c = sum(S.class_cost(u, ci) for ci in cls_g[gi])
    c += sum(S.avoid_cost_one(u, pi) for pi in av_g[gi])
    return c


def per_kick_tuples(lasers):
    """每个 kick：400 个 (v0,v1,v2,v3) 按 lasers 顺序。"""
    n = len(lasers)
    out = []
    for vals in itertools.product(range(5), repeat=n):
        ok = True
        for i, a in enumerate(lasers):
            for j, b in enumerate(lasers):
                if i >= j or abs(a - b) > GAP:
                    continue
                if abs(ALPH[vals[i]] - ALPH[vals[j]]) < SEP:
                    ok = False
                    break
            if not ok:
                break
        if ok:
            out.append(vals)
    return out


def setup_group(ks, lasers):
    pk = per_kick_tuples(lasers)
    kix = {k: i for i, k in enumerate(ks)}
    sl = tuple((kix[k - 1], i) for i, k in enumerate(ks) if k - 1 in kix)
    pr = tuple((a, b) for i, a in enumerate(lasers)
               for b in lasers[i + 1:] if abs(a - b) <= GAP)
    vi = {ki: [S.VIDX[(lasers[j], ks[ki])] for j in range(len(lasers))]
          for ki in range(len(ks))}
    return pk, sl, pr, vi


def spread_ok_vals(acc, ks, lasers, sl, pr, mki):
    n = len(lasers)
    for l_idx in range(n):
        ds = []
        for kp, ki in sl:
            if ki <= mki:
                ds.append(ALPH[acc[kp][l_idx]] - ALPH[acc[ki][l_idx]])
        if len(ds) >= 2:
            ds.sort()
            for t in range(len(ds) - 1):
                if ds[t + 1] - ds[t] < SEP:
                    return False
    for a, b in pr:
        ai, bi = lasers.index(a), lasers.index(b)
        for x, y in ((ai, bi), (bi, ai)):
            ds = []
            for kp, ki in sl:
                if ki <= mki:
                    ds.append(ALPH[acc[kp][x]] - ALPH[acc[ki][y]])
            if len(ds) >= 2:
                ds.sort()
                for t in range(len(ds) - 1):
                    if ds[t + 1] - ds[t] < SEP:
                        return False
    return True


def enumerate_group(gi, ks, lasers, cls_g, av_g, progress_every=50000):
    pk, sl, pr, vi = setup_group(ks, lasers)
    print(f"  kick赋值/kick={len(pk)}  self_links={len(sl)}  pairs={len(pr)}", flush=True)

    feas_keys = set()
    leaves = nodes = 0
    t0 = time.time()
    acc = [[0] * len(lasers) for _ in ks]

    def bt(ki):
        nonlocal leaves, nodes
        if ki == len(ks):
            leaves += 1
            u = np.zeros(len(S.VARS), dtype=np.int32)
            for kki, k in enumerate(ks):
                for j, l in enumerate(lasers):
                    u[vi[kki][j]] = ALPH[acc[kki][j]]
            if group_local_cost(u, gi, cls_g, av_g) == 0:
                key = tuple(ALPH[acc[kki][j]]
                            for kki in range(len(ks)) for j in range(len(lasers)))
                feas_keys.add(key)
            if leaves % progress_every == 0:
                print(f"    叶={leaves:,} 可行={len(feas_keys):,}  "
                      f"节点={nodes:,}  [{time.time()-t0:.0f}s]", flush=True)
            return
        for tup in pk:
            nodes += 1
            for j, ai in enumerate(tup):
                acc[ki][j] = ai
            if spread_ok_vals(acc, ks, lasers, sl, pr, ki):
                bt(ki + 1)

    bt(0)
    elapsed = time.time() - t0
    print(f"  ★ 组{lasers} 叶={leaves:,} 可行={len(feas_keys):,}  "
          f"节点={nodes:,}  {elapsed:.1f}s", flush=True)

    sols = []
    for key in feas_keys:
        d = {}
        idx = 0
        for kki, k in enumerate(ks):
            for l in lasers:
                d[(l, k)] = key[idx]
                idx += 1
        sols.append(d)
    return sols


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--group", type=int, nargs="*", default=None,
                    help="组编号 0..3，默认全部")
    ap.add_argument("--cache", default=None, help="默认 exhaust_groups_r15_24step.json 或 per-group 后缀")
    ap.add_argument("--merge-only", action="store_true",
                    help="只合并已有 cache 并跑全局验收")
    a = ap.parse_args()

    groups = [(list(ks), list(ls)) for ks, ls in F.laser_groups()]
    rebuild()
    cls_g, av_g = group_constraint_indices(groups)

    if a.merge_only:
        if not os.path.isfile(a.cache):
            print(f"无缓存 {a.cache}")
            sys.exit(1)
        print(f"缓存已有，请运行: python docs/tcode/exhaust_24step_r15.py --skip-group-enum")
        sys.exit(0)

    gids = a.group if a.group is not None else list(range(len(groups)))

    print("=" * 64)
    print(f"kick 完整回溯  alphabet={list(ALPH)}  groups={gids}")
    print("=" * 64, flush=True)

    for gi in gids:
        cache_path = a.cache or (
            CACHE if len(gids) == len(groups) else os.path.join(HERE, f"exhaust_groups_r15_g{gi}.json"))
        if len(gids) == 1 and a.cache is None:
            cache_path = os.path.join(HERE, f"exhaust_groups_r15_g{gi}.json")

        ks, ls = groups[gi]
        print(f"\n[组 {gi}] lasers={ls} kicks={tuple(ks)}  cache={cache_path}", flush=True)
        sols = enumerate_group(gi, ks, ls, cls_g, av_g)
        with open(cache_path, "w", encoding="utf-8") as f:
            json.dump({
                "group": gi,
                "lasers": ls,
                "kicks": ks,
                "alphabet": list(ALPH),
                "solutions": [list(map(list, s.items())) for s in sols],
                "count": len(sols),
            }, f, ensure_ascii=False)
        print(f"  已写 {cache_path}  count={len(sols)}", flush=True)

    # 合并 per-group 缓存
    merged = [[] for _ in groups]
    for gi in range(len(groups)):
        p = os.path.join(HERE, f"exhaust_groups_r15_g{gi}.json")
        if os.path.isfile(p):
            with open(p, encoding="utf-8") as f:
                merged[gi] = [dict((tuple(k), v) for k, v in s)
                              for s in json.load(f)["solutions"]]
    if all(merged[gi] for gi in range(len(groups))):
        with open(CACHE, "w", encoding="utf-8") as f:
            json.dump({
                "alphabet": list(ALPH),
                "groups": [[list(map(list, s.items())) for s in g] for g in merged],
                "counts": [len(g) for g in merged],
            }, f, ensure_ascii=False)
        print(f"\n★ 已合并 → {CACHE}  counts={[len(g) for g in merged]}", flush=True)
        print("下一步: python docs/tcode/exhaust_24step_r15.py --skip-group-enum", flush=True)


if __name__ == "__main__":
    main()
