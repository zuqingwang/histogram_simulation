# -*- coding: utf-8 -*-
"""
_alphabet_search.py —— 在离散字母表上搜零残留 tcode

v40 新约束：
  可选编码值只有 M 个（默认 5），均匀铺在 [0, budget]（默认 100ns）。
  档距大 → 展宽后不易糊；总跨度有界 → 鬼影贴真峰邻域。

用法：
  python docs/tcode/_alphabet_search.py --budget 100 --levels 5 --ratio 1.5 --minutes 10
"""
import argparse
import os
import sys
import time
import multiprocessing as mp

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import solve_tcode as S
import gen_tcode_figures as G

SEP, GAP = 12, 2


def make_alphabet(budget, levels):
    return np.asarray(np.round(np.linspace(0, budget, levels)), dtype=np.int32)


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


def solve_alphabet(alphabet, rng, max_steps=40000, plateau=2000):
    """min-conflicts，候选值只取自 alphabet。"""
    alphabet = np.asarray(alphabet, dtype=np.int32)
    n = len(S.VARS)
    u = alphabet[rng.integers(0, len(alphabet), n)].astype(np.int32)
    cost = S.total_cost(u)
    best_cost, stale = cost, 0

    for _ in range(max_steps):
        if cost == 0:
            return u
        bad_c = [ci for ci in range(len(S.CLASSES)) if S.class_cost(u, ci) > 0]
        bad_a = [pi for pi in range(len(S.AVOID_PAIRS)) if S.avoid_cost_one(u, pi)]
        if not bad_c and not bad_a:
            return u

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
                - sum(S.class_cost(u, ci) for ci in touched_c)
                - sum(S.avoid_cost_one(u, pi) for pi in touched_a))

        best, bestc = [], 1 << 30
        old = int(u[vi])
        for val in alphabet:
            c = (base
                 + sum(S._class_cost_with(u, ci, vi, int(val)) for ci in touched_c)
                 + sum(int(abs(int(val) - int(u[ib if ia == vi else ia])) < S.SEP)
                       for pi in touched_a
                       for ia, ib in [S.AVOID_PAIRS[pi]]))
            if c < bestc:
                best, bestc = [int(val)], c
            elif c == bestc:
                best.append(int(val))

        pick = old if (bestc == cost and rng.random() < 0.2) else int(best[rng.integers(0, len(best))])
        u[vi] = pick
        cost = (base
                + sum(S.class_cost(u, ci) for ci in touched_c)
                + sum(S.avoid_cost_one(u, pi) for pi in touched_a))
        if cost < best_cost:
            best_cost, stale = cost, 0
        else:
            stale += 1
            if stale >= plateau:
                return None
    return None


def _worker(payload):
    seeds, alphabet, ratio, seconds = payload
    sys.path.insert(0, HERE)
    import solve_tcode as Sloc
    import gen_tcode_figures as Gloc
    # rebuild locally for spawn
    Sloc.CROSSTALK_MAX_GAP = GAP
    Sloc.set_sep(SEP)
    Sloc.CLASSES = Sloc.build_classes(gap=GAP)
    Sloc.AVOID_PAIRS = Sloc.build_avoid_true(gap=GAP)
    Sloc.VAR2CLS = [[] for _ in Sloc.VARS]
    for ci, (_, prs) in enumerate(Sloc.CLASSES):
        for ia, ib in prs:
            Sloc.VAR2CLS[ia].append(ci)
            Sloc.VAR2CLS[ib].append(ci)
    Sloc.VAR2CLS = [sorted(set(c)) for c in Sloc.VAR2CLS]
    Sloc.VAR2AVOID = [[] for _ in Sloc.VARS]
    for pi, (ia, ib) in enumerate(Sloc.AVOID_PAIRS):
        Sloc.VAR2AVOID[ia].append(pi)
        Sloc.VAR2AVOID[ib].append(pi)

    alphabet = np.asarray(alphabet, dtype=np.int32)
    t0 = time.time()
    n_feas = n_try = 0
    best_ga = None

    def solve_local(rng):
        n = len(Sloc.VARS)
        u = alphabet[rng.integers(0, len(alphabet), n)].astype(np.int32)
        cost = Sloc.total_cost(u)
        best_cost, stale = cost, 0
        for _ in range(30000):
            if cost == 0:
                return u
            bad_c = [ci for ci in range(len(Sloc.CLASSES)) if Sloc.class_cost(u, ci) > 0]
            bad_a = [pi for pi in range(len(Sloc.AVOID_PAIRS)) if Sloc.avoid_cost_one(u, pi)]
            if not bad_c and not bad_a:
                return u
            if bad_a and (not bad_c or rng.random() < 0.45):
                pi = int(bad_a[rng.integers(0, len(bad_a))])
                ia, ib = Sloc.AVOID_PAIRS[pi]
                vi = int(ia if rng.random() < 0.5 else ib)
            else:
                prs = Sloc.CLASSES[int(rng.choice(bad_c))][1]
                vi = int(prs[int(rng.integers(0, len(prs)))][int(rng.integers(0, 2))])
            touched_c = Sloc.VAR2CLS[vi]
            touched_a = Sloc.VAR2AVOID[vi]
            base = (cost
                    - sum(Sloc.class_cost(u, ci) for ci in touched_c)
                    - sum(Sloc.avoid_cost_one(u, pi) for pi in touched_a))
            best, bestc = [], 1 << 30
            old = int(u[vi])
            for val in alphabet:
                av = 0
                for pi in touched_a:
                    ia, ib = Sloc.AVOID_PAIRS[pi]
                    other = int(u[ib if ia == vi else ia])
                    av += int(abs(int(val) - other) < SEP)
                c = (base
                     + sum(Sloc._class_cost_with(u, ci, vi, int(val)) for ci in touched_c)
                     + av)
                if c < bestc:
                    best, bestc = [int(val)], c
                elif c == bestc:
                    best.append(int(val))
            pick = old if (bestc == cost and rng.random() < 0.2) else int(best[rng.integers(0, len(best))])
            u[vi] = pick
            cost = (base
                    + sum(Sloc.class_cost(u, ci) for ci in touched_c)
                    + sum(Sloc.avoid_cost_one(u, pi) for pi in touched_a))
            if cost < best_cost:
                best_cost, stale = cost, 0
            else:
                stale += 1
                if stale >= 1500:
                    return None
        return None

    def resid(u, step):
        fn = Sloc.make_code_fn(u)
        fr = Gloc.build_firings(fn)
        ga = gb = kill = tb = 0
        for D in np.arange(5.0, 601.0, step):
            r = Gloc.simulate(D, fr, max_gap=GAP, ratio=ratio)
            ga += r["ga"]; gb += r["gb"]; kill += r["kill"]; tb += r["tb"]
            if ga or kill:
                return ga, gb, kill, tb, False
        return ga, gb, kill, tb, True

    for seed in seeds:
        if time.time() - t0 > seconds:
            break
        n_try += 1
        u = solve_local(np.random.default_rng(seed))
        if u is None:
            continue
        n_feas += 1
        ga, gb, kill, tb, clean = resid(u, 10.0)
        if best_ga is None or ga < best_ga:
            best_ga = ga
        if not clean:
            continue
        ga, gb, kill, tb, clean = resid(u, 2.0)
        if clean:
            return {"ok": True, "u": [int(x) for x in u], "seed": seed,
                    "feas": n_feas, "tries": n_try, "secs": time.time() - t0,
                    "gb": gb, "tb": tb, "best_ga": 0}

    return {"ok": False, "u": None, "seed": seeds[0],
            "feas": n_feas, "tries": n_try, "secs": time.time() - t0,
            "best_ga": best_ga}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--budget", type=int, default=100)
    ap.add_argument("--levels", type=int, default=5)
    ap.add_argument("--ratio", type=float, default=1.5)
    ap.add_argument("--minutes", type=float, default=8.0)
    ap.add_argument("--jobs", type=int,
                    default=max(1, (os.cpu_count() or 2) - 1))
    a = ap.parse_args()

    rebuild()
    alphabet = make_alphabet(a.budget, a.levels)
    print(f"离散字母表搜索 · ratio={a.ratio} SEP={SEP} gap={GAP}", flush=True)
    print(f"  budget={a.budget}ns  levels={a.levels}  alphabet={list(alphabet)}",
          flush=True)
    print(f"  档距={int(alphabet[1]-alphabet[0]) if len(alphabet)>1 else 0}ns  "
          f"散开类={len(S.CLASSES)} 避真={len(S.AVOID_PAIRS)}", flush=True)
    print(f"  {a.jobs} 进程 × {a.minutes:.0f} 分钟\n", flush=True)

    seconds = a.minutes * 60.0
    seeds_per = 200
    # 多轮派发，每 worker 一批 seed
    t0 = time.time()
    total_feas = 0
    best_ga = None
    found = None
    seed0 = 30000 + a.budget * 17 + a.levels * 3

    with mp.Pool(a.jobs) as pool:
        round_id = 0
        while time.time() - t0 < seconds and found is None:
            remain = seconds - (time.time() - t0)
            payloads = []
            for w in range(a.jobs):
                seeds = list(range(seed0 + w * seeds_per,
                                   seed0 + (w + 1) * seeds_per))
                payloads.append((seeds, alphabet.tolist(), a.ratio, remain))
            seed0 += a.jobs * seeds_per
            round_id += 1
            for res in pool.imap_unordered(_worker, payloads):
                total_feas += res["feas"]
                if res["best_ga"] is not None:
                    if best_ga is None or res["best_ga"] < best_ga:
                        best_ga = res["best_ga"]
                tag = "★零残留" if res["ok"] else "批次结束"
                print(f"  [{time.time()-t0:6.1f}s] {tag} seed0={res['seed']} "
                      f"可行{res['feas']}/{res['tries']} "
                      f"最好残留={res['best_ga']} ({res['secs']:.0f}s)",
                      flush=True)
                if res["ok"]:
                    found = res
                    pool.terminate()
                    break
            if found:
                break
            # 单轮若很快结束，继续下一轮
            if time.time() - t0 >= seconds:
                break

    print("\n" + "=" * 64, flush=True)
    print(f"合计可行≈{total_feas}  最好残留={best_ga}", flush=True)
    if found is None:
        print("未找到零残留解。可增大 levels 或 minutes。", flush=True)
        sys.exit(2)

    u = np.array(found["u"], dtype=np.int32)
    name = f"tcode_table_v40_r{a.ratio}_L{a.levels}_{a.budget}ns.py"
    out = os.path.join(HERE, name)
    note = (f"v40 离散字母表；levels={a.levels} alphabet={list(alphabet)} "
            f"ratio={a.ratio} gap={GAP} SEP={SEP}")
    S.dump_table(u, out, SEP, a.budget, note=note)
    # 附加字母表常量
    with open(out, "a", encoding="utf-8") as f:
        f.write(f"\nTCODE_ALPHABET = {list(map(int, alphabet))}\n")
        f.write(f"TCODE_N_LEVELS = {a.levels}\n")
    print(f"★ 写入 {out}", flush=True)
    print(f"  max={int(u.max())}  使用值={sorted(set(map(int, u)))}", flush=True)
    for l in S.LASER_IDS:
        print(f"  L{l:<2d} " + " ".join(
            f"K{k}={int(u[S.VIDX[(l,k)]]):>3d}" for k in S.KICKS_OF[l]))
    sys.exit(0)


if __name__ == "__main__":
    mp.freeze_support()
    main()
