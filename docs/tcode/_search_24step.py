# -*- coding: utf-8 -*-
"""
_search_24step.py —— 字母表 24*(0..4) = {0,24,48,72,96} 搜 ratio=1.5 零残留

用法：
  python docs/tcode/_search_24step.py --minutes 15 --jobs 8
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

SEP, GAP, RATIO = 12, 2, 1.5
ALPH = np.array([0, 24, 48, 72, 96], dtype=np.int32)


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


def solve_alphabet(alphabet, rng, max_steps=40000, plateau=2000):
    alphabet = np.asarray(alphabet, dtype=np.int32)
    u = alphabet[rng.integers(0, len(alphabet), len(S.VARS))].astype(np.int32)
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
            av = 0
            for pi in touched_a:
                ia, ib = S.AVOID_PAIRS[pi]
                other = int(u[ib if ia == vi else ia])
                av += int(abs(int(val) - other) < SEP)
            c = (base
                 + sum(S._class_cost_with(u, ci, vi, int(val)) for ci in touched_c)
                 + av)
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


def residual(u, step):
    fn = S.make_code_fn(u)
    fr = G.build_firings(fn)
    ga = gb = kill = tb = 0
    for D in np.arange(5.0, 601.0, step):
        r = G.simulate(D, fr, max_gap=GAP, ratio=RATIO)
        ga += r["ga"]; gb += r["gb"]; kill += r["kill"]; tb += r["tb"]
        if ga or kill:
            return ga, gb, kill, tb, False
    return ga, gb, kill, tb, True


def worker(payload):
    seed0, seconds = payload
    sys.path.insert(0, HERE)
    import solve_tcode as Sloc
    import gen_tcode_figures as Gloc
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

    alph = np.array([0, 24, 48, 72, 96], dtype=np.int32)
    t0 = time.time()
    n_feas = n_try = 0
    best_ga = None

    def solve_local(rng):
        u = alph[rng.integers(0, len(alph), len(Sloc.VARS))].astype(np.int32)
        cost = Sloc.total_cost(u)
        best_cost, stale = cost, 0
        for _ in range(40000):
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
            for val in alph:
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
                if stale >= 2000:
                    return None
        return None

    def resid(u, step):
        fn = Sloc.make_code_fn(u)
        fr = Gloc.build_firings(fn)
        ga = gb = kill = tb = 0
        for D in np.arange(5.0, 601.0, step):
            r = Gloc.simulate(D, fr, max_gap=GAP, ratio=RATIO)
            ga += r["ga"]; gb += r["gb"]; kill += r["kill"]; tb += r["tb"]
            if ga or kill:
                return ga, gb, kill, tb, False
        return ga, gb, kill, tb, True

    seed = seed0
    seen_feas = set()  # 避免完全相同可行解反复验收
    while time.time() - t0 < seconds:
        n_try += 1
        cur_seed = seed
        seed += 1
        u = solve_local(np.random.default_rng(cur_seed))
        if u is None:
            continue
        key = tuple(int(x) for x in u)
        if key in seen_feas:
            continue
        seen_feas.add(key)
        n_feas += 1
        ga, gb, kill, tb, clean = resid(u, 10.0)
        if best_ga is None or ga < best_ga:
            best_ga = ga
        if not clean:
            if n_feas % 100 == 0:
                print(f"  worker={seed0} 当前seed={cur_seed} 可行{n_feas}/{n_try} "
                      f"唯一可行={len(seen_feas)} 最好残留={best_ga} "
                      f"[{time.time()-t0:.0f}s]", flush=True)
            continue
        ga, gb, kill, tb, clean = resid(u, 2.0)
        if clean:
            return {"ok": True, "u": [int(x) for x in u], "seed": cur_seed,
                    "worker": seed0, "feas": n_feas, "tries": n_try,
                    "unique": len(seen_feas), "secs": time.time() - t0,
                    "gb": gb, "tb": tb}
    return {"ok": False, "worker": seed0, "seed_end": seed - 1,
            "feas": n_feas, "tries": n_try, "unique": len(seen_feas),
            "secs": time.time() - t0, "best_ga": best_ga}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--minutes", type=float, default=15.0)
    ap.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    a = ap.parse_args()

    rebuild()
    print(f"字母表 24*(0..4) = {list(ALPH)}  ratio={RATIO}  SEP={SEP}  gap={GAP}",
          flush=True)
    print(f"  散开类={len(S.CLASSES)}  避真={len(S.AVOID_PAIRS)}  "
          f"{a.jobs} 进程 × {a.minutes:.0f} 分钟\n", flush=True)

    seconds = a.minutes * 60.0
    # 每个 worker 独占 1e9 的 seed 段，进程内 seed0, seed0+1, ... 绝不重叠
    seeds = [int(1_000_000_000 * (k + 1)) for k in range(a.jobs)]
    print(f"  worker 起点 seed = {seeds}", flush=True)
    t0 = time.time()
    found = None
    total_feas = 0

    with mp.Pool(a.jobs) as pool:
        for res in pool.imap_unordered(worker, [(s, seconds) for s in seeds]):
            total_feas += res["feas"]
            tag = "★零残留" if res["ok"] else "结束"
            if res["ok"]:
                print(f"  [{time.time()-t0:6.1f}s] {tag} worker={res['worker']} "
                      f"命中seed={res['seed']} 可行{res['feas']}/{res['tries']} "
                      f"唯一={res['unique']} ({res['secs']:.0f}s)", flush=True)
            else:
                print(f"  [{time.time()-t0:6.1f}s] {tag} worker={res['worker']} "
                      f"seed→{res['seed_end']} 可行{res['feas']}/{res['tries']} "
                      f"唯一={res['unique']} 最好残留={res.get('best_ga')} "
                      f"({res['secs']:.0f}s)", flush=True)
            if res["ok"]:
                found = res
                pool.terminate()
                break

    print("\n" + "=" * 64, flush=True)
    print(f"合计可行≈{total_feas}  用时={time.time()-t0:.0f}s", flush=True)
    if found is None:
        print("本轮未找到零残留解。", flush=True)
        sys.exit(2)

    u = np.array(found["u"], dtype=np.int32)
    path = os.path.join(HERE, "tcode_table_v40_r1.5_L5_24step_96ns.py")
    S.dump_table(u, path, SEP, 96,
                 note=f"v40 alphabet 24*(0..4); ratio={RATIO} zero")
    with open(path, "a", encoding="utf-8") as f:
        f.write("\nTCODE_ALPHABET=[0, 24, 48, 72, 96]\nTCODE_N_LEVELS=5\n")
    print(f"★ 写入 {path}", flush=True)
    print(f"  seed={found['seed']}  使用值={sorted(set(map(int, u)))}", flush=True)
    for l in S.LASER_IDS:
        print(f"  L{l:<2d} " + " ".join(
            f"K{k}={int(u[S.VIDX[(l,k)]]):>2d}" for k in S.KICKS_OF[l]))
    sys.exit(0)


if __name__ == "__main__":
    mp.freeze_support()
    main()
