# -*- coding: utf-8 -*-
"""
_zero_fast.py —— 用「快搜」批量产可行解，再筛零残留

快搜来源（对话里已验证）：
  1) fast_search_v22：码差抽样 + 差分求解（24~36ns 秒级）
  2) search_v22 风格：多进程 + 大量短重启（min-conflicts 重尾，短重启才快）

判据顺序：约束可行 → 粗筛零残留 → 细扫零残留 → 再压预算。

用法：
  python docs/tcode/_zero_fast.py --ratio 1.5 --budgets 36 32 28 --seconds 60
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
import fast_search_v22 as F

SEP = 12
GAP = 2


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


def _short_restart_batch(payload):
    """worker：一批短重启，找到一个可行解就返回。"""
    seeds, budget, sep, steps, plateau = payload
    import solve_tcode as Sloc
    Sloc.set_sep(sep)
    Sloc.CROSSTALK_MAX_GAP = GAP
    # 约束表已在父进程建好并通过 fork/spawn 带过来；spawn 时需重建
    if not getattr(Sloc, "CLASSES", None) or len(Sloc.CLASSES) == 0:
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
    for seed in seeds:
        u = Sloc.solve(budget, np.random.default_rng(seed),
                       max_steps=steps, plateau=plateau)
        if u is not None:
            return [int(x) for x in u]
    return None


def hunt(B, ratio, seconds, jobs, want_zero=1):
    """在预算 B 内找零残留解。返回 (u, 说明)。"""
    rebuild()
    groups = F.laser_groups()
    quads = F.valid_quads(B, SEP)
    rng = np.random.default_rng(B * 7919 + 17)
    t0 = time.time()
    n_feas = 0
    best_ga = None
    zeros = []

    # ---- 通道 A：码差抽样（主通道，紧预算秒级）----
    while time.time() - t0 < seconds * 0.7 and len(zeros) < want_zero:
        if not quads:
            break
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
        ga, gb, kill, tb, clean = residual(u, ratio, 10.0)
        if not clean:
            if best_ga is None or ga < best_ga:
                best_ga = ga
            continue
        ga, gb, kill, tb, clean = residual(u, ratio, 2.0)
        if clean:
            zeros.append(u.copy())
            print(f"    [{time.time()-t0:5.1f}s] ★零残留#{len(zeros)} "
                  f"(码差抽样，可行{n_feas}，鬼{gb})", flush=True)

    # ---- 通道 B：多进程短重启（补漏）----
    if len(zeros) < want_zero and time.time() - t0 < seconds:
        remain = seconds - (time.time() - t0)
        seed0 = 100000 + B * 97
        per = 12
        with mp.Pool(jobs) as pool:
            while time.time() - t0 < seconds and len(zeros) < want_zero:
                payloads = []
                for w in range(jobs):
                    seeds = list(range(seed0 + w * per,
                                       seed0 + (w + 1) * per))
                    payloads.append((seeds, B, SEP, 5000, 800))
                seed0 += jobs * per
                for u_list in pool.imap_unordered(_short_restart_batch, payloads):
                    if u_list is None:
                        continue
                    u = np.array(u_list, dtype=np.int32)
                    n_feas += 1
                    ga, gb, kill, tb, clean = residual(u, ratio, 10.0)
                    if not clean:
                        if best_ga is None or ga < best_ga:
                            best_ga = ga
                        continue
                    ga, gb, kill, tb, clean = residual(u, ratio, 2.0)
                    if clean:
                        zeros.append(u.copy())
                        print(f"    [{time.time()-t0:5.1f}s] ★零残留#{len(zeros)} "
                              f"(短重启，可行{n_feas}，鬼{gb})", flush=True)
                        if len(zeros) >= want_zero:
                            break
                if time.time() - t0 > seconds:
                    break

    if zeros:
        best = min(zeros, key=lambda u: int(u.max()))
        return best, (f"可行{n_feas}，零残留{len(zeros)}组，"
                      f"{time.time()-t0:.0f}s，max={int(best.max())}")
    return None, (f"可行{n_feas}，最好残留={best_ga}，{time.time()-t0:.0f}s")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ratio", type=float, default=1.5)
    ap.add_argument("--budgets", type=int, nargs="+",
                    default=[40, 36, 32, 28, 24])
    ap.add_argument("--seconds", type=float, default=90.0)
    ap.add_argument("--jobs", type=int,
                    default=max(1, (os.cpu_count() or 2) - 1))
    a = ap.parse_args()

    rebuild()
    print(f"零残留快搜 · ratio={a.ratio} SEP={SEP} gap={GAP} "
          f"jobs={a.jobs}", flush=True)
    print(f"散开类={len(S.CLASSES)} 避真对={len(S.AVOID_PAIRS)}\n", flush=True)

    best = None
    for B in a.budgets:
        print(f"  B={B}ns ...", flush=True)
        u, msg = hunt(B, a.ratio, a.seconds, a.jobs)
        if u is not None:
            print(f"  B={B}ns ★ {msg}", flush=True)
            best = (B, u.copy())
            name = f"tcode_table_zero_r{a.ratio}_{int(u.max())}ns.py"
            S.dump_table(u, os.path.join(HERE, name), SEP, int(u.max()),
                         note=f"零残留优先再压；ratio={a.ratio} gap={GAP} SEP={SEP} "
                              f"（_zero_fast）")
        else:
            print(f"  B={B}ns 失败  {msg}", flush=True)
            # 继续往下试（有时更低档反而好抽）

    if best is None:
        print("\n未找到零残留解")
        return
    B, u = best
    print(f"\n本轮最短零残留 = {B}ns（max={int(u.max())}）")
    for l in S.LASER_IDS:
        print(f"  L{l:<2d} " + " ".join(
            f"K{k}={int(u[S.VIDX[(l,k)]]):>2d}" for k in S.KICKS_OF[l]))


if __name__ == "__main__":
    # Windows spawn 需要
    mp.freeze_support()
    main()
