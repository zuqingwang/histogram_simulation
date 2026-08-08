# -*- coding: utf-8 -*-
"""
_check_32.py —— 长时间确认 ratio=1.5 在 32ns 是否存在零残留解

策略：fast_search_v22 码差抽样（紧预算最快）；多进程并行不同 seed。
只要出现一张细扫零残留就立刻落盘并退出成功码。

用法：
  python docs/tcode/_check_32.py --minutes 45 --jobs 8
"""
import argparse
import os
import sys
import time
import multiprocessing as mp

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

BUDGET = 32
SEP = 12
GAP = 2
RATIO = 1.5


def _worker(payload):
    """一个进程：用独立 seed 狂抽，找到零残留返回 u 列表，否则返回统计。"""
    seed, seconds, report_every = payload
    sys.path.insert(0, HERE)
    import solve_tcode as S
    import gen_tcode_figures as G
    import fast_search_v22 as F

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

    groups = F.laser_groups()
    quads = F.valid_quads(BUDGET, SEP)
    rng = np.random.default_rng(seed)
    t0 = time.time()
    n_try = n_feas = 0
    best_ga = None
    hist = {}  # residual count histogram

    def residual(u, step):
        fn = S.make_code_fn(u)
        fr = G.build_firings(fn)
        ga = kill = 0
        gb = tb = 0
        for D in np.arange(5.0, 601.0, step):
            r = G.simulate(D, fr, max_gap=GAP, ratio=RATIO)
            ga += r["ga"]; gb += r["gb"]; kill += r["kill"]; tb += r["tb"]
            if ga or kill:
                return ga, gb, kill, tb, False
        return ga, gb, kill, tb, True

    while time.time() - t0 < seconds:
        n_try += 1
        codes, ok = {}, True
        for ks, ls in groups:
            sol = None
            for _ in range(80):
                sol = F.try_group(BUDGET, SEP, ks, ls, quads, rng)
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
        ga, gb, kill, tb, clean = residual(u, 10.0)
        hist[ga] = hist.get(ga, 0) + 1
        if best_ga is None or ga < best_ga:
            best_ga = ga
        if clean:
            ga2, gb2, kill2, tb2, clean2 = residual(u, 2.0)
            if clean2:
                return {
                    "ok": True, "seed": seed, "u": [int(x) for x in u],
                    "feas": n_feas, "tries": n_try, "secs": time.time() - t0,
                    "gb": gb2, "tb": tb2, "best_ga": 0, "hist": hist,
                }
        if report_every and n_feas % report_every == 0:
            # 通过文件副作用不好；imap 只能返回最终结果。此处仅累计。
            pass

    return {
        "ok": False, "seed": seed, "u": None,
        "feas": n_feas, "tries": n_try, "secs": time.time() - t0,
        "best_ga": best_ga, "hist": hist,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--minutes", type=float, default=45.0)
    ap.add_argument("--jobs", type=int,
                    default=max(1, (os.cpu_count() or 2) - 1))
    a = ap.parse_args()

    seconds = a.minutes * 60.0
    print(f"检查 32ns 零残留 · ratio={RATIO} SEP={SEP} gap={GAP}", flush=True)
    print(f"  {a.jobs} 进程 × {a.minutes:.0f} 分钟 = 合计约 "
          f"{a.jobs * a.minutes:.0f} 进程·分钟", flush=True)
    print(f"  策略：fast_search 码差抽样 → 粗筛10m → 细扫2m\n", flush=True)

    # 错开 seed，每进程独立搜满整段时间；任一命中就提前结束
    seeds = [1000 + i * 9973 for i in range(a.jobs)]
    payloads = [(s, seconds, 50) for s in seeds]
    t0 = time.time()
    found = None
    total_feas = 0
    merged_hist = {}
    best_ga = None

    with mp.Pool(a.jobs) as pool:
        # 用 imap_unordered：谁先找到零残留就先返回
        for res in pool.imap_unordered(_worker, payloads):
            total_feas += res["feas"]
            for k, v in res["hist"].items():
                merged_hist[k] = merged_hist.get(k, 0) + v
            if res["best_ga"] is not None:
                if best_ga is None or res["best_ga"] < best_ga:
                    best_ga = res["best_ga"]
            tag = "★零残留" if res["ok"] else "结束"
            print(f"  [{time.time()-t0:6.1f}s] worker seed={res['seed']} "
                  f"{tag} 可行{res['feas']}/{res['tries']} "
                  f"最好残留={res['best_ga']} ({res['secs']:.0f}s)",
                  flush=True)
            if res["ok"]:
                found = res
                pool.terminate()
                break

    print("\n" + "=" * 64, flush=True)
    print(f"合计可行解 ≈ {total_feas}，残留直方图: "
          + ", ".join(f"{k}:{v}" for k, v in sorted(merged_hist.items())),
          flush=True)
    print(f"全局最好残留 = {best_ga}", flush=True)

    if found is None:
        print(f"\n结论：在约 {a.jobs}×{a.minutes:.0f} 进程·分钟内，"
              f"32ns / ratio=1.5 **未找到**零残留解。", flush=True)
        print("（不等于严格证明无解，但强烈支持「最短零残留仍是 36ns」。）",
              flush=True)
        sys.exit(2)

    import solve_tcode as S
    S.set_sep(SEP)
    u = np.array(found["u"], dtype=np.int32)
    out = os.path.join(HERE, f"tcode_table_zero_r{RATIO}_{BUDGET}ns.py")
    S.dump_table(u, out, SEP, BUDGET,
                 note=f"32ns 零残留确认；ratio={RATIO} gap={GAP} SEP={SEP}")
    print(f"\n结论：32ns ★存在零残留解！已写入 {out}", flush=True)
    for l in S.LASER_IDS:
        print(f"  L{l:<2d} " + " ".join(
            f"K{k}={int(u[S.VIDX[(l,k)]]):>2d}" for k in S.KICKS_OF[l]))
    sys.exit(0)


if __name__ == "__main__":
    mp.freeze_support()
    main()
