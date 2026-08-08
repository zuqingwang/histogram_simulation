# -*- coding: utf-8 -*-
"""
search_v22.py —— 多进程 + 大量短重启，搜 v22 码表（散开 + 避真峰）
====================================================================
为什么要短重启：
    min-conflicts 局部搜索的耗时是【重尾分布】—— 成功的那次通常 1000~2000 步
    就出解了；一旦陷进局部极小，再磨 10 万步也是白搭。
    所以「大量短重启」远好于「少量长跑」。

用法（PowerShell）：
    python docs/tcode/search_v22.py --budget 36
    python docs/tcode/search_v22.py --budget 32 --minutes 30 --jobs 8

参数：
    --budget   码预算 [ns]，必填
    --sep      落点最小间隔 / 避真间隔 [ns]，默认 12
    --steps    每次重启的最大步数，默认 4000（短！）
    --minutes  总时限 [分钟]，默认 10
    --jobs     并行进程数，默认 = CPU 核数 − 1
    --want     找到几组解就停，默认 5（多找几组好挑残留最低的）
    --out      输出文件名，默认 tcode_table_v22_<budget>ns.py

缩写：SEP（落点最小间隔，取回波峰宽 + 裕度）、
      XM（XtalkMark，串扰标记）、TOF（Time of Flight，飞行时间）。
"""
import argparse
import os
import sys
import time
import multiprocessing as mp

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)


def _try_seeds(payload):
    """一个 worker：用一批 seed 做短重启，找到就立刻返回。"""
    seeds, budget, sep, steps, plateau = payload
    import solve_tcode as S
    S.set_sep(sep)
    for seed in seeds:
        rng = np.random.default_rng(seed)
        u = S.solve(budget, rng, max_steps=steps, plateau=plateau)
        if u is not None:
            return seed, [int(x) for x in u]
    return None, None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--budget",  type=int, required=True)
    ap.add_argument("--sep",     type=int, default=12)
    ap.add_argument("--steps",   type=int, default=5000)
    ap.add_argument("--plateau", type=int, default=800,
                    help="连续这么多步不降代价就放弃本次重启")
    ap.add_argument("--minutes", type=float, default=10.0)
    ap.add_argument("--jobs",    type=int, default=max(1, (os.cpu_count() or 2) - 1))
    ap.add_argument("--want",    type=int, default=5)
    ap.add_argument("--out",     type=str, default=None)
    a = ap.parse_args()

    out_name = a.out or f"tcode_table_v22_{a.budget}ns.py"
    out_path = os.path.join(HERE, out_name)

    print(f"预算 {a.budget} ns | 间隔 {a.sep} ns | 每次重启 {a.steps} 步 "
          f"(plateau={a.plateau}) | {a.jobs} 进程 | 时限 {a.minutes:.0f} 分钟",
          flush=True)

    t0 = time.time()
    deadline = t0 + a.minutes * 60
    found, seed_base, n_tried = [], 0, 0
    per_worker = 8                      # 每个 worker 一次吃多个 seed

    with mp.Pool(a.jobs) as pool:
        while time.time() < deadline and len(found) < a.want:
            payloads = []
            for w in range(a.jobs):
                seeds = list(range(seed_base + w * per_worker,
                                   seed_base + (w + 1) * per_worker))
                payloads.append((seeds, a.budget, a.sep, a.steps, a.plateau))
            seed_base += a.jobs * per_worker
            n_tried += a.jobs * per_worker
            for seed, u in pool.imap_unordered(_try_seeds, payloads):
                if u is not None:
                    found.append(u)
                    print(f"  [{time.time()-t0:6.1f}s] 找到第 {len(found)} 组解 "
                          f"(seed={seed}, max={max(u)} ns)", flush=True)
                    if len(found) >= a.want:
                        break
            print(f"  [{time.time()-t0:6.1f}s] 已试约 {n_tried} 次重启，"
                  f"命中 {len(found)} 组", flush=True)

    if not found:
        print(f"\n时限内没找到。这【不等于无解】—— 理论下界是 {1.5*a.sep:.0f} ns。")
        print(f"可以加时间或加步数再试：--minutes 60 --steps 8000")
        sys.exit(1)

    # ---- 挑残留最低的一组 ----
    import solve_tcode as S
    S.set_sep(a.sep)
    print(f"\n评估 {len(found)} 组解（1~600m 粗扫）...", flush=True)
    scored = []
    for u in found:
        arr = np.array(u)
        t = S.evaluate(S.make_code_fn(arr), 1.6, 20.0)
        scored.append((t["ga"], t["kill"], max(u), arr))
    scored.sort(key=lambda x: (x[0], x[1], x[2]))
    best = scored[-0 if False else 0][3]

    fn = S.make_code_fn(best)
    ok, fails = S.check_avoid_true(fn, a.sep)
    spread_bad = 0
    for _desc, prs in S.CLASSES:
        d = sorted(fn(*S.VARS[ia]) - fn(*S.VARS[ib]) for ia, ib in prs)
        gaps = [d[i + 1] - d[i] for i in range(len(d) - 1)]
        if gaps and min(gaps) < a.sep:
            spread_bad += 1

    print(f"\n最优解：max(tx) = {max(best)} ns")
    print(f"  避真峰检查：{'通过' if ok else f'失败（{len(fails)} 处）'}")
    print(f"  散开检查  ：{'通过' if spread_bad == 0 else f'失败（{spread_bad} 组）'}")
    print("  实测（1~600m，步长 2m）：")
    for ratio in (1.6, 2.5):
        t = S.evaluate(fn, ratio, step=2.0)
        print(f"    ratio={ratio}: 鬼影 {t['gb']} -> 残留 {t['ga']} "
              f"({t['ga']/max(t['gb'],1):.3%})   误杀 {t['kill']}/{t['tb']}")

    S.dump_table(best, out_path, a.sep, a.budget,
                 note=f"v22: 散开+避真峰；search_v22.py 搜得，用时 {time.time()-t0:.0f}s")
    print("\n  tx_trig_dly [ns]")
    for l in S.LASER_IDS:
        print(f"    L{l:<3d} " + "  ".join(
            f"K{k:<2d}={int(best[S.VIDX[(l,k)]]):<3d}" for k in S.KICKS_OF[l]))


if __name__ == "__main__":
    mp.freeze_support()
    main()
