# -*- coding: utf-8 -*-
"""PoD_esti v30 —— 多进程纯噪声扫描（统一 bg 步长 0.25）。

每档任务键：(N, bg_k)；仿真 noise_amb = bg / N。
不再用 AMB×N 前缀和（那会导致 N=2/4 的 bg 步长变成 0.5/1.0）。

★ v30 相对 v20 的唯一增量：额外累计 hist_std（单条 hist_add 在统计窗 152 个 bin 上的
样本 std），让 notebook 模块 10 能直接复用这批 1e6 MC，不必再单跑一次小扫描。

用法：
    python build_pod_core_v30.py
    python run_pod_v30_noise_scan.py --workers 20
    python run_pod_v30_noise_scan.py --limit 4 --n-mc 20000   # 冒烟
"""
from __future__ import annotations

import argparse
import functools
import os
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

print = functools.partial(print, flush=True)  # noqa: A001

os.environ.setdefault("POD_CORE_QUIET", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import numpy as np

import pod_esti_v30_core as core


def _job(args):
    """单个 (N, bg) 档。"""
    n_shots, k, bg_target, n_mc, chunk, seed0 = args
    nt_amb = float(bg_target) / float(n_shots)
    r_det = float(core.r_det_for_noise(nt_amb, core.N_PIX_MACRO))
    e_lam = float(core.e_lambda_for_r_det(r_det))
    p_eq = float(core.p_bin_equilibrium(r_det)[0])
    inv_tab = core.build_inv_table(r_det)

    acc = dict(
        noise_sum=0.0, noise_sumsq=0.0,
        bg_sum=0.0, bg_sumsq=0.0,
        hist_std_sum=0.0,          # ★ v30
        peak_cnt=np.zeros(core.N_PIX_MACRO * n_shots + 2, dtype=np.int64),
        nn=0,
    )
    done, part = 0, 0
    while done < n_mc:
        m = min(chunk, n_mc - done)
        rng = np.random.default_rng(seed0 + 104_729 * part)
        hist_i = core.noise_hists_per_shot(
            m, n_shots, r_det, rng, inv_tab=inv_tab)
        st = core.stats_from_hist_i(hist_i, n_shots_list=[n_shots])
        b = st[n_shots]
        acc["noise_sum"] += b["noise_sum"]
        acc["noise_sumsq"] += b["noise_sumsq"]
        acc["bg_sum"] += b["bg_sum"]
        acc["bg_sumsq"] += b["bg_sumsq"]
        acc["hist_std_sum"] += b["hist_std_sum"]
        acc["peak_cnt"] += b["peak_cnt"]
        acc["nn"] += b["n"]
        done += m
        part += 1

    nn = max(acc["nn"], 1)
    return {
        "n_shots": int(n_shots), "k": int(k),
        "bg_target": float(bg_target), "noise_amb": float(nt_amb),
        "r_det": r_det, "e_lambda": e_lam, "p_eq": p_eq,
        "noise_amb_mc": acc["noise_sum"] / nn,
        "noise_amb_std": float(np.sqrt(max(
            acc["noise_sumsq"] / nn - (acc["noise_sum"] / nn) ** 2, 0.0))),
        "bg_mc": acc["bg_sum"] / nn,
        "bg_std": float(np.sqrt(max(
            acc["bg_sumsq"] / nn - (acc["bg_sum"] / nn) ** 2, 0.0))),
        "hist_std": acc["hist_std_sum"] / nn,      # ★ v30
        "peak_cnt": acc["peak_cnt"],
        "nn": acc["nn"],
    }


def _empty_res(grid, n_mc):
    ng = len(grid)
    res = {}
    for n in core.N_SHOTS_LIST:
        n_tr = core.N_PIX_MACRO * n
        res[n] = {
            "n_shots": n, "n_tr": n_tr,
            "noise_target": grid.copy(),
            "noise_amb_target": np.round(grid / n, 6),
            "r_det": np.zeros(ng), "e_lambda": np.zeros(ng), "p_eq": np.zeros(ng),
            "noise_mc": np.zeros(ng), "noise_std": np.zeros(ng),
            "noise_amb_mc": np.zeros(ng), "noise_amb_std": np.zeros(ng),
            "hist_std": np.zeros(ng),                      # ★ v30
            "peak_cnt": np.zeros((ng, n_tr + 2), dtype=np.int64),
            "done": np.zeros(ng, dtype=bool),
        }
    return res


def _apply(res, r):
    n, k = r["n_shots"], r["k"]
    R = res[n]
    R["r_det"][k] = r["r_det"]
    R["e_lambda"][k] = r["e_lambda"]
    R["p_eq"][k] = r["p_eq"]
    R["noise_amb_mc"][k] = r["noise_amb_mc"]
    R["noise_amb_std"][k] = r["noise_amb_std"]
    R["noise_mc"][k] = r["bg_mc"]
    R["noise_std"][k] = r["bg_std"]
    R["hist_std"][k] = r["hist_std"]
    R["peak_cnt"][k] = r["peak_cnt"]
    R["done"][k] = True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=20)
    ap.add_argument("--n-mc", type=int, default=None)
    ap.add_argument("--chunk", type=int, default=None)
    ap.add_argument("--limit", type=int, default=0,
                    help="每 N 只跑前 limit 个 bg 档（冒烟）")
    ap.add_argument("--checkpoint-every", type=int, default=8)
    args = ap.parse_args()

    n_mc = int(args.n_mc or core.N_MC_NOISE)
    chunk = int(args.chunk or core.MC_CHUNK)
    grid = np.asarray(core.BG_GRID, float)
    if args.limit:
        grid = grid[: int(args.limit)]
    grid_key = np.asarray(core.BG_GRID, float)  # 缓存键始终用完整 BG_GRID

    res = core._try_load_noise_cache(core.CACHE_NOISE, grid_key)
    loaded = core.CACHE_NOISE if res is not None else None
    if res is None:
        res = core._try_load_noise_cache(core.CACHE_NOISE_CKPT, grid_key)
        loaded = core.CACHE_NOISE_CKPT if res is not None else None
    if res is None:
        res = _empty_res(np.asarray(core.BG_GRID, float), n_mc)
        print("未找到缓存，从零开始")
    else:
        for n in core.N_SHOTS_LIST:
            if n not in res:
                res = _empty_res(np.asarray(core.BG_GRID, float), n_mc)
                loaded = None
                print("缓存结构不完整，从零开始")
                break
            R = res[n]
            if "done" not in R:
                R["done"] = np.array([int(c.sum()) > 0 for c in R["peak_cnt"]], bool)
            if len(R["noise_target"]) != len(core.BG_GRID):
                res = _empty_res(np.asarray(core.BG_GRID, float), n_mc)
                loaded = None
                print("网格长度不匹配，从零开始")
                break
        if loaded:
            n_done = sum(int(np.sum(res[n]["done"])) for n in core.N_SHOTS_LIST)
            n_tot = len(core.BG_GRID) * len(core.N_SHOTS_LIST)
            print(f"从 {loaded} 载入，已完成 {n_done}/{n_tot} 档")

    full_grid = np.asarray(core.BG_GRID, float)
    todo = []
    for n in core.N_SHOTS_LIST:
        for k in range(len(full_grid)):
            if args.limit and k >= int(args.limit):
                continue
            if not bool(res[n]["done"][k]):
                todo.append((n, k, float(full_grid[k]), n_mc, chunk,
                             2000 + 1009 * n + 17 * k))

    print(f"统一 BG {len(full_grid)} 档 × N={list(core.N_SHOTS_LIST)} × {n_mc:,} 条；"
          f"待算 {len(todo)}；workers={args.workers}（进程）")
    if not todo:
        core._save_noise_cache(core.CACHE_NOISE, res, grid_key)
        print("已全部完成")
        return

    t0 = time.time()
    done_n = 0
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futs = {pool.submit(_job, j): j for j in todo}
        for fut in as_completed(futs):
            r = fut.result()
            _apply(res, r)
            done_n += 1
            pk = core.peak_stats_from_cnt(r["peak_cnt"])
            el = time.time() - t0
            eta = el / done_n * (len(todo) - done_n)
            pct = 100.0 * done_n / len(todo)
            print(f"  [{done_n}/{len(todo)} {pct:5.1f}%] N={r['n_shots']} "
                  f"bg={r['bg_target']:.2f}（amb={r['noise_amb']:.3f}）→ "
                  f"bg_mc={r['bg_mc']:.3f} peakμ={pk['mean']:.2f}  "
                  f"已用 {el/60:.1f} min，预计剩余 {eta/60:.1f} min")
            if done_n % args.checkpoint_every == 0 or done_n == len(todo):
                core._save_noise_cache(core.CACHE_NOISE_CKPT, res, grid_key)
                print(f"    …检查点已写入 {core.CACHE_NOISE_CKPT}")

    core._save_noise_cache(core.CACHE_NOISE, res, grid_key)
    if os.path.exists(core.CACHE_NOISE_CKPT):
        try:
            os.remove(core.CACHE_NOISE_CKPT)
        except OSError:
            pass
    print(f"[噪声扫描完成] → {core.CACHE_NOISE}，总用时 {(time.time()-t0)/60:.1f} min")


if __name__ == "__main__":
    main()
