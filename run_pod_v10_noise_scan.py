# -*- coding: utf-8 -*-
"""PoD_esti v10 —— 多进程纯噪声扫描（hist_i 前缀和）。

为什么 notebook 吃不满 CPU：
    noise_hists_per_shot / noise_macro_hist_fast 大量时间在 Python 循环里，
    持有 GIL（Global Interpreter Lock，全局解释器锁）。
    ThreadPoolExecutor 再多线程也只有一核在跑 Python。
    与 v05 的 run_pod_scan_v05.py 一样，必须用 ProcessPoolExecutor 才能吃满多核。

用法（PowerShell）：先停掉 notebook 里正在跑的模块 5，再执行：
    $env:PYTHONIOENCODING="utf-8"
    python build_pod_core_v10.py                 # 若改过 notebook 计算 cell
    python run_pod_v10_noise_scan.py --workers 20
    python run_pod_v10_noise_scan.py --limit 4 --n-mc 20000   # 冒烟

产物：pod_esti_v10_cache_noise.npz
然后回 notebook 重跑模块 5 cell —— 应秒级命中缓存。
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

import pod_esti_v10_core as core


def _job(args):
    """单个 AMB noise 档：仿 N_SHOTS_MAX 发 hist_i，前缀和得 N=1/2/4 统计。"""
    k, nt_amb, n_mc, chunk, seed0 = args
    r_det = float(core.r_det_for_noise(float(nt_amb), core.N_PIX_MACRO))
    e_lam = float(core.e_lambda_for_r_det(r_det))
    p_eq = float(core.p_bin_equilibrium(r_det)[0])
    inv_tab = core.build_inv_table(r_det)

    acc = {
        n: dict(
            noise_sum=0.0, noise_sumsq=0.0,
            bg_sum=0.0, bg_sumsq=0.0,
            peak_cnt=np.zeros(core.N_PIX_MACRO * n + 2, dtype=np.int64),
            nn=0,
        )
        for n in core.N_SHOTS_LIST
    }
    done, part = 0, 0
    while done < n_mc:
        m = min(chunk, n_mc - done)
        rng = np.random.default_rng(seed0 + 104_729 * part)
        hist_i = core.noise_hists_per_shot(
            m, core.N_SHOTS_MAX, r_det, rng, inv_tab=inv_tab)
        st = core.stats_from_hist_i(hist_i)
        for n in core.N_SHOTS_LIST:
            a, b = acc[n], st[n]
            a["noise_sum"] += b["noise_sum"]
            a["noise_sumsq"] += b["noise_sumsq"]
            a["bg_sum"] += b["bg_sum"]
            a["bg_sumsq"] += b["bg_sumsq"]
            a["peak_cnt"] += b["peak_cnt"]
            a["nn"] += b["n"]
        done += m
        part += 1

    out = {
        "k": k, "noise_amb": float(nt_amb),
        "r_det": r_det, "e_lambda": e_lam, "p_eq": p_eq,
        "by_N": {},
    }
    for n in core.N_SHOTS_LIST:
        a = acc[n]
        nn = max(a["nn"], 1)
        out["by_N"][n] = {
            "noise_amb_mc": a["noise_sum"] / nn,
            "noise_amb_std": float(np.sqrt(max(
                a["noise_sumsq"] / nn - (a["noise_sum"] / nn) ** 2, 0.0))),
            "bg_mc": a["bg_sum"] / nn,
            "bg_std": float(np.sqrt(max(
                a["bg_sumsq"] / nn - (a["bg_sum"] / nn) ** 2, 0.0))),
            "peak_cnt": a["peak_cnt"],
            "nn": a["nn"],
        }
    return out


def _empty_res(grid, n_mc):
    ng = len(grid)
    res = {}
    for n in core.N_SHOTS_LIST:
        n_tr = core.N_PIX_MACRO * n
        res[n] = {
            "n_shots": n, "n_tr": n_tr,
            "noise_target": np.round(grid * n, 4),
            "noise_amb_target": grid.copy(),
            "r_det": np.zeros(ng), "e_lambda": np.zeros(ng), "p_eq": np.zeros(ng),
            "noise_mc": np.zeros(ng), "noise_std": np.zeros(ng),
            "noise_amb_mc": np.zeros(ng), "noise_amb_std": np.zeros(ng),
            "peak_cnt": np.zeros((ng, n_tr + 2), dtype=np.int64),
            "done": np.zeros(ng, dtype=bool),
        }
    return res


def _apply(res, r):
    k = r["k"]
    for n in core.N_SHOTS_LIST:
        R = res[n]
        d = r["by_N"][n]
        R["r_det"][k] = r["r_det"]
        R["e_lambda"][k] = r["e_lambda"]
        R["p_eq"][k] = r["p_eq"]
        R["noise_amb_mc"][k] = d["noise_amb_mc"]
        R["noise_amb_std"][k] = d["noise_amb_std"]
        R["noise_mc"][k] = d["bg_mc"]       # 兼容 notebook：noise_mc = bg
        R["noise_std"][k] = d["bg_std"]
        R["peak_cnt"][k] = d["peak_cnt"]
        R["done"][k] = True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=20)
    ap.add_argument("--n-mc", type=int, default=None,
                    help=f"默认用 core.N_MC_NOISE={core.N_MC_NOISE}")
    ap.add_argument("--chunk", type=int, default=None,
                    help=f"默认 MC_CHUNK={core.MC_CHUNK}")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--checkpoint-every", type=int, default=4)
    args = ap.parse_args()

    n_mc = int(args.n_mc or core.N_MC_NOISE)
    chunk = int(args.chunk or core.MC_CHUNK)
    grid = np.asarray(core.NOISE_GRID_AMB, float)
    if args.limit:
        grid = grid[: int(args.limit)]
    grid_key = grid.copy()

    res = core._try_load_noise_cache(core.CACHE_NOISE, grid_key)
    loaded = core.CACHE_NOISE if res is not None else None
    if res is None:
        res = core._try_load_noise_cache(core.CACHE_NOISE_CKPT, grid_key)
        loaded = core.CACHE_NOISE_CKPT if res is not None else None
    if res is None:
        res = _empty_res(grid, n_mc)
        print("未找到缓存，从零开始")
    else:
        # 兼容：补 done / 缺字段
        for n in core.N_SHOTS_LIST:
            if n not in res:
                res = _empty_res(grid, n_mc)
                loaded = None
                print("缓存结构不完整，从零开始")
                break
            R = res[n]
            if "done" not in R:
                R["done"] = np.array([int(c.sum()) > 0 for c in R["peak_cnt"]], bool)
            if len(R["noise_target"]) != len(grid):
                res = _empty_res(grid, n_mc)
                loaded = None
                print("网格长度不匹配，从零开始")
                break
        if loaded:
            print(f"从 {loaded} 载入，已完成 "
                  f"{int(np.sum(res[core.N_SHOTS_LIST[0]]['done']))}/{len(grid)} 档")

    todo = [k for k in range(len(grid))
            if not all(bool(res[n]["done"][k]) for n in core.N_SHOTS_LIST)]
    print(f"AMB {len(grid)} 档 × N={list(core.N_SHOTS_LIST)}（前缀和）× {n_mc:,} 条；"
          f"待算 {len(todo)}；workers={args.workers}（进程）")
    if not todo:
        core._save_noise_cache(core.CACHE_NOISE, res, grid_key)
        print("已全部完成")
        return

    jobs = [(k, float(grid[k]), n_mc, chunk, 2000 + 1009 * k) for k in todo]
    t0 = time.time()
    done_n = 0
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futs = {pool.submit(_job, j): j[0] for j in jobs}
        for fut in as_completed(futs):
            r = fut.result()
            _apply(res, r)
            done_n += 1
            pk = core.peak_stats_from_cnt(r["by_N"][4]["peak_cnt"])
            el = time.time() - t0
            eta = el / done_n * (len(todo) - done_n)
            pct = 100.0 * done_n / len(todo)
            print(f"  [{done_n}/{len(todo)} {pct:5.1f}%] k={r['k']} "
                  f"noise={r['noise_amb']:.2f} → bg(N=1/2/4)="
                  f"{r['by_N'][1]['bg_mc']:.3f}/"
                  f"{r['by_N'][2]['bg_mc']:.3f}/"
                  f"{r['by_N'][4]['bg_mc']:.3f}  "
                  f"peakμ(N=4)={pk['mean']:.2f}  "
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
