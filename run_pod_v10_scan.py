# -*- coding: utf-8 -*-
"""PoD_esti v10 —— 基于 hist_i / hist_add 的全量重算扫描（不复用 v05 缓存）。

架构
  每次任务仿 N_SHOTS_MAX=4 发 → hist_i；再对 N∈{1,2,4} 做前缀和得 hist_add。
  noise = 单次 hist_i 统计窗均值；bg / peak 在 hist_add 上统计。

两个子命令
  noise   纯噪声扫 NOISE_GRID_AMB → pod_esti_v10_cache_noise.npz
  signal  固定 boost、扫 noise → pod_esti_v10_cache_signal.npz

用法（PowerShell）
  $env:PYTHONIOENCODING="utf-8"
  python run_pod_v10_scan.py noise --workers 20
  python run_pod_v10_scan.py noise --limit 4 --n-mc 20000   # 冒烟
  python run_pod_v10_scan.py signal --workers 20
  python run_pod_v10_scan.py all --workers 20
"""
from __future__ import annotations

import argparse
import functools
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

print = functools.partial(print, flush=True)  # noqa: A001

os.environ.setdefault("POD_CORE_QUIET", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import numpy as np
import pod_esti_v10_core as core

# ---- 本版扫描参数（物理量沿用 core，不改）----
N_MC_NOISE_DEFAULT = 200_000
N_MC_SIG_DEFAULT = 8_000
MC_CHUNK_NOISE = 2_500          # 20 进程时控制内存；每块内仿 4 发
MC_CHUNK_SIG = 1_000
SEED_BASE_NOISE = 710_000
SEED_BASE_SIG = 720_000

# 固定信号档（与 peak_vs_noise 同量级，避免顶满硬上限）
BOOST_LIST = [0.0, 0.004, 0.008, 0.016, 0.032]


def _atomic_savez(path, **kwargs):
    path = str(path)
    tmp = path + ".tmp.npz"
    np.savez_compressed(tmp, **kwargs)
    os.replace(tmp, path)


# ============================ 纯噪声 ============================
def _noise_job(args):
    """单档：目标单次 noise → 仿 N_SHOTS_MAX 发 → 派生 N=1/2/4 统计。"""
    k, noise_t, n_mc, chunk, seed0 = args
    # noise_target 是单次 hist_i（27 SPAD）的底噪
    r_det = float(core.r_det_for_noise(noise_t, core.N_PIX_MACRO))
    e_lam = float(core.e_lambda_for_r_det(r_det))
    p_eq = float(core.p_bin_equilibrium(r_det)[0])
    inv_tab = core.build_inv_table(r_det)

    acc = {n: {
        "noise_sum": 0.0, "noise_sumsq": 0.0,
        "bg_sum": 0.0, "bg_sumsq": 0.0,
        "peak_cnt": np.zeros(core.N_PIX_MACRO * n + 2, dtype=np.int64),
        "n": 0,
    } for n in core.N_SHOTS_LIST}

    done = 0
    part = 0
    while done < n_mc:
        m = min(chunk, n_mc - done)
        rng = np.random.default_rng(seed0 + 104_729 * part)
        hist_i = core.noise_hists_per_shot(m, core.N_SHOTS_MAX, r_det, rng, inv_tab=inv_tab)
        st = core.stats_from_hist_i(hist_i)
        for n in core.N_SHOTS_LIST:
            a = acc[n]
            b = st[n]
            a["noise_sum"] += b["noise_sum"]
            a["noise_sumsq"] += b["noise_sumsq"]
            a["bg_sum"] += b["bg_sum"]
            a["bg_sumsq"] += b["bg_sumsq"]
            a["peak_cnt"] += b["peak_cnt"]
            a["n"] += b["n"]
        done += m
        part += 1

    out = {
        "k": k, "noise_target": float(noise_t),
        "r_det": r_det, "e_lambda": e_lam, "p_eq": p_eq,
    }
    for n in core.N_SHOTS_LIST:
        a = acc[n]
        nn = max(a["n"], 1)
        out[f"n{n}"] = {
            "n": a["n"],
            "noise_mc": a["noise_sum"] / nn,
            "noise_std": float(np.sqrt(max(a["noise_sumsq"] / nn - (a["noise_sum"] / nn) ** 2, 0.0))),
            "bg_mc": a["bg_sum"] / nn,
            "bg_std": float(np.sqrt(max(a["bg_sumsq"] / nn - (a["bg_sum"] / nn) ** 2, 0.0))),
            "peak_cnt": a["peak_cnt"],
        }
    return out


def _noise_cache_key(n_mc, grid):
    return (np.asarray(grid, float), np.asarray(core.N_SHOTS_LIST, int), int(n_mc))


def _save_noise(path, rows, n_mc, grid):
    grid = np.asarray(grid, float)
    gkey, nslist, n_mc = _noise_cache_key(n_mc, grid)
    ng = len(grid)
    r_det = np.full(ng, np.nan)
    e_lam = np.full(ng, np.nan)
    p_eq = np.full(ng, np.nan)
    done = np.zeros(ng, dtype=bool)
    for k, r in rows.items():
        r_det[k] = r["r_det"]
        e_lam[k] = r["e_lambda"]
        p_eq[k] = r["p_eq"]
        done[k] = True
    payload = {
        "grid_key": gkey, "n_shots_list": nslist, "n_mc": n_mc,
        "noise_target": grid,
        "r_det": r_det, "e_lambda": e_lam, "p_eq": p_eq, "done": done,
    }
    for n in core.N_SHOTS_LIST:
        n_tr = core.N_PIX_MACRO * n
        payload[f"noise_mc_{n}"] = np.full(ng, np.nan)
        payload[f"noise_std_{n}"] = np.full(ng, np.nan)
        payload[f"bg_mc_{n}"] = np.full(ng, np.nan)
        payload[f"bg_std_{n}"] = np.full(ng, np.nan)
        payload[f"peak_cnt_{n}"] = np.zeros((ng, n_tr + 2), dtype=np.int64)
        for k, r in rows.items():
            d = r[f"n{n}"]
            payload[f"noise_mc_{n}"][k] = d["noise_mc"]
            payload[f"noise_std_{n}"][k] = d["noise_std"]
            payload[f"bg_mc_{n}"][k] = d["bg_mc"]
            payload[f"bg_std_{n}"][k] = d["bg_std"]
            payload[f"peak_cnt_{n}"][k] = d["peak_cnt"]
    _atomic_savez(path, **payload)


def _load_noise(path, n_mc, grid):
    if not os.path.exists(path):
        return {}
    z = np.load(path, allow_pickle=True)
    gkey, nslist, n_mc_k = _noise_cache_key(n_mc, grid)
    if (int(z["n_mc"]) != n_mc_k
            or not np.array_equal(z["n_shots_list"], nslist)
            or z["grid_key"].shape != gkey.shape
            or not np.allclose(z["grid_key"], gkey)):
        print(f"[noise] 缓存键不匹配，忽略 {path}")
        return {}
    rows = {}
    done = z["done"]
    for k in range(len(grid)):
        if not done[k]:
            continue
        rec = {
            "k": k, "noise_target": float(z["noise_target"][k]),
            "r_det": float(z["r_det"][k]),
            "e_lambda": float(z["e_lambda"][k]),
            "p_eq": float(z["p_eq"][k]),
        }
        for n in core.N_SHOTS_LIST:
            rec[f"n{n}"] = {
                "n": n_mc,
                "noise_mc": float(z[f"noise_mc_{n}"][k]),
                "noise_std": float(z[f"noise_std_{n}"][k]),
                "bg_mc": float(z[f"bg_mc_{n}"][k]),
                "bg_std": float(z[f"bg_std_{n}"][k]),
                "peak_cnt": np.asarray(z[f"peak_cnt_{n}"][k], dtype=np.int64),
            }
        rows[k] = rec
    print(f"[noise] 从 {path} 载入 {len(rows)}/{len(grid)} 档")
    return rows


def run_noise(n_mc, workers, limit=None, chunk=None):
    grid = np.asarray(core.NOISE_GRID_AMB, float)
    if limit is not None:
        grid = grid[: int(limit)]
    chunk = int(chunk or MC_CHUNK_NOISE)
    rows = _load_noise(core.CACHE_NOISE, n_mc, grid)
    if not rows:
        rows = _load_noise(core.CACHE_NOISE_CKPT, n_mc, grid)

    todo = [k for k in range(len(grid)) if k not in rows]
    print(f"[noise] 网格 {len(grid)} 档，已完成 {len(rows)}，待算 {len(todo)}；"
          f"N_MC={n_mc:,}，workers={workers}，chunk={chunk}")
    if not todo:
        _save_noise(core.CACHE_NOISE, rows, n_mc, grid)
        print("[noise] 已全部完成")
        return rows

    jobs = [(k, float(grid[k]), n_mc, chunk, SEED_BASE_NOISE + 1009 * k) for k in todo]
    t0 = time.time()
    done_n = 0
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futs = {pool.submit(_noise_job, j): j[0] for j in jobs}
        for fut in as_completed(futs):
            r = fut.result()
            rows[r["k"]] = r
            done_n += 1
            n1 = r["n1"]
            print(f"  [{done_n}/{len(todo)}] k={r['k']} noise_t={r['noise_target']:.2f} "
                  f"→ noise={n1['noise_mc']:.3f} bg(N=1/2/4)="
                  f"{r['n1']['bg_mc']:.3f}/{r['n2']['bg_mc']:.3f}/{r['n4']['bg_mc']:.3f} "
                  f"peak_mean(N=4)={core.peak_stats_from_cnt(r['n4']['peak_cnt'])['mean']:.2f} "
                  f"({time.time()-t0:.0f}s)")
            if done_n % 4 == 0 or done_n == len(todo):
                _save_noise(core.CACHE_NOISE_CKPT, rows, n_mc, grid)
    _save_noise(core.CACHE_NOISE, rows, n_mc, grid)
    if os.path.exists(core.CACHE_NOISE_CKPT):
        try:
            os.remove(core.CACHE_NOISE_CKPT)
        except OSError:
            pass
    print(f"[noise] 完成，写入 {core.CACHE_NOISE}，用时 {time.time()-t0:.0f}s")
    return rows


# ============================ 固定信号 × 扫 noise ============================
def _signal_job(args):
    """单档 noise：对所有 BOOST_LIST 各仿一次 N_SHOTS_MAX，派生各 N 的 peak_cnt。"""
    k, noise_t, r_det, n_mc, chunk, seed0 = args
    r_amb = float(r_det / core.PDE)
    n_boost = len(BOOST_LIST)
    # cnt[n][i_boost] = bincount
    cnt = {n: [np.zeros(core.N_PIX_MACRO * n + 2, dtype=np.int64)
               for _ in range(n_boost)] for n in core.N_SHOTS_LIST}
    bg_sum = {n: np.zeros(n_boost) for n in core.N_SHOTS_LIST}
    bg_n = {n: np.zeros(n_boost, dtype=np.int64) for n in core.N_SHOTS_LIST}

    for ib, boost in enumerate(BOOST_LIST):
        done = 0
        part = 0
        while done < n_mc:
            m = min(chunk, n_mc - done)
            rng = np.random.default_rng(seed0 + 10_007 * ib + 104_729 * part)
            hist_i = core.binary_macro_stepping_per_shot(
                m, core.F_VALS, core.N_SHOTS_MAX,
                core.R_SIG_UNIT_POD, core.TF_POD, r_amb, core.CENTERS_SIG,
                rng, boost=float(boost),
            )
            # 信号窗 peak；bg 用信号窗均值作对照（窗短，仅作记录）
            for n in core.N_SHOTS_LIST:
                hadd = core.hist_add_from_prefix(hist_i, n)
                pk = hadd.max(axis=1)
                bg = hadd.mean(axis=1)
                cnt[n][ib] += np.bincount(pk, minlength=cnt[n][ib].size)
                bg_sum[n][ib] += bg.sum()
                bg_n[n][ib] += m
            done += m
            part += 1

    out = {"k": k, "noise_target": float(noise_t), "r_det": float(r_det),
           "r_amb": r_amb}
    for n in core.N_SHOTS_LIST:
        out[f"n{n}"] = {
            "peak_cnt": np.stack(cnt[n], axis=0),  # (n_boost, n_tr+2)
            "bg_mc": bg_sum[n] / np.maximum(bg_n[n], 1),
        }
    return out


def _save_signal(path, rows, n_mc, grid, boosts):
    grid = np.asarray(grid, float)
    boosts = np.asarray(boosts, float)
    payload = {
        "grid_key": grid, "boosts": boosts, "n_mc": int(n_mc),
        "n_shots_list": np.asarray(core.N_SHOTS_LIST, int),
        "noise_target": grid,
        "done": np.array([k in rows for k in range(len(grid))], dtype=bool),
        "r_det": np.full(len(grid), np.nan),
    }
    for n in core.N_SHOTS_LIST:
        n_tr = core.N_PIX_MACRO * n
        payload[f"peak_cnt_{n}"] = np.zeros(
            (len(boosts), len(grid), n_tr + 2), dtype=np.int64)
        payload[f"bg_mc_{n}"] = np.full((len(boosts), len(grid)), np.nan)
    for k, r in rows.items():
        payload["r_det"][k] = r["r_det"]
        for n in core.N_SHOTS_LIST:
            payload[f"peak_cnt_{n}"][:, k, :] = r[f"n{n}"]["peak_cnt"]
            payload[f"bg_mc_{n}"][:, k] = r[f"n{n}"]["bg_mc"]
    _atomic_savez(path, **payload)


def _load_signal(path, n_mc, grid, boosts):
    if not os.path.exists(path):
        return {}
    z = np.load(path, allow_pickle=True)
    grid = np.asarray(grid, float)
    boosts = np.asarray(boosts, float)
    if (int(z["n_mc"]) != int(n_mc)
            or not np.allclose(z["grid_key"], grid)
            or not np.allclose(z["boosts"], boosts)
            or not np.array_equal(z["n_shots_list"], np.asarray(core.N_SHOTS_LIST))):
        print(f"[signal] 缓存键不匹配，忽略 {path}")
        return {}
    rows = {}
    for k in range(len(grid)):
        if not z["done"][k]:
            continue
        rec = {"k": k, "noise_target": float(z["noise_target"][k]),
               "r_det": float(z["r_det"][k]),
               "r_amb": float(z["r_det"][k] / core.PDE)}
        for n in core.N_SHOTS_LIST:
            rec[f"n{n}"] = {
                "peak_cnt": np.asarray(z[f"peak_cnt_{n}"][:, k, :], dtype=np.int64),
                "bg_mc": np.asarray(z[f"bg_mc_{n}"][:, k], dtype=float),
            }
        rows[k] = rec
    print(f"[signal] 从 {path} 载入 {len(rows)}/{len(grid)} 档")
    return rows


def run_signal(n_mc, workers, limit=None, chunk=None, noise_rows=None):
    grid = np.asarray(core.NOISE_GRID_AMB, float)
    if limit is not None:
        grid = grid[: int(limit)]
    chunk = int(chunk or MC_CHUNK_SIG)
    # r_det：优先用 noise 扫描实测反解；否则用理论 r_det_for_noise
    r_dets = np.array([core.r_det_for_noise(float(nt), core.N_PIX_MACRO) for nt in grid])
    if noise_rows:
        for k, r in noise_rows.items():
            if k < len(grid):
                r_dets[k] = r["r_det"]

    rows = _load_signal(core.CACHE_SIG, n_mc, grid, BOOST_LIST)
    if not rows:
        rows = _load_signal(core.CACHE_SIG_CKPT, n_mc, grid, BOOST_LIST)

    todo = [k for k in range(len(grid)) if k not in rows]
    print(f"[signal] 网格 {len(grid)} 档 × boost{BOOST_LIST}，已完成 {len(rows)}，待算 {len(todo)}；"
          f"N_MC={n_mc:,}，workers={workers}")
    if not todo:
        _save_signal(core.CACHE_SIG, rows, n_mc, grid, BOOST_LIST)
        print("[signal] 已全部完成")
        return rows

    jobs = [(k, float(grid[k]), float(r_dets[k]), n_mc, chunk,
             SEED_BASE_SIG + 1009 * k) for k in todo]
    t0 = time.time()
    done_n = 0
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futs = {pool.submit(_signal_job, j): j[0] for j in jobs}
        for fut in as_completed(futs):
            r = fut.result()
            rows[r["k"]] = r
            done_n += 1
            st = core.peak_stats_from_cnt(r["n4"]["peak_cnt"][-1])  # 最强信号
            print(f"  [{done_n}/{len(todo)}] k={r['k']} noise_t={r['noise_target']:.2f} "
                  f"peak_mean(N=4,boost={BOOST_LIST[-1]})={st['mean']:.2f} "
                  f"({time.time()-t0:.0f}s)")
            if done_n % 4 == 0 or done_n == len(todo):
                _save_signal(core.CACHE_SIG_CKPT, rows, n_mc, grid, BOOST_LIST)
    _save_signal(core.CACHE_SIG, rows, n_mc, grid, BOOST_LIST)
    if os.path.exists(core.CACHE_SIG_CKPT):
        try:
            os.remove(core.CACHE_SIG_CKPT)
        except OSError:
            pass
    print(f"[signal] 完成，写入 {core.CACHE_SIG}，用时 {time.time()-t0:.0f}s")
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["noise", "signal", "all"])
    ap.add_argument("--workers", type=int, default=20)
    ap.add_argument("--n-mc", type=int, default=None)
    ap.add_argument("--limit", type=int, default=None, help="只跑前若干 noise 档（冒烟）")
    ap.add_argument("--chunk", type=int, default=None)
    args = ap.parse_args()

    if args.mode in ("noise", "all"):
        n_mc = args.n_mc or N_MC_NOISE_DEFAULT
        rows = run_noise(n_mc, args.workers, limit=args.limit, chunk=args.chunk)
    else:
        rows = None

    if args.mode in ("signal", "all"):
        n_mc = args.n_mc or N_MC_SIG_DEFAULT
        # all 模式：signal 用自己的默认 n_mc，不被 noise 的 n_mc 带偏
        if args.mode == "all" and args.n_mc is None:
            n_mc = N_MC_SIG_DEFAULT
        if rows is None and os.path.exists(core.CACHE_NOISE):
            # 尝试载入以取 r_det
            grid = np.asarray(core.NOISE_GRID_AMB, float)
            if args.limit is not None:
                grid = grid[: args.limit]
            rows = _load_noise(core.CACHE_NOISE, N_MC_NOISE_DEFAULT
                               if args.n_mc is None else args.n_mc, grid)
        run_signal(n_mc, args.workers, limit=args.limit, chunk=args.chunk,
                   noise_rows=rows)


if __name__ == "__main__":
    main()
