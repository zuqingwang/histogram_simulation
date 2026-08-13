# -*- coding: utf-8 -*-
"""PoD_esti v30 —— 多进程 PoD 临界点扫描（统一 bg 网格）。

用法：
    python build_pod_core_v30.py
    python run_pod_v30_noise_scan.py --workers 20   # 须先有噪声缓存
    python run_pod_v30_pod_scan.py --workers 20
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
import pod_esti_v30_core as core


def _job(args):
    n_shots, k, seed = args
    return core.solve_pod_noise(n_shots, k, seed)


def _progress_msg(value):
    """★ v30：critical 里只有 POD_FARS 那几个 tag，不能再写死 100ppm
    （写死会让每一档都误报「无有效交点」，而其实结果是好的）。"""
    crit = value.get("critical", {})
    if value.get("invalid"):
        return f"跳过：{value['invalid']}"
    tag = next((core.FAR_TAG[f] for f in core.POD_FARS if f in core.FAR_TAG
                and core.FAR_TAG[f] in crit), None)
    p90 = crit.get(tag, {}).get("0.90") if tag else None
    if not p90:
        return "无有效交点"
    return (f"{core.FAR_TAG_TO_LABEL[tag]} E90="
            f"{p90['boost']*core.E_PULSE_BASE*1e9:.3g} nJ，"
            f"验证PoD={p90['pod']:.3f}，peak均值={p90['peak_mean']:.2f}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=None)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--checkpoint-every", type=int, default=8)
    args = ap.parse_args()

    if not core._noise_cache_complete(core.NOISE_RES):
        print("缺少完整噪声缓存 / THRESH。请先：")
        print("  python run_pod_v30_noise_scan.py --workers 20")
        return 2
    if not core.THRESH:
        print("THRESH 为空，无法求 PoD")
        return 2

    n_cpu = os.cpu_count() or 1
    workers = args.workers or min(getattr(core, "N_WORKERS", 20), n_cpu)

    grid_key = np.concatenate([core.NOISE_GRID[n] for n in core.N_SHOTS_LIST])
    expected = {(ns, float(nt))
                for ns in core.N_SHOTS_LIST for nt in core.NOISE_GRID[ns]}

    res = None
    loaded_from = None
    for cand in [core.CACHE_POD, core.CACHE_POD_CKPT]:
        res = core._try_load_pod_cache(cand, grid_key)
        if res is not None:
            loaded_from = cand
            break
    if res is None:
        res = {}
        print("未找到 PoD 缓存，从零开始")
    else:
        print(f"从 {loaded_from} 载入 {len(res)} 档，断点续跑")

    pending = []
    for n_shots in core.N_SHOTS_LIST:
        for k in range(len(core.NOISE_GRID[n_shots])):
            key = (n_shots, float(core.NOISE_GRID[n_shots][k]))
            if key in res:
                continue
            pending.append((n_shots, k, 7000 + n_shots * 1_000_000 + k * 20_000))
    if args.limit:
        pending = pending[: args.limit]

    n_total = len(expected)
    if not pending:
        print(f"全部 {n_total} 档均已完成")
        core._save_pod_cache(core.CACHE_POD, res, grid_key)
        return 0

    print(f"PoD 扫描：待跑 {len(pending)} / {n_total} 档（统一 bg）；"
          f"FAR×{len(core.FAR_TAGS)}；ProcessPool workers={workers}")
    os.environ["POD_CORE_QUIET"] = "1"

    t0 = time.time()
    done_n = 0
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futs = {pool.submit(_job, spec): spec for spec in pending}
        for fut in as_completed(futs):
            key, value = fut.result()
            res[key] = value
            done_n += 1
            if done_n % args.checkpoint_every == 0:
                core._save_pod_cache(core.CACHE_POD_CKPT, res, grid_key)
                print(f"    …检查点已写入 {core.CACHE_POD_CKPT}")
            el = time.time() - t0
            eta = el / done_n * (len(pending) - done_n)
            pct = 100.0 * done_n / len(pending)
            print(f"  [{done_n}/{len(pending)} {pct:5.1f}%] N={key[0]} "
                  f"bg={key[1]:.2f}：{_progress_msg(value)}；"
                  f"已用 {el/60:.1f} min，预计剩余 {eta/60:.1f} min")

    core._save_pod_cache(core.CACHE_POD, res, grid_key)
    print(f"[PoD 扫描完成] {(time.time()-t0)/60:.1f} min → {core.CACHE_POD}")
    if len(res) >= n_total and os.path.exists(core.CACHE_POD_CKPT):
        try:
            os.remove(core.CACHE_POD_CKPT)
        except OSError:
            pass
    return 0


if __name__ == "__main__":
    sys.exit(main())
