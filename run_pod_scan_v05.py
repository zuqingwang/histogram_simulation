# -*- coding: utf-8 -*-
"""PoD_esti v05 —— 多进程 PoD 临界点全量扫描。

为什么要用进程而不是线程：
    binary_macro_stepping 是「Python 层逐细网格步循环 + 小数组 NumPy」的结构，
    绝大部分时间持有 GIL（Global Interpreter Lock，全局解释器锁），
    所以 ThreadPoolExecutor 无论开多少路，CPU 总占用都只有十几到三十几个百分点。
    改成 ProcessPoolExecutor 后每个进程有独立解释器，才能真正吃满多核。

用法（PowerShell）：
    $env:PYTHONIOENCODING="utf-8"
    python run_pod_scan_v05.py                 # 全量 208 档
    python run_pod_scan_v05.py --limit 8       # 冒烟测试：只跑 8 档
    python run_pod_scan_v05.py --workers 16    # 手动指定进程数

产物：pod_esti_v05_cache_pod.npz（增量检查点 pod_esti_v05_cache_pod.partial.npz）
跑完后在 PoD_esti_v05.ipynb 里直接命中缓存，无需重算。
"""
import argparse
import functools
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

# 长时间后台任务必须实时看到进度，否则 stdout 缓冲会让终端一片空白
print = functools.partial(print, flush=True)  # noqa: A001

import numpy as np

# 模块级 import：spawn 出来的子进程会 import 本文件，从而共享同一份内核构建流程。
import pod_esti_v05_core as core


def _job(args):
    """子进程任务：求解单个 noise 档的全部 FAR × PoD 临界点。"""
    n_shots, k, seed = args
    return core.solve_pod_noise(n_shots, k, seed)


def _progress_msg(value):
    c100 = value.get("critical", {}).get("100ppm", {})
    p90 = c100.get("0.90")
    if not p90:
        return "无有效交点"
    return (f"E90={p90['boost']*core.E_PULSE_BASE*1e9:.3g} nJ，"
            f"验证PoD={p90['pod']:.3f}，peak均值={p90['peak_mean']:.2f}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=None,
                    help="进程数，默认取 N_WORKERS 与 CPU 数的较小值")
    ap.add_argument("--limit", type=int, default=0,
                    help="只跑前 N 档（冒烟测试用），0 表示全量")
    ap.add_argument("--checkpoint-every", type=int, default=8,
                    help="每完成多少档写一次增量检查点")
    args = ap.parse_args()

    n_cpu = os.cpu_count() or 1
    workers = args.workers or min(getattr(core, "N_WORKERS", 20), n_cpu)

    grid_key = np.concatenate([core.NOISE_GRID[n] for n in core.N_SHOTS_LIST])
    expected = {(ns, float(nt))
                for ns in core.N_SHOTS_LIST for nt in core.NOISE_GRID[ns]}

    # ---- 断点续跑：主缓存 → 检查点 ----
    res = None
    loaded_from = None
    for cand in [core.CACHE_POD, core.CACHE_POD_CKPT]:
        res = core._try_load_pod_cache(cand, grid_key)
        if res is not None:
            loaded_from = cand
            break
    if res is None:
        res = {}
        print("未找到可用的 PoD 缓存，从零开始扫描")
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
        pending = pending[:args.limit]

    n_total = len(expected)
    if not pending:
        print(f"全部 {n_total} 档均已完成，无需计算")
        core._save_pod_cache(core.CACHE_POD, res, grid_key)
        return

    print(f"PoD 临界点扫描：待跑 {len(pending)} / 共 {n_total} 档；"
          f"每档求 {len(core.FAR_TAGS)} 档 FAR × PoD50/90")
    print(f"并行方式：ProcessPoolExecutor，进程数 {workers}（本机 CPU {n_cpu}）")

    # 子进程 import 内核时静音，否则 20 份启动日志会淹没进度
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
            if done_n == 1 or done_n % 5 == 0 or done_n == len(pending):
                el = time.time() - t0
                rate = el / done_n
                eta = rate * (len(pending) - done_n)
                print(f"  [{done_n}/{len(pending)}] N_shots={key[0]} "
                      f"noise={key[1]:.2f}：{_progress_msg(value)}；"
                      f"累计 {el/60:.1f} min，预计剩余 {eta/60:.1f} min")

    core._save_pod_cache(core.CACHE_POD, res, grid_key)
    print(f"扫描完成，总用时 {(time.time()-t0)/60:.1f} min；已写入 {core.CACHE_POD}")
    if len(res) >= n_total and os.path.exists(core.CACHE_POD_CKPT):
        try:
            os.remove(core.CACHE_POD_CKPT)
        except OSError:
            pass


if __name__ == "__main__":
    sys.exit(main())
