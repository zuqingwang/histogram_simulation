# -*- coding: utf-8 -*-
"""peak_distribution v01 —— 信号强度 ×2 时，peak 分布统计量是否也 ×2？

要回答的问题
    固定环境噪声，信号强度（boost）加倍时：
      · peak 分布的众数 / 中位数 / 均值 / p90 是否也加倍？
      · PoD50 / PoD90 临界能量处，上述统计量如何缩放？
    PoD（Probability of Detection，检测概率）= P(peak ≥ T)，T 取自 PoD_esti_v05 的 FAR 阈值。

设计
  1. 物理内核 import pod_esti_v05_core，不复制参数。
  2. boost 用对数网格，使大量 (b, 2b) 对同时落在网格上，便于算缩放比。
  3. 缓存存 peak 的 bincount（充分统计量）。
  4. ProcessPoolExecutor（GIL 限制，线程吃不满 CPU）。

用法（PowerShell）
    $env:PYTHONIOENCODING="utf-8"
    python peak_distribution_scan.py
    python peak_distribution_scan.py --limit 4 --n-mc 2000   # 冒烟
"""
import argparse
import functools
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

print = functools.partial(print, flush=True)  # noqa: A001

import numpy as np
import pod_esti_v05_core as core

# ============================ 本工作参数 ============================
# 对数网格：相邻点近似 ×√2，因此每隔 2 档正好是 ×2，便于缩放对比
_BOOST_POS = np.round(np.logspace(np.log10(0.001), np.log10(0.08), 25), 6)
BOOST_LIST = np.unique(np.concatenate([[0.0], _BOOST_POS])).tolist()

# 代表噪声档（必须落在 core.NOISE_GRID 上）
NOISE_PICK = {
    1: [0.5, 2.0, 5.0, 8.0],
    4: [1.0, 5.0, 15.0, 30.0],
}

N_MC_DEFAULT = 8000
MC_CHUNK = 2000
SEED_BASE = 77_000

CACHE = "peak_distribution_v01_cache.npz"
CACHE_CKPT = "peak_distribution_v01_cache.partial.npz"

# PoD 用的 FAR 标签（与 core.FAR_TAGS 对齐）
POD_FAR_TAGS = ["100ppm", "10ppm"]


def _noise_index(n_shots, noise_t):
    grid = core.NOISE_GRID[n_shots]
    k = int(np.argmin(np.abs(grid - noise_t)))
    if abs(float(grid[k]) - float(noise_t)) > 1e-6:
        raise ValueError(f"noise={noise_t} 不在 NOISE_GRID[{n_shots}] 上，最近={grid[k]}")
    return k


def _peak_bincount(boost, n_shots, r_amb, n_mc, seed, n_tr):
    cnt = np.zeros(n_tr + 2, dtype=np.int64)
    for s in range(0, n_mc, MC_CHUNK):
        m = min(MC_CHUNK, n_mc - s)
        pk = core._peaks_chunk(boost, n_shots, r_amb, m, seed + 104_729 * s)
        cnt += np.bincount(pk, minlength=n_tr + 2)[: n_tr + 2]
    return cnt


def _job(args):
    """单个 (n_shots, noise) 任务：扫完整 BOOST_LIST。"""
    n_shots, noise_t, n_mc, seed0 = args
    n_tr = 27 * n_shots
    k = _noise_index(n_shots, noise_t)
    R = core.NOISE_RES[n_shots]
    r_amb = float(R["r_det"][k] / core.PDE)
    noise_mc = float(R["noise_mc"][k])
    e_lam = float(R["e_lambda"][k])

    cnt = np.zeros((len(BOOST_LIST), n_tr + 2), dtype=np.int64)
    for i, boost in enumerate(BOOST_LIST):
        cnt[i] = _peak_bincount(boost, n_shots, r_amb, n_mc,
                                seed0 + 1009 * i, n_tr)

    # 该 noise 档的 FAR 阈值
    Tr = core.THRESH[n_shots]
    T_map = {tag: int(Tr["T" + tag][k]) for tag in core.FAR_TAGS}

    return (n_shots, float(noise_t)), {
        "cnt": cnt, "noise_target": float(noise_t),
        "noise_mc": noise_mc, "e_lambda": e_lam, "r_amb": r_amb,
        "T_map": T_map, "k": k,
    }


def _cache_key(n_mc):
    return (np.asarray(BOOST_LIST, float),
            {n: np.asarray(NOISE_PICK[n], float) for n in (1, 4)},
            int(n_mc))


def _save(path, res, n_mc):
    boosts, picks, n_mc = _cache_key(n_mc)
    out = {
        "boosts": boosts, "n_mc": n_mc,
        "far_tags": np.array(core.FAR_TAGS),
        "pod_far_tags": np.array(POD_FAR_TAGS),
        "e_pulse_base": float(core.E_PULSE_BASE),
    }
    for n_shots in (1, 4):
        noises = picks[n_shots]
        n_tr = 27 * n_shots
        nb, ng = len(BOOST_LIST), len(noises)
        cnt = np.zeros((nb, ng, n_tr + 2), dtype=np.int64)
        done = np.zeros(ng, dtype=bool)
        noise_mc = np.full(ng, np.nan)
        e_lam = np.full(ng, np.nan)
        T_arr = np.full((len(core.FAR_TAGS), ng), -1, dtype=int)
        for j, nt in enumerate(noises):
            r = res.get((n_shots, float(nt)))
            if r is None:
                continue
            cnt[:, j, :] = r["cnt"]
            done[j] = True
            noise_mc[j] = r["noise_mc"]
            e_lam[j] = r["e_lambda"]
            for ti, tag in enumerate(core.FAR_TAGS):
                T_arr[ti, j] = r["T_map"][tag]
        out[f"noise_{n_shots}"] = noises
        out[f"cnt_{n_shots}"] = cnt
        out[f"done_{n_shots}"] = done
        out[f"noisemc_{n_shots}"] = noise_mc
        out[f"elam_{n_shots}"] = e_lam
        out[f"T_{n_shots}"] = T_arr
    core._atomic_savez(path, **out)


def _load(path, n_mc):
    if not os.path.exists(path):
        return None
    boosts, picks, n_mc = _cache_key(n_mc)
    z = np.load(path, allow_pickle=True)
    if "boosts" not in z.files:
        return None
    if not (np.allclose(z["boosts"], boosts) and int(z["n_mc"]) == n_mc):
        return None
    for n in (1, 4):
        if f"noise_{n}" not in z.files or not np.allclose(z[f"noise_{n}"], picks[n]):
            return None
    res = {}
    for n_shots in (1, 4):
        noises = z[f"noise_{n_shots}"]
        cnt = z[f"cnt_{n_shots}"]
        done = z[f"done_{n_shots}"]
        for j, nt in enumerate(noises):
            if not done[j]:
                continue
            T_map = {tag: int(z[f"T_{n_shots}"][ti, j])
                     for ti, tag in enumerate(z["far_tags"])}
            res[(n_shots, float(nt))] = {
                "cnt": cnt[:, j, :], "noise_target": float(nt),
                "noise_mc": float(z[f"noisemc_{n_shots}"][j]),
                "e_lambda": float(z[f"elam_{n_shots}"][j]),
                "T_map": T_map,
            }
    return res


def _mean_from_cnt(c):
    v = np.arange(c.size, dtype=float)
    return float((v * c).sum() / max(c.sum(), 1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=None)
    ap.add_argument("--n-mc", type=int, default=N_MC_DEFAULT)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--checkpoint-every", type=int, default=2)
    args = ap.parse_args()

    n_cpu = os.cpu_count() or 1
    workers = args.workers or min(getattr(core, "N_WORKERS", 20), n_cpu)
    n_mc = args.n_mc

    # 校验 noise 档都在网格上
    for n, lst in NOISE_PICK.items():
        for nt in lst:
            _noise_index(n, nt)

    res = None
    for cand in (CACHE, CACHE_CKPT):
        res = _load(cand, n_mc)
        if res is not None:
            print(f"从 {cand} 载入 {len(res)} 档，断点续跑")
            break
    if res is None:
        res = {}
        print("未找到可用缓存，从零开始扫描")

    pending = []
    for n_shots in (1, 4):
        for nt in NOISE_PICK[n_shots]:
            key = (n_shots, float(nt))
            if key in res:
                continue
            pending.append((n_shots, float(nt), n_mc,
                            SEED_BASE + n_shots * 1_000_000 + int(nt * 1000)))
    if args.limit:
        pending = pending[: args.limit]

    n_total = sum(len(v) for v in NOISE_PICK.values())
    if not pending:
        print(f"全部 {n_total} 档均已完成")
        _save(CACHE, res, n_mc)
        return

    print(f"peak_distribution 扫描：待跑 {len(pending)} / 共 {n_total} 档")
    print(f"每档 {len(BOOST_LIST)} 个 boost × {n_mc:,} 条 MC")
    print(f"boost = {BOOST_LIST[:3]} … {BOOST_LIST[-2:]}  （共 {len(BOOST_LIST)}）")
    print(f"E [nJ] ≈ {[round(b * core.E_PULSE_BASE * 1e9, 3) for b in BOOST_LIST[:3]]}"
          f" … {[round(b * core.E_PULSE_BASE * 1e9, 2) for b in BOOST_LIST[-2:]]}")
    print(f"并行：ProcessPoolExecutor {workers} 进程（CPU {n_cpu}）")

    os.environ["POD_CORE_QUIET"] = "1"
    t0, done_n = time.time(), 0
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futs = {pool.submit(_job, spec): spec for spec in pending}
        for fut in as_completed(futs):
            key, value = fut.result()
            res[key] = value
            done_n += 1
            if done_n % args.checkpoint_every == 0:
                _save(CACHE_CKPT, res, n_mc)
            el = time.time() - t0
            eta = el / done_n * (len(pending) - done_n)
            m0 = _mean_from_cnt(value["cnt"][0])
            m1 = _mean_from_cnt(value["cnt"][-1])
            print(f"  [{done_n}/{len(pending)}] N={key[0]} noise={key[1]:.2f}："
                  f"peak均值 boost0={m0:.2f} / boost最大={m1:.2f}；"
                  f"累计 {el/60:.1f} min，剩余 {eta/60:.1f} min")

    _save(CACHE, res, n_mc)
    print(f"扫描完成，总用时 {(time.time()-t0)/60:.1f} min；已写入 {CACHE}")
    if len(res) >= n_total and os.path.exists(CACHE_CKPT):
        try:
            os.remove(CACHE_CKPT)
        except OSError:
            pass


if __name__ == "__main__":
    sys.exit(main())
