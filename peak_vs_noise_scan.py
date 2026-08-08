# -*- coding: utf-8 -*-
"""peak_vs_noise v01 —— 固定信号强度下，peak 分布随环境噪声增强的演化扫描。

要回答的问题
    给定信号强度（固定 boost），环境噪声 noise 逐渐增强时，
    信号窗内 peak（宏像元 1 ns bin 的最大累加计数）的分布怎么变？是线性的吗？

noise 的口径与 PoD_esti 全项目一致：宏像元（9×3 = 27 个 SPAD）在 1 ns bin 上的平均累加计数。

设计要点
  1. 物理内核直接 import pod_esti_v05_core，不复制任何物理参数，保证与 PoD_esti_v05 一致。
  2. boost = 0 的那一档就是**纯噪声基线**，用来把 peak 拆成「噪声本底 + 信号净增量」。
  3. 额外算 noise = 0（r_amb = 0）的**无噪声纯信号**参考档，用于可加性检验：
     若 peak(信号+噪声) < peak(纯信号) + peak(纯噪声)，即为次可加（sub-additive）→ 非线性。
  4. peak 取值域只有 0…n_tr（≤108），所以 **bincount 就是充分统计量**：
     整个扫描的缓存只有约 1 MB，却保留了完整分布，事后可算任意分位数、画任意分布图。
  5. 用 ProcessPoolExecutor 而不是线程：binary_macro_stepping 是
     「Python 层逐细网格步循环 + 小数组 NumPy」，几乎全程持有 GIL
     （Global Interpreter Lock，全局解释器锁），多线程吃不满 CPU。

用法（PowerShell）
    $env:PYTHONIOENCODING="utf-8"
    python peak_vs_noise_scan.py                # 全量 208 档
    python peak_vs_noise_scan.py --limit 6      # 冒烟测试
    python peak_vs_noise_scan.py --workers 16 --n-mc 4000

产物：peak_vs_noise_v01_cache.npz（增量检查点 peak_vs_noise_v01_cache.partial.npz）
"""
import argparse
import functools
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

# 长任务必须实时刷新，否则 stdout 缓冲会让终端一片空白
print = functools.partial(print, flush=True)  # noqa: A001

import numpy as np

# 模块级 import：spawn 出来的子进程会 import 本文件，从而共享同一份内核构建流程
import pod_esti_v05_core as core

# ============================ 本工作的参数 ============================
# 固定信号强度档（boost = 相对 E_PULSE_BASE = 799.4 nJ 的倍数），
# 对应单脉冲能量 0 / 3.20 / 6.40 / 12.79 / 25.58 nJ。
# 第 0 档必须是 0.0 —— 那是纯噪声基线。
# 取值依据见 worklog_peak_vs_noise.md：boost ≳ 0.05 时 peak 已顶到二值硬上限，
# 噪声的影响会被硬上限压掉，看不出规律，所以只取 0.004…0.032。
BOOST_LIST = [0.0, 0.004, 0.008, 0.016, 0.032]

N_MC_DEFAULT = 8000     # 每个 (n_shots, noise, boost) 的 MC 条数
MC_CHUNK     = 2000     # 单次 binary_macro_stepping 的条数，控制峰值内存
SEED_BASE    = 91_000

CACHE      = "peak_vs_noise_v01_cache.npz"
CACHE_CKPT = "peak_vs_noise_v01_cache.partial.npz"

K_NOFLOOR = -1          # 特殊档号：noise = 0（r_amb = 0），无噪声纯信号参考


# ============================ 子进程任务 ============================
def _peak_bincount(boost, n_shots, r_amb, n_mc, seed, n_tr):
    """跑 n_mc 条 MC，返回 peak 的 bincount（长度 n_tr + 2）。"""
    cnt = np.zeros(n_tr + 2, dtype=np.int64)
    for s in range(0, n_mc, MC_CHUNK):
        m = min(MC_CHUNK, n_mc - s)
        pk = core._peaks_chunk(boost, n_shots, r_amb, m, seed + 104_729 * s)
        cnt += np.bincount(pk, minlength=n_tr + 2)[: n_tr + 2]
    return cnt


def _job(args):
    """单个 (n_shots, noise 档) 任务：把该档下所有 boost 一次算完。

    k = K_NOFLOOR 时表示 noise = 0（r_amb = 0）的无噪声纯信号参考档。
    """
    n_shots, k, n_mc, seed0 = args
    n_tr = 27 * n_shots

    if k == K_NOFLOOR:
        r_amb, noise_t, noise_mc, e_lam = 0.0, 0.0, 0.0, 0.0
    else:
        R = core.NOISE_RES[n_shots]
        noise_t = float(R["noise_target"][k])
        noise_mc = float(R["noise_mc"][k])
        e_lam = float(R["e_lambda"][k])
        r_amb = float(R["r_det"][k] / core.PDE)

    cnt = np.zeros((len(BOOST_LIST), n_tr + 2), dtype=np.int64)
    for i, boost in enumerate(BOOST_LIST):
        cnt[i] = _peak_bincount(boost, n_shots, r_amb, n_mc,
                                seed0 + 1009 * i, n_tr)
    return (n_shots, k), {
        "cnt": cnt, "noise_target": noise_t,
        "noise_mc": noise_mc, "e_lambda": e_lam, "r_amb": r_amb,
    }


# ============================ 缓存读写 ============================
def _cache_key(n_mc):
    """缓存键：noise 网格 + 信号档 + MC 条数。任一项变化都不得复用旧缓存。"""
    return (np.concatenate([core.NOISE_GRID[n] for n in core.N_SHOTS_LIST]),
            np.asarray(BOOST_LIST, float), int(n_mc))


def _save(path, res, n_mc):
    grid_key, boosts, n_mc = _cache_key(n_mc)
    out = {"grid_key": grid_key, "boosts": boosts, "n_mc": n_mc,
           "far_tags": np.array(core.FAR_TAGS)}
    for n_shots in core.N_SHOTS_LIST:
        grid = core.NOISE_GRID[n_shots]
        n_tr = 27 * n_shots
        nb, ng = len(BOOST_LIST), len(grid)
        cnt = np.zeros((nb, ng, n_tr + 2), dtype=np.int64)
        done = np.zeros(ng, dtype=bool)
        noise_mc = np.full(ng, np.nan)
        e_lam = np.full(ng, np.nan)
        for k in range(ng):
            r = res.get((n_shots, k))
            if r is None:
                continue
            cnt[:, k, :] = r["cnt"]
            done[k] = True
            noise_mc[k] = r["noise_mc"]
            e_lam[k] = r["e_lambda"]
        out[f"noise_{n_shots}"] = grid
        out[f"cnt_{n_shots}"] = cnt
        out[f"done_{n_shots}"] = done
        out[f"noisemc_{n_shots}"] = noise_mc
        out[f"elam_{n_shots}"] = e_lam
        # 无噪声纯信号参考档
        r0 = res.get((n_shots, K_NOFLOOR))
        out[f"cnt0_{n_shots}"] = (r0["cnt"] if r0 is not None
                                 else np.zeros((nb, n_tr + 2), dtype=np.int64))
        out[f"done0_{n_shots}"] = np.array(r0 is not None)
        # 各档 FAR 阈值，随缓存一起存，让 notebook 自包含
        Tr = core.THRESH[n_shots]
        out[f"T_{n_shots}"] = np.stack(
            [Tr["T" + tag] for tag in core.FAR_TAGS]).astype(int)
    core._atomic_savez(path, **out)


def _load(path, n_mc):
    """只有网格、信号档、MC 条数全部一致才接受缓存。"""
    if not os.path.exists(path):
        return None
    grid_key, boosts, n_mc = _cache_key(n_mc)
    z = np.load(path, allow_pickle=True)
    if "boosts" not in z.files:
        return None
    if not (np.array_equal(z["grid_key"], grid_key)
            and np.allclose(z["boosts"], boosts)
            and int(z["n_mc"]) == n_mc):
        return None
    res = {}
    for n_shots in core.N_SHOTS_LIST:
        cnt = z[f"cnt_{n_shots}"]
        done = z[f"done_{n_shots}"]
        grid = z[f"noise_{n_shots}"]
        for k in range(len(grid)):
            if not done[k]:
                continue
            res[(n_shots, k)] = {
                "cnt": cnt[:, k, :], "noise_target": float(grid[k]),
                "noise_mc": float(z[f"noisemc_{n_shots}"][k]),
                "e_lambda": float(z[f"elam_{n_shots}"][k]), "r_amb": np.nan,
            }
        if bool(z[f"done0_{n_shots}"]):
            res[(n_shots, K_NOFLOOR)] = {
                "cnt": z[f"cnt0_{n_shots}"], "noise_target": 0.0,
                "noise_mc": 0.0, "e_lambda": 0.0, "r_amb": 0.0,
            }
    return res


# ============================ 主流程 ============================
def _mean_from_cnt(c):
    v = np.arange(c.size, dtype=float)
    return (v * c).sum() / max(c.sum(), 1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=None)
    ap.add_argument("--n-mc", type=int, default=N_MC_DEFAULT)
    ap.add_argument("--limit", type=int, default=0,
                    help="只跑前 N 档（冒烟测试），0 = 全量")
    ap.add_argument("--checkpoint-every", type=int, default=10)
    args = ap.parse_args()

    n_cpu = os.cpu_count() or 1
    workers = args.workers or min(getattr(core, "N_WORKERS", 20), n_cpu)
    n_mc = args.n_mc

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
    for n_shots in core.N_SHOTS_LIST:
        # 无噪声纯信号参考档排在最前，便于尽早看到基准
        for k in [K_NOFLOOR] + list(range(len(core.NOISE_GRID[n_shots]))):
            if (n_shots, k) in res:
                continue
            pending.append((n_shots, k, n_mc,
                            SEED_BASE + n_shots * 1_000_000 + (k + 2) * 20_000))
    if args.limit:
        pending = pending[: args.limit]

    n_total = sum(len(core.NOISE_GRID[n]) + 1 for n in core.N_SHOTS_LIST)
    if not pending:
        print(f"全部 {n_total} 档均已完成，无需计算")
        _save(CACHE, res, n_mc)
        return

    print(f"peak 分布扫描：待跑 {len(pending)} / 共 {n_total} 档")
    print(f"每档 {len(BOOST_LIST)} 个信号强度 × {n_mc:,} 条 MC "
          f"= {len(BOOST_LIST)*n_mc:,} 条")
    print(f"信号档 boost = {BOOST_LIST}")
    print(f"      即 E = {[round(b*core.E_PULSE_BASE*1e9, 2) for b in BOOST_LIST]} nJ"
          f"（boost=0 为纯噪声基线）")
    print(f"并行：ProcessPoolExecutor {workers} 进程（本机 CPU {n_cpu}）")

    os.environ["POD_CORE_QUIET"] = "1"   # 子进程 import 内核时静音

    t0, done_n = time.time(), 0
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futs = {pool.submit(_job, spec): spec for spec in pending}
        for fut in as_completed(futs):
            key, value = fut.result()
            res[key] = value
            done_n += 1
            if done_n % args.checkpoint_every == 0:
                _save(CACHE_CKPT, res, n_mc)
            if done_n == 1 or done_n % 10 == 0 or done_n == len(pending):
                el = time.time() - t0
                eta = el / done_n * (len(pending) - done_n)
                nt = value["noise_target"]
                m_noise = _mean_from_cnt(value["cnt"][0])
                m_top = _mean_from_cnt(value["cnt"][-1])
                print(f"  [{done_n}/{len(pending)}] N_shots={key[0]} "
                      f"noise={nt:.2f}：peak均值 纯噪声={m_noise:.2f}，"
                      f"最强信号档={m_top:.2f}；"
                      f"累计 {el/60:.1f} min，预计剩余 {eta/60:.1f} min")

    _save(CACHE, res, n_mc)
    print(f"扫描完成，总用时 {(time.time()-t0)/60:.1f} min；已写入 {CACHE}")
    if len(res) >= n_total and os.path.exists(CACHE_CKPT):
        try:
            os.remove(CACHE_CKPT)
        except OSError:
            pass


if __name__ == "__main__":
    sys.exit(main())
