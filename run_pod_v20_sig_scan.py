# -*- coding: utf-8 -*-
"""PoD_esti v20 模块 9.3 / 15 的信号扫描（固定信号强度 × 统一 bg 网格），多进程 + 缓存。

在 v11 里这一步是写在 notebook cell 里【串行】跑的：
    48 bg × 3 个 N × 9 个 boost × 8000 MC
单线程要跑很久，而且中断就全丢。本脚本把它拆成 (N, bg, boost) 三元组任务，
丢进 ProcessPool（规则三：默认 20 进程），每完成若干档写一次 .partial.npz 检查点。

产出缓存 `pod_esti_v20_cache_signal.npz`，字段与 notebook 的 CACHE_SIG 完全兼容：
    grid_key    : BG_GRID（48 档）
    boosts      : BOOST_LIST（9 档，含 boost=0 作为"无信号"基准）
    n_mc        : 每档 MC 条数
    n_shots_list: [1,2,4]
    peak_cnt_<N>: (n_boost, n_bg, 27*N+2) 的 peak bincount —— **完整分布**，
                  模块 15 要用它看"同信号、不同 bg 时 peak 分布怎么变"
    done_<N>    : (n_boost, n_bg) 布尔完成标记（本脚本新增，notebook 读取时忽略也无妨）

peak 的口径与模块 7 的 PoD 完全一致：在信号子窗 CENTERS_SIG 上取 hist_add 的最大值。
boost=0 那一行就是"同一个信号窗、同一个 bg、但没有信号"的对照组。

用法：
    $env:PYTHONIOENCODING="utf-8"
    python run_pod_v20_sig_scan.py --workers 20
    python run_pod_v20_sig_scan.py --limit 3 --n-mc 500     # 冒烟
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

try:                                  # v20 内核优先；缺失时退回 v11（两者物理部分逐字相同）
    import pod_esti_v20_core as core
except ImportError:
    import pod_esti_v11_core as core

CACHE = "pod_esti_v20_cache_signal.npz"
CACHE_CKPT = "pod_esti_v20_cache_signal.partial.npz"

N_LIST = [1, 2, 4]
# 与 v11 模块 9.3 相同的 boost 档位；boost=0 是无信号基准，模块 15 要用
BOOST_LIST = np.round(np.arange(0.0, 0.032 + 1e-12, 0.004), 6)


def _job(a):
    n_shots, k, bg_target, ib, boost, n_mc, chunk, seed0 = a
    nt_amb = float(bg_target) / float(n_shots)
    r_det = float(core.r_det_for_noise(nt_amb, core.N_PIX_MACRO))
    r_amb = r_det / core.PDE
    n_tr = core.N_PIX_MACRO * n_shots

    cnt = np.zeros(n_tr + 2, dtype=np.int64)
    done, part = 0, 0
    while done < n_mc:
        m = min(chunk, n_mc - done)
        rng = np.random.default_rng(seed0 + 7919 * part)
        hi = core.binary_macro_stepping_per_shot(
            m, core.F_VALS, n_shots, core.R_SIG_UNIT_POD, core.TF_POD,
            r_amb, core.CENTERS_SIG, rng, boost=float(boost),
        )
        pk = hi.sum(axis=1).max(axis=1).astype(np.int64)
        cnt += np.bincount(pk, minlength=cnt.size)
        done += m
        part += 1

    v = np.arange(cnt.size, dtype=float)
    tot = max(int(cnt.sum()), 1)
    mu = float((v * cnt).sum() / tot)
    sd = float(np.sqrt(max((v * v * cnt).sum() / tot - mu * mu, 0.0)))
    return dict(n_shots=int(n_shots), k=int(k), ib=int(ib),
                bg=float(bg_target), boost=float(boost),
                cnt=cnt, mean=mu, std=sd)


def _empty(nb, ng):
    return {n: {"peak_cnt": np.zeros((nb, ng, core.N_PIX_MACRO * n + 2),
                                     dtype=np.int64),
                "done": np.zeros((nb, ng), dtype=bool)} for n in N_LIST}


def _save(path, res, grid, n_mc):
    tmp = path + ".tmp.npz"
    payload = {"grid_key": grid, "boosts": np.asarray(BOOST_LIST, float),
               "n_mc": int(n_mc), "n_shots_list": np.asarray(N_LIST)}
    for n in N_LIST:
        payload[f"peak_cnt_{n}"] = res[n]["peak_cnt"]
        payload[f"done_{n}"] = res[n]["done"]
    np.savez_compressed(tmp, **payload)
    os.replace(tmp, path)


def _load(path, grid, n_mc):
    if not os.path.exists(path):
        return None
    try:
        z = np.load(path)
        if (int(z["n_mc"]) != int(n_mc)
                or not np.array_equal(z["n_shots_list"], np.asarray(N_LIST))
                or z["grid_key"].shape != grid.shape
                or not np.allclose(z["grid_key"], grid)
                or z["boosts"].shape != BOOST_LIST.shape
                or not np.allclose(z["boosts"], BOOST_LIST)):
            return None
        out = {}
        for n in N_LIST:
            cnt = np.array(z[f"peak_cnt_{n}"])
            if f"done_{n}" in z.files:
                dn = np.array(z[f"done_{n}"])
            else:
                dn = cnt.sum(axis=2) > 0          # 兼容 notebook 写出的旧缓存
            out[n] = {"peak_cnt": cnt, "done": dn}
        return out
    except Exception:
        return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=20)
    ap.add_argument("--n-mc", type=int, default=8_000)
    ap.add_argument("--chunk", type=int, default=1_000)
    ap.add_argument("--limit", type=int, default=0, help="只跑前 N 个 bg 档（冒烟用）")
    ap.add_argument("--checkpoint-every", type=int, default=20)
    args = ap.parse_args()

    grid = np.asarray(core.BG_GRID, float)
    ng, nb = len(grid), len(BOOST_LIST)

    print("=" * 84)
    print(f"v20 信号扫描：{ng} bg × {nb} boost × N={N_LIST} = {ng*nb*len(N_LIST)} 档，"
          f"每档 {args.n_mc:,} MC")
    print(f"boost 档位：{list(BOOST_LIST)}（boost=0 = 同窗口无信号基准）")
    print(f"信号子窗 {core.CENTERS_SIG.size} 个 bin，细网格 {core.TF_POD.size} 步")
    print("=" * 84)

    res = _load(CACHE, grid, args.n_mc) or _load(CACHE_CKPT, grid, args.n_mc)
    if res is None:
        res = _empty(nb, ng)
        print("未找到缓存，从零开始")
    else:
        nd = sum(int(res[n]["done"].sum()) for n in N_LIST)
        print(f"命中缓存，已完成 {nd}/{ng*nb*len(N_LIST)} 档")

    todo = [(n, k, float(grid[k]), ib, float(BOOST_LIST[ib]),
             args.n_mc, args.chunk, 31000 + 1009 * n + 37 * k + 7919 * ib)
            for n in N_LIST for k in range(ng) for ib in range(nb)
            if not res[n]["done"][ib, k] and not (args.limit and k >= args.limit)]
    # 单档耗时 ≈ n_tr × bg，重活先派发可缩短尾部空转
    todo.sort(key=lambda j: -(j[0] * max(j[2], 0.25)))
    print(f"待算 {len(todo)}；workers={args.workers}")

    if todo:
        t0 = time.time()
        dn = 0
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            futs = [pool.submit(_job, j) for j in todo]
            for fut in as_completed(futs):
                r = fut.result()
                n, k, ib = r["n_shots"], r["k"], r["ib"]
                res[n]["peak_cnt"][ib, k] = r["cnt"]
                res[n]["done"][ib, k] = True
                dn += 1
                el = time.time() - t0
                eta = el / dn * (len(todo) - dn)
                print(f"  [{dn}/{len(todo)} {100*dn/len(todo):5.1f}%] "
                      f"N={n} bg={r['bg']:5.2f} boost={r['boost']:.3f} → "
                      f"peakμ={r['mean']:6.2f} peakσ={r['std']:5.2f}"
                      f"　已用 {el/60:.1f} min，剩 {eta/60:.1f} min")
                if dn % args.checkpoint_every == 0 or dn == len(todo):
                    _save(CACHE_CKPT, res, grid, args.n_mc)
        _save(CACHE, res, grid, args.n_mc)
        if os.path.exists(CACHE_CKPT):
            try:
                os.remove(CACHE_CKPT)
            except OSError:
                pass
        print(f"[信号扫描完成] → {CACHE}，{(time.time()-t0)/60:.1f} min")
    else:
        _save(CACHE, res, grid, args.n_mc)
        print(f"全部命中缓存 → {CACHE}")


if __name__ == "__main__":
    main()
