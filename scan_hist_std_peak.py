# -*- coding: utf-8 -*-
"""扫描三个量随噪声的变化（N=1/2/4），出三联图。

三个量（都在 hist_add 的统计窗内，纯环境光、无信号）：
  ① within-hist std：**单条**直方图内 152 个 bin 的样本标准差，
     再对所有 MC 条数取平均。注意这是「一条 hist 内部的起伏」，
     不是「同一个 bin 在多条 hist 之间的起伏」——理想独立时两者相等。
  ② peak 均值
  ③ peak 标准差（跨 MC 条数）

扫描轴：统一 BG_GRID（步长 0.25）。每档令单发底 noise = bg / N，
        所以同一横坐标 bg 上三个 N 严格可比；另出一版横轴为 noise 的图。

缓存 + 多进程（项目规则三）。用法：
    $env:PYTHONIOENCODING="utf-8"
    python scan_hist_std_peak.py --workers 20 --n-mc 100000
    python scan_hist_std_peak.py --limit 4 --n-mc 5000     # 冒烟
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

try:                                  # v20 起优先用 v20 内核；两者物理部分逐字相同
    import pod_esti_v20_core as core
except ImportError:
    import pod_esti_v11_core as core

CACHE = "scan_hist_std_peak_cache.npz"
CACHE_CKPT = "scan_hist_std_peak_cache.partial.npz"
N_LIST = [1, 2, 4]
COLORS = {1: "tab:blue", 2: "tab:green", 4: "tab:red"}


def _job(args):
    """一个 (N, bg) 档：返回 within-hist std 的均值 + peak 的 bincount。"""
    n_shots, k, bg_target, n_mc, chunk, seed0 = args
    i0, i1 = core.I_STAT0, core.I_STAT1
    noise = float(bg_target) / float(n_shots)
    r_det = float(core.r_det_for_noise(noise, core.N_PIX_MACRO))
    inv_tab = core.build_inv_table(r_det)

    std_sum = 0.0
    bg_sum = 0.0
    peak_cnt = np.zeros(core.N_PIX_MACRO * n_shots + 2, dtype=np.int64)
    nn = 0
    done, part = 0, 0
    while done < n_mc:
        m = min(chunk, n_mc - done)
        rng = np.random.default_rng(seed0 + 7919 * part)
        hi = core.noise_hists_per_shot(m, n_shots, r_det, rng, inv_tab=inv_tab)
        add = hi[:, :, i0:i1].sum(axis=1).astype(np.float64)   # (m, N_STAT)
        std_sum += float(add.std(axis=1, ddof=1).sum())        # 每条 hist 各自的 std
        bg_sum += float(add.mean(axis=1).sum())
        peak_cnt += np.bincount(add.max(axis=1).astype(np.int64),
                                minlength=peak_cnt.size)
        nn += m
        done += m
        part += 1

    return {
        "n_shots": int(n_shots), "k": int(k),
        "bg_target": float(bg_target), "noise": noise,
        "hist_std_mean": std_sum / max(nn, 1),
        "bg_mc": bg_sum / max(nn, 1),
        "peak_cnt": peak_cnt, "nn": nn,
    }


def _peak_mean_std(cnt):
    cnt = np.asarray(cnt, float)
    n = cnt.sum()
    if n <= 0:
        return np.nan, np.nan
    t = np.arange(cnt.size, dtype=float)
    mu = float((cnt * t).sum() / n)
    var = float((cnt * t * t).sum() / n - mu * mu)
    return mu, float(np.sqrt(max(var, 0.0)))


def _save(path, res, grid, n_mc):
    tmp = path + ".tmp.npz"
    np.savez_compressed(
        tmp, grid=grid, n_mc=n_mc, n_list=np.asarray(N_LIST),
        **{f"{key}_{n}": res[n][key] for n in N_LIST
           for key in ("hist_std", "bg_mc", "peak_mean", "peak_std", "done")})
    os.replace(tmp, path)


def _load(path, grid, n_mc):
    if not os.path.exists(path):
        return None
    try:
        z = np.load(path)
        if int(z["n_mc"]) != int(n_mc) or not np.allclose(z["grid"], grid):
            return None
        return {n: {key: np.array(z[f"{key}_{n}"]) for key in
                    ("hist_std", "bg_mc", "peak_mean", "peak_std", "done")}
                for n in N_LIST}
    except Exception:
        return None


def _empty(ng):
    return {n: {"hist_std": np.zeros(ng), "bg_mc": np.zeros(ng),
                "peak_mean": np.zeros(ng), "peak_std": np.zeros(ng),
                "done": np.zeros(ng, dtype=bool)} for n in N_LIST}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=20)
    ap.add_argument("--n-mc", type=int, default=100_000)
    ap.add_argument("--chunk", type=int, default=5_000)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--checkpoint-every", type=int, default=8)
    args = ap.parse_args()

    grid = np.asarray(core.BG_GRID, float)
    ng = len(grid)
    res = _load(CACHE, grid, args.n_mc) or _load(CACHE_CKPT, grid, args.n_mc)
    if res is None:
        res = _empty(ng)
        print("未找到缓存，从零开始")
    else:
        nd = sum(int(res[n]["done"].sum()) for n in N_LIST)
        print(f"命中缓存，已完成 {nd}/{ng*len(N_LIST)} 档")

    todo = [(n, k, float(grid[k]), args.n_mc, args.chunk, 4200 + 1009 * n + 31 * k)
            for n in N_LIST for k in range(ng)
            if not res[n]["done"][k] and not (args.limit and k >= args.limit)]

    print(f"BG {ng} 档 × N={N_LIST} × {args.n_mc:,} MC；待算 {len(todo)}；"
          f"workers={args.workers}")

    if todo:
        t0 = time.time()
        dn = 0
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            futs = [pool.submit(_job, j) for j in todo]
            for fut in as_completed(futs):
                r = fut.result()
                n, k = r["n_shots"], r["k"]
                mu, sd = _peak_mean_std(r["peak_cnt"])
                res[n]["hist_std"][k] = r["hist_std_mean"]
                res[n]["bg_mc"][k] = r["bg_mc"]
                res[n]["peak_mean"][k] = mu
                res[n]["peak_std"][k] = sd
                res[n]["done"][k] = True
                dn += 1
                el = time.time() - t0
                eta = el / dn * (len(todo) - dn)
                print(f"  [{dn}/{len(todo)} {100*dn/len(todo):5.1f}%] "
                      f"N={n} bg={r['bg_target']:.2f}（noise={r['noise']:.3f}）→ "
                      f"histσ={r['hist_std_mean']:.3f} peakμ={mu:.3f} peakσ={sd:.3f}"
                      f"　已用 {el/60:.1f} min，剩 {eta/60:.1f} min")
                if dn % args.checkpoint_every == 0 or dn == len(todo):
                    _save(CACHE_CKPT, res, grid, args.n_mc)
        _save(CACHE, res, grid, args.n_mc)
        if os.path.exists(CACHE_CKPT):
            try:
                os.remove(CACHE_CKPT)
            except OSError:
                pass
        print(f"[扫描完成] → {CACHE}，{(time.time()-t0)/60:.1f} min")

    # ---------------- 作图 ----------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import rcParams
    rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
    rcParams["axes.unicode_minus"] = False

    ok = {n: res[n]["done"] for n in N_LIST}

    def draw(xkey, fname, xlabel, suptitle):
        fig, ax = plt.subplots(1, 3, figsize=(16.5, 4.9))

        for n in N_LIST:
            msk = ok[n]
            x = grid[msk] if xkey == "bg" else grid[msk] / n
            ax[0].plot(x, res[n]["hist_std"][msk], "-o", color=COLORS[n],
                       ms=2.6, lw=1.8, label=f"N={n}　MC 实测")
            bgv = grid[msk]
            ax[0].plot(x, np.sqrt(bgv * (1 - bgv / (27 * n))), "--",
                       color=COLORS[n], lw=1.1, alpha=0.75,
                       label=f"N={n}　解析 √(bg(1−bg/27N))")
        if xkey == "bg":
            ax[0].plot(grid, np.sqrt(grid), "k:", lw=1.4, label="纯泊松 √bg")
        ax[0].set_xlabel(xlabel)
        ax[0].set_ylabel("单条 hist 内 152 bin 的 std（对 MC 取平均）")
        ax[0].set_title("① 每条 hist 自身的起伏")
        ax[0].legend(fontsize=7.2)
        ax[0].grid(alpha=0.3)

        for n in N_LIST:
            msk = ok[n]
            x = grid[msk] if xkey == "bg" else grid[msk] / n
            ax[1].plot(x, res[n]["peak_mean"][msk], "-o", color=COLORS[n],
                       ms=2.6, lw=1.8, label=f"N={n}")
            ax[1].plot(x, res[n]["bg_mc"][msk], ":", color=COLORS[n], lw=1.0,
                       alpha=0.7)
        ax[1].set_xlabel(xlabel)
        ax[1].set_ylabel("peak 均值 [计数]")
        ax[1].set_title("② peak 均值（点线 = 实测 bg，作参照）")
        ax[1].legend(fontsize=8)
        ax[1].grid(alpha=0.3)

        for n in N_LIST:
            msk = ok[n]
            x = grid[msk] if xkey == "bg" else grid[msk] / n
            ax[2].plot(x, res[n]["peak_std"][msk], "-o", color=COLORS[n],
                       ms=2.6, lw=1.8, label=f"N={n}")
        ax[2].axvline(np.nan)
        ax[2].set_xlabel(xlabel)
        ax[2].set_ylabel("peak 标准差 [计数]")
        ax[2].set_title("③ peak 标准差（跨 MC 条数）")
        ax[2].legend(fontsize=8)
        ax[2].grid(alpha=0.3)

        fig.suptitle(suptitle, fontsize=12.5)
        fig.tight_layout()
        fig.savefig(fname, dpi=130, bbox_inches="tight")
        print(f"图已保存 → {fname}")

    draw("bg", "scan_hist_std_peak_vs_bg.png", "bg（hist_add 统计窗均值）",
         f"三联图 vs bg —— 同 bg 下比 N（noise=bg/N，每档 {args.n_mc:,} 次 MC）")
    draw("noise", "scan_hist_std_peak_vs_noise.png", "noise（单发 hist_i 统计窗均值）",
         f"三联图 vs noise —— 单发底相同时比 N（bg=N·noise，每档 {args.n_mc:,} 次 MC）")

    # ---------------- 抽样表 ----------------
    print("\n抽样数值（每 8 档一行）")
    print(f"{'bg':>6} | " + " | ".join(
        f"N={n}: histσ  peakμ  peakσ" for n in N_LIST))
    for k in range(0, ng, 8):
        if not all(ok[n][k] for n in N_LIST):
            continue
        row = " | ".join(
            f"{res[n]['hist_std'][k]:6.3f} {res[n]['peak_mean'][k]:6.2f} "
            f"{res[n]['peak_std'][k]:6.3f}" for n in N_LIST)
        print(f"{grid[k]:6.2f} | {row}")


if __name__ == "__main__":
    main()
