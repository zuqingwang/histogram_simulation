# -*- coding: utf-8 -*-
"""测量统计窗内 bin 与 bin 的相关性，定出有效独立 bin 数 M_eff。

背景：引擎里「bin 被点亮 ⟺ 它之前最近一次雪崩距它 < T_OVER(≈8 ns)」，
所以一次雪崩会点亮**一连串**相邻 bin（直到下一次雪崩截断）。
这会造成相邻 bin 之间的【正相关】，而不是死时间式的负相关。

输出：
  - 单 bin 边缘方差（跨 MC）vs 二项预言 bg(1−bg/27N)
  - 单条 hist 内 152 bin 的样本方差均值
  - 自相关函数 ACF(lag)
  - 由 E[s²]/σ² 与由 peak 反推的 M_eff
"""
from __future__ import annotations

import argparse
import os

os.environ.setdefault("POD_CORE_QUIET", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np
from scipy.stats import binom

import pod_esti_v11_core as core

I0, I1 = core.I_STAT0, core.I_STAT1
M = core.N_STAT


def peak_moments_indep(n_tr, p, m):
    t = np.arange(0, n_tr + 1)
    F = binom.cdf(t, n_tr, p)
    tail = 1.0 - F ** m
    e1 = tail.sum()
    e2 = ((2 * t + 1) * tail).sum()
    return e1, float(np.sqrt(max(e2 - e1 * e1, 0.0)))


def m_eff_from_peak(n_tr, p, target_peak_mean, lo=1.0, hi=400.0):
    """二分：让独立模型的 peak 均值等于实测值，反推有效独立 bin 数。"""
    for _ in range(60):
        mid = np.sqrt(lo * hi)
        if peak_moments_indep(n_tr, p, mid)[0] < target_peak_mean:
            lo = mid
        else:
            hi = mid
    return float(np.sqrt(lo * hi))


def run(bg, n_shots, n_mc, chunk):
    noise = bg / n_shots
    r_det = float(core.r_det_for_noise(noise, core.N_PIX_MACRO))
    inv_tab = core.build_inv_table(r_det)
    n_tr = core.N_PIX_MACRO * n_shots
    p_eq = float(core.p_bin_equilibrium(r_det)[0])

    lags = np.arange(0, 21)
    acf_num = np.zeros(lags.size)
    tot_n = 0
    s2_sum = 0.0
    x_sum = 0.0
    x2_sum = 0.0
    pk = []
    done, part = 0, 0
    while done < n_mc:
        m = min(chunk, n_mc - done)
        rng = np.random.default_rng(31337 + 7919 * part)
        hi = core.noise_hists_per_shot(m, n_shots, r_det, rng, inv_tab=inv_tab)
        a = hi[:, :, I0:I1].sum(axis=1).astype(np.float64)
        s2_sum += float(a.var(axis=1, ddof=1).sum())
        x_sum += float(a.sum())
        x2_sum += float((a * a).sum())
        pk.append(a.max(axis=1))
        c = a - a.mean(axis=1, keepdims=True)
        for j, L in enumerate(lags):
            acf_num[j] += float((c[:, :M - L] * c[:, L:]).sum())
        tot_n += m
        done += m
        part += 1

    pk = np.concatenate(pk)
    mean_bin = x_sum / (tot_n * M)
    var_bin = x2_sum / (tot_n * M) - mean_bin ** 2      # 单 bin 边缘方差
    s2_mean = s2_sum / tot_n                             # 单条 hist 内样本方差
    acf = acf_num / np.array([tot_n * (M - L) for L in lags])
    acf = acf / acf[0]

    var_theory = bg * (1 - bg / n_tr)
    return dict(
        bg=bg, n_shots=n_shots, n_tr=n_tr, p_eq=p_eq,
        mean_bin=mean_bin, var_bin=var_bin, var_theory=var_theory,
        s2_mean=s2_mean, acf=acf, lags=lags,
        peak_mean=float(pk.mean()), peak_std=float(pk.std()),
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-mc", type=int, default=20000)
    ap.add_argument("--chunk", type=int, default=5000)
    args = ap.parse_args()

    print("=" * 92)
    print(f"bin 间相关诊断（统计窗 M={M} bins，T_OVER={core.T_OVER*1e9:.2f} ns）")
    print("=" * 92)

    for bg, n_shots in [(4.0, 1), (4.0, 4), (12.0, 1), (12.0, 4), (1.0, 1), (1.0, 4)]:
        r = run(bg, n_shots, args.n_mc, args.chunk)
        m_pk = m_eff_from_peak(r["n_tr"], r["p_eq"], r["peak_mean"])
        pk_indep = peak_moments_indep(r["n_tr"], r["p_eq"], M)[0]
        # 由 E[s²]=σ²[1 − Σρ/(M(M−1))] 反推平均相关
        rho_bar = (1 - r["s2_mean"] / r["var_bin"])
        print(f"\n--- bg={bg:g}, N={n_shots}  (noise={bg/n_shots:.3f}, "
              f"p_eq={r['p_eq']:.4f}) ---")
        print(f"  单 bin 均值        实测 {r['mean_bin']:.4f}   目标 {bg:g}")
        print(f"  单 bin 边缘方差    实测 {r['var_bin']:.4f}   二项预言 "
              f"{r['var_theory']:.4f}   比 {r['var_bin']/r['var_theory']:.4f}")
        print(f"  单条 hist 内样本方差 {r['s2_mean']:.4f}   "
              f"比边缘方差 {r['s2_mean']/r['var_bin']:.4f}")
        print(f"  ⇒ 平均成对相关 ρ̄ ≈ {rho_bar/(1):.5f}"
              f"（>0 表示正相关；乘 M(M−1) 即 Σρ）")
        print(f"  ACF: " + "  ".join(
            f"L{L}={v:+.3f}" for L, v in zip(r["lags"][:11], r["acf"][:11])))
        print(f"  peak 均值 实测 {r['peak_mean']:.3f}   "
              f"独立模型(M=152) {pk_indep:.3f}")
        print(f"  ⇒ 由 peak 反推有效独立 bin 数 M_eff ≈ {m_pk:.1f}"
              f"（名义 152，比值 {m_pk/M:.3f}）")


if __name__ == "__main__":
    main()
