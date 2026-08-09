# -*- coding: utf-8 -*-
"""同 bg 下不同 N_shots 的 peak / 阈值差异 —— 解析模型数值验证。

模型：每个统计窗 bin 的计数 X ~ Binomial(n_tr=27N, p=bg/(27N))
      peak = max over M=152 个 bin（先按独立处理，再讨论相关修正）

产物：theory_peak_bg_multishot_fig.png + 控制台表格（供 markdown 文档引用）
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import brentq
from scipy.stats import binom, norm

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
rcParams["axes.unicode_minus"] = False

N_PIX = 27          # 宏像元 SPAD 数
M_BINS = 152        # 统计窗 bin 数（0–200 ns 掐头去尾 24 ns）
N_LIST = [1, 2, 4]
BG_GRID = np.round(np.arange(0.25, 12.0 + 1e-9, 0.25), 4)
FARS = [(0.01, "1%"), (100e-6, "100ppm"), (10e-6, "10ppm")]


def peak_moments(n_tr, p, m=M_BINS):
    """独立 Binomial 的 M 重极值：E[peak], std[peak]。

    E[max] = Σ_{t>=0} (1 - F(t)^M)；二阶矩用 E[X²]=Σ (2t+1)(1-F(t)^M)。
    """
    t = np.arange(0, n_tr + 1)
    F = binom.cdf(t, n_tr, p)
    Fm = F ** m
    tail = 1.0 - Fm                      # P(peak > t)
    e1 = tail.sum()
    e2 = ((2 * t + 1) * tail).sum()
    return e1, float(np.sqrt(max(e2 - e1 * e1, 0.0)))


def thr_indep(n_tr, p, far, m=M_BINS):
    """满足 P(peak >= T) <= far 的最小整数 T。"""
    a_bin = 1.0 - (1.0 - far) ** (1.0 / m)
    t = np.arange(0, n_tr + 2)
    sf = binom.sf(t - 1, n_tr, p)        # P(X >= t)
    ok = np.where(sf <= a_bin)[0]
    return int(t[ok[0]]) if ok.size else n_tr + 1


def thr_continuous(n_tr, p, far, m=M_BINS):
    """连续阈值：在 ln SF(T) 上线性插值求 SF=α_bin 的实数 T。

    整数 T 会让 ρ=T_N/T_1 变成台阶，掩盖真实趋势；连续阈值用于看趋势。
    """
    a_bin = 1.0 - (1.0 - far) ** (1.0 / m)
    t = np.arange(0, n_tr + 2)
    sf = binom.sf(t - 1, n_tr, p)
    ls = np.log(np.maximum(sf, 1e-300))
    tgt = np.log(a_bin)
    idx = np.where(ls <= tgt)[0]
    if idx.size == 0 or idx[0] == 0:
        return np.nan
    j = idx[0]
    return float(t[j - 1] + (ls[j - 1] - tgt) / (ls[j - 1] - ls[j]))


def thr_poisson(bg, far, m=M_BINS):
    """N→∞ 极限（纯泊松）的连续阈值：解 T·ln(T/bg) − T + bg = L。"""
    L = -np.log(1.0 - (1.0 - far) ** (1.0 / m))
    f = lambda T: T * np.log(T / bg) - T + bg - L
    return float(brentq(f, bg * (1 + 1e-9) + 1e-9, bg + 400.0))


def thr_ld(bg, n_shots, far, m=M_BINS):
    """大偏差一阶修正：T_N ≈ T∞ − (T∞−bg)² / (2·27N·ln(T∞/bg))。"""
    t_inf = thr_poisson(bg, far, m)
    dt = (t_inf - bg) ** 2 / (2.0 * N_PIX * n_shots * np.log(t_inf / bg))
    return t_inf - dt, t_inf, dt


def thr_gauss(bg, n_tr, far, m=M_BINS):
    """高斯极值闭式：T ≈ bg + z·σ，σ=√(bg(1-bg/n_tr))。"""
    a_bin = 1.0 - (1.0 - far) ** (1.0 / m)
    z = norm.ppf(1.0 - a_bin)
    sd = np.sqrt(bg * (1.0 - bg / n_tr))
    return bg + z * sd, z, sd


print("=" * 78)
print("同 bg 下 peak / 阈值的解析模型（独立 Binomial 极值，M=152 bin）")
print("=" * 78)
print(f"n_tr = 27·N；同 bg 时 p = bg/(27N)；μ_bin = bg 与 N 无关")
print(f"Var_bin = bg·(1 − bg/(27N))  ⇒ N 越大方差越大（欠离散越弱）\n")

# ---------------- 表 1：方差比 ----------------
print("表 1　单 bin 标准差 σ_bin = √(bg(1−bg/27N))")
print(f"{'bg':>6} | " + " | ".join(f"N={n}:σ" for n in N_LIST) + " | σ4/σ1")
for bg in [1.0, 3.0, 6.0, 9.0, 12.0]:
    sds = [np.sqrt(bg * (1 - bg / (N_PIX * n))) for n in N_LIST]
    print(f"{bg:6.2f} | " + " | ".join(f"{s:8.4f}" for s in sds)
          + f" | {sds[-1]/sds[0]:.4f}")

# ---------------- 表 2：阈值与倍数 ----------------
results = {}
for far, flab in FARS:
    for n in N_LIST:
        n_tr = N_PIX * n
        results[(far, n)] = np.array(
            [thr_indep(n_tr, bg / n_tr, far) for bg in BG_GRID], float)

print("\n表 2　阈值 T 与倍数 ρ=T_N/T_1（精确二项极值）")
for far, flab in FARS:
    print(f"\n  FAR={flab}")
    print(f"  {'bg':>6} | {'T(N=1)':>7} {'T(N=2)':>7} {'T(N=4)':>7} | "
          f"{'ρ2':>6} {'ρ4':>6}")
    for bg in [0.5, 1.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0]:
        k = int(np.argmin(np.abs(BG_GRID - bg)))
        T = [results[(far, n)][k] for n in N_LIST]
        print(f"  {BG_GRID[k]:6.2f} | {T[0]:7.0f} {T[1]:7.0f} {T[2]:7.0f} | "
              f"{T[1]/T[0]:6.3f} {T[2]/T[0]:6.3f}")

# ---------------- 表 3：peak 均值 / std ----------------
print("\n表 3　peak 均值与标准差（M=152 独立 bin 极值）")
print(f"  {'bg':>6} | " + " ".join(f"{'μ'+str(n):>7}" for n in N_LIST)
      + " | " + " ".join(f"{'s'+str(n):>6}" for n in N_LIST)
      + " | Δμ(4−1) rel%")
pm = {n: [] for n in N_LIST}
ps = {n: [] for n in N_LIST}
for bg in BG_GRID:
    for n in N_LIST:
        n_tr = N_PIX * n
        e, s = peak_moments(n_tr, bg / n_tr)
        pm[n].append(e)
        ps[n].append(s)
for n in N_LIST:
    pm[n] = np.array(pm[n])
    ps[n] = np.array(ps[n])
for bg in [0.5, 1.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0]:
    k = int(np.argmin(np.abs(BG_GRID - bg)))
    d = pm[4][k] - pm[1][k]
    print(f"  {BG_GRID[k]:6.2f} | " + " ".join(f"{pm[n][k]:7.3f}" for n in N_LIST)
          + " | " + " ".join(f"{ps[n][k]:6.3f}" for n in N_LIST)
          + f" | {d:+7.3f} {100*d/pm[1][k]:+6.2f}%")

# ---------------- 高斯闭式对照 ----------------
print("\n表 4　高斯闭式 T≈bg+z·σ 与精确二项对比（FAR=1%）")
far = 0.01
print(f"  {'bg':>6} | {'z':>5} | {'T1_g':>6} {'T1_ex':>6} | {'T4_g':>6} {'T4_ex':>6} | "
      f"{'ρ4_g':>6} {'ρ4_ex':>6}")
for bg in [1.0, 3.0, 6.0, 9.0, 12.0]:
    k = int(np.argmin(np.abs(BG_GRID - bg)))
    t1g, z, _ = thr_gauss(bg, N_PIX * 1, far)
    t4g, _, _ = thr_gauss(bg, N_PIX * 4, far)
    t1e, t4e = results[(far, 1)][k], results[(far, 4)][k]
    print(f"  {bg:6.2f} | {z:5.2f} | {t1g:6.2f} {t1e:6.0f} | {t4g:6.2f} {t4e:6.0f} | "
          f"{t4g/t1g:6.3f} {t4e/t1e:6.3f}")

# ---------------- 有效 bin 数敏感性 ----------------
print("\n表 5　有效独立 bin 数 M_eff 的影响（FAR=1%，bg=6）")
bg = 6.0
for m_eff in [152, 76, 38, 19]:
    t1 = thr_indep(27, bg / 27, far, m=m_eff)
    t4 = thr_indep(108, bg / 108, far, m=m_eff)
    print(f"  M_eff={m_eff:4d} → T1={t1:3d}  T4={t4:3d}  ρ4={t4/t1:.3f}")

# ---------------- 表 6：连续阈值（去掉整数台阶）----------------
cont = {}
for far, flab in FARS:
    for n in N_LIST:
        n_tr = N_PIX * n
        cont[(far, n)] = np.array(
            [thr_continuous(n_tr, bg / n_tr, far) for bg in BG_GRID], float)

print("\n表 6　连续阈值与平滑倍数（整数量化已去除）")
for far, flab in FARS:
    print(f"\n  FAR={flab}")
    print(f"  {'bg':>6} | {'Tc1':>6} {'Tc2':>6} {'Tc4':>6} | {'ρ2':>6} {'ρ4':>6} "
          f"| {'T4−T1':>6}")
    for bg in [0.5, 1.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0]:
        k = int(np.argmin(np.abs(BG_GRID - bg)))
        T = [cont[(far, n)][k] for n in N_LIST]
        print(f"  {BG_GRID[k]:6.2f} | {T[0]:6.2f} {T[1]:6.2f} {T[2]:6.2f} | "
              f"{T[1]/T[0]:6.3f} {T[2]/T[0]:6.3f} | {T[2]-T[0]:+6.2f}")

# ---------------- 表 7：大偏差闭式 vs 精确 ----------------
print("\n表 7　大偏差闭式 T_N ≈ T∞ − (T∞−bg)²/(2·27N·ln(T∞/bg))　（FAR=1%）")
far = 0.01
print(f"  {'bg':>6} | {'T∞':>6} | {'T1_LD':>6} {'T1_ex':>6} | "
      f"{'T4_LD':>6} {'T4_ex':>6} | {'ΔT_LD':>6} {'ΔT_ex':>6}")
for bg in [1.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0]:
    k = int(np.argmin(np.abs(BG_GRID - bg)))
    t1l, t_inf, _ = thr_ld(bg, 1, far)
    t4l, _, _ = thr_ld(bg, 4, far)
    t1e, t4e = cont[(far, 1)][k], cont[(far, 4)][k]
    print(f"  {bg:6.2f} | {t_inf:6.2f} | {t1l:6.2f} {t1e:6.2f} | "
          f"{t4l:6.2f} {t4e:6.2f} | {t4l-t1l:+6.2f} {t4e-t1e:+6.2f}")

# ---------------- 表 8：二值硬上限风险 ----------------
print("\n表 8　N=1 的二值硬上限（cap=27）风险：T 距离 cap 还有多少")
for far, flab in FARS:
    row = []
    for bg in [8.0, 10.0, 11.0, 12.0]:
        k = int(np.argmin(np.abs(BG_GRID - bg)))
        row.append(f"bg={bg:>4.1f}:T1={results[(far, 1)][k]:>3.0f}(余{27-results[(far,1)][k]:>3.0f})")
    print(f"  FAR={flab:>7}  " + "  ".join(row))

# ---------------- 作图 ----------------
colors = {1: "tab:blue", 2: "tab:green", 4: "tab:red"}
fig, ax = plt.subplots(2, 2, figsize=(12.5, 9.0))

a = ax[0, 0]
for n in N_LIST:
    a.plot(BG_GRID, np.sqrt(BG_GRID * (1 - BG_GRID / (N_PIX * n))),
           color=colors[n], lw=2, label=f"N={n}")
a.plot(BG_GRID, np.sqrt(BG_GRID), "k:", lw=1.5, label="纯泊松 √bg")
a.set_xlabel("bg"); a.set_ylabel(r"$\sigma_{bin}$")
a.set_title("① 单 bin 标准差：N 越大越接近泊松")
a.grid(alpha=0.3); a.legend()

a = ax[0, 1]
for n in N_LIST:
    a.plot(BG_GRID, results[(0.01, n)], color=colors[n], lw=2, label=f"N={n}")
a.set_xlabel("bg"); a.set_ylabel("T @ FAR=1%")
a.set_title("② 阈值曲线（精确二项极值）")
a.grid(alpha=0.3); a.legend()

a = ax[1, 0]
for far, flab in FARS:
    for n in [2, 4]:
        ls = "-" if n == 4 else "--"
        a.plot(BG_GRID, cont[(far, n)] / cont[(far, 1)],
               ls, lw=1.8, label=f"N={n}, FAR={flab}")
a.plot(BG_GRID, results[(0.01, 4)] / results[(0.01, 1)], color="0.6", lw=1.0,
       alpha=0.8, label="N=4, FAR=1%（整数台阶）")
a.axhline(1.0, color="k", ls=":", lw=1.2)
a.set_xlabel("bg"); a.set_ylabel(r"$\rho_{N/1}=T_N/T_1$")
a.set_title("③ 阈值倍数（连续阈值）：不是常数，随 bg 单调上升")
a.grid(alpha=0.3); a.legend(fontsize=8)

a = ax[1, 1]
for n in N_LIST:
    a.plot(BG_GRID, pm[n], color=colors[n], lw=2, label=f"N={n} 均值")
    a.fill_between(BG_GRID, pm[n] - ps[n], pm[n] + ps[n],
                   color=colors[n], alpha=0.15)
a.set_xlabel("bg"); a.set_ylabel("peak")
a.set_title("④ peak 均值 ±1σ（阴影）")
a.grid(alpha=0.3); a.legend()

fig.suptitle("同 bg 下 N_shots 对 peak 分布与阈值的影响（解析模型）", fontsize=13)
fig.tight_layout()
fig.savefig("theory_peak_bg_multishot_fig.png", dpi=130, bbox_inches="tight")
print("\n图已保存 → theory_peak_bg_multishot_fig.png")

# ---------------- 供文档引用的汇总 ----------------
print("\n汇总（FAR=1%）：")
r4 = results[(0.01, 4)] / results[(0.01, 1)]
r2 = results[(0.01, 2)] / results[(0.01, 1)]
print(f"  ρ2 范围 {r2.min():.3f} – {r2.max():.3f}")
print(f"  ρ4 范围 {r4.min():.3f} – {r4.max():.3f}")
print(f"  bg≤2 时 ρ4 ≈ {r4[BG_GRID <= 2].mean():.3f}；"
      f"bg≥9 时 ρ4 ≈ {r4[BG_GRID >= 9].mean():.3f}")
dmu = (pm[4] - pm[1])
print(f"  Δpeak_mean(4−1) 范围 {dmu.min():+.2f} – {dmu.max():+.2f} "
      f"（相对 {100*(dmu/pm[1]).min():+.1f}% – {100*(dmu/pm[1]).max():+.1f}%）")
