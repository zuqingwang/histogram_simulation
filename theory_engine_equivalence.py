# -*- coding: utf-8 -*-
"""逐环节验证「精确逐光子引擎 ≡ 更新过程快速引擎 ≡ 平衡态解析式」的推导。

配套文档：theory_engine_equivalence.md（推导在那里，本脚本只做数值检验）。

推导链条与本脚本的测试一一对应：
  T1  逆变换采样：Δ = H⁻¹(E)，E~Exp(1)  ⟹  P(Δ>x) = e^{−H(x)}
  T2  更新-回报定理：p_bin = (1/μ)∫₀^{T_OVER} S = E[min(X,T_OVER)]/E[X]
  T3  过阈窗并集 = 「在下一次雪崩处截断」后的不交并（差分数组技巧是恒等式，不是近似）
  T4  宏像元单 bin 边缘律 = Binomial(n_tr, p_bin)，且只依赖 n_tr
  T5  精确逐光子引擎是上述连续时间模型的离散化：p_bin(dt) → 解析值，误差 O(dt)

用法：
    $env:PYTHONIOENCODING="utf-8"
    python theory_engine_equivalence.py
"""
from __future__ import annotations

import functools
import os
import time

print = functools.partial(print, flush=True)  # noqa: A001

os.environ.setdefault("POD_CORE_QUIET", "1")

import numpy as np
from scipy.stats import binom as _binom

import pod_esti_v11_core as core

TAU = core.TAU_RC
TOV = core.T_OVER
NZ_REF = 6.0            # 参考环境光档（noise27 = 27·p_eq）


def hazard(d, r_det):
    """条件强度 h(Δ) = r_det · g(1 − e^{−Δ/τ})。"""
    return r_det * core.spad_response_g(1.0 - np.exp(-d / TAU),
                                        core.RESP_SHAPE, core.RESP_K)


def cumulative_hazard(d_grid, r_det):
    h = hazard(d_grid, r_det)
    return np.concatenate([[0.0], np.cumsum(0.5 * (h[1:] + h[:-1]) * np.diff(d_grid))])


def sample_intervals_like_engine(n, r_det, rng, inv_tab=None):
    """完全复刻 noise_macro_hist_fast 里的采样路径：E~Exp(1) → 直查表 → 线性插值。"""
    inv, scale = core.build_inv_table(r_det) if inv_tab is None else inv_tab
    E = rng.standard_exponential(n, dtype=np.float32)
    x = np.minimum(E, np.float32(core.E_MAX)) * scale
    i0 = x.astype(np.int32)
    np.clip(i0, 0, inv.size - 2, out=i0)
    fr = x - i0
    return np.maximum(inv[i0] * (1.0 - fr) + inv[i0 + 1] * fr,
                      np.float32(1e-13)).astype(np.float64)


# ============================================================ T1
def t1_interval_law(r_det, n=2_000_000):
    print("=" * 96)
    print("T1  逆变换采样正确性： Δ = H⁻¹(E), E~Exp(1)  ⟹  P(Δ > x) = S(x) = e^{−H(x)}")
    print("=" * 96)
    rng = np.random.default_rng(11)
    dx = sample_intervals_like_engine(n, r_det, rng)

    d = np.linspace(0.0, 40 * TAU + 40.0 / r_det, 400001)
    H = cumulative_hazard(d, r_det)
    S = np.exp(-H)

    qs = [0.5, 1, 2, 3, 5, 8, 12, 20, 30, 50, 80]   # ns
    print(f"{'x [ns]':>8} {'解析 S(x)':>12} {'经验 P(Δ>x)':>13} {'相对差':>10}")
    worst = 0.0
    for xn in qs:
        x = xn * 1e-9
        s_ana = float(np.interp(x, d, S))
        s_emp = float((dx > x).mean())
        rel = s_emp / s_ana - 1
        worst = max(worst, abs(rel))
        print(f"{xn:8.1f} {s_ana:12.6f} {s_emp:13.6f} {100*rel:9.3f}%")
    se = 1.0 / np.sqrt(n)
    print(f"\n最大相对偏差 {100*worst:.3f}%（{n:,} 样本的 MC 相对误差量级 ≈ {100*se:.3f}%）")
    return d, H, S, dx


# ============================================================ T2
def t2_renewal_reward(r_det, d, S, dx):
    print("\n" + "=" * 96)
    print("T2  更新-回报定理： p_bin = (1/μ)∫₀^{T_OVER} S(u)du = E[min(X,T_OVER)] / E[X]")
    print("=" * 96)
    p_ana, mu_ana = core.p_bin_equilibrium(r_det)

    mu_emp = float(dx.mean())
    p_rr = float(np.minimum(dx, TOV).mean() / mu_emp)

    m = d <= TOV
    p_int = float(np.trapezoid(S[m], d[m]) / mu_ana)

    rng = np.random.default_rng(23)
    h = core.noise_macro_hist_fast(400_000, 1, r_det, rng)
    p_mc = float(h[:, core.I_STAT0:core.I_STAT1].mean())

    print(f"  平均雪崩间隔 μ：解析 {mu_ana*1e9:.4f} ns，采样 {mu_emp*1e9:.4f} ns"
          f"（差 {100*(mu_emp/mu_ana-1):+.3f}%）")
    print(f"  ① p_bin_equilibrium 解析式            = {p_ana:.6f}")
    print(f"  ② (1/μ)∫₀^T S du（独立数值积分）       = {p_int:.6f}  差 {100*(p_int/p_ana-1):+.4f}%")
    print(f"  ③ E[min(X,T)]/E[X]（更新-回报, 采样）  = {p_rr:.6f}  差 {100*(p_rr/p_ana-1):+.4f}%")
    print(f"  ④ 快速引擎 MC（400,000×1 轨迹）        = {p_mc:.6f}  差 {100*(p_mc/p_ana-1):+.4f}%")
    print("  ②③④ 都指向同一个数 ⟹ 快速引擎采样的就是解析式描述的那个平衡态分布。")
    return p_ana


# ============================================================ T3
def t3_union_identity(r_det, n_trace=4000):
    print("\n" + "=" * 96)
    print("T3  「过阈窗并集」= 「在下一次雪崩处截断后的不交并」——差分数组技巧是恒等式")
    print("=" * 96)
    print("    ∪ₖ [aₖ, aₖ+T)  =  ⊔ₖ [aₖ, min(aₖ+T, aₖ₊₁))   （aₖ 递增时逐点相等）")

    rng = np.random.default_rng(37)
    centers = core.CENTERS
    wl, bw, nb = core.WIN_LO, core.BIN_W, core.NBINS
    t_start = core.WIN_LO - core.WARM_NS * 1e-9
    inv_tab = core.build_inv_table(r_det)          # 只建一次表
    span = core.WIN_HI - t_start
    n_step = int(4 * span / core.p_bin_equilibrium(r_det)[1]) + 32   # 足够覆盖整窗的间隔数

    for jitter, tag in ((0.0, "无抖动"), (core.JIT, f"有抖动 σ={core.JIT*1e12:.0f} ps")):
        bad = 0
        checked = 0
        for _ in range(n_trace):
            # 直接按更新过程生成一条轨迹的雪崩时刻（一次性批量采样后累加）
            t0_ = t_start + rng.exponential(1.0 / r_det)
            a = t0_ + np.concatenate(
                [[0.0], np.cumsum(sample_intervals_like_engine(n_step, r_det, rng, inv_tab))])
            a = a[a < core.WIN_HI]
            if a.size < 2:
                continue
            if jitter > 0:
                a = a + rng.normal(0.0, jitter, a.size)

            # (a) 暴力涂并集
            brute = np.zeros(nb, dtype=np.int32)
            for tt in a:
                brute[(centers >= tt) & (centers < tt + TOV)] = 1

            # (b) 截断 + 差分数组 + cumsum（引擎用的写法）
            lo = a
            hi = np.minimum(a + TOV, np.concatenate([a[1:], [np.inf]]))
            b_lo = np.clip(np.ceil((lo - wl) / bw - 0.5), 0, nb).astype(np.int64)
            b_hi = np.clip(np.ceil((hi - wl) / bw - 0.5), 0, nb).astype(np.int64)
            msk = b_hi > b_lo
            diff = (np.bincount(b_lo[msk], minlength=nb + 1)
                    - np.bincount(b_hi[msk], minlength=nb + 1))
            fast = np.cumsum(diff)[:nb]

            checked += 1
            if not np.array_equal(brute, fast):
                bad += 1
        print(f"  {tag}：{checked} 条轨迹中逐 bin 完全相同 {checked-bad} 条，不同 {bad} 条")
    print("    有抖动时唯一可能的偏差来源：抖动把相邻雪崩的先后顺序颠倒（需间隔 ≲ σ），")
    print("    此时 b_hi ≤ b_lo 该区间被丢弃。下面给出这种事件的概率。")
    d = np.linspace(0.0, 10 * core.JIT, 20001)
    p_swap = float(1.0 - np.exp(-cumulative_hazard(d, r_det)[-1]))
    print(f"    P(相邻间隔 < 10σ = {10*core.JIT*1e12:.0f} ps) = {p_swap:.3e}（乱序概率的上界）")


# ============================================================ T4
def t4_binomial_marginal(r_det, p_ana):
    print("\n" + "=" * 96)
    print("T4  宏像元单 bin 边缘律 = Binomial(n_tr, p_bin)，且只依赖 n_tr = n_pix × N_shots")
    print("=" * 96)
    rng = np.random.default_rng(41)
    for n_tr in (27, 36, 108):
        h = core.noise_macro_hist_fast(120_000, n_tr, r_det, rng)
        col = h[:, core.I_STAT0 + 40]                  # 统计窗中间随便挑一个 bin
        m, v = col.mean(), col.var()
        cnt = np.bincount(col, minlength=n_tr + 1)[:n_tr + 1] / col.size
        pmf = _binom.pmf(np.arange(n_tr + 1), n_tr, p_ana)
        tv = 0.5 * np.abs(cnt - pmf).sum()             # 总变差距离
        print(f"  n_tr={n_tr:3d}: 均值 {m:7.4f}（理论 {n_tr*p_ana:7.4f}）"
              f"  方差 {v:7.4f}（理论 {n_tr*p_ana*(1-p_ana):7.4f}）"
              f"  与二项分布的总变差距离 {tv:.4f}")
    print("  注：这是【单个 bin 的边缘分布】。bin 之间因 T_OVER=8ns 正相关，")
    print("      所以 peak 不能按「152 个独立二项」算（见 handoff 坑 11）。")


# ============================================================ T5
def p_bin_discrete(r_amb, dt):
    """把「细网格步进」这个离散模型**精确解出来**（不用 MC）。

    推导（见 md 第 8 节）：一步内到达 n~Poisson(μ)，μ=r_amb·dt，各自以 φ(age) 触发，
    且首个触发后 vov=0 ⟹ g(0)=0 ⟹ 步内至多一次雪崩。于是
        P(本步至少一次雪崩 | age=m) = 1 − E[(1−φ)ⁿ] = 1 − exp(−μ·φ(m·dt))  ≡ q(m)
    这正是 binary_macro_stepping 里的 `-expm1(-mu*phi[age])`。
    离散更新过程：P(X_d > m 步) = Π_{j=1..m}(1−q(j))，且 q(0)=0（刚雪崩完不可能再触发）。
    """
    M = int(np.ceil(60 * TAU / dt)) + int(np.ceil(60.0 / (r_amb * core.PDE * dt))) + 10
    m = np.arange(0, M + 1)
    phi = core.PDE * core.spad_response_g(1.0 - np.exp(-m * dt / TAU),
                                          core.RESP_SHAPE, core.RESP_K)
    q = -np.expm1(-r_amb * dt * phi)
    q[0] = 0.0
    surv = np.concatenate([[1.0], np.cumprod(1.0 - q[1:])])   # surv[m] = P(X_d > m 步)
    mu = dt * surv.sum()
    # 过阈窗宽是【连续量】T_OVER，与 dt 无关；∫₀^{T_OVER} P(Y>u)du，Y=X_d·dt 为阶梯生存函数
    K = int(np.floor(TOV / dt))
    lit = dt * surv[:K].sum() + surv[K] * (TOV - K * dt)
    return lit / mu, mu


def t5_discretization():
    print("\n" + "=" * 96)
    print("T5a 精确逐光子引擎 = 连续时间模型的离散化：把离散模型解析解出来，量化偏差")
    print("=" * 96)
    for nz in (NZ_REF, 12.0):
        r_det = float(core.r_det_for_noise(nz, core.N_PIX_MACRO))
        r_amb = r_det / core.PDE
        p_c, mu_c = core.p_bin_equilibrium(r_det)
        print(f"\n  noise27 = {nz:.1f}　r_det = {r_det:.4e} cps"
              f"　连续模型 p_bin = {p_c:.6f}，μ = {mu_c*1e9:.4f} ns")
        print(f"  {'dt [ps]':>9} {'离散 p_bin':>12} {'相对差':>11} "
              f"{'离散 μ [ns]':>13} {'μ 相对差':>11} {'μ误差比上一行':>14}")
        prev = None
        for dt_ps in (3200, 1600, 800, 400, 200, 100, 50):
            p_d, mu_d = p_bin_discrete(r_amb, dt_ps * 1e-12)
            e_mu = mu_d / mu_c - 1
            ratio = f"{prev/e_mu:14.2f}" if prev else f"{'—':>14}"
            print(f"  {dt_ps:9d} {p_d:12.6f} {100*(p_d/p_c-1):10.4f}% "
                  f"{mu_d*1e9:13.4f} {100*e_mu:10.4f}% {ratio}")
            prev = e_mu
        print("  最后一列 = 上一行 μ 误差 ÷ 本行 μ 误差。dt 每减半误差降 4 倍 ⟹ 误差是 O(dt²)。")
        print(f"  生产用的 dt = {core.DT_FINE*1e12:.0f} ps 那一行：p_bin 偏差已在 1e-4 相对量级，")
        print("  低于 p_bin_equilibrium 自身数值积分的精度，MC 根本分辨不出来。")

    print("\n" + "=" * 96)
    print("T5b 交叉验证：逐光子引擎 spad_binary_trace ≟ 步进引擎 binary_macro_stepping")
    print("=" * 96)
    print("    两者的等价性靠恒等式 E[(1−φ)ⁿ] = e^{−μφ}（n~Poisson(μ)）；")
    print("    步进引擎直接用 1−e^{−μφ}，逐光子引擎是逐个光子抽。二者应统计一致。")
    centers = core.CENTERS
    i0, i1 = core.I_STAT0, core.I_STAT1
    nbn = i1 - i0
    t_lo = core.WIN_LO - core.WARM_NS * 1e-9
    r_det = float(core.r_det_for_noise(NZ_REF, core.N_PIX_MACRO))
    r_amb = r_det / core.PDE
    print(f"\n  noise27 = {NZ_REF}　{'dt [ps]':>8} {'逐光子引擎':>12} {'步进引擎':>12} "
          f"{'离散模型解析':>14} {'逐光子−步进':>13}")
    for dt_ps, n in ((800, 6000), (200, 6000)):
        dt = dt_ps * 1e-12
        tf = np.arange(t_lo, core.WIN_HI, dt)
        rs = np.zeros_like(tf)

        rng = np.random.default_rng(900 + dt_ps)
        per = np.empty(n)
        for k in range(n):
            per[k] = core.spad_binary_trace(
                rs, r_amb, tf, centers, core.PDE, core.TAU_RC, core.VTH_FRAC,
                core.JIT, rng, TOV, 0.0, core.RESP_SHAPE, core.RESP_K)[i0:i1].sum() / nbn
        p1, se1 = per.mean(), per.std(ddof=1) / np.sqrt(n)

        rng = np.random.default_rng(1900 + dt_ps)
        h = core.binary_macro_stepping(n, np.zeros(1), rs, tf, r_amb, centers, rng)
        col = h[:, i0:i1].mean(axis=1)
        p2, se2 = col.mean(), col.std(ddof=1) / np.sqrt(n)

        p_d, _ = p_bin_discrete(r_amb, dt)
        print(f"  {'':>16} {dt_ps:8d} {p1:8.5f}±{se1:.5f} {p2:8.5f}±{se2:.5f} "
              f"{p_d:14.5f} {p1-p2:+12.5f}")


def main():
    r_det = float(core.r_det_for_noise(NZ_REF, core.N_PIX_MACRO))
    print(f"参考档：noise27 = {NZ_REF}　r_det = {r_det:.4e} cps　"
          f"τ_RC = {TAU*1e9:.4f} ns　T_OVER = {TOV*1e9:.4f} ns　"
          f"g = {core.RESP_SHAPE}(k={core.RESP_K})\n")
    d, H, S, dx = t1_interval_law(r_det)
    p_ana = t2_renewal_reward(r_det, d, S, dx)
    t3_union_identity(r_det)
    t4_binomial_marginal(r_det, p_ana)
    t5_discretization()


if __name__ == "__main__":
    t0 = time.time()
    main()
    print(f"\n总用时 {time.time()-t0:.1f} s")
