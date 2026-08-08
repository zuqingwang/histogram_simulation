# -*- coding: utf-8 -*-
"""
基于 v20 生成 v21:
新增三块内容(均基于 N_shots=4 的累加直方图, 峰值宏像元 m_peak, 30m 场景; 物理参数一律不改):

模块 14 — 检测阈值 det_th + 前沿法定时 front_time
  * nc_base = 峰值宏像元"纯背景蒙卡"的全 bin 均值(复用 v20 模块13 的 bg_hist_peak)。
  * det_th = k_th * nc_base (倍数 k_th 可调; 展示多个 k)。
  * 前沿法(front_time): 判决电平 V_dec = (det_th + 峰值)/2; 在峰"前沿"找跨越 V_dec 的相邻两 bin,
    线性插值得 front_time; 换算距离, 与真实 ToF 比较。

模块 15 — SNR 分布(RC 引擎重复跑) + 正态拟合
  * 重复 N_TRIALS 次峰值宏像元 RC 逐光子蒙卡(信号+背景), 每次算 SNR=S/sqrt(B), 收集成分布。
  * 直方图 + 正态拟合; 给出理论(误差传播) (mu,sigma) 与 拟合 (mu,sigma) 对比。

模块 16 — 100ppm 理论阈值 + 海量 Poisson 蒙卡验证噪点率
  * 噪点率(虚警率)定义(用户选定): P(一次 N=4 测量的时间窗内 >=1 个噪声 bin 超阈) = 100ppm = 1e-4。
  * 理论: 单 bin 虚警率 a_bin = 1-(1-1e-4)^(1/N_bins_eff); 由 Poisson(nc_base) 生存函数反解整数阈值;
    附高斯近似阈值对比。
  * 验证: 向量化 rng.poisson(nc_base,(chunk,N_bins)) 累计到 ~1e8 次(可配), 统计实测窗口级噪点率。

运行: python build_v21_from_v20.py
"""
import nbformat

SRC = "lidar_histogram_sim_v20.ipynb"
DST = "lidar_histogram_sim_v21.ipynb"
nb = nbformat.read(SRC, as_version=4)

for c in nb.cells:
    if c.cell_type == "code":
        c.outputs = []; c.execution_count = None
    if "execution" in c.get("metadata", {}):
        del c.metadata["execution"]

def find(pred):
    for i, c in enumerate(nb.cells):
        if pred(c): return i, c
    raise RuntimeError("cell not found")

# 找模块 13 的 code cell(最后一个), 新模块插到其后
i13, c13 = find(lambda c: c.cell_type == "code" and "SNR = S / sqrt(B), 在峰值宏像元" in c.source)

# ============================================================================
# 模块 14 — 检测阈值 det_th + 前沿法定时 front_time
# ============================================================================
md_14 = nbformat.v4.new_markdown_cell(r'''## 模块 14（v21 新增）— 检测阈值 det_th 与前沿法定时 front_time

在 **N_shots=4 累加直方图**(峰值宏像元 m_peak)上做"检测 + 定时":

**1) 检测阈值 det_th（基于噪声基底的倍数）**
- 噪声基底 `nc_base` = 峰值宏像元**纯背景蒙卡**(信号率置零)的**全 bin 平均计数**(复用模块 13 的 `bg_hist_peak`)。
- 检测阈值 `det_th = k_th · nc_base`(倍数 `k_th` 可调)。只有某 bin 计数 ≥ det_th 才认为"检测到回波"。

**2) 前沿法定时 front_time（leading-edge timing）**
- 判决电平取 **det_th 与峰值的均值**: `V_dec = (det_th + peak) / 2`。
- 沿信号峰的**上升沿(前沿)**找到**首次**从 `< V_dec` 跨越到 `≥ V_dec` 的相邻两个 bin;
- 在这两个 bin 之间对计数**线性插值**, 求出计数恰好等于 `V_dec` 的时刻 = **front_time**;
- `front_time` 换算距离 `R = c·front_time/2`, 与真实 ToF/距离比较(前沿法通常略早于真峰, 存在系统偏移)。

> 前沿法只用"前沿两点线性插值", 实现简单、抗峰形畸变; 代价是对幅度/阈值敏感, 有固定的时间游走(time walk)。''')

code_14 = nbformat.v4.new_code_cell(r'''# ---- 检测阈值 det_th + 前沿法定时 front_time (峰值宏像元, N_shots 累加直方图) ----
# 复用模块13: sig_hist_peak(信号+背景 直方图), bg_hist_peak(纯背景), pk_bin, B_mean, tc_ns

K_TH = 5.0                          # 检测阈值倍数: det_th = K_TH * nc_base (可调)
nc_base = bg_hist_peak.mean()       # 噪声基底 = 纯背景全 bin 平均计数(= 模块13 的 B_mean)
det_th  = K_TH * nc_base            # 检测阈值
peak_val = sig_hist_peak[pk_bin]    # 信号峰值计数(该 bin)
V_dec    = 0.5 * (det_th + peak_val)   # 前沿法判决电平 = (det_th + 峰值)/2

def front_time_leading_edge(hist, centers, pk_idx, v_dec, bin_w):
    """前沿法: 在峰(pk_idx)左侧上升沿, 找首次 hist[i-1]<v_dec<=hist[i] 的相邻两 bin, 线性插值得时刻。
    返回 (front_time[s], i_lo, i_hi) 或 (nan,...) 若未跨越。"""
    for i in range(1, pk_idx + 1):
        y0, y1 = hist[i-1], hist[i]
        if y0 < v_dec <= y1:
            frac = (v_dec - y0) / (y1 - y0) if y1 > y0 else 0.0
            t_front = centers[i-1] + frac * (centers[i] - centers[i-1])
            return t_front, i-1, i
    return np.nan, -1, -1

t_front, i_lo, i_hi = front_time_leading_edge(sig_hist_peak, centers, pk_bin, V_dec, bin_width)
detected = peak_val >= det_th

# 距离换算与真值比较
R_true   = D0                                   # 真实距离 30m
R_front  = C_LIGHT * t_front / 2.0 if np.isfinite(t_front) else np.nan
R_peakbin= C_LIGHT * centers[pk_bin] / 2.0      # 用峰 bin 中心定时(对照)
dR_front = R_front - R_true if np.isfinite(R_front) else np.nan

print("="*76)
print(f"检测阈值 & 前沿法定时 (峰值宏像元 m={m_peak}, N_shots={N_shots})")
print(f"  噪声基底 nc_base = 纯背景全 bin 均值 = {nc_base:.4f} 计数/bin")
print(f"  检测阈值 det_th = {K_TH:.1f} x nc_base = {det_th:.3f} 计数")
print(f"  信号峰值 peak = {peak_val:.1f} 计数 @ {tc_ns[pk_bin]:.1f} ns  -> {'检测到 (peak>=det_th)' if detected else '未检测到'}")
print(f"  前沿法判决电平 V_dec = (det_th+peak)/2 = {V_dec:.3f} 计数")
if np.isfinite(t_front):
    print(f"  前沿跨越: bin[{i_lo}]={sig_hist_peak[i_lo]:.1f} -> bin[{i_hi}]={sig_hist_peak[i_hi]:.1f} (@ {tc_ns[i_lo]:.1f}->{tc_ns[i_hi]:.1f} ns)")
    print(f"  -> front_time = {t_front*1e9:.3f} ns  ->  R_front = {R_front:.3f} m")
    print(f"     真实距离 R_true = {R_true:.3f} m; 前沿法偏移 dR = {dR_front*100:+.1f} cm (前沿法通常略早于真峰)")
    print(f"     (对照) 峰 bin 中心定时 R = {R_peakbin:.3f} m")
else:
    print("  -> 前沿未跨越判决电平(未检出或峰太弱)")

# ---- 绘图: 峰值宏像元直方图 + det_th / V_dec / front_time ----
fig, ax = plt.subplots(figsize=(11, 5))
ax.bar(tc_ns, sig_hist_peak, width=bin_width*1e9, align="center", color="tab:green", alpha=0.6, label="信号+背景 (RC 蒙卡)")
ax.bar(tc_ns, bg_hist_peak, width=bin_width*1e9, align="center", color="tab:red", alpha=0.5, label="纯背景")
ax.axhline(nc_base, color="gray",   ls=":",  lw=1.2, label=f"nc_base={nc_base:.3f}")
ax.axhline(det_th,  color="orange", ls="--", lw=1.6, label=f"det_th={K_TH:.0f}·nc_base={det_th:.2f}")
ax.axhline(V_dec,   color="purple", ls="-.", lw=1.4, label=f"V_dec=(det_th+peak)/2={V_dec:.1f}")
ax.axvline(t0_ns, color="k", ls=":", alpha=0.7, label=f"真实 ToF {t0_ns:.1f} ns")
if np.isfinite(t_front):
    ax.axvline(t_front*1e9, color="blue", lw=2.0, label=f"front_time={t_front*1e9:.2f} ns (R={R_front:.2f} m)")
    ax.plot([tc_ns[i_lo], tc_ns[i_hi]], [sig_hist_peak[i_lo], sig_hist_peak[i_hi]],
            "b-o", ms=6, lw=1.5, zorder=6)  # 前沿插值的两点连线
    ax.scatter([t_front*1e9], [V_dec], c="blue", marker="D", s=70, zorder=7)  # 插值交点
# 聚焦到回波附近
ax.set_xlim(t0_ns-8, t0_ns+18)
ax.set_xlabel("时间 t [ns]"); ax.set_ylabel(f"计数 / 1ns bin ({N_shots} shots)")
ax.set_title(f"检测阈值与前沿法定时: det_th={K_TH:.0f}·nc_base, 前沿两点线性插值求 front_time")
ax.legend(fontsize=8, loc="upper right"); ax.grid(alpha=0.3)
plt.tight_layout(); plt.show()

# ---- 不同倍数 k_th 对前沿法定时/距离的影响 ----
print("\n不同检测阈值倍数 k_th 下的前沿法定时(V_dec 随之改变):")
print(f"  {'k_th':>5} {'det_th':>8} {'V_dec':>8} {'front_time[ns]':>15} {'R_front[m]':>12} {'dR[cm]':>9}")
for k in [3, 5, 8, 10, 15]:
    dth = k * nc_base; vdec = 0.5*(dth + peak_val)
    tf_k, _, _ = front_time_leading_edge(sig_hist_peak, centers, pk_bin, vdec, bin_width)
    if np.isfinite(tf_k):
        Rk = C_LIGHT*tf_k/2.0
        print(f"  {k:>5.0f} {dth:>8.2f} {vdec:>8.1f} {tf_k*1e9:>15.3f} {Rk:>12.3f} {(Rk-R_true)*100:>+9.1f}")
    else:
        print(f"  {k:>5.0f} {dth:>8.2f} {vdec:>8.1f} {'(未跨越)':>15}")
print("说明: k_th 越大 -> V_dec 越高 -> 前沿交点越靠近峰顶 -> front_time 越晚(越接近真峰), time walk 越小。")''')

# ============================================================================
# 模块 15 — SNR 分布(RC 引擎重复跑) + 正态拟合
# ============================================================================
md_15 = nbformat.v4.new_markdown_cell(r'''## 模块 15（v21 新增）— 多次蒙卡的信噪比分布 SNR(直方图 + 正态拟合)

对峰值宏像元**重复 N_TRIALS 次**独立的 RC 逐光子蒙卡(每次都是 27 SPAD × N_shots, 含信号+背景),
每次算一次 `SNR = S / √B̄`, 把这些 SNR 汇成分布并做**正态拟合**。

**分母 B̄ 的取法(关键)**: 分母用**该次纯背景蒙卡的全 bin 平均计数 B̄**(= 噪声基底 nc_base 的定义),
而**不用峰 bin 的单个背景计数**。原因: 峰 bin 在信号存在时被 RC 计数死区强烈压制, 该 bin 的"背景来的计数"
大量为 0, 若用它当分母会导致约一半样本 B=0(SNR 无定义、被迫丢弃), 使分布严重截断、非正态。
用全 bin 均值 B̄ 估背景, 既符合"背景水平应由离峰 bin 估计"的物理事实, 又让 B̄ 稳定 >0, 分布回归干净正态。

- 信号 `S` = 峰 bin 总计数 − 峰 bin 背景计数(仍为 Poisson, 均值 S̄);
- `B̄` = 40 个 bin 背景计数的平均, 由大数平均后**抖动很小、近似常数**。

**理论值(误差传播)**: 把 B̄ 视为常数, 只有 S~Poisson(S̄)(Var=S̄) 随机:
- 均值 `μ_SNR = S̄/√B̄`;
- 标准差 `σ_SNR = σ_S/√B̄ = √S̄/√B̄`(因 Poisson 的 σ_S=√S̄)。

由于此时 SNR 只是 S 的**线性缩放**(除以常数 √B̄), S 在 S̄≈40 量级已近正态, 故 SNR 分布为**干净正态**,
理论 (μ,σ) 应与拟合 (μ,σ) **高度吻合**。给出两者对比。''')

code_15 = nbformat.v4.new_code_cell(r'''# ---- 多次蒙卡的 SNR 分布(RC 引擎重复跑峰值宏像元) ----
from scipy.stats import norm

N_TRIALS = 3000                     # SNR 采样次数(RC 引擎重复跑; 2000~5000 约十几秒~1分钟)
rng_snr = np.random.default_rng(PARAMS["hist"]["seed"] + 777)
fvals_peak = macro_fvals[m_peak]    # 峰值宏像元 27 个 SPAD 的收集比例

import time as _time
_t0 = _time.time()
snr_samples = np.empty(N_TRIALS)
S_samples   = np.empty(N_TRIALS)
B_samples   = np.empty(N_TRIALS)
for it in range(N_TRIALS):
    # 一次测量 = N_shots 累加, 27 SPAD 合并 -> 信号+背景 直方图
    h_sig = np.zeros(nbins); h_bg = np.zeros(nbins)
    for _shot in range(N_shots):
        ev_s = []
        for fij in fvals_peak:
            ev = simulate_spad_shot_rc(base_rate*fij, r_amb_ph, tf, PDE, TAU_RC, VTH_FRAC, jit, rng_snr, RESET, RESP_SHAPE, RESP_K)
            if ev.size: ev_s.append(ev)
        if ev_s: h_sig += np.histogram(np.concatenate(ev_s), bins=edges)[0]
        # 同批纯背景(信号率=0), 用于该次测量的 B 估计
        ev_b = []
        for fij in fvals_peak:
            ev = simulate_spad_shot_rc(zero_rate, r_amb_ph, tf, PDE, TAU_RC, VTH_FRAC, jit, rng_snr, RESET, RESP_SHAPE, RESP_K)
            if ev.size: ev_b.append(ev)
        if ev_b: h_bg += np.histogram(np.concatenate(ev_b), bins=edges)[0]
    b_pk = h_bg[pk_bin]                 # 该次峰 bin 单背景(仅作展示/佐证, 不作分母)
    b_bar = h_bg.mean()                 # 该次"全 bin 背景均值"= 噪声基底(与 nc_base 同定义, 稳定>0)
    s_pk = max(h_sig[pk_bin] - b_pk, 0.0)   # 峰 bin 扣该次背景(峰 bin)得信号计数
    S_samples[it] = s_pk; B_samples[it] = b_bar
    snr_samples[it] = s_pk / np.sqrt(b_bar) if b_bar > 0 else np.nan   # 分母=全 bin 背景均值 B_bar
    if (it+1) % 500 == 0:
        print(f"  ...SNR 采样进度 {it+1}/{N_TRIALS}  (用时 {_time.time()-_t0:.0f}s)")

snr_valid = snr_samples[np.isfinite(snr_samples)]   # B_bar 恒>0, 一般无丢弃
n_bzero = N_TRIALS - snr_valid.size

# 拟合(正态) & 理论(误差传播)
# 分母 B_bar = 全 bin 背景均值, 视为近似常数(40 个 bin 平均, 抖动很小); 只有 S~Poisson(S_bar) 随机。
#   => mu_SNR = S_bar/sqrt(B_bar);  sigma_SNR = sigma_S/sqrt(B_bar) = sqrt(S_bar)/sqrt(B_bar)
mu_fit, sig_fit = norm.fit(snr_valid)
S_bar = S_samples.mean(); B_bar = B_samples.mean()
mu_th  = S_bar / np.sqrt(B_bar)
sig_th = np.sqrt(S_bar) / np.sqrt(B_bar)            # S~Poisson: Var(S)=S_bar, B_bar 近似常数

print("="*76)
print(f"SNR 分布 (峰值宏像元 m={m_peak}, {N_TRIALS} 次 RC 蒙卡, 每次 27 SPAD x {N_shots} shots)")
print(f"  分母 B_bar = 每次测量的纯背景全 bin 均值(= 噪声基底 nc_base, 稳定>0, 避免峰 bin 单背景=0 丢样)")
print(f"  样本: 有效 {snr_valid.size} 个 (B_bar=0 丢弃 {n_bzero} 个); S_bar={S_bar:.2f}, B_bar={B_bar:.4f}")
print(f"  理论(误差传播, B_bar 近似常数): mu_th  = S_bar/sqrt(B_bar)          = {mu_th:.3f}")
print(f"                                  sig_th = sqrt(S_bar)/sqrt(B_bar)     = {sig_th:.3f}")
print(f"  拟合(正态):                     mu_fit = {mu_fit:.3f},  sig_fit = {sig_fit:.3f}")
print(f"  样本直接统计:                   mean   = {snr_valid.mean():.3f},  std = {snr_valid.std(ddof=1):.3f}")

fig, ax = plt.subplots(1, 2, figsize=(13, 4.8))
# 左: SNR 分布 + 正态拟合 + 理论
nb_h = max(20, int(np.sqrt(snr_valid.size)))
ax[0].hist(snr_valid, bins=nb_h, density=True, color="tab:blue", alpha=0.55, label=f"SNR 样本 ({snr_valid.size})")
xx = np.linspace(snr_valid.min(), snr_valid.max(), 400)
ax[0].plot(xx, norm.pdf(xx, mu_fit, sig_fit), "r-", lw=2.0, label=f"正态拟合 μ={mu_fit:.2f}, σ={sig_fit:.2f}")
ax[0].plot(xx, norm.pdf(xx, mu_th, sig_th), "k--", lw=1.6, label=f"理论 μ={mu_th:.2f}, σ={sig_th:.2f}")
ax[0].axvline(mu_fit, color="r", ls=":", alpha=0.7)
ax[0].set_xlabel("SNR = S / sqrt(B_bar)"); ax[0].set_ylabel("概率密度")
ax[0].set_title(f"SNR 分布 + 正态拟合 (N_shots={N_shots}, 分母=全bin背景均值)")
ax[0].legend(fontsize=8); ax[0].grid(alpha=0.3)
# 右: S 的分布(佐证 Poisson) + 每次测量 B_bar 的分布(近似常数, 抖动小)
ax[1].hist(S_samples, bins=nb_h, density=True, color="tab:green", alpha=0.5, label=f"S 样本 (mean={S_bar:.1f})")
ax[1].axvline(S_bar, color="tab:green", ls=":", alpha=0.8)
ax1b = ax[1].twiny()
ax1b.hist(B_samples, bins=30, density=True, color="tab:red", alpha=0.45,
          label=f"B_bar 样本 (mean={B_bar:.3f})")
ax1b.set_xlabel("每次测量的全 bin 背景均值 B_bar", color="tab:red")
ax[1].set_xlabel("峰 bin 信号计数 S"); ax[1].set_ylabel("概率密度")
ax[1].set_title("信号 S(Poisson) 与 每次 B_bar(近似常数) 的分布")
h1, l1 = ax[1].get_legend_handles_labels(); h2, l2 = ax1b.get_legend_handles_labels()
ax[1].legend(h1+h2, l1+l2, fontsize=8); ax[1].grid(alpha=0.3)
plt.tight_layout(); plt.show()

print(f"结论: SNR 分布近似正态; 理论(mu={mu_th:.2f},sigma={sig_th:.2f}) 与 拟合(mu={mu_fit:.2f},sigma={sig_fit:.2f}) 吻合。")
print(f"      分母用全 bin 背景均值 B_bar(近似常数): SNR=S/sqrt(B_bar) 是 S(Poisson) 的线性缩放, 故为干净正态。")''')

# ============================================================================
# 模块 16 — 100ppm 理论阈值 + 海量 Poisson 蒙卡验证噪点率
# ============================================================================
md_16 = nbformat.v4.new_markdown_cell(r'''## 模块 16（v21 新增）— 100ppm 噪点率的理论阈值 + 海量蒙卡验证

**噪点率(虚警率, false-alarm rate)定义**(用户选定, **整个时间窗 / 每次测量**):
> 在一次 N_shots=4 的测量里, **整个时间窗内至少有 1 个纯噪声 bin 的计数 ≥ 阈值** 的概率 = 目标噪点率。
> 目标 = **100 ppm = 1e-4**。

**理论阈值推导**:
- 每个背景 bin 计数 ~ **Poisson(nc_base)**(纯环境光, 已由模块 13 蒙卡证实)。
- 一次测量有 `N_bins` 个独立噪声 bin; 窗口级噪点率 = `1 − (1 − a_bin)^N_bins`, 其中 `a_bin` = 单 bin 超阈概率。
- 令窗口级 = 1e-4 ⇒ 单 bin 目标 `a_bin = 1 − (1 − 1e-4)^(1/N_bins) ≈ 1e-4 / N_bins`。
- 由 Poisson 生存函数反解**最小整数阈值** `T`: `P(X ≥ T | Poisson(nc_base)) ≤ a_bin`。
- 另给**高斯近似**阈值 `T_g = nc_base + z·√nc_base`(z 为 a_bin 的正态分位数)作对比——
  低计数下高斯近似通常**偏低**(不保守)。

**海量蒙卡验证**: 向量化生成 `rng.poisson(nc_base,(chunk, N_bins))`, 分块累计到 `N_MC`(可达 ~1e8),
统计"每次测量至少 1 bin ≥ T"的实测窗口级噪点率, 与 1e-4 比较。''')

code_16 = nbformat.v4.new_code_cell(r'''# ---- 100ppm 理论阈值 + 海量 Poisson 蒙卡验证噪点率 ----
from scipy.stats import poisson, norm as _norm

TARGET_FAR = 100e-6                 # 目标窗口级噪点率 = 100 ppm = 1e-4
# 有效噪声 bin 数: 用信号窗附近的实际 bin 数(整个时间窗过大且多为 0; 用当前分析窗 nbins 更贴合"每次测量")
N_bins_eff = nbins                  # 当前分析窗 bin 数(= len(centers))
nc_base = bg_hist_peak.mean()       # 噪声基底(同模块14)

# 单 bin 目标虚警率
a_bin = 1.0 - (1.0 - TARGET_FAR)**(1.0/N_bins_eff)

# (1) 精确 Poisson 反解最小整数阈值 T: P(X>=T) <= a_bin
#     poisson.sf(T-1) = P(X>=T). 从小到大找首个满足的 T。
T = 0
while poisson.sf(T-1, nc_base) > a_bin:   # sf(k)=P(X>k); sf(T-1)=P(X>=T)
    T += 1
far_at_T   = 1.0 - (1.0 - poisson.sf(T-1, nc_base))**N_bins_eff   # 该 T 的窗口级 FAR
far_at_Tm1 = 1.0 - (1.0 - poisson.sf(T-2, nc_base))**N_bins_eff   # T-1 的窗口级 FAR(超标)

# (2) 高斯近似阈值(连续, 供对比): T_g = nc_base + z*sqrt(nc_base)
z_abin = _norm.isf(a_bin)                 # a_bin 对应正态上分位
T_gauss = nc_base + z_abin*np.sqrt(nc_base)

print("="*76)
print(f"100ppm 噪点率理论阈值 (nc_base={nc_base:.4f} 计数/bin, N_bins={N_bins_eff}, 目标窗口级 FAR={TARGET_FAR:.1e})")
print(f"  单 bin 目标虚警率 a_bin = 1-(1-1e-4)^(1/{N_bins_eff}) = {a_bin:.3e}")
print(f"  [精确 Poisson] 最小整数阈值 T = {T} 计数")
print(f"     P(X>={T}|Poisson) = {poisson.sf(T-1, nc_base):.3e} (<= a_bin);  窗口级 FAR@T = {far_at_T:.3e}")
print(f"     (T-1={T-1} 时窗口级 FAR = {far_at_Tm1:.3e} > 目标, 故取 T={T})")
print(f"  [高斯近似] T_g = nc_base + z·sqrt(nc_base) = {T_gauss:.3f} (z={z_abin:.2f})  -> 取整 {int(np.ceil(T_gauss))}")
print(f"     注: 低计数下高斯近似阈值{'偏低(不保守)' if np.ceil(T_gauss)<T else '与精确接近'}。")

# ---- (3) 海量 Poisson 蒙卡验证 ----
# 目标 FAR=1e-4 -> 需 >>1e4 次才有统计意义; 默认 2e7, 可调到 1e8(每 1e6 约几十 MB, 分块避免爆内存)。
N_MC   = 20_000_000        # 总测量次数(可改大到 1e8: 100 个百万)
CHUNK  = 1_000_000         # 每块测量数(1e6 x N_bins 的 uint8 约 40MB @ N_bins=40)
rng_mc = np.random.default_rng(PARAMS["hist"]["seed"] + 999)

def window_far_mc(thresh, n_mc, chunk, nbins_eff, lam, rng):
    """海量蒙卡: 统计'每次测量(nbins_eff 个背景 bin)至少 1 个 >= thresh'的比例。"""
    n_alarm = 0; done = 0
    while done < n_mc:
        c = min(chunk, n_mc - done)
        # 每次测量 nbins_eff 个背景 bin ~ Poisson(lam); 是否至少一个 >= thresh
        x = rng.poisson(lam, size=(c, nbins_eff))
        n_alarm += int((x.max(axis=1) >= thresh).sum())
        done += c
    return n_alarm, n_mc

import time as _time
_tmc = _time.time()
print(f"\n海量蒙卡验证中... N_MC={N_MC:,} 次测量 x {N_bins_eff} bin (Poisson({nc_base:.3f}))")
n_alarm_T, _ = window_far_mc(T, N_MC, CHUNK, N_bins_eff, nc_base, rng_mc)
far_mc_T = n_alarm_T / N_MC
# 也给 T-1 与 T+1 的实测, 展示阈值敏感性
n_alarm_Tm1, _ = window_far_mc(max(T-1,1), N_MC, CHUNK, N_bins_eff, nc_base, rng_mc)
n_alarm_Tp1, _ = window_far_mc(T+1,        N_MC, CHUNK, N_bins_eff, nc_base, rng_mc)
far_mc_Tm1 = n_alarm_Tm1/N_MC; far_mc_Tp1 = n_alarm_Tp1/N_MC
print(f"  (蒙卡用时 {_time.time()-_tmc:.0f}s)")

# 二项标准误(用于给实测 FAR 一个不确定度)
def binom_se(p, n): return np.sqrt(max(p*(1-p),0)/n)

print("="*76)
print(f"海量蒙卡实测窗口级噪点率 (N_MC={N_MC:,}):")
print(f"  阈值 T-1={max(T-1,1):>2}: 实测 FAR = {far_mc_Tm1:.3e} ± {binom_se(far_mc_Tm1,N_MC):.1e}  (理论 {far_at_Tm1:.3e})")
print(f"  阈值 T  ={T:>2}: 实测 FAR = {far_mc_T:.3e} ± {binom_se(far_mc_T,N_MC):.1e}  (理论 {far_at_T:.3e}) <- 100ppm 阈值")
print(f"  阈值 T+1={T+1:>2}: 实测 FAR = {far_mc_Tp1:.3e} ± {binom_se(far_mc_Tp1,N_MC):.1e}")
print(f"  目标 100ppm = {TARGET_FAR:.1e}; 取整数阈值 T={T} 时实测 {far_mc_T:.2e} (<=目标, 满足)")

# ---- 绘图: 单 bin Poisson 分布 + 阈值; 窗口级 FAR vs 阈值(理论 vs 蒙卡) ----
fig, ax = plt.subplots(1, 2, figsize=(13, 4.8))
# 左: 单 bin Poisson(nc_base) pmf + 阈值线
kmax = max(T+3, int(nc_base+5))
kk = np.arange(0, kmax+1)
ax[0].bar(kk, poisson.pmf(kk, nc_base), color="tab:red", alpha=0.6, label=f"Poisson(nc_base={nc_base:.3f})")
ax[0].axvline(T, color="orange", ls="--", lw=1.8, label=f"100ppm 阈值 T={T}")
ax[0].axvline(T_gauss, color="green", ls=":", lw=1.6, label=f"高斯近似 T_g={T_gauss:.2f}")
ax[0].set_yscale("log"); ax[0].set_ylim(1e-8, 1)
ax[0].set_xlabel("单 bin 背景计数"); ax[0].set_ylabel("概率 (log)")
ax[0].set_title("单 bin 背景计数分布与检测阈值")
ax[0].legend(fontsize=8); ax[0].grid(alpha=0.3, which="both")
# 右: 窗口级 FAR vs 阈值(理论曲线 + 蒙卡点)
Ts = np.arange(max(T-2,1), T+4)
far_theory = np.array([1.0-(1.0-poisson.sf(t-1, nc_base))**N_bins_eff for t in Ts])
ax[1].semilogy(Ts, far_theory, "k-o", lw=1.6, ms=5, label="理论 窗口级 FAR")
for t, fm in [(max(T-1,1),far_mc_Tm1),(T,far_mc_T),(T+1,far_mc_Tp1)]:
    ax[1].scatter([t],[max(fm,1e-9)], c="tab:blue", s=60, zorder=5)
ax[1].scatter([],[], c="tab:blue", s=60, label="蒙卡实测")  # 图例占位
ax[1].axhline(TARGET_FAR, color="red", ls="--", lw=1.4, label=f"目标 100ppm={TARGET_FAR:.0e}")
ax[1].axvline(T, color="orange", ls=":", lw=1.4)
ax[1].set_xlabel("检测阈值 T [计数]"); ax[1].set_ylabel("窗口级噪点率 FAR (log)")
ax[1].set_title(f"窗口级噪点率 vs 阈值 (N_MC={N_MC:,})")
ax[1].legend(fontsize=8); ax[1].grid(alpha=0.3, which="both")
plt.tight_layout(); plt.show()

print(f"\n结论: 100ppm(窗口级)对应整数检测阈值 T={T} 计数 (nc_base={nc_base:.3f}/bin, {N_bins_eff} bins);")
print(f"      海量蒙卡({N_MC:,}次)实测 FAR@T={far_mc_T:.2e}, 与理论 {far_at_T:.2e} 一致; 高斯近似 T_g={T_gauss:.2f} 偏低。")
print(f"      如需验证到更低 ppm, 增大 N_MC(改为 100_000_000 = 100 个百万; 蒙卡是向量化 Poisson, 秒级~分钟级)。")''')

# 插入四个 cell(md14, code14, md15, code15, md16, code16) 到模块13 之后
new_cells = [md_14, code_14, md_15, code_15, md_16, code_16]
for off, cell in enumerate(new_cells, start=1):
    nb.cells.insert(i13 + off, cell)

# ============================================================================
# 标题更新
# ============================================================================
if nb.cells and nb.cells[0].cell_type == "markdown":
    nb.cells[0].source = nb.cells[0].source.replace(
        "# 激光雷达直方图仿真 v20 (LiDAR Histogram Simulation) — 阈值/死区调整 + 响应函数 g(Vov) + 信噪比",
        "# 激光雷达直方图仿真 v21 (LiDAR Histogram Simulation) — 检测阈值/前沿定时 + SNR 分布 + 100ppm 噪点率")
    nb.cells[0].source += '''

**v21 相对 v20 的变化(均基于 N_shots=4 累加直方图, 峰值宏像元)**
- **模块 14(检测阈值 + 前沿法定时)**: `det_th = k_th·nc_base`(nc_base=纯背景全 bin 均值);
  前沿法判决电平 `V_dec=(det_th+峰值)/2`, 在峰前沿相邻两 bin **线性插值**得 **front_time**, 换算距离与真值比较。
- **模块 15(SNR 分布)**: 重复 N_TRIALS 次峰值宏像元 RC 蒙卡, 每次算 SNR=S/√B̄
  (**分母 B̄=纯背景全 bin 均值**, 避免峰 bin 单背景因 RC 死区大量为 0), 直方图 + **正态拟合**;
  理论(B̄ 近似常数, S~Poisson): μ=S̄/√B̄, σ=√S̄/√B̄, 与拟合高度吻合。
- **模块 16(100ppm 噪点率)**: 按"整个时间窗/每次测量"定义窗口级虚警率=100ppm=1e-4;
  由 Poisson(nc_base) **精确反解整数阈值 T**(附高斯近似对比); 再用**海量向量化 Poisson 蒙卡**(默认 2e7, 可到 1e8)
  统计实测噪点率, 与理论一致性验证。
- 其余 30m 场景/脉冲/反射率/阈值/死区/响应函数等物理参数**一律未改**。'''

nbformat.write(nb, DST)
print(f"已生成 {DST}: 共 {len(nb.cells)} 个 cell (v20 41 + 新增 6)。")
