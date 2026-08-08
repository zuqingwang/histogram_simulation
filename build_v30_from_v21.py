# -*- coding: utf-8 -*-
# =============================================================================
# build_v30_from_v21.py
#   由 lidar_histogram_sim_v21_clean.ipynb 派生出 lidar_histogram_sim_v30.ipynb
#   —— 不覆盖 v21_clean，只读它、另存为 v30。
#
# v30 与 v21 的本质差异（物理统计方式）：
#   v21：SPAD 记录每个光子的 timestamp（时间戳），再 histogram 成直方图。
#   v30：系统【不记录 timestamp】。激光发射后以【1ns 时钟】逐点采样，
#        每个采样时刻只判断"该 SPAD 输出电压是否越过 60% 阈值"——
#        过了记 1，没过记 0。单个光子引发一次雪崩 → 输出被拉高一个
#        【8ns 宽、高度=1】的过阈值窗（窗宽由 τ_RC 与 Vth 推出，见下）。
#
#   两大现象由同一机制解释：
#     · 小信号展宽：单雪崩把回波在时间轴上摊成 8ns 宽的 1，且雪崩触发
#       时刻随机 → 直方图被抹宽。
#     · 大信号饱和：二值下单个 SPAD 每个 1ns bin、每次 shot 至多贡献 1；
#       雪崩密集时 8ns 窗互相顺延、连成一片全 1 → 峰 bin 触顶。
#       宏像元 27 个 SPAD 求和、再累加 N_shots ⇒ 峰 bin 硬上限 = 27×N_shots。
#
#   延迟游标（统一向后 delay，用于亚 ns 测距）：
#     系统计时零点 = start；激光【实际发光时刻】相对 start 向后延迟
#       t_laser = (tx_trig_dly + tx_trig_tcode)·1ns + delta_dly·(1/12)ns
#     1ns 采样格【固定不动】（钉在 start 起的整 ns）；改变 t_laser 让回波
#     相对采样格平移，扫描 delta_dly=0..11 即得 1/12 ns(≈83.3ps) 的测距分辨。
#
#   本脚本做法（最小侵入、复用 v21 物理链路）：
#     · 模块 0–9（PARAMS、激光/光学/信道/目标、光链路、RC 恢复机制展示）原样保留；
#       其中模块 8b 的 Vov(t) RC 恢复曲线正是 8ns 过阈窗的物理来源，保留作教学。
#     · 新增【模块 9b】：定义 v30 二值采样引擎 spad_binary_trace() 与过阈窗宽
#       T_OVER = −τ_RC·ln(1−Vth_frac) ≈ 8ns；给 PARAMS 补充 timing 三参数；
#       并画"单 SPAD 二值波形"演示（弱信号展宽 vs 强信号饱和）。
#     · 模块 11 / 13 / 15：改用二值引擎产出 macro_hist / bg_hist_peak（契约不变）。
#     · 模块 12：保留热图(图A)+直方图(图B)，删去 timestamp eventplot(图C)。
#     · 新增【模块 12b】：delta_dly 游标扫描 → 亚 ns 测距演示。
#     · 模块 14 / 16：仅依赖契约变量(macro_hist/bg_hist_peak/pk_bin/…)，无缝复用。
#     · 最前面插入 v30 说明标题 cell。
# =============================================================================
import nbformat

SRC = "lidar_histogram_sim_v21_clean.ipynb"
DST = "lidar_histogram_sim_v30.ipynb"

nb = nbformat.read(SRC, as_version=4)

# 清空所有输出与执行计数，得到干净的可分发 notebook
for c in nb.cells:
    if c.cell_type == "code":
        c.outputs = []
        c.execution_count = None
    if "execution" in c.get("metadata", {}):
        del c.metadata["execution"]


def find(pred):
    """返回首个满足 pred 的 (index, cell)，找不到抛错。"""
    for i, c in enumerate(nb.cells):
        if pred(c):
            return i, c
    raise RuntimeError("cell not found")


def find_code_startswith(marker):
    """按 code cell 源码是否包含 marker 定位（marker 取该 cell 的独特首行片段）。"""
    return find(lambda c: c.cell_type == "code" and marker in c.source)


# =============================================================================
# 1) 新增【模块 9b】二值采样引擎（定义 + 单 SPAD 波形演示）
#    插入位置：模块 9（CELL 源码含 'thr = fpix0.max() * 0.01'）之后。
# =============================================================================
md_9b = nbformat.v4.new_markdown_cell(r'''## 模块 9b（v30 核心）— 二值采样引擎：8ns 过阈值窗 + 1ns 时钟 0/1 采样

**这是 v30 与 v21 的分水岭。** v21 记录每个光子的 timestamp 再做直方图；v30 **不记录 timestamp**，而是：

1. 激光发射后，系统以 **1ns 时钟**在每个整 ns 采样一次；
2. 每个采样点只判断该 SPAD 输出电压**是否越过 60% 阈值**（`Vth_frac`）——过则记 **1**，否则记 **0**；
3. 一次雪崩（单光子被探测）把输出拉高成一个 **8ns 宽、高度=1** 的"过阈值窗"。窗宽由 RC 恢复推出：

$$t_\mathrm{over} = -\,\tau_\mathrm{RC}\cdot\ln\!\left(1-V_\mathrm{th,frac}\right)
= -\,8.7315\,\mathrm{ns}\times\ln(1-0.60)\approx 8.0\,\mathrm{ns}$$

> 即：雪崩后 Vov 跌到 0，按 τ_RC 指数充回；充回到 60% 阈值前，读出端一直判"过阈=1"，历时恰好 ≈8ns。**不新增独立参数**，直接由 τ_RC 与 Vth_frac 导出。

**两大现象，同一机制：**

```
弱信号(几个分立雪崩)          强信号(雪崩密集 → 8ns 窗顺延堆积)
 阈值─ ┌──8ns──┐   ┌──8ns──┐    阈值─ ┌────────连成一片──────────┐
   1 ─┤        ├───┤       ├      1 ─┤   (每个新雪崩把窗口再顺延)  ├
   0 ─┘        └───┘       └      0 ─┘                            └
     └ 单雪崩被摊成 8ns 宽 →展宽      └ 单 SPAD 每 bin 每 shot 至多 1 →饱和
```

**宏像元合并 & 饱和上限：** 每个 SPAD 对每个 1ns bin 独立给 0/1；宏像元 27 个 SPAD 相加（0–27），再累加 N_shots ⇒ **峰 bin 硬上限 = 27 × N_shots**。

**延迟游标（亚 ns 测距）：** 1ns 采样格固定钉在 start 起的整 ns；激光实际发光相对 start 向后延迟

$$t_\mathrm{laser} = (\texttt{tx\_trig\_dly}+\texttt{tx\_trig\_tcode})\cdot 1\,\mathrm{ns} + \texttt{delta\_dly}\cdot\tfrac{1}{12}\,\mathrm{ns}$$

改变 `t_laser` 让回波相对采样格平移；扫描 `delta_dly = 0..11` 即得 **1/12 ns ≈ 83.3 ps** 的测距分辨（见模块 12b）。默认三参数=0（`t_laser=0`），此时回波峰位与 v21 一致，仅统计方式变为二值。''')

code_9b = nbformat.v4.new_code_cell(r'''# ---- v30 二值采样引擎：定义 + 过阈窗宽 + timing 三参数 + 单 SPAD 波形演示 ----

# (1) 过阈值窗宽 T_OVER：由 RC 恢复与阈值推出，≈8ns（不新增独立参数）
#     雪崩后 Vov 从 0 按 τ_RC 指数充回，充到 Vth_frac 之前一直判"过阈=1"。
T_OVER = -TAU_RC * np.log(1.0 - VTH_FRAC)     # = -8.7315e-9 * ln(0.4) ≈ 8.0e-9 s
print(f"过阈值窗宽 T_OVER = -τ_RC·ln(1-Vth_frac) = {T_OVER*1e9:.3f} ns  "
      f"(τ_RC={TAU_RC*1e9:.4f}ns, Vth_frac={VTH_FRAC:.0%})")

# (2) timing 三参数：统一把"激光实际发光时刻"相对 start 向后延迟。
#     采样格(1ns 时钟)固定不动，靠 t_laser 让回波相对采样格平移 → 亚 ns 测距。
#     默认全 0 => t_laser=0 => 回波峰位与 v21 一致（仅统计方式变二值）。
#     注：此处仅【新增】timing 键，不改动 PARAMS 中任何既有物理参数值。
PARAMS.setdefault("timing", {
    "tx_trig_dly":   0,      # 单位 1ns（整数 ns 粗延迟之一）
    "tx_trig_tcode": 0,      # 单位 1ns（整数 ns 粗延迟之二）
    "delta_dly":     0,      # 单位 1/12 ns（≈83.3ps 精延迟游标；扫 0..11 测距）
})

def laser_delay(timing):
    """由 timing 三参数算激光相对 start 的向后延迟 t_laser [s]。"""
    return ((timing["tx_trig_dly"] + timing["tx_trig_tcode"]) * 1e-9
            + timing["delta_dly"] * (1e-9 / 12.0))

T_LASER = laser_delay(PARAMS["timing"])
print(f"激光延迟 t_laser = (tx_trig_dly+tx_trig_tcode)·1ns + delta_dly·(1/12)ns = {T_LASER*1e9:.4f} ns "
      f"(默认三参数={PARAMS['timing']})")


def spad_binary_trace(r_sig_fine, r_amb_ph, tf, centers, PDE_max, tau_rc, Vth_frac,
                      jitter_sigma, rng, t_over, t_laser=0.0,
                      resp_shape="linear", resp_k=3.0):
    """v30 单 SPAD、单次 shot 的【二值采样】。

    返回长度 = len(centers) 的 int8 数组：每个 1ns 采样点的 0/1（是否过阈值）。

    机制（与 v21 RC 引擎共用"雪崩触发"物理，但输出方式不同）：
      1) 光子按 Poisson(r·dt) 到达细网格 tf（信号 r_sig_fine + 背景 r_amb_ph）；
      2) 每个光子按 触发概率 = PDE_max·g(vov_frac) 判定是否引发雪崩；
         vov_frac 按 RC 恢复 1-exp(-(t-last)/τ) 计算（同 v21，恢复期灵敏度渐增）；
      3) 每次雪崩把输出在 [t_av, t_av + t_over] 内拉高为 1（可顺延/堆积）；
      4) 采样点 = centers + t_laser（激光延迟把回波相对采样格整体后移）；
         采样点落入任一拉高区间 → 该 bin 记 1。（v21 是把 timestamp 塞进 histogram）
    """
    dt = tf[1] - tf[0]
    mu = (r_sig_fine + r_amb_ph) * dt
    n_ph = rng.poisson(mu)                      # 各细网格点到达光子数 ~ Poisson
    nbn = len(centers)
    out = np.zeros(nbn, dtype=np.int8)
    if n_ph.sum() == 0:
        return out
    t_arr = np.repeat(tf, n_ph)                 # 展开成逐光子到达时刻
    u = rng.random(t_arr.size)
    # --- 逐光子判定雪崩（Vov 按 RC 恢复；触发即把 Vov 拉回 0 重新恢复）---
    last = -1e30
    inv_tau = 1.0 / tau_rc
    av = []                                     # 记录雪崩发生时刻
    for k in range(t_arr.size):
        t = t_arr[k]
        d = (t - last) * inv_tau
        vov_frac = 1.0 - np.exp(-d) if d < 700 else 1.0
        p_fire = PDE_max * spad_response_g(vov_frac, resp_shape, resp_k)
        if u[k] < p_fire:
            av.append(t)
            last = t
    if not av:
        return out
    av = np.asarray(av)
    if jitter_sigma > 0:                        # IRF 抖动作用在雪崩(→过阈窗)时刻
        av = av + rng.normal(0.0, jitter_sigma, av.size)
    # --- 二值采样：采样点固定钉在 start 起整 ns(centers 不动)；
    #     激光延迟 t_laser 把回波(雪崩→8ns 过阈窗)【整体向后移】：
    #     窗口 = [t_av + t_laser, t_av + t_laser + t_over]，采样点落入即置 1。
    #     (t_laser=0 时窗口=[t_av, t_av+t_over]，与 v21 峰位一致；
    #      t_laser 增大 → 回波右移 → 质心随之增大，斜率+1，见模块 12b。)
    for tt in av:
        lo = tt + t_laser
        out[(centers >= lo) & (centers < lo + t_over)] = 1
    return out


def over_waveform(av_times, t_grid, t_over):
    """辅助(仅演示用)：在细网格 t_grid 上，把每次雪崩的 [t_av, t_av+t_over] 置 1，
    得到连续的"过阈值方波"，用于画 8ns 窗的堆积效果。"""
    w = np.zeros_like(t_grid)
    for tt in np.atleast_1d(av_times):
        w[(t_grid >= tt) & (t_grid < tt + t_over)] = 1.0
    return w


# ---- 单 SPAD 二值波形演示：弱信号(展宽) vs 强信号(饱和) ----
# 说明：这里为可视化对比，人为设置"弱/强"两档信号率(强档= r_sig 放大若干倍)，
#       仅用于演示二值机制，不改动 PARAMS 中的物理场景参数。
rng_demo = np.random.default_rng(PARAMS["hist"]["seed"] + 96)
tgrid = np.arange(t_lo, t_hi, dt_fine)                 # 细网格(画连续方波用)
tgrid_ns = tgrid * 1e9

# 直接指定几处"雪崩时刻"以清晰展示 8ns 窗（弱：分立；强：密集堆积）
av_weak   = np.array([t0_ns - 2.0, t0_ns + 12.0]) * 1e-9          # 两次分立雪崩
av_strong = (t0_ns + np.array([-3, -1, 1, 2, 4, 6, 9, 11.0])) * 1e-9  # 密集雪崩

w_weak   = over_waveform(av_weak,   tgrid, T_OVER)
w_strong = over_waveform(av_strong, tgrid, T_OVER)

# 1ns 采样点上的 0/1（用 centers 作采样时刻）
def sample_binary(av_times, centers, t_over, t_laser=0.0):
    # 采样点固定(centers)，激光延迟把回波窗 [av+t_laser, av+t_laser+t_over] 右移
    b = np.zeros(len(centers), dtype=int)
    for tt in np.atleast_1d(av_times):
        lo = tt + t_laser
        b[(centers >= lo) & (centers < lo + t_over)] = 1
    return b

b_weak   = sample_binary(av_weak,   centers, T_OVER)
b_strong = sample_binary(av_strong, centers, T_OVER)

fig, ax = plt.subplots(1, 2, figsize=(13, 4.4), sharey=True)
for a, (w, b, av, ttl) in zip(
        ax,
        [(w_weak,   b_weak,   av_weak,   "弱信号：单雪崩→8ns 窗(展宽)"),
         (w_strong, b_strong, av_strong, "强信号：8ns 窗顺延堆积→连成一片(饱和)")]):
    a.plot(tgrid_ns, w, color="tab:blue", lw=1.4, label="过阈值(连续) over(t)")
    a.axhline(1.0, color="0.7", ls="-", lw=0.6)
    # 1ns 采样点 0/1
    a.plot(tc_ns, b, "o", color="tab:green", ms=5, label="1ns 采样 0/1")
    a.vlines(tc_ns, 0, b, color="tab:green", lw=1.0, alpha=0.5)
    # 雪崩时刻标记
    for j, tt in enumerate(av * 1e9):
        a.axvline(tt, color="tab:red", ls=":", lw=1.0,
                  label="雪崩时刻" if j == 0 else None)
    a.axvline(t0_ns, color="k", ls="--", alpha=0.5, lw=1.0, label=f"ToF {t0_ns:.1f} ns")
    a.set_xlabel("时间 t [ns]"); a.set_title(ttl, fontsize=10)
    a.set_ylim(-0.08, 1.15); a.grid(alpha=0.3); a.legend(fontsize=8, loc="upper right")
ax[0].set_ylabel("过阈值电平 / 采样值")
plt.suptitle(f"v30 二值采样引擎：每雪崩→{T_OVER*1e9:.1f}ns 过阈窗，1ns 时钟采 0/1", fontsize=11)
plt.tight_layout(); plt.show()

print("="*76)
print("v30 二值采样引擎已就绪：")
print(f"  · 过阈窗宽 T_OVER = {T_OVER*1e9:.2f} ns（由 τ_RC 与 Vth 导出，非独立参数）")
print(f"  · 弱信号：{av_weak.size} 次分立雪崩 → 采到 {int(b_weak.sum())} 个 '1' bin（每雪崩摊成~8ns 宽）")
print(f"  · 强信号：{av_strong.size} 次密集雪崩 → 采到 {int(b_strong.sum())} 个 '1' bin（窗堆积、趋于连续=饱和前兆）")
print(f"  · 单 SPAD 每 bin 每 shot 至多 1 ⇒ 宏像元(27 SPAD)×N_shots 后峰 bin 硬上限 = {27*PARAMS['hist']['N_shots']}")''')

i9, _ = find_code_startswith("thr = fpix0.max() * 0.01")
nb.cells.insert(i9 + 1, md_9b)
nb.cells.insert(i9 + 2, code_9b)


# =============================================================================
# 2) 重写【模块 11】：逐-SPAD 二值采样 → macro_hist（契约不变）
# =============================================================================
code_11_new = r'''nbins = len(centers)
base_rate = signal_photon_rate_fine(echo0, 1.0, tf)   # 单位收集比例(f=1)的信号率, 后续乘 f_ij
N_shots = PARAMS["hist"]["N_shots"]
macro_fvals = [fpix0[:, m*By_m:(m+1)*By_m].ravel() for m in range(n_macro)]

# 激光延迟(默认 0)：把回波相对 1ns 采样格整体后移 t_laser
T_LASER = laser_delay(PARAMS["timing"])

# ---- v30 二值模型: 逐-SPAD 二值采样, 27 SPAD 求和, 累加 N_shots ----
#   每个 SPAD 每 bin 每 shot ∈ {0,1}; 宏像元求和后 ∈ [0,27]; 再 ×N_shots 累加。
rng = np.random.default_rng(PARAMS["hist"]["seed"] + 11)
macro_hist = np.zeros((n_macro, nbins))
for _shot in range(N_shots):
    for m in range(n_macro):
        acc = np.zeros(nbins, dtype=np.int32)          # 本 shot 本宏像元的 27-SPAD 求和
        for fij in macro_fvals[m]:
            acc += spad_binary_trace(base_rate*fij, r_amb_ph, tf, centers, PDE, TAU_RC,
                                     VTH_FRAC, jit, rng, T_OVER, T_LASER, RESP_SHAPE, RESP_K)
        macro_hist[m] += acc

macro_cap = n_pix_macro * N_shots                       # 二值硬上限(峰 bin 触顶值)

# ---- 泊松期望计数(未二值封顶)作参考虚线: 若不受"每bin至多1"限制会有多少 ----
#   与 v21 相同的解析上限, 展示强信号处二值饱和相对泊松期望的"削顶"差异。
irf_k = gaussian_kernel(jit, dt_fine)
bin_idx = np.clip(((tf - t_lo)/bin_width).astype(int), 0, nbins-1)
macro_lam = np.zeros((n_macro, nbins))
for m in range(n_macro):
    r_det = (base_rate*macro_fsum[m] + n_pix_macro*r_amb_ph) * PDE
    r_det = np.convolve(r_det, irf_k, mode="same") * dt_fine
    macro_lam[m] = N_shots * np.bincount(bin_idx, weights=r_det*dt_fine, minlength=nbins)

tot_bin, tot_la = macro_hist.sum(), macro_lam.sum()
pk_peakbin = macro_hist[m_peak].max()
n_saturated = int((macro_hist[m_peak] >= macro_cap - 1e-9).sum())
print("="*76)
print(f"v30 二值采样: {n_macro} 宏 × {n_pix_macro} SPAD × {N_shots} shots "
      f"(过阈窗 T_OVER={T_OVER*1e9:.1f}ns, Vth={VTH_FRAC*100:.0f}%)")
print(f"  峰 bin 二值硬上限 macro_cap = 27×N_shots = {macro_cap}")
print(f"  说明: 二值总计数(bin 级 0/1 之和)与'泊松光子期望'量纲不同, 不宜直接比'保留率';")
print(f"        单光子被摊成 {T_OVER*1e9:.0f}ns 宽的多个 '1' bin, 故二值总和可高于光子期望。")
print(f"  全体二值总计数 = {tot_bin:.0f}  (泊松光子期望参考 {tot_la:.0f}, 仅作虚线对照)")
print(f"峰值宏像元 m={m_peak}: 峰 bin={pk_peakbin:.0f}/{macro_cap} "
      f"({'已触顶饱和' if pk_peakbin>=macro_cap-1e-9 else '未触顶'}), "
      f"触顶 bin 数={n_saturated}, 总计数={macro_hist[m_peak].sum():.0f}; "
      f"边缘 m=0 总计数={macro_hist[0].sum():.0f}")'''

i11, c11 = find_code_startswith("base_rate = signal_photon_rate_fine(echo0, 1.0, tf)")
c11.source = code_11_new


# =============================================================================
# 3) 重写【模块 12】：热图(图A) + 直方图(图B, 含二值上限线)，删去 timestamp eventplot(图C)
# =============================================================================
code_12_new = r'''# ---- 图 A: 宏像元-时间 二值计数热图 ----
fig, ax = plt.subplots(figsize=(11, 5))
im = ax.imshow(macro_hist, origin="upper", aspect="auto", cmap="inferno",
               extent=[tc_ns[0]-0.5, tc_ns[-1]+0.5, n_macro-0.5, -0.5],
               vmin=0, vmax=macro_cap)
ax.axvline(t0_ns, color="cyan", ls=":", lw=1.2, label=f"真实 ToF {t0_ns:.1f} ns")
ax.axhline(m_peak, color="lime", ls=":", lw=1.0, alpha=0.7, label=f"峰值宏像元 m={m_peak}")
ax.set_xlabel("时间 t [ns]"); ax.set_ylabel("宏像元序号 m (0=顶部, 沿长边 y)")
ax.set_title(f"宏像元(9×3)二值直方图热图 (v30: 8ns 过阈窗, N_shots={N_shots}, 上限={macro_cap})")
ax.legend(fontsize=9, loc="upper right")
plt.colorbar(im, ax=ax, label=f"计数 / (宏像元, 1ns bin)  [0..{macro_cap}]")
plt.tight_layout(); plt.show()

# ---- 图 B: 代表性宏像元直方图 (点线=二值实测; 虚线=泊松期望未封顶; 水平线=二值硬上限) ----
fig, ax = plt.subplots(figsize=(11, 4.6))
reps_m = sorted(set([0, max(0, m_peak-6), m_peak, min(n_macro-1, m_peak+6)]))
colors = ["tab:gray", "tab:green", "tab:red", "tab:orange", "tab:purple"]
for m, c in zip(reps_m, colors):
    ax.plot(tc_ns, macro_hist[m], color=c, lw=1.4, marker="o", ms=4,
            label=f"m={m} 二值 (Σf={macro_fsum[m]:.3f})")          # 二值实测: 标记+连线
    ax.plot(tc_ns, macro_lam[m], color=c, lw=1.1, ls="--", alpha=0.6)  # 泊松期望(未封顶)参考
ax.axhline(macro_cap, color="k", ls="-.", lw=1.2, alpha=0.8,
           label=f"二值硬上限 27×N_shots={macro_cap}")             # 饱和封顶线
ax.axvline(t0_ns, color="k", ls=":", alpha=0.6, label=f"ToF {t0_ns:.1f} ns")
ax.set_xlabel("时间 t [ns]"); ax.set_ylabel(f"计数 / 1ns bin ({N_shots} shots)")
ax.set_title("代表性宏像元直方图 (点线=二值实测, 虚线=泊松期望未封顶, 点划线=饱和上限)")
ax.legend(fontsize=8); ax.grid(alpha=0.3)
plt.tight_layout(); plt.show()

# 说明: 强信号宏像元(如 m_peak)峰区二值实测被"削顶"到 27×N_shots, 而泊松期望虚线会更高
#       —— 这正是 v30 二值模型的【饱和】; 弱信号宏像元二者接近(未触顶)。
#       v21 的"27 SPAD 单次 timestamp eventplot(图C)"已按 v30 设计删去。
print(f"图A/B 完成: 峰值宏像元 m={m_peak} 峰 bin={macro_hist[m_peak].max():.0f}/{macro_cap} "
      f"({'饱和' if macro_hist[m_peak].max()>=macro_cap-1e-9 else '未饱和'}); "
      f"泊松期望峰 bin≈{macro_lam[m_peak].max():.0f}(未封顶参考)。")'''

i12, c12 = find_code_startswith("# ---- 图 A: 热图 ----")
c12.source = code_12_new


# =============================================================================
# 4) 新增【模块 12b】：delta_dly 游标扫描 → 亚 ns 测距演示
#    插入位置：模块 12 之后。
# =============================================================================
md_12b = nbformat.v4.new_markdown_cell(r'''## 模块 12b（v30 新增）— 延迟游标 `delta_dly` 扫描：亚 ns（1/12 ns）测距

1ns 采样格是**固定**的，本身只能给出 1ns（≈15cm）的量化。v30 靠**移动激光发光时刻**突破它：

$$t_\mathrm{laser}=(\texttt{tx\_trig\_dly}+\texttt{tx\_trig\_tcode})\cdot 1\,\mathrm{ns}+\texttt{delta\_dly}\cdot\tfrac{1}{12}\,\mathrm{ns}$$

固定粗延迟，令 `delta_dly = 0,1,…,11`，每档把回波（8ns 过阈窗）**整体向后移** 1/12 ns ≈ 83.3 ps（采样格不动）。对每档重复测量取平均直方图，用**前沿半高插值时刻**（leading-edge，连续量）定位回波；该时刻随 `t_laser` 严格**线性移动、斜率≈1** ⇒ 等效获得 **1/12 ns** 的测距分辨（12 档合成 ≈1.25cm 级）。

> 为何不用质心？8ns 过阈窗是"平顶"，质心对亚 ns 位移不灵敏，且受 1ns 量化台阶与固定统计窗边界牵引，斜率会系统性偏低。前沿半高插值只依赖上升沿、是连续量，对整体平移天然线性——正是模块 14 前沿定时法的思想。图中同时画出两者作对比。''')

code_12b = nbformat.v4.new_code_cell(r'''# ---- delta_dly 游标扫描: 亚 ns 测距演示 (峰值宏像元, 二值累加) ----
#   固定 tx_trig_dly / tx_trig_tcode, 只扫 delta_dly=0..11; 每档 t_laser 递增 1/12 ns。
#   每档重复 N_REP_DLY 次测量取平均直方图(压噪), 再用【前沿半高插值】定位回波时刻,
#   同时给出【峰区质心】作对比。前沿法对整体平移严格线性(斜率≈1)。
base_timing = dict(PARAMS["timing"])                 # 备份, 演示后原样还原(不改用户设置)
fvals_pk = macro_fvals[m_peak]
delta_list = list(range(12))
N_REP_DLY = 200                                      # 每档重复次数(平均直方图压噪; 越大越平滑)
front_ns = np.empty(len(delta_list))                 # 各档前沿半高时刻 [ns]
cent_ns  = np.empty(len(delta_list))                 # 各档峰区质心   [ns] (对比用)
tl_ns    = np.empty(len(delta_list))                 # 各档 t_laser    [ns]

# 分析窗: 以默认峰 bin 为中心 ±6 bin(质心用); 前沿法在整段 centers 上找上升沿
pk0 = int(np.argmax(macro_hist[m_peak]))
lo_w, hi_w = max(0, pk0-6), min(nbins, pk0+7)

def _run_avg_hist(tl, seed0):
    """峰值宏像元: 重复 N_REP_DLY 次(每次 27 SPAD × N_shots 二值累加)取平均直方图。"""
    h = np.zeros(nbins)
    for rep in range(N_REP_DLY):
        rng_d = np.random.default_rng(seed0 + rep)
        for _shot in range(N_shots):
            acc = np.zeros(nbins, dtype=np.int32)
            for fij in fvals_pk:
                acc += spad_binary_trace(base_rate*fij, r_amb_ph, tf, centers, PDE, TAU_RC,
                                         VTH_FRAC, jit, rng_d, T_OVER, tl, RESP_SHAPE, RESP_K)
            h += acc
    return h / N_REP_DLY                              # 每次测量的平均直方图

def _leading_edge_time(h):
    """前沿半高插值: 上升沿上首次跨越 (baseline+peak)/2 的线性插值时刻 [ns]。"""
    base = h[:max(pk0-6,1)].mean() if pk0 > 6 else h.min()   # 峰左侧底噪基线
    pk_i = int(np.argmax(h)); v_half = 0.5*(base + h[pk_i])
    for i in range(1, pk_i+1):
        if h[i-1] < v_half <= h[i]:
            frac = (v_half - h[i-1]) / (h[i] - h[i-1]) if h[i] > h[i-1] else 0.0
            return tc_ns[i-1] + frac*(tc_ns[i]-tc_ns[i-1])
    return np.nan

for q, dd in enumerate(delta_list):
    tm = dict(base_timing); tm["delta_dly"] = dd
    tl = laser_delay(tm); tl_ns[q] = tl*1e9
    h = _run_avg_hist(tl, PARAMS["hist"]["seed"] + 500 + dd*1000)
    front_ns[q] = _leading_edge_time(h)
    ww = h[lo_w:hi_w]; xx = tc_ns[lo_w:hi_w]
    cent_ns[q] = (xx*ww).sum()/ww.sum() if ww.sum() > 0 else np.nan

# 线性拟合(理论斜率应≈1: 激光每后移 Δ, 回波定位时刻同步后移 Δ)
okf = np.isfinite(front_ns); okc = np.isfinite(cent_ns)
slope_f, icpt_f = np.polyfit(tl_ns[okf], front_ns[okf], 1)
slope_c, icpt_c = np.polyfit(tl_ns[okc], cent_ns[okc], 1)

fig, ax = plt.subplots(1, 2, figsize=(13, 4.6))
# 左: 若干档的峰区平均直方图, 看回波随 delta_dly 平移
for dd in [0, 3, 6, 9, 11]:
    tm = dict(base_timing); tm["delta_dly"] = dd; tl = laser_delay(tm)
    h = _run_avg_hist(tl, PARAMS["hist"]["seed"] + 500 + dd*1000)
    ax[0].plot(tc_ns, h, marker="o", ms=3, lw=1.2, label=f"delta_dly={dd} (t_laser={tl*1e9:.2f}ns)")
ax[0].axvline(t0_ns, color="k", ls=":", alpha=0.6, label=f"ToF {t0_ns:.1f} ns")
ax[0].set_xlim(tc_ns[max(lo_w-1,0)], tc_ns[min(hi_w,nbins-1)])
ax[0].set_xlabel("时间 t [ns]"); ax[0].set_ylabel(f"平均二值计数 / bin ({N_shots} shots)")
ax[0].set_title(f"峰区平均直方图随 delta_dly 平移 (每档 1/12 ns, 平均 {N_REP_DLY} 次)")
ax[0].legend(fontsize=8); ax[0].grid(alpha=0.3)
# 右: 前沿半高时刻 & 质心 vs t_laser (前沿法斜率≈1)
ax[1].plot(tl_ns, front_ns, "o-", color="tab:blue", ms=5, label="前沿半高时刻(连续)")
ax[1].plot(tl_ns, slope_f*tl_ns+icpt_f, "b--", lw=1.0, label=f"前沿拟合 斜率={slope_f:.3f}")
ax[1].plot(tl_ns, cent_ns, "s-", color="tab:purple", ms=4, alpha=0.7, label="峰区质心(对比)")
ax[1].plot(tl_ns, slope_c*tl_ns+icpt_c, color="tab:purple", ls=":", lw=1.0,
           label=f"质心拟合 斜率={slope_c:.3f}")
ax[1].set_xlabel("激光延迟 t_laser [ns] (delta_dly·1/12 ns)")
ax[1].set_ylabel("回波定位时刻 [ns]")
ax[1].set_title("回波时刻随游标线性移动 → 亚 ns(1/12 ns) 分辨")
ax[1].legend(fontsize=8); ax[1].grid(alpha=0.3)
plt.tight_layout(); plt.show()

print("="*76)
print(f"delta_dly 扫描(0..11, 步进 1/12 ns={1000/12:.1f}ps, 每档平均 {N_REP_DLY} 次):")
print(f"  前沿半高法: 拟合斜率 = {slope_f:.3f} (理论≈1: 激光后移 Δ → 回波前沿同步后移 Δ)")
print(f"  峰区质心法: 拟合斜率 = {slope_c:.3f} (受 8ns 平顶/1ns 量化牵引, 系统性偏低, 仅作对比)")
print(f"  ⇒ 1ns 固定采样格 + 12 档游标 = 等效 1/12 ns 时间分辨 ≈ {C_LIGHT*(1e-9/12)/2*100:.2f} cm 测距步进")
# 还原用户原始 timing 设置(本演示不改动 PARAMS)
PARAMS["timing"] = base_timing''')

i12b, _ = find_code_startswith("# ---- 图 A: 宏像元-时间 二值计数热图 ----")
nb.cells.insert(i12b + 1, md_12b)
nb.cells.insert(i12b + 2, code_12b)


# =============================================================================
# 5) 重写【模块 13】背景蒙卡：改用二值引擎产 bg_hist_peak（其余 SNR 逻辑/绘图不变）
# =============================================================================
code_13_new = r'''# ---- 信噪比 SNR = S / sqrt(B), 在峰值宏像元 m_peak 的信号峰 bin (v30 二值模型) ----
# 复用模块11 的 macro_hist(信号+背景, 二值采样 N_shots)。
# 另跑"纯环境光"二值采样(信号率=0)估计每 bin 背景计数 B。

# 1) 纯背景二值采样: 峰值宏像元 27 SPAD, 信号率置 0, 仅环境光, 累计 N_shots
zero_rate = np.zeros_like(base_rate)
T_LASER = laser_delay(PARAMS["timing"])
rng_bg = np.random.default_rng(PARAMS["hist"]["seed"] + 313)
bg_hist_peak = np.zeros(nbins)
for _shot in range(N_shots):
    acc = np.zeros(nbins, dtype=np.int32)
    for fij in macro_fvals[m_peak]:          # fij 不影响背景(背景与 f_pix 无关), 但保持 27 次一致
        acc += spad_binary_trace(zero_rate, r_amb_ph, tf, centers, PDE, TAU_RC,
                                 VTH_FRAC, jit, rng_bg, T_OVER, T_LASER, RESP_SHAPE, RESP_K)
    bg_hist_peak += acc

# 2) 信号峰 bin: 取峰值宏像元 macro_hist 的最大 bin 作为信号峰位置
sig_hist_peak = macro_hist[m_peak]
pk_bin = int(np.argmax(sig_hist_peak))
tot_peak = sig_hist_peak[pk_bin]             # 峰 bin 总计数(信号+背景)
B_peak = bg_hist_peak[pk_bin]                # 峰 bin 背景计数(纯环境光二值)
S_peak = max(tot_peak - B_peak, 0.0)         # 扣背景后信号计数

# 背景平均(用全 bin 纯背景二值均值, 更稳健地代表底噪水平)
B_mean = bg_hist_peak.mean()
SNR_sqrtB   = S_peak / np.sqrt(B_peak) if B_peak > 0 else np.inf
SNR_sqrtB_m = S_peak / np.sqrt(B_mean) if B_mean > 0 else np.inf
SNR_sqrtSB  = S_peak / np.sqrt(S_peak + B_peak) if (S_peak + B_peak) > 0 else 0.0

print("="*76)
print(f"信噪比 SNR (峰值宏像元 m={m_peak}, 信号峰 bin @ {tc_ns[pk_bin]:.0f} ns, N_shots={N_shots}, v30 二值)")
print(f"  峰 bin 总计数(信号+背景) = {tot_peak:.1f}  (二值硬上限 {macro_cap})")
print(f"  背景 B(纯环境光二值, 该 bin) = {B_peak:.3f};  背景均值(全 bin) = {B_mean:.3f}")
print(f"  信号 S = 峰 bin 总计数 - B = {S_peak:.1f}")
print(f"  -> SNR = S/sqrt(B)      = {SNR_sqrtB:.2f}   [主定义, 背景受限]")
print(f"     SNR = S/sqrt(B_mean) = {SNR_sqrtB_m:.2f}   (用全 bin 背景均值作 B)")
print(f"     SNR = S/sqrt(S+B)    = {SNR_sqrtSB:.2f}   (含信号散粒噪声, 供参考)")

# 3) 各宏像元峰 bin 的 SNR 分布(用统一背景均值 B_mean 作分母, 快速给全局图景)
snr_per_macro = np.zeros(n_macro)
for m in range(n_macro):
    s_tot = macro_hist[m][pk_bin]
    s_sig = max(s_tot - B_mean, 0.0)
    snr_per_macro[m] = s_sig / np.sqrt(B_mean) if B_mean > 0 else 0.0

fig, ax = plt.subplots(1, 2, figsize=(13, 4.6))
# 左: 峰值宏像元直方图 + 标注 S / B / 峰 bin
ax[0].bar(tc_ns, sig_hist_peak, width=bin_width*1e9, align="center", color="tab:green", alpha=0.7,
          label="信号+背景 (二值采样)")
ax[0].bar(tc_ns, bg_hist_peak, width=bin_width*1e9, align="center", color="tab:red", alpha=0.6,
          label="纯背景 B (环境光二值)")
ax[0].axhline(macro_cap, color="k", ls="-.", lw=1.0, alpha=0.7, label=f"二值上限={macro_cap}")
ax[0].axvline(t0_ns, color="k", ls=":", alpha=0.7, label=f"真实 ToF {t0_ns:.1f} ns")
ax[0].annotate(f"峰 bin\nS={S_peak:.0f}, B={B_peak:.2f}\nSNR=S/sqrt(B)={SNR_sqrtB:.1f}",
               (tc_ns[pk_bin], tot_peak), xytext=(tc_ns[pk_bin]+4, tot_peak*0.8), fontsize=8,
               arrowprops=dict(arrowstyle="->", lw=0.8))
ax[0].set_xlabel("时间 t [ns]"); ax[0].set_ylabel(f"计数 / 1ns bin ({N_shots} shots)")
ax[0].set_title(f"峰值宏像元 m={m_peak}: 信号峰 S 与背景 B (v30 二值)")
ax[0].legend(fontsize=8); ax[0].grid(alpha=0.3)

# 右: 各宏像元峰 bin 的 SNR
ax[1].plot(np.arange(n_macro), snr_per_macro, lw=1.6, marker="o", ms=3, color="tab:purple")
ax[1].axvline(m_peak, color="lime", ls=":", lw=1.2, label=f"峰值宏像元 m={m_peak} (SNR={snr_per_macro[m_peak]:.1f})")
ax[1].set_xlabel("宏像元序号 m (沿长边 y)"); ax[1].set_ylabel("SNR = S/sqrt(B) @ 峰 bin")
ax[1].set_title(f"各宏像元信号峰 bin 的 SNR (背景 B_mean={B_mean:.2f})")
ax[1].legend(fontsize=9); ax[1].grid(alpha=0.3)
plt.tight_layout(); plt.show()

print(f"峰值宏像元 SNR(S/sqrt(B)) = {SNR_sqrtB:.2f}; 边缘宏像元 m=0 SNR = {snr_per_macro[0]:.2f}")
print("(定义: SNR=S/sqrt(B), S=扣背景信号计数, B=峰 bin 背景计数; 背景由纯环境光二值采样估计)")'''

i13, c13 = find_code_startswith("# ---- 信噪比 SNR = S / sqrt(B), 在峰值宏像元 m_peak 的信号峰 bin ----")
c13.source = code_13_new


# =============================================================================
# 6) 重写【模块 15】SNR 分布：改用二值引擎重复采样（统计/拟合/绘图逻辑不变）
# =============================================================================
code_15_new = r'''# ---- 多次蒙卡的 SNR 分布(v30 二值引擎重复跑峰值宏像元) ----
from scipy.stats import norm

N_TRIALS = 3000                     # SNR 采样次数(二值引擎重复跑; 2000~5000 约十几秒~1分钟)
rng_snr = np.random.default_rng(PARAMS["hist"]["seed"] + 777)
fvals_peak = macro_fvals[m_peak]    # 峰值宏像元 27 个 SPAD 的收集比例
T_LASER = laser_delay(PARAMS["timing"])

import time as _time
_t0 = _time.time()
snr_samples = np.empty(N_TRIALS)
S_samples   = np.empty(N_TRIALS)
B_samples   = np.empty(N_TRIALS)
for it in range(N_TRIALS):
    # 一次测量 = N_shots 累加, 27 SPAD 二值求和 -> 信号+背景 直方图
    h_sig = np.zeros(nbins); h_bg = np.zeros(nbins)
    for _shot in range(N_shots):
        acc_s = np.zeros(nbins, dtype=np.int32); acc_b = np.zeros(nbins, dtype=np.int32)
        for fij in fvals_peak:
            acc_s += spad_binary_trace(base_rate*fij, r_amb_ph, tf, centers, PDE, TAU_RC,
                                       VTH_FRAC, jit, rng_snr, T_OVER, T_LASER, RESP_SHAPE, RESP_K)
            # 同批纯背景(信号率=0), 用于该次测量的 B 估计
            acc_b += spad_binary_trace(zero_rate, r_amb_ph, tf, centers, PDE, TAU_RC,
                                       VTH_FRAC, jit, rng_snr, T_OVER, T_LASER, RESP_SHAPE, RESP_K)
        h_sig += acc_s; h_bg += acc_b
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
# 分母 B_bar = 全 bin 背景均值, 视为近似常数(nbins 个 bin 平均, 抖动很小); 只有 S~Poisson(S_bar) 随机。
#   => mu_SNR = S_bar/sqrt(B_bar);  sigma_SNR = sigma_S/sqrt(B_bar) = sqrt(S_bar)/sqrt(B_bar)
mu_fit, sig_fit = norm.fit(snr_valid)
S_bar = S_samples.mean(); B_bar = B_samples.mean()
mu_th  = S_bar / np.sqrt(B_bar)
sig_th = np.sqrt(S_bar) / np.sqrt(B_bar)            # S~Poisson: Var(S)=S_bar, B_bar 近似常数

print("="*76)
print(f"SNR 分布 (峰值宏像元 m={m_peak}, {N_TRIALS} 次二值蒙卡, 每次 27 SPAD x {N_shots} shots)")
print(f"  注: v30 二值下 S 受 27×N_shots={macro_cap} 封顶; 若峰未饱和, S 近似 Poisson, SNR 近似正态。")
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
ax[0].set_title(f"SNR 分布 + 正态拟合 (v30 二值, N_shots={N_shots})")
ax[0].legend(fontsize=8); ax[0].grid(alpha=0.3)
# 右: S 的分布 + 每次测量 B_bar 的分布(近似常数, 抖动小)
ax[1].hist(S_samples, bins=nb_h, density=True, color="tab:green", alpha=0.5, label=f"S 样本 (mean={S_bar:.1f})")
ax[1].axvline(S_bar, color="tab:green", ls=":", alpha=0.8)
ax1b = ax[1].twiny()
ax1b.hist(B_samples, bins=30, density=True, color="tab:red", alpha=0.45,
          label=f"B_bar 样本 (mean={B_bar:.3f})")
ax1b.set_xlabel("每次测量的全 bin 背景均值 B_bar", color="tab:red")
ax[1].set_xlabel("峰 bin 信号计数 S"); ax[1].set_ylabel("概率密度")
ax[1].set_title("信号 S 与 每次 B_bar(近似常数) 的分布")
h1, l1 = ax[1].get_legend_handles_labels(); h2, l2 = ax1b.get_legend_handles_labels()
ax[1].legend(h1+h2, l1+l2, fontsize=8); ax[1].grid(alpha=0.3)
plt.tight_layout(); plt.show()

print(f"结论: v30 二值下 SNR 分布近似正态(峰未饱和时); 理论(mu={mu_th:.2f},sigma={sig_th:.2f}) 与 拟合(mu={mu_fit:.2f},sigma={sig_fit:.2f}) 吻合。")
print(f"      分母用全 bin 背景均值 B_bar(近似常数): SNR=S/sqrt(B_bar) 是 S 的线性缩放。")'''

i15, c15 = find_code_startswith("# ---- 多次蒙卡的 SNR 分布(RC 引擎重复跑峰值宏像元) ----")
c15.source = code_15_new


# =============================================================================
# 7) 最前面插入 v30 说明标题 cell（不改动 v21 原标题 cell）
# =============================================================================
md_title = nbformat.v4.new_markdown_cell(r'''# LiDAR 直方图仿真 v30 —— 二值采样模型（8ns 过阈窗 + 1ns 时钟 0/1）

> **本版由 v21 派生。** 与 v21 最大的差异在"SPAD 统计方式"：
>
> | | v21（旧） | **v30（本版）** |
> |---|---|---|
> | 记录量 | 每个光子的 **timestamp** | **不记 timestamp**；1ns 时钟采 **0/1** |
> | 单光子响应 | 一个事件（δ） | **8ns 宽、高度=1** 的过阈值窗 |
> | 直方图 | timestamp 落 bin 计数 | 采样点是否过阈（每 SPAD 每 bin 每 shot ∈ {0,1}） |
> | 峰 bin 上限 | 无硬上限 | **27 × N_shots**（饱和） |
> | 展宽/饱和 | — | 由 8ns 窗（展宽）与二值封顶+窗堆积（饱和）统一解释 |
> | 亚 ns 测距 | — | `delta_dly` 游标（1/12 ns）扫描（模块 12b） |
>
> **保留复用 v21：** 模块 0–9 的物理链路（激光脉冲、发射/接收光学、信道、朗伯目标、
> 椭圆高斯像斑、光链路、环境光、RC 恢复机制）原样不变；其中模块 8b 的 Vov(t) RC
> 恢复曲线正是 **8ns 过阈窗宽 `T_OVER = −τ_RC·ln(1−Vth_frac) ≈ 8ns`** 的物理来源。
>
> **v30 核心新增/改写：** 模块 9b（二值引擎 + 波形演示）、模块 11/12/13/15（改用二值
> 采样产 `macro_hist`/`bg_hist_peak`，图12 删去 timestamp eventplot、加饱和上限线）、
> 模块 12b（`delta_dly` 亚 ns 测距）。下游模块 14/16 依赖契约变量无缝复用。
>
> ⚠️ 物理场景参数（目标距离/反射率/脉冲/bin 宽/N_shots 等）与 v21 **完全一致**，
> 仅新增 `PARAMS["timing"]` 三参数（默认全 0，即 `t_laser=0`，回波峰位与 v21 一致）。''')
nb.cells.insert(0, md_title)


# =============================================================================
# 写出 v30 notebook
# =============================================================================
nbformat.write(nb, DST)
print(f"[OK] 已生成 {DST}，共 {len(nb.cells)} 个 cell（源 {SRC} 未改动）。")
