# -*- coding: utf-8 -*-
"""
基于 v12 生成 v13:
- 保留 v12 模块 0-7(参数 + 物理链路 + 硬死时间蒙卡引擎 simulate_spad_shot), 清空输出;
- 新增模块 7b: RC 恢复 SPAD 蒙卡引擎 simulate_spad_shot_rc(过电压 Vov 按 RC 指数恢复,
  PDE∝Vov, 计数需 Vov>=Vth, 仅计数事件复位);
- 模块 8/9/11/12 改用 RC 引擎, 并与 v12 硬死时间(14ns)对比。
参数(用户确认): τ_RC=6ns, Vth=10% Vov_max, reset='count'(仅计数事件复位)。
运行: python build_v13_from_v12.py
"""
import nbformat

SRC = "lidar_histogram_sim_v12.ipynb"
DST = "lidar_histogram_sim_v13.ipynb"

nb = nbformat.read(SRC, as_version=4)

# --- 保留 v12 中直到(不含)"## 模块 8"的所有 cell(模块 0-7, 含硬死时间引擎) ---
kept = []
for c in nb.cells:
    if c.cell_type == "markdown" and c.source.lstrip().startswith("## 模块 8"):
        break
    kept.append(c)

for c in kept:
    if c.cell_type == "code":
        c.outputs = []
        c.execution_count = None
    if "execution" in c.get("metadata", {}):
        del c.metadata["execution"]

# 替换首个标题 cell 为 v13 标题
if kept and kept[0].cell_type == "markdown":
    kept[0].source = r'''# 激光雷达直方图仿真 v13 (LiDAR Histogram Simulation) — RC 恢复 SPAD 模型

链路(激光→TX→信道→目标→RX)与 **v4 完全一致**; 逐光子蒙特卡洛(Monte Carlo, 蒙卡)、
宏像元(Macro Pixel)统计沿用 **v10/v11/v12**(模块 0–7 原样复用, 含 v12 的硬死时间引擎作对照)。

**v13 相对 v12 的变化 —— 用 RC 恢复模型取代硬死时间**
- v12 用**硬死时间(hard dead time)**: 雪崩后固定 14ns 内完全不响应, 之后瞬间满灵敏。
  这抹掉了 SPAD 恢复期的**部分灵敏**——真实器件恢复期仍能以较低概率触发较弱雪崩。
- v13 新增 **RC 恢复模型**(模块 7b, 更接近物理): 雪崩后过电压
  **Vov(OverVoltage, 过电压 = V_bias − V_br)** 通过淬灭电阻 R 与结电容 C_J 按 **RC 指数**恢复:
  $$V_{ov}(\Delta t)=V_{ov,max}\,(1-e^{-\Delta t/\tau}),\qquad \tau=R\,C_J$$
  - 恢复期来光子的触发概率 **PDE(Δt)=PDE_max · Vov(Δt)/Vov_max**(线性 PDE∝Vov);
  - 触发后**只有脉冲幅度(∝当前 Vov)≥ 阈值 Vth 才被读出电路计一次数**;
  - **仅计数事件复位**: 只有够阈值的计数事件才重置 Vov、重新开始 RC 恢复(亚阈弱雪崩不复位)。
- **参数**(本版设定): τ_RC = **6 ns**, Vth = **10% · Vov_max**, Vov_max = 3.3 V, reset = **仅计数事件**。
- 模块 8/9/11/12 全部改用 RC 引擎, 并**与 v12 硬死时间(14ns)对比**, 展示两模型差异。

**为什么 RC 更可靠(尤其本场景)**
- 本场景中心 SPAD 单 shot 入射约 64 光子、脉宽约 6ns(远短于恢复时间), 光子极密集。
- 硬死时间把恢复期一刀切死 ⇒ 严重低估计数(峰值区仅剩~6%); RC 允许恢复途中被较弱触发 ⇒ 更真实(~18%)。

> ⚠️ 运行方式: **Kernel → Restart & Run All**(从上到下顺序执行)。
> 宏像元方向锁定: **3 沿长边 y, 9 沿短边 x**(9×3, 沿长边 40 个, 每个 27 SPAD)。'''

# ---------------- v13 新增/替换 cells ----------------
NEW = []
def md(s): NEW.append(("md", s))
def code(s): NEW.append(("code", s))

# ===== 模块 7b: RC 引擎 =====
md(r'''## 模块 7b（v13 新增）— RC 恢复 SPAD 蒙特卡洛引擎

雪崩后过电压 Vov 不是"突变恢复", 而是按 RC 指数充回(见 SPAD 等效电路: 淬灭电阻 R + 结电容 C_J):
$$V_{ov}(\Delta t)=V_{ov,max}\,(1-e^{-\Delta t/\tau}),\qquad \tau=R\,C_J$$
其中 Δt = 距上次**计数事件**经过的时间(仅计数事件复位)。

单次 shot 逐光子处理:
1. 每精细 bin 入射光子 ~Poisson(r·dt), 展开为逐光子到达时刻。
2. 光子到达时, 当前 `vov_frac = 1 - exp(-Δt/τ)`; 触发概率 `PDE_max · vov_frac`(线性 PDE∝Vov)。
3. 触发后: 若 `vov_frac ≥ Vth_frac` ⇒ **计一次数**, 记录时间戳(加 IRF 抖动), 并复位(Δt 归零);
   否则(亚阈)⇒ 不计数; `reset='all'` 时也复位, `reset='count'` 时不复位。

> 与硬死时间的关键区别: RC 恢复期是**渐变灵敏**(概率随 Vov 连续上升), 而非"非 0 即 1"的硬开关。''')

code(r'''def simulate_spad_shot_rc(r_sig_fine, r_amb_ph, tf, PDE_max, tau_rc, Vth_frac,
                          jitter_sigma, rng, reset_mode="count", Vov_max=3.3):
    """RC 恢复 SPAD 单次 shot 蒙卡, 返回被计数事件的记录时间戳数组 [s]。
    r_sig_fine: 信号光子率(不含 PDE)在精细网格 tf 上; r_amb_ph: 环境光光子率(标量, 不含 PDE)。
    tau_rc: RC 恢复时间常数 = R·C_J; Vth_frac: 计数所需最小 Vov 占 Vov_max 比例;
    reset_mode: 'count'=仅计数事件复位 Vov; 'all'=任何雪崩(含亚阈)都复位。"""
    dt = tf[1] - tf[0]
    mu = (r_sig_fine + r_amb_ph) * dt
    n_ph = rng.poisson(mu)
    if n_ph.sum() == 0:
        return np.empty(0)
    t_arr = np.repeat(tf, n_ph)            # 逐光子到达时刻(升序)
    u = rng.random(t_arr.size)             # 每光子随机数(判触发)
    det = []
    last = -1e30                            # 上次复位时刻
    inv_tau = 1.0 / tau_rc
    for k in range(t_arr.size):
        t = t_arr[k]
        d = (t - last) * inv_tau
        vov_frac = 1.0 - np.exp(-d) if d < 700 else 1.0   # 当前 Vov/Vov_max
        if u[k] < PDE_max * vov_frac:                     # 雪崩触发(概率∝Vov)
            if vov_frac >= Vth_frac:                      # 幅度够阈 -> 计数
                det.append(t)
                last = t                                  # 计数事件复位 Vov
            elif reset_mode == "all":
                last = t                                  # 亚阈雪崩也复位(仅 all 模式)
    det = np.asarray(det)
    if jitter_sigma > 0 and det.size:
        det = det + rng.normal(0.0, jitter_sigma, det.size)   # IRF 抖动只加在记录时间戳
    return det

# ---- RC 参数(本版设定) ----
TAU_RC   = 6e-9        # RC 恢复时间常数 = R·C_J
VTH_FRAC = 0.10        # 计数所需最小 Vov (占 Vov_max)
VOV_MAX  = 3.3         # 满过电压 [V]
RESET    = "count"     # 仅计数事件复位

# 等效"硬死区"(恢复到 Vth 前完全不计数)与"渐变灵敏"分界
t_deadzone = -np.log(1 - VTH_FRAC) * TAU_RC
print(f"RC 引擎就绪: τ_RC={TAU_RC*1e9:.1f} ns, Vth={VTH_FRAC*100:.0f}%·Vov_max, reset='{RESET}'")
print(f"  Vov 恢复曲线: 1τ->{100*(1-np.exp(-1)):.0f}%, 2.3τ->{100*(1-np.exp(-2.3)):.0f}%, 5τ->{100*(1-np.exp(-5)):.1f}%")
print(f"  低于 Vth 的'硬死区'≈{t_deadzone*1e9:.2f} ns, 之后为渐变灵敏(与硬死时间 14ns 的一刀切不同)")''')

# ===== 模块 8(替换: RC vs 硬死时间)=====
md(r'''## 模块 8（v13 更新）— 单 SPAD: RC 恢复模型 vs 硬死时间（timestamp + 1ns 直方图）

取像斑中心 SPAD, 只在 30m 目标 ToF 附近仿真。**同一批入射光子**下对比两种模型:
- **图 1（timestamp）**: 单次 shot, RC 模型 vs 硬死时间(14ns) 各自探测到的光子时间戳(raster)。
- **图 2（1ns 直方图 + 均值）**: 多 shot 蒙卡均值, RC vs 硬死时间 vs 解析上限(无死时间)。

> RC 恢复期渐变灵敏 ⇒ 峰值区计数应**高于**硬死时间、**低于**无死时间上限。''')

code(r'''# ---- 时间窗(与 v12 相同, 参数不变) ----
t0 = time_of_flight(D0)
pre, post = 10e-9, 30e-9
dt_fine  = PARAMS["hist"]["dt_fine"]
bin_width = PARAMS["hist"]["bin_width"]
t_lo, t_hi = t0 - pre, t0 + post
tf = np.arange(t_lo, t_hi, dt_fine)

PDE = PARAMS["spad"]["PDE"]; jit = PARAMS["spad"]["jitter_sigma"]
t_dead = 14e-9                                # v12 硬死时间(作对照)

f_ij = fpix0[i0, j0]
r_sig = signal_photon_rate_fine(echo0, f_ij, tf)

# ---- (1) 单次 shot: RC vs 硬死时间(同 seed) ----
rng_rc  = np.random.default_rng(PARAMS["hist"]["seed"])
ev1_rc  = simulate_spad_shot_rc(r_sig, r_amb_ph, tf, PDE, TAU_RC, VTH_FRAC, jit, rng_rc, RESET)
rng_hd  = np.random.default_rng(PARAMS["hist"]["seed"])
ev1_hd  = simulate_spad_shot(r_sig, r_amb_ph, tf, PDE, t_dead, jit, rng_hd)
ev1_rc_ns = np.sort(ev1_rc)*1e9; ev1_hd_ns = np.sort(ev1_hd)*1e9
centers, c1_rc, edges = events_to_hist(ev1_rc, t_lo, t_hi, bin_width)
tc_ns = centers*1e9; t0_ns = t0*1e9

print("="*76)
print(f"单 SPAD (i0,j0)=({i0},{j0}), f_pix={f_ij:.3e}; 30m ToF={t0_ns:.1f} ns")
print(f"单次 shot: RC 模型探测 {ev1_rc.size} 个, 硬死时间(14ns)探测 {ev1_hd.size} 个")
print(f"  RC timestamp[ns]: " + (", ".join(f"{t:.2f}" for t in ev1_rc_ns) if ev1_rc.size else "(无)"))
print(f"  硬死 timestamp[ns]: " + (", ".join(f"{t:.2f}" for t in ev1_hd_ns) if ev1_hd.size else "(无)"))

# ---- (2) 多 shot 均值: RC / 硬死 / 无死上限 ----
Nrep = 3000
rng_a = np.random.default_rng(PARAMS["hist"]["seed"]+1)
rng_b = np.random.default_rng(PARAMS["hist"]["seed"]+1)
acc_rc = np.zeros(len(centers)); acc_hd = np.zeros(len(centers))
for _ in range(Nrep):
    e = simulate_spad_shot_rc(r_sig, r_amb_ph, tf, PDE, TAU_RC, VTH_FRAC, jit, rng_a, RESET)
    acc_rc += np.histogram(e, bins=edges)[0]
for _ in range(Nrep):
    e = simulate_spad_shot(r_sig, r_amb_ph, tf, PDE, t_dead, jit, rng_b)
    acc_hd += np.histogram(e, bins=edges)[0]
mc_rc = acc_rc/Nrep; mc_hd = acc_hd/Nrep

# 解析无死时间上限
r_det = (r_sig + r_amb_ph) * PDE
r_det = np.convolve(r_det, gaussian_kernel(jit, dt_fine), mode="same") * dt_fine
bin_idx = np.clip(((tf - t_lo)/bin_width).astype(int), 0, len(centers)-1)
lam = np.bincount(bin_idx, weights=r_det*dt_fine, minlength=len(centers))

print(f"窗口每 shot 总计数: 无死上限={lam.sum():.3f}, RC={mc_rc.sum():.3f}, 硬死14ns={mc_hd.sum():.3f}")
print(f"  -> RC 保留 {100*mc_rc.sum()/lam.sum():.1f}%, 硬死时间保留 {100*mc_hd.sum()/lam.sum():.1f}% (RC 更接近真实)")
pk = int(lam.argmax())
print(f"峰值 bin @ {tc_ns[pk]:.0f} ns: 无死={lam[pk]:.3f}, RC={mc_rc[pk]:.3f}, 硬死={mc_hd[pk]:.3f}")

# ---- 绘图 1: timestamp raster ----
fig, ax = plt.subplots(figsize=(11, 3.0))
ax.eventplot([ev1_rc_ns, ev1_hd_ns], lineoffsets=[2,1], linelengths=0.8,
             colors=["tab:green","tab:red"])
ax.axvline(t0_ns, color="k", ls=":", label=f"真实 ToF {t0_ns:.1f} ns")
ax.set_yticks([1,2]); ax.set_yticklabels(["硬死 14ns","RC 恢复"])
ax.set_xlim(t_lo*1e9, t_hi*1e9); ax.set_xlabel("时间 t [ns]")
ax.set_title(f"单个 SPAD · 单次 shot 探测 timestamp: RC 恢复(τ={TAU_RC*1e9:.0f}ns) vs 硬死时间(14ns)")
ax.legend(fontsize=9, loc="upper right"); ax.grid(alpha=0.3, axis="x")
plt.tight_layout(); plt.show()

# ---- 绘图 2: 直方图 + 均值对比 ----
fig, ax = plt.subplots(2, 1, figsize=(11, 7), sharex=True)
ax[0].bar(tc_ns, c1_rc, width=bin_width*1e9, align="center", color="tab:green", alpha=0.8)
ax[0].axvline(t0_ns, color="r", ls=":", label=f"真实 ToF {t0_ns:.1f} ns")
ax[0].set_ylabel("计数 / 1ns bin")
ax[0].set_title(f"单 SPAD · 单次 shot 直方图 (RC 模型, τ={TAU_RC*1e9:.0f}ns, Vth={VTH_FRAC*100:.0f}%)")
ax[0].legend(fontsize=9); ax[0].grid(alpha=0.3)

ax[1].plot(tc_ns, lam, "k--", lw=1.5, label="解析上限(无死时间)")
ax[1].plot(tc_ns, mc_rc, color="tab:green", lw=1.8, marker="o", ms=3, label=f"RC 恢复 τ={TAU_RC*1e9:.0f}ns")
ax[1].plot(tc_ns, mc_hd, color="tab:red", lw=1.8, marker="s", ms=3, label="硬死时间 14ns")
ax[1].axvline(t0_ns, color="k", ls=":", alpha=0.6)
ax[1].set_xlabel("时间 t [ns]"); ax[1].set_ylabel(f"每 shot 平均计数 / bin ({Nrep} shots)")
ax[1].set_title("多 shot 均值: RC 介于无死上限与硬死时间之间 -> RC 更真实")
ax[1].legend(fontsize=9); ax[1].grid(alpha=0.3)
plt.tight_layout(); plt.show()''')

# ===== 模块 9(RC 版)=====
md(r'''## 模块 9（v13 更新）— 阵列内不同 SPAD 的响应差异（RC 恢复模型）

像斑椭圆高斯 ⇒ 每个 SPAD 收集比例不同。用 RC 模型比较沿长边 y 不同偏移的 SPAD。
- **左**: 像元收集比例 f_pix 热图(聚焦被照区)。
- **右**: 沿长边 y 取中心/中段/边缘代表 SPAD, RC 模型下各自多 shot 平均计数-时间。''')

code(r'''thr = fpix0.max() * 0.01
xs = np.where(fx0 > fx0.max()*0.01)[0]; ys = np.where(fy0 > fy0.max()*0.01)[0]
ix = slice(max(0, xs.min()-1), min(fpix0.shape[0], xs.max()+2))
iy = slice(max(0, ys.min()-1), min(fpix0.shape[1], ys.max()+2))
sub = fpix0[ix, iy]

offs = [0, 10, 25, 45]
offs = [d for d in offs if 0 <= j0 + d < fpix0.shape[1]]
reps = 800
rng9 = np.random.default_rng(PARAMS["hist"]["seed"] + 7)

fig, ax = plt.subplots(1, 2, figsize=(13, 4.6))
im = ax[0].imshow(sub.T, origin="lower", aspect="auto", cmap="viridis",
                  extent=[ix.start, ix.stop-1, iy.start, iy.stop-1])
ax[0].scatter([i0], [j0], c="r", marker="x", s=80, label="中心 SPAD")
ax[0].set_xlabel("像元 x 序号 (短边)"); ax[0].set_ylabel("像元 y 序号 (长边)")
ax[0].set_title("像元空间收集比例 f_pix (聚焦被照区)")
ax[0].legend(fontsize=8); plt.colorbar(im, ax=ax[0], label="f_pix")

for d in offs:
    j = j0 + d
    f_ij_d = fpix0[i0, j]
    r_sig_d = signal_photon_rate_fine(echo0, f_ij_d, tf)
    acc = np.zeros(len(centers))
    for _ in range(reps):
        ev = simulate_spad_shot_rc(r_sig_d, r_amb_ph, tf, PDE, TAU_RC, VTH_FRAC, jit, rng9, RESET)
        acc += np.histogram(ev, bins=edges)[0]
    ax[1].plot(tc_ns, acc/reps, lw=1.6, marker="o", ms=3,
               label=f"Δy={d} 像元 (f={f_ij_d:.2e})")
ax[1].axvline(t0_ns, color="k", ls=":", alpha=0.6, label=f"ToF {t0_ns:.1f} ns")
ax[1].set_xlabel("时间 t [ns]"); ax[1].set_ylabel(f"每 shot 平均计数 / bin ({reps} shots)")
ax[1].set_title("不同 SPAD(沿长边 y 偏移)的计数-时间 (RC 模型): 越偏离中心收集越少")
ax[1].legend(fontsize=8); ax[1].grid(alpha=0.3)
plt.tight_layout(); plt.show()

print("说明: RC 模型下, 中心 SPAD 因光子密集、恢复期部分灵敏, 计数比硬死时间高; 边缘 SPAD 光子稀疏, 两模型接近。")''')

# ===== 模块 10(沿用 v12: 去图留字)=====
md(r'''## 模块 10（沿用 v12）— Macro Pixel（宏像元）定义 9×3（仅文字说明与数值）

**定义**: 短边 x 的 9 个 SPAD 全取, 长边 y 每 3 个为一组。
阵列 `9×120` ⇒ 沿长边共 `120/3 = 40` 个宏像元, 每个 `9×3 = 27` 个 SPAD。
宏像元 `m` 覆盖全部 x、y 索引 `[3m, 3m+3)`。像斑长边 σ≈200µm ⇒ 中间宏像元信号多、两端仅底噪。''')

code(r'''Bx_m, By_m = 9, 3
a = PARAMS["spad_array"]
assert a["Nx"] == Bx_m and a["Ny"] % By_m == 0
n_macro = a["Ny"] // By_m
n_pix_macro = Bx_m * By_m
macro_fsum = np.array([fpix0[:, m*By_m:(m+1)*By_m].sum() for m in range(n_macro)])
m_peak = int(macro_fsum.argmax())
print(f"宏像元 = {Bx_m}×{By_m}; 共 {n_macro} 个, 每个 {n_pix_macro} SPAD; 峰值宏像元 m_peak={m_peak} (Σf={macro_fsum[m_peak]:.3f})")
print(f"全部 Σf_pix={macro_fsum.sum():.3f}; 边缘 m=0 的 Σf={macro_fsum[0]:.3e} (近乎纯底噪)")''')

# ===== 模块 11(RC 版宏像元)=====
md(r'''## 模块 11（v13 更新）— 逐-SPAD 蒙卡统计每个宏像元直方图（RC 恢复模型）

每个宏像元 27 个 SPAD 各自独立跑 **RC 引擎**(每 SPAD 各自 Vov 恢复), 合并事件后 1ns 直方图, 累加 N_shots。
同时给出 v12 硬死时间(14ns)的全体总计数做对比, 展示 RC 模型对计数损失的修正。''')

code(r'''nbins = len(centers)
base_rate = signal_photon_rate_fine(echo0, 1.0, tf)
N_shots = PARAMS["hist"]["N_shots"]
macro_fvals = [fpix0[:, m*By_m:(m+1)*By_m].ravel() for m in range(n_macro)]

# ---- RC 模型: 逐-SPAD 蒙卡 ----
rng = np.random.default_rng(PARAMS["hist"]["seed"] + 11)
macro_hist = np.zeros((n_macro, nbins))
for _shot in range(N_shots):
    for m in range(n_macro):
        ev_all = []
        for fij in macro_fvals[m]:
            ev = simulate_spad_shot_rc(base_rate*fij, r_amb_ph, tf, PDE, TAU_RC, VTH_FRAC, jit, rng, RESET)
            if ev.size: ev_all.append(ev)
        if ev_all:
            macro_hist[m] += np.histogram(np.concatenate(ev_all), bins=edges)[0]

# ---- 硬死时间(14ns)对比: 只算全体总计数 ----
rng2 = np.random.default_rng(PARAMS["hist"]["seed"] + 11)
tot_hd = 0
for _shot in range(N_shots):
    for m in range(n_macro):
        for fij in macro_fvals[m]:
            ev = simulate_spad_shot(base_rate*fij, r_amb_ph, tf, PDE, t_dead, jit, rng2)
            tot_hd += ev.size

# ---- 解析上限(无死时间) ----
irf_k = gaussian_kernel(jit, dt_fine)
bin_idx = np.clip(((tf - t_lo)/bin_width).astype(int), 0, nbins-1)
macro_lam = np.zeros((n_macro, nbins))
for m in range(n_macro):
    r_det = (base_rate*macro_fsum[m] + n_pix_macro*r_amb_ph) * PDE
    r_det = np.convolve(r_det, irf_k, mode="same") * dt_fine
    macro_lam[m] = N_shots * np.bincount(bin_idx, weights=r_det*dt_fine, minlength=nbins)

tot_rc, tot_la = macro_hist.sum(), macro_lam.sum()
print("="*76)
print(f"逐-SPAD 蒙卡: {n_macro} 宏 × {n_pix_macro} SPAD × {N_shots} shots (RC: τ={TAU_RC*1e9:.0f}ns, Vth={VTH_FRAC*100:.0f}%)")
print(f"全体总计数: 无死上限={tot_la:.0f}, RC={tot_rc:.0f}, 硬死14ns={tot_hd:.0f}")
print(f"  -> 保留率: RC={100*tot_rc/tot_la:.1f}%, 硬死时间={100*tot_hd/tot_la:.1f}% (RC 高出 {tot_rc-tot_hd:.0f} 计数)")
print(f"峰值宏像元 m={m_peak}: RC 峰值 bin={macro_hist[m_peak].max():.0f}, 总计数={macro_hist[m_peak].sum():.0f}; "
      f"边缘 m=0 总计数={macro_hist[0].sum():.0f}")''')

# ===== 模块 12(RC 版可视化)=====
md(r'''## 模块 12（v13 更新）— 宏像元-时间热图（RC 模型）+ 峰值宏像元 27 SPAD timestamp

- **图 A**: 40 宏像元 × 时间(1ns bin) RC 模型计数热图。
- **图 B**: 代表性宏像元直方图(柱=RC 蒙卡, 虚线=无死上限参考)。
- **图 C**: 峰值宏像元单次 shot, 27 个 SPAD 各自 RC 探测 timestamp(每行一 SPAD)。''')

code(r'''# ---- 图 A: 热图 ----
fig, ax = plt.subplots(figsize=(11, 5))
im = ax.imshow(macro_hist, origin="upper", aspect="auto", cmap="inferno",
               extent=[tc_ns[0]-0.5, tc_ns[-1]+0.5, n_macro-0.5, -0.5])
ax.axvline(t0_ns, color="cyan", ls=":", lw=1.2, label=f"真实 ToF {t0_ns:.1f} ns")
ax.axhline(m_peak, color="lime", ls=":", lw=1.0, alpha=0.7, label=f"峰值宏像元 m={m_peak}")
ax.set_xlabel("时间 t [ns]"); ax.set_ylabel("宏像元序号 m (0=顶部, 沿长边 y)")
ax.set_title(f"宏像元(9×3)直方图热图 (RC 模型 τ={TAU_RC*1e9:.0f}ns, N_shots={N_shots}, PDE={PDE})")
ax.legend(fontsize=9, loc="upper right"); plt.colorbar(im, ax=ax, label="计数 / (宏像元, 1ns bin)")
plt.tight_layout(); plt.show()

# ---- 图 B: 代表性宏像元直方图(RC 蒙卡: 标记+连线; 上限: 虚线) ----
fig, ax = plt.subplots(figsize=(11, 4.6))
reps_m = sorted(set([0, max(0, m_peak-6), m_peak, min(n_macro-1, m_peak+6)]))
colors = ["tab:gray", "tab:green", "tab:red", "tab:orange", "tab:purple"]
for m, c in zip(reps_m, colors):
    ax.plot(tc_ns, macro_hist[m], color=c, lw=1.4, marker="o", ms=4,
            label=f"m={m} RC 蒙卡 (Σf={macro_fsum[m]:.3f})")            # 蒙卡数据: 标记+连线(不画柱状)
    ax.plot(tc_ns, macro_lam[m], color=c, lw=1.1, ls="--", alpha=0.6)    # 无死时间上限: 虚线参考
ax.axvline(t0_ns, color="k", ls=":", alpha=0.6, label=f"ToF {t0_ns:.1f} ns")
ax.set_xlabel("时间 t [ns]"); ax.set_ylabel(f"计数 / 1ns bin ({N_shots} shots)")
ax.set_title("代表性宏像元直方图 (点线=RC 蒙卡, 虚线=无死时间上限参考)")
ax.legend(fontsize=8); ax.grid(alpha=0.3)
plt.tight_layout(); plt.show()

# ---- 图 C: 峰值宏像元 27 SPAD RC timestamp ----
rng_ts = np.random.default_rng(PARAMS["hist"]["seed"] + 200)
sp_events = []
for fij in macro_fvals[m_peak]:
    ev = simulate_spad_shot_rc(base_rate*fij, r_amb_ph, tf, PDE, TAU_RC, VTH_FRAC, jit, rng_ts, RESET)
    sp_events.append(np.sort(ev)*1e9)
n_ev_total = sum(len(e) for e in sp_events)

fig, ax = plt.subplots(figsize=(11, 5))
ax.eventplot(sp_events, lineoffsets=np.arange(1, n_pix_macro+1), linelengths=0.8, colors="tab:green")
ax.axvline(t0_ns, color="r", ls=":", label=f"真实 ToF {t0_ns:.1f} ns")
ax.set_xlim(t_lo*1e9, t_hi*1e9); ax.set_ylim(0.5, n_pix_macro+0.5)
ax.set_xlabel("时间 t [ns]"); ax.set_ylabel("宏像元内 SPAD 序号 (1..27)")
ax.set_title(f"峰值宏像元 m={m_peak} · 单次 shot · 27 SPAD 各自 RC 探测 timestamp (τ={TAU_RC*1e9:.0f}ns)")
ax.legend(fontsize=9); ax.grid(alpha=0.3, axis="x")
plt.tight_layout(); plt.show()

print(f"峰值宏像元 m={m_peak} 单次 shot(RC): 27 个 SPAD 共探测 {n_ev_total} 个光子。")
print(f"每个 SPAD 各自 Vov 按 RC 恢复(τ={TAU_RC*1e9:.0f}ns), 恢复期渐变灵敏, 合并后 1ns 直方图。")''')

# --- 组装 & 写出 ---
for typ, src in NEW:
    kept.append(nbformat.v4.new_markdown_cell(src) if typ == "md"
                else nbformat.v4.new_code_cell(src))
nb.cells = kept
nbformat.write(nb, DST)
print(f"已生成 {DST}: 共 {len(kept)} 个 cell (保留 v12 模块0-7 + 新增 {len(NEW)})。")
