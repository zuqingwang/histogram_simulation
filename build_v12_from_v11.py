# -*- coding: utf-8 -*-
"""
基于 v11 生成 v12:
- 保留 v11 模块 0-7(参数 + 物理链路 + 单 SPAD 蒙卡引擎, 引擎已支持 dead time>0), 清空输出;
- 模块 8/9/11/12 启用 SPAD dead time = 14 ns(不再是 0);
- 模块 8: 增加"单 SPAD 探测光子 timestamp"输出(单次 shot + 多次 shot raster), 1ns 直方图保留;
- 模块 10: 去掉收集比例图, 只保留文字说明与数值(计算供 11/12 复用);
- 模块 11: dead time=14ns 下逐-SPAD 蒙卡, 校验改为量化 dead time 计数损失;
- 模块 12: 删"各宏像元总计数柱状图"(图 B 左), 增加"峰值宏像元 27 个 SPAD 各自 timestamp"。
运行: python build_v12_from_v11.py
"""
import nbformat

SRC = "lidar_histogram_sim_v11.ipynb"
DST = "lidar_histogram_sim_v12.ipynb"

nb = nbformat.read(SRC, as_version=4)

# --- 保留 v11 中直到(不含)"## 模块 8"的所有 cell(即模块 0-7) ---
kept = []
for c in nb.cells:
    if c.cell_type == "markdown" and c.source.lstrip().startswith("## 模块 8"):
        break
    kept.append(c)

# 清空保留 cell 的输出/执行序号
for c in kept:
    if c.cell_type == "code":
        c.outputs = []
        c.execution_count = None
    if "execution" in c.get("metadata", {}):
        del c.metadata["execution"]

# 替换首个标题 cell 为 v12 标题
if kept and kept[0].cell_type == "markdown":
    kept[0].source = r'''# 激光雷达直方图仿真 v12 (LiDAR Histogram Simulation) — dead time 实装 + 探测 timestamp

链路(激光→TX→信道→目标→RX)与 **v4 完全一致**; 单 SPAD(Single-Photon Avalanche Diode,
单光子雪崩二极管)逐光子蒙特卡洛(Monte Carlo, 蒙卡)引擎与 **v10 完全一致**;
宏像元(Macro Pixel)统计与 **v11 完全一致**(模块 0–7 原样复用)。

**v12 相对 v11 的变化**
- **实装 dead time(死时间)= 14 ns**(non-paralyzable, 非瘫痪型): 一次雪崩后 14 ns 内到达的光子被忽略,
  且不延长死时间; 死时间从**雪崩时刻**算起, IRF(Instrument Response Function, 仪器响应函数)抖动只加在记录时间戳上。
  这会引入 **pile-up(脉冲堆积效应: 前面的光子占住探测器, 后到光子被系统性漏计, 直方图前沿被抬高)**。
- **模块 8**: 增加"单个 SPAD 探测光子 **timestamp**(时间戳, 含 dead time + IRF 抖动)"的直接输出——
  单次 shot 的事件时刻列表 + 多次 shot 的 raster(每行一次 shot); **1ns bin 直方图保留**。
- **模块 10**: 只保留宏像元定义的文字说明与数值(供模块 11/12 复用), **不再绘制收集比例图**。
- **模块 12**: 删除"各宏像元总计数柱状图"(原图 B 左); 增加"峰值宏像元 27 个 SPAD 各自 timestamp"raster。

**关于噪声(信号也是 Poisson)**
- 信号光子与环境光子到达探测器都服从 **Poisson 过程(散粒噪声 shot noise, 光子离散性的量子涨落)**;
  即便脉冲能量完全确定, 单次 shot 实际探测光子数仍随机涨落(相干态测光子数严格为 Poisson)。
- 我们计算的率函数 r(t) 只是**期望**; 每个精细 bin 入射光子数 ~Poisson(r·dt), 再逐光子按
  PDE(Photon Detection Efficiency, 光子探测效率)判定是否触发雪崩。**这正是"逐光子蒙卡"的前提。**

> ⚠️ 运行方式: **Kernel → Restart & Run All**(从上到下顺序执行)。
> 说明: 宏像元方向锁定为 **3 沿长边 y, 9 沿短边 x**(9×3, 沿长边 40 个, 每个 27 SPAD)。'''

# ---------------- v12 新增/替换 cells ----------------
NEW = []
def md(s): NEW.append(("md", s))
def code(s): NEW.append(("code", s))

# ===== 模块 8(替换)=====
md(r'''## 模块 8（v12 更新）— 单 SPAD 探测 timestamp（含 dead time=14ns + IRF 抖动）+ 1ns 直方图

只在 30m 目标 ToF(Time of Flight, 飞行时间)附近的时间窗内仿真。取像斑中心 SPAD:
- **图 1（timestamp）**: 直接给出"某一次、某个 SPAD、探测到的光子**时间戳**"——
  上=单次 shot 的事件竖线(raster), 下=多次 shot 堆叠(每行一次 shot), 展示 shot 间的随机与稀疏。
- **图 2（1ns 直方图 + 验证）**: 上=单次 shot 的 1ns bin 直方图; 下=多 shot 蒙卡均值(含 14ns dead time)
  vs 解析期望 λ(dead time=0 参考)。二者差值即 **dead time 造成的计数损失(pile-up)**。

> dead time 从**雪崩时刻**算起(未加抖动); IRF 抖动(σ=100ps)只叠加在被记录的 timestamp 上。''')

code(r'''# ---- 时间窗: 只在 30m 目标 ToF 附近(与 v10/v11 相同, 参数不变) ----
t0 = time_of_flight(D0)
pre, post = 10e-9, 30e-9                      # ToF 前 10ns / 后 30ns
dt_fine  = PARAMS["hist"]["dt_fine"]          # 10 ps 精细网格
bin_width = PARAMS["hist"]["bin_width"]       # 1 ns 直方图 bin
t_lo, t_hi = t0 - pre, t0 + post
tf = np.arange(t_lo, t_hi, dt_fine)

PDE = PARAMS["spad"]["PDE"]; jit = PARAMS["spad"]["jitter_sigma"]
t_dead = 14e-9                                # ★ v12: SPAD dead time = 14 ns(不再是 0)

# ---- 中心 SPAD 的信号光子率(不含 PDE) ----
f_ij = fpix0[i0, j0]
r_sig = signal_photon_rate_fine(echo0, f_ij, tf)
tc = None  # 占位, 下面定义 centers 后赋值

# ---- (1) 单次 shot 的探测 timestamp(含 dead time + IRF 抖动) ----
rng_demo = np.random.default_rng(PARAMS["hist"]["seed"])
ev1 = simulate_spad_shot(r_sig, r_amb_ph, tf, PDE, t_dead, jit, rng_demo)
ev1_ns = np.sort(ev1) * 1e9
centers, c1, edges = events_to_hist(ev1, t_lo, t_hi, bin_width)
tc_ns = centers * 1e9; t0_ns = t0 * 1e9

print("="*76)
print(f"单 SPAD (i0,j0)=({i0},{j0}), f_pix={f_ij:.3e}; dead time={t_dead*1e9:.0f} ns; 30m ToF={t0_ns:.1f} ns")
N_sig_inc = float(np.trapezoid(r_sig, tf)); N_amb_inc = r_amb_ph*(t_hi-t_lo)
print(f"单 shot 入射(期望): 信号≈{N_sig_inc:.3f} ph, 环境≈{N_amb_inc:.3f} ph "
      f"(注: 实际每 shot 入射数 ~Poisson, 会涨落)")
print(f"本次单 shot 探测到 {ev1.size} 个光子, timestamp [ns]:")
print("  " + (", ".join(f"{t:.3f}" for t in ev1_ns) if ev1.size else "(本次无探测)"))

# ---- (2) 多次 shot 的 timestamp(展示 shot 间随机性) ----
N_raster = 30
rng_ras = np.random.default_rng(PARAMS["hist"]["seed"] + 100)
ras = [np.sort(simulate_spad_shot(r_sig, r_amb_ph, tf, PDE, t_dead, jit, rng_ras))*1e9
       for _ in range(N_raster)]
n_ev_ras = [len(e) for e in ras]
print(f"{N_raster} 次 shot 各自探测光子数: min={min(n_ev_ras)}, max={max(n_ev_ras)}, "
      f"均值={np.mean(n_ev_ras):.2f} (体现散粒噪声下的 shot 间涨落)")

# ---- (3) 多 shot 蒙卡均值(含 dead time) vs 解析 λ(dead time=0 参考) ----
Nrep = 3000
rng_val = np.random.default_rng(PARAMS["hist"]["seed"] + 1)
acc = np.zeros(len(centers))
for _ in range(Nrep):
    ev = simulate_spad_shot(r_sig, r_amb_ph, tf, PDE, t_dead, jit, rng_val)
    cc, _ = np.histogram(ev, bins=edges)
    acc += cc
mc_mean = acc / Nrep                          # 每 shot 平均计数/bin(含 14ns dead time)

# 解析 λ = (r_sig+r_amb)·PDE ⊗ IRF, dead time=0(理论上限参考)
r_det = (r_sig + r_amb_ph) * PDE
r_det = np.convolve(r_det, gaussian_kernel(jit, dt_fine), mode="same") * dt_fine
bin_idx = np.clip(((tf - t_lo)/bin_width).astype(int), 0, len(centers)-1)
lam = np.bincount(bin_idx, weights=r_det*dt_fine, minlength=len(centers))

tot_mc, tot_la = mc_mean.sum(), lam.sum()
pk = int(lam.argmax())
print(f"峰值 bin @ {tc_ns[pk]:.0f} ns: 蒙卡(含 {t_dead*1e9:.0f}ns dead time)={mc_mean[pk]:.4f} "
      f"vs 解析(dead time=0)={lam[pk]:.4f} (峰值区 pile-up 最强 -> 蒙卡偏低)")
print(f"窗口每 shot 总计数: 蒙卡={tot_mc:.4f} vs 解析(dead time=0)={tot_la:.4f} "
      f"-> dead time 计数损失 {100*(tot_la-tot_mc)/max(tot_la,1e-12):.2f}%")

# ---- 绘图 1: 探测 timestamp ----
fig, ax = plt.subplots(2, 1, figsize=(11, 5.6))
ax[0].eventplot([ev1_ns], lineoffsets=1, linelengths=0.8, colors="steelblue")
ax[0].axvline(t0_ns, color="r", ls=":", label=f"真实 ToF {t0_ns:.1f} ns")
ax[0].set_xlim(t_lo*1e9, t_hi*1e9); ax[0].set_yticks([]); ax[0].set_ylabel("单次 shot")
ax[0].set_title(f"单个 SPAD · 单次 shot 探测光子 timestamp (dead time={t_dead*1e9:.0f}ns, PDE={PDE}, 含 IRF 抖动)")
ax[0].legend(fontsize=9, loc="upper right"); ax[0].grid(alpha=0.3, axis="x")

ax[1].eventplot(ras, lineoffsets=np.arange(1, N_raster+1), linelengths=0.8, colors="steelblue")
ax[1].axvline(t0_ns, color="r", ls=":")
ax[1].set_xlim(t_lo*1e9, t_hi*1e9)
ax[1].set_xlabel("时间 t [ns]"); ax[1].set_ylabel("shot 序号")
ax[1].set_title(f"{N_raster} 次 shot 的探测 timestamp(每行一次 shot): 探测时刻随机、稀疏")
ax[1].grid(alpha=0.3, axis="x")
plt.tight_layout(); plt.show()

# ---- 绘图 2: 1ns 直方图(保留) + 验证 ----
fig, ax = plt.subplots(2, 1, figsize=(11, 7), sharex=True)
ax[0].bar(tc_ns, c1, width=bin_width*1e9, align="center", color="steelblue")
ax[0].axvline(t0_ns, color="r", ls=":", label=f"真实 ToF {t0_ns:.1f} ns")
ax[0].set_ylabel("计数 / 1ns bin")
ax[0].set_title(f"单个 SPAD · 单次 shot 计数-时间直方图 (1ns bin, dead time={t_dead*1e9:.0f}ns, seed={PARAMS['hist']['seed']})")
ax[0].legend(fontsize=9); ax[0].grid(alpha=0.3)

ax[1].bar(tc_ns, mc_mean, width=bin_width*1e9, align="center", color="lightsteelblue",
          label=f"蒙卡均值 ({Nrep} shots, dead time={t_dead*1e9:.0f}ns)")
ax[1].plot(tc_ns, lam, "r-", lw=1.8, marker="o", ms=3, label="解析 λ = PDE·r_ph 卷积 IRF (dead time=0 参考)")
ax[1].axvline(t0_ns, color="k", ls=":", alpha=0.6)
ax[1].set_xlabel("时间 t [ns]"); ax[1].set_ylabel("每 shot 平均计数 / bin")
ax[1].set_title("多 shot 蒙卡均值(含 dead time) vs 解析期望(dead time=0): 差值即 dead time 损失")
ax[1].legend(fontsize=9); ax[1].grid(alpha=0.3)
plt.tight_layout(); plt.show()''')

# ===== 模块 9(保留 v11 原逻辑, dead time 自动继承 14ns)=====
md(r'''## 模块 9（沿用）— 阵列内不同 SPAD 的响应差异（现含 dead time=14ns）

像斑是椭圆高斯 ⇒ **每个 SPAD 收集比例不同**, 计数率也不同(中心高、边缘低)。dead time 沿用模块 8 的 14ns。
- **左**: 像元收集比例 f_pix 热图(聚焦被照区域)。
- **右**: 沿长边 y 取中心/中段/边缘几个代表 SPAD, 各自多 shot 平均计数-时间对比。''')

code(r'''# 聚焦被照区域(收集比例 > 峰值 1% 的像元范围)做热图
thr = fpix0.max() * 0.01
xs = np.where(fx0 > fx0.max()*0.01)[0]; ys = np.where(fy0 > fy0.max()*0.01)[0]
ix = slice(max(0, xs.min()-1), min(fpix0.shape[0], xs.max()+2))
iy = slice(max(0, ys.min()-1), min(fpix0.shape[1], ys.max()+2))
sub = fpix0[ix, iy]

# 沿长边 y 取代表像元(相对中心 j0 的偏移, 单位=像元)
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
        ev = simulate_spad_shot(r_sig_d, r_amb_ph, tf, PDE, t_dead, jit, rng9)
        cc, _ = np.histogram(ev, bins=edges)
        acc += cc
    ax[1].plot(tc_ns, acc/reps, lw=1.6, marker="o", ms=3,
               label=f"Δy={d} 像元 (f={f_ij_d:.2e})")
ax[1].axvline(t0_ns, color="k", ls=":", alpha=0.6, label=f"ToF {t0_ns:.1f} ns")
ax[1].set_xlabel("时间 t [ns]"); ax[1].set_ylabel(f"每 shot 平均计数 / bin ({reps} shots)")
ax[1].set_title("不同 SPAD(沿长边 y 偏移)的计数-时间: 越偏离中心收集越少")
ax[1].legend(fontsize=8); ax[1].grid(alpha=0.3)
plt.tight_layout(); plt.show()

print("说明: y 方向 σ 很大(像斑长边≈800µm), 故沿 y 偏移几十个像元收集比例才明显下降; "
      "这正体现'对每一个 SPAD 单独做时间响应推导'。")''')

# ===== 模块 10(替换: 去图, 只留文字与数值)=====
md(r'''## 模块 10（v12 更新）— Macro Pixel（宏像元）定义 9×3（仅文字说明与数值，不绘图）

**定义**: 短边 x 的 9 个 SPAD 全取, 长边 y 每 3 个为一组。
阵列 `Nx×Ny = 9×120` ⇒ 沿长边(y)自上而下共 `120/3 = 40` 个宏像元, 每个含 `9×3 = 27` 个 SPAD。
- 宏像元 `m` (m=0 在顶部) 覆盖: 全部 x, y 索引 `[3m, 3m+3)`。
- 像斑为椭圆高斯(中心≈阵列中心, 长边 y 方向 σ≈200µm ≈ 20 像元) ⇒ 中间的宏像元收集信号多、
  两端的宏像元几乎只有环境光底噪。

> 本版按要求**不再绘制**收集比例图, 仅保留定义、数值输出与后续模块所需的量(`n_macro`, `macro_fsum`, `m_peak`)。''')

code(r'''# ---- 宏像元定义: 9(短边x全部) × 3(长边y) ----
Bx_m, By_m = 9, 3
a = PARAMS["spad_array"]
assert a["Nx"] == Bx_m, f"短边 SPAD 数={a['Nx']} 应等于 Bx_m={Bx_m}"
assert a["Ny"] % By_m == 0, f"长边 {a['Ny']} 不能被 {By_m} 整除"
n_macro = a["Ny"] // By_m                       # = 40
n_pix_macro = Bx_m * By_m                        # = 27

# 每个宏像元的信号空间收集比例(其 27 个 SPAD 的 f_pix 之和), 用主目标(30m)的像斑
macro_fsum = np.array([fpix0[:, m*By_m:(m+1)*By_m].sum() for m in range(n_macro)])
m_peak = int(macro_fsum.argmax())               # 收集最强的宏像元(≈像斑中心所在)

print(f"宏像元 = {Bx_m}(短边x全部) × {By_m}(长边y);  阵列 {a['Nx']}×{a['Ny']} → 共 {n_macro} 个宏像元, 每个 {n_pix_macro} 个 SPAD")
print(f"像斑中心 SPAD 在 y={j0} → 落在宏像元 m={j0//By_m}; 收集最强宏像元 m_peak={m_peak} (Σf_pix={macro_fsum[m_peak]:.3f})")
print(f"全部宏像元 Σf_pix 合计 = {macro_fsum.sum():.3f} (应≈整像斑落片比例 {fpix0.sum():.3f})")
print(f"边缘宏像元 m=0 的 Σf_pix = {macro_fsum[0]:.3e} (≈0, 近乎纯环境光底噪)")''')

# ===== 模块 11(替换: dead time=14ns, 校验改为损失量化)=====
md(r'''## 模块 11（v12 更新）— 逐-SPAD 蒙卡统计每个宏像元的直方图（dead time=14ns）

对每个宏像元, 让其 **27 个 SPAD 各自独立**跑蒙卡引擎(逐光子 → PDE 判定 → dead time=14ns → IRF 抖动),
把 27 个 SPAD 的事件**合并**成该宏像元的事件流, 再直方图化(1ns bin), 累加 `N_shots` 次。

- 复用模块 8 的时间窗、精细网格与 `t_dead`(=14ns)。
- **每个 SPAD 各自独立死时间**(合并前判定), 这正是逐-SPAD 结构的意义。
- 解析 λ(=`PDE·(base·Σf + 27·r_amb) ⊗ IRF`, **dead time=0**)作为**理论上限参考**:
  蒙卡因 dead time 的 pile-up 应**不高于**该参考; 二者差值即 dead time 计数损失。''')

code(r'''# 复用模块 8 定义的: tf, t_lo, t_hi, bin_width, centers, edges, tc_ns, t0_ns, dt_fine, PDE, jit, t_dead, r_amb_ph
nbins = len(centers)
base_rate = signal_photon_rate_fine(echo0, 1.0, tf)     # f_pix=1 的信号率形状(30m, tilt=0 无展宽)
N_shots = PARAMS["hist"]["N_shots"]

# 预存每个宏像元 27 个 SPAD 的空间收集比例
macro_fvals = [fpix0[:, m*By_m:(m+1)*By_m].ravel() for m in range(n_macro)]   # 每个: 长度27

# ---- 逐-SPAD 蒙卡(每 SPAD 各自 14ns dead time), 累加 N_shots ----
rng = np.random.default_rng(PARAMS["hist"]["seed"] + 11)
macro_hist = np.zeros((n_macro, nbins))
for _shot in range(N_shots):
    for m in range(n_macro):
        ev_all = []
        for fij in macro_fvals[m]:
            ev = simulate_spad_shot(base_rate * fij, r_amb_ph, tf, PDE, t_dead, jit, rng)
            if ev.size:
                ev_all.append(ev)
        if ev_all:
            cc, _ = np.histogram(np.concatenate(ev_all), bins=edges)
            macro_hist[m] += cc

# ---- 解析上限参考(dead time=0, 累加 N_shots): λ = PDE·(base·Σf + 27·r_amb) ⊗ IRF ----
irf_k = gaussian_kernel(jit, dt_fine)
bin_idx = np.clip(((tf - t_lo) / bin_width).astype(int), 0, nbins - 1)
macro_lam = np.zeros((n_macro, nbins))
for m in range(n_macro):
    r_det = (base_rate * macro_fsum[m] + n_pix_macro * r_amb_ph) * PDE
    r_det = np.convolve(r_det, irf_k, mode="same") * dt_fine
    macro_lam[m] = N_shots * np.bincount(bin_idx, weights=r_det * dt_fine, minlength=nbins)

# ---- 校验 / dead time 损失量化 ----
tot_mc, tot_la = macro_hist.sum(), macro_lam.sum()
pk_mc = macro_hist[m_peak].max(); pk_la = macro_lam[m_peak].max()
print("="*76)
print(f"逐-SPAD 蒙卡完成: {n_macro} 宏像元 × {n_pix_macro} SPAD × {N_shots} shots "
      f"= {n_macro*n_pix_macro*N_shots} 次单-SPAD 仿真 (dead time={t_dead*1e9:.0f}ns)")
print(f"时间窗 [{t_lo*1e9:.1f}, {t_hi*1e9:.1f}] ns, 30m ToF={t0_ns:.1f} ns, bin=1ns, {nbins} bins")
print(f"全体总计数: 蒙卡(含 dead time)={tot_mc:.0f} vs 解析(dead time=0 上限)={tot_la:.1f} "
      f"-> dead time 计数损失 {100*(tot_la-tot_mc)/max(tot_la,1e-9):.2f}%")
print(f"峰值宏像元 m={m_peak} 峰值 bin: 蒙卡={pk_mc:.0f} vs 解析(dead time=0)={pk_la:.1f} "
      f"(峰值区 pile-up 最强 -> 损失 {100*(pk_la-pk_mc)/max(pk_la,1e-9):.1f}%)")
print(f"信号最强宏像元 m={m_peak}: 峰值 bin={macro_hist[m_peak].max():.0f}, 总计数={macro_hist[m_peak].sum():.0f}; "
      f"边缘宏像元 m=0: 总计数={macro_hist[0].sum():.0f}(近乎纯底噪)")''')

# ===== 模块 12(替换: 删图B左, 加 timestamp)=====
md(r'''## 模块 12（v12 更新）— 宏像元-时间热图 + 峰值宏像元 27 个 SPAD 的 timestamp

- **图 A**: 40 个宏像元 × 时间(1ns bin)的**计数热图**(含 14ns dead time)。
- **图 B**: 峰值宏像元与几个偏离宏像元的**直方图**(蒙卡柱 + 解析线, 线为 dead time=0 上限参考)。
- **图 C（v12 新增, 替代原总计数柱状图）**: 峰值宏像元 `m_peak` 在**单次 shot** 里,
  其 **27 个 SPAD 各自探测到的光子 timestamp**(每行一个 SPAD)——直观展示"宏像元 = 27 个逐-SPAD 事件合并"。''')

code(r'''# ---- 图 A: 宏像元-时间 热图 ----
fig, ax = plt.subplots(figsize=(11, 5))
im = ax.imshow(macro_hist, origin="upper", aspect="auto", cmap="inferno",
               extent=[tc_ns[0]-0.5, tc_ns[-1]+0.5, n_macro-0.5, -0.5])
ax.axvline(t0_ns, color="cyan", ls=":", lw=1.2, label=f"真实 ToF {t0_ns:.1f} ns")
ax.axhline(m_peak, color="lime", ls=":", lw=1.0, alpha=0.7, label=f"峰值宏像元 m={m_peak}")
ax.set_xlabel("时间 t [ns]"); ax.set_ylabel("宏像元序号 m (0=顶部, 沿长边 y)")
ax.set_title(f"宏像元(9×3)直方图热图: 每宏像元 27 个 SPAD 逐-蒙卡求和 (N_shots={N_shots}, PDE={PDE}, dead time={t_dead*1e9:.0f}ns)")
ax.legend(fontsize=9, loc="upper right"); plt.colorbar(im, ax=ax, label="计数 / (宏像元, 1ns bin)")
plt.tight_layout(); plt.show()

# ---- 图 B: 代表性宏像元直方图(柱=蒙卡含 dead time, 线=解析 dead time=0 参考) ----
fig, ax = plt.subplots(figsize=(11, 4.6))
reps_m = sorted(set([0, max(0, m_peak-6), m_peak, min(n_macro-1, m_peak+6)]))
colors = ["tab:gray", "tab:green", "tab:red", "tab:orange", "tab:purple"]
for m, c in zip(reps_m, colors):
    ax.bar(tc_ns, macro_hist[m], width=bin_width*1e9, align="center", alpha=0.35, color=c)
    ax.plot(tc_ns, macro_lam[m], color=c, lw=1.6, label=f"m={m} (Σf={macro_fsum[m]:.3f})")
ax.axvline(t0_ns, color="k", ls=":", alpha=0.6, label=f"ToF {t0_ns:.1f} ns")
ax.set_xlabel("时间 t [ns]"); ax.set_ylabel(f"计数 / 1ns bin ({N_shots} shots)")
ax.set_title("代表性宏像元直方图(柱=蒙卡含 dead time, 线=解析 dead time=0 参考)")
ax.legend(fontsize=8); ax.grid(alpha=0.3)
plt.tight_layout(); plt.show()

# ---- 图 C(v12 新增): 峰值宏像元单次 shot, 27 个 SPAD 各自 timestamp ----
rng_ts = np.random.default_rng(PARAMS["hist"]["seed"] + 200)
sp_events = []
for fij in macro_fvals[m_peak]:
    ev = simulate_spad_shot(base_rate * fij, r_amb_ph, tf, PDE, t_dead, jit, rng_ts)
    sp_events.append(np.sort(ev) * 1e9)
n_ev_total = sum(len(e) for e in sp_events)
n_ev_each = [len(e) for e in sp_events]

fig, ax = plt.subplots(figsize=(11, 5))
ax.eventplot(sp_events, lineoffsets=np.arange(1, n_pix_macro+1), linelengths=0.8, colors="steelblue")
ax.axvline(t0_ns, color="r", ls=":", label=f"真实 ToF {t0_ns:.1f} ns")
ax.set_xlim(t_lo*1e9, t_hi*1e9); ax.set_ylim(0.5, n_pix_macro+0.5)
ax.set_xlabel("时间 t [ns]"); ax.set_ylabel("宏像元内 SPAD 序号 (1..27)")
ax.set_title(f"峰值宏像元 m={m_peak} · 单次 shot · 27 个 SPAD 各自探测 timestamp (dead time={t_dead*1e9:.0f}ns)")
ax.legend(fontsize=9); ax.grid(alpha=0.3, axis="x")
plt.tight_layout(); plt.show()

print(f"峰值宏像元 m={m_peak} 单次 shot: 27 个 SPAD 共探测 {n_ev_total} 个光子 "
      f"(各 SPAD: min={min(n_ev_each)}, max={max(n_ev_each)})。")
print(f"每个 SPAD 各自 14ns dead time, 合并成宏像元事件流后再直方图化。")
print(f"底噪水平: 边缘宏像元 m=0 总计数≈{macro_hist[0].sum():.0f} ({N_shots} shots, 27 SPAD × 环境光)。")''')

# --- 组装 & 写出 ---
for typ, src in NEW:
    kept.append(nbformat.v4.new_markdown_cell(src) if typ == "md"
                else nbformat.v4.new_code_cell(src))
nb.cells = kept
nbformat.write(nb, DST)
print(f"已生成 {DST}: 共 {len(kept)} 个 cell (保留 v11 模块0-7 + 新增 {len(NEW)})。")
