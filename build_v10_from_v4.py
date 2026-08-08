# -*- coding: utf-8 -*-
"""
基于 v4 生成 v10:
- 原样保留 v4 的"模块 0-5"(参数 + 全部物理链路函数), 保证物理一致、可与 v4 对照;
- 追加 v10 新模块 6-9: 单 SPAD 蒙特卡洛(逐光子, dead time=0, 预留 dead time 接口)。
运行: python build_v10_from_v4.py
"""
import nbformat

SRC = "lidar_histogram_sim_v4.ipynb"
DST = "lidar_histogram_sim_v10.ipynb"

nb = nbformat.read(SRC, as_version=4)

# --- 保留 v4 中直到(不含)"## 模块 6"的所有 cell ---
kept = []
for c in nb.cells:
    if c.cell_type == "markdown" and c.source.lstrip().startswith("## 模块 6"):
        break
    kept.append(c)

# 清空保留 cell 的输出/执行序号, 让文件干净(用户会 Restart & Run All)
for c in kept:
    if c.cell_type == "code":
        c.outputs = []
        c.execution_count = None

# 把第一个标题 cell(v4 大标题)替换为 v10 标题
if kept and kept[0].cell_type == "markdown":
    kept[0].source = MD_TITLE = r'''# 激光雷达直方图仿真 v10 (LiDAR Histogram Simulation) — 单 SPAD 蒙特卡洛

链路(激光→TX→信道→目标→RX)与 **v4 完全一致**(模块 0–5 原样复用, 便于对照)。
v10 的新东西在**模块 6–9**: 把"求期望 λ 再取 Poisson"的做法换成**逐光子蒙特卡洛**
(Monte Carlo, 蒙卡), 面向**单个 SPAD**(Single-Photon Avalanche Diode, 单光子雪崩二极管)
做时间上的计数响应推导。

**v10 相对 v4 的变化**
- 新增单 SPAD 光子到达率(不含 PDE, Photon Detection Efficiency, 光子探测效率)的空间分配。
- 新增逐光子蒙卡探测引擎: 入射光子 ~Poisson → 每光子以 PDE 概率触发雪崩 → **n 个光子全未触发则无计数**。
- IRF(Instrument Response Function, 仪器响应函数)高斯抖动逐事件叠加在时间戳上。
- **预留 dead time(死时间)接口**: 当前 t_dead=0 不生效; 探测一次后 t_dead 内不再响应, 将来直接改非零即可。
- 仿真只在 30m 目标 ToF(Time of Flight, 飞行时间)附近的时间窗内进行(省算力)。

**尚未实装(按约定留待后续)**
1. 自上而下每 3×9 为一个 macro pixel(宏像元), 共 40 个, 统计各宏像元内 SPAD 的 1ns bin 直方图。
2. dead time 的实际影响(接口已预留, 现在恒为 0)。

> ⚠️ 运行方式: **Kernel → Restart & Run All**(从上到下顺序执行)。'''

# ---------------- v10 新增 cells ----------------
NEW_CELLS = []
def md(s): NEW_CELLS.append(("md", s))
def code(s): NEW_CELLS.append(("code", s))

# ===== 模块 6 =====
md(r'''## 模块 6（v10 新增）— 单 SPAD 光子到达率（不含 PDE）

把到达像面的回波**总功率**按**椭圆高斯像斑**分配到每个物理像元(用误差函数 erf 对像元窗口积分),
得到**每个 SPAD 各自的**信号光子到达率 r_ph(t)(单位 photons/s, **不含 PDE**)。
环境光按单像元均匀铺底, 也给出不含 PDE 的光子率。

- 像斑 1/e² 全宽 s(来自模块 5) → 高斯 σ = s/4; 像元窗口 [c−pitch/2, c+pitch/2] 内高斯积分 = erf 差分。
- 像斑中心默认对准阵列中心; 每个像元收集比例不同(中心高、边缘低) ⇒ 每个 SPAD 的率不同。
- 倾角展宽(几何真实, 非 IRF)并入信号率(对 30m、tilt=0 无影响, 代码通用保留)。

> 关键: 这里用**未乘 PDE** 的光子率; PDE 留到模块 7 蒙卡里逐光子判定。''')

code(r'''from scipy.special import erf

def gaussian_kernel(sigma, dt, n_sigma=5):
    """归一化高斯卷积核(∫k·dt=1); 供 IRF/倾角展宽卷积用(与 v4 定义一致)。"""
    if sigma <= 0:
        return np.array([1.0 / dt])
    half = max(1, int(np.ceil(n_sigma * sigma / dt)))
    tk = np.arange(-half, half + 1) * dt
    k = np.exp(-0.5 * (tk / sigma)**2)
    return k / (k.sum() * dt)

def ambient_photon_rate_per_pixel(p=PARAMS):
    """单像元环境光【光子到达率】[photons/s], 不含 PDE
    (与 v4 的 ambient_count_rate_per_pixel 相差一个 PDE 因子)。"""
    if not p["ambient"]["enable"]:
        return 0.0
    E = p["ambient"]["E_lambda"] * (p["rx"]["filter_bw"] * 1e9)   # 带内辐照 [W/m²]
    L = p["ambient"]["surface_rho"] * E / np.pi                   # 辐亮度 [W/m²/sr]
    iFOV = p["spad_array"]["pitch"] / p["rx"]["f_RX"]
    P_amb = L * iFOV**2 * rx_area(p)                              # 单像元带内功率 [W]
    return P_amb / E_PHOTON * p["rx"]["T_RX"] * p["rx"]["T_filter"]   # 不含 PDE

def pixel_grid(p=PARAMS):
    """阵列各像元中心坐标 [m](以阵列中心为原点); 像斑中心默认对准阵列中心。"""
    a = p["spad_array"]; pitch = a["pitch"]
    xi = (np.arange(a["Nx"]) - (a["Nx"] - 1) / 2.0) * pitch   # x=短边
    yj = (np.arange(a["Ny"]) - (a["Ny"] - 1) / 2.0) * pitch   # y=长边
    return xi, yj

def pixel_collection_matrix(D, p=PARAMS):
    """每个像元在椭圆高斯像斑上的空间收集比例 f_pix[i,j] (∑≤1, 其余漏到阵列外)。
    返回 (f_pix[Nx,Ny], fx[Nx], fy[Ny])。"""
    sx, sy = rx_image_spot_size(D, p)          # 1/e² 全宽 (x=短边, y=长边)
    sig_x, sig_y = sx / 4.0, sy / 4.0          # 高斯 σ
    xi, yj = pixel_grid(p); pitch = p["spad_array"]["pitch"]
    def _frac(centers, sig):                   # 各像元窗口内一维高斯积分(erf 差分)
        lo = (centers - pitch / 2.0) / (np.sqrt(2) * sig)
        hi = (centers + pitch / 2.0) / (np.sqrt(2) * sig)
        return 0.5 * (erf(hi) - erf(lo))
    fx = _frac(xi, sig_x); fy = _frac(yj, sig_y)
    return np.outer(fx, fy), fx, fy

def signal_photon_rate_fine(echo, f_pix_ij, tf, p=PARAMS):
    """单 SPAD 信号【光子到达率】(不含 PDE)在精细网格 tf 上。
    已含倾角几何展宽, 未含 IRF; f_pix_ij 为该像元空间收集比例(标量)。"""
    t0 = time_of_flight(echo["D"])
    r = pulse_temporal(tf - t0, p) * link_factor(echo, p) / E_PHOTON * f_pix_ij
    sig_b = echo_range_broadening_sigma(echo["D"], echo["tilt_deg"], p)   # 几何展宽(非 IRF)
    if sig_b > 0:
        dt_fine = tf[1] - tf[0]
        r = np.convolve(r, gaussian_kernel(sig_b, dt_fine), mode="same") * dt_fine
    return r

# --- 自检: 主目标(frac 最大, 本例 30m) 的像元收集分布 ---
echo0 = max(PARAMS["target"]["echoes"], key=lambda e: e["frac"])
D0 = echo0["D"]
fpix0, fx0, fy0 = pixel_collection_matrix(D0)
i0, j0 = int(np.argmax(fx0)), int(np.argmax(fy0))     # 收集最强的像元(≈阵列中心)
sx0, sy0 = rx_image_spot_size(D0)
r_amb_ph = ambient_photon_rate_per_pixel()
print(f"主目标 D={D0:.0f} m  (frac={echo0['frac']}, ρ={echo0['rho']}, tilt={echo0['tilt_deg']}°)")
print(f"像斑 1/e² 全宽: x(短边)={sx0*1e6:.1f} µm, y(长边)={sy0*1e6:.1f} µm  → σ_x={sx0/4*1e6:.1f}, σ_y={sy0/4*1e6:.1f} µm")
print(f"阵列 {PARAMS['spad_array']['Nx']}×{PARAMS['spad_array']['Ny']}, pitch={PARAMS['spad_array']['pitch']*1e6:.0f} µm; "
      f"像斑总落片比例 ∑f_pix={fpix0.sum():.3f}")
print(f"最强像元 (i0,j0)=({i0},{j0}): f_pix={fpix0[i0,j0]:.3e}  (fx={fx0[i0]:.3f}, fy={fy0[j0]:.3f})")
print(f"单像元环境光光子率(不含 PDE) r_amb = {r_amb_ph:.3e} ph/s "
      f"(×PDE={PARAMS['spad']['PDE']} → 计数率 {r_amb_ph*PARAMS['spad']['PDE']:.3e} cps, 应与 v4 一致)")''')

# ===== 模块 7 =====
md(r'''## 模块 7（v10 新增）— 单 SPAD 蒙特卡洛探测引擎（逐光子，预留 dead time）

一次 shot 的处理:
1. 每精细 bin 入射光子数 ~ **Poisson(r_ph·dt)**, 展开成逐光子到达时刻(=精细 bin 中心)。
2. **有限 PDE**: 每个入射光子独立以概率 PDE 触发雪崩; **n 个光子全未触发 ⇒ 无计数**(蒙卡本质)。
3. **IRF** 高斯抖动只加在被记录事件的时间戳上(雪崩时刻本身留给 dead time 判定)。
4. **dead time 已预留**(当前 t_dead=0 不生效): 探测一次后 t_dead 内不再响应。

返回该 SPAD 本次 shot 的事件时间戳数组 → 可 rebin 成 1ns 直方图。''')

code(r'''def simulate_spad_shot(r_sig_fine, r_amb_ph, tf, PDE, t_dead, jitter_sigma, rng):
    """单 SPAD 单次 shot 蒙卡, 返回被记录事件时间戳数组 [s]。
    r_sig_fine: 信号光子率(不含 PDE)在精细网格 tf 上; r_amb_ph: 环境光光子率(标量, 不含 PDE)。
    dead time=0 时每光子独立(向量化); t_dead>0 走逐事件扫描(将来用)。"""
    dt = tf[1] - tf[0]
    mu = (r_sig_fine + r_amb_ph) * dt                 # 每精细 bin 期望入射光子数
    n_ph = rng.poisson(mu)                            # 每 bin 实际入射光子数
    if n_ph.sum() == 0:
        return np.empty(0)
    t_arr = np.repeat(tf, n_ph)                       # 逐光子到达时刻(bin 中心), 已升序
    fired = rng.random(t_arr.size) < PDE              # 每个入射光子是否触发雪崩(PDE)

    if t_dead <= 0:                                   # —— dead time=0: 每光子独立, 向量化 ——
        t_fire = t_arr[fired]
        if jitter_sigma > 0 and t_fire.size:
            t_fire = t_fire + rng.normal(0.0, jitter_sigma, t_fire.size)   # IRF 抖动(只影响记录时间戳)
        return t_fire

    # —— dead time>0: 逐事件扫描(接口预留, 现在不会走到) ——
    det = []; last = -np.inf
    for t, f in zip(t_arr, fired):
        if not f:
            continue
        if (t - last) < t_dead:                       # 死时间内: 吞掉, 不响应
            continue
        last = t                                      # 记录雪崩时刻(未加 IRF)
        det.append(t + (rng.normal(0.0, jitter_sigma) if jitter_sigma > 0 else 0.0))
    return np.asarray(det)

def events_to_hist(events, t_lo, t_hi, bin_width):
    """事件时间戳 → 直方图。返回 (bin 中心, 计数, bin 边界)。"""
    nb = int(round((t_hi - t_lo) / bin_width))
    edges = t_lo + np.arange(nb + 1) * bin_width
    counts, _ = np.histogram(events, bins=edges)
    centers = 0.5 * (edges[:-1] + edges[1:])
    return centers, counts, edges

print("蒙卡引擎就绪: simulate_spad_shot() 逐光子(PDE 判定 + IRF 抖动 + dead time 预留); events_to_hist() 直方图化。")''')

# ===== 模块 8 =====
md(r'''## 模块 8（v10 新增）— 单次 shot 演示 + 多 shot 验证蒙卡无偏

只在 30m 目标 ToF 附近的时间窗内仿真。选像斑中心 SPAD:
- **上图**: 单个 SPAD **单次 shot** 的计数-时间直方图 —— 即"某一次、某个 SPAD、计数随时间"。
- **下图**: 大量 shot 的**蒙卡均值** vs **解析期望 λ = PDE·r_ph ⊗ IRF**。二者吻合 ⇒ 蒙卡实现无偏。''')

code(r'''# ---- 时间窗: 只在 30m 目标 ToF 附近 ----
t0 = time_of_flight(D0)
pre, post = 10e-9, 30e-9                      # ToF 前 10ns / 后 30ns(双指数尾 + 裕量)
dt_fine  = PARAMS["hist"]["dt_fine"]          # 10 ps 精细网格
bin_width = PARAMS["hist"]["bin_width"]       # 1 ns 直方图 bin
t_lo, t_hi = t0 - pre, t0 + post
tf = np.arange(t_lo, t_hi, dt_fine)

PDE = PARAMS["spad"]["PDE"]; jit = PARAMS["spad"]["jitter_sigma"]; t_dead = 0.0

# ---- 中心 SPAD 的信号光子率(不含 PDE) ----
f_ij = fpix0[i0, j0]
r_sig = signal_photon_rate_fine(echo0, f_ij, tf)

# ---- (1) 单次 shot ----
rng_demo = np.random.default_rng(PARAMS["hist"]["seed"])
ev1 = simulate_spad_shot(r_sig, r_amb_ph, tf, PDE, t_dead, jit, rng_demo)
centers, c1, edges = events_to_hist(ev1, t_lo, t_hi, bin_width)

# ---- (2) 多 shot 蒙卡均值 ----
Nrep = 3000
rng_val = np.random.default_rng(PARAMS["hist"]["seed"] + 1)
acc = np.zeros(len(centers))
for _ in range(Nrep):
    ev = simulate_spad_shot(r_sig, r_amb_ph, tf, PDE, t_dead, jit, rng_val)
    cc, _ = np.histogram(ev, bins=edges)
    acc += cc
mc_mean = acc / Nrep                          # 每 shot 平均计数/bin

# ---- (3) 解析期望 λ = (r_sig+r_amb)·PDE ⊗ IRF, rebin 到 1ns(每 shot) ----
r_det = (r_sig + r_amb_ph) * PDE
r_det = np.convolve(r_det, gaussian_kernel(jit, dt_fine), mode="same") * dt_fine
bin_idx = np.clip(((tf - t_lo) / bin_width).astype(int), 0, len(centers) - 1)
lam = np.bincount(bin_idx, weights=r_det * dt_fine, minlength=len(centers))

# ---- 数字汇总 ----
tc_ns = centers * 1e9; t0_ns = t0 * 1e9
N_sig_inc = float(np.trapezoid(r_sig, tf))               # 单 shot 入射信号光子(到该像元)
N_amb_inc = r_amb_ph * (t_hi - t_lo)                     # 单 shot 入射环境光子(窗口内)
print("="*76)
print(f"单 SPAD (i0,j0)=({i0},{j0}), f_pix={f_ij:.3e}; 时间窗 [{t_lo*1e9:.1f}, {t_hi*1e9:.1f}] ns, 30m ToF={t0_ns:.1f} ns")
print(f"峰值信号光子率(不含 PDE)={r_sig.max():.3e} ph/s;  环境光子率={r_amb_ph:.3e} ph/s")
print(f"单 shot 入射: 信号≈{N_sig_inc:.2f} ph, 环境≈{N_amb_inc:.3f} ph  → 探测≈信号{N_sig_inc*PDE:.2f}+环境{N_amb_inc*PDE:.3f} ph")
print(f"本次单 shot 记录事件数={ev1.size};  峰值 bin 计数={c1.max()}")
peak = lam.argmax()
print(f"验证(峰值 bin @ {tc_ns[peak]:.0f} ns): 蒙卡均值={mc_mean[peak]:.4f} vs 解析 λ={lam[peak]:.4f} "
      f"(相对差 {100*abs(mc_mean[peak]-lam[peak])/max(lam[peak],1e-12):.1f}%)")
tot_mc, tot_la = mc_mean.sum(), lam.sum()
print(f"验证(窗口内每 shot 总计数): 蒙卡={tot_mc:.4f} vs 解析={tot_la:.4f} (相对差 {100*abs(tot_mc-tot_la)/tot_la:.2f}%)")

# ---- 绘图 ----
fig, ax = plt.subplots(2, 1, figsize=(11, 7), sharex=True)
ax[0].bar(tc_ns, c1, width=bin_width*1e9, align="center", color="steelblue")
ax[0].axvline(t0_ns, color="r", ls=":", label=f"真实 ToF {t0_ns:.1f} ns")
ax[0].set_ylabel("计数 / 1ns bin")
ax[0].set_title(f"单个 SPAD · 单次 shot 计数-时间 (30m, PDE={PDE}, dead time=0, seed={PARAMS['hist']['seed']})")
ax[0].legend(fontsize=9); ax[0].grid(alpha=0.3)

ax[1].bar(tc_ns, mc_mean, width=bin_width*1e9, align="center", color="lightsteelblue",
          label=f"蒙卡均值 ({Nrep} shots)")
ax[1].plot(tc_ns, lam, "r-", lw=1.8, marker="o", ms=3, label="解析期望 λ = PDE·r_ph 卷积 IRF")
ax[1].axvline(t0_ns, color="k", ls=":", alpha=0.6)
ax[1].set_xlabel("时间 t [ns]"); ax[1].set_ylabel("每 shot 平均计数 / bin")
ax[1].set_title("多 shot 蒙卡均值 vs 解析期望 (吻合 -> 蒙卡实现无偏)")
ax[1].legend(fontsize=9); ax[1].grid(alpha=0.3)
plt.tight_layout(); plt.show()''')

# ===== 模块 9 =====
md(r'''## 模块 9（v10 新增）— 阵列内不同 SPAD 的响应差异

像斑是椭圆高斯 ⇒ **每个 SPAD 收集比例不同**, 计数率也不同(中心高、边缘低)。
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

# --- 组装并写出 ---
for typ, src in NEW_CELLS:
    kept.append(nbformat.v4.new_markdown_cell(src) if typ == "md"
                else nbformat.v4.new_code_cell(src))

nb.cells = kept
nbformat.write(nb, DST)
print(f"已生成 {DST}: 共 {len(kept)} 个 cell (保留 v4 前置 + 新增 {len(NEW_CELLS)})。")
