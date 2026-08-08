# -*- coding: utf-8 -*-
"""
基于 v10 生成 v11:
- 保留 v10 全部模块 0-9(参数 + 物理链路 + 单 SPAD 蒙卡引擎), 清空输出;
- 追加 v11 新模块 10-12: Macro Pixel(宏像元) = 9(短边x全部)×3(长边y),
  沿长边共 40 个; 用逐-SPAD 蒙卡统计每个宏像元的 1ns 直方图(仅 30m 目标, ToF 附近)。
运行: python build_v11_from_v10.py
"""
import nbformat

SRC = "lidar_histogram_sim_v10.ipynb"
DST = "lidar_histogram_sim_v11.ipynb"

nb = nbformat.read(SRC, as_version=4)

# 清空所有代码 cell 的输出/执行序号(用户会 Restart & Run All)
for c in nb.cells:
    if c.cell_type == "code":
        c.outputs = []
        c.execution_count = None

# 替换首个标题 cell 为 v11 标题
if nb.cells and nb.cells[0].cell_type == "markdown":
    nb.cells[0].source = r'''# 激光雷达直方图仿真 v11 (LiDAR Histogram Simulation) — 宏像元 × 单 SPAD 蒙特卡洛

链路(激光→TX→信道→目标→RX)与 **v4 完全一致**; 单 SPAD(Single-Photon Avalanche Diode,
单光子雪崩二极管)逐光子蒙特卡洛(Monte Carlo, 蒙卡)引擎与 **v10 完全一致**(模块 0–9 原样复用)。

**v11 相对 v10 的变化(新模块 10–12)**
- 定义 **Macro Pixel(宏像元)= 9 × 3**: 短边 x 的 9 个 SPAD **全取**, 长边 y 每 **3** 个为一组;
  阵列 9×120 → 沿长边自上而下共 **120/3 = 40 个宏像元**, 每个含 9×3 = **27 个 SPAD**。
- 用**逐-SPAD 蒙卡**(复用 v10 引擎, dead time=0)统计**每个宏像元**的 **1ns bin** 直方图:
  宏像元直方图 = 其 27 个 SPAD 各自独立蒙卡事件之**和**。
- 仅针对 **30m 目标**, 只在其 ToF(Time of Flight, 飞行时间)附近的时间窗内仿真。
- 逐-SPAD 结构为将来 **dead time(死时间)** 预留: dead time=0 时等价于"率相加的单一泊松过程"(用作无偏校验);
  加 dead time 后每个 SPAD 各自死时间, 结构无需重写。

**尚未实装(按约定留待后续)**
- dead time 的实际影响(v10 引擎已预留接口, 现恒为 0)。

> ⚠️ 运行方式: **Kernel → Restart & Run All**(从上到下顺序执行)。
> 说明: 上一版口头总结曾把"3 的方向"说反, 现更正并锁定: **3 沿长边 y, 9 沿短边 x**。'''

# ---------------- v11 新增 cells ----------------
NEW = []
def md(s): NEW.append(("md", s))
def code(s): NEW.append(("code", s))

# ===== 模块 10 =====
md(r'''## 模块 10（v11 新增）— Macro Pixel（宏像元）定义 9×3 与收集比例分布

**定义**: 短边 x 的 9 个 SPAD 全取, 长边 y 每 3 个为一组。
阵列 `Nx×Ny = 9×120` ⇒ 沿长边(y)自上而下共 `120/3 = 40` 个宏像元, 每个含 `9×3 = 27` 个 SPAD。
- 宏像元 `m` (m=0 在顶部) 覆盖: 全部 x, y 索引 `[3m, 3m+3)`。
- 像斑为椭圆高斯(中心≈阵列中心, 长边 y 方向 σ≈200µm ≈ 20 像元) ⇒ 中间的宏像元收集信号多、
  两端的宏像元几乎只有环境光底噪。下图给出每个宏像元的**信号空间收集比例** `Σf_pix`(其 27 个 SPAD 之和)。''')

code(r'''# ---- 宏像元定义: 9(短边x全部) × 3(长边y) ----
Bx_m, By_m = 9, 3
a = PARAMS["spad_array"]
assert a["Nx"] == Bx_m, f"短边 SPAD 数={a['Nx']} 应等于 Bx_m={Bx_m}"
assert a["Ny"] % By_m == 0, f"长边 {a['Ny']} 不能被 {By_m} 整除"
n_macro = a["Ny"] // By_m                      # = 40
n_pix_macro = Bx_m * By_m                       # = 27

# 每个宏像元的信号空间收集比例(其 27 个 SPAD 的 f_pix 之和), 用主目标(30m)的像斑
macro_fsum = np.array([fpix0[:, m*By_m:(m+1)*By_m].sum() for m in range(n_macro)])
m_peak = int(macro_fsum.argmax())              # 收集最强的宏像元(≈像斑中心所在)

print(f"宏像元 = {Bx_m}(短边x全部) × {By_m}(长边y);  阵列 {a['Nx']}×{a['Ny']} → 共 {n_macro} 个宏像元, 每个 {n_pix_macro} 个 SPAD")
print(f"像斑中心 SPAD 在 y={j0} → 落在宏像元 m={j0//By_m}; 收集最强宏像元 m_peak={m_peak} (Σf_pix={macro_fsum[m_peak]:.3f})")
print(f"全部宏像元 Σf_pix 合计 = {macro_fsum.sum():.3f} (应≈整像斑落片比例 {fpix0.sum():.3f})")

fig, ax = plt.subplots(figsize=(10, 3.6))
ax.bar(np.arange(n_macro), macro_fsum, color="tab:blue")
ax.axvline(m_peak, color="r", ls=":", label=f"收集最强 m={m_peak}")
ax.set_xlabel("宏像元序号 m (0=顶部, 沿长边 y 自上而下)")
ax.set_ylabel("信号收集比例 Σf_pix")
ax.set_title(f"各宏像元(9×3)的信号空间收集比例: 中间强、两端弱 (30m 像斑长边 σ≈{sy0/4*1e6:.0f}µm)")
ax.legend(); ax.grid(alpha=0.3, axis="y")
plt.tight_layout(); plt.show()''')

# ===== 模块 11 =====
md(r'''## 模块 11（v11 新增）— 逐-SPAD 蒙卡统计每个宏像元的直方图（dead time=0）

对每个宏像元, 让其 **27 个 SPAD 各自独立**跑 v10 蒙卡引擎(逐光子 → PDE 判定 → IRF 抖动),
把 27 个 SPAD 的事件**合并**成该宏像元的事件流, 再直方图化(1ns bin), 累加 `N_shots` 次。

- 复用模块 8 的时间窗与网格(30m ToF 附近, `dt_fine`=10ps 精细网格, 1ns 直方图 bin)。
- 环境光对 27 个 SPAD 都存在 ⇒ 每个宏像元底噪 ∝ 27。
- **无偏校验**: dead time=0 时, 逐-SPAD 求和等价于率相加的单一泊松过程,
  解析期望 `λ_m = PDE·(base·Σf_pix + 27·r_amb) ⊗ IRF`。蒙卡累加应与之吻合。''')

code(r'''# 复用模块 8 定义的: tf, t_lo, t_hi, bin_width, centers, edges, tc_ns, t0_ns, dt_fine, PDE, jit, t_dead, r_amb_ph
nbins = len(centers)
base_rate = signal_photon_rate_fine(echo0, 1.0, tf)     # f_pix=1 的信号率形状(30m, tilt=0 无展宽)
N_shots = PARAMS["hist"]["N_shots"]

# 预存每个宏像元 27 个 SPAD 的空间收集比例
macro_fvals = [fpix0[:, m*By_m:(m+1)*By_m].ravel() for m in range(n_macro)]   # 每个: 长度27

# ---- 逐-SPAD 蒙卡, 累加 N_shots ----
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

# ---- 解析期望(每宏像元, 累加 N_shots): λ = PDE·(base·Σf + 27·r_amb) ⊗ IRF ----
irf_k = gaussian_kernel(jit, dt_fine)
bin_idx = np.clip(((tf - t_lo) / bin_width).astype(int), 0, nbins - 1)
macro_lam = np.zeros((n_macro, nbins))
for m in range(n_macro):
    r_det = (base_rate * macro_fsum[m] + n_pix_macro * r_amb_ph) * PDE
    r_det = np.convolve(r_det, irf_k, mode="same") * dt_fine
    macro_lam[m] = N_shots * np.bincount(bin_idx, weights=r_det * dt_fine, minlength=nbins)

# ---- 校验 ----
tot_mc, tot_la = macro_hist.sum(), macro_lam.sum()
pk_mc = macro_hist[m_peak].max(); pk_la = macro_lam[m_peak].max()
print("="*76)
print(f"逐-SPAD 蒙卡完成: {n_macro} 宏像元 × {n_pix_macro} SPAD × {N_shots} shots "
      f"= {n_macro*n_pix_macro*N_shots} 次单-SPAD 仿真")
print(f"时间窗 [{t_lo*1e9:.1f}, {t_hi*1e9:.1f}] ns, 30m ToF={t0_ns:.1f} ns, bin=1ns, {nbins} bins")
print(f"校验(全体总计数): 蒙卡={tot_mc:.0f} vs 解析={tot_la:.1f} (相对差 {100*abs(tot_mc-tot_la)/tot_la:.2f}%)")
print(f"校验(峰值宏像元 m={m_peak} 峰值 bin): 蒙卡={pk_mc:.0f} vs 解析={pk_la:.1f} "
      f"(相对差 {100*abs(pk_mc-pk_la)/max(pk_la,1e-9):.1f}%)")
print(f"信号最强宏像元 m={m_peak}: 峰值 bin 计数={macro_hist[m_peak].max():.0f}, "
      f"总计数={macro_hist[m_peak].sum():.0f}; 边缘宏像元 m=0: 总计数={macro_hist[0].sum():.0f}(近乎纯底噪)")''')

# ===== 模块 12 =====
md(r'''## 模块 12（v11 新增）— 可视化：宏像元-时间热图 + 代表性宏像元直方图

- **图 A**: 40 个宏像元 × 时间(1ns bin)的**计数热图** —— 直观展示信号集中在中间宏像元、ToF≈200ns 处。
- **图 B(左)**: 每个宏像元的**总计数**随宏像元序号变化(信号+底噪), 峰在像斑中心所在宏像元。
- **图 B(右)**: 峰值宏像元与几个偏离宏像元的**直方图**(蒙卡柱 + 解析线), 展示"越偏离中心信号越弱"。''')

code(r'''# ---- 图 A: 宏像元-时间 热图 ----
fig, ax = plt.subplots(figsize=(11, 5))
im = ax.imshow(macro_hist, origin="upper", aspect="auto", cmap="inferno",
               extent=[tc_ns[0]-0.5, tc_ns[-1]+0.5, n_macro-0.5, -0.5])
ax.axvline(t0_ns, color="cyan", ls=":", lw=1.2, label=f"真实 ToF {t0_ns:.1f} ns")
ax.axhline(m_peak, color="lime", ls=":", lw=1.0, alpha=0.7, label=f"峰值宏像元 m={m_peak}")
ax.set_xlabel("时间 t [ns]"); ax.set_ylabel("宏像元序号 m (0=顶部, 沿长边 y)")
ax.set_title(f"宏像元(9×3)直方图热图: 每宏像元 27 个 SPAD 逐-蒙卡求和 (N_shots={N_shots}, PDE={PDE}, dead time=0)")
ax.legend(fontsize=9, loc="upper right"); plt.colorbar(im, ax=ax, label="计数 / (宏像元, 1ns bin)")
plt.tight_layout(); plt.show()

# ---- 图 B: 总计数分布 + 代表性宏像元直方图 ----
fig, ax = plt.subplots(1, 2, figsize=(13, 4.6))
macro_tot = macro_hist.sum(axis=1)
ax[0].bar(np.arange(n_macro), macro_tot, color="tab:blue")
ax[0].axvline(m_peak, color="r", ls=":", label=f"峰值 m={m_peak}")
ax[0].set_xlabel("宏像元序号 m (0=顶部)"); ax[0].set_ylabel(f"总计数 ({N_shots} shots)")
ax[0].set_title("各宏像元总计数(信号+底噪): 中间强、两端仅底噪")
ax[0].legend(); ax[0].grid(alpha=0.3, axis="y")

reps = sorted(set([0, max(0, m_peak-6), m_peak, min(n_macro-1, m_peak+6)]))
colors = ["tab:gray", "tab:green", "tab:red", "tab:orange", "tab:purple"]
for m, c in zip(reps, colors):
    ax[1].bar(tc_ns, macro_hist[m], width=bin_width*1e9, align="center", alpha=0.35, color=c)
    ax[1].plot(tc_ns, macro_lam[m], color=c, lw=1.6,
               label=f"m={m} (Σf={macro_fsum[m]:.3f})")
ax[1].axvline(t0_ns, color="k", ls=":", alpha=0.6, label=f"ToF {t0_ns:.1f} ns")
ax[1].set_xlabel("时间 t [ns]"); ax[1].set_ylabel(f"计数 / 1ns bin ({N_shots} shots)")
ax[1].set_title("代表性宏像元直方图(柱=蒙卡, 线=解析期望)")
ax[1].legend(fontsize=8); ax[1].grid(alpha=0.3)
plt.tight_layout(); plt.show()

print(f"峰值宏像元 m={m_peak} 位于像斑中心(y≈{j0}); 越向两端(m→0 或 m→{n_macro-1})信号越弱, 最终只剩环境光底噪。")
print(f"底噪水平: 每宏像元 27 个 SPAD × 环境光, 边缘宏像元 m=0 总计数≈{macro_hist[0].sum():.0f} ({N_shots} shots)。")''')

# --- 组装 & 写出 ---
for typ, src in NEW:
    nb.cells.append(nbformat.v4.new_markdown_cell(src) if typ == "md"
                    else nbformat.v4.new_code_cell(src))
nbformat.write(nb, DST)
print(f"已生成 {DST}: 共 {len(nb.cells)} 个 cell (v10 全部 + 新增 {len(NEW)})。")
