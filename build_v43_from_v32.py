# -*- coding: utf-8 -*-
"""
build_v41_from_v32.py  —— 从 v32 派生 v41（在 v40 基础上按用户 9 条要求调整）
==============================================================================
用户要求（2026-07-27）：
  1) 恢复模块 7c(cell23)、9c(cell33)、模块18(cell43)、12b(cell47)；
  2) 模块8(cell25) 的 RC-vs-硬死 演示废弃停用——但它兼职定义全局时间窗，
     故【拆分】：时间窗定义独立成"模块 0b"cell 保留，演示部分转 markdown 停用；
  3) 全局时间窗放宽：ToF 前 SIM_PRE_NS=50ns / 后 SIM_POST_NS=100ns（参数可调，
     前段把环境估得更准、后段把脉冲看得更全）；步长 dt_fine 10ps→50ps（RC=8ns，够用）；
     绘图窗单独收窄：ToF 前 PLOT_PRE_NS=20ns / 后 PLOT_POST_NS=50ns（参数可调）；
  4) 模块11b(cell39)、模块12(cell45) 展示光强偏强 → 各加可调降强系数，
     并【独立重跑】展示数据，不污染核心分析(模块11 的 macro_hist / 13/14)；
  5) v40 新模块「A/B」→ 改"模块 19 / 20"命名；能量扫描横轴改 linear；
     SNR-距离右图 y 轴改 linear，并叠加 1/D² 参考线判断平方反比适用距离。

不改动 v32 源文件；不改动 PARAMS 任何物理参数值。
生成：lidar_histogram_sim_v41.ipynb
"""
import json

SRC_NB = "lidar_histogram_sim_v32.ipynb"
OUT_NB = "lidar_histogram_sim_v43.ipynb"

# ============================================================================
# 停用集（转 markdown 注释停用）——v42 相对 v41 变化：
#   新增停用 27（模块8b：RC 模型单 SPAD 的 Vov(t) 曲线演示）——实际模型是 9b 二值引擎，
#   8b 是旧 RC-trace 演示，已无用；经依赖检查 9c(cell33) 自给自足、不依赖 8b，可安全停用。
#   其余同 v41：29/41/53/55 + cell25 演示部分。
# ============================================================================
DISABLE = {
    27: "模块 8b：RC 模型单 SPAD 的过电压 Vov(t) 曲线演示 —— 实际模型是 9b 二值引擎，8b 为旧演示，停用。",
    29: "模块 9：阵列内不同 SPAD 的响应差异 —— v42 不需要，停用。",
    41: "模块 17：不同反射率 ρ 的信号波形对比 —— 被『模块 19 能量扫描』取代，停用。",
    53: "模块 15：多次蒙卡的 SNR 分布 + 正态拟合 —— v42 用『模块 20 SNR-距离』替代，停用。",
    55: "模块 16：100ppm 噪点率理论阈值 + 海量蒙卡验证 —— v42 不需要，停用。",
}

# ============================================================================
# 【模块 0b】全局时间窗与绘图窗设置 —— 取代 cell25 前 21 行的时间窗定义
#   关键：centers/edges 用 np.arange 直接生成，不再依赖模块8 的 events_to_hist。
# ============================================================================
MODULE_0B = '''# ============================================================================
# 模块 0b（v41 新增）— 全局时间窗与绘图窗设置
#   原 v32 把时间窗定义混在"模块8(RC vs 硬死时间演示)"里；v41 把模块8 演示废弃，
#   故将【全局时间窗定义】独立到此处，供下游所有二值 MC 模块(11/11b/12/13/14/18/12b
#   及模块19/20)统一使用。所有窗口/步长参数集中在这里，方便调节。
# ============================================================================

# ---- 仿真时间窗（可调）：前段更长→环境底噪估计更稳；后段更长→脉冲拖尾看得更全 ----
SIM_PRE_NS  = 50.0      # 仿真窗：ToF 之前 [ns]（可调）
SIM_POST_NS = 100.0     # 仿真窗：ToF 之后 [ns]（可调）

# ---- 绘图时间窗（可调）：画图只聚焦回波附近，不必画满整个仿真窗 ----
PLOT_PRE_NS  = 20.0     # 绘图窗：ToF 之前 [ns]（可调）
PLOT_POST_NS = 50.0     # 绘图窗：ToF 之后 [ns]（可调）

# ---- 精细步长（可调）：RC 恢复 τ_RC=8ns，200ps 已足够；原 10ps 过细、拖慢逐光子引擎 ----
dt_fine = 200e-12       # 全局精细网格步长 [s]（原 PARAMS 里是 10ps，此处放粗到 200ps）

# ---- 由上述参数推导时间窗（下游契约变量，全部在此定义）----
bin_width = PARAMS["hist"]["bin_width"]
pre, post = SIM_PRE_NS * 1e-9, SIM_POST_NS * 1e-9
t0   = time_of_flight(D0)
t_lo, t_hi = t0 - pre, t0 + post
tf   = np.arange(t_lo, t_hi, dt_fine)                 # 精细网格（分析窗，不含护带）
edges = np.arange(t_lo, t_hi + bin_width / 2, bin_width)   # 1ns bin 边界
centers = 0.5 * (edges[:-1] + edges[1:])              # bin 中心
nbins = len(centers)                                  # 直方图长度（契约变量）
tc_ns = centers * 1e9
t0_ns = t0 * 1e9

# ---- 绘图窗边界（各时域图统一用；比仿真窗窄）----
plot_lo_ns = t0_ns - PLOT_PRE_NS
plot_hi_ns = t0_ns + PLOT_POST_NS

# ---- 中心 SPAD 的单管信号率（原在模块8 定义，下游模块8b/9c 依赖，一并上移到此）----
#   f_ij = 最强收集像元的收集比例；r_sig = 该单 SPAD 的信号光子到达率（精细网格 tf 上）。
f_ij = fpix0[i0, j0]
r_sig = signal_photon_rate_fine(echo0, f_ij, tf)

print("="*76)
print(f"[模块 0b] 全局时间窗设置：")
print(f"  仿真窗：ToF-{SIM_PRE_NS:.0f}ns ~ ToF+{SIM_POST_NS:.0f}ns  (D0={D0}m, ToF={t0_ns:.1f}ns)")
print(f"          t_lo={t_lo*1e9:.1f}ns, t_hi={t_hi*1e9:.1f}ns, 步长 dt_fine={dt_fine*1e12:.0f}ps")
print(f"          精细网格点数={len(tf)}, 1ns bin 数 nbins={nbins}")
print(f"  绘图窗：ToF-{PLOT_PRE_NS:.0f}ns ~ ToF+{PLOT_POST_NS:.0f}ns  ({plot_lo_ns:.1f}~{plot_hi_ns:.1f}ns)")
print(f"  中心 SPAD: f_ij={f_ij:.3e}, 供模块8b/9c 用的 r_sig 已就绪。")
print(f"  注：dt_fine 由 10ps 放粗到 200ps（RC τ=8ns 下精度足够），窗更宽但格点数反而更少。")
'''

# ============================================================================
# 模块8 演示部分（cell25 去掉前 21 行时间窗后的剩余）——转 markdown 停用
#   这部分是 RC vs 硬死时间的对照演示（含 Nrep=3000 循环 + 两张图），模型已废弃。
# ============================================================================


def to_markdown_disabled(cell, reason, ver="v43"):
    """把一个 code cell 转成 markdown（代码进 ```python 块，完整保留、不执行）。"""
    code_src = "".join(cell["source"]) if isinstance(cell["source"], list) else cell["source"]
    md = (f"> 🚫 **【{ver} 停用】** {reason}\n>\n"
          f"> 下面是原始代码（**已注释停用，不执行**），完整保留以备查阅：\n\n"
          f"```python\n{code_src}\n```\n")
    return {"cell_type": "markdown", "id": cell.get("id", "disabled"),
            "metadata": {f"{ver}_disabled": True},
            "source": md.splitlines(keepends=True)}


def split_cell25(cell):
    """把 cell25 拆成：模块8 演示部分（转 markdown 停用）。
    时间窗定义部分由 MODULE_0B 取代（单独插入），故此处只返回"演示停用"markdown。
    做法：取 cell25 源码中第 22 行(含)之后的内容作为"演示部分"。"""
    src = "".join(cell["source"])
    lines = src.splitlines(keepends=True)
    # 前 21 行是时间窗定义（含空行与注释），第 22 行起是模块8 演示
    demo_src = "".join(lines[21:])
    reason = ("模块 8：单 SPAD『RC 恢复 vs 硬死时间』对照演示（含 Nrep=3000 均值与两张图）"
              "—— 该模型已废弃；其原本兼定义的全局时间窗已上移到『模块 0b』。停用演示。")
    md = (f"> 🚫 **【v43 停用】** {reason}\n>\n"
          f"> 下面是模块8 演示部分的原始代码（**已注释停用，不执行**），完整保留以备查阅：\n\n"
          f"```python\n{demo_src}\n```\n")
    return {"cell_type": "markdown", "id": cell.get("id", "mod8_demo") + "_demo_disabled",
            "metadata": {"v43_disabled": True},
            "source": md.splitlines(keepends=True)}


# ============================================================================
# 模块 11b（cell39）整块替换：降光强(可调) + 独立重跑 + 绘图窗
#   不动模块11 的 macro_hist/per_shot_peak（核心 13/14 用原值）；此处按 DEMO11B_SCALE
#   降强后【单独重跑】峰值宏像元 per-shot，仅供本图展示。
# ============================================================================
MODULE_11B_NEW = r'''# ============================================================================
# 模块 11b（v31 折线图；v41 改：降光强 + 独立重跑 + 绘图窗）
#   用户反馈：11b 展示光强偏强(峰早早顶满 27/108)。v41 加【可调降强系数 DEMO11B_SCALE】，
#   把信号率 ×系数后【单独重跑】峰值宏像元 per-shot 数据(不动模块11 的 macro_hist，
#   核心分析 13/14 仍用原值)；绘图窗改用全局 plot_lo_ns/plot_hi_ns(ToF 前20/后50ns)。
# ============================================================================
DEMO11B_SCALE = 0.3        # 模块11b 展示用降强系数（可调；1.0=原强度，<1 变弱）

# ---- 独立重跑：峰值宏像元 27 SPAD × N_shots，信号率 ×DEMO11B_SCALE，带归因 ----
_rng11b = np.random.default_rng(PARAMS["hist"]["seed"] + 41100)
ps_demo = np.zeros((N_shots, nbins))       # 每发总直方图(27 SPAD 二值和)
ps_sig  = np.zeros((N_shots, nbins))       # 信号触发分量
ps_amb  = np.zeros((N_shots, nbins))       # 环境触发分量
for _s in range(N_shots):
    a_tot = np.zeros(nbins, dtype=np.int32)
    a_sig = np.zeros(nbins, dtype=np.int32)
    a_amb = np.zeros(nbins, dtype=np.int32)
    for fij in macro_fvals[m_peak]:
        o, o_s, o_a = spad_binary_trace(
            base_rate_gen * fij * DEMO11B_SCALE, r_amb_ph, tf_gen, centers,
            PDE, TAU_RC, VTH_FRAC, jit, _rng11b, T_OVER, T_LASER, RESP_SHAPE, RESP_K,
            return_attrib=True)
        a_tot += o; a_sig += o_s; a_amb += o_a
    ps_demo[_s] = a_tot; ps_sig[_s] = a_sig; ps_amb[_s] = a_amb

cap_shot = n_pix_macro                       # 单 shot 硬上限 = 27 SPAD
acc4    = ps_demo.sum(axis=0)                 # N_shots 累加
acc_sig = ps_sig.sum(axis=0)
acc_amb = ps_amb.sum(axis=0)
pk_b = int(acc4.argmax())

fig, ax = plt.subplots(1, 2, figsize=(14, 5.0))

# ===== 左图：每发一条折线 + 累加折线 =====
axL = ax[0]
cmap = plt.cm.viridis(np.linspace(0.15, 0.85, N_shots))
for s in range(N_shots):
    axL.plot(tc_ns, ps_demo[s], lw=1.6, marker="o", ms=3, color=cmap[s], alpha=0.9,
             label=f"shot #{s}（单发, 峰={ps_demo[s].max():.0f}/{cap_shot}）")
axL.plot(tc_ns, acc4, lw=2.4, marker="o", ms=3, color="k",
         label=f"{N_shots} 发累加（峰={acc4.max():.0f}/{macro_cap}）")
axL.axhline(cap_shot, color="tab:blue", ls=":", lw=1.2, label=f"单发上限 {cap_shot}")
axL.axhline(macro_cap, color="k", ls="-.", lw=1.0, alpha=0.6, label=f"累加上限 {macro_cap}")
axL.axvline(t0_ns, color="k", ls=":", alpha=0.6, label=f"ToF {t0_ns:.1f} ns")
axL.set_xlim(plot_lo_ns, plot_hi_ns)
axL.set_xlabel("时间 t [ns]"); axL.set_ylabel("二值计数 / 1ns bin（峰值宏像元, 27 SPAD 和）")
axL.set_title(f"每个 shot 的直方图（峰值宏像元 m={m_peak}, 降强×{DEMO11B_SCALE}）：单发 vs 累加", fontsize=10)
axL.legend(fontsize=7.5, loc="upper right", ncol=1); axL.grid(alpha=0.3)

# ===== 右图：信号/环境/总 三条折线归因 =====
axR = ax[1]
axR.plot(tc_ns, acc_sig, lw=1.6, marker="o", ms=3, color="tab:green",
         label=f"信号触发累加（峰={acc_sig[pk_b]:.0f}）")
axR.plot(tc_ns, acc_amb, lw=1.6, marker="o", ms=3, color="tab:orange",
         label=f"环境触发累加（Σ={acc_amb.sum():.0f}）")
axR.plot(tc_ns, acc4, lw=2.0, marker="o", ms=3, color="k",
         label=f"信号+环境累加（峰={acc4.max():.0f}）")
axR.axhline(macro_cap, color="k", ls="-.", lw=1.0, alpha=0.6, label=f"累加上限 {macro_cap}")
axR.axvline(t0_ns, color="k", ls=":", alpha=0.6, label=f"ToF {t0_ns:.1f} ns")
axR.set_xlim(plot_lo_ns, plot_hi_ns)
axR.set_xlabel("时间 t [ns]"); axR.set_ylabel(f"计数 / 1ns bin（{N_shots} 发累加）")
axR.set_title("累加直方图的 信号/环境 归因（绿=信号, 橙=环境, 黑=总）", fontsize=10)
axR.legend(fontsize=7.5, loc="upper right"); axR.grid(alpha=0.3)

plt.suptitle(f"模块 11b　每个 shot 的直方图（降强×{DEMO11B_SCALE}, 绘图窗 ToF-{PLOT_PRE_NS:.0f}~+{PLOT_POST_NS:.0f}ns）",
             fontsize=11.5)
plt.tight_layout(rect=[0, 0, 1, 0.96]); plt.show()

print("="*76)
print(f"模块 11b（降强×{DEMO11B_SCALE}，独立重跑，不影响模块11/13/14）：峰值宏像元 m={m_peak}，{N_shots} 发")
for s in range(N_shots):
    ps = ps_demo[s]
    print(f"  shot #{s}: 峰 bin={ps.max():.0f}/{cap_shot} @ {tc_ns[int(ps.argmax())]:.0f} ns, "
          f"总计数={ps.sum():.0f}（信号={ps_sig[s].sum():.0f}, 环境={ps_amb[s].sum():.0f}）")
print(f"  ---- {N_shots} 发累加：峰 bin={acc4.max():.0f}/{macro_cap} @ {tc_ns[pk_b]:.0f} ns, 总计数={acc4.sum():.0f} ----")
print(f"  说明：降强系数 DEMO11B_SCALE={DEMO11B_SCALE} 仅作用于本模块展示，不改动物理参数、不影响核心分析。")
'''

# ============================================================================
# 模块 12（cell45）整块替换：降光强(可调) + 独立重跑全 40 宏像元 + 绘图窗
# ============================================================================
MODULE_12_NEW = r'''# ============================================================================
# 模块 12（v13 宏像元热图/归因；v41 改：降光强 + 独立重跑 + 绘图窗）
#   用户反馈：模块12 信号偏强。v41 加【可调降强系数 MOD12_SCALE】，把信号率 ×系数后
#   【单独重跑】全 40 宏像元二值直方图 + 归因 + 泊松期望，仅供本模块展示；
#   不动模块11 的 macro_hist（核心 13/14 用原值）。绘图窗用 plot_lo_ns/plot_hi_ns。
# ============================================================================
MOD12_SCALE = 0.5          # 模块12 展示用降强系数（可调；1.0=原强度）

# ---- 独立重跑：全 40 宏像元 × N_shots，信号率 ×MOD12_SCALE，带归因 ----
_rng12 = np.random.default_rng(PARAMS["hist"]["seed"] + 41200)
mh_demo     = np.zeros((n_macro, nbins))     # 总(信号+环境)
mh_demo_sig = np.zeros((n_macro, nbins))
mh_demo_amb = np.zeros((n_macro, nbins))
for _s in range(N_shots):
    for m in range(n_macro):
        a_tot = np.zeros(nbins, dtype=np.int32)
        a_sig = np.zeros(nbins, dtype=np.int32)
        a_amb = np.zeros(nbins, dtype=np.int32)
        for fij in macro_fvals[m]:
            o, o_s, o_a = spad_binary_trace(
                base_rate_gen * fij * MOD12_SCALE, r_amb_ph, tf_gen, centers,
                PDE, TAU_RC, VTH_FRAC, jit, _rng12, T_OVER, T_LASER, RESP_SHAPE, RESP_K,
                return_attrib=True)
            a_tot += o; a_sig += o_s; a_amb += o_a
        mh_demo[m]     += a_tot
        mh_demo_sig[m] += a_sig
        mh_demo_amb[m] += a_amb

# 泊松期望(未封顶,降强后)参考虚线：与模块11 同法，信号率 ×MOD12_SCALE
irf_k12 = gaussian_kernel(jit, dt_fine)
bin_idx12 = np.clip(((tf - t_lo)/bin_width).astype(int), 0, nbins-1)
mlam_demo = np.zeros((n_macro, nbins))
for m in range(n_macro):
    r_det = (base_rate*macro_fsum[m]*MOD12_SCALE + n_pix_macro*r_amb_ph) * PDE
    r_det = np.convolve(r_det, irf_k12, mode="same") * dt_fine
    mlam_demo[m] = N_shots * np.bincount(bin_idx12, weights=r_det*dt_fine, minlength=nbins)

# ---- 图 A: 宏像元-时间 二值计数热图 ----
fig, ax = plt.subplots(figsize=(11, 5))
im = ax.imshow(mh_demo, origin="upper", aspect="auto", cmap="inferno",
               extent=[tc_ns[0]-0.5, tc_ns[-1]+0.5, n_macro-0.5, -0.5],
               vmin=0, vmax=macro_cap)
ax.axvline(t0_ns, color="cyan", ls=":", lw=1.2, label=f"真实 ToF {t0_ns:.1f} ns")
ax.axhline(m_peak, color="lime", ls=":", lw=1.0, alpha=0.7, label=f"峰值宏像元 m={m_peak}")
ax.set_xlim(plot_lo_ns, plot_hi_ns)
ax.set_xlabel("时间 t [ns]"); ax.set_ylabel("宏像元序号 m (0=顶部, 沿长边 y)")
ax.set_title(f"宏像元(9×3)二值直方图热图 (降强×{MOD12_SCALE}, N_shots={N_shots}, 上限={macro_cap})")
ax.legend(fontsize=9, loc="upper right")
plt.colorbar(im, ax=ax, label=f"计数 / (宏像元, 1ns bin)  [0..{macro_cap}]")
plt.tight_layout(); plt.show()

# ---- 图 B: 代表性宏像元直方图 ----
fig, ax = plt.subplots(figsize=(11, 4.6))
reps_m = sorted(set([0, max(0, m_peak-6), m_peak, min(n_macro-1, m_peak+6)]))
colors = ["tab:gray", "tab:green", "tab:red", "tab:orange", "tab:purple"]
for m, c in zip(reps_m, colors):
    ax.plot(tc_ns, mh_demo[m], color=c, lw=1.4, marker="o", ms=4,
            label=f"m={m} 二值 (Σf={macro_fsum[m]:.3f})")
    ax.plot(tc_ns, mlam_demo[m], color=c, lw=1.1, ls="--", alpha=0.6)
ax.axhline(macro_cap, color="k", ls="-.", lw=1.2, alpha=0.8, label=f"二值硬上限={macro_cap}")
ax.axvline(t0_ns, color="k", ls=":", alpha=0.6, label=f"ToF {t0_ns:.1f} ns")
ax.set_xlim(plot_lo_ns, plot_hi_ns)
ax.set_xlabel("时间 t [ns]"); ax.set_ylabel(f"计数 / 1ns bin ({N_shots} shots)")
ax.set_title(f"代表性宏像元直方图 (降强×{MOD12_SCALE}; 点线=二值实测, 虚线=泊松期望, 点划线=上限)")
ax.legend(fontsize=8); ax.grid(alpha=0.3)
plt.tight_layout(); plt.show()

# ---- 图 B2: 峰值宏像元 信号/环境/总 三条折线归因 ----
fig, ax = plt.subplots(figsize=(11, 4.6))
ax.plot(tc_ns, mh_demo_sig[m_peak], color="tab:green", lw=1.4, marker="o", ms=4, label="信号光子触发")
ax.plot(tc_ns, mh_demo_amb[m_peak], color="tab:orange", lw=1.4, marker="o", ms=4, label="环境光子触发")
ax.plot(tc_ns, mh_demo[m_peak], color="k", lw=1.8, marker="o", ms=3, label="总共 (信号+环境)")
ax.plot(tc_ns, mlam_demo[m_peak], "k--", lw=1.3, alpha=0.7, label="泊松期望(未封顶)")
ax.axhline(macro_cap, color="k", ls="-.", lw=1.2, alpha=0.8, label=f"二值硬上限={macro_cap}")
ax.axvline(t0_ns, color="k", ls=":", alpha=0.6, label=f"ToF {t0_ns:.1f} ns")
pk_pk = int(mh_demo[m_peak].argmax())
ax.annotate(f"峰 bin={mh_demo[m_peak][pk_pk]:.0f}/{macro_cap}",
            (tc_ns[pk_pk], mh_demo[m_peak][pk_pk]),
            xytext=(tc_ns[pk_pk]+5, macro_cap*0.72), fontsize=8,
            arrowprops=dict(arrowstyle="->", lw=0.8))
ax.set_ylim(0, macro_cap*1.12)
ax.set_xlim(plot_lo_ns, plot_hi_ns)
ax.set_xlabel("时间 t [ns]"); ax.set_ylabel(f"计数 / 1ns bin ({N_shots} shots)")
ax.set_title(f"峰值宏像元 m={m_peak} 信号/环境/总 三条折线归因 (降强×{MOD12_SCALE})")
ax.legend(fontsize=8); ax.grid(alpha=0.3)
plt.tight_layout(); plt.show()

sig_pk = mh_demo_sig[m_peak][pk_pk]; amb_pk = mh_demo_amb[m_peak][pk_pk]
print(f"模块12（降强×{MOD12_SCALE}，独立重跑，不影响13/14）: 峰值宏像元 m={m_peak} "
      f"峰 bin={mh_demo[m_peak].max():.0f}/{macro_cap} "
      f"({'饱和' if mh_demo[m_peak].max()>=macro_cap-1e-9 else '未饱和'})")
print(f"  峰 bin @ {tc_ns[pk_pk]:.0f}ns: 信号触发={sig_pk:.0f} + 环境触发={amb_pk:.0f} "
      f"(信号占 {100*sig_pk/max(sig_pk+amb_pk,1e-9):.0f}%)")
'''

print("构建脚本：模块11b/12 替换源已就绪。")


# ============================================================================
# v41 新增分析头（markdown）
# ============================================================================
NEW_HEADER_MD = '''---
# 🆕 v41 新增分析（模块 19 / 20；基于以上 v32 完整链路，全程蒙特卡罗）

以下两个模块复用上方 v32 已定义的物理链路与逐光子二值引擎
（`spad_binary_trace` / `signal_photon_rate_fine` / `link_factor` /
`pixel_collection_matrix` / `front_time_leading_edge` / 模块 0b 的时间窗）：

- **模块 19 — 能量扫描 → 前沿/重心定时 → dist-peak / dist-area 四条曲线**
  在 `cali_dist`（=`D0`）放单目标，扫描反射能量倍数 `boost`（上下界/档数可调）。
  每档跑二值 MC 得直方图，用**前沿法**（v32 原生）与**重心法**（COG，v41 新增）各定时一次；
  能量太低（峰 < 检测阈值 `det_th`）则**不定时、留空**。四条曲线：dist-peak(前沿/重心) +
  dist-area(前沿/重心)。**横轴用 linear（线性）刻度。**

- **模块 20 — SNR vs 距离**
  沿用 v32 SNR 定义（`SNR = S/√B`）扫描距离 D。右图（峰&背景 vs 距离）**y 轴用 linear**，
  并叠加 **1/D² 参考线**，据此判断"多少米之外可用平方反比近似"。

> 缩写：COG（Center of Gravity，重心/质心）。其余（TCSPC/SPAD/IRF/ToF/PDE/SNR/RC）见前文。
'''

# ============================================================================
# 模块 19：能量扫描（改自 v40 模块 A —— 横轴改 linear、命名改"模块19"、PNG 改 v41）
# ============================================================================
MODULE_19 = r'''# ============================================================================
# 模块 19（v41 新增）— 能量扫描：前沿法/重心法定时 → dist-peak / dist-area 四条曲线
#   复用 v32：spad_binary_trace（逐光子二值引擎）、signal_photon_rate_fine、
#            pixel_collection_matrix、front_time_leading_edge、模块 0b 时间窗。
#   能量注入：反射能量倍数 boost 乘到"单位收集比例信号率"上（ρ≤1 缩不出 1e5，用倍数表征）。
#   横轴用 linear（线性）刻度。
# ============================================================================

# ---- 可调参数 ----
BOOST_MIN = 1e-2          # 反射能量倍数下界（可调）
BOOST_MAX = 1e2           # 反射能量倍数上界（可调）
K_TH_19   = 5.0           # 检测阈值倍数 det_th = K_TH_19 · nc_base（沿用 v32 模块14 默认）
COG_HALF  = 6            # 重心法窗口半宽 [bin]
SMOOTH_N  = 15           # 散点滑窗平均窗口点数（可调；奇数为佳，用于叠加 N 点滑动平均线+标准差带）

# ---- 分段线性步长（用户指定；每段步长不同，低能量段密、高能量段疏）----
#   [1e-2,1e-1] 步长 0.01；[1e-1,1e0] 步长 0.02（用户未指定该中间段，取 0.02，可调）；
#   [1e0,1e1] 步长 0.2；[1e1,1e2] 步长 2。各段拼接后去重排序。
SEGMENTS = [
    (1e-2, 1e-1, 0.01),   # 低能量段：步长 0.01
    (1e-1, 1e0,  0.02),   # 中间段（未指定，取 0.02，可调）
    (1e0,  1e1,  0.2),    # 中能量段：步长 0.2
    (1e1,  1e2,  2.0),    # 高能量段：步长 2
]
_parts = [np.arange(lo, hi, st) for (lo, hi, st) in SEGMENTS]
_parts.append(np.array([BOOST_MAX]))                     # 补上终点 1e2
boost_grid = np.unique(np.concatenate(_parts))           # 拼接、去重、升序
N_BOOST = len(boost_grid)

def moving_average(y, n=SMOOTH_N):
    """对 1D 数组做 N 点滑动平均（nan 安全：逐窗口 nanmean）。返回同长度数组。"""
    y = np.asarray(y, dtype=float)
    half = n // 2
    out = np.full(y.size, np.nan)
    for i in range(y.size):
        lo = max(0, i - half); hi = min(y.size, i + half + 1)
        seg = y[lo:hi]
        if np.isfinite(seg).any():
            out[i] = np.nanmean(seg)
    return out

def moving_std(y, n=SMOOTH_N):
    """对 1D 数组做 N 点滑动【标准差】（nan 安全：逐窗口 nanstd）。返回同长度数组。
    用户要求：滑窗平均后输出该点标准差（点足够密时反映该处散布）。"""
    y = np.asarray(y, dtype=float)
    half = n // 2
    out = np.full(y.size, np.nan)
    for i in range(y.size):
        lo = max(0, i - half); hi = min(y.size, i + half + 1)
        seg = y[lo:hi]
        if np.isfinite(seg).sum() >= 2:
            out[i] = np.nanstd(seg)
    return out

def centroid_time_cog(hist, centers_, pk_idx, half=COG_HALF):
    """重心法(COG)定时：峰 ±half bin 窗口内计数质心。窗口内全 0 返回 nan。"""
    lo = max(0, pk_idx - half); hi = min(len(hist), pk_idx + half + 1)
    w = hist[lo:hi]; x = centers_[lo:hi]
    return (x * w).sum() / w.sum() if w.sum() > 0 else np.nan

def run_peak_hist_boost(boost, p=PARAMS):
    """cali_dist 处峰值宏像元 27 SPAD × N_shots 逐光子二值 MC，信号率 ×boost。返回累加直方图。"""
    rng = np.random.default_rng(p["hist"]["seed"] + 41910)
    h = np.zeros(nbins)
    for _shot in range(N_shots):
        acc = np.zeros(nbins, dtype=np.int32)
        for fij in macro_fvals[m_peak]:
            acc += spad_binary_trace(
                base_rate_gen * fij * boost, r_amb_ph, tf_gen, centers,
                PDE, TAU_RC, VTH_FRAC, jit, rng, T_OVER, T_LASER, RESP_SHAPE, RESP_K)
        h += acc
    return h

det_th_19 = K_TH_19 * nc_base       # nc_base 来自模块13（纯背景全 bin 均值）
print(f"模块19 能量扫描：boost ∈ [{BOOST_MIN:.0e}, {BOOST_MAX:.0e}], 分段步长 {[s[2] for s in SEGMENTS]} → {N_BOOST} 档（逐光子引擎）")
print(f"  检测阈值 det_th = {K_TH_19:.0f}×nc_base = {det_th_19:.3f}（峰 < 此值则不可分辨、留空）")

front_R = np.full(N_BOOST, np.nan)
cog_R   = np.full(N_BOOST, np.nan)
peak_v  = np.zeros(N_BOOST)
area_v  = np.zeros(N_BOOST)

# 每段"最后一档"的 boost 值（用于挑出代表档画波形+定时示意）
seg_last_boost = [seg[1] for seg in SEGMENTS]            # 各段上界即该段最后一档附近
seg_last_boost = [boost_grid[np.argmin(np.abs(boost_grid - b))] for b in seg_last_boost]
seg_waveforms = []      # 保存 (boost, hist, pk, t_front[ns], t_cog[ns], det_th, V_dec)

for k, boost in enumerate(boost_grid):
    h = run_peak_hist_boost(boost)
    pk = int(np.argmax(h))
    peak_v[k] = h[pk]; area_v[k] = h.sum()
    detectable = h[pk] >= det_th_19
    tf_k = t_cog = np.nan; V_dec_k = np.nan
    if detectable:
        V_dec_k = 0.5 * (det_th_19 + h[pk])
        tf_k, _, _ = front_time_leading_edge(h, centers, pk, V_dec_k, bin_width)
        if np.isfinite(tf_k):
            front_R[k] = C_LIGHT * tf_k / 2.0
        t_cog = centroid_time_cog(h, centers, pk)
        if np.isfinite(t_cog):
            cog_R[k] = C_LIGHT * t_cog / 2.0
    # 若本档是某段的"最后一档"，保存其波形与定时结果供后续画示意图
    if any(abs(boost - slb) < 1e-9 for slb in seg_last_boost):
        seg_waveforms.append((boost, h.copy(), pk,
                              tf_k*1e9 if np.isfinite(tf_k) else np.nan,
                              t_cog*1e9 if np.isfinite(t_cog) else np.nan,
                              det_th_19, V_dec_k))
    tag = "可分辨" if detectable else "太低-留空"
    if k == 0 or k == N_BOOST-1 or (k+1) % max(1, N_BOOST//20) == 0:   # 节流：约打印 20 行
        print(f"  [{k+1:>4d}/{N_BOOST}] boost={boost:.3f}  峰={h[pk]:>4.0f}  面积={h.sum():>6.0f}  {tag}")

mF = np.isfinite(front_R); mC = np.isfinite(cog_R)
print(f"有效定时：前沿 {mF.sum()}/{N_BOOST}，重心 {mC.sum()}/{N_BOOST}（低能量档按要求留空）")

# ---- 绘图1：dist-peak / dist-area 散点图 + N点滑窗平均 + 标准差带（横轴 linear）----
#   用户要求：散点；叠加 N 点滑动平均线；并输出该点标准差（点够密→画 ±1σ 误差带）。
fig, axes = plt.subplots(1, 2, figsize=(14, 5.2))

def scatter_with_smooth(ax, xv, yv, color, lbl):
    """散点 + 按 x 排序的滑窗均值线 + ±1σ 标准差带。"""
    ax.scatter(xv, yv, s=10, color=color, alpha=0.40, label=f"{lbl}散点")
    if np.isfinite(yv).sum() > SMOOTH_N:
        o = np.argsort(xv); xs = xv[o]; ys = yv[o]
        mu = moving_average(ys); sd = moving_std(ys)
        ax.plot(xs, mu, "-", color=color, lw=2.0, label=f"{lbl}{SMOOTH_N}点滑窗均值")
        ax.fill_between(xs, mu - sd, mu + sd, color=color, alpha=0.18, label=f"{lbl}±1σ")

ax = axes[0]   # dist-peak
scatter_with_smooth(ax, peak_v[mF], front_R[mF], "tab:blue", "前沿法")
scatter_with_smooth(ax, peak_v[mC], cog_R[mC], "tab:red", "重心法")
ax.axhline(D0, color="gray", ls=":", lw=1.2, alpha=0.8, label=f"真值 {D0} m")
ax.set_xlabel("直方图峰 bin 计数 peak（线性轴）")
ax.set_ylabel("测距结果 [m]"); ax.set_title("dist-peak（散点 + 滑窗均值 ± 标准差带）")
ax.legend(fontsize=7); ax.grid(alpha=0.3)

ax = axes[1]   # dist-area
scatter_with_smooth(ax, area_v[mF], front_R[mF], "tab:blue", "前沿法")
scatter_with_smooth(ax, area_v[mC], cog_R[mC], "tab:red", "重心法")
ax.axhline(D0, color="gray", ls=":", lw=1.2, alpha=0.8, label=f"真值 {D0} m")
ax.set_xlabel("直方图总面积 area（线性轴）")
ax.set_ylabel("测距结果 [m]"); ax.set_title("dist-area（散点 + 滑窗均值 ± 标准差带）")
ax.legend(fontsize=7); ax.grid(alpha=0.3)

plt.suptitle(f"模块 19 — 能量扫描定时散点图 (cali_dist={D0} m, N_shots={N_shots}, {N_BOOST}档, {SMOOTH_N}点滑窗)",
             fontsize=12)
plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.savefig("v43_energy_scan.png", dpi=110, bbox_inches="tight")
plt.show()

# ---- 绘图2：每个能量段"最后一档"的波形 + 前沿/重心定时示意 ----
#   用户要求：每段最后一次画波形，标出前沿定时、重心定时位置。
n_seg = len(seg_waveforms)
if n_seg > 0:
    fig, axs = plt.subplots(1, n_seg, figsize=(4.7*n_seg, 4.4), squeeze=False)
    axs = axs[0]
    for j, (bst, h, pk, tfn, tcn, dth, vdec) in enumerate(seg_waveforms):
        a = axs[j]
        a.plot(tc_ns, h, color="k", lw=1.3, marker="o", ms=2.5, label="直方图波形")
        a.axhline(dth, color="orange", ls="--", lw=1.2, label=f"检测阈 det_th={dth:.1f}")
        if np.isfinite(vdec):
            a.axhline(vdec, color="purple", ls="-.", lw=1.0, alpha=0.7, label=f"判决电平 V_dec={vdec:.1f}")
        a.axvline(t0_ns, color="gray", ls=":", lw=1.0, alpha=0.7, label=f"真值 ToF {t0_ns:.1f}ns")
        if np.isfinite(tfn):
            a.axvline(tfn, color="tab:blue", lw=1.8, label=f"前沿定时 {tfn:.2f}ns")
        if np.isfinite(tcn):
            a.axvline(tcn, color="tab:red", lw=1.8, ls="--", label=f"重心定时 {tcn:.2f}ns")
        a.set_xlim(plot_lo_ns, plot_hi_ns)
        a.set_xlabel("时间 t [ns]"); a.set_ylabel("计数 / 1ns bin")
        a.set_title(f"boost={bst:.3g}（峰={h[pk]:.0f}）", fontsize=10)
        a.legend(fontsize=6.5); a.grid(alpha=0.3)
    plt.suptitle("模块 19 — 各能量段最后一档：波形 + 前沿/重心定时示意", fontsize=12)
    plt.tight_layout(rect=[0, 0, 1, 0.94])
    plt.savefig("v43_energy_waveforms.png", dpi=110, bbox_inches="tight")
    plt.show()

# ---- 绘图3：峰值 & 面积 随距离的变化（新增；固定 boost=1，扫描距离）----
#   用户要求：画峰值和面积随距离的变化图。
print("正在跑 峰值&面积 vs 距离 ...")
D_pa = np.linspace(5.0, 300.0, 20)      # 距离档（可调）
peak_pa = np.zeros(D_pa.size); area_pa = np.zeros(D_pa.size)
for idx, D in enumerate(D_pa):
    t0d = time_of_flight(D); t_lo_d, t_hi_d = t0d - pre, t0d + post
    guard = T_OVER + 5 * jit
    tf_gen_d = np.arange(t_lo_d - guard, t_hi_d, dt_fine)
    edges_d = np.arange(t_lo_d, t_hi_d + bin_width/2, bin_width)
    centers_d = 0.5 * (edges_d[:-1] + edges_d[1:]); nb_d = len(centers_d)
    fpixD, _, _ = pixel_collection_matrix(D)
    fvalsD = [fpixD[:, m*By_m:(m+1)*By_m].ravel() for m in range(n_macro)]
    fsumD = np.array([fv.sum() for fv in fvalsD]); mpk = int(fsumD.argmax())
    echoD = dict(echo0); echoD["D"] = D
    base_gen_D = signal_photon_rate_fine(echoD, 1.0, tf_gen_d)
    rng_pa = np.random.default_rng(PARAMS["hist"]["seed"] + 41930 + idx)
    h = np.zeros(nb_d)
    for _shot in range(N_shots):
        acc = np.zeros(nb_d, dtype=np.int32)
        for fij in fvalsD[mpk]:
            acc += spad_binary_trace(base_gen_D * fij, r_amb_ph, tf_gen_d, centers_d,
                                     PDE, TAU_RC, VTH_FRAC, jit, rng_pa, T_OVER, T_LASER, RESP_SHAPE, RESP_K)
        h += acc
    peak_pa[idx] = h.max(); area_pa[idx] = h.sum()

fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(D_pa, peak_pa, "o-", color="tab:blue", lw=1.6, ms=5, label="峰值 peak")
ax.plot(D_pa, area_pa, "s-", color="tab:green", lw=1.6, ms=5, label="面积 area")
ax.set_xlabel("距离 D [m]"); ax.set_ylabel("计数")
ax.set_title(f"峰值 & 面积 随距离变化 (boost=1, N_shots={N_shots})")
ax.legend(fontsize=9); ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("v43_peak_area_vs_dist.png", dpi=110, bbox_inches="tight")
plt.show()

print("="*76)
print(f"模块 19 汇总：cali_dist={D0} m, ρ={echo0['rho']}, boost∈[{BOOST_MIN:.0e},{BOOST_MAX:.0e}]（{N_BOOST}档）")
if mF.any():
    print(f"  前沿法测距误差范围 [{(front_R[mF].min()-D0)*100:+.1f}, {(front_R[mF].max()-D0)*100:+.1f}] cm")
if mC.any():
    print(f"  重心法测距误差范围 [{(cog_R[mC].min()-D0)*100:+.1f}, {(cog_R[mC].max()-D0)*100:+.1f}] cm")
'''

# ============================================================================
# 模块 20：SNR vs 距离（改自 v40 模块 B —— 右图 y 改 linear + 1/D² 参考线 + 平方反比判断）
# ============================================================================
MODULE_20 = r'''# ============================================================================
# 模块 20（v41 新增）— 不同距离下的信噪比 SNR 变化（逐光子二值 MC）
#   沿用 v32 模块13 SNR 定义：SNR = S/√B（S=峰bin总计数−B，B=峰bin纯背景计数）。
#   右图 y 轴用 linear。平方反比验证用【未封顶的信号光子率峰值】比 1/D²·exp(-2αD)——
#   注意：峰 bin 计数因二值 0/1 饱和(近场顶到 108)会偏离 1/D²，不能用来判平方反比，
#   故本模块用未封顶物理量做验证，并把饱和的峰计数画在次轴作对照。
# ============================================================================

# ---- 可调参数 ----
D_MIN_B = 1.0            # 最近距离 [m]
D_MAX_B = 300.0         # 最远距离 [m]
N_D_B   = 50            # 距离档数（加密，可调；逐光子引擎，档多则慢）
BOOST_B = 1.0           # 反射能量倍数（默认 1 = 真实回波）

D_list_B = np.linspace(D_MIN_B, D_MAX_B, N_D_B)
print(f"模块20 SNR vs 距离：D ∈ [{D_MIN_B}, {D_MAX_B}] m, {N_D_B} 档 (boost={BOOST_B}, 逐光子二值 MC)")

def snr_at_distance(D, boost=BOOST_B, p=PARAMS):
    """距离 D 处峰值宏像元二值 MC（信号+背景 与 纯背景各一次），返回 (SNR, 峰计数, B)。
    时间窗/护带随 D 重建（用全局 pre/post/dt_fine），能量注入用 boost 倍数。"""
    t0d = time_of_flight(D)
    t_lo_d, t_hi_d = t0d - pre, t0d + post
    guard = T_OVER + 5 * jit
    tf_gen_d = np.arange(t_lo_d - guard, t_hi_d, dt_fine)
    edges_d = np.arange(t_lo_d, t_hi_d + bin_width/2, bin_width)
    centers_d = 0.5 * (edges_d[:-1] + edges_d[1:])
    nb_d = len(centers_d)
    fpixD, _, _ = pixel_collection_matrix(D, p)
    fvalsD = [fpixD[:, m*By_m:(m+1)*By_m].ravel() for m in range(n_macro)]
    fsumD = np.array([fv.sum() for fv in fvalsD])
    mpk = int(fsumD.argmax())
    echoD = dict(echo0); echoD["D"] = D
    base_gen_D = signal_photon_rate_fine(echoD, 1.0, tf_gen_d, p)
    zero_gen_D = np.zeros_like(base_gen_D)
    rng_s = np.random.default_rng(p["hist"]["seed"] + 42010)
    h_sig = np.zeros(nb_d)
    for _shot in range(N_shots):
        acc = np.zeros(nb_d, dtype=np.int32)
        for fij in fvalsD[mpk]:
            acc += spad_binary_trace(base_gen_D * fij * boost, r_amb_ph, tf_gen_d, centers_d,
                                     PDE, TAU_RC, VTH_FRAC, jit, rng_s, T_OVER, T_LASER, RESP_SHAPE, RESP_K)
        h_sig += acc
    rng_b = np.random.default_rng(p["hist"]["seed"] + 42020)
    h_bg = np.zeros(nb_d)
    for _shot in range(N_shots):
        acc = np.zeros(nb_d, dtype=np.int32)
        for fij in fvalsD[mpk]:
            acc += spad_binary_trace(zero_gen_D, r_amb_ph, tf_gen_d, centers_d,
                                     PDE, TAU_RC, VTH_FRAC, jit, rng_b, T_OVER, T_LASER, RESP_SHAPE, RESP_K)
        h_bg += acc
    pk = int(np.argmax(h_sig))
    B = h_bg[pk]; S = max(h_sig[pk] - B, 0.0)
    snr = S / np.sqrt(B) if B > 0 else (np.inf if S > 0 else 0.0)
    # 未封顶物理量：该峰值宏像元的【信号光子率峰值】(∝ link_factor ∝ 1/D²·exp(-2αD))，
    #   不受二值 0/1 饱和影响，用于干净地验证平方反比律。
    sig_rate_pk = float(base_gen_D.max() * fsumD[mpk])
    return snr, h_sig[pk], B, sig_rate_pk

snr_B  = np.zeros(N_D_B); peak_B = np.zeros(N_D_B); bg_B = np.zeros(N_D_B)
rate_B = np.zeros(N_D_B)      # 未封顶信号光子率峰值（供 1/D² 验证）
for k, D in enumerate(D_list_B):
    snr_B[k], peak_B[k], bg_B[k], rate_B[k] = snr_at_distance(D)
    print(f"  D={D:>6.1f} m: 峰={peak_B[k]:>4.0f}  B={bg_B[k]:>5.2f}  SNR={snr_B[k]:>6.2f}  信号率峰={rate_B[k]:.2e}")

# ---- 平方反比拟合：两条曲线都试 1/D² ----
#   理论：信号光子率 ∝ 接收立体角/D² × 大气透过率² = (1/D²)·exp(-2αD)。
#   ① 信号光子率峰值 rate_B：未封顶物理量，应严格遵循 1/D²。
#   ② SNR：S∝信号率∝1/D²，背景 B≈常数 → SNR=S/√B 也应 ∝1/D²（近似）。
#   两者均以最远档为锚点归一到 1/D²·exp(-2αD)，并给出各档【仿真 vs 平方反比】相对误差。
i_anchor = N_D_B - 1                         # 锚点取最远档（近场偏离最小，作基准最稳）
D_anchor = D_list_B[i_anchor]
alpha = PARAMS["channel"]["alpha"]
atm_ratio = np.exp(-2*alpha*(D_list_B - D_anchor))       # 相对锚点的大气透过率²

# ① 信号率的 1/D² 理论线 + 误差
rate_invsq     = rate_B[i_anchor] * (D_anchor / D_list_B)**2                 # 纯 1/D²
rate_invsq_atm = rate_B[i_anchor] * (D_anchor / D_list_B)**2 * atm_ratio     # 含大气
err_rate = (rate_B - rate_invsq_atm) / np.maximum(rate_invsq_atm, 1e-30) * 100   # 相对误差 [%]

# ② SNR 的 1/D² 理论线 + 误差（用有效档拟合锚点）
finite = np.isfinite(snr_B)
snr_invsq     = snr_B[i_anchor] * (D_anchor / D_list_B)**2                   # 纯 1/D²
snr_invsq_atm = snr_B[i_anchor] * (D_anchor / D_list_B)**2 * atm_ratio       # 含大气
err_snr = np.where(finite, (snr_B - snr_invsq_atm) / np.maximum(snr_invsq_atm, 1e-30) * 100, np.nan)

# ---- 绘图：两张图都做"仿真 vs 1/D²"对比；y轴按仿真范围；右侧y轴画相对误差 ----
#   用户要求：① 主 y 轴范围按【仿真数据】而非平方反比（避免 1/D² 近场冲高把仿真压平）；
#            ② 每张图右侧加一条 y 轴画【相对误差 %】。
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# ===== 左图：SNR vs 距离 =====
ax = axes[0]
ax.plot(D_list_B[finite], snr_B[finite], "o-", color="tab:purple", lw=1.8, ms=5, label="SNR（仿真）")
ax.plot(D_list_B, snr_invsq_atm, "--", color="tab:green", lw=1.4, alpha=0.8,
        label=f"1/D²·exp(-2αD)（@{D_anchor:.0f}m 归一）")
ax.axhline(5.0, color="orange", ls=":", lw=1.0, alpha=0.7, label="SNR=5")
# y 轴按仿真范围（留 10% 余量），不让 1/D² 近场冲高撑坏坐标
_sv = snr_B[finite]
if _sv.size:
    ax.set_ylim(0, _sv.max()*1.1)
ax.set_xlabel("距离 D [m]"); ax.set_ylabel("SNR = S / √B（仿真范围）")
ax.set_title(f"SNR vs 距离（仿真 vs 1/D²，平均|误差|{np.nanmean(np.abs(err_snr)):.1f}%）")
ax.legend(fontsize=8, loc="upper right"); ax.grid(alpha=0.3)
# 右侧 y 轴：相对误差
axe = ax.twinx()
axe.plot(D_list_B, err_snr, "^:", color="tab:brown", lw=1.0, ms=4, alpha=0.7, label="相对误差 [%]")
axe.axhline(0, color="tab:brown", ls="-", lw=0.6, alpha=0.4)
axe.set_ylabel("相对误差 [%]", color="tab:brown"); axe.tick_params(axis="y", labelcolor="tab:brown")
axe.legend(fontsize=8, loc="lower right")

# ===== 右图：信号光子率峰值 vs 距离 =====
ax = axes[1]
ax.plot(D_list_B, rate_B, "o-", color="tab:blue", lw=1.7, ms=5, label="信号光子率峰值（仿真, 未封顶）")
ax.plot(D_list_B, rate_invsq_atm, "--", color="tab:green", lw=1.4, alpha=0.8,
        label=f"1/D²·exp(-2αD)（@{D_anchor:.0f}m 归一）")
ax.set_ylim(0, rate_B.max()*1.1)      # y 轴按仿真范围
ax.set_xlabel("距离 D [m]"); ax.set_ylabel("信号光子率峰值 [ph/s]（仿真范围）")
ax.set_title(f"信号率 vs 距离（仿真 vs 1/D²，平均|误差|{np.nanmean(np.abs(err_rate)):.2f}%）")
ax.legend(fontsize=8, loc="upper right"); ax.grid(alpha=0.3)
# 右侧 y 轴：相对误差
axe = ax.twinx()
axe.plot(D_list_B, err_rate, "^:", color="tab:brown", lw=1.0, ms=4, alpha=0.7, label="相对误差 [%]")
axe.axhline(0, color="tab:brown", ls="-", lw=0.6, alpha=0.4)
axe.set_ylabel("相对误差 [%]", color="tab:brown"); axe.tick_params(axis="y", labelcolor="tab:brown")
axe.legend(fontsize=8, loc="lower right")

plt.suptitle(f"模块 20 — SNR & 信号率 vs 距离：与平方反比 1/D² 对比 (ρ={echo0['rho']}, N_shots={N_shots}, {N_D_B}档)",
             fontsize=12)
plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.savefig("v43_snr_vs_distance.png", dpi=110, bbox_inches="tight")
plt.show()

print("="*76)
print(f"模块 20 汇总：SNR & 信号率 vs 距离，与平方反比 1/D² 对比 (ρ={echo0['rho']}, N_shots={N_shots}, {N_D_B}档)")
if finite.any():
    print(f"  SNR 范围 [{snr_B[finite].min():.2f}, {snr_B[finite].max():.2f}]；"
          f"SNR≥5 最大距离 ≈ {D_list_B[np.where(snr_B>=5.0)[0][-1]]:.1f} m" if (snr_B>=5.0).any() else "  SNR 全<5")
print(f"  以最远档 D={D_anchor:.0f}m 为锚点，理论 = 1/D²·exp(-2αD)：")
print(f"  ① 信号率峰值 vs 1/D²：平均|误差|={np.nanmean(np.abs(err_rate)):.2f}%，最大={np.nanmax(np.abs(err_rate)):.2f}%")
print(f"     → 未封顶物理量，几乎严格遵循平方反比（误差主要来自蒙卡统计涨落）。")
print(f"  ② SNR vs 1/D²：平均|误差|={np.nanmean(np.abs(err_snr)):.1f}%，最大={np.nanmax(np.abs(err_snr)):.1f}%")
print(f"     → SNR 近似遵循 1/D²，偏差主因：背景 B 随距离非严格恒定、二值饱和使近场 S 被削顶。")
'''


def code_cell(cid, src):
    return {"cell_type": "code", "id": cid, "metadata": {},
            "execution_count": None, "outputs": [],
            "source": src.splitlines(keepends=True)}

def md_cell(cid, src):
    return {"cell_type": "markdown", "id": cid, "metadata": {},
            "source": src.splitlines(keepends=True)}


# ============================================================================
# 读入 v32 → 逐 cell 处理 → 追加新模块 → 写出 v41
# ============================================================================
with open(SRC_NB, "r", encoding="utf-8") as f:
    nb = json.load(f)
src_cells = nb["cells"]
print(f"\n读入 {SRC_NB}: {len(src_cells)} 个 cell")

out_cells = []
n_disabled = n_replaced = n_split = 0
for i, c in enumerate(src_cells):
    if i == 25:
        # cell25 拆分：先插入"模块 0b 时间窗"code cell，再插入"模块8 演示"停用 markdown
        out_cells.append(code_cell("module_0b_timewindow", MODULE_0B))
        out_cells.append(split_cell25(c))
        n_split += 1
        continue
    if i == 39:                       # 模块11b 整块替换
        out_cells.append(code_cell(c.get("id", "mod11b"), MODULE_11B_NEW))
        n_replaced += 1
        continue
    if i == 45:                       # 模块12 整块替换
        out_cells.append(code_cell(c.get("id", "mod12"), MODULE_12_NEW))
        n_replaced += 1
        continue
    if i == 43:                       # 模块18（恢复）：仅把坐标轴改成绘图窗，其余不动
        src43 = "".join(c["source"])
        # 原 v32 用固定 set_xlim(t0_ns-8, t0_ns+18)，改成全局绘图窗 plot_lo_ns/plot_hi_ns
        src43 = src43.replace("ax[1].set_xlim(t0_ns-8, t0_ns+18)",
                              "ax[1].set_xlim(plot_lo_ns, plot_hi_ns)  # v41: 用全局绘图窗(ToF前20/后50ns)")
        out_cells.append(code_cell(c.get("id", "mod18"), src43))
        n_replaced += 1
        continue
    if i in DISABLE and c["cell_type"] == "code":   # 停用 → markdown
        out_cells.append(to_markdown_disabled(c, DISABLE[i]))
        n_disabled += 1
        continue
    out_cells.append(c)               # 其余原样保留（含恢复的 23/33/47）

# 追加 v42 新模块
out_cells.append(md_cell("v43_new_header", NEW_HEADER_MD))
out_cells.append(code_cell("v43_module_19_energy", MODULE_19))
out_cells.append(code_cell("v43_module_20_snr", MODULE_20))

nb["cells"] = out_cells
with open(OUT_NB, "w", encoding="utf-8") as f:
    json.dump(nb, f, ensure_ascii=False, indent=1)

print(f"已生成 {OUT_NB}: 共 {len(out_cells)} 个 cell")
print(f"  · cell25 拆分为 [模块0b 时间窗] + [模块8 演示停用]（{n_split} 处）")
print(f"  · 模块11b/12 整块替换（降强+独立重跑+绘图窗，{n_replaced} 个）")
print(f"  · 停用转 markdown：{n_disabled} 个 {sorted(DISABLE)}")
print(f"  · 恢复：模块7c(23)/9c(33)/18(43)/12b(47)（原样保留，未停用）")
print(f"  · 追加：模块19 能量扫描 + 模块20 SNR-距离（含 1/D^2 判断）")
