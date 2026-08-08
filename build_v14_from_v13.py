# -*- coding: utf-8 -*-
"""
基于 v13 生成 v14:
任务1: 新增模块 6b —— SPAD 阵列俯视图(网格 + 光斑 1/e² 椭圆轮廓 + 按 f_pix 涂色)。
任务2: 新增模块 8b —— RC 情况下单个 SPAD 某次蒙卡采样的过电压 Vov(t) 曲线, 标阈值。
任务3: RC 参数/死时间/阈值等移入模块0 PARAMS["spad"], 模块7b/8 改为从字典读取。
运行: python build_v14_from_v13.py
"""
import nbformat

SRC = "lidar_histogram_sim_v13.ipynb"
DST = "lidar_histogram_sim_v14.ipynb"

nb = nbformat.read(SRC, as_version=4)

# 清空输出/执行序号 + 去执行时间戳
for c in nb.cells:
    if c.cell_type == "code":
        c.outputs = []
        c.execution_count = None
    if "execution" in c.get("metadata", {}):
        del c.metadata["execution"]

def get_cell(pred):
    for i, c in enumerate(nb.cells):
        if pred(c):
            return i, c
    raise RuntimeError("cell not found")

# ================= 任务3a: 模块0 参数字典新增 SPAD 参数 =================
i0c, c0 = get_cell(lambda c: c.cell_type=="code" and '"spad"' in c.source and "spad_array" in c.source)
old_spad = '''    # ---------- SPAD 器件 ----------
    "spad": {
        "PDE": 0.30,            # 光子探测效率
        "DCR": 0.0e3,           # 单像元暗计数率 [cps]
        "jitter_sigma": 100e-12,# IRF 高斯 sigma [s]
    },'''
new_spad = '''    # ---------- SPAD 器件 ----------
    "spad": {
        "PDE": 0.30,            # 光子探测效率(=PDE_max, 满过电压时的峰值 PDE)
        "DCR": 0.0e3,           # 单像元暗计数率 [cps]
        "jitter_sigma": 100e-12,# IRF 高斯 sigma [s]
        # --- 死时间 / 恢复模型参数(v14: 从此处统一配置) ---
        "t_dead": 14e-9,        # 硬死时间 [s] (hard dead time 对照模型, non-paralyzable)
        "tau_rc": 6e-9,         # RC 恢复时间常数 τ=R·C_J [s] (RC 模型)
        "Vov_max": 3.3,         # 满过电压 Vov_max [V] (=V_bias − V_br)
        "Vth_frac": 0.10,       # 计数所需最小过电压占比 (Vth = Vth_frac·Vov_max)
        "reset_mode": "count",  # RC 复位方式: "count"=仅计数事件复位 / "all"=任何雪崩都复位
    },'''
assert old_spad in c0.source, "模块0 spad 字典锚点未匹配"
c0.source = c0.source.replace(old_spad, new_spad)

# ================= 任务3b: 模块7b 改为从字典读取参数 =================
i7b, c7b = get_cell(lambda c: c.cell_type=="code" and "def simulate_spad_shot_rc(" in c.source)
old_rcparam = '''# ---- RC 参数(本版设定) ----
TAU_RC   = 6e-9        # RC 恢复时间常数 = R·C_J
VTH_FRAC = 0.10        # 计数所需最小 Vov (占 Vov_max)
VOV_MAX  = 3.3         # 满过电压 [V]
RESET    = "count"     # 仅计数事件复位

# 等效"硬死区"(恢复到 Vth 前完全不计数)与"渐变灵敏"分界
t_deadzone = -np.log(1 - VTH_FRAC) * TAU_RC
print(f"RC 引擎就绪: τ_RC={TAU_RC*1e9:.1f} ns, Vth={VTH_FRAC*100:.0f}%·Vov_max, reset='{RESET}'")
print(f"  Vov 恢复曲线: 1τ->{100*(1-np.exp(-1)):.0f}%, 2.3τ->{100*(1-np.exp(-2.3)):.0f}%, 5τ->{100*(1-np.exp(-5)):.1f}%")
print(f"  低于 Vth 的'硬死区'≈{t_deadzone*1e9:.2f} ns, 之后为渐变灵敏(与硬死时间 14ns 的一刀切不同)")'''

new_rcparam = '''def simulate_spad_shot_rc_trace(r_sig_fine, r_amb_ph, tf, PDE_max, tau_rc, Vth_frac,
                                jitter_sigma, rng, reset_mode="count", Vov_max=3.3):
    """与 simulate_spad_shot_rc 同逻辑, 但额外记录 Vov(t) 轨迹供可视化(模块 8b 用)。
    返回 dict: t_photon(所有入射光子时刻), fired(是否触发雪崩), counted(是否计数),
    vov_at_photon(光子到达瞬间的 Vov 比例), t_count(计数事件时刻), t_grid/vov_grid(连续 Vov 曲线)。"""
    dt = tf[1] - tf[0]
    mu = (r_sig_fine + r_amb_ph) * dt
    n_ph = rng.poisson(mu)
    t_arr = np.repeat(tf, n_ph) if n_ph.sum() else np.empty(0)
    u = rng.random(t_arr.size)
    fired = np.zeros(t_arr.size, dtype=bool)
    counted = np.zeros(t_arr.size, dtype=bool)
    vov_at = np.zeros(t_arr.size)
    resets = []                                  # 复位时刻(用于重建连续曲线)
    last = -1e30; inv_tau = 1.0 / tau_rc
    for k in range(t_arr.size):
        t = t_arr[k]; d = (t - last) * inv_tau
        vf = 1.0 - np.exp(-d) if d < 700 else 1.0
        vov_at[k] = vf
        if u[k] < PDE_max * vf:
            fired[k] = True
            if vf >= Vth_frac:
                counted[k] = True; last = t; resets.append(t)
            elif reset_mode == "all":
                last = t; resets.append(t)
    # 连续 Vov(t) 曲线: 在精细网格上按"距上次复位"重建
    tg = tf.copy()
    last = -1e30; vg = np.ones(tg.size)
    ri = 0; resets = np.asarray(resets)
    for idx in range(tg.size):
        t = tg[idx]
        while ri < resets.size and resets[ri] <= t:
            last = resets[ri]; ri += 1
        d = (t - last) * inv_tau
        vg[idx] = 1.0 - np.exp(-d) if d < 700 else 1.0
    t_count = t_arr[counted]
    return {"t_photon": t_arr, "fired": fired, "counted": counted, "vov_at_photon": vov_at,
            "t_count": t_count, "t_grid": tg, "vov_grid": vg}

# ---- RC / 死时间参数: 统一从 PARAMS["spad"] 读取(v14) ----
_sp = PARAMS["spad"]
PDE      = _sp["PDE"]            # = PDE_max
JIT      = _sp["jitter_sigma"]
T_DEAD   = _sp["t_dead"]         # 硬死时间(对照)
TAU_RC   = _sp["tau_rc"]         # RC 恢复时间常数
VOV_MAX  = _sp["Vov_max"]
VTH_FRAC = _sp["Vth_frac"]
RESET    = _sp["reset_mode"]
jit = JIT; t_dead = T_DEAD       # 小写别名(下游模块 8/9/11/12 沿用旧命名)

# 等效"硬死区"(恢复到 Vth 前完全不计数)与"渐变灵敏"分界
t_deadzone = -np.log(1 - VTH_FRAC) * TAU_RC
print(f"参数来自 PARAMS['spad']: PDE_max={PDE}, 硬死时间={T_DEAD*1e9:.1f}ns, τ_RC={TAU_RC*1e9:.1f}ns, "
      f"Vov_max={VOV_MAX}V, Vth={VTH_FRAC*100:.0f}%, reset='{RESET}'")
print(f"  Vov 恢复曲线: 1τ->{100*(1-np.exp(-1)):.0f}%, 2.3τ->{100*(1-np.exp(-2.3)):.0f}%, 5τ->{100*(1-np.exp(-5)):.1f}%")
print(f"  低于 Vth 的'硬死区'≈{t_deadzone*1e9:.2f} ns, 之后为渐变灵敏(与硬死时间一刀切不同)")'''

assert old_rcparam in c7b.source, "模块7b RC 参数锚点未匹配"
c7b.source = c7b.source.replace(old_rcparam, new_rcparam)

# 模块7b markdown 里补一句 trace 函数说明
i7bmd, c7bmd = get_cell(lambda c: c.cell_type=="markdown" and c.source.lstrip().startswith("## 模块 7b"))
c7bmd.source += '''

> v14: 参数(PDE、硬死时间、τ_RC、Vov_max、Vth、reset)统一从 `PARAMS["spad"]` 读取;
> 另附 `simulate_spad_shot_rc_trace()` 记录 Vov(t) 轨迹, 供模块 8b 画过电压曲线。'''

# ================= 任务3c: 模块8 改为从字典读取(删去硬编码 t_dead/jit) =================
i8, c8 = get_cell(lambda c: c.cell_type=="code" and "t_dead = 14e-9" in c.source)
old8 = '''PDE = PARAMS["spad"]["PDE"]; jit = PARAMS["spad"]["jitter_sigma"]
t_dead = 14e-9                                # v12 硬死时间(作对照)'''
new8 = '''# 参数已在模块 7b 从 PARAMS["spad"] 读入: PDE, jit, t_dead, TAU_RC, VTH_FRAC, VOV_MAX, RESET'''
assert old8 in c8.source, "模块8 参数锚点未匹配"
c8.source = c8.source.replace(old8, new8)

# 模块8/9/11/12 里若把 TAU_RC*1e9=6 硬写在标题, 保持不变(它们引用变量, 无碍)

# ================= 任务1: 模块6b —— 阵列俯视图(网格+光斑+涂色)=================
i6, c6 = get_cell(lambda c: c.cell_type=="code" and "from scipy.special import erf" in c.source)

md_6b = nbformat.v4.new_markdown_cell(r'''## 模块 6b（v14 新增）— SPAD 阵列俯视图：网格 + 光斑轮廓 + 按收集比例涂色

在像面(SPAD 阵列平面)上直观展示回波像斑落在阵列上的样子:
- **方格**: 9×120 物理 SPAD 阵列(pitch=10µm), 每格一个 SPAD;
- **光斑轮廓**: 30m 目标回波的椭圆高斯像斑, 画出 1/e²(全宽 s)与 1σ(=s/4)两条椭圆;
- **涂色**: 每个 SPAD 按其空间收集比例 `f_pix` 着色(颜色越亮=接收信号越多)。

> x=短边(9 个 SPAD), y=长边(120 个 SPAD); 像斑长边 σ_y≈200µm≈20 像元, 短边 σ_x≈11µm≈1 像元。''')

code_6b = nbformat.v4.new_code_cell(r'''from matplotlib.patches import Ellipse, Rectangle
from matplotlib.collections import PatchCollection

# 阵列几何(以阵列中心为原点, 单位 µm)
a = PARAMS["spad_array"]; pitch_um = a["pitch"]*1e6
Nx, Ny = a["Nx"], a["Ny"]
xi_um = (np.arange(Nx) - (Nx-1)/2.0) * pitch_um     # 各 SPAD 中心 x
yj_um = (np.arange(Ny) - (Ny-1)/2.0) * pitch_um     # 各 SPAD 中心 y

# 像斑参数(30m 主目标)
sx_um, sy_um = sx0*1e6, sy0*1e6                       # 1/e² 全宽 [µm]
sig_x_um, sig_y_um = sx_um/4.0, sy_um/4.0            # 1σ [µm]

# 涂色: 归一化 f_pix
fnorm = fpix0 / fpix0.max()
cmap = plt.cm.viridis

fig, axes = plt.subplots(1, 2, figsize=(13, 6.2),
                         gridspec_kw={"width_ratios":[1, 2.6]})

# ---- 左: 全阵列(y 方向很长) ----
def draw_array(ax, show_all_grid):
    rects = []; colors = []
    for ii in range(Nx):
        for jj in range(Ny):
            x = xi_um[ii]-pitch_um/2; y = yj_um[jj]-pitch_um/2
            rects.append(Rectangle((x, y), pitch_um, pitch_um))
            colors.append(cmap(fnorm[ii, jj]))
    pc = PatchCollection(rects, facecolor=colors, edgecolor=(0,0,0,0.15), linewidth=0.2)
    ax.add_collection(pc)
    # 光斑椭圆(中心在原点)
    for k, ls, lab in [(1.0,"-","光斑 1/e² 全宽"), (0.5,"--","光斑 1σ")]:
        # 1/e² 全宽 = 2*(s/2)=s => 半轴 = s/2; 1σ 半轴 = sig
        wx = (sx_um if k==1.0 else 2*sig_x_um); wy = (sy_um if k==1.0 else 2*sig_y_um)
        ax.add_patch(Ellipse((0,0), wx, wy, fill=False, edgecolor="red",
                             ls=ls, lw=1.6, label=lab))
    ax.scatter([xi_um[i0]], [yj_um[j0]], c="white", edgecolor="red", marker="o", s=60,
               zorder=5, label=f"中心 SPAD ({i0},{j0})")

draw_array(axes[0], True)
axes[0].set_xlim(xi_um[0]-pitch_um, xi_um[-1]+pitch_um)
axes[0].set_ylim(yj_um[0]-pitch_um, yj_um[-1]+pitch_um)
axes[0].set_aspect("equal")
axes[0].set_xlabel("x [µm] (短边, 9 SPAD)"); axes[0].set_ylabel("y [µm] (长边, 120 SPAD)")
axes[0].set_title("全阵列 9×120\n(y 方向真实比例)")
axes[0].legend(fontsize=7, loc="upper right")

# ---- 右: 放大被照区(y 方向裁剪到光斑附近) ----
draw_array(axes[1], False)
y_half = 2.2*sig_y_um
axes[1].set_xlim(xi_um[0]-pitch_um, xi_um[-1]+pitch_um)
axes[1].set_ylim(-y_half, y_half)
axes[1].set_aspect("equal")
axes[1].set_xlabel("x [µm] (短边)"); axes[1].set_ylabel("y [µm] (长边)")
axes[1].set_title(f"被照区放大 (±2.2σ_y)  像斑 1/e²: {sx_um:.0f}×{sy_um:.0f}µm, σ: {sig_x_um:.0f}×{sig_y_um:.0f}µm")
axes[1].legend(fontsize=8, loc="upper right")

sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(0, fpix0.max()))
plt.colorbar(sm, ax=axes[1], label="SPAD 收集比例 f_pix", fraction=0.046, pad=0.04)
plt.tight_layout(); plt.show()

n_lit = int((fpix0 > fpix0.max()*0.01).sum())
print(f"光斑覆盖: 被照 SPAD(f_pix>1%峰值) 约 {n_lit} 个 / 共 {Nx*Ny} 个")
print(f"像斑长边 σ_y={sig_y_um:.0f}µm ≈ {sig_y_um/pitch_um:.0f} 像元; 短边 σ_x={sig_x_um:.0f}µm ≈ {sig_x_um/pitch_um:.0f} 像元")
print(f"中心 SPAD ({i0},{j0}) f_pix={fpix0[i0,j0]:.3e}; 涂色越亮=收集越多。")''')

# 插入到模块6 code cell(i6)之后
nb.cells.insert(i6+1, md_6b)
nb.cells.insert(i6+2, code_6b)

# ================= 任务2: 模块8b —— 单 SPAD 某次蒙卡的 Vov(t) 曲线 =================
# 找模块8 code cell(此时索引已因插入 6b 而+2)
i8, c8 = get_cell(lambda c: c.cell_type=="code" and "ev1_rc  = simulate_spad_shot_rc(" in c.source)

md_8b = nbformat.v4.new_markdown_cell(r'''## 模块 8b（v14 新增）— RC 模型：单个 SPAD 某次蒙卡的过电压 Vov(t) 曲线

画出**中心 SPAD** 在**某一次蒙卡采样**里, 过电压 `Vov(t)` 随时间的演化(这正是 RC 恢复模型的核心内部状态):
- **蓝线**: Vov(t)/Vov_max —— 每次计数事件后跌到 0, 再按 `1−e^(−Δt/τ)` 指数充回;
- **红色虚线**: 计数阈值 `Vth = Vth_frac·Vov_max`(只有 Vov 恢复到此线以上, 光子触发才被计数);
- **绿色 ↑**: 被**计数**的光子(触发且 Vov≥Vth); **灰色 ×**: 触发但**亚阈**(不计数); **淡点**: 未触发的入射光子。

> 展示"恢复期渐变灵敏": 阈值线以下也可能触发(灰×), 但幅度不够不计数; 阈值线以上才计数(绿↑)。''')

code_8b = nbformat.v4.new_code_cell(r'''# 选一次"事件较丰富"的蒙卡采样(多试几个 seed, 取计数事件多的, 便于观察)
best = None
for s in range(20):
    rng_t = np.random.default_rng(1000 + s)
    tr = simulate_spad_shot_rc_trace(r_sig, r_amb_ph, tf, PDE, TAU_RC, VTH_FRAC,
                                     0.0, rng_t, RESET, VOV_MAX)   # jitter=0: Vov 曲线用雪崩时刻
    nc = int(tr["counted"].sum())
    if best is None or nc > best[0]:
        best = (nc, s, tr)
n_cnt, seed_used, tr = best

tph_ns = tr["t_photon"]*1e9
fired = tr["fired"]; counted = tr["counted"]; vov_ph = tr["vov_at_photon"]
tg_ns = tr["t_grid"]*1e9; vg = tr["vov_grid"]
Vth = VTH_FRAC

fig, ax = plt.subplots(figsize=(11.5, 5.2))
# Vov(t) 连续曲线
ax.plot(tg_ns, vg, color="tab:blue", lw=1.4, label="Vov(t)/Vov_max (RC 恢复)")
# 阈值线
ax.axhline(Vth, color="red", ls="--", lw=1.4, label=f"计数阈值 Vth = {VTH_FRAC:.0%}·Vov_max = {VTH_FRAC*VOV_MAX:.2f} V")
# ToF
ax.axvline(t0_ns, color="k", ls=":", alpha=0.6, label=f"真实 ToF {t0_ns:.1f} ns")

# 光子标记: 分三类
m_cnt = counted
m_sub = fired & (~counted)             # 触发但亚阈
m_no  = (~fired)                       # 未触发
ax.scatter(tph_ns[m_no],  vov_ph[m_no],  s=14, c="0.6", alpha=0.35, marker=".", label="入射光子(未触发)")
ax.scatter(tph_ns[m_sub], vov_ph[m_sub], s=55, c="0.35",
           marker="x", linewidth=1.4, label="触发但亚阈(不计数)")
ax.scatter(tph_ns[m_cnt], vov_ph[m_cnt], s=90, c="tab:green", marker="^",
           edgecolor="k", linewidth=0.5, zorder=6, label="计数事件(触发且 Vov≥Vth)")

ax.set_ylim(-0.03, 1.05)
ax.set_xlim(t_lo*1e9, t_hi*1e9)
ax.set_xlabel("时间 t [ns]"); ax.set_ylabel("Vov / Vov_max")
ax.set_title(f"单个 SPAD 某次蒙卡采样的过电压 Vov(t)  (中心 SPAD, τ_RC={TAU_RC*1e9:.0f}ns, "
             f"Vth={VTH_FRAC:.0%}, 本次计数 {n_cnt} 个)")
ax.legend(fontsize=8, loc="center right", ncol=1)
ax.grid(alpha=0.3)
plt.tight_layout(); plt.show()

n_ph_tot = tr["t_photon"].size; n_fire = int(fired.sum())
print(f"本次蒙卡采样(seed={1000+seed_used}): 入射光子 {n_ph_tot} 个, 触发雪崩 {n_fire} 个, 计数 {n_cnt} 个")
print(f"  其中 {n_fire-n_cnt} 个是'触发但亚阈'(发生在 Vov 未恢复到 {VTH_FRAC:.0%} 时, 幅度不够不计数)")
print(f"  Vth = {VTH_FRAC:.0%}·Vov_max = {VTH_FRAC*VOV_MAX:.2f} V; 每次计数后 Vov 跌 0, 按 τ={TAU_RC*1e9:.0f}ns 指数充回")
print(f"  注: 此图 Vov 曲线用雪崩时刻(未加 IRF 抖动), 抖动只影响记录的 timestamp。")''')

# 插入到模块8 code cell 之后(即模块9 之前)
nb.cells.insert(i8+1, md_8b)
nb.cells.insert(i8+2, code_8b)

# ================= 标题更新 =================
if nb.cells and nb.cells[0].cell_type == "markdown":
    nb.cells[0].source = nb.cells[0].source.replace(
        "# 激光雷达直方图仿真 v13 (LiDAR Histogram Simulation) — RC 恢复 SPAD 模型",
        "# 激光雷达直方图仿真 v14 (LiDAR Histogram Simulation) — 阵列俯视图 + Vov 曲线 + 参数化")
    nb.cells[0].source += '''

**v14 相对 v13 的变化**
- 新增**模块 6b**: SPAD 阵列俯视图(网格 + 光斑 1/e²/1σ 椭圆轮廓 + 每个 SPAD 按 f_pix 涂色)。
- 新增**模块 8b**: RC 模型下单个 SPAD 某次蒙卡采样的**过电压 Vov(t) 曲线**, 标出计数阈值 Vth。
- **参数化**: 硬死时间、τ_RC、Vov_max、Vth、reset 全部移入 `PARAMS["spad"]`, 模块 7b/8 从字典读取。'''

nbformat.write(nb, DST)
print(f"已生成 {DST}: 共 {len(nb.cells)} 个 cell。")
