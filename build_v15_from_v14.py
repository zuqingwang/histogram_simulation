# -*- coding: utf-8 -*-
"""
基于 v14 生成 v15:
新增"目标深度展宽"建模(SPAD 架构下 10ns 上升沿的唯一物理来源):
  1) 连续深度分布 Δz: 回波 = 激光脉冲 ⊗ 深度核(均匀/高斯) ⊗ IRF;
  2) 离散多子目标: 深度范围内多个物体, 回波 = 各子目标(各自延迟+相对反射率)叠加。
明确不建"电子学带宽展宽"(SPAD+TDC 只记录雪崩时刻, 后级带宽不影响时间戳)。

实现方式:
  - PARAMS 保持 30m 主目标不变(不擅改用户场景); 新增字段仅作可选、缺省即 v14 行为;
  - 模块 6c 新增: 深度→时间波形 的推导说明 + 深度展宽核函数 depth_broadening_kernel();
  - 扩展 signal_photon_rate_fine: 支持 echo["depth_span"]/["depth_profile"]/["sub_targets"];
  - 模块 8c 新增: 演示不同 Δz / 多子目标下回波上升沿变化(局部构造 echo, 不改全局参数)。
运行: python build_v15_from_v14.py
"""
import nbformat

SRC = "lidar_histogram_sim_v14.ipynb"
DST = "lidar_histogram_sim_v15.ipynb"
nb = nbformat.read(SRC, as_version=4)

for c in nb.cells:
    if c.cell_type == "code":
        c.outputs = []; c.execution_count = None
    if "execution" in c.get("metadata", {}):
        del c.metadata["execution"]

def find(pred):
    for i, c in enumerate(nb.cells):
        if pred(c): return i, c
    raise RuntimeError("not found")

# ========== 1. 参数字典: 给 30m echo 增加"深度展宽"可选字段(缺省不改变行为) ==========
i0c, c0 = find(lambda c: c.cell_type=="code" and '"echoes"' in c.source)
old_echo = '''            {"D": 30.0,  "rho": 0.10, "frac": 1.00, "tilt_deg": 0.0,  "name": "近处物体(边缘拦截)"},'''
new_echo = '''            {"D": 30.0,  "rho": 0.10, "frac": 1.00, "tilt_deg": 0.0,  "name": "近处物体(边缘拦截)",
             # v15 深度展宽(可选, 缺省=0 即无展宽, 与 v14 一致):
             "depth_span": 0.0,        # 连续深度范围 Δz [m] (沿视线); 0=不展宽
             "depth_profile": "uniform",  # "uniform" 均匀 / "gauss" 高斯
             "sub_targets": []},       # 离散子目标: [{"dD":相对主距偏移[m], "rho_rel":相对反射率}, ...]'''
assert old_echo in c0.source, "echo 锚点未匹配"
c0.source = c0.source.replace(old_echo, new_echo)

# ========== 2. 模块 6c: 深度展宽核函数 + 推导说明(插在模块6之后, 模块6b之前) ==========
# 找模块6 code cell(erf 定义处)
i6, _ = find(lambda c: c.cell_type=="code" and "from scipy.special import erf" in c.source)

md_6c = nbformat.v4.new_markdown_cell(r'''## 模块 6c（v15 新增）— 目标深度展宽：回波上升沿的物理来源

**问题**: 实际 SPAD 激光雷达回波上升沿可宽达 ~10ns, 而激光脉冲固有上升沿仅 ~0.6ns, 差一个数量级。

**盘点时间展宽机制**(探测回波 = 各机制波形的卷积):
| 机制 | 来源 | 量级 | 模型 |
|---|---|---|---|
| 激光脉冲固有波形 | 双指数上升 τ_r=0.7ns | 上升沿 0.62ns | ✅ |
| SPAD IRF 抖动 | 雪崩建立统计涨落 | FWHM 0.24ns | ✅ |
| **目标深度展宽** | 光斑内目标沿视线的深度分布 Δz | **2Δz/c**: 1.5m→10ns | ★ v15 新增 |
| ~~电子学带宽~~ | ~~TIA/放大器带宽~~ | ~~t_rise≈0.35/BW~~ | ✗ 不建* |

> *SPAD+TDC(Time-to-Digital Converter, 时间数字转换器)是**数字式**探测: TDC 只记录**雪崩过阈时刻**,
> 后级模拟带宽把电脉冲幅度波形拖长, 但**不改变记录的时间戳**。故 SPAD 架构下电子学带宽**不产生**上升沿展宽,
> 时间分辨只受 IRF 限制。因此本模型**唯一**能产生 10ns 上升沿的物理机制是**目标深度展宽**。

**深度→时间映射**: 深度 z 处散射体回波延迟 `t = 2z/c`(双程)。
深度分布 g(z) 按 `dt = 2dz/c` 映射为时间展宽核; **均匀深度 Δz 的时间全宽 = 2Δz/c**。
- Δz=15cm → 1.0ns;  Δz=75cm → 5.0ns;  **Δz=150cm → 10.0ns**(与"10ns↔150cm"一致)。''')

code_6c = nbformat.v4.new_code_cell(r'''def depth_broadening_kernel(depth_span, profile, dt_fine):
    """目标深度展宽的时间核。depth_span=沿视线深度范围 Δz [m]; profile='uniform'/'gauss'。
    返回归一化核(积分×dt=1), 与信号率卷积时用 mode='same' 再 ×dt。"""
    if depth_span <= 0:
        return np.array([1.0/dt_fine])          # δ 核(不展宽)
    T_span = 2.0 * depth_span / C_LIGHT          # 时间全宽 = 2Δz/c(双程)
    if profile == "gauss":
        # 把 depth_span 视为高斯 1σ 对应的深度 -> 时间 σ
        sig_t = T_span                            # 用 2Δz/c 作为 1σ(可调约定)
        return gaussian_kernel(sig_t, dt_fine)
    # 默认 uniform: 时间上宽 T_span 的矩形
    n = max(1, int(round(T_span / dt_fine)))
    k = np.ones(n)
    return k / (k.sum() * dt_fine)

# 快速自检: 均匀深度核的时间全宽
for dz in [0.15, 0.75, 1.5]:
    ker = depth_broadening_kernel(dz, "uniform", PARAMS["hist"]["dt_fine"])
    W = ker.size * PARAMS["hist"]["dt_fine"] * 1e9
    print(f"深度 Δz={dz*100:.0f}cm -> 均匀时间核全宽 = {W:.2f} ns (= 2Δz/c)")
print("(深度展宽核将在 signal_photon_rate_fine 中与激光脉冲卷积, 形成宽上升沿)")''')

nb.cells.insert(i6+1, md_6c)
nb.cells.insert(i6+2, code_6c)

# ========== 3. 扩展 signal_photon_rate_fine: 深度核 + 多子目标 ==========
i_sig, c_sig = find(lambda c: c.cell_type=="code" and "def signal_photon_rate_fine" in c.source)
old_sig = '''def signal_photon_rate_fine(echo, f_pix_ij, tf, p=PARAMS):
    """单 SPAD 信号【光子到达率】(不含 PDE)在精细网格 tf 上。
    已含倾角几何展宽, 未含 IRF; f_pix_ij 为该像元空间收集比例(标量)。"""
    t0 = time_of_flight(echo["D"])
    r = pulse_temporal(tf - t0, p) * link_factor(echo, p) / E_PHOTON * f_pix_ij
    sig_b = echo_range_broadening_sigma(echo["D"], echo["tilt_deg"], p)   # 几何展宽(非 IRF)
    if sig_b > 0:
        dt_fine = tf[1] - tf[0]
        r = np.convolve(r, gaussian_kernel(sig_b, dt_fine), mode="same") * dt_fine
    return r'''

new_sig = '''def _single_echo_rate(echo, f_pix_ij, tf, p, link_scale=1.0):
    """单个(子)回波的光子率: 双指数脉冲 × 链路 × 收集比例, 已含倾角几何展宽(高斯)。"""
    t0 = time_of_flight(echo["D"])
    r = pulse_temporal(tf - t0, p) * link_factor(echo, p) * link_scale / E_PHOTON * f_pix_ij
    sig_b = echo_range_broadening_sigma(echo["D"], echo.get("tilt_deg", 0.0), p)  # 倾角几何展宽
    if sig_b > 0:
        dt_fine = tf[1] - tf[0]
        r = np.convolve(r, gaussian_kernel(sig_b, dt_fine), mode="same") * dt_fine
    return r

def signal_photon_rate_fine(echo, f_pix_ij, tf, p=PARAMS):
    """单 SPAD 信号【光子到达率】(不含 PDE)在精细网格 tf 上。
    v15: 在 v14(倾角几何展宽)基础上, 叠加两类目标深度展宽——
      (a) 连续深度分布 echo["depth_span"](Δz)/["depth_profile"]: 与深度核卷积;
      (b) 离散子目标 echo["sub_targets"]: 各自延迟(dD)+相对反射率(rho_rel)叠加。
    所有新字段缺省时(depth_span=0, sub_targets=[]) 行为与 v14 完全一致。"""
    dt_fine = tf[1] - tf[0]
    # 主回波
    r = _single_echo_rate(echo, f_pix_ij, tf, p)
    # (b) 离散子目标叠加(相对主距离偏移 dD, 相对反射率 rho_rel)
    for sub in echo.get("sub_targets", []):
        sub_echo = dict(echo)
        sub_echo["D"] = echo["D"] + sub["dD"]
        sub_echo["sub_targets"] = []; sub_echo["depth_span"] = 0.0
        r = r + _single_echo_rate(sub_echo, f_pix_ij, tf, p, link_scale=sub.get("rho_rel", 1.0))
    # (a) 连续深度展宽核卷积
    dz = echo.get("depth_span", 0.0)
    if dz > 0:
        ker = depth_broadening_kernel(dz, echo.get("depth_profile", "uniform"), dt_fine)
        if ker.size > 1:
            r = np.convolve(r, ker, mode="same") * dt_fine
    return r'''
assert old_sig in c_sig.source, "signal_photon_rate_fine 锚点未匹配"
c_sig.source = c_sig.source.replace(old_sig, new_sig)

# ========== 4. 模块 8c: 深度展宽对回波上升沿的演示(局部 echo, 不改全局) ==========
# 插在模块 8b code cell 之后
i8b, _ = find(lambda c: c.cell_type=="code" and "simulate_spad_shot_rc_trace(" in c.source and "best" in c.source)

md_8c = nbformat.v4.new_markdown_cell(r'''## 模块 8c（v15 新增）— 目标深度展宽对回波上升沿的影响

用**局部构造的 echo**(不改全局 30m 场景参数)演示深度展宽如何把回波上升沿从 ~0.6ns 拉宽到 ~10ns:
- **图 1（连续深度 Δz）**: 均匀深度 Δz=0/15/75/150cm 时的回波波形(信号光子率), 标注上升沿(10-90%);
- **图 2（离散多子目标）**: 1.5m 深度内 3 个物体(不同距离/反射率)合成的复合回波。

> Δz=150cm 对应时间全宽 2Δz/c=10ns —— 即"1.5m 深度内有其他物体, 回波就展宽到 ~10ns"。''')

code_8c = nbformat.v4.new_code_cell(r'''import copy
# 精细网格覆盖 30m 回波附近, 稍宽以容纳深度展宽
t0_30 = time_of_flight(30.0)
tf_d = np.arange(t0_30 - 5e-9, t0_30 + 25e-9, PARAMS["hist"]["dt_fine"])
tfd_ns = tf_d*1e9; t0_30_ns = t0_30*1e9
f_demo = fpix0[i0, j0]

def rise_1090(r, t):
    if r.max() <= 0: return np.nan, np.nan, np.nan
    rn = r/r.max(); pk = rn.argmax()
    tlo = t[:pk+1][np.searchsorted(rn[:pk+1], 0.1)]
    thi = t[:pk+1][np.searchsorted(rn[:pk+1], 0.9)]
    return (thi-tlo)*1e9, tlo*1e9, thi*1e9

base_echo = next(e for e in PARAMS["target"]["echoes"] if e["D"]==30.0)

# ---- 图 1: 连续深度 Δz ----
fig, ax = plt.subplots(1, 2, figsize=(13, 4.8))
for dz_cm in [0, 15, 75, 150]:
    e = copy.deepcopy(base_echo); e["depth_span"] = dz_cm/100.0; e["depth_profile"]="uniform"; e["sub_targets"]=[]
    r = signal_photon_rate_fine(e, f_demo, tf_d)
    rise, _, _ = rise_1090(r, tf_d)
    ax[0].plot(tfd_ns, r/max(r.max(),1e-30), lw=1.6,
               label=f"Δz={dz_cm}cm (2Δz/c={2*(dz_cm/100)/C_LIGHT*1e9:.1f}ns, 上升沿{rise:.1f}ns)")
ax[0].axvline(t0_30_ns, color="k", ls=":", alpha=0.6, label=f"主目标 ToF {t0_30_ns:.1f} ns")
ax[0].set_xlabel("时间 t [ns]"); ax[0].set_ylabel("归一化信号光子率")
ax[0].set_title("连续深度展宽: Δz 越大, 上升沿越宽")
ax[0].legend(fontsize=8); ax[0].grid(alpha=0.3)

# ---- 图 2: 离散多子目标 ----
e2 = copy.deepcopy(base_echo); e2["depth_span"]=0.0
e2["sub_targets"] = [{"dD": 0.75, "rho_rel": 0.6}, {"dD": 1.5, "rho_rel": 0.8}]  # 30.75m/31.5m
r_multi = signal_photon_rate_fine(e2, f_demo, tf_d)
e0 = copy.deepcopy(base_echo); e0["depth_span"]=0.0; e0["sub_targets"]=[]
r_single = signal_photon_rate_fine(e0, f_demo, tf_d)
ax[1].plot(tfd_ns, r_single, lw=1.4, color="tab:blue", label="单目标(30m)")
ax[1].plot(tfd_ns, r_multi, lw=1.8, color="tab:red", label="3子目标(30/30.75/31.5m)")
for dD, rr in [(0,1.0),(0.75,0.6),(1.5,0.8)]:
    ax[1].axvline(time_of_flight(30.0+dD)*1e9, color="gray", ls="--", alpha=0.5)
ax[1].set_xlabel("时间 t [ns]"); ax[1].set_ylabel("信号光子率 (未归一)")
ax[1].set_title("离散多子目标: 1.5m 深度内 3 物体合成复合回波")
ax[1].legend(fontsize=8); ax[1].grid(alpha=0.3)
plt.tight_layout(); plt.show()

# 定量输出
print("连续深度展宽 -> 回波上升沿(10-90%):")
for dz_cm in [0, 15, 75, 150]:
    e = copy.deepcopy(base_echo); e["depth_span"]=dz_cm/100.0; e["sub_targets"]=[]
    r = signal_photon_rate_fine(e, f_demo, tf_d); rise,_,_ = rise_1090(r, tf_d)
    print(f"  Δz={dz_cm:3d}cm: 时间全宽 2Δz/c={2*(dz_cm/100)/C_LIGHT*1e9:5.1f}ns -> 上升沿 {rise:.2f} ns")
print("结论: SPAD 架构下, 目标深度是 10ns 级上升沿的物理来源; Δz=1.5m 即达 ~10ns 时间全宽。")
print("(注: 全局 30m 场景参数未改, 本模块用局部 echo 演示; 电子学带宽按物理事实不建模)")''')

nb.cells.insert(i8b+1, md_8c)
nb.cells.insert(i8b+2, code_8c)

# ========== 5. 标题更新 ==========
if nb.cells and nb.cells[0].cell_type=="markdown":
    nb.cells[0].source = nb.cells[0].source.replace(
        "# 激光雷达直方图仿真 v14 (LiDAR Histogram Simulation) — 阵列俯视图 + Vov 曲线 + 参数化",
        "# 激光雷达直方图仿真 v15 (LiDAR Histogram Simulation) — 目标深度展宽（回波上升沿）")
    nb.cells[0].source += '''

**v15 相对 v14 的变化 — 目标深度展宽(回波上升沿的物理来源)**
- 新增**模块 6c**: 盘点时间展宽机制, 深度→时间映射, 深度展宽核函数 `depth_broadening_kernel()`。
- 扩展 `signal_photon_rate_fine`: 支持 (a) 连续深度分布 Δz 卷积; (b) 离散多子目标叠加。
  新字段缺省时(depth_span=0, sub_targets=[])行为与 v14 完全一致。
- 新增**模块 8c**: 演示 Δz=0/15/75/150cm 及多子目标下回波上升沿(用局部 echo, 不改全局 30m 场景)。
- **明确不建电子学带宽展宽**: SPAD+TDC 数字式探测只记录雪崩时刻, 后级带宽不影响时间戳,
  故 SPAD 架构下 10ns 上升沿的唯一物理来源是**目标深度**(Δz=1.5m ↔ 10ns)。'''

nbformat.write(nb, DST)
print(f"已生成 {DST}: 共 {len(nb.cells)} 个 cell。")
