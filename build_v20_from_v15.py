# -*- coding: utf-8 -*-
"""
基于 v15 生成 v20:
任务1: SPAD 计数阈值 Vth_frac 0.10 -> 0.60; 调 RC 时间常数使"可再次探测的最小间隔"(死区)=8ns。
       死区 t_deadzone = -tau_rc * ln(1 - Vth_frac); 令其=8ns 且 Vth_frac=0.60:
       tau_rc = 8ns / (-ln(1-0.60)) = 8ns / (-ln 0.40) = 8.7315 ns。
任务2: SPAD 对光子的响应 g(Vov) 成可选函数关系(linear / exp), 触发概率 = PDE_max * g(vov_frac);
       g(0)=0(Vov=0 完全不响应), g(1)=1(满 Vov 触发概率=PDE_max); exp 用凹形(快起饱和)。
       新增模块 7c: 画 g(Vov) 响应曲线(linear vs exp)。
任务3: 新增模块 13: 信噪比 SNR = S/sqrt(B) (信号峰 bin; S=扣背景信号计数, B=该 bin 背景计数)。

只改用户明确要求处; 30m 场景/脉冲/反射率/bin 宽等物理参数一律不动。
运行: python build_v20_from_v15.py
"""
import nbformat

SRC = "lidar_histogram_sim_v15.ipynb"
DST = "lidar_histogram_sim_v20.ipynb"
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

# ============================================================================
# 任务1 + 任务2a: 模块0 PARAMS["spad"] —— 改阈值/τ_RC, 新增响应函数字段
# ============================================================================
i0, c0 = find(lambda c: c.cell_type == "code" and '"tau_rc": 6e-9' in c.source)
old_spad = '''        "t_dead": 14e-9,        # 硬死时间 [s] (hard dead time 对照模型, non-paralyzable)
        "tau_rc": 6e-9,         # RC 恢复时间常数 τ=R·C_J [s] (RC 模型)
        "Vov_max": 3.3,         # 满过电压 Vov_max [V] (=V_bias − V_br)
        "Vth_frac": 0.10,       # 计数所需最小过电压占比 (Vth = Vth_frac·Vov_max)
        "reset_mode": "count",  # RC 复位方式: "count"=仅计数事件复位 / "all"=任何雪崩都复位
    },'''
new_spad = '''        "t_dead": 14e-9,        # 硬死时间 [s] (hard dead time 对照模型, non-paralyzable)
        # v20: Vth_frac 0.10->0.60, 并令 τ_RC 使"可再次探测的最小间隔"(死区)=8ns。
        #   死区 t_deadzone = -τ_RC·ln(1-Vth_frac); 令=8ns 且 Vth_frac=0.60 ->
        #   τ_RC = 8e-9 / (-ln(1-0.60)) = 8e-9 / 0.916291 = 8.7315 ns。
        "tau_rc": 8.7315e-9,    # RC 恢复时间常数 τ=R·C_J [s] (v20: 使死区=8ns, 见上)
        "Vov_max": 3.3,         # 满过电压 Vov_max [V] (=V_bias − V_br)
        "Vth_frac": 0.60,       # 计数所需最小过电压占比 (v20: 0.10->0.60; Vth=Vth_frac·Vov_max)
        "reset_mode": "count",  # RC 复位方式: "count"=仅计数事件复位 / "all"=任何雪崩都复位
        # --- v20 新增: SPAD 对光子的响应函数 g(vov_frac), 触发概率=PDE_max·g(vov_frac) ---
        #   约束: g(0)=0 (Vov=0 完全不响应), g(1)=1 (满 Vov 触发概率=PDE_max)。
        "resp_shape": "exp",    # "linear": g(x)=x (与 v15 一致) / "exp": 凹形快起饱和
        "resp_k": 3.0,          # exp 曲率(仅 resp_shape="exp" 生效): g(x)=(1-e^(-k·x))/(1-e^(-k)), k>0 凹
    },'''
assert old_spad in c0.source, "模块0 spad 锚点未匹配"
c0.source = c0.source.replace(old_spad, new_spad)

# ============================================================================
# 任务2b: 模块7b —— 定义响应函数 g; RC 引擎(两个函数)触发概率改为 PDE_max·g(vov_frac)
# ============================================================================
i7b, c7b = find(lambda c: c.cell_type == "code" and "def simulate_spad_shot_rc(" in c.source)

# (1) 在文件顶部插入响应函数 spad_response_g 的定义
old_head = '''def simulate_spad_shot_rc(r_sig_fine, r_amb_ph, tf, PDE_max, tau_rc, Vth_frac,
                          jitter_sigma, rng, reset_mode="count", Vov_max=3.3):'''
new_head = '''def spad_response_g(vov_frac, shape="linear", k=3.0):
    """SPAD 对光子的响应函数 g(vov_frac) ∈ [0,1]: 触发概率 = PDE_max · g(vov_frac)。
    约束: g(0)=0 (Vov=0 完全不响应, 对光子完全不触发), g(1)=1 (满 Vov 触发概率=PDE_max)。
      shape="linear": g(x)=x            (线性, PDE∝Vov; 与 v15 完全一致)
      shape="exp"   : g(x)=(1-e^(-k·x))/(1-e^(-k))  (k>0 凹形: 低 Vov 即有明显响应, 随后饱和)
    vov_frac 可为标量或数组。"""
    x = np.clip(vov_frac, 0.0, 1.0)
    if shape == "exp":
        return (1.0 - np.exp(-k * x)) / (1.0 - np.exp(-k))
    return x   # linear (默认)

def simulate_spad_shot_rc(r_sig_fine, r_amb_ph, tf, PDE_max, tau_rc, Vth_frac,
                          jitter_sigma, rng, reset_mode="count", Vov_max=3.3,
                          resp_shape="linear", resp_k=3.0):'''
assert old_head in c7b.source, "模块7b RC 函数头锚点未匹配"
c7b.source = c7b.source.replace(old_head, new_head)

# (2) simulate_spad_shot_rc 主循环触发判据: PDE_max*vov_frac -> PDE_max*g(vov_frac)
old_fire1 = '''        vov_frac = 1.0 - np.exp(-d) if d < 700 else 1.0   # 当前 Vov/Vov_max
        if u[k] < PDE_max * vov_frac:                     # 雪崩触发(概率∝Vov)'''
new_fire1 = '''        vov_frac = 1.0 - np.exp(-d) if d < 700 else 1.0   # 当前 Vov/Vov_max
        p_fire = PDE_max * spad_response_g(vov_frac, resp_shape, resp_k)  # 触发概率=PDE_max·g(Vov)
        if u[k] < p_fire:                                 # 雪崩触发(概率∝g(Vov))'''
assert old_fire1 in c7b.source, "模块7b 主引擎触发判据锚点未匹配"
c7b.source = c7b.source.replace(old_fire1, new_fire1)

# (3) trace 函数签名与触发判据同步加 g
old_trace_head = '''def simulate_spad_shot_rc_trace(r_sig_fine, r_amb_ph, tf, PDE_max, tau_rc, Vth_frac,
                                jitter_sigma, rng, reset_mode="count", Vov_max=3.3):'''
new_trace_head = '''def simulate_spad_shot_rc_trace(r_sig_fine, r_amb_ph, tf, PDE_max, tau_rc, Vth_frac,
                                jitter_sigma, rng, reset_mode="count", Vov_max=3.3,
                                resp_shape="linear", resp_k=3.0):'''
assert old_trace_head in c7b.source, "模块7b trace 函数头锚点未匹配"
c7b.source = c7b.source.replace(old_trace_head, new_trace_head)

old_fire2 = '''        vf = 1.0 - np.exp(-d) if d < 700 else 1.0
        vov_at[k] = vf
        if u[k] < PDE_max * vf:'''
new_fire2 = '''        vf = 1.0 - np.exp(-d) if d < 700 else 1.0
        vov_at[k] = vf
        if u[k] < PDE_max * spad_response_g(vf, resp_shape, resp_k):'''
assert old_fire2 in c7b.source, "模块7b trace 触发判据锚点未匹配"
c7b.source = c7b.source.replace(old_fire2, new_fire2)

# (4) 参数读取块: 增加 RESP_SHAPE / RESP_K, 打印同步更新(死区/响应)
old_read = '''VTH_FRAC = _sp["Vth_frac"]
RESET    = _sp["reset_mode"]
jit = JIT; t_dead = T_DEAD       # 小写别名(下游模块 8/9/11/12 沿用旧命名)

# 等效"硬死区"(恢复到 Vth 前完全不计数)与"渐变灵敏"分界
t_deadzone = -np.log(1 - VTH_FRAC) * TAU_RC
print(f"参数来自 PARAMS['spad']: PDE_max={PDE}, 硬死时间={T_DEAD*1e9:.1f}ns, τ_RC={TAU_RC*1e9:.1f}ns, "
      f"Vov_max={VOV_MAX}V, Vth={VTH_FRAC*100:.0f}%, reset='{RESET}'")
print(f"  Vov 恢复曲线: 1τ->{100*(1-np.exp(-1)):.0f}%, 2.3τ->{100*(1-np.exp(-2.3)):.0f}%, 5τ->{100*(1-np.exp(-5)):.1f}%")
print(f"  低于 Vth 的'硬死区'≈{t_deadzone*1e9:.2f} ns, 之后为渐变灵敏(与硬死时间一刀切不同)")'''
new_read = '''VTH_FRAC = _sp["Vth_frac"]
RESET    = _sp["reset_mode"]
RESP_SHAPE = _sp.get("resp_shape", "linear")   # v20: 响应函数形状
RESP_K     = _sp.get("resp_k", 3.0)            # v20: exp 曲率
jit = JIT; t_dead = T_DEAD       # 小写别名(下游模块 8/9/11/12 沿用旧命名)

# 计数死区: Vov 需恢复到 Vth 才能再次计数 -> 最小间隔 t_deadzone = -τ·ln(1-Vth_frac)
t_deadzone = -np.log(1 - VTH_FRAC) * TAU_RC
print(f"参数来自 PARAMS['spad']: PDE_max={PDE}, 硬死时间={T_DEAD*1e9:.1f}ns, tau_RC={TAU_RC*1e9:.3f}ns, "
      f"Vov_max={VOV_MAX}V, Vth={VTH_FRAC*100:.0f}%, reset='{RESET}'")
print(f"  响应函数 g(Vov): shape='{RESP_SHAPE}'" + (f", k={RESP_K}" if RESP_SHAPE=="exp" else "")
      + f"  (触发概率=PDE_max*g(vov_frac); g(0)=0, g(1)=1)")
print(f"  Vov 恢复曲线: 1tau->{100*(1-np.exp(-1)):.0f}%, 2.3tau->{100*(1-np.exp(-2.3)):.0f}%, 5tau->{100*(1-np.exp(-5)):.1f}%")
print(f"  [v20] 计数死区(可再次探测最小间隔) = -tau*ln(1-Vth) = {t_deadzone*1e9:.2f} ns  (目标 8ns)")'''
assert old_read in c7b.source, "模块7b 参数读取锚点未匹配"
c7b.source = c7b.source.replace(old_read, new_read)

# 模块7b markdown 补一句 v20 说明
i7bmd, c7bmd = find(lambda c: c.cell_type == "markdown" and c.source.lstrip().startswith("## 模块 7b"))
c7bmd.source += '''

> **v20 变化**: (1) 触发概率由 `PDE_max·vov_frac` 改为 **`PDE_max·g(vov_frac)`**, 响应函数
> `g` 可选 `linear`(=旧版) 或 `exp`(凹形快起饱和), 满足 g(0)=0、g(1)=1(见 `spad_response_g`);
> (2) 阈值 `Vth_frac` 提到 **0.60**, τ_RC 调到 **8.73ns** 使计数死区(可再次探测最小间隔)=**8ns**。'''

# ============================================================================
# 任务2c: 新增模块 7c —— 响应函数 g(Vov) 曲线 (插在模块7b code cell 之后)
# ============================================================================
md_7c = nbformat.v4.new_markdown_cell(r'''## 模块 7c（v20 新增）— SPAD 对光子的响应函数 g(Vov)

单个入射光子的**触发概率 = PDE_max · g(vov_frac)**, 其中 `vov_frac = Vov/Vov_max ∈ [0,1]`。
响应函数 `g` 描述"过电压越高、越容易触发雪崩", 满足两个物理边界:
- **g(0)=0**: 过电压为 0 时对光子**完全不响应**(刚雪崩淬灭、尚未恢复);
- **g(1)=1**: 过电压满值时触发概率 = **PDE_max**(=0.30, 峰值 PDE)。

两种可选形状(`PARAMS["spad"]["resp_shape"]`):
- **linear**: g(x)=x —— PDE 正比于 Vov(与 v13–v15 一致);
- **exp(凹形)**: g(x)=(1−e^(−k·x))/(1−e^(−k)), k>0 —— 低 Vov 即有明显响应、随后饱和, 更贴近真实 SPAD 的 PDE–Vov 曲线。

> 注意: 响应函数 g 决定"**能否触发雪崩**"; 计数阈值 `Vth_frac` 决定"**触发后幅度够不够被计数**"。二者独立。
> 本版全局采用 `shape='{}'`。'''.format("exp"))

code_7c = nbformat.v4.new_code_cell(r'''# 画响应函数 g(vov_frac): linear vs exp(几种 k), 并标注全局所用形状
xx = np.linspace(0, 1, 400)
fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))

# 左: g(x) 本身
ax[0].plot(xx, spad_response_g(xx, "linear"), lw=2.0, color="tab:blue", label="linear: g(x)=x")
for kk, c in zip([2, 3, 5], ["tab:orange", "tab:red", "tab:purple"]):
    ax[0].plot(xx, spad_response_g(xx, "exp", kk), lw=1.8, ls="--", color=c,
               label=f"exp k={kk}: (1-e^(-kx))/(1-e^(-k))")
ax[0].scatter([0, 1], [0, 1], c="k", zorder=5)
ax[0].annotate("g(0)=0\n(Vov=0 完全不响应)", (0, 0), xytext=(0.06, 0.16), fontsize=8,
               arrowprops=dict(arrowstyle="->", lw=0.8))
ax[0].annotate("g(1)=1\n(满 Vov -> 触发概率=PDE_max)", (1, 1), xytext=(0.40, 0.82), fontsize=8,
               arrowprops=dict(arrowstyle="->", lw=0.8))
ax[0].set_xlabel("vov_frac = Vov / Vov_max"); ax[0].set_ylabel("g(vov_frac)")
ax[0].set_title("响应函数 g(Vov): linear vs exp(凹形)")
ax[0].legend(fontsize=8, loc="upper left"); ax[0].grid(alpha=0.3)

# 右: 实际触发概率 = PDE_max·g, 叠加阈值线(Vth 左侧即便触发也不计数)
ax[1].plot(xx, PDE * spad_response_g(xx, "linear"), lw=1.6, color="tab:blue", alpha=0.7, label="linear")
ax[1].plot(xx, PDE * spad_response_g(xx, RESP_SHAPE, RESP_K), lw=2.2, color="tab:red",
           label=f"本版 shape='{RESP_SHAPE}'" + (f" k={RESP_K}" if RESP_SHAPE=="exp" else ""))
ax[1].axhline(PDE, color="gray", ls=":", lw=1.0, label=f"PDE_max={PDE}")
ax[1].axvline(VTH_FRAC, color="green", ls="--", lw=1.4, label=f"计数阈值 Vth_frac={VTH_FRAC:.2f}")
ax[1].axvspan(0, VTH_FRAC, color="green", alpha=0.06)
ax[1].set_xlabel("vov_frac = Vov / Vov_max"); ax[1].set_ylabel("单光子触发概率 = PDE_max·g")
ax[1].set_title("触发概率与计数阈值(阈值左侧: 触发也不计数)")
ax[1].legend(fontsize=8, loc="center right"); ax[1].grid(alpha=0.3)
plt.tight_layout(); plt.show()

print(f"全局响应函数: shape='{RESP_SHAPE}'" + (f", k={RESP_K}" if RESP_SHAPE=="exp" else ""))
print(f"  g(0)={spad_response_g(0.0, RESP_SHAPE, RESP_K):.3f}, g(0.5)={spad_response_g(0.5, RESP_SHAPE, RESP_K):.3f}, "
      f"g(1)={spad_response_g(1.0, RESP_SHAPE, RESP_K):.3f}  (触发概率 = {PDE}*g)")
print(f"  满 Vov 触发概率 = PDE_max*g(1) = {PDE*spad_response_g(1.0, RESP_SHAPE, RESP_K):.3f} (=PDE_max)")
if RESP_SHAPE == "exp":
    print(f"  exp 凹形: 在 vov_frac=Vth={VTH_FRAC:.2f} 处 g={spad_response_g(VTH_FRAC, RESP_SHAPE, RESP_K):.3f} "
          f"(线性同点 g={VTH_FRAC:.2f}); 低 Vov 响应被抬高。")''')

nb.cells.insert(i7b + 1, md_7c)
nb.cells.insert(i7b + 2, code_7c)

# ============================================================================
# 全局: 所有 simulate_spad_shot_rc(...) 调用补上 RESP_SHAPE, RESP_K 实参
#   原调用统一以 "..., RESET)" 结尾 -> 改为 "..., RESET, RESP_SHAPE, RESP_K)"
#   trace 调用以 "..., RESET, VOV_MAX)" 结尾 -> 改为 "..., RESET, VOV_MAX, RESP_SHAPE, RESP_K)"
# ============================================================================
for c in nb.cells:
    if c.cell_type != "code":
        continue
    if "simulate_spad_shot_rc" not in c.source:
        continue
    # trace 调用(先处理, 避免与下面的通用替换冲突)
    c.source = c.source.replace(
        "VTH_FRAC,\n                                     0.0, rng_t, RESET, VOV_MAX)",
        "VTH_FRAC,\n                                     0.0, rng_t, RESET, VOV_MAX, RESP_SHAPE, RESP_K)")
    # 普通 RC 调用: 结尾 ", jit, rng*, RESET)" 补两个响应参数
    import re as _re
    c.source = _re.sub(
        r"(simulate_spad_shot_rc\([^\n]*?, jit, (?:rng\w*|rng_rc|rng_ts|rng9|rng_a), RESET)\)",
        r"\1, RESP_SHAPE, RESP_K)", c.source)

# ============================================================================
# 任务3: 新增模块 13 —— 信噪比 SNR = S / sqrt(B)
#   在信号峰 bin: S = 扣背景后信号计数, B = 该 bin 背景(环境光)计数。
#   峰值宏像元 m_peak 用逐-SPAD RC 蒙卡累计 N_shots; 背景 B 由纯环境光另跑蒙卡估计。
# ============================================================================
# 找模块12 code cell(最后一个), 插到其后
i12, c12 = find(lambda c: c.cell_type == "code" and "峰值宏像元 m={m_peak} 单次 shot(RC)" in c.source)

md_13 = nbformat.v4.new_markdown_cell(r'''## 模块 13（v20 新增）— 信噪比 SNR = S / √B

在**峰值宏像元 m_peak** 的**信号峰 bin** 处计算信噪比(Signal-to-Noise Ratio, 信噪比):
$$\mathrm{SNR}=\frac{S}{\sqrt{B}}$$
- **S**(信号计数): 峰 bin 的总计数**扣除**该 bin 的背景 = `hist_peak − B`;
- **B**(背景计数): 仅环境光(+暗计数)在该 bin 的期望计数, 由**纯环境光蒙卡**(信号率置零)估计;
- 分母 **√B** 为背景散粒噪声(Poisson)标准差。

> 说明: 采用 **S/√B**(背景受限定义)。同一峰 bin 也给出 √(S+B) 供参考;
> 全部基于 **RC 模型 + 响应函数 g** 的逐-SPAD 蒙卡累计(N_shots), 与前面模块一致。''')

code_13 = nbformat.v4.new_code_cell(r'''# ---- 信噪比 SNR = S / sqrt(B), 在峰值宏像元 m_peak 的信号峰 bin ----
# 复用模块11 的 macro_hist(信号+背景, RC 逐-SPAD 蒙卡 N_shots)。
# 另跑"纯环境光"蒙卡(信号率=0)估计每 bin 背景计数 B。

# 1) 纯背景蒙卡: 峰值宏像元 27 SPAD, 信号率置 0, 仅环境光, 累计 N_shots
zero_rate = np.zeros_like(base_rate)
rng_bg = np.random.default_rng(PARAMS["hist"]["seed"] + 313)
bg_hist_peak = np.zeros(nbins)
for _shot in range(N_shots):
    ev_all = []
    for fij in macro_fvals[m_peak]:          # fij 不影响背景(背景与 f_pix 无关), 但保持 27 次一致
        ev = simulate_spad_shot_rc(zero_rate, r_amb_ph, tf, PDE, TAU_RC, VTH_FRAC, jit, rng_bg, RESET, RESP_SHAPE, RESP_K)
        if ev.size: ev_all.append(ev)
    if ev_all:
        bg_hist_peak += np.histogram(np.concatenate(ev_all), bins=edges)[0]

# 2) 信号峰 bin: 取峰值宏像元 macro_hist 的最大 bin 作为信号峰位置
sig_hist_peak = macro_hist[m_peak]
pk_bin = int(np.argmax(sig_hist_peak))
tot_peak = sig_hist_peak[pk_bin]             # 峰 bin 总计数(信号+背景)
B_peak = bg_hist_peak[pk_bin]                # 峰 bin 背景计数(纯环境光蒙卡)
S_peak = max(tot_peak - B_peak, 0.0)         # 扣背景后信号计数

# 背景平均(用离峰 bin 的纯背景蒙卡均值, 更稳健地代表底噪水平)
B_mean = bg_hist_peak.mean()
SNR_sqrtB   = S_peak / np.sqrt(B_peak) if B_peak > 0 else np.inf
SNR_sqrtB_m = S_peak / np.sqrt(B_mean) if B_mean > 0 else np.inf
SNR_sqrtSB  = S_peak / np.sqrt(S_peak + B_peak) if (S_peak + B_peak) > 0 else 0.0

print("="*76)
print(f"信噪比 SNR (峰值宏像元 m={m_peak}, 信号峰 bin @ {tc_ns[pk_bin]:.0f} ns, N_shots={N_shots})")
print(f"  峰 bin 总计数(信号+背景) = {tot_peak:.1f}")
print(f"  背景 B(纯环境光蒙卡, 该 bin) = {B_peak:.3f};  背景均值(全 bin) = {B_mean:.3f}")
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
          label="信号+背景 (RC 蒙卡)")
ax[0].bar(tc_ns, bg_hist_peak, width=bin_width*1e9, align="center", color="tab:red", alpha=0.6,
          label="纯背景 B (环境光蒙卡)")
ax[0].axvline(t0_ns, color="k", ls=":", alpha=0.7, label=f"真实 ToF {t0_ns:.1f} ns")
ax[0].annotate(f"峰 bin\nS={S_peak:.0f}, B={B_peak:.2f}\nSNR=S/sqrt(B)={SNR_sqrtB:.1f}",
               (tc_ns[pk_bin], tot_peak), xytext=(tc_ns[pk_bin]+4, tot_peak*0.8), fontsize=8,
               arrowprops=dict(arrowstyle="->", lw=0.8))
ax[0].set_xlabel("时间 t [ns]"); ax[0].set_ylabel(f"计数 / 1ns bin ({N_shots} shots)")
ax[0].set_title(f"峰值宏像元 m={m_peak}: 信号峰 S 与背景 B")
ax[0].legend(fontsize=8); ax[0].grid(alpha=0.3)

# 右: 各宏像元峰 bin 的 SNR
ax[1].plot(np.arange(n_macro), snr_per_macro, lw=1.6, marker="o", ms=3, color="tab:purple")
ax[1].axvline(m_peak, color="lime", ls=":", lw=1.2, label=f"峰值宏像元 m={m_peak} (SNR={snr_per_macro[m_peak]:.1f})")
ax[1].set_xlabel("宏像元序号 m (沿长边 y)"); ax[1].set_ylabel("SNR = S/sqrt(B) @ 峰 bin")
ax[1].set_title(f"各宏像元信号峰 bin 的 SNR (背景 B_mean={B_mean:.2f})")
ax[1].legend(fontsize=9); ax[1].grid(alpha=0.3)
plt.tight_layout(); plt.show()

print(f"峰值宏像元 SNR(S/sqrt(B)) = {SNR_sqrtB:.2f}; 边缘宏像元 m=0 SNR = {snr_per_macro[0]:.2f}")
print("(定义: SNR=S/sqrt(B), S=扣背景信号计数, B=峰 bin 背景计数; 背景由纯环境光 RC 蒙卡估计)")''')

nb.cells.insert(i12 + 1, md_13)
nb.cells.insert(i12 + 2, code_13)

# ============================================================================
# 标题更新
# ============================================================================
if nb.cells and nb.cells[0].cell_type == "markdown":
    nb.cells[0].source = nb.cells[0].source.replace(
        "# 激光雷达直方图仿真 v15 (LiDAR Histogram Simulation) — 目标深度展宽（回波上升沿）",
        "# 激光雷达直方图仿真 v20 (LiDAR Histogram Simulation) — 阈值/死区调整 + 响应函数 g(Vov) + 信噪比")
    nb.cells[0].source += '''

**v20 相对 v15 的变化**
- **任务1(阈值/死区)**: 计数阈值 `Vth_frac` 0.10 → **0.60**; τ_RC 由 6ns 调到 **8.73ns**,
  使"可再次探测的最小间隔"(计数死区 = −τ·ln(1−Vth_frac)) = **8ns**。
- **任务2(响应函数)**: SPAD 单光子触发概率由 `PDE_max·vov_frac` 改为 **`PDE_max·g(vov_frac)`**;
  `g` 可选 **linear**(=旧版) 或 **exp**(凹形快起饱和), 满足 g(0)=0(Vov=0 完全不响应)、g(1)=1(满 Vov→PDE_max)。
  新增**模块 7c** 画 g(Vov) 响应曲线; 本版全局用 `shape='exp'`(k=3)。
- **任务3(信噪比)**: 新增**模块 13** 计算 **SNR = S/√B**(信号峰 bin; S=扣背景信号计数, B=纯环境光蒙卡背景计数)。
- 其余 30m 场景/脉冲/反射率/bin 宽等物理参数**一律未改**。'''

nbformat.write(nb, DST)
print(f"已生成 {DST}: 共 {len(nb.cells)} 个 cell。")
