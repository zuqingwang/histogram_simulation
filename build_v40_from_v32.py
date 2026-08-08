# -*- coding: utf-8 -*-
"""
build_v40_from_v32.py  （第二版：真正从 v32 派生）
====================================================
用户明确要求：忘掉上一版被旧 v40 带偏的做法，重新【从 v32 完整派生】v40。
经与用户确认的两项关键取舍：
  · notebook 结构 = 【完整 v32 + 把不必要模块注释停用】（不是精简重建）；
  · 蒙特卡罗二值引擎 = 【照搬 v32 逐光子 spad_binary_trace】（不做等价加速改写）。

因此本脚本的做法是：
  1) 逐字节读入 v32 的全部 56 个 cell，【原样保留】，一个字都不改；
  2) 把 8 个与 v40 主题无关的分析模块 code cell【转成 markdown】（代码放进 ```python 块，
     完整可见但不执行）——这就是"不必要的代码直接注释"，且经 AST 依赖分析确认它们是
     安全叶子（不被任何保留 cell 引用），停用不破坏执行链；
  3) 末尾【追加】3 个新 cell：
       · 新增说明（markdown）
       · 模块 A：能量扫描 → 前沿法/重心法定时 → dist-peak / dist-area 四条曲线；
       · 模块 B：SNR vs 距离。
     两个新模块复用 v32 的物理链路与【逐光子二值引擎】，全程蒙特卡罗。

只生成 lidar_histogram_sim_v40.ipynb，不改动 v32 及其它文件。

停用清单（8 个 code cell 及其模块，均与 v40 主题无关，经依赖分析为安全叶子）：
  cell23 模块7c  SPAD 响应函数 g(Vov) 单独对比图
  cell29 模块9   阵列内不同 SPAD 响应差异
  cell33 模块9c  单个 SPAD 二值时域波形演示
  cell41 模块17  不同反射率 ρ 的信号波形对比（被 v40 能量扫描取代）
  cell43 模块18  信号强度倍数扫描（正是 v40 能量扫描要做的，旧版重做）
  cell47 模块12b delta_dly 亚 ns 游标测距扫描
  cell53 模块15  多次蒙卡的 SNR 分布 + 正态拟合
  cell55 模块16  100ppm 噪点率理论阈值 + 海量蒙卡验证
"""
import json

SRC_NB = "lidar_histogram_sim_v32.ipynb"
OUT_NB = "lidar_histogram_sim_v40.ipynb"

# 要停用（转 markdown）的 code cell 索引，及一句停用理由
DISABLE = {
    23: "模块 7c：SPAD 响应函数 g(Vov) 单独对比图 —— v40 不需要，停用。",
    29: "模块 9：阵列内不同 SPAD 的响应差异 —— v40 不需要，停用。",
    33: "模块 9c：单个 SPAD 的二值时域波形演示 —— v40 不需要，停用。",
    41: "模块 17：不同反射率 ρ 的信号波形对比 —— 被 v40『能量扫描』取代，停用。",
    43: "模块 18：信号强度倍数扫描 —— 正是 v40『能量扫描』要做的，旧模块停用、下方新模块重做。",
    47: "模块 12b：delta_dly 亚 ns 游标测距扫描 —— v40 不需要，停用。",
    53: "模块 15：多次蒙卡的 SNR 分布 + 正态拟合 —— v40 用『SNR vs 距离』替代，停用。",
    55: "模块 16：100ppm 噪点率理论阈值 + 海量蒙卡验证 —— v40 不需要，停用。",
}

# 读入 v32
with open(SRC_NB, "r", encoding="utf-8") as f:
    nb = json.load(f)
cells = nb["cells"]
print(f"读入 {SRC_NB}: {len(cells)} 个 cell")

# ---- 把停用 code cell 转成 markdown（代码进 ```python 块，完整保留、不执行）----
n_disabled = 0
for i, reason in DISABLE.items():
    c = cells[i]
    assert c["cell_type"] == "code", f"cell{i} 不是 code，停用清单有误"
    code_src = "".join(c["source"])
    md = (f"> 🚫 **【v40 停用】** {reason}\n>\n"
          f"> 下面是该模块的原始代码（**已注释停用，不执行**），完整保留以备查阅：\n\n"
          f"```python\n{code_src}\n```\n")
    cells[i] = {
        "cell_type": "markdown",
        "id": c.get("id", f"disabled_{i}"),
        "metadata": {"v40_disabled": True},
        "source": md.splitlines(keepends=True),
    }
    n_disabled += 1
print(f"已停用（转 markdown）{n_disabled} 个 code cell: {sorted(DISABLE)}")

# ============================================================================
# 追加新 cell —— 新模块说明 + 模块 A（能量扫描）+ 模块 B（SNR vs 距离）
# ============================================================================

NEW_HEADER_MD = '''---
# 🆕 v40 新增分析（基于以上 v32 完整链路，全程蒙特卡罗）

以下两个模块是 v40 相对 v32 的新增内容，**直接复用上方 v32 已定义的物理链路与逐光子二值引擎**
（`spad_binary_trace` / `signal_photon_rate_fine` / `link_factor` / `pixel_collection_matrix` / `front_time_leading_edge` 等）：

- **模块 A — 能量扫描 → 前沿/重心定时 → dist-peak / dist-area 四条曲线**
  在 `cali_dist`（=`D0`=30 m）放单目标，扫描反射能量倍数 `boost ∈ [1e-5, 1e5]`（上下界/步长可调）。
  每档跑二值 MC 得直方图，用**前沿法**（v32 原生 `front_time_leading_edge`）与**重心法**（v40 新增，峰邻域质心）
  各定时一次；能量太低（峰 < 检测阈值 `det_th`）则**不定时、留空**。画四条曲线：
  dist-peak（前沿/重心各一条）+ dist-area（前沿/重心各一条）。

- **模块 B — SNR vs 距离**
  沿用 v32 的 SNR 定义（`SNR = S/√B`，B=峰 bin 纯背景计数），扫描距离 D，二值 MC 逐点估计。

> 缩写：COG（Center of Gravity，重心/质心）。其余缩写（TCSPC/SPAD/IRF/ToF/PDE/SNR/RC）见前文。
'''

# ---- 模块 A：能量扫描（复用 v32 逐光子引擎；重心法为 v40 新增）----
MODULE_A = r'''# ============================================================================
# 模块 A（v40 新增）— 能量扫描：前沿法/重心法定时 → dist-peak / dist-area 四条曲线
#
# 复用 v32：spad_binary_trace（逐光子二值引擎，原样）、signal_photon_rate_fine、
#          pixel_collection_matrix、front_time_leading_edge、laser_delay、护带机制。
# 能量注入方式：反射能量倍数 boost 乘到"单位收集比例信号率"上（ρ≤1 无法缩放出 1e5，
#              故用倍数表征反射能量，物理上等价于回波光子数的整体缩放）。
# ============================================================================

# ---- 可调参数 ----
BOOST_MIN = 1e-5          # 反射能量倍数下界（可调）
BOOST_MAX = 1e3           # 反射能量倍数上界（可调；验证阶段先用 1e3，确认无 bug 后可调回 1e5）
N_BOOST   = 12            # 能量档数（对数等分步长，可调；照搬逐光子引擎，高档较慢，勿设过大）
K_TH_A    = 5.0           # 检测阈值倍数 det_th = K_TH_A · nc_base（沿用 v32 模块14 默认 5.0）
COG_HALF  = 6            # 重心法窗口半宽 [bin]（峰 ±COG_HALF 内算质心）

boost_grid = np.logspace(np.log10(BOOST_MIN), np.log10(BOOST_MAX), N_BOOST)

# ---- 重心法（COG，v40 新增；v32 原生只有前沿法）----
def centroid_time_cog(hist, centers, pk_idx, half=COG_HALF):
    """重心法定时：取峰 pk_idx ±half 个 bin 的窗口算计数质心，返回 t_cog[s]。窗口内全 0 返回 nan。"""
    lo = max(0, pk_idx - half)
    hi = min(len(hist), pk_idx + half + 1)
    w = hist[lo:hi]; x = centers[lo:hi]
    return (x * w).sum() / w.sum() if w.sum() > 0 else np.nan

# ---- 单能量档：跑峰值宏像元二值 MC，返回该档直方图（完全复现 v32 模块11 的峰值宏像元流程）----
def run_peak_hist_boost(boost, p=PARAMS):
    """在 cali_dist(=D0) 处，对峰值宏像元 27 SPAD 跑 N_shots 的逐光子二值 MC，
    信号率整体 ×boost。返回该宏像元累加直方图（长度 nbins，与全局 centers 对齐）。
    —— 时间窗/护带/网格全部沿用 v32 上方已定义的 t_lo/t_hi/tf_gen/centers/nbins。"""
    rng = np.random.default_rng(p["hist"]["seed"] + 40010)   # v40 专用种子，不动 v32 其它 rng
    h = np.zeros(nbins)
    for _shot in range(N_shots):
        acc = np.zeros(nbins, dtype=np.int32)
        for fij in macro_fvals[m_peak]:                      # 峰值宏像元的 27 个 SPAD
            acc += spad_binary_trace(
                base_rate_gen * fij * boost, r_amb_ph, tf_gen, centers,
                PDE, TAU_RC, VTH_FRAC, jit, rng, T_OVER, T_LASER, RESP_SHAPE, RESP_K)
        h += acc
    return h

# ---- 背景基底（复用 v32 模块13 已算好的 bg_hist_peak / nc_base）----
det_th_A = K_TH_A * nc_base      # nc_base 来自模块13（纯背景全 bin 均值）
print(f"能量扫描：boost ∈ [{BOOST_MIN:.0e}, {BOOST_MAX:.0e}], {N_BOOST} 档（逐光子引擎，高能量档较慢）")
print(f"  检测阈值 det_th = {K_TH_A:.0f}×nc_base = {det_th_A:.3f}（峰 < 此值则不可分辨、留空）")

front_R = np.full(N_BOOST, np.nan)   # 前沿法测距 [m]
cog_R   = np.full(N_BOOST, np.nan)   # 重心法测距 [m]
peak_v  = np.zeros(N_BOOST)          # 峰 bin 计数
area_v  = np.zeros(N_BOOST)          # 直方图总面积（窗内总计数）

for k, boost in enumerate(boost_grid):
    h = run_peak_hist_boost(boost)
    pk = int(np.argmax(h))
    peak_v[k] = h[pk]
    area_v[k] = h.sum()
    detectable = h[pk] >= det_th_A          # 与 v32 前沿法同一可检测判据
    if detectable:
        V_dec_k = 0.5 * (det_th_A + h[pk])                      # 前沿判决电平(同 v32)
        tf_k, _, _ = front_time_leading_edge(h, centers, pk, V_dec_k, bin_width)
        if np.isfinite(tf_k):
            front_R[k] = C_LIGHT * tf_k / 2.0
        t_cog = centroid_time_cog(h, centers, pk)               # 重心法(v40 新增)
        if np.isfinite(t_cog):
            cog_R[k] = C_LIGHT * t_cog / 2.0
    tag = "可分辨" if detectable else "太低-留空"
    print(f"  [{k+1:>2d}/{N_BOOST}] boost={boost:.2e}  峰={h[pk]:>4.0f}  面积={h.sum():>5.0f}  {tag}")

mF = np.isfinite(front_R); mC = np.isfinite(cog_R)
print(f"有效定时：前沿 {mF.sum()}/{N_BOOST}，重心 {mC.sum()}/{N_BOOST}（低能量档按要求留空）")

# ---- 绘图：dist-peak / dist-area 四条曲线 ----
fig, axes = plt.subplots(1, 2, figsize=(14, 5.2))

# (左) dist-peak：横轴=峰计数，纵轴=测距，前沿/重心各一条
ax = axes[0]
ax.plot(peak_v[mF], front_R[mF], "o-", color="tab:blue", lw=1.6, ms=6,
        label=f"前沿法 (有效 {mF.sum()}/{N_BOOST})")
ax.plot(peak_v[mC], cog_R[mC], "s-", color="tab:red", lw=1.6, ms=6,
        label=f"重心法 (有效 {mC.sum()}/{N_BOOST})")
ax.axhline(D0, color="gray", ls=":", lw=1.2, alpha=0.8, label=f"真值 {D0} m")
ax.set_xscale("log"); ax.set_xlabel("直方图峰 bin 计数 peak")
ax.set_ylabel("测距结果 [m]"); ax.set_title("dist-peak（前沿法 & 重心法）")
ax.legend(fontsize=9); ax.grid(alpha=0.3, which="both")

# (右) dist-area：横轴=总面积，纵轴=测距，前沿/重心各一条
ax = axes[1]
ax.plot(area_v[mF], front_R[mF], "o-", color="tab:blue", lw=1.6, ms=6, label="前沿法")
ax.plot(area_v[mC], cog_R[mC], "s-", color="tab:red", lw=1.6, ms=6, label="重心法")
ax.axhline(D0, color="gray", ls=":", lw=1.2, alpha=0.8, label=f"真值 {D0} m")
ax.set_xscale("log"); ax.set_xlabel("直方图总面积 area（窗内总计数）")
ax.set_ylabel("测距结果 [m]"); ax.set_title("dist-area（前沿法 & 重心法）")
ax.legend(fontsize=9); ax.grid(alpha=0.3, which="both")

plt.suptitle(f"v40 模块 A — 能量扫描定时 (cali_dist={D0} m, N_shots={N_shots}, K_th={K_TH_A:.0f}, 逐光子二值 MC)",
             fontsize=12)
plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.savefig("v40_energy_scan.png", dpi=110, bbox_inches="tight")
plt.show()

print("="*76)
print(f"模块 A 汇总：cali_dist={D0} m, ρ={echo0['rho']}, boost∈[{BOOST_MIN:.0e},{BOOST_MAX:.0e}]（{N_BOOST}档）")
print(f"  四条曲线：dist-peak(前沿/重心) + dist-area(前沿/重心)")
if mF.any():
    print(f"  前沿法测距误差范围 [{(front_R[mF].min()-D0)*100:+.1f}, {(front_R[mF].max()-D0)*100:+.1f}] cm")
if mC.any():
    print(f"  重心法测距误差范围 [{(cog_R[mC].min()-D0)*100:+.1f}, {(cog_R[mC].max()-D0)*100:+.1f}] cm")
'''

# ---- 模块 B：SNR vs 距离（复用 v32 的 S/√B 定义与逐光子引擎）----
MODULE_B = r'''# ============================================================================
# 模块 B（v40 新增）— 不同距离下的信噪比 SNR 变化（逐光子二值 MC）
#
# 沿用 v32 模块13 的 SNR 定义：SNR = S/√B，其中
#   S = 峰 bin 总计数 − B（扣背景后信号），B = 峰 bin 处纯背景计数（纯环境光二值采样）。
# 对每个距离 D：重建时间窗/护带（同 v32 方式），跑"信号+背景"与"纯背景"各一次二值 MC。
# ============================================================================

# ---- 可调参数 ----
D_MIN_B = 5.0            # 最近距离 [m]
D_MAX_B = 250.0         # 最远距离 [m]
N_D_B   = 16            # 距离档数（逐光子引擎，勿设过大）
BOOST_B = 1.0           # 该扫描的反射能量倍数（默认 1 = 真实回波）

D_list_B = np.linspace(D_MIN_B, D_MAX_B, N_D_B)
print(f"SNR vs 距离：D ∈ [{D_MIN_B}, {D_MAX_B}] m, {N_D_B} 档 (boost={BOOST_B}, 逐光子二值 MC)")

def snr_at_distance(D, boost=BOOST_B, p=PARAMS):
    """在距离 D 处跑峰值宏像元二值 MC（信号+背景 与 纯背景各一次），返回 (SNR, 峰计数, B)。
    时间窗/护带/网格按 v32 方式随 D 重建；能量注入用 boost 倍数。"""
    # --- 随 D 重建时间窗与护带（复现 v32 模块8/11 的构造）---
    t0d = time_of_flight(D)
    t_lo_d, t_hi_d = t0d - pre, t0d + post
    tf_d = np.arange(t_lo_d, t_hi_d, dt_fine)
    guard = T_OVER + 5 * jit
    tf_gen_d = np.arange(t_lo_d - guard, t_hi_d, dt_fine)
    nb_d = len(tf_d)                                        # 采样点数(=nbins)
    edges_d = np.arange(t_lo_d, t_hi_d + bin_width/2, bin_width)
    centers_d = 0.5 * (edges_d[:-1] + edges_d[1:])
    nb_d = len(centers_d)
    # --- 该 D 的像元收集矩阵 & 峰值宏像元（随 D 变化，重算）---
    fpixD, _, _ = pixel_collection_matrix(D, p)
    fvalsD = [fpixD[:, m*By_m:(m+1)*By_m].ravel() for m in range(n_macro)]
    fsumD = np.array([fv.sum() for fv in fvalsD])
    mpk = int(fsumD.argmax())
    # --- 单位收集比例信号率（护带网格上），用当前 echo0 的 ρ/frac/tilt ---
    echoD = dict(echo0); echoD["D"] = D                     # 复用 echo0 的反射率等，仅改距离
    base_gen_D = signal_photon_rate_fine(echoD, 1.0, tf_gen_d, p)
    zero_gen_D = np.zeros_like(base_gen_D)
    # --- 信号+背景 ---
    rng_s = np.random.default_rng(p["hist"]["seed"] + 40100)
    h_sig = np.zeros(nb_d)
    for _shot in range(N_shots):
        acc = np.zeros(nb_d, dtype=np.int32)
        for fij in fvalsD[mpk]:
            acc += spad_binary_trace(base_gen_D * fij * boost, r_amb_ph, tf_gen_d, centers_d,
                                     PDE, TAU_RC, VTH_FRAC, jit, rng_s, T_OVER, T_LASER, RESP_SHAPE, RESP_K)
        h_sig += acc
    # --- 纯背景（信号率=0）---
    rng_b = np.random.default_rng(p["hist"]["seed"] + 40200)
    h_bg = np.zeros(nb_d)
    for _shot in range(N_shots):
        acc = np.zeros(nb_d, dtype=np.int32)
        for fij in fvalsD[mpk]:
            acc += spad_binary_trace(zero_gen_D, r_amb_ph, tf_gen_d, centers_d,
                                     PDE, TAU_RC, VTH_FRAC, jit, rng_b, T_OVER, T_LASER, RESP_SHAPE, RESP_K)
        h_bg += acc
    pk = int(np.argmax(h_sig))
    B = h_bg[pk]
    S = max(h_sig[pk] - B, 0.0)
    snr = S / np.sqrt(B) if B > 0 else (np.inf if S > 0 else 0.0)
    return snr, h_sig[pk], B

snr_B  = np.zeros(N_D_B)
peak_B = np.zeros(N_D_B)
bg_B   = np.zeros(N_D_B)
for k, D in enumerate(D_list_B):
    snr_B[k], peak_B[k], bg_B[k] = snr_at_distance(D)
    print(f"  D={D:>6.1f} m: 峰={peak_B[k]:>4.0f}  B={bg_B[k]:>5.2f}  SNR={snr_B[k]:>6.2f}")

# ---- 绘图 ----
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
finite = np.isfinite(snr_B)

ax = axes[0]   # SNR vs 距离
ax.plot(D_list_B[finite], snr_B[finite], "o-", color="tab:purple", lw=1.8, ms=6)
ax.axhline(5.0, color="orange", ls=":", lw=1.2, alpha=0.8, label="SNR=5（常用检测门限）")
ax.axhline(1.0, color="gray", ls=":", lw=1.0, alpha=0.7, label="SNR=1")
ax.set_xlabel("距离 D [m]"); ax.set_ylabel("SNR = S / √B")
ax.set_title(f"信噪比 vs 距离 (ρ={echo0['rho']}, N_shots={N_shots}, 逐光子二值 MC)")
ax.legend(fontsize=9); ax.grid(alpha=0.3)

ax = axes[1]   # 峰 & 背景 vs 距离
ax.semilogy(D_list_B, np.clip(peak_B, 0.1, None), "o-", color="tab:blue", lw=1.6, ms=5, label="峰 bin 计数")
ax.semilogy(D_list_B, np.clip(bg_B, 0.1, None), "s-", color="tab:red", lw=1.6, ms=5, label="峰 bin 背景 B")
ax.set_xlabel("距离 D [m]"); ax.set_ylabel("计数 / bin (log)")
ax.set_title(f"峰 bin 计数 & 背景 vs 距离 (N_shots={N_shots})")
ax.legend(fontsize=9); ax.grid(alpha=0.3, which="both")

plt.tight_layout()
plt.savefig("v40_snr_vs_distance.png", dpi=110, bbox_inches="tight")
plt.show()

print("="*76)
print(f"模块 B 汇总：SNR vs 距离 (ρ={echo0['rho']}, N_shots={N_shots}, boost={BOOST_B})")
if finite.any():
    print(f"  SNR 范围 [{snr_B[finite].min():.2f}, {snr_B[finite].max():.2f}]")
ge5 = np.where(snr_B >= 5.0)[0]
print(f"  SNR≥5 的最大距离 ≈ {D_list_B[ge5[-1]]:.1f} m" if len(ge5) else "  所有距离 SNR < 5")
'''


def code_cell(cid, src):
    return {"cell_type": "code", "id": cid, "metadata": {},
            "execution_count": None, "outputs": [],
            "source": src.splitlines(keepends=True)}

def md_cell(cid, src):
    return {"cell_type": "markdown", "id": cid, "metadata": {},
            "source": src.splitlines(keepends=True)}

cells.append(md_cell("v40_new_header", NEW_HEADER_MD))
cells.append(code_cell("v40_module_a_energy", MODULE_A))
cells.append(code_cell("v40_module_b_snr", MODULE_B))

nb["cells"] = cells
with open(OUT_NB, "w", encoding="utf-8") as f:
    json.dump(nb, f, ensure_ascii=False, indent=1)
print(f"已生成 {OUT_NB}: 共 {len(cells)} 个 cell "
      f"(v32 原 56 个，其中 {n_disabled} 个转 markdown 停用；末尾追加 3 个新 cell)")
