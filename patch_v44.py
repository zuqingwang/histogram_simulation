# -*- coding: utf-8 -*-
"""
patch_v44.py —— 从 v43 notebook 派生 v44（直接改 JSON，保留用户在 v43 里改过的参数）
用户三项要求（只改代码、不运行）：
  1) 模块12b 左图：正常整数采样点画不同 delta_dly 的波形，强调采样点固定在整数 ns、
     看不同 delta_dly 在整数采样点上的取值变化（加茎线更直观）。
  2) 模块19 峰值&面积 vs 距离：峰值用左 y 轴、面积用右 y 轴（twinx）。
  3) 模块20：改为 c/D² 拟合（c 为拟合参数）；已知足够远处必是平方反比，故只拟合"后 N 个点"，
     从小到大扫描 N（等价于起始距离从远towards近推进），找到误差首次超过阈值(如5%)的位置，
     判定"从几米开始符合平方反比"。
"""
import json

NB = "lidar_histogram_sim_v44.ipynb"
with open(NB, "r", encoding="utf-8") as f:
    nb = json.load(f)
cells = nb["cells"]


def find_cell(pred):
    for i, c in enumerate(cells):
        if pred(c):
            return i
    return -1


def set_src(i, text):
    cells[i]["source"] = text.splitlines(keepends=True)
    cells[i]["outputs"] = []           # 清掉旧输出（本次不运行，输出留着会误导）
    cells[i]["execution_count"] = None


# ============================================================================
# 改动 1：模块 12b（cell id=dc64ca18）——左图强调整数采样点 + 说明
#   只改左图 ax[0] 的绘制：原 plot(折线+o) → 折线 + 整数采样点茎线(vlines)+散点，
#   并在标题/注释里点明"采样点固定在整数 ns，delta_dly 只移动回波"。其余（右图、扫描、
#   拟合、打印、还原 timing）一律不动。
# ============================================================================
i12b = find_cell(lambda c: c.get("id") == "dc64ca18")
assert i12b >= 0, "未找到模块12b (dc64ca18)"
src12b = "".join(cells[i12b]["source"])

# 原左图循环块（精确匹配 v43 内容）
OLD_LEFT = '''# 左: 若干档的峰区平均直方图, 看回波随 delta_dly 平移
for dd in [0, 3, 6, 9, 11]:
    tm = dict(base_timing); tm["delta_dly"] = dd; tl = laser_delay(tm)
    h = _run_avg_hist(tl, PARAMS["hist"]["seed"] + 500 + dd*1000)
    ax[0].plot(tc_ns, h, marker="o", ms=3, lw=1.2, label=f"delta_dly={dd} (t_laser={tl*1e9:.2f}ns)")
ax[0].axvline(t0_ns, color="k", ls=":", alpha=0.6, label=f"ToF {t0_ns:.1f} ns")
ax[0].set_xlim(tc_ns[max(lo_w-1,0)], tc_ns[min(hi_w,nbins-1)])
ax[0].set_xlabel("时间 t [ns]"); ax[0].set_ylabel(f"平均二值计数 / bin ({N_shots} shots)")
ax[0].set_title(f"峰区平均直方图随 delta_dly 平移 (每档 1/12 ns, 平均 {N_REP_DLY} 次)")
ax[0].legend(fontsize=8); ax[0].grid(alpha=0.3)'''

NEW_LEFT = '''# 左: 若干档的峰区平均直方图, 看回波随 delta_dly 平移。
#   注意：采样点(HDC 时钟)始终【固定在整数 ns】(tc_ns 均为整数)，delta_dly 只让回波(雪崩
#   过阈窗)相对采样格平移 k/12 ns；因此这里就是"正常整数采样点采样"，观察不同 delta_dly
#   下【同一批整数采样点】上的取值如何变化（等价采样时钟相对回波偏移 2/12、3/12、6/12 ns…）。
_dd_show = [0, 3, 6, 9, 11]
_cmap_dd = plt.cm.viridis(np.linspace(0.1, 0.85, len(_dd_show)))
for _ci, dd in enumerate(_dd_show):
    tm = dict(base_timing); tm["delta_dly"] = dd; tl = laser_delay(tm)
    h = _run_avg_hist(tl, PARAMS["hist"]["seed"] + 500 + dd*1000)
    col = _cmap_dd[_ci]
    # 折线连接 + 整数采样点散点（强调采样落在整数 ns 格点上）
    ax[0].plot(tc_ns, h, "-", lw=1.0, alpha=0.55, color=col)
    ax[0].plot(tc_ns, h, "o", ms=5, color=col,
               label=f"delta_dly={dd} ({dd}/12 ns, t_laser={tl*1e9:.2f}ns)")
ax[0].axvline(t0_ns, color="k", ls=":", alpha=0.6, label=f"ToF {t0_ns:.1f} ns")
ax[0].set_xlim(tc_ns[max(lo_w-1,0)], tc_ns[min(hi_w,nbins-1)])
# 标出整数 ns 采样格点（竖直细线），直观说明采样时钟固定在整数 ns
for _xg in tc_ns[max(lo_w-1,0):min(hi_w,nbins-1)+1]:
    ax[0].axvline(_xg, color="0.85", lw=0.6, zorder=0)
ax[0].set_xlabel("时间 t [ns]（灰竖线=整数 ns 采样格点）"); ax[0].set_ylabel(f"平均二值计数 / bin ({N_shots} shots)")
ax[0].set_title(f"整数采样点上、不同 delta_dly 的波形 (每档 {1000/12:.1f} ps, 平均 {N_REP_DLY} 次)")
ax[0].legend(fontsize=8); ax[0].grid(alpha=0.3)'''

assert OLD_LEFT in src12b, "模块12b 左图未精确匹配（v43 内容可能已改）"
src12b_new = src12b.replace(OLD_LEFT, NEW_LEFT)
set_src(i12b, src12b_new)
print(f"[1/3] 模块12b(cell{i12b}) 左图已改：整数采样格点 + 不同 delta_dly 波形")


# ============================================================================
# 改动 2：模块 19（cell id=v43_module_19_energy → 改 id 为 v44_）——峰值&面积左右双 y 轴
# ============================================================================
i19 = find_cell(lambda c: c.get("id") == "v43_module_19_energy")
assert i19 >= 0, "未找到模块19"
src19 = "".join(cells[i19]["source"])

OLD_PA = '''fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(D_pa, peak_pa, "o-", color="tab:blue", lw=1.6, ms=5, label="峰值 peak")
ax.plot(D_pa, area_pa, "s-", color="tab:green", lw=1.6, ms=5, label="面积 area")
ax.set_xlabel("距离 D [m]"); ax.set_ylabel("计数")
ax.set_title(f"峰值 & 面积 随距离变化 (boost=1, N_shots={N_shots})")
ax.legend(fontsize=9); ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("v43_peak_area_vs_dist.png", dpi=110, bbox_inches="tight")
plt.show()'''

NEW_PA = '''#   用户要求：峰值用左 y 轴、面积用右 y 轴（twinx）。
fig, ax = plt.subplots(figsize=(8, 5))
ln1 = ax.plot(D_pa, peak_pa, "o-", color="tab:blue", lw=1.6, ms=5, label="峰值 peak（左轴）")
ax.set_xlabel("距离 D [m]")
ax.set_ylabel("峰值 peak [计数]", color="tab:blue"); ax.tick_params(axis="y", labelcolor="tab:blue")
ax.grid(alpha=0.3)
ax2 = ax.twinx()                                    # 右 y 轴画面积
ln2 = ax2.plot(D_pa, area_pa, "s-", color="tab:green", lw=1.6, ms=5, label="面积 area（右轴）")
ax2.set_ylabel("面积 area [计数]", color="tab:green"); ax2.tick_params(axis="y", labelcolor="tab:green")
ax.set_title(f"峰值(左轴) & 面积(右轴) 随距离变化 (boost=1, N_shots={N_shots})")
lns = ln1 + ln2; ax.legend(lns, [l.get_label() for l in lns], fontsize=9, loc="upper right")
plt.tight_layout()
plt.savefig("v44_peak_area_vs_dist.png", dpi=110, bbox_inches="tight")
plt.show()'''

assert OLD_PA in src19, "模块19 峰值面积图未精确匹配"
src19_new = src19.replace(OLD_PA, NEW_PA)
# 同步把模块19 其它 PNG 名 v43→v44
src19_new = src19_new.replace('savefig("v43_energy_scan.png"', 'savefig("v44_energy_scan.png"')
src19_new = src19_new.replace('savefig("v43_energy_waveforms.png"', 'savefig("v44_energy_waveforms.png"')
set_src(i19, src19_new)
cells[i19]["id"] = "v44_module_19_energy"
print(f"[2/3] 模块19(cell{i19}) 峰值&面积已改：左右双 y 轴；PNG→v44")


# ============================================================================
# 改动 3：模块 20（cell id=v43_module_20_snr → v44_）——c/D² 拟合 + 扫描起始距离
#   替换原"以最远档归一"的做法为：对"后 N 个点"用最小二乘拟合 c（rate ≈ c/D²·exp(-2αD)），
#   从小 N 到大 N 扫描（起点距离逐步向近推进），找误差首次 > TOL 的 N，判定平方反比起始距离。
# ============================================================================
i20 = find_cell(lambda c: c.get("id") == "v43_module_20_snr")
assert i20 >= 0, "未找到模块20"
src20 = "".join(cells[i20]["source"])

# 替换从"# ---- 平方反比拟合"到 cell 末尾的整段（拟合 + 绘图 + 打印）
_marker = "# ---- 平方反比拟合"
head20 = src20[:src20.index(_marker)]     # 保留数据采集部分（snr_at_distance + 扫描循环）

NEW_FIT_AND_PLOT = '''# ---- c/D² 拟合 + 平方反比起始距离扫描 ----
#   模型：rate(D) ≈ c / D² · exp(-2αD)，c 为拟合参数（含发射能量/口径/反射率等常数）。
#   已知足够远处必然是平方反比；近场因光斑/收集比例/饱和可能偏离。
#   做法：对【最远的 N 个点】做最小二乘拟合 c，再算这 N 个点的最大相对误差；
#         从小 N 到大 N 扫描（等价起始距离由远向近推进），找到误差首次 > TOL 的 N，
#         该 N 对应的起始距离即"从几米开始符合平方反比"。
alpha = PARAMS["channel"]["alpha"]
TOL_FIT = 0.05                      # 平方反比判定阈值（相对误差 5%，可调）
basis = (1.0 / D_list_B**2) * np.exp(-2*alpha*D_list_B)   # c 的基函数 1/D²·exp(-2αD)

def fit_c_lastN(y, N):
    """用最远 N 个点最小二乘拟合 c（模型 y≈c·basis），返回 (c, 该N点最大相对误差, 起始距离)。"""
    yl = y[-N:]; bl = basis[-N:]
    ok = np.isfinite(yl) & (bl > 0)
    if ok.sum() < 2:
        return np.nan, np.inf, D_list_B[-N]
    c = np.sum(bl[ok]*yl[ok]) / np.sum(bl[ok]**2)         # 最小二乘闭式解
    pred = c * bl
    rel = np.abs(yl - pred) / np.maximum(np.abs(pred), 1e-30)
    return c, np.nanmax(rel[ok]), D_list_B[-N]

def scan_start_distance(y):
    """从小 N 到大 N 扫描（起点从远向近推进），返回 (最优 N, 起始距离, c, 各N的最大误差数组)。
    判定：随 N 增大纳入更近的点，最大误差首次 > TOL_FIT 的前一个 N 即平方反比适用的最近边界。"""
    Ns = np.arange(3, len(D_list_B)+1)          # 至少 3 点起拟
    maxerr = np.full(Ns.size, np.nan); cs = np.full(Ns.size, np.nan)
    for j, N in enumerate(Ns):
        cs[j], maxerr[j], _ = fit_c_lastN(y, N)
    good = maxerr <= TOL_FIT
    if good.any():
        last_ok = np.where(good)[0].max()       # 仍满足阈值的最大 N
        N_best = Ns[last_ok]; D_start = D_list_B[-N_best]; c_best = cs[last_ok]
    else:
        N_best = Ns[0]; D_start = D_list_B[-N_best]; c_best = cs[0]
    return N_best, D_start, c_best, Ns, maxerr

# 分别对 信号率峰值 rate_B 与 SNR 做拟合与扫描
N_rate, Dstart_rate, c_rate, Ns_r, maxerr_r = scan_start_distance(rate_B)
finite = np.isfinite(snr_B)
snr_for_fit = np.where(finite, snr_B, np.nan)
N_snr, Dstart_snr, c_snr, Ns_s, maxerr_s = scan_start_distance(snr_for_fit)

# 拟合曲线（用各自最优 c）
rate_fit = c_rate * basis
snr_fit  = c_snr  * basis
err_rate = (rate_B - rate_fit) / np.maximum(np.abs(rate_fit), 1e-30) * 100
err_snr  = np.where(finite, (snr_B - snr_fit) / np.maximum(np.abs(snr_fit), 1e-30) * 100, np.nan)

# ---- 绘图：两图都做 c/D² 拟合；y轴按仿真范围；右侧y轴画相对误差 ----
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# ===== 左图：SNR vs 距离 =====
ax = axes[0]
ax.plot(D_list_B[finite], snr_B[finite], "o-", color="tab:purple", lw=1.6, ms=4, label="SNR（仿真）")
ax.plot(D_list_B, snr_fit, "--", color="tab:green", lw=1.5, alpha=0.85,
        label=f"c/D²·exp(-2αD) 拟合（后{N_snr}点, c={c_snr:.2e}）")
ax.axvline(Dstart_snr, color="tab:red", ls="-.", lw=1.3, alpha=0.8,
           label=f"平方反比起始 ≈ {Dstart_snr:.1f} m (误差<{TOL_FIT*100:.0f}%)")
ax.axhline(5.0, color="orange", ls=":", lw=1.0, alpha=0.7, label="SNR=5")
_sv = snr_B[finite]
if _sv.size: ax.set_ylim(0, _sv.max()*1.1)      # y 轴按仿真范围
ax.set_xlabel("距离 D [m]"); ax.set_ylabel("SNR = S / √B（仿真范围）")
ax.set_title(f"SNR vs 距离：c/D² 拟合 → 平方反比起始 ≈ {Dstart_snr:.1f} m")
ax.legend(fontsize=7.5, loc="upper right"); ax.grid(alpha=0.3)
axe = ax.twinx()
axe.plot(D_list_B, err_snr, "^:", color="tab:brown", lw=1.0, ms=3, alpha=0.6, label="相对误差 [%]")
axe.axhline(TOL_FIT*100, color="tab:red", ls=":", lw=0.8, alpha=0.6)
axe.axhline(-TOL_FIT*100, color="tab:red", ls=":", lw=0.8, alpha=0.6)
axe.set_ylabel("相对误差 [%]", color="tab:brown"); axe.tick_params(axis="y", labelcolor="tab:brown")
axe.legend(fontsize=8, loc="lower right")

# ===== 右图：信号光子率峰值 vs 距离 =====
ax = axes[1]
ax.plot(D_list_B, rate_B, "o-", color="tab:blue", lw=1.6, ms=4, label="信号光子率峰值（仿真, 未封顶）")
ax.plot(D_list_B, rate_fit, "--", color="tab:green", lw=1.5, alpha=0.85,
        label=f"c/D²·exp(-2αD) 拟合（后{N_rate}点, c={c_rate:.2e}）")
ax.axvline(Dstart_rate, color="tab:red", ls="-.", lw=1.3, alpha=0.8,
           label=f"平方反比起始 ≈ {Dstart_rate:.1f} m (误差<{TOL_FIT*100:.0f}%)")
ax.set_ylim(0, rate_B.max()*1.1)                # y 轴按仿真范围
ax.set_xlabel("距离 D [m]"); ax.set_ylabel("信号光子率峰值 [ph/s]（仿真范围）")
ax.set_title(f"信号率 vs 距离：c/D² 拟合 → 平方反比起始 ≈ {Dstart_rate:.1f} m")
ax.legend(fontsize=7.5, loc="upper right"); ax.grid(alpha=0.3)
axe = ax.twinx()
axe.plot(D_list_B, err_rate, "^:", color="tab:brown", lw=1.0, ms=3, alpha=0.6, label="相对误差 [%]")
axe.axhline(TOL_FIT*100, color="tab:red", ls=":", lw=0.8, alpha=0.6)
axe.axhline(-TOL_FIT*100, color="tab:red", ls=":", lw=0.8, alpha=0.6)
axe.set_ylabel("相对误差 [%]", color="tab:brown"); axe.tick_params(axis="y", labelcolor="tab:brown")
axe.legend(fontsize=8, loc="lower right")

plt.suptitle(f"模块 20 — SNR & 信号率 vs 距离：c/D² 拟合 + 平方反比起始距离扫描 "
             f"(ρ={echo0['rho']}, N_shots={N_shots}, {N_D_B}档, 阈值{TOL_FIT*100:.0f}%)", fontsize=11.5)
plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.savefig("v44_snr_vs_distance.png", dpi=110, bbox_inches="tight")
plt.show()

print("="*76)
print(f"模块 20 汇总：c/D²·exp(-2αD) 拟合（c 为拟合参数）+ 平方反比起始距离扫描（阈值 {TOL_FIT*100:.0f}%）")
print(f"  ① 信号率峰值：拟合 c={c_rate:.3e}；后 {N_rate} 点满足 → 平方反比【从 ≈ {Dstart_rate:.1f} m 起成立】")
print(f"  ② SNR       ：拟合 c={c_snr:.3e}；后 {N_snr} 点满足 → 平方反比【从 ≈ {Dstart_snr:.1f} m 起成立】")
print(f"  说明：从最远点起纳入越来越近的点做 c/D² 拟合，一旦某更近点使最大相对误差 >{TOL_FIT*100:.0f}%，")
print(f"        即认为平方反比在该距离以内不再成立。信号率(未封顶)通常比 SNR 更早满足平方反比。")
'''

src20_new = head20 + NEW_FIT_AND_PLOT
set_src(i20, src20_new)
cells[i20]["id"] = "v44_module_20_snr"
print(f"[3/3] 模块20(cell{i20}) 已改：c/D^2 拟合 + 扫描起始距离；PNG→v44")


# 顺带把 v44 新模块的 header id 同步（若存在）
ihdr = find_cell(lambda c: c.get("id") == "v43_new_header")
if ihdr >= 0:
    cells[ihdr]["id"] = "v44_new_header"

with open(NB, "w", encoding="utf-8") as f:
    json.dump(nb, f, ensure_ascii=False, indent=1)
print(f"\n已写出 {NB}（共 {len(cells)} cell）")
