# -*- coding: utf-8 -*-
"""生成 `peak_vs_energy_v01.ipynb`（工作名 `peak_vs_energy`）。

改 cell 内容请改本文件再重跑，不要直接改 .ipynb —— 否则下次重建会被覆盖。

结构沿用 PoD_esti_v30 的三段式：
    模块说明（markdown，只留基本描述与缩写）
  → 画图参数 cell（每个参数注明它控制哪根轴）
  → 画图 cell（只读缓存，不重算）
"""
from __future__ import annotations

import ast
import json
import sys

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

NB_PATH = "peak_vs_energy_v01.ipynb"

CELLS: list[tuple[str, str]] = []


def md(src: str) -> None:
    CELLS.append(("markdown", src.strip("\n")))


def code(src: str) -> None:
    CELLS.append(("code", src.strip("\n")))


# =====================================================================
md(r"""
# 信号能量扫描：二值 SPAD 累加直方图的 peak 与半高全宽

**这个文件在做什么**：把回波信号能量从「平均不到一个光子」一路扫到「深度饱和」，
在**没有环境光（bg = 0）**的条件下，看宏像元 N 发累加直方图的
**峰值计数 peak** 和**半高全宽 FWHM** 怎么随能量变化——包括它们的分布、均值和标准差。

物理内核整个复用 `pod_esti_v30_core.py`，本文件不定义任何新物理；
蒙特卡洛由 `run_peak_energy_scan.py` 离线跑好落盘，这里只读缓存画图。

**像斑照度口径：`f_pix` 均匀分布。** 宏像元内 27 个 SPAD 的空间收集比例取成**完全相同**
（总量 $\sum f_{pix}$ 与真实像斑一致，只是把它平摊到 27 个 SPAD 上）。
于是同一个 bin 里各条轨迹的点亮概率 $p_t$ 严格相等，单 bin 计数是**标准二项分布**，
本文件的所有结论都可以直接和二项公式对照。
真实像斑（`f_pix` 最大/最小相差约 241 倍）的对照数据在
`peak_vs_energy_cache_real.npz`，切换方式见 `run_peak_energy_scan.py` 的 `F_PIX_MODE`。

**缩写与口径**

| 记号 | 含义 |
|---|---|
| SPAD | Single-Photon Avalanche Diode，单光子雪崩二极管，1 bit 器件，一个 bin 内至多记 1 |
| ToF | Time of Flight，飞行时间 |
| PDE | Photon Detection Efficiency，光子探测效率 |
| MC | Monte Carlo，蒙特卡洛 |
| FWHM | Full Width at Half Maximum，半高全宽 |
| `boost` | **信号能量倍率**，本文件的横轴。`boost = 1` 等于默认 ρ=0.1 场景的回波强度；只有比例有意义 |
| `bg` | `hist_add` 统计窗内每 bin 的平均计数。本文件固定 **bg = 0**（无环境光） |
| `hist_add` | N 发累加直方图，取值 0…`n_tr` |
| `n_tr` | 参与累加的轨迹数 = 27 SPAD × N_shots，即二值硬上限（N=1/2/4 → 27/54/108） |
| `peak` | `hist_add` 在统计窗 152 个 bin 内的最大计数 |
| `T_OVER` | 过阈窗宽，约 8 ns。一次雪崩会把连续约 8 个 bin 点亮 |

**FWHM 的口径**：对**单次实现的 `hist_add` 波形**测半高全宽（取包含峰位的那一段连续过半高区间，
两端线性插值到亚 bin，1 bin = 1 ns），再对多次蒙卡求平均。
`peak` 低于 4 计数的实现不计入（此时半高只有 1–2 个计数，FWHM 没有意义）。
""")

# =====================================================================
md(r"""
## 模块 0　载入缓存与派生统计量

读 `peak_vs_energy_cache_<F_PIX_MODE>.npz`（由 `run_peak_energy_scan.py` 生成），
从落盘的充分统计量派生 peak 均值/标准差、FWHM 均值/标准差、平均波形。
本 cell 不做任何蒙特卡洛。

缓存文件名带模式后缀，`uniform`（本文件默认，27 个 SPAD 照度均匀）与 `real`（真实像斑轮廓）
各存各的，互不覆盖。载入后会打印当前是哪一种，以及 `f_pix` 的实际离散度。

缓存里存的是分布本身（`peak_cnt` 是 peak 的完整计数分布，`fwhm_cnt` 是 FWHM 的分布），
不是原始样本，所以文件很小但信息完整。
""")

code(r"""
import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

# 中文字体：必须给成【列表】才有回退链。只写一个字体名时，
# 该字体缺哪个字形就直接画成方框（项目里旧脚本的老毛病）。
plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

os.environ.setdefault("POD_CORE_QUIET", "1")
import pod_esti_v30_core as core
import run_peak_energy_scan as scan

CACHE = scan.CACHE
if not os.path.exists(CACHE):
    raise FileNotFoundError(
        f"找不到 {CACHE}。先在终端跑：python -u run_peak_energy_scan.py --workers 20 --n-mc 20000")

Z       = np.load(CACHE, allow_pickle=False)
BOOSTS  = Z["boosts"]
N_MC    = int(Z["n_mc"])
BG      = float(Z["bg"])
FPIX    = Z["f_pix"]
FPIX_MODE = str(Z["f_pix_mode"])
FPIX_CV = float(FPIX.std() / FPIX.mean())      # 离散度：0 = 照度完全均匀
N_LIST  = list(scan.N_SHOTS_LIST)
TC_NS   = core.TC_NS
FW_EDGE = float(Z["fwhm_bin_ns"]) * np.arange(Z[f"fwhm_cnt_{N_LIST[0]}"].shape[1] + 1)
FW_CTR  = 0.5 * (FW_EDGE[1:] + FW_EDGE[:-1])

_COLOR_N = {1: "#1f77b4", 2: "#d62728", 4: "#2ca02c"}

LIN_FIT_PEAK_MAX = 2.0    # 拟合「低能线性参考斜率」时只用 peak 均值低于此值的档
DEEP_BOOST_MIN   = 20.0   # 拟合「FWHM 每十倍能量增宽」时只用 boost 高于此值的档（已完全封顶）


def autoscale_y(a, series, pad=0.08, floor_zero=True):
    # 按【当前 xlim 之内】的数据定纵轴。matplotlib 的自动缩放看的是全部数据，
    # 只设 xlim 不设 ylim 时，远处的数据（比如深饱和段的 peak=106）会把放大图的纵轴撑开，
    # 低能段被压成贴地的一条平线。series 是若干 (x, y) 对。
    lo, hi = a.get_xlim()
    vals = [y[(x >= lo) & (x <= hi) & np.isfinite(y)] for x, y in series]
    vals = [v for v in vals if v.size]
    if not vals:
        return
    v = np.concatenate(vals)
    y0, y1 = float(v.min()), float(v.max())
    d = max(y1 - y0, 1e-9) * pad
    a.set_ylim(0.0 if (floor_zero and y0 >= 0) else y0 - d, y1 + d)


def _mu_sd_from_counts(cnt, values):
    tot = np.maximum(cnt.sum(axis=1).astype(float), 1.0)
    mu = (cnt * values).sum(axis=1) / tot
    m2 = (cnt * values * values).sum(axis=1) / tot
    return mu, np.sqrt(np.maximum(m2 - mu * mu, 0.0))


STAT = {}
for n in N_LIST:
    cnt = Z[f"peak_cnt_{n}"]
    n_tr = core.N_PIX_MACRO * n
    pmu, psd = _mu_sd_from_counts(cnt, np.arange(cnt.shape[1])[None, :])

    nval = Z[f"fwhm_nval_{n}"].astype(float)
    with np.errstate(invalid="ignore", divide="ignore"):
        fmu = Z[f"fwhm_sum_{n}"] / nval
        fsd = np.sqrt(np.maximum(Z[f"fwhm_sumsq_{n}"] / nval - fmu * fmu, 0.0))
    ok = nval > 0

    S = dict(
        n_tr=n_tr,
        peak_pmf=cnt / float(N_MC),
        peak_mu=pmu, peak_sd=psd,
        fwhm_mu=np.where(ok, fmu, np.nan),
        fwhm_sd=np.where(ok, fsd, np.nan),
        fwhm_valid=nval / float(N_MC),
        fwhm_pmf=Z[f"fwhm_cnt_{n}"] / np.maximum(nval, 1.0)[:, None],
        wave=Z[f"wave_sum_{n}"] / float(N_MC),
    )

    # 低能线性参考：过原点最小二乘。peak 很小时 E[peak] ≈ 期望总计数，本该严格正比于能量
    m = (BOOSTS > 0) & (pmu <= LIN_FIT_PEAK_MAX)
    S["lin_slope"] = float((pmu[m] * BOOSTS[m]).sum() / (BOOSTS[m] ** 2).sum()) if m.sum() >= 2 else np.nan

    # 深饱和段 FWHM 的对数增长率：FWHM = a + b·log10(boost)
    d = (BOOSTS >= DEEP_BOOST_MIN) & np.isfinite(S["fwhm_mu"])
    if d.sum() >= 3:
        b_, a_ = np.polyfit(np.log10(BOOSTS[d]), S["fwhm_mu"][d], 1)
        S["fwhm_per_decade"], S["fwhm_a"] = float(b_), float(a_)
    else:
        S["fwhm_per_decade"], S["fwhm_a"] = np.nan, np.nan
    STAT[n] = S

print("=" * 96)
print(f"载入 {CACHE}：{BOOSTS.size} 个能量档，每档 {N_MC:,} 次 MC，bg = {BG:g}（无环境光）")
print(f"  boost 范围 {BOOSTS[BOOSTS > 0].min():.4g} → {BOOSTS.max():.4g}，非均匀采样（坐标轴仍用线性）")
print(f"  宏像元 {core.N_PIX_MACRO} SPAD，统计窗 bin [{core.I_STAT0}, {core.I_STAT1})，共 {core.N_STAT} 个 bin")
print(f"  像斑照度 f_pix = '{FPIX_MODE}'：Σf_pix = {FPIX.sum():.5f}，"
      f"每 SPAD {FPIX.min():.3e} … {FPIX.max():.3e}，离散度 std/mean = {FPIX_CV:.3f}"
      + ("（完全均匀 ⇒ 各轨迹 p_t 相等 ⇒ 单 bin 服从二项分布）" if FPIX_CV < 1e-9 else
         "（不均匀 ⇒ 各轨迹 p_t 不同 ⇒ 单 bin 服从泊松二项分布）"))
print(f"  FWHM 有效门槛：peak >= {int(Z['fwhm_min_peak'])} 计数")
print("=" * 96)
for n in N_LIST:
    S = STAT[n]
    i_sat = int(np.argmax(S["peak_mu"] >= 0.95 * S["peak_mu"].max()))
    print(f"  N_shots={n}（n_tr={S['n_tr']:>3}）: peak 均值封顶于 {S['peak_mu'].max():6.2f} "
          f"（{100 * S['peak_mu'].max() / S['n_tr']:.1f}% 上限，boost≈{BOOSTS[i_sat]:.3g} 起）；"
          f"peak 标准差最大 {np.nanmax(S['peak_sd']):.2f}；"
          f"FWHM 深饱和段每十倍能量增宽 {S['fwhm_per_decade']:.2f} ns")
""")

# =====================================================================
md(r"""
## 模块 1　能量扫描总览：平均波形怎么从单光子长成饱和平台

沿用 `lidar_histogram_sim_v45.ipynb` 模块 18 的画法：一条曲线一个能量档，颜色按 `boost` 对数映射，
用色条代替十几条图例。

- **左图**：平均 `hist_add` 波形（线性纵轴）。看削顶饱和平台是怎么先在峰区出现、再向两侧展宽的。
- **右图**：同样的波形**各自归一化到自身峰值**。峰高被抹掉后，只剩形状变化——
  这就是 FWHM 增长的直观来源。灰色横线是半高 0.5，它与曲线的两个交点之间就是 FWHM。
""")

code(r"""
# ---------------- 模块 1 画图参数 ----------------
M1_N_SHOW    = 4            # 【画哪一档】只画一个 N_shots，默认 N=4（n_tr=108）
M1_B_MIN     = 1e-3         # 【选哪些能量档】下界：再弱的波形几乎全零，画上去只是贴地一条线
M1_B_MAX     = 1e4          # 【选哪些能量档】上界
M1_N_CURVES  = 18           # 【画几条曲线】在上面区间里按 log10(boost) 等距挑这么多档
M1_XLIM_NS   = (94, 146)    # 【两图横轴】时间 [ns]。ToF = 100 ns，最强档的过半高区间要到约 131 ns，
                            #   右端留够余量才不会把最宽那几条的下降沿切掉
M1_YLIM_L    = None         # 【左图纵轴】计数范围；None = 自动
M1_SHOW_CAP  = True         # 【左图】是否画二值硬上限 n_tr 的水平线
M1_FIGSIZE   = (14.6, 5.3)
""")

code(r"""
_n = M1_N_SHOW
_S = STAT[_n]

_sel = np.where((BOOSTS >= M1_B_MIN) & (BOOSTS <= M1_B_MAX))[0]
_lb = np.log10(BOOSTS[_sel])
_pick = _sel[np.unique([int(np.argmin(np.abs(_lb - t)))
                        for t in np.linspace(_lb.min(), _lb.max(), M1_N_CURVES)])]

_norm = matplotlib.colors.LogNorm(vmin=BOOSTS[_pick].min(), vmax=BOOSTS[_pick].max())
_cmap = plt.cm.viridis

fig, ax = plt.subplots(1, 2, figsize=M1_FIGSIZE, constrained_layout=True)
for i in _pick:
    c = _cmap(_norm(BOOSTS[i]))
    w = _S["wave"][i]
    ax[0].plot(TC_NS, w, color=c, lw=1.35)
    if w.max() > 0:
        ax[1].plot(TC_NS, w / w.max(), color=c, lw=1.35)

if M1_SHOW_CAP:
    ax[0].axhline(_S["n_tr"], color="k", ls="-.", lw=1.2, alpha=0.85,
                  label=f"二值硬上限 n_tr = {_S['n_tr']}")
    ax[0].legend(fontsize=8.5, loc="upper right")
ax[0].set_ylabel(f"平均计数 / 1 ns bin（N_shots={_n} 累加，{N_MC:,} 次 MC 平均）")
ax[0].set_title("左：平均累加波形——弱信号贴地，强信号削顶成平台", fontsize=10.5)
if M1_YLIM_L is not None:
    ax[0].set_ylim(*M1_YLIM_L)

ax[1].axhline(0.5, color="0.35", ls="--", lw=1.4, label="半高 0.5（FWHM 由此量取）")
ax[1].set_ylabel("归一化到各自峰值")
ax[1].set_ylim(0, 1.06)
ax[1].set_title("右：各自归一化——峰高抹掉后只剩形状，能量越大越宽", fontsize=10.5)
ax[1].legend(fontsize=8.5, loc="upper right")

for a in ax:
    a.set_xlabel("时间 t [ns]")
    a.set_xlim(*M1_XLIM_NS)
    a.grid(alpha=0.3)

_sm = plt.cm.ScalarMappable(cmap=_cmap, norm=_norm); _sm.set_array([])
_cb = fig.colorbar(_sm, ax=ax, fraction=0.035, pad=0.015)
_cb.set_label("信号能量倍率 boost（对数刻度）")

fig.suptitle(f"模块 1　能量扫描总览（bg = {BG:g}，N_shots = {_n}，"
             f"boost {BOOSTS[_pick].min():.3g} → {BOOSTS[_pick].max():.3g}）", fontsize=12.5)
plt.show()
""")

# =====================================================================
md(r"""
## 模块 2　peak 的分布怎么随能量变化

热图的每一竖列是一个能量档上 `peak` 的完整概率分布（由 20,000 次 MC 的计数分布直接给出，不是拟合）。
白线是均值，虚线是均值 ± 1 标准差。

上排是低能段放大，下排是全局。三列分别是 N_shots = 1 / 2 / 4，注意三者的二值硬上限不同（27 / 54 / 108）。
""")

code(r"""
# ---------------- 模块 2 画图参数 ----------------
M2_XLIM_ZOOM   = (0.0, 0.005)   # 【上排横轴】低能段 boost 范围
M2_XLIM_GLOBAL = (0.0, 0.5)     # 【下排横轴】全局 boost 范围。照度均匀时全部 SPAD 同步饱和，
                                #   boost≈0.3 就封顶，横轴铺到 0.5 已经走完整个压缩过程
M2_YFRAC       = 1.06           # 【纵轴】画到 n_tr 的多少倍
M2_YFRAC_ZOOM  = None           # 【上排纵轴】单独指定 peak 上界；None = 自动贴合该段数据
M2_CMAP        = "magma"        # 【配色】概率密度
M2_NORM_COLUMN = True           # 【配色口径】True = 每一竖列各自除以本列最大值。
                                #   低能档的概率几乎全压在 peak=0 上（≈1），不归一化会把色标拉满、
                                #   让其他列全黑。定量数值请看叠加的均值线与模块 6 的表
M2_VMAX_Q      = 0.99           # 【颜色上限】取（归一化后）概率的这个分位，防止个别极窄分布吃掉色标
M2_SHOW_BAND   = True           # 【是否叠】均值线与 ±1 标准差虚线
M2_FIGSIZE     = (16.8, 8.6)
""")

code(r"""
def _edges(x):
    # 由非均匀采样点生成 pcolormesh 需要的格子边界
    m = 0.5 * (x[1:] + x[:-1])
    return np.concatenate([[x[0] - (m[0] - x[0])], m, [x[-1] + (x[-1] - m[-1])]])


_xe = _edges(BOOSTS)
fig, axes = plt.subplots(2, len(N_LIST), figsize=M2_FIGSIZE)

for row, xlim in enumerate([M2_XLIM_ZOOM, M2_XLIM_GLOBAL]):
    _in = (BOOSTS >= xlim[0]) & (BOOSTS <= xlim[1])
    for col, n in enumerate(N_LIST):
        a = axes[row, col]
        S = STAT[n]
        pmf = S["peak_pmf"]
        ye = np.arange(pmf.shape[1] + 1) - 0.5

        img = pmf / np.maximum(pmf.max(axis=1, keepdims=True), 1e-12) if M2_NORM_COLUMN else pmf
        _v = img[_in]
        vmax = np.quantile(_v[_v > 0], M2_VMAX_Q) if (_v > 0).any() else 1.0
        im = a.pcolormesh(_xe, ye, img.T, cmap=M2_CMAP, vmin=0, vmax=vmax,
                          shading="flat", rasterized=True)

        if M2_SHOW_BAND:
            a.plot(BOOSTS, S["peak_mu"], "-", color="w", lw=1.7, label="peak 均值")
            a.plot(BOOSTS, S["peak_mu"] + S["peak_sd"], "--", color="w", lw=1.0, alpha=0.85,
                   label="均值 ± 1 标准差")
            a.plot(BOOSTS, np.maximum(S["peak_mu"] - S["peak_sd"], 0), "--", color="w",
                   lw=1.0, alpha=0.85)

        a.set_xlim(*xlim)
        if row == 0:
            top = (M2_YFRAC_ZOOM if M2_YFRAC_ZOOM is not None
                   else max(2.0, 1.35 * (S["peak_mu"][_in] + 2 * S["peak_sd"][_in]).max()))
            a.set_ylim(-0.5, top)
        else:
            a.set_ylim(-0.5, S["n_tr"] * M2_YFRAC)
            a.axhline(S["n_tr"], color="c", ls="-.", lw=1.2, alpha=0.9)
        a.set_xlabel("信号能量倍率 boost")
        if col == 0:
            a.set_ylabel(("低能段放大\n" if row == 0 else "全局\n") + "peak [计数/bin]")
        a.set_title(f"N_shots={n}（n_tr={S['n_tr']}）", fontsize=10.5)
        if row == 0 and col == 0:
            # 白字必须配深底：这一格左上角常常是亮黄色，透明图例会看不清
            a.legend(fontsize=7.8, loc="upper left", labelcolor="w",
                     facecolor="0.15", edgecolor="0.5", framealpha=0.8)
        fig.colorbar(im, ax=a, fraction=0.045, pad=0.02).set_label(
            "列内归一化概率" if M2_NORM_COLUMN else "概率", fontsize=8)

fig.suptitle(f"模块 2　peak 的概率分布随信号能量的变化（bg = {BG:g}，每档 {N_MC:,} 次 MC）",
             fontsize=12.5)
fig.tight_layout(rect=[0, 0, 1, 0.95])
plt.show()
""")

# =====================================================================
md(r"""
## 模块 3　图 i：peak 均值 vs 信号能量

三个横轴范围看同一条曲线（全部线性坐标轴）：

- **左：低能段** —— 检验 peak 与能量是否成正比。虚线是过原点的线性外推参考，
  斜率由 peak 均值低于 2 计数的那些档拟合得到。曲线什么时候离开这条虚线，就是饱和开始的地方。
- **中：全局** —— 一路到约 95% 二值硬上限，看完整的压缩过程。
- **右：穿过饱和** —— peak 彻底封顶后就是一条平线（真正还在变的是 FWHM，见模块 5）。

**另附三张体检图**：

1. **线性拟合双轴图**（`boost ∈ [0, 0.1]`）：左 = peak + 过原点直线；右 = 相对误差。
2. **饱和指数拟合双轴图**（同一 boost 窗）：形式
   $\mathrm{peak} \approx A\,(1-e^{-\alpha\cdot\mathrm{boost}})$
   （不是无界的 $A e^{bx}$——peak 有硬上限 $n_{tr}$）。
   低能极限 $A\alpha\cdot\mathrm{boost}$ 应接近线性斜率。左 = peak + 拟合曲线；右 = 相对误差。
3. **5% 线性区图**：相对低能线性斜率的误差，标出连续满足 `|err|≤5%` 的 peak 上限。
""")

code(r"""
# ---------------- 模块 3 画图参数 ----------------
M3_XLIM_ZOOM      = (0.0, 0.005)  # 【左图横轴】低能段 boost 范围
M3_XLIM_GLOBAL    = (0.0, 0.5)    # 【中图横轴】全局 boost 范围（照度均匀时 boost≈0.3 即封顶）
M3_XLIM_DEEP      = (0.0, 200.0)  # 【右图横轴】穿过饱和的 boost 范围
M3_YLIM_ZOOM      = None          # 【左图纵轴】peak 计数；None = 自动
M3_YLIM_GLOBAL    = None          # 【中图纵轴】None = 自动（0 → 最大 n_tr）
M3_YLIM_DEEP      = None          # 【右图纵轴】None = 自动
M3_SHOW_LINREF    = True          # 【左图】是否叠过原点的线性外推参考虚线
M3_SHOW_CAP       = True          # 【中/右图】是否画各档的二值硬上限 n_tr
M3_USE_LOG_X_DEEP = False         # 【右图横轴】True = 改成对数轴。线性轴上 4 个数量级会压扁，
                                  #   想看完整深饱和段可以打开；定量结论另见模块 6 的表
M3_MARKER_MS      = 3.4
M3_FIGSIZE        = (16.8, 5.2)
# ---- 拟合双轴图（boost 窗；左 peak / 右相对误差）----
M3_FIT_XLIM       = (0.0, 0.1)    # 【拟合图横轴】信号能量倍率 boost 范围（线性/指数共用）
M3_FIT_YLIM_PEAK  = None          # 【左纵轴】peak 均值；None = 按窗内数据自动
M3_FIT_YLIM_ERR   = None          # 【右纵轴】相对误差；None = 按窗内 err 自动（含正负）
M3_FIT_SLOPE_MODE = "low_energy"  # 【线性斜率口径】"low_energy" / "full_window"
M3_FIT_SHOW_TOL   = False         # 【线性拟合右纵轴】是否叠 ±M3_LIN_TOL 灰带（默认关掉）
M3_FIT_FIGSIZE    = (15.6, 4.8)
# ---- 饱和指数拟合 ----
# 形式 peak ≈ A * (1 - exp(-α * boost))；低能极限斜率 = A*α
M3_EXP_ENABLE     = True          # 【开关】是否画饱和指数拟合双轴图
M3_EXP_A_MODE     = "free"        # 【渐近线】"free" = 拟合 A；"n_tr" = 固定 A = n_tr（二值硬上限）
M3_EXP_SHOW_LIN   = True          # 【左纵轴】是否叠低能线性参考（便于对比）
M3_EXP_SHOW_TOL   = False         # 【指数拟合右纵轴】是否叠 ±M3_LIN_TOL 灰带（默认关掉）
M3_EXP_YLIM_ERR   = None          # 【右纵轴】相对误差范围；None = 自动
# ---- 5% 线性区附图 ----
M3_LIN_TOL        = 0.05          # 【线性容差】|peak/(斜率×boost) − 1| ≤ 此值 视为仍线性
M3_LIN_XLIM_PEAK  = None          # 【附图横轴】peak 均值范围；None = 自动到略超过 5% 越界点
M3_LIN_YLIM_ERR   = (-0.20, 0.05) # 【附图纵轴】相对误差；饱和后误差朝负向走（实测 < 线性外推）
M3_LIN_FIGSIZE    = (8.4, 4.8)
""")

code(r"""
fig, ax = plt.subplots(1, 3, figsize=M3_FIGSIZE)
_specs = [(M3_XLIM_ZOOM,   M3_YLIM_ZOOM,   "左：低能段——peak 与能量是否成正比"),
          (M3_XLIM_GLOBAL, M3_YLIM_GLOBAL, "中：全局——一路压缩到二值硬上限"),
          (M3_XLIM_DEEP,   M3_YLIM_DEEP,   "右：穿过饱和——peak 已封顶不再变")]

for k, (xlim, ylim, ttl) in enumerate(_specs):
    a = ax[k]
    for n in N_LIST:
        S = STAT[n]
        a.plot(BOOSTS, S["peak_mu"], "-o", color=_COLOR_N[n], lw=1.6, ms=M3_MARKER_MS,
               label=f"N_shots={n}（n_tr={S['n_tr']}）")
        if M3_SHOW_CAP and k > 0:
            a.axhline(S["n_tr"], color=_COLOR_N[n], ls="-.", lw=1.0, alpha=0.55)

    if k == 0 and M3_SHOW_LINREF:
        xr = np.linspace(0, xlim[1], 50)
        for n in N_LIST:
            a.plot(xr, STAT[n]["lin_slope"] * xr, ":", color=_COLOR_N[n], lw=1.5, alpha=0.9)
        a.plot([], [], ":", color="0.3", lw=1.5, label="过原点线性外推（低能段拟合）")

    if k == 2 and M3_USE_LOG_X_DEEP:
        a.set_xscale("log")
        a.set_xlim(max(BOOSTS[BOOSTS > 0].min(), 1e-3), BOOSTS.max())
    else:
        a.set_xlim(*xlim)
    if ylim is not None:
        a.set_ylim(*ylim)
    else:
        autoscale_y(a, [(BOOSTS, STAT[n]["peak_mu"]) for n in N_LIST])
    a.set_xlabel("信号能量倍率 boost（比例量，单位无关）")
    a.set_ylabel("peak 均值 [计数/bin]")
    a.set_title(ttl, fontsize=10.5)
    a.grid(alpha=0.3)
    a.legend(fontsize=8, loc="best")

fig.suptitle(f"模块 3　图 i：peak 均值 vs 信号能量（bg = {BG:g}，每档 {N_MC:,} 次 MC，线性坐标轴）",
             fontsize=12.5)
fig.tight_layout(rect=[0, 0, 1, 0.93])
plt.show()

# ---- 拟合双轴图：boost 窗；左 peak / 右相对误差 ----
# 斜率口径由 M3_FIT_SLOPE_MODE 决定（见参数 cell）。
FIT_WIN = {}
_b_hi = float(M3_FIT_XLIM[1])
_mode = str(M3_FIT_SLOPE_MODE).strip().lower()
for n in N_LIST:
    S = STAT[n]
    m_show = (BOOSTS > 0) & (BOOSTS <= _b_hi)
    b, p = BOOSTS[m_show], S["peak_mu"][m_show]
    if _mode == "full_window":
        m_fit = m_show
    else:
        # 默认：只拿真正还线性的点定斜率，再把参考线铺满整个显示窗
        m_fit = m_show & (S["peak_mu"] <= LIN_FIT_PEAK_MAX)
    b_f, p_f = BOOSTS[m_fit], S["peak_mu"][m_fit]
    if m_fit.sum() >= 2:
        slope = float((p_f * b_f).sum() / (b_f ** 2).sum())
    else:
        slope = float(S["lin_slope"])  # 回退到模块 0 的低能斜率
    with np.errstate(invalid="ignore", divide="ignore"):
        err = np.where(b > 0, p / (slope * b) - 1.0, np.nan)
    FIT_WIN[n] = dict(mask=m_show, boost=b, peak=p, err=err, slope=slope,
                      n_fit=int(m_fit.sum()))

fig, axes = plt.subplots(1, len(N_LIST), figsize=M3_FIT_FIGSIZE, sharex=True)
if len(N_LIST) == 1:
    axes = [axes]
for ax, n in zip(axes, N_LIST):
    S, F = STAT[n], FIT_WIN[n]
    c = _COLOR_N[n]
    ax.plot(F["boost"], F["peak"], "-o", color=c, lw=1.6, ms=M3_MARKER_MS,
            label="peak 均值（MC）")
    if np.isfinite(F["slope"]):
        xr = np.linspace(0.0, _b_hi, 80)
        ax.plot(xr, F["slope"] * xr, ":", color="0.25", lw=1.8,
                label=f"过原点拟合  slope={F['slope']:.1f}")
    ax.set_xlim(*M3_FIT_XLIM)
    if M3_FIT_YLIM_PEAK is not None:
        ax.set_ylim(*M3_FIT_YLIM_PEAK)
    else:
        _ymax = float(np.nanmax(F["peak"])) if F["peak"].size else 1.0
        ax.set_ylim(0.0, max(1.0, 1.08 * _ymax))
    ax.set_xlabel("信号能量倍率 boost")
    ax.set_ylabel("peak 均值 [计数/bin]", color=c)
    ax.tick_params(axis="y", labelcolor=c)
    ax.set_title(f"N_shots={n}（n_tr={S['n_tr']}）", fontsize=10.5)
    ax.grid(alpha=0.3)

    ax2 = ax.twinx()
    if M3_FIT_SHOW_TOL:
        ax2.axhspan(-M3_LIN_TOL, M3_LIN_TOL, color="0.85", alpha=0.65, zorder=0,
                    label=f"|err|≤{M3_LIN_TOL:.0%} 参考带")
    ax2.axhline(0.0, color="0.4", ls=":", lw=1.0)
    ax2.plot(F["boost"], F["err"], "-s", color="0.2", lw=1.3, ms=M3_MARKER_MS - 0.4,
             alpha=0.9, label="相对误差")
    if M3_FIT_YLIM_ERR is not None:
        ax2.set_ylim(*M3_FIT_YLIM_ERR)
    else:
        _e = F["err"][np.isfinite(F["err"])]
        if _e.size:
            _lo, _hi = float(_e.min()), float(_e.max())
            _pad = 0.08 * max(0.2, _hi - _lo)
            ax2.set_ylim(_lo - _pad, _hi + _pad)
    ax2.set_ylabel("相对误差 = peak/(斜率×boost) − 1", color="0.2")
    ax2.tick_params(axis="y", labelcolor="0.2")

    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, fontsize=7.5, loc="upper left")

_mode_cn = ("低能点定斜率" if _mode != "full_window" else "整窗全部点定斜率")
fig.suptitle(f"模块 3　附图：boost∈[{M3_FIT_XLIM[0]:g},{M3_FIT_XLIM[1]:g}] 线性拟合"
             f"（{_mode_cn}；左=peak，右=相对误差；bg={BG:g}）", fontsize=12.0)
fig.tight_layout(rect=[0, 0, 1, 0.93])
plt.show()

print("=" * 96)
print(f"窗内线性拟合（显示 boost ≤ {M3_FIT_XLIM[1]:g}；斜率模式 = {M3_FIT_SLOPE_MODE!r}）")
print(f"  {'N':>3} {'n_tr':>5} {'斜率':>12} {'拟合点数':>8} {'窗末peak':>10} {'窗末相对误差':>14}")
for n in N_LIST:
    S, F = STAT[n], FIT_WIN[n]
    e_last = float(F["err"][-1]) if F["err"].size and np.isfinite(F["err"][-1]) else np.nan
    p_last = float(F["peak"][-1]) if F["peak"].size else np.nan
    print(f"  {n:>3} {S['n_tr']:>5} {F['slope']:>12.2f} {F['n_fit']:>8d} "
          f"{p_last:>10.2f} {e_last:>+14.2%}")

# ---- 饱和指数拟合双轴图：peak ≈ A*(1-exp(-α*boost)) ----
# 物理动机：二值轨迹点亮概率 ~ 1-e^{-κ·能量}；低能极限 A*α 应接近线性斜率。
# 不要用无界的 A*exp(b*x)——peak 有硬上限 n_tr。
FIT_EXP = {}
if M3_EXP_ENABLE:
    from scipy.optimize import curve_fit as _curve_fit

    def _sat_exp(b, A, alpha):
        return A * (1.0 - np.exp(-np.asarray(b, dtype=float) * alpha))

    _b_hi = float(M3_FIT_XLIM[1])
    _a_mode = str(M3_EXP_A_MODE).strip().lower()
    for n in N_LIST:
        S = STAT[n]
        m = (BOOSTS > 0) & (BOOSTS <= _b_hi)
        b, p = BOOSTS[m].astype(float), S["peak_mu"][m].astype(float)
        n_tr = float(S["n_tr"])
        slope0 = float(S["lin_slope"]) if np.isfinite(S["lin_slope"]) else float(p.max() / max(b.max(), 1e-12))
        alpha0 = max(slope0 / max(n_tr, 1.0), 1e-6)
        ok = False
        A_hat, alpha_hat = np.nan, np.nan
        try:
            if _a_mode == "n_tr":
                def _f1(bb, alpha):
                    return n_tr * (1.0 - np.exp(-bb * alpha))
                (alpha_hat,), _ = _curve_fit(
                    _f1, b, p, p0=[alpha0],
                    bounds=(1e-8, np.inf), maxfev=20000)
                A_hat = n_tr
            else:
                (A_hat, alpha_hat), _ = _curve_fit(
                    _sat_exp, b, p, p0=[n_tr, alpha0],
                    bounds=([1e-6, 1e-8], [np.inf, np.inf]), maxfev=20000)
            ok = True
        except Exception as _e:
            print(f"  [警告] N={n} 饱和指数拟合失败：{_e}")
        pred = _sat_exp(b, A_hat, alpha_hat) if ok else np.full_like(p, np.nan)
        with np.errstate(invalid="ignore", divide="ignore"):
            err = np.where(pred > 0, p / pred - 1.0, np.nan)
        rmse = float(np.sqrt(np.nanmean((p - pred) ** 2))) if ok else np.nan
        # 与同窗线性外推（低能斜率）比一下 RMSE，量化「指数是否更好」
        lin_pred = slope0 * b
        rmse_lin = float(np.sqrt(np.nanmean((p - lin_pred) ** 2)))
        FIT_EXP[n] = dict(
            boost=b, peak=p, pred=pred, err=err, ok=ok,
            A=float(A_hat), alpha=float(alpha_hat),
            slope_low=A_hat * alpha_hat if ok else np.nan,
            lin_slope=slope0, rmse=rmse, rmse_lin=rmse_lin, n_fit=int(m.sum()),
        )

    fig, axes = plt.subplots(1, len(N_LIST), figsize=M3_FIT_FIGSIZE, sharex=True)
    if len(N_LIST) == 1:
        axes = [axes]
    for ax, n in zip(axes, N_LIST):
        S, F = STAT[n], FIT_EXP[n]
        c = _COLOR_N[n]
        ax.plot(F["boost"], F["peak"], "-o", color=c, lw=1.6, ms=M3_MARKER_MS,
                label="peak 均值（MC）")
        if F["ok"]:
            xr = np.linspace(0.0, _b_hi, 120)
            ax.plot(xr, _sat_exp(xr, F["A"], F["alpha"]), "-", color="0.15", lw=1.9,
                    label=f"A(1-e^(-αb))  A={F['A']:.1f}, α={F['alpha']:.2f}")
        if M3_EXP_SHOW_LIN and np.isfinite(F["lin_slope"]):
            xr = np.linspace(0.0, _b_hi, 80)
            ax.plot(xr, F["lin_slope"] * xr, ":", color="0.45", lw=1.5,
                    label=f"低能线性  slope={F['lin_slope']:.1f}")
        ax.set_xlim(*M3_FIT_XLIM)
        if M3_FIT_YLIM_PEAK is not None:
            ax.set_ylim(*M3_FIT_YLIM_PEAK)
        else:
            _ymax = float(np.nanmax(F["peak"])) if F["peak"].size else 1.0
            ax.set_ylim(0.0, max(1.0, 1.08 * _ymax))
        ax.set_xlabel("信号能量倍率 boost")
        ax.set_ylabel("peak 均值 [计数/bin]", color=c)
        ax.tick_params(axis="y", labelcolor=c)
        ax.set_title(f"N_shots={n}（n_tr={S['n_tr']}）", fontsize=10.5)
        ax.grid(alpha=0.3)

        ax2 = ax.twinx()
        if M3_EXP_SHOW_TOL:
            ax2.axhspan(-M3_LIN_TOL, M3_LIN_TOL, color="0.85", alpha=0.65, zorder=0,
                        label=f"|err|≤{M3_LIN_TOL:.0%} 参考带")
        ax2.axhline(0.0, color="0.4", ls=":", lw=1.0)
        ax2.plot(F["boost"], F["err"], "-s", color="0.2", lw=1.3, ms=M3_MARKER_MS - 0.4,
                 alpha=0.9, label="相对误差（对指数）")
        _ylim_e = M3_EXP_YLIM_ERR if M3_EXP_YLIM_ERR is not None else M3_FIT_YLIM_ERR
        if _ylim_e is not None:
            ax2.set_ylim(*_ylim_e)
        else:
            _e = F["err"][np.isfinite(F["err"])]
            if _e.size:
                _lo, _hi = float(_e.min()), float(_e.max())
                _pad = 0.08 * max(0.05, _hi - _lo)
                ax2.set_ylim(_lo - _pad, _hi + _pad)
        ax2.set_ylabel("相对误差 = peak / 指数拟合 − 1", color="0.2")
        ax2.tick_params(axis="y", labelcolor="0.2")
        h1, l1 = ax.get_legend_handles_labels()
        h2, l2 = ax2.get_legend_handles_labels()
        ax.legend(h1 + h2, l1 + l2, fontsize=7.0, loc="upper left")

    _a_cn = ("A 自由拟合" if _a_mode != "n_tr" else "A 固定为 n_tr")
    fig.suptitle(f"模块 3　附图：饱和指数拟合 peak≈A(1-e^(-α·boost))"
                 f"（{_a_cn}；boost∈[{M3_FIT_XLIM[0]:g},{M3_FIT_XLIM[1]:g}]；bg={BG:g}）",
                 fontsize=12.0)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    plt.show()

    print("=" * 96)
    print(f"饱和指数拟合 peak≈A(1-exp(-α·boost))（boost≤{M3_FIT_XLIM[1]:g}；A 模式={M3_EXP_A_MODE!r}）")
    print(f"  {'N':>3} {'n_tr':>5} {'A':>10} {'α':>10} {'A·α(低能斜率)':>14} "
          f"{'低能线性斜率':>12} {'RMSE_指数':>10} {'RMSE_线性':>10} {'窗末相对误差':>14}")
    for n in N_LIST:
        F = FIT_EXP[n]
        e_last = float(F["err"][-1]) if F["err"].size and np.isfinite(F["err"][-1]) else np.nan
        print(f"  {n:>3} {STAT[n]['n_tr']:>5} {F['A']:>10.2f} {F['alpha']:>10.3f} "
              f"{F['slope_low']:>14.2f} {F['lin_slope']:>12.2f} "
              f"{F['rmse']:>10.3f} {F['rmse_lin']:>10.3f} {e_last:>+14.2%}")

# ---- 附图：相对线性外推的误差 vs peak，标出 |err|≤5% 的连续区间上限 ----
# 口径：从低能往上，连续满足 |peak/(斜率×boost)−1| ≤ M3_LIN_TOL 的最后一档；
# 一旦越界就停，后面即使偶然再落回容差内也不算（避免饱和回折误判）。
# 注意：这里的斜率仍用模块 0 低能段（peak≤LIN_FIT_PEAK_MAX）拟合的 STAT['lin_slope']，
# 与上方「整窗拟合」不同——整窗会被饱和段压低斜率，不适合用来划「可当线性」上限。
LIN_OK = {}
for n in N_LIST:
    S = STAT[n]
    ref = S["lin_slope"] * BOOSTS
    with np.errstate(invalid="ignore", divide="ignore"):
        err = np.where(BOOSTS > 0, S["peak_mu"] / ref - 1.0, np.nan)
    ok = np.isfinite(err) & (np.abs(err) <= M3_LIN_TOL)
    # 从第一个有效点起连续 True 的末尾
    i_last = -1
    for i in range(err.size):
        if BOOSTS[i] <= 0 or not np.isfinite(err[i]):
            continue
        if ok[i]:
            i_last = i
        else:
            break
    LIN_OK[n] = dict(
        err=err,
        i_last=i_last,
        peak_max=float(S["peak_mu"][i_last]) if i_last >= 0 else np.nan,
        boost_max=float(BOOSTS[i_last]) if i_last >= 0 else np.nan,
        frac=float(S["peak_mu"][i_last] / S["n_tr"]) if i_last >= 0 else np.nan,
    )

fig, ax = plt.subplots(figsize=M3_LIN_FIGSIZE)
ax.axhspan(-M3_LIN_TOL, M3_LIN_TOL, color="0.85", alpha=0.7, zorder=0,
           label=f"|相对误差| ≤ {M3_LIN_TOL:.0%} 容差带")
ax.axhline(0.0, color="0.35", ls=":", lw=1.2)
for n in N_LIST:
    S = STAT[n]
    L = LIN_OK[n]
    m = BOOSTS > 0
    ax.plot(S["peak_mu"][m], L["err"][m], "-o", color=_COLOR_N[n],
            lw=1.6, ms=M3_MARKER_MS,
            label=f"N_shots={n}（n_tr={S['n_tr']}）")
    if L["i_last"] >= 0:
        ax.axvline(L["peak_max"], color=_COLOR_N[n], ls="--", lw=1.3, alpha=0.85)
        ax.plot(L["peak_max"], L["err"][L["i_last"]], "*", color=_COLOR_N[n],
                ms=14, mec="k", mew=0.6, zorder=5)
ax.set_xlabel("peak 均值 [计数/bin]")
ax.set_ylabel("相对误差 = peak / (斜率×boost) − 1")
ax.set_title(f"模块 3　附图：相对低能斜率的误差（容差 {M3_LIN_TOL:.0%}；"
             f"虚线+★ = 连续满足容差的 peak 上限）", fontsize=10.5)
if M3_LIN_XLIM_PEAK is not None:
    ax.set_xlim(*M3_LIN_XLIM_PEAK)
else:
    # 横轴铺到略超过最远的 5% 越界点，让上限位置清楚；别被饱和段拖到 n_tr
    _xmax = max((L["peak_max"] for L in LIN_OK.values() if np.isfinite(L["peak_max"])),
                default=10.0)
    ax.set_xlim(0.0, max(2.0, 1.6 * _xmax))
if M3_LIN_YLIM_ERR is not None:
    ax.set_ylim(*M3_LIN_YLIM_ERR)
ax.grid(alpha=0.3)
ax.legend(fontsize=8, loc="best")
fig.tight_layout()
plt.show()

print("=" * 96)
print("低能段线性度体检：peak 均值 / (斜率 × boost)，等于 1 就是严格正比")
print(f"  {'boost':>10}" + "".join(f"{'N=' + str(n):>22}" for n in N_LIST))
for i, b in enumerate(BOOSTS):
    if b <= 0 or b > M3_XLIM_ZOOM[1]:
        continue
    row = f"  {b:>10.5f}"
    for n in N_LIST:
        S = STAT[n]
        mu, ref = S["peak_mu"][i], S["lin_slope"] * b
        row += f"{mu:>12.3f}（{mu / ref if ref > 0 else np.nan:>5.3f}）"
    print(row)

print("\n" + "=" * 96)
print(f"|相对误差| ≤ {M3_LIN_TOL:.0%} 的连续线性区（从低能往上，一旦越界即停；"
      f"斜率=模块0低能段拟合）")
print(f"  {'N':>3} {'n_tr':>5} {'斜率':>10} {'peak上限':>10} {'占硬上限':>10} "
      f"{'对应boost':>12} {'该档相对误差':>14}")
for n in N_LIST:
    S, L = STAT[n], LIN_OK[n]
    if L["i_last"] < 0:
        print(f"  {n:>3} {S['n_tr']:>5} {S['lin_slope']:>10.2f} {'—':>10}")
        continue
    e = L["err"][L["i_last"]]
    print(f"  {n:>3} {S['n_tr']:>5} {S['lin_slope']:>10.2f} {L['peak_max']:>10.2f} "
          f"{L['frac']:>10.1%} {L['boost_max']:>12.5g} {e:>+14.2%}")
""")

# =====================================================================
md(r"""
## 模块 4　图 ii 与图 iii：peak 的标准差

`peak` 的标准差随能量**不是单调的**，而是拱形：弱信号端被「几乎不可能亮」压扁，
强信号端被「几乎必定亮」压扁，中间涨落最大。

**照度均匀时，这条拱就是教科书里的二项曲线。** 27 个 SPAD 的收集比例相同 ⇒
同一个 bin 里 $n_{tr}$ 条轨迹的点亮概率 $p_t$ 严格相等 ⇒ 单 bin 计数服从二项分布 $B(n_{tr},p)$：

$$\sigma = \sqrt{n_{tr}\,p\,(1-p)},\qquad \max_p \sigma = \tfrac12\sqrt{n_{tr}}\ \text{ 于 }\ p=0.5$$

实测拱顶落在占比 **0.53**、拱顶幅度与 $\tfrac12\sqrt{n_{tr}}$ 的比值是 **0.993 / 1.004 / 1.002**
（N = 1 / 2 / 4），逐档比值在占比 0.05–0.95 全程落在 **0.99–1.005**（见下方打印表）。
左图里实测曲线整条压在二项参考的粗浅光晕上，看不出偏差。
这和纯噪声侧 PoD_esti 用的 $\sigma^2 = bg\,(1-bg/n_{tr})$ 是**同一条规律**——
因为两边现在都满足「各轨迹概率相同」这个前提。

**`peak` 是 152 个 bin 的最大值，不是单个 bin 的计数**，所以严格说它不服从二项分布。
但在 **bg = 0** 这个前提下，差别小到测不出来：对照实验（`check_uniform_vs_real_fpix.py`）
直接比 `peak 的 σ / 最亮单 bin 的 σ`，**全能量段都是 0.978–1.001**，均匀与真实像斑都一样。
道理很直白：没有环境光时统计窗里只有一处有信号，其余 bin 恒为 0，
**根本没有第二个候选者来竞争最大值**，于是 `peak` 就等于最亮那个 bin 的计数。
把 bg 开成非零之后这一条不再成立，届时需要重新评估。

注意占比 > 0.95 之后不要再用「实测 σ / 由 peak 均值代入的二项 σ」这个比值判断——
那里 peak 均值已经贴死 $n_{tr}$，把它当二项的 $p$ 代回去，分母趋于 0，比值失去意义。

**换成真实像斑会怎样**（对照缓存 `peak_vs_energy_cache_real.npz`）：`f_pix` 最大/最小相差约 241 倍，
$p_t$ 互不相同，单 bin 变成**泊松二项分布**，

$$\sigma^2=\sum_t p_t(1-p_t)=n_{tr}\,\bar p\,(1-\bar p)-n_{tr}\operatorname{Var}_t(p_t)$$

右边第二项恒为正 ⇒ 方差恒低于二项，拱顶被推到占比 0.2–0.4 那一段。
**本文件用的是均匀照度，不受这一项影响。**

- **左（图 ii）**：标准差 vs peak 均值，★ 标出拱顶，阴影标出 σ 在峰值 97% 以内的平台段，
  点线是纯二项参考 $\sqrt{n_{tr}p(1-p)}$。
- **中 / 右（图 iii）**：标准差 vs 能量，分别看低能段与全局。
""")

code(r"""
# ---------------- 模块 4 画图参数 ----------------
M4_XLIM_PEAKMU  = None          # 【左图横轴】peak 均值范围；None = 自动（0 → 最大 n_tr）
M4_XLIM_ZOOM    = (0.0, 0.05)   # 【中图横轴】低能段 boost 范围。照度均匀时拱顶（占比 0.5）
                                #   落在 boost≈0.03，正好在这一段中间
M4_XLIM_GLOBAL  = (0.0, 0.5)    # 【右图横轴】全局 boost 范围
M4_YLIM_SD      = None          # 【三张图共用纵轴】peak 标准差范围；None = 自动
M4_NORMALIZE_X  = False         # 【左图横轴】True = 把 peak 均值除以 n_tr 归一化成占比，
                                #   便于验证「拱顶在 p≈0.5」；False = 直接用计数
M4_SHOW_BINOM   = True          # 【左图】是否叠纯二项参考曲线 sqrt(n_tr·p(1-p))
M4_BINOM_LW     = 5.0           # 【左图参考线】线宽。均匀照度下实测与参考重合，参考线画粗、
M4_BINOM_ALPHA  = 0.28          #   画浅、压在底层当光晕，实测细线盖在上面才看得出是两条
M4_SHOW_PLATEAU = True          # 【左图】是否用阴影标出「拱顶平台」（σ ≥ 峰值的 M4_PLATEAU_FRAC）。
                                #   真实 f_pix 下拱很平，★ 的位置基本由 MC 噪声决定，不要单看 ★
M4_PLATEAU_FRAC = 0.97          # 【左图阴影】平台判定阈值。20,000 次 MC 下 σ 自身的统计误差
                                #   约 σ/√(2N) ≈ 0.5%，取 97% 大致框出「实际难以区分」的一段
M4_MARKER_MS    = 3.4
M4_FIGSIZE      = (16.8, 5.2)
""")

code(r"""
fig, ax = plt.subplots(1, 3, figsize=M4_FIGSIZE)

# ---- 左：图 ii　标准差 vs peak 均值 ----
for n in N_LIST:
    S = STAT[n]
    x = S["peak_mu"] / S["n_tr"] if M4_NORMALIZE_X else S["peak_mu"]
    ax[0].plot(x, S["peak_sd"], "-o", color=_COLOR_N[n], lw=1.6, ms=M4_MARKER_MS,
               label=f"N_shots={n}（n_tr={S['n_tr']}）")
    i_top = int(np.nanargmax(S["peak_sd"]))
    ax[0].plot(x[i_top], S["peak_sd"][i_top], "*", color=_COLOR_N[n], ms=15,
               mec="k", mew=0.6, zorder=5)
    # 拱顶其实是一段平台：σ 在峰值 99% 以内的档，argmax 落在哪一档基本由 MC 噪声决定。
    # 画出平台范围，免得 ★ 被当成一个精确的位置。
    if M4_SHOW_PLATEAU:
        pl = np.where(S["peak_sd"] >= M4_PLATEAU_FRAC * S["peak_sd"][i_top])[0]
        if pl.size > 1:
            ax[0].axvspan(x[pl].min(), x[pl].max(), color=_COLOR_N[n], alpha=0.10, lw=0)
if M4_SHOW_BINOM:
    # 照度均匀时实测曲线与二项参考几乎完全重合，画成同粗细的两条线就只看得到一条。
    # 参考线改成压在底下的粗浅色「光晕」，实测细线盖在上面，重合与否一眼可辨。
    for n in N_LIST:
        S = STAT[n]
        p = np.linspace(0, 1, 400)
        xr = p if M4_NORMALIZE_X else p * S["n_tr"]
        ax[0].plot(xr, np.sqrt(S["n_tr"] * p * (1 - p)), "-", color=_COLOR_N[n],
                   lw=M4_BINOM_LW, alpha=M4_BINOM_ALPHA, zorder=1, solid_capstyle="round")
    ax[0].plot([], [], "-", color="0.45", lw=M4_BINOM_LW, alpha=M4_BINOM_ALPHA,
               label=r"纯二项 $\sqrt{n_{tr}p(1-p)}$（底层粗浅线）")
ax[0].set_xlabel("peak 均值占硬上限的比例" if M4_NORMALIZE_X else "peak 均值 [计数/bin]")
ax[0].set_title("左（图 ii）：peak 标准差 vs peak 均值（★ 为拱顶）", fontsize=10.5)
if M4_XLIM_PEAKMU is not None:
    ax[0].set_xlim(*M4_XLIM_PEAKMU)

# ---- 中 / 右：图 iii　标准差 vs 能量 ----
for k, (xlim, ttl) in enumerate([(M4_XLIM_ZOOM, "中（图 iii）：标准差 vs 能量——低能段"),
                                 (M4_XLIM_GLOBAL, "右（图 iii）：标准差 vs 能量——全局")]):
    a = ax[k + 1]
    for n in N_LIST:
        a.plot(BOOSTS, STAT[n]["peak_sd"], "-o", color=_COLOR_N[n], lw=1.6,
               ms=M4_MARKER_MS, label=f"N_shots={n}")
    a.set_xlim(*xlim)
    a.set_xlabel("信号能量倍率 boost（比例量，单位无关）")
    a.set_title(ttl, fontsize=10.5)
    if M4_YLIM_SD is None:
        autoscale_y(a, [(BOOSTS, STAT[n]["peak_sd"]) for n in N_LIST])

for a in ax:
    a.set_ylabel("peak 标准差 [计数]")
    if M4_YLIM_SD is not None:
        a.set_ylim(*M4_YLIM_SD)
    a.grid(alpha=0.3)
    a.legend(fontsize=8, loc="best")

fig.suptitle(f"模块 4　图 ii / iii：peak 标准差（bg = {BG:g}，每档 {N_MC:,} 次 MC，线性坐标轴）",
             fontsize=12.5)
fig.tight_layout(rect=[0, 0, 1, 0.93])
plt.show()

print("=" * 96)
print("peak 标准差的拱顶 —— 与纯二项 √(n_tr·p(1−p)) 对照")
print(f"  当前照度模式 f_pix = '{FPIX_MODE}'，离散度 std/mean = {FPIX_CV:.3f}"
      + ("（均匀 ⇒ 各轨迹 p_t 相等 ⇒ 单 bin 服从二项分布，拱顶应在占比 0.5、峰值 √(n_tr)/2）"
         if FPIX_CV < 1e-9 else
         "（不均匀 ⇒ 泊松二项，方差被压低 n_tr·Var_t(p_t)，拱顶左移）"))
print(f"  {'N':>3} {'n_tr':>5} {'拱顶σ':>8} {'纯二项峰值':>11} {'实测/纯二项':>12} "
      f"{'argmax占比':>11} {'平台占比区间':>16}")
for n in N_LIST:
    S = STAT[n]
    i = int(np.nanargmax(S["peak_sd"]))
    sd_binom = np.sqrt(S["n_tr"] * 0.25)
    frac = S["peak_mu"] / S["n_tr"]
    pl = np.where(S["peak_sd"] >= M4_PLATEAU_FRAC * S["peak_sd"][i])[0]
    span = f"{frac[pl].min():.2f} – {frac[pl].max():.2f}" if pl.size > 1 else "—"
    print(f"  {n:>3} {S['n_tr']:>5} {S['peak_sd'][i]:>8.2f} {sd_binom:>11.2f} "
          f"{S['peak_sd'][i] / sd_binom:>12.3f} {frac[i]:>11.3f} {span:>16}")
print(f"\n  注意：拱顶是一段平台，不是一个点——末列区间内 σ 只变化 {1 - M4_PLATEAU_FRAC:.0%} 以内，")
print("  argmax 的具体落点对 MC 噪声与能量网格都敏感，引用时不要报三位有效数字。")
print("  另：三档 N 的数值高度一致不是独立证据——它们由同一批 4 发实现前缀和派生，随机性共享。")
print("\n  逐档核对「实测 σ / 二项 σ」（只列占比 0.05–0.95 之间的档，两端 σ 太小、比值没意义）：")
_n = N_LIST[-1]
_S = STAT[_n]
_fr = _S["peak_mu"] / _S["n_tr"]
_m = (_fr > 0.05) & (_fr < 0.95)
print(f"  N_shots={_n}: " + "  ".join(
    f"{_fr[i]:.2f}→{_S['peak_sd'][i] / np.sqrt(_S['n_tr'] * _fr[i] * (1 - _fr[i])):.3f}"
    for i in np.where(_m)[0][::2]))
print("\n  全程贴着 1 ⇒ 照度均匀时 peak 的涨落就是二项涨落。")
print("  占比 > 0.95 的档没有列：那里 peak 均值已贴死 n_tr，拿它当二项的 p 代回去分母趋于 0，")
print("  这个比值失去意义。")
print("  另有对照实验直接比 peak 与最亮单 bin（check_uniform_vs_real_fpix.py）：")
print("  σ_peak/σ_单bin 全能量段 0.978–1.001 —— bg=0 时窗内只有一处有信号、没有第二个候选者，")
print("  所以「152 bin 取最大值」这一步自始至终不起作用。开了 bg 之后需重新评估。")
""")

# =====================================================================
md(r"""
## 模块 5　图 iv 与图 v：半高全宽

半高全宽有两个和直觉不同的地方：

**一是有地板，约等于 `T_OVER` ≈ 8 ns。** 一次雪崩就会把连续约 8 个 bin 点亮，
所以哪怕信号只有一个光子，累加波形的半高全宽也降不到 8 ns 以下。
低能段量到的宽度是**器件的过阈窗宽**，不是激光脉冲宽度。

**二是不会饱和。** 照度均匀时 27 个 SPAD 同步逼近饱和，`peak` 在 `boost ≈ 0.3` 就彻底封顶
（真实像斑要等最暗的 SPAD 点亮，得拖到 `boost ≈ 20`），但 FWHM 之后仍然按
「每十倍能量增宽约固定 ns 数」一直长下去——因为峰高被硬上限锁死后，
脉冲的前后沿被继续抬高到半高线以上，过半高的时间跨度就随能量的对数增长。

**三是几乎与 N_shots 无关。** N = 1 / 2 / 4 三条曲线在左、中两图上基本重合——
累加更多发只把峰抬高，不改变波形的相对形状，所以半高全宽只由能量决定。
正因为重合，三条线必须用不同粗细与线型才看得出来（见 `M5_STYLE_N`）。

灰色阴影标出 FWHM 有效样本占比不足 50% 的能量区间（该处多数实现的 peak 还够不到 4 计数），
这一段的 FWHM 只由少数较亮的实现贡献，有选择偏倚，不要当结论用。
""")

code(r"""
# ---------------- 模块 5 画图参数 ----------------
M5_XLIM_GLOBAL    = (0.0, 0.5)    # 【左图横轴】全局 boost 范围（照度均匀时 boost≈0.3 即封顶）
M5_XLIM_DEEP      = (0.0, 200.0)  # 【中图横轴】穿过饱和的 boost 范围
M5_XLIM_PEAKMU    = None          # 【右图横轴】peak 均值范围；None = 自动
M5_YLIM_FWHM      = None          # 【三张图共用纵轴】FWHM [ns]；None = 自动
M5_SHOW_BAND      = True          # 【是否叠】FWHM 的 ±1 标准差阴影带（逐次实现之间的涨落）
M5_STYLE_N = {1: dict(ls="-",  lw=3.0, ms=0.0, alpha=0.9),    # 【线型】三档 N_shots 的 FWHM 曲线
              2: dict(ls="--", lw=2.0, ms=0.0, alpha=0.95),   #   几乎完全重合，必须靠粗细+线型区分，
              4: dict(ls="-",  lw=1.1, ms=3.2, alpha=1.0)}    #   否则后画的会把先画的盖死
M5_SHOW_TOVER     = True          # 【是否画】T_OVER 地板水平线
M5_MIN_VALID      = 0.5           # 【灰色阴影】有效样本占比低于此值的能量段标为不可信
M5_SHOW_DEEP_FIT  = True          # 【中图】是否叠深饱和段的对数拟合直线并标注每十倍增宽
M5_USE_LOG_X_DEEP = False         # 【中图横轴】True = 改成对数轴，可一眼看到对数直线
M5_MARKER_MS      = 3.4
M5_FIGSIZE        = (16.8, 5.2)
""")

code(r"""
fig, ax = plt.subplots(1, 3, figsize=M5_FIGSIZE)
_TOVER_NS = core.T_OVER * 1e9


def _shade_invalid(a, xlim):
    # 把 FWHM 有效样本占比不足的能量段涂灰，提醒该段有选择偏倚
    v = STAT[N_LIST[-1]]["fwhm_valid"]
    bad = BOOSTS[(v < M5_MIN_VALID) & (BOOSTS <= xlim[1])]
    if bad.size:
        a.axvspan(xlim[0], min(bad.max(), xlim[1]), color="0.85", alpha=0.55, zorder=0,
                  label=f"有效样本 < {M5_MIN_VALID:.0%}，勿引用")


# ---- 左（图 iv）：FWHM vs 能量，全局 ----
# ---- 中（图 iv）：FWHM vs 能量，穿过饱和 ----
for k, (xlim, ttl) in enumerate([(M5_XLIM_GLOBAL, "左（图 iv）：FWHM vs 能量——全局"),
                                 (M5_XLIM_DEEP, "中（图 iv）：FWHM vs 能量——穿过饱和后仍在增宽")]):
    a = ax[k]
    _shade_invalid(a, xlim)
    for n in N_LIST:
        S = STAT[n]
        a.plot(BOOSTS, S["fwhm_mu"], marker="o", color=_COLOR_N[n],
               label=f"N_shots={n}", **M5_STYLE_N[n])
        if M5_SHOW_BAND:
            a.fill_between(BOOSTS, S["fwhm_mu"] - S["fwhm_sd"], S["fwhm_mu"] + S["fwhm_sd"],
                           color=_COLOR_N[n], alpha=0.14, lw=0)
    if k == 1 and M5_SHOW_DEEP_FIT:
        n = N_LIST[-1]
        S = STAT[n]
        d = BOOSTS >= DEEP_BOOST_MIN
        xr = np.linspace(DEEP_BOOST_MIN, max(xlim[1], BOOSTS.max()), 200)
        a.plot(xr, S["fwhm_a"] + S["fwhm_per_decade"] * np.log10(xr), "--", color="k",
               lw=1.4, label=f"深饱和拟合：每十倍能量 +{S['fwhm_per_decade']:.2f} ns")
    if k == 1 and M5_USE_LOG_X_DEEP:
        a.set_xscale("log")
        a.set_xlim(max(BOOSTS[BOOSTS > 0].min(), 1e-3), BOOSTS.max())
    else:
        a.set_xlim(*xlim)
    a.set_xlabel("信号能量倍率 boost（比例量，单位无关）")
    a.set_title(ttl, fontsize=10.5)
    if M5_YLIM_FWHM is None:
        autoscale_y(a, [(BOOSTS, STAT[n]["fwhm_mu"]) for n in N_LIST]
                    + [(BOOSTS, np.full_like(BOOSTS, _TOVER_NS))], floor_zero=False)

# ---- 右（图 v）：FWHM vs peak 均值 ----
a = ax[2]
for n in N_LIST:
    S = STAT[n]
    a.plot(S["peak_mu"], S["fwhm_mu"], marker="o", color=_COLOR_N[n],
           label=f"N_shots={n}（n_tr={S['n_tr']}）", **M5_STYLE_N[n])
    if M5_SHOW_BAND:
        a.fill_between(S["peak_mu"], S["fwhm_mu"] - S["fwhm_sd"], S["fwhm_mu"] + S["fwhm_sd"],
                       color=_COLOR_N[n], alpha=0.14, lw=0)
    a.axvline(S["n_tr"], color=_COLOR_N[n], ls="-.", lw=1.0, alpha=0.5)
if M5_XLIM_PEAKMU is not None:
    a.set_xlim(*M5_XLIM_PEAKMU)
a.set_xlabel("peak 均值 [计数/bin]")
a.set_title("右（图 v）：FWHM vs peak 均值——末端竖直上翘即封顶后继续增宽", fontsize=10.5)

for a in ax:
    if M5_SHOW_TOVER:
        a.axhline(_TOVER_NS, color="0.35", ls=":", lw=1.6,
                  label=f"T_OVER 地板 = {_TOVER_NS:.0f} ns")
    a.set_ylabel("半高全宽 FWHM [ns]")
    if M5_YLIM_FWHM is not None:
        a.set_ylim(*M5_YLIM_FWHM)
    a.grid(alpha=0.3)
    a.legend(fontsize=7.8, loc="best")

fig.suptitle(f"模块 5　图 iv / v：半高全宽（bg = {BG:g}，逐次实现测量后对 {N_MC:,} 次 MC 求平均，"
             f"线性坐标轴）", fontsize=12.5)
fig.tight_layout(rect=[0, 0, 1, 0.93])
plt.show()
""")

# =====================================================================
md(r"""
## 模块 6　关键数值汇总

把全部能量档的 peak 与 FWHM 打成表。**深饱和段线性横轴画不下（4 个数量级），
所以那一段的定量结论以本表为准。**
""")

code(r"""
print("=" * 108)
print(f"信号能量扫描汇总　bg = {BG:g}，每档 {N_MC:,} 次 MC，宏像元 {core.N_PIX_MACRO} SPAD")
print("=" * 108)
for n in N_LIST:
    S = STAT[n]
    print(f"\n--- N_shots = {n}（二值硬上限 n_tr = {S['n_tr']}）---")
    print(f"  {'boost':>11} {'peak均值':>9} {'占上限':>8} {'peak标准差':>11} "
          f"{'FWHM均值[ns]':>13} {'FWHM标准差':>11} {'FWHM有效%':>10}")
    for i, b in enumerate(BOOSTS):
        if b <= 0:
            continue
        print(f"  {b:>11.5g} {S['peak_mu'][i]:>9.2f} {S['peak_mu'][i] / S['n_tr']:>8.3f} "
              f"{S['peak_sd'][i]:>11.2f} {S['fwhm_mu'][i]:>13.2f} {S['fwhm_sd'][i]:>11.2f} "
              f"{100 * S['fwhm_valid'][i]:>9.1f}%")

print("\n" + "=" * 108)
print("结论速览")
print("=" * 108)
for n in N_LIST:
    S = STAT[n]
    i_top = int(np.nanargmax(S["peak_sd"]))
    i_lin = np.where((BOOSTS > 0) & (S["peak_mu"] > 0.9 * S["lin_slope"] * BOOSTS))[0]
    b_lin = BOOSTS[i_lin].max() if i_lin.size else np.nan
    d = np.isfinite(S["fwhm_mu"]) & (S["fwhm_valid"] >= 0.5)
    print(f"  N_shots={n}（n_tr={S['n_tr']}）：")
    print(f"    peak 与能量保持 10% 以内线性：boost 到 {b_lin:.4g}（peak 均值约 "
          f"{S['lin_slope'] * b_lin:.1f} 计数，占上限 {S['lin_slope'] * b_lin / S['n_tr']:.1%}）")
    print(f"    peak 封顶值 {S['peak_mu'].max():.2f}（占上限 {S['peak_mu'].max() / S['n_tr']:.1%}）；"
          f"标准差拱顶 {S['peak_sd'][i_top]:.2f} 出现在占上限 {S['peak_mu'][i_top] / S['n_tr']:.2f} 处")
    print(f"    FWHM 从 {np.nanmin(S['fwhm_mu'][d]):.2f} ns（T_OVER 地板）涨到 "
          f"{np.nanmax(S['fwhm_mu'][d]):.2f} ns，深饱和段每十倍能量 +{S['fwhm_per_decade']:.2f} ns")
""")


# =====================================================================
def main() -> None:
    for i, (kind, src) in enumerate(CELLS):
        if kind == "code":
            try:
                ast.parse(src)
            except SyntaxError as e:
                raise SystemExit(f"cell {i} 语法错误：{e}") from e

    nb = {
        "cells": [
            {"cell_type": k, "metadata": {},
             **({"outputs": [], "execution_count": None} if k == "code" else {}),
             "source": (s + "\n").splitlines(keepends=True)}
            for k, s in CELLS
        ],
        "metadata": {
            "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
            "language_info": {"name": "python", "version": "3.11"},
        },
        "nbformat": 4, "nbformat_minor": 5,
    }
    with open(NB_PATH, "w", encoding="utf-8") as f:
        json.dump(nb, f, ensure_ascii=False, indent=1)

    nmd = sum(1 for k, _ in CELLS if k == "markdown")
    print(f"[已生成] {NB_PATH}：{len(CELLS)} cells（markdown {nmd} / code {len(CELLS) - nmd}）")
    for i, (k, s) in enumerate(CELLS):
        head = s.lstrip("#").strip().splitlines()[0][:64]
        print(f"  {i:>3} {k:<8} {len(s.splitlines()):>4}L  {head}")


if __name__ == "__main__":
    main()
