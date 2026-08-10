# -*- coding: utf-8 -*-
"""由 PoD_esti_v11.ipynb 生成 PoD_esti_v20.ipynb。

v20 = v11 的**全部内容一个不删**，外加 6 个新模块（11–16）和 1 个理论汇总（17）：

  模块 11  每条 hist 内 152 bin 的 std 均值 / peak 均值 / peak 标准差 随 bg
  模块 12  连续（实数）阈值曲线 —— 折线，不再是整数阶梯
  模块 13  FAR=5% / 1% / 100 ppm 下 PoD50 与 PoD90 所需信号的均值
  模块 14  平方反比测远（纯 1/D² 与含大气衰减两种口径）
  模块 15  同信号强度、不同 bg 时 peak 分布怎么变（均值是否只是加 bg、std 怎么变）
  模块 16  宏像元 3×9 vs 3×6 阈值对比（本轮旁支分析并入正本）
  模块 17  理论汇总（阈值倍数模型 + 引擎一致性模型 + 两种口径的说明）

同时对 v11 主体做的**最小改动**（不改物理参数）：
  * 缓存名换成 pod_esti_v20_*，并把 v11 的缓存登记为 fallback（规则三：读到即同步写回 v20 主名），
    所以已经跑完的 v11 噪声/PoD 结果不会重算。
  * 模块 9.3 的信号扫描从 notebook 内【串行】改为调用多进程脚本 run_pod_v20_sig_scan.py。
  * 自动开跑的脚本名 v11 → v20；输出图名 pod_v11_* → pod_v20_*。

用法：python upgrade_pod_esti_v20_from_v11.py
"""
from __future__ import annotations

import json
import shutil
from pathlib import Path

SRC = Path("PoD_esti_v11.ipynb")
DST = Path("PoD_esti_v20.ipynb")


_ID = [0]


def _next_id() -> str:
    """nbformat 4.5 要求每个 cell 有唯一 id。"""
    _ID[0] += 1
    return f"v20cell{_ID[0]:03d}"


def md(src: str) -> dict:
    return {"cell_type": "markdown", "id": _next_id(), "metadata": {},
            "source": src.strip("\n").splitlines(keepends=True)}


def code(src: str) -> dict:
    return {"cell_type": "code", "id": _next_id(), "execution_count": None,
            "metadata": {}, "outputs": [],
            "source": src.strip("\n").splitlines(keepends=True)}


def cell_src(nb, i) -> str:
    return "".join(nb["cells"][i]["source"])


def set_src(nb, i, s: str) -> None:
    nb["cells"][i]["source"] = s.splitlines(keepends=True)


def sub(nb, i, old, new, required=True):
    s = cell_src(nb, i)
    if old not in s:
        if required:
            raise SystemExit(f"cell {i} 找不到待替换文本：{old[:70]!r}")
        return
    set_src(nb, i, s.replace(old, new))


# ============================================================ 新模块源码
M11_MD = r"""
## 模块 11 — ★ v20：每条 hist 的 std、peak 均值、peak 标准差 随 bg

这一模块补齐"**对每个 bg 档，除了阈值以外还要看的三个统计量**"：

| # | 量 | 精确定义 |
|---|---|---|
| ① | **单条 hist 内的 std** | 取**一条** `hist_add` 在统计窗 152 个 bin 上的样本标准差（`ddof=1`），再对所有 MC 条数求平均。注意这是"**一条直方图内部**的起伏"，不是"同一个 bin 在多条直方图之间的起伏" —— 两者只有在 bin 之间独立时才相等 |
| ② | **peak 均值** | 统计窗内最大 bin 的均值 |
| ③ | **peak 标准差** | 统计窗内最大 bin 跨 MC 条数的标准差 |

### ① 单条 hist 内 std 的解析对照（两条线）

* **纯泊松**：$\sigma=\sqrt{\mathrm{bg}}$
* **二值饱和（Binomial）**：每个 bin 是 $\mathrm{Bin}(n_{tr},p_{eq})$，$n_{tr}=27N$，$p_{eq}=\mathrm{bg}/n_{tr}$，故
  $$\sigma_{bin}=\sqrt{n_{tr}\,p_{eq}(1-p_{eq})}=\sqrt{\mathrm{bg}\left(1-\frac{\mathrm{bg}}{27N}\right)}$$
  这就是 **Fano 因子 $F=1-\mathrm{bg}/(27N)<1$** 的来源：SPAD 是 1 bit 的，一个 bin 最多点亮一次，
  所以噪声比泊松**更小**（亚泊松）。$N$ 越大 $n_{tr}$ 越大，压缩越弱，$\sigma_{bin}$ 越接近 $\sqrt{\mathrm{bg}}$。

### ③ peak 标准差的解析对照（Gumbel 极值公式）

peak 是统计窗内约 $M_{eff}$ 个近高斯 bin（每个 $\approx\mathcal N(\mathrm{bg},\sigma_{bin}^2)$）取极大值。
高斯极大值服从 **Gumbel** 极限，标准差有闭式：

$$
\sigma_{peak}\;\approx\;\frac{\pi}{\sqrt6}\,\frac{\sigma_{bin}}{z_M},
\qquad z_M=\sqrt{2\ln M_{eff}}\;\approx\;\frac{\mu_{peak}-\mathrm{bg}}{\sigma_{bin}}
$$

这里用**实测的** $\mu_{peak}$ 反推 $z_M$（峰位在均值上方几个 $\sigma_{bin}$），**不含拟合参数**，
所以它是一个真正的检验：「已知峰位，峰宽是不是极值理论给的宽度」。
它把三张图串起来 —— ①给 $\sigma_{bin}$、②给 $\mu_{peak}$、③预测 $\sigma_{peak}$。

> **适用范围**：Gumbel 是 $M\to\infty$ 的渐近结果。大 bg（$\mathrm{bg}\gtrsim6$）计数多、单 bin 近高斯，
> 解析与实测贴合到几个百分点；小 bg 时计数少、分布离散且强右偏，解析会系统性偏低，属正常。

数据来自 `scan_hist_std_peak.py`（缓存 `scan_hist_std_peak_cache.npz`，48 bg × N∈{1,2,4} × 100,000 MC）。
横轴给两个版本：**vs bg**（同底噪比 N）与 **vs noise**（同单发底噪比 N）。
"""

M11_CODE = r'''
# ---- 模块 11：单条 hist 内 std / peak 均值 / peak std 随 bg ----
_COLORS_N   = {1: "tab:blue", 2: "tab:green", 4: "tab:red"}   # v20 各新模块统一配色
HSP_CACHE   = "scan_hist_std_peak_cache.npz"
HSP_SCRIPT  = "scan_hist_std_peak.py"
N_MC_HSP    = 100_000


def _load_hist_std_cache(path, grid, n_mc, n_list):
    if not os.path.exists(path):
        return None
    try:
        z = np.load(path)
        if int(z["n_mc"]) != int(n_mc) or z["grid"].shape != grid.shape \
                or not np.allclose(z["grid"], grid):
            return None
        out = {n: {k: np.array(z[f"{k}_{n}"]) for k in
                   ("hist_std", "bg_mc", "peak_mean", "peak_std", "done")}
               for n in n_list}
        return out if all(np.all(out[n]["done"]) for n in n_list) else None
    except Exception:
        return None


HSP = _load_hist_std_cache(HSP_CACHE, np.asarray(BG_GRID, float),
                           N_MC_HSP, N_SHOTS_LIST)
if HSP is None:
    import sys
    print("=" * 72)
    print(f"未找到完整 {HSP_CACHE} → 自动调用多进程扫描")
    print("=" * 72)
    _rc = _run_cmd_stream([sys.executable, HSP_SCRIPT,
                           "--workers", str(int(N_WORKERS)),
                           "--n-mc", str(int(N_MC_HSP))])
    if _rc != 0:
        raise RuntimeError(f"{HSP_SCRIPT} 失败，请查看上方输出")
    HSP = _load_hist_std_cache(HSP_CACHE, np.asarray(BG_GRID, float),
                               N_MC_HSP, N_SHOTS_LIST)
    if HSP is None:
        raise RuntimeError(f"扫描结束仍无法载入 {HSP_CACHE}")
print(f"模块 11 数据就绪：{len(BG_GRID)} bg × N={N_SHOTS_LIST} × {N_MC_HSP:,} MC")

_bgg = np.asarray(BG_GRID, float)


def _draw_m11(xmode, fname, xlabel, suptitle):
    """xmode='bg' 横轴用 bg；'noise' 横轴用 bg/N（=单发底噪）。"""
    fig, ax = plt.subplots(1, 3, figsize=(17.0, 5.0))

    # ① 单条 hist 内的 std
    for n in N_SHOTS_LIST:
        x = _bgg if xmode == "bg" else _bgg / n
        ax[0].plot(x, HSP[n]["hist_std"], "-", color=_COLORS_N[n], lw=2.0,
                   label=f"N={n}　MC 实测")
        ax[0].plot(x, np.sqrt(_bgg * (1.0 - _bgg / (N_PIX_MACRO * n))), "--",
                   color=_COLORS_N[n], lw=1.1, alpha=0.8,
                   label=f"N={n}　解析 √(bg(1−bg/{N_PIX_MACRO*n}))")
    if xmode == "bg":
        ax[0].plot(_bgg, np.sqrt(_bgg), "k:", lw=1.5, label="纯泊松 √bg（无饱和）")
    ax[0].set_xlabel(xlabel)
    ax[0].set_ylabel("单条 hist 内 152 个 bin 的 std（对 MC 取平均）")
    ax[0].set_title("① 每条直方图自身的起伏（亚泊松）", fontsize=11)
    ax[0].legend(fontsize=7.4); ax[0].grid(alpha=0.3)

    # ② peak 均值
    for n in N_SHOTS_LIST:
        x = _bgg if xmode == "bg" else _bgg / n
        ax[1].plot(x, HSP[n]["peak_mean"], "-", color=_COLORS_N[n], lw=2.0,
                   label=f"N={n} peak 均值")
        ax[1].plot(x, HSP[n]["bg_mc"], ":", color=_COLORS_N[n], lw=1.1, alpha=0.7)
    ax[1].set_xlabel(xlabel); ax[1].set_ylabel("peak 均值 [计数/bin]")
    ax[1].set_title("② peak 均值（点线 = 实测 bg，作参照）", fontsize=11)
    ax[1].legend(fontsize=8); ax[1].grid(alpha=0.3)

    # ③ peak 标准差：实测 vs Gumbel 极值解析
    #    peak = max of ~M_eff 个近高斯 bin，Gumbel 极限下
    #    σ_peak ≈ (π/√6)·σ_bin/z_M，其中 z_M=√(2 ln M_eff)≈(peakμ−bg)/σ_bin。
    #    用实测 peakμ 反出 z_M（不含拟合参数），是「给定峰位、峰宽是否为极值宽度」的检验。
    for n in N_SHOTS_LIST:
        x = _bgg if xmode == "bg" else _bgg / n
        ax[2].plot(x, HSP[n]["peak_std"], "-", color=_COLORS_N[n], lw=2.0,
                   label=f"N={n}　MC 实测")
        sig_bin = np.sqrt(_bgg * (1.0 - _bgg / (N_PIX_MACRO * n)))
        z_M = (HSP[n]["peak_mean"] - _bgg) / np.maximum(sig_bin, 1e-9)
        sig_evt = (np.pi / np.sqrt(6.0)) * sig_bin / np.maximum(z_M, 1e-9)
        ax[2].plot(x, sig_evt, "--", color=_COLORS_N[n], lw=1.2, alpha=0.85,
                   label=f"N={n}　Gumbel 解析")
    ax[2].set_xlabel(xlabel); ax[2].set_ylabel("peak 标准差 [计数/bin]")
    ax[2].set_title("③ peak 标准差：实测 vs Gumbel 极值解析", fontsize=11)
    ax[2].legend(fontsize=7.0, ncol=2); ax[2].grid(alpha=0.3)

    fig.suptitle(suptitle, fontsize=12.5)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(fname, dpi=120, bbox_inches="tight")
    plt.show()


_draw_m11("bg", "pod_v20_m11_vs_bg.png", "bg（hist_add 统计窗每 bin 均值）",
          f"模块 11　三个统计量 vs bg —— 同底噪比 N（noise=bg/N，每档 {N_MC_HSP:,} MC）")
_draw_m11("noise", "pod_v20_m11_vs_noise.png", "noise（单发 hist_i 统计窗均值）",
          f"模块 11　三个统计量 vs noise —— 同单发底噪比 N（bg=N·noise，每档 {N_MC_HSP:,} MC）")

print("=" * 118)
print("histσ=单条 hist 内 std；解析σ=√(bg(1−bg/27N))；peakσ=实测 peak std；σ_EVT=Gumbel 极值解析")
print(f"{'bg':>6} | " + " | ".join(
    f"N={n}: histσ  解析σ  peakμ  peakσ  σ_EVT" for n in N_SHOTS_LIST))
for k in range(0, len(_bgg), 6):
    cells = []
    for n in N_SHOTS_LIST:
        bgv = _bgg[k]
        ana = np.sqrt(bgv * (1.0 - bgv / (N_PIX_MACRO * n)))
        z_M = (HSP[n]['peak_mean'][k] - bgv) / max(ana, 1e-9)
        s_evt = (np.pi / np.sqrt(6.0)) * ana / max(z_M, 1e-9)
        cells.append(f"{HSP[n]['hist_std'][k]:6.3f} {ana:6.3f} "
                     f"{HSP[n]['peak_mean'][k]:6.2f} {HSP[n]['peak_std'][k]:6.3f} "
                     f"{s_evt:6.3f}")
    print(f"{_bgg[k]:6.2f} | " + " | ".join(cells))

print("\n【模块 11 读图要点】")
print("  ① 实测 histσ 应贴合 √(bg(1−bg/27N))，明显低于纯泊松 √bg —— 这是 1 bit SPAD 的亚泊松压缩；")
print("     N 越大 n_tr=27N 越大，压缩越弱，三条线在小 bg 处几乎重合、大 bg 处分开。")
print("  ② peak 均值远高于 bg（152 个 bin 取极大值），且 peak−bg 随 bg 增长后趋于平缓。")
print("  ③ peak 标准差对比 Gumbel 极值解析 σ_EVT=(π/√6)·σ_bin/z_M（z_M=(peakμ−bg)/σ_bin）：")
print("     大 bg（bg≳6）两者贴合到几个百分点；小 bg 时 EVT 偏低，因为计数少、分布离散且强右偏，")
print("     还没进入极值定理成立的渐近区。这条曲线把 ①（σ_bin）②（peakμ）③（peakσ）三者串起来。")
'''

M12_MD = r"""
## 模块 12 — ★ v20：阈值曲线用折线（点与点直线相连），并给连续阈值

### 需求澄清

"不要画阶梯"指的是**画法**：把每个 bg 档算出的阈值当成一个数据点，
**点与点之间用直线连成折线**，而不是用阶梯函数（`step`，带竖直立边）去连。
本模块图 A 的实线+点就是这个折线。模块 6 的 noise–threshold 曲线本来就是折线；
这里补一张以 bg 为横轴、6 条 FAR 都画全的折线版，并额外叠一条更平滑的连续阈值。

> 注意：阈值本身仍是**整数**（硬件只能比较整数计数），折线只是把这些整数点连起来看趋势。
> 下面的"连续阈值 $T_c$"是另一件事 —— 它把整数量化也抹掉，方便比较和拟合。

### 连续阈值的定义

生存函数 $S(T)=P(\mathrm{peak}\ge T)$ 只在整数上有定义。在 $\log S$ 上做线性插值，
定义**连续阈值** $T_c$ 为使插值后的生存函数恰好等于目标 FAR 的那个实数：

$$
T_c = j-1 + \frac{\ln S(j\!-\!1) - \ln \mathrm{FAR}}{\ln S(j\!-\!1) - \ln S(j)},
\qquad j=\min\{T: S(T)<\mathrm{FAR}\}
$$

* $T_c$ 与整数阈值的关系是 $T_{\text{整数}} = \lceil T_c \rceil$，两者**不矛盾**；
* 用 $\log S$ 而不是 $S$ 插值，是因为尾部近似指数衰减，$\log S$ 才接近线性；
* **硬件仍然只能用整数阈值**。$T_c$ 是用来看趋势、做比较、做拟合的量。

### 分辨极限

MC 条数 $M=10^6$ 时，生存函数最小只能分辨到 $1/M=10^{-6}$。
若某档在目标 FAR 处已经落到 $S=0$，$T_c$ 记为 `NaN`（图上断开），不做外推。
10 ppm 这一档在小 bg 时最容易碰到这个极限。
"""

M12_CODE = r'''
# ---- 模块 12：连续阈值 ----
def far_threshold_continuous(cnt, target_far):
    """连续（实数）阈值：在 log(生存函数) 上线性插值，使 S(T_c)=target_far。

    返回 NaN 表示 MC 分辨不到（目标 FAR 处生存函数已经为 0）。
    整数阈值 = ceil(T_c)，与 far_threshold_from_cnt 一致。
    """
    cnt = np.asarray(cnt)
    n = int(cnt.sum())
    if n <= 0:
        return np.nan
    n_ge = np.concatenate([[n], n - np.cumsum(cnt)[:-1]]).astype(float)
    sf = n_ge / n
    idx = np.where(sf < target_far)[0]
    if idx.size == 0:
        return np.nan                      # 整条曲线都在目标之上，MC 不够
    j = int(idx[0])
    if j == 0:
        return 0.0
    s_hi, s_lo = sf[j - 1], sf[j]
    if s_lo <= 0.0 or s_hi <= 0.0:
        return np.nan                      # 掉进 MC 分辨极限，不外推
    if s_hi <= target_far:
        return float(j - 1)
    w = (np.log(s_hi) - np.log(target_far)) / (np.log(s_hi) - np.log(s_lo))
    return float((j - 1) + w)


THRESH_C = {}
for n in N_SHOTS_LIST:
    R = NOISE_RES[n]
    rec = {"bg": np.asarray(R["noise_mc"], float)}
    st = [peak_stats_from_cnt(c) for c in R["peak_cnt"]]
    rec["peak_mean"] = np.array([s["mean"] for s in st])
    rec["peak_std"] = np.array([s["std"] for s in st])
    for far in TARGET_FARS:
        tag = FAR_TAG[far]
        rec["Tc" + tag] = np.array(
            [far_threshold_continuous(c, far) for c in R["peak_cnt"]])
    THRESH_C[n] = rec

_nan_report = []
for n in N_SHOTS_LIST:
    for far in TARGET_FARS:
        bad = int(np.sum(~np.isfinite(THRESH_C[n]["Tc" + FAR_TAG[far]])))
        if bad:
            _nan_report.append(f"N={n} {FAR_LABEL[far]}: {bad}/{len(BG_GRID)} 档超出 MC 分辨")
print("连续阈值计算完成。" + ("　".join(_nan_report) if _nan_report else "全部档位均在 MC 分辨范围内。"))

# --- 图 A：阈值 vs bg。整数阈值的采样点用【直线】连成折线（不是阶梯函数）；
#          连续阈值 Tc 作为更平滑的趋势线叠上去 ---
fig, axes = plt.subplots(1, len(N_SHOTS_LIST),
                         figsize=(5.7 * len(N_SHOTS_LIST), 5.0), sharex=True)
_far_cols = plt.cm.viridis(np.linspace(0.05, 0.92, len(TARGET_FARS)))
for a, n in zip(np.atleast_1d(axes), N_SHOTS_LIST):
    bg = THRESH_C[n]["bg"]
    for far, c in zip(TARGET_FARS, _far_cols):
        tag = FAR_TAG[far]
        # 硬件用的整数阈值：每个 bg 档一个点，点与点之间用直线连（折线）
        a.plot(bg, THRESH[n]["T" + tag], "-", marker=".", ms=3.5, color=c,
               lw=1.5, label=f"{FAR_LABEL[far]}")
        # 连续（实数）阈值：更平滑的趋势线，虚线叠加
        a.plot(bg, THRESH_C[n]["Tc" + tag], "--", color=c, lw=1.0, alpha=0.7)
    a.plot(bg, bg, ":", color="0.45", lw=1.3, label="参考 T=bg")
    a.axhline(NOISE_RES[n]["n_tr"], color="k", ls="-.", lw=1.1, alpha=0.7,
              label=f"二值硬上限 {NOISE_RES[n]['n_tr']}")
    a.set_xlabel("bg [计数/bin]"); a.set_ylabel("阈值 T [计数]")
    a.set_title(f"N_shots={n}（实线+点=整数阈值折线，虚线=连续 $T_c$）", fontsize=10.5)
    a.legend(fontsize=7.6, ncol=2); a.grid(alpha=0.3)
fig.suptitle("模块 12A　阈值 vs bg —— 采样点用直线相连的折线（6 条 FAR）", fontsize=12.5)
fig.tight_layout(rect=[0, 0, 1, 0.93])
fig.savefig("pod_v20_m12_Tc_vs_bg.png", dpi=120, bbox_inches="tight")
plt.show()

# --- 图 B：把三个 N 叠在同一张图上比（每条 FAR 一个 panel）---
_FAR_M12 = [0.05, 0.01, 100e-6]
fig, axes = plt.subplots(2, len(_FAR_M12), figsize=(5.6 * len(_FAR_M12), 8.4))
for j, far in enumerate(_FAR_M12):
    tag = FAR_TAG[far]
    for n in N_SHOTS_LIST:
        bg = THRESH_C[n]["bg"]
        axes[0, j].plot(bg, THRESH_C[n]["Tc" + tag], "-", color=_COLORS_N[n],
                        lw=2.0, label=f"N={n}")
        axes[1, j].plot(bg, THRESH_C[n]["Tc" + tag] - bg, "-", color=_COLORS_N[n],
                        lw=2.0, label=f"N={n}")
    axes[0, j].set_title(f"$T_c$ vs bg　FAR={FAR_LABEL[far]}", fontsize=11)
    axes[0, j].set_ylabel("$T_c$ [计数]")
    axes[1, j].set_title(f"$T_c$ − bg　FAR={FAR_LABEL[far]}", fontsize=11)
    axes[1, j].set_ylabel("$T_c$ − bg [计数]")
    for r in (0, 1):
        axes[r, j].set_xlabel("bg [计数/bin]")
        axes[r, j].grid(alpha=0.3); axes[r, j].legend(fontsize=8)
fig.suptitle("模块 12B　同 bg 下 N=1/2/4 的连续阈值与阈值余量", fontsize=12.5)
fig.tight_layout(rect=[0, 0, 1, 0.94])
fig.savefig("pod_v20_m12_Tc_compare.png", dpi=120, bbox_inches="tight")
plt.show()

# --- 图 C：有效 z 值 (Tc − peak_mean)/peak_std ---
fig, axes = plt.subplots(1, len(_FAR_M12), figsize=(5.4 * len(_FAR_M12), 4.4),
                         sharey=True)
for a, far in zip(np.atleast_1d(axes), _FAR_M12):
    tag = FAR_TAG[far]
    for n in N_SHOTS_LIST:
        rec = THRESH_C[n]
        z = (rec["Tc" + tag] - rec["peak_mean"]) / np.maximum(rec["peak_std"], 1e-9)
        a.plot(rec["bg"], z, "-", color=_COLORS_N[n], lw=1.9, label=f"N={n}")
    a.set_xlabel("bg"); a.set_ylabel(r"$z=(T_c-\mu_{peak})/\sigma_{peak}$")
    a.set_title(f"FAR={FAR_LABEL[far]}", fontsize=11)
    a.grid(alpha=0.3); a.legend(fontsize=8)
fig.suptitle("模块 12C　有效 z 值：阈值离 peak 均值几个 peak 标准差", fontsize=12.5)
fig.tight_layout(rect=[0, 0, 1, 0.92])
fig.savefig("pod_v20_m12_zeff.png", dpi=120, bbox_inches="tight")
plt.show()

print("=" * 112)
print("连续阈值抽样（括号内为硬件用的整数阈值 = ceil(Tc)）")
_hdr = f"{'bg':>6}"
for far in _FAR_M12:
    for n in N_SHOTS_LIST:
        _hdr += f"{'N'+str(n)+'@'+FAR_LABEL[far]:>16}"
print(_hdr)
for k in range(0, len(BG_GRID), 6):
    row = f"{BG_GRID[k]:6.2f}"
    for far in _FAR_M12:
        tag = FAR_TAG[far]
        for n in N_SHOTS_LIST:
            tc = THRESH_C[n]["Tc" + tag][k]
            ti = THRESH[n]["T" + tag][k]
            row += f"{tc:11.2f}({ti:3d})" if np.isfinite(tc) else f"{'nan':>11}({ti:3d})"
    print(row)
'''

M13_MD = r"""
## 模块 13 — ★ v20：FAR = 5% / 1% / 100 ppm 下，PoD50 与 PoD90 所需信号的均值

对这三条阈值曲线，逐 bg、逐 N_shots 给出**刚好达到 50% / 90% 探测概率所需的信号**。
"所需信号"给四种口径，避免歧义：

| 口径 | 定义 | 说明 |
|---|---|---|
| `peak_mean` | 临界点上 `hist_add` 峰值的均值 | 含底噪，直接从 MC 统计 |
| **`S_net = peak_mean − bg`** | 净峰高 | 扣掉平均底噪后的峰高 |
| **`ΔS = peak_mean − peak_mean(boost=0)`** | 相对**同窗口无信号**基线的增量 | 最干净的"信号带来了多少" |
| `E_crit` | 临界发射脉冲能量 [nJ] | 可直接与激光器指标对照 |

> `S_net` 与 `ΔS` 不相等，且 `ΔS < S_net`。原因：即使没有信号，
> 在 15 个 bin 的信号窗里取最大值本身也已经高于 bg（极值统计的抬升 ≈ 2σ）。
> 所以 `S_net` 会**高估**信号贡献，`ΔS` 才是"信号净增量"。
> `ΔS` 需要模块 9.3 / 15 的 `boost=0` 基线（同一子窗、同一 bg、同一 N），本模块自动引用。

数据来源：模块 7 的 `POD_RES`（每个临界点都做过 `N_MC_POD_VERIFY` 次独立验证）。
"""

M13_CODE = r'''
# ---- 模块 13：三条 FAR × PoD50/90 的所需信号 ----
FAR_M13 = [0.05, 0.01, 100e-6]           # 5% / 1% / 100 ppm
LEVELS_M13 = [0.50, 0.90]
_LS_LEVEL = {0.50: "--", 0.90: "-"}


def _sig_baseline_peak(n, k):
    """同一信号子窗、同一 bg、同一 N，boost=0 时的 peak 均值（无信号基线）。"""
    try:
        return float(SIG_M9[n]["peak_mean"][0, k])
    except Exception:
        return np.nan


def collect_pod_signal(n_shots, far_tag, level):
    """返回按 bg 排好的 (bg, boost, pod, peak_mean, peak_std, T, S_net, dS, E_nJ, Nph)。"""
    rows = []
    for k, nt in enumerate(NOISE_GRID[n_shots]):
        r = POD_RES.get((n_shots, float(nt)))
        if not r or "critical" not in r:
            continue
        rec = r["critical"].get(far_tag, {}).get(f"{level:.2f}")
        if not rec:
            continue
        bg = float(r["noise"])
        base = _sig_baseline_peak(n_shots, k)
        rows.append((
            bg, rec["boost"], rec["pod"], rec["peak_mean"], rec["peak_std"],
            r["T_map"][far_tag], rec["peak_mean"] - bg,
            rec["peak_mean"] - base,
            rec["boost"] * E_PULSE_BASE * 1e9,
            rec["boost"] * _NPH_BASE,
        ))
    a = np.asarray(rows, float)
    return a[np.argsort(a[:, 0])] if a.size else a


SIGREQ = {(n, FAR_TAG[f], lv): collect_pod_signal(n, FAR_TAG[f], lv)
          for n in N_SHOTS_LIST for f in FAR_M13 for lv in LEVELS_M13}
_ncov = {k: (0 if v.size == 0 else v.shape[0]) for k, v in SIGREQ.items()}
print(f"模块 13 覆盖度：每条曲线取到 {min(_ncov.values())}–{max(_ncov.values())} / "
      f"{len(BG_GRID)} 个 bg 档")
if not np.isfinite([_sig_baseline_peak(N_SHOTS_LIST[0], 0)]).all():
    print("  注意：未取到 boost=0 基线（模块 9.3 缓存缺失），ΔS 一列将为 NaN")

# --- 图 A：净峰高 S_net 与 相对无信号基线的增量 ΔS ---
fig, axes = plt.subplots(2, len(FAR_M13), figsize=(5.8 * len(FAR_M13), 8.8))
for j, far in enumerate(FAR_M13):
    tag = FAR_TAG[far]
    for n in N_SHOTS_LIST:
        for lv in LEVELS_M13:
            a = SIGREQ[(n, tag, lv)]
            if not a.size:
                continue
            axes[0, j].plot(a[:, 0], a[:, 6], _LS_LEVEL[lv], color=_COLORS_N[n],
                            lw=1.9, label=f"N={n} PoD{int(lv*100)}")
            axes[1, j].plot(a[:, 0], a[:, 7], _LS_LEVEL[lv], color=_COLORS_N[n],
                            lw=1.9, label=f"N={n} PoD{int(lv*100)}")
    axes[0, j].set_title(f"净峰高 $S_{{net}}$ = peak−bg　FAR={FAR_LABEL[far]}", fontsize=11)
    axes[0, j].set_ylabel("$S_{net}$ [计数]")
    axes[1, j].set_title(f"净增量 $\\Delta S$（相对无信号基线）　FAR={FAR_LABEL[far]}",
                         fontsize=11)
    axes[1, j].set_ylabel("$\\Delta S$ [计数]")
    for r in (0, 1):
        axes[r, j].set_xlabel("bg [计数/bin]")
        axes[r, j].grid(alpha=0.3); axes[r, j].legend(fontsize=7.4, ncol=2)
fig.suptitle("模块 13A　达到 PoD50 / PoD90 所需的信号峰高（实线 PoD90，虚线 PoD50）",
             fontsize=12.5)
fig.tight_layout(rect=[0, 0, 1, 0.94])
fig.savefig("pod_v20_m13_signal_counts.png", dpi=120, bbox_inches="tight")
plt.show()

# --- 图 B：临界发射能量与等效信号光子数 ---
fig, axes = plt.subplots(2, len(FAR_M13), figsize=(5.8 * len(FAR_M13), 8.8))
for j, far in enumerate(FAR_M13):
    tag = FAR_TAG[far]
    for n in N_SHOTS_LIST:
        for lv in LEVELS_M13:
            a = SIGREQ[(n, tag, lv)]
            if not a.size:
                continue
            axes[0, j].semilogy(a[:, 0], a[:, 8], _LS_LEVEL[lv], color=_COLORS_N[n],
                                lw=1.9, label=f"N={n} PoD{int(lv*100)}")
            axes[1, j].semilogy(a[:, 0], a[:, 9], _LS_LEVEL[lv], color=_COLORS_N[n],
                                lw=1.9, label=f"N={n} PoD{int(lv*100)}")
    axes[0, j].set_title(f"临界发射能量　FAR={FAR_LABEL[far]}", fontsize=11)
    axes[0, j].set_ylabel("$E_{crit}$ [nJ]")
    axes[1, j].set_title(f"等效信号光子数　FAR={FAR_LABEL[far]}", fontsize=11)
    axes[1, j].set_ylabel("到达宏像元的信号光子数")
    for r in (0, 1):
        axes[r, j].set_xlabel("bg [计数/bin]")
        axes[r, j].grid(alpha=0.3, which="both"); axes[r, j].legend(fontsize=7.4, ncol=2)
fig.suptitle("模块 13B　达到 PoD50 / PoD90 所需的发射能量与信号光子数", fontsize=12.5)
fig.tight_layout(rect=[0, 0, 1, 0.94])
fig.savefig("pod_v20_m13_signal_energy.png", dpi=120, bbox_inches="tight")
plt.show()

print("=" * 126)
print("模块 13　抽样表（FAR=1%）")
print(f"{'N':>3}{'bg':>7}{'PoD':>6}{'T':>5}{'peakμ':>9}{'peakσ':>8}"
      f"{'S_net':>9}{'ΔS':>9}{'E[nJ]':>11}{'Nph':>10}{'验证PoD':>9}")
for n in N_SHOTS_LIST:
    a_by_lv = {lv: SIGREQ[(n, FAR_TAG[0.01], lv)] for lv in LEVELS_M13}
    for lv in LEVELS_M13:
        a = a_by_lv[lv]
        if not a.size:
            continue
        for i in range(0, a.shape[0], max(1, a.shape[0] // 8)):
            print(f"{n:>3d}{a[i,0]:>7.2f}{int(lv*100):>6d}{a[i,5]:>5.0f}"
                  f"{a[i,3]:>9.2f}{a[i,4]:>8.2f}{a[i,6]:>9.2f}{a[i,7]:>9.2f}"
                  f"{a[i,8]:>11.4g}{a[i,9]:>10.4g}{a[i,2]:>9.3f}")
'''

M14_MD = r"""
## 模块 14 — ★ v20：平方反比下的测远估计

### 两种口径

把临界能量倍数 `boost` 折算成距离时，本仓库一直用的是**含大气衰减**的版本
（`equiv_distance()`，模块 8 已有）：

$$\mathrm{boost}(D)=\frac{D_{ref}^2}{D^2}\,e^{-2\alpha (D-D_{ref})},\qquad \alpha=0.1\ \mathrm{km^{-1}}$$

本模块**补上用户要求的纯平方反比口径**，即忽略大气吸收、只保留 $1/D^2$：

$$\mathrm{boost}=\frac{D_{ref}^2}{D^2}\ \Longrightarrow\ \boxed{D_{1/r^2}=\frac{D_{ref}}{\sqrt{\mathrm{boost}}}}$$

其中 $D_{ref}=15$ m 是仿真里的基准目标距离，`boost=1` 对应基准回波。

两者的关系：纯平方反比是**上界**，含衰减的结果一定更近。
在 $\alpha=0.1$ /km 下，100 m 处衰减因子 $e^{-2\alpha(100-15)}=e^{-0.017}\approx 0.983$，
只差 1.7%；到 500 m 才差到 $e^{-0.097}\approx 0.91$。
也就是说在本仓库关心的距离段，**大气衰减不是主导，平方反比就是主要规律**。

### 前提条件（很重要）

$1/D^2$ 成立的前提是**目标是充满视场的朗伯扩展面**（回波功率 $\propto \rho A_{RX}/D^2$）。
若目标比光斑小（点目标 / 远处小物体），实际会更陡。本模块的结论只在扩展面假设下成立，
与 `equiv_distance()` 的假设一致。

### 曲线上的抖动是哪来的

测距曲线不是光滑的，逐点会有几米的上下跳动，有两个来源，**都不是物理效应**：

1. **整数阈值的阶梯**。$T$ 只能取整数，bg 连续变化时 $T$ 一格一格跳，
   临界能量跟着跳（模块 12 的连续阈值就是为了看清这一点）。
2. **PoD 临界点求解的 MC 噪声**。每个临界点用 `N_MC_POD_VERIFY=5000` 次独立验证，
   收敛容差 `POD_VERIFY_TOL=0.02`，所以 PoD 落在 0.88–0.92 之间就算收敛，
   对应的能量还有几个百分点的不确定度；$D\propto 1/\sqrt{E}$ 又把它折一半。

**看趋势，不要抠单点。** 需要更平滑的曲线就加大 `N_MC_POD_VERIFY` 并收紧 `POD_VERIFY_TOL`。
"""

M14_CODE = r'''
# ---- 模块 14：平方反比测远 ----
def range_inverse_square(boost, D_ref=D_TARGET):
    """纯平方反比：boost = (D_ref/D)^2 → D = D_ref / sqrt(boost)。"""
    b = np.asarray(boost, float)
    out = np.full(b.shape, np.nan)
    m = np.isfinite(b) & (b > 0)
    out[m] = D_ref / np.sqrt(b[m])
    return out if out.ndim else float(out)


RANGE_M14 = {}
for n in N_SHOTS_LIST:
    for far in FAR_M13:
        tag = FAR_TAG[far]
        for lv in LEVELS_M13:
            a = SIGREQ[(n, tag, lv)]
            if not a.size:
                RANGE_M14[(n, tag, lv)] = a
                continue
            d_sq = range_inverse_square(a[:, 1])
            d_at = np.array([equiv_distance(b) for b in a[:, 1]])
            RANGE_M14[(n, tag, lv)] = np.column_stack([a[:, 0], d_sq, d_at, a[:, 1]])

# --- 图 A：测距 vs bg（实线 = 纯 1/D²，虚线 = 含大气衰减）---
fig, axes = plt.subplots(len(LEVELS_M13), len(FAR_M13),
                         figsize=(5.8 * len(FAR_M13), 4.6 * len(LEVELS_M13)))
axes = np.atleast_2d(axes)
for i, lv in enumerate(LEVELS_M13):
    for j, far in enumerate(FAR_M13):
        tag = FAR_TAG[far]
        a_ = axes[i, j]
        for n in N_SHOTS_LIST:
            d = RANGE_M14[(n, tag, lv)]
            if not d.size:
                continue
            a_.plot(d[:, 0], d[:, 1], "-", color=_COLORS_N[n], lw=2.0,
                    label=f"N={n} 纯 1/D²")
            a_.plot(d[:, 0], d[:, 2], "--", color=_COLORS_N[n], lw=1.2, alpha=0.8,
                    label=f"N={n} 含大气衰减")
        a_.set_xlabel("bg [计数/bin]"); a_.set_ylabel("等效探测距离 [m]")
        a_.set_title(f"PoD{int(lv*100)}　FAR={FAR_LABEL[far]}", fontsize=11)
        a_.grid(alpha=0.3); a_.legend(fontsize=7.2, ncol=2)
fig.suptitle(f"模块 14　平方反比测远（基准 D_ref={D_TARGET:g} m，反射率 {RHO_TARGET:.0%}，"
             f"发射能量固定为临界值）", fontsize=12.5)
fig.tight_layout(rect=[0, 0, 1, 0.94])
fig.savefig("pod_v20_m14_range_vs_bg.png", dpi=120, bbox_inches="tight")
plt.show()

# --- 图 B：两种口径的相对差 ---
fig, ax = plt.subplots(figsize=(7.6, 4.8))
for n in N_SHOTS_LIST:
    d = RANGE_M14[(n, FAR_TAG[0.01], 0.90)]
    if not d.size:
        continue
    ax.plot(d[:, 1], 100 * (d[:, 1] - d[:, 2]) / np.maximum(d[:, 1], 1e-9),
            "-o", ms=3, color=_COLORS_N[n], lw=1.7, label=f"N={n}")
ax.set_xlabel("纯 1/D² 测距 [m]")
ax.set_ylabel("大气衰减让测距缩短的百分比 [%]")
ax.set_title("模块 14B　大气衰减相对平方反比的修正量（FAR=1%, PoD90）", fontsize=11.5)
ax.grid(alpha=0.3); ax.legend(fontsize=9)
fig.tight_layout()
fig.savefig("pod_v20_m14_atten_correction.png", dpi=120, bbox_inches="tight")
plt.show()

print("=" * 108)
print("模块 14　平方反比测距抽样（PoD90）")
print(f"{'N':>3}{'bg':>7}" + "".join(
    f"{FAR_LABEL[f]+' 1/D²':>14}{FAR_LABEL[f]+' +衰减':>14}" for f in FAR_M13))
for n in N_SHOTS_LIST:
    ref = RANGE_M14[(n, FAR_TAG[0.01], 0.90)]
    if not ref.size:
        continue
    for i in range(0, ref.shape[0], max(1, ref.shape[0] // 8)):
        bg = ref[i, 0]
        row = f"{n:>3d}{bg:>7.2f}"
        for far in FAR_M13:
            d = RANGE_M14[(n, FAR_TAG[far], 0.90)]
            if not d.size:
                row += f"{'-':>14}{'-':>14}"
                continue
            k = int(np.argmin(np.abs(d[:, 0] - bg)))
            row += f"{d[k,1]:>14.1f}{d[k,2]:>14.1f}"
        print(row)
'''

M15_MD = r"""
## 模块 15 — ★ v20：同一信号强度、不同 bg 时 peak 分布怎么变

这一模块回答三个具体问题（数据来自模块 9.3 的完整 `peak_cnt` 分布，
每个 `(N, bg, boost)` 都有 8,000 条 MC 的**完整直方图**，不只是均值方差）：

1. **分布形状怎么变**：把同一个 `boost` 在不同 bg 下的 peak 概率质量函数（PMF）叠在一起看。
2. **peak 均值是不是"只是加上 bg"**：定义
   $$\Delta\mu(b,\mathrm{bg}) = \mu_{peak}(b,\mathrm{bg}) - \mu_{peak}(0,\mathrm{bg})$$
   如果"加信号 = 在原来的基础上平移"，那么 $\Delta\mu$ 应该**与 bg 无关**（水平线）。
3. **peak 标准差怎么变**：对比 $\sigma_{peak}(b,\mathrm{bg})$ 与 $\sigma_{peak}(0,\mathrm{bg})$。

另外给出**偏度**（三阶中心矩 / $\sigma^3$）随 bg 的变化，用来判断形状是不是真的只是平移。

### "抢占"模型：为什么 $\Delta\mu$ 会随 bg 掉下去

SPAD 是 1 bit 的：**一个 bin 在一发里只能亮一次**。如果环境光已经把这个 SPAD 点亮了，
信号光子再来也不产生额外计数。设单个 SPAD 单发被环境光点亮的概率是 $p_{eq}$、
被信号点亮的概率是 $q$，两者独立竞争，则该 bin 亮的概率是
$1-(1-p_{eq})(1-q)$，**信号带来的净增量**是

$$
\Delta p = \big[1-(1-p_{eq})(1-q)\big]-p_{eq} = q\,(1-p_{eq})
$$

于是有一个不含拟合参数的预言：

$$
\boxed{\ \frac{\Delta\mu(\mathrm{bg})}{\Delta\mu(\mathrm{bg}_{\min})}\approx\frac{1-p_{eq}(\mathrm{bg})}{1-p_{eq}(\mathrm{bg}_{\min})}\ }
$$

图 15B 下排会把这条**黑色虚线**画上去与实测对照。
两者的差距就是"抢占"之外的第二个机制：**极值竞争** ——
bg 高时无信号基线 $\mu(0)$ 本身已被 15 个 bin 取极大值抬高，信号 bin 未必总是全窗最大值，
所以实测会掉得比"抢占"模型更快。

### 实测结论（8,000 MC/档，`boost=0.016`）

| 现象 | 数值 |
|---|---|
| $\Delta\mu$ **不是常数** | N=1 时从 bg=0.25 的 **6.77** 掉到 bg=10.25 的 **2.69**（掉 60%）；N=4 从 27.8 掉到 22.5 |
| $\sigma_{peak}$ 加信号后**变大而不是变小** | $\sigma(b)/\sigma(0)$ 从低 bg 的 2.8–5.3 降到高 bg 的 1.04–1.64，但**始终 ≥ 1** |
| 形状确实在变 | 纯噪声 peak 右偏（偏度 0.2–1.0，极值分布特征）；加信号后偏度掉到 **≈0.05**，接近对称 |

所以对"**peak 的均值是否简单地只是加上 bg**"这个问题，答案是：

> **不是。** 低 bg 时近似成立，但 bg 越高，同样的信号能顶上去的高度越少。
> 主因是 1 bit SPAD 的"抢占"（环境光先点亮就轮不到信号），次因是极值竞争。
> 而且分布形状从右偏的极值分布变成近似对称的二项分布，
> **用"均值 + 标准差"两个数已经不足以描述它**，必须看完整分布（图 15A）。
"""

M15_CODE = r'''
# ---- 模块 15：同信号强度、不同 bg 下 peak 分布的变化 ----
def _moments_from_cnt(cnt):
    """由 bincount 求 (均值, 标准差, 偏度)。"""
    cnt = np.asarray(cnt, float)
    n = cnt.sum()
    if n <= 0:
        return np.nan, np.nan, np.nan
    v = np.arange(cnt.size, dtype=float)
    m1 = (v * cnt).sum() / n
    c2 = (((v - m1) ** 2) * cnt).sum() / n
    c3 = (((v - m1) ** 3) * cnt).sum() / n
    sd = np.sqrt(max(c2, 0.0))
    sk = c3 / sd ** 3 if sd > 1e-12 else np.nan
    return float(m1), float(sd), float(sk)


_BOOSTS = np.asarray(BOOST_LIST_M9, float)
_BG_SHOW = [0.5, 3.0, 6.0, 9.0, 12.0]                # 用于分布叠加的 bg 档
_BOOST_SHOW = float(_BOOSTS[len(_BOOSTS) // 2])      # 中间那档信号
_ib_show = int(np.argmin(np.abs(_BOOSTS - _BOOST_SHOW)))
_k_show = [int(np.argmin(np.abs(np.asarray(BG_GRID, float) - b))) for b in _BG_SHOW]

M15 = {}
for n in N_SHOTS_LIST:
    mu = np.zeros((len(_BOOSTS), len(BG_GRID)))
    sd = np.zeros_like(mu); sk = np.zeros_like(mu)
    for ib in range(len(_BOOSTS)):
        for k in range(len(BG_GRID)):
            mu[ib, k], sd[ib, k], sk[ib, k] = _moments_from_cnt(
                SIG_M9[n]["peak_cnt"][ib, k])
    M15[n] = dict(mu=mu, sd=sd, sk=sk)
print(f"模块 15：{len(_BOOSTS)} boost × {len(BG_GRID)} bg × N={N_SHOTS_LIST} 的完整分布已就绪；"
      f"分布叠加用 boost={_BOOST_SHOW:g}")

# --- 图 A：固定 boost，不同 bg 的 peak PMF 叠加 ---
fig, axes = plt.subplots(2, len(N_SHOTS_LIST),
                         figsize=(5.6 * len(N_SHOTS_LIST), 8.2))
axes = np.atleast_2d(axes)
_cbg = plt.cm.plasma(np.linspace(0.05, 0.85, len(_k_show)))
for j, n in enumerate(N_SHOTS_LIST):
    for c, k in zip(_cbg, _k_show):
        cnt = SIG_M9[n]["peak_cnt"][_ib_show, k].astype(float)
        tot = max(cnt.sum(), 1.0)
        x = np.arange(cnt.size)
        m = cnt > 0
        axes[0, j].plot(x[m], (cnt / tot)[m], "-", color=c, lw=1.6,
                        label=f"bg={BG_GRID[k]:g}")
        mu_, sd_, _ = _moments_from_cnt(cnt)
        axes[1, j].plot(((x - mu_) / max(sd_, 1e-9))[m], (cnt / tot)[m] * sd_, "-",
                        color=c, lw=1.6, label=f"bg={BG_GRID[k]:g}")
    axes[0, j].set_title(f"N={n}　boost={_BOOST_SHOW:g}　peak 分布", fontsize=11)
    axes[0, j].set_xlabel("peak [计数]"); axes[0, j].set_ylabel("概率")
    axes[1, j].set_title(f"N={n}　标准化后 $(peak-\\mu)/\\sigma$", fontsize=11)
    axes[1, j].set_xlabel(r"$(peak-\mu)/\sigma$"); axes[1, j].set_ylabel("概率密度")
    for r in (0, 1):
        axes[r, j].grid(alpha=0.3); axes[r, j].legend(fontsize=7.6)
fig.suptitle("模块 15A　同一信号强度下，peak 分布随 bg 的变化（下排标准化后看形状是否不变）",
             fontsize=12.5)
fig.tight_layout(rect=[0, 0, 1, 0.94])
fig.savefig("pod_v20_m15_peak_pmf.png", dpi=120, bbox_inches="tight")
plt.show()

# --- 图 B：Δμ 是否与 bg 无关（"是不是只是加上 bg"）---
fig, axes = plt.subplots(2, len(N_SHOTS_LIST),
                         figsize=(5.6 * len(N_SHOTS_LIST), 8.2))
axes = np.atleast_2d(axes)
_cb = plt.cm.viridis(np.linspace(0.05, 0.92, len(_BOOSTS) - 1))
for j, n in enumerate(N_SHOTS_LIST):
    base = M15[n]["mu"][0]
    for c, ib in zip(_cb, range(1, len(_BOOSTS))):
        d = M15[n]["mu"][ib] - base
        axes[0, j].plot(BG_GRID, d, "-", color=c, lw=1.6,
                        label=f"b={_BOOSTS[ib]:g}")
        axes[1, j].plot(BG_GRID, d / max(d[0], 1e-9), "-", color=c, lw=1.6)
    # 「抢占」模型：信号净增量 ∝ (1−p_eq)，不含任何拟合参数
    _peq = np.asarray(NOISE_RES[n]["p_eq"], float)
    axes[1, j].plot(BG_GRID, (1.0 - _peq) / (1.0 - _peq[0]), "k--", lw=1.8,
                    label=r"抢占模型 $(1-p_{eq})/(1-p_{eq}^{min})$")
    axes[0, j].set_title(f"N={n}　$\\Delta\\mu=\\mu(b)-\\mu(0)$", fontsize=11)
    axes[0, j].set_ylabel(r"$\Delta\mu$ [计数]")
    axes[1, j].set_title(f"N={n}　归一到最小 bg（水平=只是平移）", fontsize=11)
    axes[1, j].set_ylabel(r"$\Delta\mu(bg)\,/\,\Delta\mu(bg_{min})$")
    axes[1, j].axhline(1.0, color="0.4", ls=":", lw=1.1)
    for r in (0, 1):
        axes[r, j].set_xlabel("bg [计数/bin]")
        axes[r, j].grid(alpha=0.3)
    axes[0, j].legend(fontsize=7.0, ncol=2)
    axes[1, j].legend(fontsize=7.6)
fig.suptitle("模块 15B　peak 均值是不是「只是加上 bg」：$\\Delta\\mu$ 随 bg 的变化",
             fontsize=12.5)
fig.tight_layout(rect=[0, 0, 1, 0.94])
fig.savefig("pod_v20_m15_shift_test.png", dpi=120, bbox_inches="tight")
plt.show()

# --- 图 C：peak 标准差与偏度 ---
fig, axes = plt.subplots(2, len(N_SHOTS_LIST),
                         figsize=(5.6 * len(N_SHOTS_LIST), 8.2))
axes = np.atleast_2d(axes)
for j, n in enumerate(N_SHOTS_LIST):
    for c, ib in zip(plt.cm.viridis(np.linspace(0.05, 0.92, len(_BOOSTS))),
                     range(len(_BOOSTS))):
        lab = f"b={_BOOSTS[ib]:g}" + ("（无信号）" if ib == 0 else "")
        axes[0, j].plot(BG_GRID, M15[n]["sd"][ib], "-", color=c, lw=1.6, label=lab)
        axes[1, j].plot(BG_GRID, M15[n]["sk"][ib], "-", color=c, lw=1.6)
    axes[0, j].set_title(f"N={n}　peak 标准差", fontsize=11)
    axes[0, j].set_ylabel(r"$\sigma_{peak}$ [计数]")
    axes[1, j].set_title(f"N={n}　peak 偏度", fontsize=11)
    axes[1, j].set_ylabel("偏度")
    axes[1, j].axhline(0.0, color="k", ls=":", lw=1.1)
    for r in (0, 1):
        axes[r, j].set_xlabel("bg [计数/bin]")
        axes[r, j].grid(alpha=0.3)
    axes[0, j].legend(fontsize=7.0, ncol=2)
fig.suptitle("模块 15C　peak 标准差与偏度随 bg（不同信号强度）", fontsize=12.5)
fig.tight_layout(rect=[0, 0, 1, 0.94])
fig.savefig("pod_v20_m15_std_skew.png", dpi=120, bbox_inches="tight")
plt.show()

print("=" * 128)
print(f"模块 15　定量摘要（boost={_BOOST_SHOW:g}）　Δμ_rel=Δμ(bg)/Δμ(bg_min)，"
      f"抢占模型=(1−p_eq)/(1−p_eq_min)")
print(f"{'N':>3}{'bg':>7}{'μ(0)':>9}{'σ(0)':>8}{'偏度(0)':>9}"
      f"{'μ(b)':>9}{'σ(b)':>8}{'偏度(b)':>9}{'Δμ':>8}{'Δμ_rel':>9}"
      f"{'抢占模型':>10}{'σ(b)/σ(0)':>11}")
for n in N_SHOTS_LIST:
    d0 = M15[n]["mu"][_ib_show][0] - M15[n]["mu"][0][0]
    peq = np.asarray(NOISE_RES[n]["p_eq"], float)
    for k in range(0, len(BG_GRID), 8):
        m0, s0, k0 = M15[n]["mu"][0][k], M15[n]["sd"][0][k], M15[n]["sk"][0][k]
        mb, sb, kb = (M15[n]["mu"][_ib_show][k], M15[n]["sd"][_ib_show][k],
                      M15[n]["sk"][_ib_show][k])
        print(f"{n:>3d}{BG_GRID[k]:>7.2f}{m0:>9.2f}{s0:>8.2f}{k0:>9.2f}"
              f"{mb:>9.2f}{sb:>8.2f}{kb:>9.2f}{mb-m0:>8.2f}"
              f"{(mb-m0)/max(d0,1e-9):>9.3f}"
              f"{(1-peq[k])/(1-peq[0]):>10.3f}{sb/max(s0,1e-9):>11.3f}")

print("\n【模块 15 结论】")
print("  1) peak 均值【不是】简单地加上 bg。Δμ 随 bg 单调下降：同样的信号，bg 越高顶起来的高度越少。")
print("     主因是 1 bit SPAD 的「抢占」——环境光先点亮就轮不到信号，净增量 ∝ (1−p_eq)；")
print("     实测比抢占模型掉得更快，多出来的那部分来自极值竞争（高 bg 时无信号基线本身已被抬高）。")
print("  2) peak 标准差加信号后【变大】，不是变小；σ(b)/σ(0) 始终 ≥ 1，只是随 bg 升高趋近 1。")
print("     信号自己带来的二项涨落叠在噪声之上，并没有把 peak「钉死」。")
print("  3) 分布形状确实在变：纯噪声 peak 右偏（极值分布），加信号后偏度掉到 ≈0.05、接近对称。")
print("     所以只报「均值 + std」不足以描述 peak，需要看完整分布（图 15A）。")
'''

M16_MD = r"""
## 模块 16 — ★ v20：宏像元 3×9 vs 3×6 阈值对比

这一节把此前独立脚本 `compare_macro_3x9_vs_3x6.py` 的结论并入正本。
完整版（6 联图 + 全部数值表）见 `compare_macro_3x9_vs_3x6.png` 与 `compare_macro_3x9_vs_3x6_log.txt`，
这里画一个浓缩版并复述关键结论。

### 口径（不写清会得到相反的结论）

* **噪声在每个 SPAD 上均匀**：环境光只由 $p_{eq}$ 一个数刻画
  （单个 SPAD、单发、单个 1 ns bin 被点亮的**平衡态概率**），与宏像元多大、累加多少发都无关。
  于是宏像元每 bin 的底噪 $\mathrm{bg}=n_{tr}\cdot p_{eq}$，$n_{tr}=n_{pix}\times N_{shots}$。
* **信号也按每 SPAD 均匀**处理（用户明确要求；若按像斑加权，结论会反过来）。
  于是宏像元收到的信号 $\propto n_{tr}$，灵敏度判据取
  $$q_{req}=\frac{T-\mathrm{bg}}{n_{tr}}\quad(\text{每 SPAD 每发需额外贡献的点亮概率，越小越灵敏})$$

### 两种横轴口径，回答的是不同问题

| 口径 | 含义 |
|---|---|
| **同一片天光**（横轴 = 照度 klux） | 各配置 bg 不同，$\mathrm{bg}\propto n_{tr}$ |
| **同一个 bg**（横轴 = 各配置自身 bg） | 同一个 bg 意味着各配置处在**不同**天光下 |

### 结论

1. **纯噪声阈值只取决于 $n_{tr}=n_{pix}\times N_{shots}$**。引擎把"SPAD 数"和"shot 数"折进同一个
   轨迹数维度，MC 实证：3×6@N=6 与 3×9@N=4（都是 $n_{tr}=108$）在 24 档上 **T@1% 最大差 0 计数**。
2. 同一片天光下 3×6@N=2 vs 3×9@N=4：绝对阈值前者低一半以上（底噪只有 1/3），
   但折算成 $q_{req}$ 前者要差 **1.8–2.0 倍**。
3. 同一个 bg 下这一对差距放大到 **2.3–2.7 倍**（把 3×6@N=2 拉到同样 bg，等于让它处在 3 倍强的天光里）。
4. 设计取舍：**宏像元缩小 1.5 倍，用 1.5 倍发数可以精确换回同样的噪声性能**（$n_{tr}$ 是不变量），
   代价是帧率；换来的是 x 方向 1.5 倍角分辨。
"""

M16_CODE = r'''
# ---- 模块 16：宏像元 3×9 vs 3×6（读独立脚本的缓存，画浓缩版）----
MACRO_CACHE = "compare_macro_3x9_vs_3x6_cache.npz"
MACRO_CFGS = [(9, 3, 1), (9, 3, 2), (9, 3, 4),
              (6, 3, 1), (6, 3, 2), (6, 3, 4), (6, 3, 6)]
MACRO_COLOR = {(9, 3, 1): "#9ecae1", (9, 3, 2): "#4292c6", (9, 3, 4): "#08519c",
               (6, 3, 1): "#fcae91", (6, 3, 2): "#fb6a4a", (6, 3, 4): "#cb181d",
               (6, 3, 6): "#67000d"}
MACRO_LS = {9: "-", 6: "--"}
N_PIX_REF_M16 = 27
E_LAMBDA_100KLUX = 0.68


def _macro_label(cfg):
    nx, ny, n = cfg
    return f"{ny}×{nx}（{nx*ny} SPAD） N={n}"


def _macro_ntr(cfg):
    nx, ny, n = cfg
    return float(nx * ny * n)


if not os.path.exists(MACRO_CACHE):
    print(f"未找到 {MACRO_CACHE}。要生成请在命令行跑：")
    print('  $env:PYTHONIOENCODING="utf-8"')
    print("  python compare_macro_3x9_vs_3x6.py --workers 20 --n-mc 200000")
else:
    _mz = np.load(MACRO_CACHE)
    _mfars = list(np.asarray(_mz["fars"], float))
    _mim = int(np.argmin(np.abs(np.asarray(_mfars) - 0.01)))   # FAR=1%
    _mT = np.asarray(_mz[f"thr_{_mim}"], float)
    _mbg = np.asarray(_mz["bg_mc"], float)
    _mok = np.asarray(_mz["done"], bool)
    _mkl = np.asarray(_mz["e_lambda"], float).max(axis=0) / E_LAMBDA_100KLUX * 100
    print(f"模块 16 载入 {MACRO_CACHE}：{_mT.shape[1]} 档环境光 × {_mT.shape[0]} 种配置，"
          f"每档 {int(_mz['n_mc']):,} MC，FAR={_mfars[_mim]:.0%}")

    def _macro_curve(cfg):
        """换成以【该配置自己的 bg】为自变量，返回 (bg, q_req, T)，bg 升序。"""
        ci = MACRO_CFGS.index(cfg)
        m = _mok[ci]
        x = _mbg[ci][m]
        q = (_mT[ci][m] - x) / _macro_ntr(cfg)
        o = np.argsort(x)
        return x[o], q[o], _mT[ci][m][o]

    fig, ax = plt.subplots(1, 4, figsize=(22.0, 5.0))

    # ① 同一片天光下的 bg
    for cfg in MACRO_CFGS:
        ci = MACRO_CFGS.index(cfg); m = _mok[ci]
        ax[0].plot(_mkl[m], _mbg[ci][m], MACRO_LS[cfg[0]], color=MACRO_COLOR[cfg],
                   lw=1.7, label=_macro_label(cfg))
    ax[0].set_xlabel("环境光照度 [klux]（等价于 p_eq）")
    ax[0].set_ylabel("bg [计数/bin]")
    ax[0].set_title("① 同一片天光下的底噪 bg = n_tr·p_eq", fontsize=10.5)
    ax[0].grid(alpha=0.3); ax[0].legend(fontsize=6.8)

    # ②③ 以自身 bg 为横轴
    for idx, (yf, ylab, ttl) in enumerate([
        (lambda cfg: _macro_curve(cfg)[2], "T@FAR=1% [计数]",
         "② 阈值 T vs 自身 bg（7 条几乎重合）"),
        (lambda cfg: _macro_curve(cfg)[1], "q_req = (T−bg)/n_tr",
         "③ 所需信号 q_req vs 自身 bg"),
    ], start=1):
        for cfg in MACRO_CFGS:
            x, q, T = _macro_curve(cfg)
            ax[idx].plot(x, yf(cfg), MACRO_LS[cfg[0]], color=MACRO_COLOR[cfg],
                         lw=1.7, label=_macro_label(cfg))
        ax[idx].set_xscale("log")
        ax[idx].set_xlabel("bg = 该配置 hist_add 统计窗每 bin 均值 [计数]")
        ax[idx].set_ylabel(ylab); ax[idx].set_title(ttl, fontsize=10.5)
        ax[idx].grid(alpha=0.3, which="both"); ax[idx].legend(fontsize=6.8)

    # ④ 成对比值（插到公共 bg 网格）
    for (ca, cb), col in zip(
            [((6, 3, 2), (9, 3, 4)), ((6, 3, 4), (9, 3, 4)),
             ((6, 3, 2), (9, 3, 2)), ((6, 3, 6), (9, 3, 4))],
            ["tab:red", "tab:orange", "tab:green", "tab:purple"]):
        bga, qa, _ = _macro_curve(ca)
        bgb, qb, _ = _macro_curve(cb)
        lo, hi = max(bga.min(), bgb.min()), min(bga.max(), bgb.max())
        if not (hi > lo):
            continue
        xs = np.geomspace(lo, hi, 80)
        ax[3].plot(xs, np.interp(xs, bga, qa) / np.interp(xs, bgb, qb), "-", lw=1.9,
                   color=col, label=f"{_macro_label(ca)} ÷ {_macro_label(cb)}")
    ax[3].axhline(1.0, color="k", lw=1.0, ls=":")
    ax[3].set_xscale("log")
    ax[3].set_xlabel("公共 bg（两个配置都工作在这个底噪上）")
    ax[3].set_ylabel("所需信号强度之比（<1 表示前者更灵敏）")
    ax[3].set_title("④ 同一 bg 下成对比较", fontsize=10.5)
    ax[3].grid(alpha=0.3, which="both"); ax[3].legend(fontsize=6.8)

    fig.suptitle("模块 16　宏像元 3×9（27 SPAD）vs 3×6（18 SPAD）　FAR=1%，"
                 "噪声与信号均按每 SPAD 均匀；阈值只取决于 n_tr = n_pix×N_shots",
                 fontsize=12.5)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    fig.savefig("pod_v20_m16_macro_compare.png", dpi=120, bbox_inches="tight")
    plt.show()

    # 等价性核对 + 同 bg 对比表
    i96, i94 = MACRO_CFGS.index((6, 3, 6)), MACRO_CFGS.index((9, 3, 4))
    mm = _mok[i96] & _mok[i94]
    print(f"[等价性核对] n_tr=108 的两种实现（3×6@N=6 vs 3×9@N=4）："
          f"T@1% 最大差 {np.abs(_mT[i96][mm]-_mT[i94][mm]).max():.0f} 计数，"
          f"bg 最大差 {np.abs(_mbg[i96][mm]-_mbg[i94][mm]).max():.4f}")
    print("\n同一个 bg 下 3×6@N=2 vs 3×9@N=4（FAR=1%）")
    print(f"{'bg':>6}{'T(3×6,N=2)':>12}{'T(3×9,N=4)':>12}"
          f"{'q_req 比':>12}")
    bga, qa, Ta = _macro_curve((6, 3, 2))
    bgb, qb, Tb = _macro_curve((9, 3, 4))
    for x in (2, 4, 6, 8, 10, 12, 14, 16):
        if not (bga.min() <= x <= bga.max() and bgb.min() <= x <= bgb.max()):
            continue
        print(f"{x:6.1f}{np.interp(x,bga,Ta):12.2f}{np.interp(x,bgb,Tb):12.2f}"
              f"{np.interp(x,bga,qa)/np.interp(x,bgb,qb):12.3f}")
'''

M17_MD = r"""
## 模块 17 — ★ v20：理论汇总（为什么这些曲线长这样）

本节把 v11 期间做的三份理论工作收进正本，作为前面所有图的解释框架。
完整推导见仓库根目录的两份 Markdown：`theory_peak_bg_multishot.md` 与 `theory_engine_equivalence.md`。

---

### 17.1　为什么同 bg 下不同 N 的阈值曲线不一样

**单 bin 分布是二项而不是泊松。** SPAD 是 1 bit 器件，一个 bin 在一发里最多点亮一次，
$N$ 发累加后每个 bin 服从 $\mathrm{Bin}(n_{tr},p_{eq})$，$n_{tr}=27N$，$\mathrm{bg}=n_{tr}p_{eq}$。于是

$$\sigma^2 = \mathrm{bg}\left(1-\frac{\mathrm{bg}}{n_{tr}}\right),\qquad
F=\frac{\sigma^2}{\mu}=1-\frac{\mathrm{bg}}{n_{tr}}<1$$

这就是**亚泊松**（Fano 因子小于 1）。**同一个 bg 下，$N$ 越大 $n_{tr}$ 越大，方差压缩越弱、$\sigma$ 越大**，
所以阈值 $T\approx\mu+z\sigma$ 也越高。这直接解释了模块 12B 里三条曲线的排序。

> **反直觉但重要**：曾经的直觉是"$N=4$、每发 noise=1 时，各发的峰不会落在同一个 bin，
> 所以叠加后的 peak 应该比 $N=1$、noise=4 时更小"。
> **这是错的。** `check_same_bg_two_ways.py` 的定向 MC 表明，同 bg 下 $N$ 大的那一边 peak **更大**。
> 原因就是上面的方差公式：$N=1$ 时 $n_{tr}=27$，$p_{eq}$ 被逼到很高，方差被压得很狠；
> $N=4$ 时 $n_{tr}=108$，同样的 bg 只需要 1/4 的 $p_{eq}$，压缩弱得多。
> "峰不对齐"确实存在，但它影响的是**信号**（相干累加），对**噪声底**的极值统计不起主导作用。

**peak 是 152 个 bin 的极大值**，要用极值统计。但 bin 之间**不独立**：过阈窗 $T_{OVER}\approx 8$ ns
让一次雪崩同时点亮相邻 8 个采样点，造成**强正相关**（`check_bin_correlation.py` 实测 ACF 在
lag 1–7 都显著为正）。因此有效独立 bin 数 $M_{eff}$ 远小于名义的 152。

> **更正记录**：早期文档里写成"死时间造成相邻 bin 负相关"，方向是反的。实测是**正相关**。

**三层阈值模型**（精度递增）：高斯极值近似 → 大偏差（LD）近似 → 精确二项分位。
在关心的 FAR 范围内 LD 近似已经足够，细节见 `theory_peak_bg_multishot.md` 第 5–8 节。

**结论**：阈值比值 $\rho_{N/1}=T_N/T_1$ **不是常数**。因为
$T\approx a\,\mathrm{bg}+z\sqrt{b\,\mathrm{bg}}$ 里 $\sqrt{\mathrm{bg}}$ 项的权重随 bg 变化，
只有 bg 很大时 $\rho$ 才趋于 $a_N/a_1$。模块 10 的残差图量化了这一点。

---

### 17.2　为什么快速引擎和 v45 的逐光子 RC 引擎是一致的

用户曾担心 PoD 把 SPAD 简化成了"8 ns 硬死时间"。**并没有。**
`check_engine_vs_v45.py` 做了三级核对：源码逐行、同种子比特级（60/60 条轨迹逐 bin 完全相同）、统计级。
RC 恢复是显式建模的：$V_{ov}$ 占比 $=1-e^{-\Delta t/\tau_{RC}}$，触发概率 $\mathrm{PDE}_{max}\cdot g(V_{ov})$，
每次雪崩把 $V_{ov}$ 打回 0。

理论上的等价链条（`theory_engine_equivalence.md`）：

1. **泊松稀释**：入射光子是泊松流，探测是独立稀释 → 探测事件仍是泊松流，强度 $r_{\det}$。
2. **更新过程**：每次雪崩把器件状态复位，下一次雪崩只依赖"距上次雪崩多久" → 更新过程，
   风险函数 $h(t)=r_{\det}\,g\!\left(1-e^{-t/\tau_{RC}}\right)$。
3. **逆变换抽样**：间隔分布 $F(t)=1-e^{-\int_0^t h}$ 可以预先做成 $H^{-1}$ 查找表，$O(1)$ 抽样。
4. **并集恒等式**：每次雪崩点亮一段长 $T_{OVER}$ 的窗；这些窗的并集用差分数组 $O(n)$ 求出，
   与逐点判断等价。
5. **更新-回报定理**：平衡态下一个 bin 被点亮的概率
   $p_{bin}=T_{OVER}/\mathbb{E}[X]$ 的连续时间形式，给出解析的 $p_{eq}$。
6. **$n_{tr}$ 折叠**：各 SPAD、各 shot 独立同分布，纯噪声下只以 $n_{tr}=n_{pix}\times N_{shots}$ 出现
   （模块 16 的实证：3×6@N=6 与 3×9@N=4 阈值完全相同）。
7. **离散化误差**：逐光子引擎用 $dt$ 步进，误差是 $O(dt^2)$；实测 $dt$ 减半误差降到约 1/4
   （比值 3.89 / 3.97 / 3.99 / 4.00，干净的二阶收敛）。

**待查项**：`theory_engine_equivalence.py` 的 T5b 里，步进引擎 `binary_macro_stepping` 在
$dt=800/200$ ps 两档都比解析值**偏低**（2.2σ / 0.7σ），怀疑是"先出 bin、再处理本步雪崩"的半步对齐。
纯噪声扫描用的是快速引擎，不受影响；但 **PoD 信号支路用的正是这个步进引擎**，建议加大样本复核。

---

### 17.3　两个必须先讲清的口径

1. **同天光 vs 同 bg**（模块 16）：同一片天光下 $n_{tr}$ 小的配置 bg 天然更低、占便宜；
   拉到同一个 bg 才是公平的探测能力比较，但这意味着两者处在不同强度的环境光里。
   **引用任何"差几倍"的数字时都必须注明是哪一种口径。**
2. **信号均匀 vs 按像斑加权**（模块 16）：按像斑
   `FX=[0.0014,0.0152,0.084,0.234,0.330,…]` 加权时，9 列砍到 6 列只丢 1.8% 信号 → 3×6 更灵敏；
   按均匀处理时 3×6 只收到 2/3 信号 → 结论反转成 3×6 差 1.2 倍。
   **两个结论都对，只是口径不同。** 本仓库当前统一采用**均匀**口径。
3. **净峰高 $S_{net}$ vs 净增量 $\Delta S$**（模块 13）：即使没有信号，在 15 个 bin 的信号窗里取最大值
   也已经高于 bg（极值抬升约 2σ），所以 $S_{net}=\mathrm{peak}-\mathrm{bg}$ 会**高估**信号贡献。
"""


# ============================================================ 主流程
def main():
    if not SRC.exists():
        raise SystemExit(f"找不到 {SRC}")
    shutil.copyfile(SRC, DST)
    nb = json.loads(DST.read_text(encoding="utf-8"))

    # ---------- 1. 标题 ----------
    s0 = cell_src(nb, 0)
    s0 = s0.replace("# PoD_esti v11", "# PoD_esti v20")
    s0 += r"""

---

## ★ v20 相对 v11 的变化

**v11 的内容一条都没有删**，只在后面追加了 7 个模块，并把两处工程问题修掉：

| 模块 | 内容 | 对应需求 |
|---|---|---|
| 11 | 每条 hist 内 152 bin 的 std 均值、peak 均值、peak 标准差 随 bg | 需求 1 后半 |
| 12 | **连续（实数）阈值曲线** —— 折线，不再是整数阶梯 | 需求 1「不要画阶梯」 |
| 13 | FAR=5% / 1% / 100 ppm 下 PoD50 与 PoD90 所需信号的均值 | 需求 2 |
| 14 | 平方反比测远（纯 $1/D^2$ 与含大气衰减两种口径） | 需求 3 |
| 15 | 同信号强度、不同 bg 时 peak 分布怎么变（均值是否只是加 bg、std 怎么变） | 需求 4 |
| 16 | 宏像元 3×9 vs 3×6 阈值对比 | 本轮旁支分析并入正本 |
| 17 | 理论汇总：阈值倍数模型 + 引擎一致性模型 + 口径说明 | 同上 |

工程改动（不动任何物理参数）：

* 缓存名换成 `pod_esti_v20_*`，同时把 v11 的缓存登记为 **fallback**，
  载入成功后同步写回 v20 主名（规则三）。**已跑完的 v11 噪声/PoD 结果不会重算。**
* 模块 9.3 的信号扫描从 notebook 内**串行**改成调用多进程脚本 `run_pod_v20_sig_scan.py`
  （48 bg × 9 boost × 3 个 N，ProcessPool 20 进程 + 断点续跑）。
"""
    set_src(nb, 0, s0)

    # ---------- 2. cell 2：缓存名与 fallback ----------
    sub(nb, 2, 'CACHE_NOISE = "pod_esti_v11_cache_noise.npz"',
        'CACHE_NOISE = "pod_esti_v20_cache_noise.npz"')
    sub(nb, 2, 'CACHE_POD   = "pod_esti_v11_cache_pod.npz"',
        'CACHE_POD   = "pod_esti_v20_cache_pod.npz"')
    sub(nb, 2, 'CACHE_SIG   = "pod_esti_v11_cache_signal.npz"  # 模块 9 固定信号',
        'CACHE_SIG   = "pod_esti_v20_cache_signal.npz"  # 模块 9.3 / 15 固定信号')
    sub(nb, 2,
        'CACHE_NOISE_FALLBACK = []  # ★ v11 禁止复用 v10/旧缓存，全量重算\n'
        'CACHE_POD_FALLBACK   = []',
        '# ★ v20：物理内核与网格与 v11 逐字相同 → v11 缓存可直接复用（读到后同步写回 v20 主名）\n'
        'CACHE_NOISE_FALLBACK = ["pod_esti_v11_cache_noise.npz"]\n'
        'CACHE_POD_FALLBACK   = ["pod_esti_v11_cache_pod.npz"]')
    sub(nb, 2, 'CACHE_NOISE_CKPT = "pod_esti_v11_cache_noise.partial.npz"',
        'CACHE_NOISE_CKPT = "pod_esti_v20_cache_noise.partial.npz"')
    sub(nb, 2, 'CACHE_POD_CKPT   = "pod_esti_v11_cache_pod.partial.npz"',
        'CACHE_POD_CKPT   = "pod_esti_v20_cache_pod.partial.npz"')
    sub(nb, 2, 'CACHE_SIG_CKPT   = "pod_esti_v11_cache_signal.partial.npz"',
        'CACHE_SIG_CKPT   = "pod_esti_v20_cache_signal.partial.npz"')
    sub(nb, 2, '  ★ v11 统一 BG_GRID', '  ★ v20 统一 BG_GRID')
    sub(nb, 2, '★ 噪声可从 v04 fallback 迁入；PoD 因 FAR 变更必须用 v05 新缓存重跑',
        '★ v11 缓存已登记为 fallback：读到即同步写回 v20 主名，不会重算')

    # ---------- 3. 自动开跑的脚本名 ----------
    s17 = cell_src(nb, 17).replace("run_pod_v11_noise_scan.py",
                                   "run_pod_v20_noise_scan.py")
    s17 = s17.replace("# ---- ★ v11：统一 BG_GRID", "# ---- ★ v20：统一 BG_GRID")
    s17 = s17.replace("build_pod_core_v11.py", "build_pod_core_v20.py")
    s17 = s17.replace("v11 噪声扫描", "v20 噪声扫描")
    s17 = s17.replace("完整 v11 噪声缓存", "完整 v20 噪声缓存")
    s17 = s17.replace("★ v11：按统一 bg 网格扫", "★ v20：按统一 bg 网格扫")
    s17 = s17.replace("run_noise_scan_v11_bg", "run_noise_scan_v20_bg")
    s17 = s17.replace("run_noise_scan_v11_amb", "run_noise_scan_v20_amb")
    s17 = s17.replace("完整，v11 统一 bg", "完整，v20 统一 bg")
    set_src(nb, 17, s17)

    s25 = cell_src(nb, 25).replace("run_pod_v11_pod_scan.py", "run_pod_v20_pod_scan.py")
    s25 = s25.replace("★ v11：按当前 N 仿", "★ v20：按当前 N 仿")
    s25 = s25.replace("# ★ v11：缺 PoD 缓存时", "# ★ v20：缺 PoD 缓存时")
    set_src(nb, 25, s25)

    # ---------- 4. 模块 9.3：改用多进程脚本 ----------
    set_src(nb, 35, M93_CODE.strip("\n") + "\n")
    sub(nb, 34, "数据由下方 cell 现场计算（或读 `CACHE_SIG`）；**不复用** v05/旧 peak_vs_noise 缓存。",
        "★ v20：数据由多进程脚本 `run_pod_v20_sig_scan.py` 产出（缓存 `pod_esti_v20_cache_signal.npz`，\n"
        "48 bg × 9 boost × N∈{1,2,4} × 8,000 MC，ProcessPool 20 进程 + 断点续跑）。\n"
        "缓存里存的是**完整 peak 分布**（bincount），模块 15 会用它看分布形状怎么变。")

    # ---------- 5. 全局图名 v11 → v20 ----------
    for i, c in enumerate(nb["cells"]):
        s = "".join(c["source"])
        if "pod_v11_" in s:
            set_src(nb, i, s.replace("pod_v11_", "pod_v20_"))

    # ---------- 6. 追加新模块 ----------
    nb["cells"] += [
        md(M11_MD), code(M11_CODE),
        md(M12_MD), code(M12_CODE),
        md(M13_MD), code(M13_CODE),
        md(M14_MD), code(M14_CODE),
        md(M15_MD), code(M15_CODE),
        md(M16_MD), code(M16_CODE),
        md(M17_MD),
    ]

    DST.write_text(json.dumps(nb, ensure_ascii=False, indent=1), encoding="utf-8")
    print(f"wrote {DST}（{len(nb['cells'])} cells）")


# 模块 9.3 的新实现（改用多进程脚本 + v20 缓存）
M93_CODE = r'''
# ---- 模块 9.3：固定信号 × 统一 bg 网格（★ v20：改用多进程脚本 run_pod_v20_sig_scan.py）----
BOOST_LIST_M9 = np.round(np.arange(0.0, 0.032 + 1e-12, 0.004), 6)
N_MC_SIG_M9 = 8000
SIG_SCRIPT = "run_pod_v20_sig_scan.py"


def _try_load_sig_m9(path):
    """载入信号扫描缓存；要求 boost 网格、bg 网格、MC 条数、N 列表全部对得上，且全部档位已完成。"""
    if not (USE_CACHE and os.path.exists(path)):
        return None
    try:
        z = np.load(path)
        if (int(z["n_mc"]) != N_MC_SIG_M9
                or not np.allclose(z["grid_key"], BG_GRID)
                or not np.array_equal(z["n_shots_list"], np.asarray(N_SHOTS_LIST))
                or not np.allclose(z["boosts"], BOOST_LIST_M9)):
            return None
        for n in N_SHOTS_LIST:
            cnt = z[f"peak_cnt_{n}"]
            done = z[f"done_{n}"] if f"done_{n}" in z.files else (cnt.sum(axis=2) > 0)
            if not np.all(done):
                return None
        return z
    except Exception:
        return None


_zsig = None
for _cand in [CACHE_SIG, CACHE_SIG_CKPT]:
    _zsig = _try_load_sig_m9(_cand)
    if _zsig is not None:
        print(f"模块 9.3 命中缓存 {_cand}")
        break

if _zsig is None:
    import sys
    print("=" * 72)
    print("未找到完整 v20 信号缓存 → 自动调用多进程扫描（ProcessPool，吃满 CPU）")
    print(f"规模：{len(BG_GRID)} bg × {len(BOOST_LIST_M9)} boost × N={N_SHOTS_LIST}"
          f" × {N_MC_SIG_M9:,} MC")
    print("=" * 72)
    _rc = _run_cmd_stream([sys.executable, SIG_SCRIPT,
                           "--workers", str(int(N_WORKERS)),
                           "--n-mc", str(int(N_MC_SIG_M9))])
    if _rc != 0:
        raise RuntimeError(f"{SIG_SCRIPT} 失败，请查看上方进度输出")
    _zsig = _try_load_sig_m9(CACHE_SIG)
    if _zsig is None:
        raise RuntimeError(f"多进程信号扫描结束但仍无法载入完整缓存 {CACHE_SIG}")
    print(f"多进程信号扫描完成，已载入 {CACHE_SIG}")

SIG_M9 = {}
for n in N_SHOTS_LIST:
    cnt = np.asarray(_zsig[f"peak_cnt_{n}"])
    mu = np.zeros((len(BOOST_LIST_M9), len(BG_GRID)))
    sd = np.zeros_like(mu)
    for i in range(len(BOOST_LIST_M9)):
        for k in range(len(BG_GRID)):
            s = peak_stats_from_cnt(cnt[i, k])
            mu[i, k] = s["mean"]; sd[i, k] = s["std"]
    SIG_M9[n] = dict(peak_cnt=cnt, peak_mean=mu, peak_std=sd)
print(f"模块 9.3 数据就绪：{len(BOOST_LIST_M9)} boost × {len(BG_GRID)} bg × "
      f"N={N_SHOTS_LIST}，每档 {N_MC_SIG_M9:,} MC（含完整 peak 分布）")

_bg_t = BG_GRID
fig, ax = plt.subplots(2, len(N_SHOTS_LIST), figsize=(5.2 * len(N_SHOTS_LIST), 8.0),
                       sharex=True)
for j, n in enumerate(N_SHOTS_LIST):
    for i, b in enumerate(BOOST_LIST_M9):
        ax[0, j].plot(_bg_t, SIG_M9[n]["peak_mean"][i], lw=1.2, label=f"b={b:g}")
        ax[1, j].plot(_bg_t, SIG_M9[n]["peak_std"][i], lw=1.2)
    ax[0, j].set_title(f"N={n} peak均值"); ax[0, j].set_ylabel("peak mean")
    ax[1, j].set_title(f"N={n} peak std"); ax[1, j].set_xlabel("bg")
    ax[1, j].set_ylabel("peak std")
    ax[0, j].legend(fontsize=7, ncol=2); ax[0, j].grid(True, alpha=0.3)
    ax[1, j].grid(True, alpha=0.3)
fig.suptitle("模块 9.3a　固定信号：peak 均值/std 随 bg（v20 统一步长）", fontsize=12)
fig.tight_layout()
fig.savefig("pod_v20_m9_sig_vs_bg.png", dpi=120, bbox_inches="tight")
plt.show()

# 平移检验：相对 boost=0 的均值增量
fig, ax = plt.subplots(1, len(N_SHOTS_LIST), figsize=(5.0 * len(N_SHOTS_LIST), 4.2))
if len(N_SHOTS_LIST) == 1:
    ax = [ax]
for j, n in enumerate(N_SHOTS_LIST):
    base = SIG_M9[n]["peak_mean"][0]
    for i, b in enumerate(BOOST_LIST_M9[1:], start=1):
        ax[j].plot(_bg_t, SIG_M9[n]["peak_mean"][i] - base, lw=1.3, label=f"b={b:g}")
    ax[j].set_xlabel("bg"); ax[j].set_ylabel("Δpeak_mean"); ax[j].set_title(f"N={n}")
    ax[j].legend(fontsize=7); ax[j].grid(True, alpha=0.3)
fig.suptitle("模块 9.3b　分布平移检验", fontsize=12)
fig.tight_layout()
fig.savefig("pod_v20_m9_shift.png", dpi=120, bbox_inches="tight")
plt.show()

print("模块 9.3 线性斜率摘要（Δpeak / boost，对 bg 平均）：")
for n in N_SHOTS_LIST:
    base = SIG_M9[n]["peak_mean"][0]
    for i, b in enumerate(BOOST_LIST_M9[1:], start=1):
        if b <= 0:
            continue
        slope = np.mean((SIG_M9[n]["peak_mean"][i] - base) / b)
        print(f"  N={n} boost={b:g} → 平均斜率 {slope:.2f}")
'''


if __name__ == "__main__":
    main()
