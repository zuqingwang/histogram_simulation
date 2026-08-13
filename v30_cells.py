# -*- coding: utf-8 -*-
"""PoD_esti v30 的**新分析层** cell 源码。

被 build_pod_esti_v30.py 引用。物理内核不在这里，由 build 脚本从
PoD_esti_v20.ipynb 逐字抽取（保证物理一条不改）。

v30 的结构约定（用户要求）：
  每个模块 = [计算/载入缓存 cell] → [绘图参数 cell] → [绘图 cell]
  绘图 cell 只读缓存里的结果，绝不重算。
  markdown 只留基本描述 + 缩写；理论与公式推导移到 theory_PoD_esti_v30.md。
"""

# ============================================================ 模块 5：纯噪声 bg 扫描

M5_MD = r"""
## 模块 5 — 纯噪声 bg 扫描（全项目唯一的大重算）

对统一 bg 网格（0.25 → 12，步长 0.25，共 48 档）逐档跑纯环境光 MC，
得到每档的 **peak 分布（bincount）**、**单条 hist 内 152 bin 的 std**，
并由 peak 分布反解各 FAR 目标下的**检测阈值 T**。

本模块的结果 `NOISE_RES` / `THRESH` 是模块 9 / 10 / 11 的**唯一数据源**，后面不再重算。

- **bg**：`hist_add`（累加 N 发后的直方图）在统计窗内每个 1 ns bin 的**平均计数**。全项目横轴统一用它。
- **noise**：单发 `hist_i` 每 bin 的平均计数，`noise = bg / N_shots`。**只用于描述单次直方图的底噪**，不再作为横轴。
- **peak**：一条 `hist_add` 在统计窗 152 个 bin 内的**最大计数**。
- **FAR**（False Alarm Rate，虚警率）：`P(peak ≥ T)`，即一个统计窗出一个噪点的概率。

> **阈值 T 为什么是整数**：`hist_add` 是 27×N 条二值（0/1）轨迹求和，取值只能是整数
> 0…27N，判定就是 `peak >= T`。所以 `T=10.25` 与 `T=11` 是同一个判定，
> 0.25 粒度的阈值在计数轴上不存在。曲线看着像阶梯，是因为 **bg 以 0.25 连续步进、而 T 只能按 1 计数跳**。
"""

M5_PARAM = r'''
# ==================== 模块 5 绘图参数 ====================
# 只影响画图，不影响任何计算；改完重跑本 cell 与下一个绘图 cell 即可。

M5_FARS        = [0.05, 0.01]   # 【叠加曲线】bg–peak 图上只叠这两条阈值（5% / 1%）
M5_XLIM        = (0.0, 12.25)   # 【上下两排横轴】bg 范围 [计数/bin]
M5_YLIM_PEAK   = None           # 【上排纵轴】peak 计数范围；None = 自动贴合数据
M5_SHOW_CAP    = False          # 【上排】是否画二值硬上限 n_tr。N=4 时上限 108 离数据很远，
                                #   画上去会把曲线压扁，所以默认关掉
M5_BAND_LO     = 0.01           # 【上排阴影带】下分位（0.01 = 1% 分位）
M5_BAND_HI     = 0.99           # 【上排阴影带】上分位
M5_STRIP_YLIM  = None           # 【下排纵轴】密度条带的 peak 范围；None = 自动按数据裁
M5_STRIP_QCUT  = 0.9999         # 【下排纵轴自动裁剪】保留到该分位，压掉长尾空白
M5_STRIP_VMAXQ = 0.98           # 【下排色标上限】取 pmf 的该分位，避免个别高峰吃掉全部色阶
M5_STRIP_CMAP  = "magma"        # 【下排配色】
M5_FIGSIZE     = (5.6 * len(N_SHOTS_LIST), 8.6)
M5_TABLE_STEP  = 6              # 【文字表】每隔几档 bg 打一行
'''

M5_PLOT = r'''
# ==================== 模块 5 绘图 ====================
# 只读 NOISE_RES / THRESH，不重算。
_bg = np.asarray(BG_GRID, float)
fig, axes = plt.subplots(2, len(N_SHOTS_LIST), figsize=M5_FIGSIZE, sharex=True)
axes = np.atleast_2d(axes)

for j, n in enumerate(N_SHOTS_LIST):
    R, Tr = NOISE_RES[n], THRESH[n]
    st = [peak_stats_from_cnt(c) for c in R["peak_cnt"]]
    pk_mean = np.array([s["mean"] for s in st])
    pk_lo = np.array([_quantile_from_cnt(c, M5_BAND_LO) for c in R["peak_cnt"]])
    pk_hi = np.array([_quantile_from_cnt(c, M5_BAND_HI) for c in R["peak_cnt"]])

    # ---- 上排：bg–peak 曲线 + 1%/5% 阈值 ----
    a = axes[0, j]
    a.fill_between(_bg, pk_lo, pk_hi, color=_COLORS_N[n], alpha=0.18,
                   label=f"peak {M5_BAND_LO:.0%}–{M5_BAND_HI:.0%} 分位带")
    a.plot(_bg, pk_mean, "-", color=_COLORS_N[n], lw=2.0, label="peak 均值")
    a.plot(_bg, _bg, ":", color="0.45", lw=1.3, label="参考 T=bg")
    for far, ls in zip(M5_FARS, ["--", "-."]):
        a.plot(_bg, Tr["T" + FAR_TAG[far]], ls, color="k", lw=1.5,
               label=f"阈值 {FAR_LABEL[far]}")
    if M5_SHOW_CAP:
        a.axhline(R["n_tr"], color="0.3", ls=(0, (1, 1)), lw=1.1,
                  label=f"二值硬上限 {R['n_tr']}")
    a.set_title(f"N_shots={n}　bg–peak（阈值为整数计数，故呈阶梯）", fontsize=10.5)
    a.set_ylabel("peak / 阈值 [计数/bin]")
    if M5_YLIM_PEAK is not None:
        a.set_ylim(*M5_YLIM_PEAK)
    else:   # 贴着数据自动定范围，别让远处的硬上限把曲线压扁
        a.set_ylim(0, max(pk_hi.max(), Tr["T" + FAR_TAG[M5_FARS[0]]].max()) * 1.28)
    a.set_xlim(*M5_XLIM)
    a.grid(alpha=0.3); a.legend(fontsize=7.2, ncol=2, loc="upper left")

    # ---- 下排：peak 密度条带 ----
    b = axes[1, j]
    pmf = np.array([c / max(c.sum(), 1) for c in R["peak_cnt"]])      # [bg, peak]
    if M5_STRIP_YLIM is not None:
        v_hi = int(M5_STRIP_YLIM[1])
    else:
        v_hi = int(max(_quantile_from_cnt(c, M5_STRIP_QCUT)
                       for c in R["peak_cnt"])) + 2
    v_hi = min(v_hi, pmf.shape[1] - 1)
    sub = pmf[:, :v_hi + 1]
    vmax = np.quantile(sub[sub > 0], M5_STRIP_VMAXQ) if np.any(sub > 0) else 1.0
    mesh = b.pcolormesh(_bg, np.arange(v_hi + 1), sub.T, cmap=M5_STRIP_CMAP,
                        vmin=0.0, vmax=vmax, shading="nearest")
    fig.colorbar(mesh, ax=b, pad=0.01).set_label("P(peak = v)", fontsize=8)
    for far, ls in zip(M5_FARS, ["--", "-."]):
        b.plot(_bg, Tr["T" + FAR_TAG[far]], ls, color="w", lw=1.4,
               label=f"阈值 {FAR_LABEL[far]}")
    b.set_title(f"N_shots={n}　peak 概率密度条带", fontsize=10.5)
    b.set_xlabel("bg [计数/bin]"); b.set_ylabel("peak [计数]")
    if M5_STRIP_YLIM is not None:
        b.set_ylim(*M5_STRIP_YLIM)
    b.set_xlim(*M5_XLIM)
    b.legend(fontsize=7.2, loc="upper left")

fig.suptitle(f"模块 5　纯噪声 bg 扫描：peak 分布与检测阈值（每档 {N_MC_NOISE:,} MC）",
             fontsize=12.5)
fig.tight_layout(rect=[0, 0, 1, 0.94])
fig.savefig("pod_v30_m5_bg_peak.png", dpi=120, bbox_inches="tight")
plt.show()

print("=" * 108)
print("模块 5 抽样表：实测 bg / peak 均值 / 单条 hist 内 std / 1% 与 5% 阈值")
_hdr = f"{'bg':>6}"
for n in N_SHOTS_LIST:
    _hdr += f" | N={n}: bg实测 peakμ histσ  T5%  T1%"
print(_hdr)
for k in range(0, len(_bg), M5_TABLE_STEP):
    row = f"{_bg[k]:6.2f}"
    for n in N_SHOTS_LIST:
        R, Tr = NOISE_RES[n], THRESH[n]
        pm = peak_stats_from_cnt(R["peak_cnt"][k])["mean"]
        row += (f" | {R['noise_mc'][k]:6.3f} {pm:6.2f} {R['hist_std'][k]:6.3f}"
                f" {Tr['T5pct'][k]:4d} {Tr['T1pct'][k]:4d}")
    print(row)
'''


# ============================================================ 模块 6：PoD 临界能量

M6_MD = r"""
## 模块 6 — PoD 临界能量（逐 bg 档求 50% / 90% 交点）

对每个 bg 档、每个 `N_shots`，扫描回波能量倍率 `boost`，用模块 5 的阈值 `T` 判定，
求 **PoD**（Probability of Detection，探测概率）达到 50% / 90% 的临界能量。

只对 `POD_FARS` 这几条 FAR 求临界能量（0.5% / 1% / 5% / 10%）；
更严的 ppm 档阈值仍然算，只是不在这里解 PoD 交点（省一大半机时）。
"""

M6_PARAM = r'''
# ==================== 模块 6 绘图参数 ====================
M6_SHOW_N      = 4               # 【画哪个 N_shots】的 PoD 曲线
M6_SHOW_BG     = [1.0, 6.0, 12.0]  # 【画哪几个 bg 档】
M6_FAR         = 0.01            # 【用哪条 FAR 的阈值】判定
M6_POD_XPAD    = 0.30            # 【上排横轴】交点两侧各留多少 decade（放大交点区）
M6_POD_YLIM    = (-0.02, 1.04)   # 【上排纵轴】PoD 范围
M6_DIST_STYLE  = "fill"          # 【下排画法】"fill"=实心直方图 / "line"=折线
M6_DIST_QLO    = 0.0005          # 【下排横轴】左端取分布的这个分位
M6_DIST_QHI    = 0.9995          # 【下排横轴】右端取分布的这个分位（保证右尾不被切掉）
M6_FIGSIZE     = (5.4 * 3, 8.0)
'''

M6_PLOT = r'''
# ==================== 模块 6 绘图 ====================
_tag6 = FAR_TAG[M6_FAR]
_ks6 = [int(np.argmin(np.abs(np.asarray(BG_GRID, float) - b))) for b in M6_SHOW_BG]
fig, axes = plt.subplots(2, len(_ks6), figsize=M6_FIGSIZE)
axes = np.atleast_2d(axes)

for j, k in enumerate(_ks6):
    key = (M6_SHOW_N, float(NOISE_GRID[M6_SHOW_N][k]))
    r = POD_RES.get(key)
    a_pod, a_dist = axes[0, j], axes[1, j]
    if not r or _tag6 not in r.get("curve", {}):
        a_pod.text(0.5, 0.5, "该档无 PoD 结果", ha="center", va="center")
        continue
    cur = r["curve"][_tag6]
    T = int(r["T_map"][_tag6])
    bgv = float(r["noise"])

    # ---- 上排：PoD–能量曲线，横轴放大到交点附近 ----
    a_pod.semilogx(cur["boost"], cur["pod"], "o-", color=_COLORS_N[M6_SHOW_N],
                   ms=3.2, lw=1.5, label="MC 实测 PoD")
    xr = []
    for lv, c in zip(POD_LEVELS, ["tab:orange", "tab:purple"]):
        rec = r["critical"].get(_tag6, {}).get(f"{lv:.2f}")
        if rec:
            a_pod.axvline(rec["boost"], color=c, ls="--", lw=1.3,
                          label=f"PoD{int(lv*100)} @ boost={rec['boost']:.3g}")
            a_pod.axhline(lv, color=c, ls=":", lw=0.9)
            xr.append(rec["boost"])
    if xr:  # 只画交点前后 M6_POD_XPAD 个 decade，掐掉两头没用的平台
        a_pod.set_xlim(10 ** (np.log10(min(xr)) - M6_POD_XPAD),
                       10 ** (np.log10(max(xr)) + M6_POD_XPAD))
    a_pod.set_ylim(*M6_POD_YLIM)
    a_pod.set_title(f"bg={bgv:.2f}　T={T}（{FAR_LABEL[M6_FAR]}）", fontsize=10.5)
    a_pod.set_xlabel("回波能量倍率 boost"); a_pod.set_ylabel("PoD")
    a_pod.grid(alpha=0.3, which="both"); a_pod.legend(fontsize=7.4)

    # ---- 下排：临界点上的 peak 分布（实心/折线，不再用空心直方图）----
    rec90 = r["critical"].get(_tag6, {}).get("0.90")
    if rec90 and rec90.get("peak_cnt") is not None:
        cnt = np.asarray(rec90["peak_cnt"], float)
        pmf = cnt / max(cnt.sum(), 1.0)
        x = np.arange(pmf.size)
        m = pmf > 0
        if M6_DIST_STYLE == "line":
            a_dist.plot(x[m], pmf[m], "-", color=_COLORS_N[M6_SHOW_N], lw=1.8)
        else:
            a_dist.bar(x[m], pmf[m], width=1.0, color=_COLORS_N[M6_SHOW_N],
                       alpha=0.75, edgecolor="none")
        a_dist.axvline(T, color="k", ls="--", lw=1.6, label=f"阈值 T={T}")
        a_dist.axvline(bgv, color="0.5", ls=":", lw=1.3, label=f"bg={bgv:.2f}")
        # 横轴按分布自身的分位取，并保证阈值 T 一定在视野里（固定 ±N 计数会切掉右尾）
        cum = np.cumsum(pmf)
        lo = float(np.searchsorted(cum, M6_DIST_QLO))
        hi = float(np.searchsorted(cum, M6_DIST_QHI))
        a_dist.set_xlim(min(lo, T) - 1.5, max(hi, T) + 1.5)
        a_dist.set_title(f"PoD90 临界点的 peak 分布（{rec90['pod']:.1%}）", fontsize=10)
        a_dist.set_xlabel("peak [计数]"); a_dist.set_ylabel("概率")
        a_dist.grid(alpha=0.3); a_dist.legend(fontsize=7.4)

fig.suptitle(f"模块 6　PoD–能量交点与临界 peak 分布（N_shots={M6_SHOW_N}）", fontsize=12.5)
fig.tight_layout(rect=[0, 0, 1, 0.94])
fig.savefig("pod_v30_m6_pod_critical.png", dpi=120, bbox_inches="tight")
plt.show()
'''


# ============================================================ 模块 7：全 bg 汇总

M7_MD = r"""
## 模块 7 — 全 bg 档汇总：临界 peak、临界能量与等效测距

把模块 6 的临界点按 bg 排开，给出 **PoD50 / PoD90** 下的临界 peak、临界发射能量，
以及在平方反比下等效到多远（`equiv_distance`）。按用户要求只画 **1% 与 5%** 两条 FAR。
"""

M7_PARAM = r'''
# ==================== 模块 7 绘图参数 ====================
M7_FARS    = [0.05, 0.01]     # 【画哪几条 FAR】只画 5% 与 1%
M7_LEVELS  = [0.50, 0.90]     # 【画哪几个 PoD 水平】
M7_XLIM    = (0.0, 12.25)     # 【横轴】bg 范围
M7_YLIM_E  = None             # 【能量子图纵轴】None = 自动（对数轴）
M7_YLIM_D  = None             # 【测距子图纵轴】None = 自动
M7_FIGSIZE = (16.5, 4.6)
'''

M7_PLOT = r'''
# ==================== 模块 7 绘图 ====================
_LS_LV = {0.50: "--", 0.90: "-"}
fig, ax = plt.subplots(1, 3, figsize=M7_FIGSIZE)
for n in N_SHOTS_LIST:
    for far in M7_FARS:
        tag = FAR_TAG[far]
        for lv in M7_LEVELS:
            arr = _collect_critical(n, tag, lv)
            if not arr.size:
                continue
            lab = f"N={n} {FAR_LABEL[far]} PoD{int(lv*100)}"
            ax[0].plot(arr[:, 0], arr[:, 2], _LS_LV[lv], color=_COLORS_N[n],
                       lw=1.6, alpha=0.9 if far == 0.01 else 0.5, label=lab)
            ax[1].semilogy(arr[:, 0], arr[:, 3], _LS_LV[lv], color=_COLORS_N[n],
                           lw=1.6, alpha=0.9 if far == 0.01 else 0.5, label=lab)
            ax[2].plot(arr[:, 0], arr[:, 4], _LS_LV[lv], color=_COLORS_N[n],
                       lw=1.6, alpha=0.9 if far == 0.01 else 0.5, label=lab)
ax[0].set_ylabel("临界 peak 均值 [计数]"); ax[0].set_title("① 临界 peak", fontsize=11)
ax[1].set_ylabel("$E_{crit}$ [nJ]");      ax[1].set_title("② 临界发射能量", fontsize=11)
ax[2].set_ylabel("等效距离 [m]");          ax[2].set_title("③ 平方反比等效测距", fontsize=11)
for i, a in enumerate(ax):
    a.set_xlabel("bg [计数/bin]"); a.set_xlim(*M7_XLIM)
    a.grid(alpha=0.3, which="both")
if M7_YLIM_E is not None:
    ax[1].set_ylim(*M7_YLIM_E)
if M7_YLIM_D is not None:
    ax[2].set_ylim(*M7_YLIM_D)
# 12 条曲线，图例放进画框会盖住数据 → 全图共用一个图例，摆在下方
_h, _l = ax[0].get_legend_handles_labels()
fig.legend(_h, _l, loc="lower center", ncol=6, fontsize=7.4,
           frameon=False, bbox_to_anchor=(0.5, -0.02))
fig.suptitle("模块 7　全 bg 档汇总（实线 PoD90，虚线 PoD50）", fontsize=12.5)
fig.tight_layout(rect=[0, 0.10, 1, 0.92])
fig.savefig("pod_v30_m7_summary.png", dpi=120, bbox_inches="tight")
plt.show()
'''


# ============================================================ 模块 9：阈值倍数 rho

M9_MD = r"""
## 模块 9 — 同 bg 下不同 N 的阈值倍数 ρ

**问的是什么**：把 N=1 / 2 / 4 拉到**同一个 bg**（同样的底噪水平）上，
它们需要的检测阈值差几倍？定义

`ρ_N = T_N / T_1`（同一个 bg、同一条 FAR）

如果 ρ 是一个与 bg 无关的常数，那么"多发累加要把阈值抬高几倍"就可以用一个数概括；
本模块就是检验这件事，并给出偏离常数有多大。

数据**全部来自模块 5 的 `THRESH`**，本模块不重算、也不重画阈值曲线本身。
"""

M9_PARAM = r'''
# ==================== 模块 9 绘图参数 ====================
M9_FARS    = [0.05, 0.01]    # 【画哪几条 FAR】
M9_XLIM    = (0.0, 12.25)    # 【横轴】bg 范围
M9_YLIM_R  = None            # 【左图纵轴】ρ 范围；None = 自动
M9_YLIM_RES= (-15, 15)       # 【右图纵轴】相对残差 [%]
M9_BG_MIN  = 1.0             # 【统计 ρ̄ 时忽略多小的 bg】小 bg 处 T 只有几个计数，量化噪声极大
M9_FIGSIZE = (13.0, 4.6)
'''

M9_PLOT = r'''
# ==================== 模块 9 绘图 ====================
_bg = np.asarray(BG_GRID, float)
_msk9 = _bg >= M9_BG_MIN
fig, ax = plt.subplots(1, 2, figsize=M9_FIGSIZE)
_rows9 = []
for far in M9_FARS:
    tag = FAR_TAG[far]
    T1 = THRESH[1]["T" + tag].astype(float)
    for n in N_SHOTS_LIST:
        if n == 1:
            continue
        rho = THRESH[n]["T" + tag].astype(float) / np.maximum(T1, 1e-9)
        rho_bar = float(np.nanmean(rho[_msk9]))
        resid = (rho - rho_bar) / rho_bar * 100.0
        ls = "-" if far == 0.01 else "--"
        ax[0].plot(_bg, rho, ls, color=_COLORS_N[n], lw=1.7,
                   label=f"N={n}/1　{FAR_LABEL[far]}")
        ax[0].axhline(rho_bar, color=_COLORS_N[n], ls=":", lw=1.0, alpha=0.7)
        ax[1].plot(_bg, resid, ls, color=_COLORS_N[n], lw=1.7,
                   label=f"N={n}/1　{FAR_LABEL[far]}")
        _rows9.append((FAR_LABEL[far], n, rho_bar,
                       float(np.nanmin(rho[_msk9])), float(np.nanmax(rho[_msk9])),
                       float(np.sqrt(np.nanmean(resid[_msk9] ** 2)))))
ax[0].axhline(1.0, color="0.4", ls=":", lw=1.2, label=r"$\rho$=1（多发不抬阈值）")
ax[0].set_xlabel("bg [计数/bin]"); ax[0].set_ylabel(r"阈值倍数 $\rho_N=T_N/T_1$")
ax[0].set_title("① 同 bg 下的阈值倍数（点线 = 各自均值）", fontsize=11)
if M9_YLIM_R is not None:
    ax[0].set_ylim(*M9_YLIM_R)
ax[1].axhline(0.0, color="k", ls=":", lw=1.2)
ax[1].set_xlabel("bg [计数/bin]"); ax[1].set_ylabel("相对残差 [%]")
ax[1].set_title(r"② 相对「$\rho$ = 常数」假设的残差", fontsize=11)
ax[1].set_ylim(*M9_YLIM_RES)
for a in ax:
    a.set_xlim(*M9_XLIM); a.grid(alpha=0.3); a.legend(fontsize=7.4)
fig.suptitle("模块 9　同 bg 下不同 N 的阈值倍数（数据来自模块 5，未重算）", fontsize=12.5)
fig.tight_layout(rect=[0, 0, 1, 0.92])
fig.savefig("pod_v30_m9_rho.png", dpi=120, bbox_inches="tight")
plt.show()

print("=" * 88)
print(f"阈值倍数汇总（只统计 bg ≥ {M9_BG_MIN:g} 的档，小 bg 处 T 量化噪声太大）")
print(f"{'FAR':>10}{'N/1':>6}{'ρ均值':>10}{'ρ最小':>9}{'ρ最大':>9}{'残差rms[%]':>12}")
for lab, n, rb, rlo, rhi, rms in _rows9:
    print(f"{lab:>10}{n:>6d}{rb:>10.4f}{rlo:>9.4f}{rhi:>9.4f}{rms:>12.2f}")
print("\n【怎么读】残差 rms 越小，说明「阈值差一个固定倍数」这个说法越站得住。")
print("  ρ 随 bg 上翘 = 大 bg 时多发累加要额外多抬阈值；根因见 theory_PoD_esti_v30.md。")
'''


# ============================================================ 模块 10：hist std / peak mu sigma

M10_MD = r"""
## 模块 10 — 单条 hist 内的 std、peak 均值、peak 标准差 随 bg

三个统计量，全部来自模块 5 的同一批 MC（每档 1e6 条），不重算：

| # | 量 | 定义 |
|---|---|---|
| ① | **单条 hist 内 std** | 一条 `hist_add` 在统计窗 152 个 bin 上的样本标准差，再对所有 MC 条取平均。注意是"**一条直方图内部**的起伏"，不是"同一个 bin 在多条直方图之间的起伏" |
| ② | **peak 均值** | 统计窗内最大 bin 的均值 |
| ③ | **peak 标准差** | 统计窗内最大 bin 跨 MC 条的标准差 |

①③ 各叠一条解析对照线，公式与推导见 `theory_PoD_esti_v30.md`。
"""

M10_PARAM = r'''
# ==================== 模块 10 绘图参数 ====================
M10_XLIM      = (0.0, 12.25)  # 【三张图横轴】bg 范围
M10_YLIM_STD  = None          # 【① 纵轴】hist 内 std；None = 自动
M10_YLIM_MEAN = None          # 【② 纵轴】peak 均值；None = 自动
M10_YLIM_PSTD = None          # 【③ 纵轴】peak 标准差；None = 自动
M10_SHOW_POISSON = True       # 【①】是否叠纯泊松 √bg 参考线
M10_SHOW_GUMBEL  = True       # 【③】是否叠 Gumbel 极值解析线
M10_FIGSIZE   = (16.5, 4.8)
'''

M10_PLOT = r'''
# ==================== 模块 10 绘图 ====================
_bg = np.asarray(BG_GRID, float)
fig, ax = plt.subplots(1, 3, figsize=M10_FIGSIZE)

for n in N_SHOTS_LIST:
    R = NOISE_RES[n]
    st = [peak_stats_from_cnt(c) for c in R["peak_cnt"]]
    pk_mean = np.array([s["mean"] for s in st])
    pk_std = np.array([s["std"] for s in st])
    sig_bin = np.sqrt(_bg * (1.0 - _bg / (N_PIX_MACRO * n)))   # 二项饱和解析

    ax[0].plot(_bg, R["hist_std"], "-", color=_COLORS_N[n], lw=2.0, label=f"N={n} MC 实测")
    ax[0].plot(_bg, sig_bin, "--", color=_COLORS_N[n], lw=1.1, alpha=0.85,
               label=f"N={n} 二项解析")
    ax[1].plot(_bg, pk_mean, "-", color=_COLORS_N[n], lw=2.0, label=f"N={n}")
    ax[2].plot(_bg, pk_std, "-", color=_COLORS_N[n], lw=2.0, label=f"N={n} MC 实测")
    if M10_SHOW_GUMBEL:
        z_M = (pk_mean - _bg) / np.maximum(sig_bin, 1e-9)
        sig_evt = (np.pi / np.sqrt(6.0)) * sig_bin / np.maximum(z_M, 1e-9)
        ax[2].plot(_bg, sig_evt, "--", color=_COLORS_N[n], lw=1.1, alpha=0.85,
                   label=f"N={n} Gumbel 解析")

if M10_SHOW_POISSON:
    ax[0].plot(_bg, np.sqrt(_bg), ":", color="k", lw=1.4, label="纯泊松 √bg（无饱和）")
ax[1].plot(_bg, _bg, ":", color="0.45", lw=1.3, label="参考 y=bg")

ax[0].set_ylabel("单条 hist 内 152 bin 的 std")
ax[0].set_title("① 每条直方图自身的起伏（亚泊松）", fontsize=11)
ax[1].set_ylabel("peak 均值 [计数/bin]")
ax[1].set_title("② peak 均值", fontsize=11)
ax[2].set_ylabel("peak 标准差 [计数/bin]")
ax[2].set_title("③ peak 标准差：实测 vs Gumbel 极值解析", fontsize=11)
for a, yl in zip(ax, [M10_YLIM_STD, M10_YLIM_MEAN, M10_YLIM_PSTD]):
    a.set_xlabel("bg [计数/bin]"); a.set_xlim(*M10_XLIM)
    if yl is not None:
        a.set_ylim(*yl)
    a.grid(alpha=0.3); a.legend(fontsize=7.0, ncol=2)
fig.suptitle(f"模块 10　三个统计量 vs bg（复用模块 5 的 {N_MC_NOISE:,} MC/档）", fontsize=12.5)
fig.tight_layout(rect=[0, 0, 1, 0.92])
fig.savefig("pod_v30_m10_stats.png", dpi=120, bbox_inches="tight")
plt.show()
'''


# ============================================================ 模块 11：有效 z 值

M11_MD = r"""
## 模块 11 — 有效 z 值：阈值离 peak 均值几个 peak 标准差

`z = (T − μ_peak) / σ_peak`，其中 `μ_peak` / `σ_peak` 是**纯噪声** peak 的均值与标准差
（都来自模块 5）。它回答"这条 FAR 的阈值，在噪声 peak 分布上站在第几个 σ"，
是把不同 bg、不同 N 拉到同一把尺子上比较阈值裕量的最直接方式。
"""

M11_PARAM = r'''
# ==================== 模块 11 绘图参数 ====================
M11_FARS   = [0.05, 0.01]    # 【画哪几条 FAR】
M11_XLIM   = (0.0, 12.25)    # 【横轴】bg 范围
M11_YLIM   = None            # 【纵轴】z 范围；None = 自动
M11_BG_MIN = 1.0             # 【统计均值时】忽略比它小的 bg
M11_FIGSIZE= (5.6 * 2, 4.6)
'''

M11_PLOT = r'''
# ==================== 模块 11 绘图 ====================
_bg = np.asarray(BG_GRID, float)
_msk11 = _bg >= M11_BG_MIN
fig, axes = plt.subplots(1, len(M11_FARS), figsize=M11_FIGSIZE, sharey=True)
axes = np.atleast_1d(axes)
for a, far in zip(axes, M11_FARS):
    tag = FAR_TAG[far]
    for n in N_SHOTS_LIST:
        R = NOISE_RES[n]
        st = [peak_stats_from_cnt(c) for c in R["peak_cnt"]]
        mu = np.array([s["mean"] for s in st])
        sd = np.array([s["std"] for s in st])
        z = (THRESH[n]["T" + tag].astype(float) - mu) / np.maximum(sd, 1e-9)
        zbar = float(np.nanmean(z[_msk11]))
        a.plot(_bg, z, "-", color=_COLORS_N[n], lw=1.4, alpha=0.75,
               label=f"N={n}（均值 {zbar:.2f}）")
        # 锯齿来自整数阈值（T 只能跳 1，bg 连续步进），均值线才是可引用的结论
        a.axhline(zbar, color=_COLORS_N[n], ls=":", lw=1.6)
    a.set_xlabel("bg [计数/bin]")
    a.set_title(f"FAR={FAR_LABEL[far]}", fontsize=11)
    a.set_xlim(*M11_XLIM)
    if M11_YLIM is not None:
        a.set_ylim(*M11_YLIM)
    a.grid(alpha=0.3); a.legend(fontsize=8)
axes[0].set_ylabel(r"$z=(T-\mu_{peak})/\sigma_{peak}$")
fig.suptitle("模块 11　有效 z 值：阈值站在噪声 peak 分布的第几个 σ", fontsize=12.5)
fig.tight_layout(rect=[0, 0, 1, 0.92])
fig.savefig("pod_v30_m11_zeff.png", dpi=120, bbox_inches="tight")
plt.show()
'''


# ============================================================ 模块 12：PoD50/90 所需信号

M12_MD = r"""
## 模块 12 — 各 FAR 下 PoD50 / PoD90 所需信号的均值

复用模块 6 的 `POD_RES`，不重算。"所需信号"给三种口径，避免歧义：

| 口径 | 定义 |
|---|---|
| `peak_mean` | 临界点上 `hist_add` 峰值的均值（含底噪） |
| `S_net = peak_mean − bg` | 净峰高 |
| `E_crit` | 临界发射脉冲能量 [nJ]，可直接对照激光器指标 |

> `S_net` 会**高估**信号贡献：即使没有信号，在信号窗里取最大值本身也已高于 bg（极值抬升）。
"""

M12_PARAM = r'''
# ==================== 模块 12 绘图参数 ====================
M12_FARS   = [0.05, 0.01]   # 【画哪几条 FAR】
M12_LEVELS = [0.50, 0.90]   # 【画哪几个 PoD 水平】
M12_XLIM   = (0.0, 12.25)   # 【横轴】bg 范围
M12_YLIM_S = None           # 【① 纵轴】净峰高；None = 自动
M12_YLIM_E = None           # 【② 纵轴】能量（对数）；None = 自动
M12_FIGSIZE= (5.8 * len(M12_FARS), 8.4)
M12_TABLE_FAR = 0.01        # 【文字表】用哪条 FAR
M12_TABLE_STEP = 8          # 【文字表】每隔几档 bg 打一行
'''

M12_PLOT = r'''
# ==================== 模块 12 绘图 ====================
_LS_LV = {0.50: "--", 0.90: "-"}
fig, axes = plt.subplots(2, len(M12_FARS), figsize=M12_FIGSIZE)
axes = np.atleast_2d(axes)
for j, far in enumerate(M12_FARS):
    tag = FAR_TAG[far]
    for n in N_SHOTS_LIST:
        for lv in M12_LEVELS:
            arr = _collect_critical(n, tag, lv)
            if not arr.size:
                continue
            lab = f"N={n} PoD{int(lv*100)}"
            axes[0, j].plot(arr[:, 0], arr[:, 2] - arr[:, 0], _LS_LV[lv],
                            color=_COLORS_N[n], lw=1.8, label=lab)
            axes[1, j].semilogy(arr[:, 0], arr[:, 3], _LS_LV[lv],
                                color=_COLORS_N[n], lw=1.8, label=lab)
    axes[0, j].set_title(f"净峰高 $S_{{net}}$　FAR={FAR_LABEL[far]}", fontsize=11)
    axes[0, j].set_ylabel("$S_{net}$ = peak − bg [计数]")
    axes[1, j].set_title(f"临界发射能量　FAR={FAR_LABEL[far]}", fontsize=11)
    axes[1, j].set_ylabel("$E_{crit}$ [nJ]")
    for r in (0, 1):
        axes[r, j].set_xlabel("bg [计数/bin]"); axes[r, j].set_xlim(*M12_XLIM)
        axes[r, j].grid(alpha=0.3, which="both"); axes[r, j].legend(fontsize=7.4, ncol=2)
if M12_YLIM_S is not None:
    for j in range(len(M12_FARS)):
        axes[0, j].set_ylim(*M12_YLIM_S)
if M12_YLIM_E is not None:
    for j in range(len(M12_FARS)):
        axes[1, j].set_ylim(*M12_YLIM_E)
fig.suptitle("模块 12　达到 PoD50 / PoD90 所需的信号（实线 PoD90，虚线 PoD50）", fontsize=12.5)
fig.tight_layout(rect=[0, 0, 1, 0.94])
fig.savefig("pod_v30_m12_signal_req.png", dpi=120, bbox_inches="tight")
plt.show()

print("=" * 96)
print(f"模块 12 抽样表（FAR={FAR_LABEL[M12_TABLE_FAR]}）")
print(f"{'N':>3}{'bg':>7}{'PoD':>6}{'T':>5}{'peakμ':>9}{'S_net':>9}{'E[nJ]':>11}")
_tag12 = FAR_TAG[M12_TABLE_FAR]
for n in N_SHOTS_LIST:
    for lv in M12_LEVELS:
        arr = _collect_critical(n, _tag12, lv)
        for i in range(0, len(arr), M12_TABLE_STEP):
            bgv, pod, pm, e_nj = arr[i, 0], arr[i, 1], arr[i, 2], arr[i, 3]
            Tv = int(arr[i, 5])
            print(f"{n:>3d}{bgv:>7.2f}{pod:>6.2f}{Tv:>5d}{pm:>9.2f}"
                  f"{pm-bgv:>9.2f}{e_nj:>11.3f}")
'''


# ============================================================ 模块 13：平方反比测远

M13_MD = r"""
## 模块 13 — 平方反比下的测距能力

信号回波功率 ∝ `1/D²`。把模块 6 的临界能量换算成"在给定发射能量下能测多远"：
给出**纯平方反比**与**含大气衰减** `exp(−2αD)` 两种口径。复用 `POD_RES`，不重算。
"""

M13_PARAM = r'''
# ==================== 模块 13 绘图参数 ====================
M13_FAR      = 0.01          # 【用哪条 FAR】
M13_LEVEL    = 0.90          # 【用哪个 PoD 水平】
M13_E_BUDGET = [1.0, 10.0, 100.0]   # 【发射能量预算 nJ】每个值一条曲线
M13_XLIM     = (0.0, 12.25)  # 【横轴】bg 范围
M13_YLIM     = None          # 【纵轴】距离 [m]；None = 自动
M13_FIGSIZE  = (12.6, 4.6)
'''

M13_PLOT = r'''
# ==================== 模块 13 绘图 ====================
# 记号：boost = 发射能量倍率；接收端衰减比 ratio = boost_crit / boost_avail，
#       满足 (D_ref/D)²·exp(−2α(D−D_ref)) = ratio。equiv_distance() 就是解这个方程。
_tag13 = FAR_TAG[M13_FAR]
fig, ax = plt.subplots(1, 2, figsize=M13_FIGSIZE)
for n in N_SHOTS_LIST:
    arr = _collect_critical(n, _tag13, M13_LEVEL)
    if not arr.size:
        continue
    bgv, boost_crit = arr[:, 0], arr[:, 6]
    for e_bud, ls in zip(M13_E_BUDGET, ["-", "--", ":"]):
        boost_avail = e_bud * 1e-9 / E_PULSE_BASE
        ratio = boost_crit / max(boost_avail, 1e-30)
        d_pure = D_TARGET / np.sqrt(np.maximum(ratio, 1e-30))   # 纯 1/D²
        d_att = np.array([equiv_distance(r) for r in ratio])    # 含 exp(−2αD)
        ax[0].plot(bgv, d_pure, ls, color=_COLORS_N[n], lw=1.7,
                   label=f"N={n} E={e_bud:g} nJ")
        ax[1].plot(bgv, d_att, ls, color=_COLORS_N[n], lw=1.7,
                   label=f"N={n} E={e_bud:g} nJ")
ax[0].set_title("① 纯平方反比 $1/D^2$", fontsize=11)
ax[1].set_title(f"② 含大气衰减 $e^{{-2\\alpha D}}$"
                f"（α={PARAMS['channel']['alpha']*1e3:.2f}/km）", fontsize=11)
for a in ax:
    a.set_xlabel("bg [计数/bin]"); a.set_ylabel("可测距离 [m]")
    a.set_xlim(*M13_XLIM)
    if M13_YLIM is not None:
        a.set_ylim(*M13_YLIM)
    a.grid(alpha=0.3); a.legend(fontsize=6.8, ncol=2)
fig.suptitle(f"模块 13　平方反比测距（FAR={FAR_LABEL[M13_FAR]}，PoD{int(M13_LEVEL*100)}）",
             fontsize=12.5)
fig.tight_layout(rect=[0, 0, 1, 0.92])
fig.savefig("pod_v30_m13_range.png", dpi=120, bbox_inches="tight")
plt.show()
'''


# ============================================================ 模块 14：同信号不同 bg 的 peak 分布

M14_MD = r"""
## 模块 14 — 同一信号强度、不同 bg 时 peak 分布怎么变

用模块 8 的固定信号扫描（每个 `(N, bg, boost)` 都存了**完整 peak 分布**），回答三件事：

1. **分布形状怎么变**：同一个 `boost` 在不同 bg 下的 peak 概率质量函数（PMF）叠一起看。
2. **peak 均值是不是"只是加上 bg"**：叠一条**理想加法线**
   `μ_ideal(bg) = μ(b, bg_min) + (bg − bg_min)`（斜率 = 1）。贴合则成立，偏离即证伪。
3. **peak 标准差怎么变**：对比有信号与无信号的 `σ_peak`，并看偏度。
"""

M14_PARAM = r'''
# ==================== 模块 14 绘图参数 ====================
M14_BG_SHOW  = [0.5, 3.0, 6.0, 9.0, 12.0]  # 【① PMF 叠加】画哪几个 bg 档
M14_BOOST_I  = None        # 【用第几档 boost】None = 取中间那档
M14_XLIM     = (0.0, 12.25)   # 【②③ 横轴】bg 范围
M14_PMF_XLIM = None        # 【① 横轴】peak 计数范围；None = 自动
M14_YLIM_MU  = None        # 【② 纵轴】peak 均值；None = 自动
M14_YLIM_SD  = None        # 【③ 纵轴】peak 标准差；None = 自动
M14_FIGSIZE  = (5.6 * len(N_SHOTS_LIST), 12.0)
'''

M14_PLOT = r'''
# ==================== 模块 14 绘图 ====================
_bg = np.asarray(BG_GRID, float)
_bo = np.asarray(BOOST_LIST_SIG, float)
_ib = int(len(_bo) // 2) if M14_BOOST_I is None else int(M14_BOOST_I)
_ks = [int(np.argmin(np.abs(_bg - b))) for b in M14_BG_SHOW]
_k0 = int(np.argmin(_bg))
_cbg = plt.cm.plasma(np.linspace(0.05, 0.85, len(_ks)))

fig, axes = plt.subplots(3, len(N_SHOTS_LIST), figsize=M14_FIGSIZE)
axes = np.atleast_2d(axes)
for j, n in enumerate(N_SHOTS_LIST):
    mu_b = SIG_RES[n]["mu"][_ib]
    sd_b = SIG_RES[n]["sd"][_ib]
    mu_0 = SIG_RES[n]["mu"][0]
    sd_0 = SIG_RES[n]["sd"][0]

    # ① 不同 bg 下的 peak PMF
    for c, k in zip(_cbg, _ks):
        cnt = np.asarray(SIG_RES[n]["cnt"][_ib, k], float)
        tot = max(cnt.sum(), 1.0)
        x = np.arange(cnt.size); m = cnt > 0
        axes[0, j].plot(x[m], (cnt / tot)[m], "-", color=c, lw=1.6,
                        label=f"bg={_bg[k]:g}")
    axes[0, j].set_title(f"N={n}　boost={_bo[_ib]:g}　peak 分布", fontsize=11)
    axes[0, j].set_xlabel("peak [计数]"); axes[0, j].set_ylabel("概率")
    if M14_PMF_XLIM is not None:
        axes[0, j].set_xlim(*M14_PMF_XLIM)
    axes[0, j].grid(alpha=0.3); axes[0, j].legend(fontsize=7.6)

    # ② 「peak = 信号峰 + bg」加法理想线
    ideal = mu_b[_k0] + (_bg - _bg[_k0])
    axes[1, j].plot(_bg, mu_b, "-", color=_COLORS_N[n], lw=2.2, label="实测 peak 均值")
    axes[1, j].plot(_bg, ideal, "k--", lw=1.9, label="理想：信号峰 + bg（斜率 1）")
    axes[1, j].plot(_bg, mu_0, ":", color="0.5", lw=1.5, label="无信号 peak 均值")
    axes[1, j].set_title(f"N={n}　加法假设检验", fontsize=11)
    axes[1, j].set_xlabel("bg [计数/bin]"); axes[1, j].set_ylabel("peak 均值 [计数]")
    axes[1, j].set_xlim(*M14_XLIM)
    if M14_YLIM_MU is not None:
        axes[1, j].set_ylim(*M14_YLIM_MU)
    axes[1, j].grid(alpha=0.3); axes[1, j].legend(fontsize=7.6)

    # ③ peak 标准差：有信号 vs 无信号
    axes[2, j].plot(_bg, sd_b, "-", color=_COLORS_N[n], lw=2.0,
                    label=f"有信号 b={_bo[_ib]:g}")
    axes[2, j].plot(_bg, sd_0, "--", color=_COLORS_N[n], lw=1.4, alpha=0.8,
                    label="无信号 b=0")
    axes[2, j].set_title(f"N={n}　peak 标准差", fontsize=11)
    axes[2, j].set_xlabel("bg [计数/bin]"); axes[2, j].set_ylabel(r"$\sigma_{peak}$ [计数]")
    axes[2, j].set_xlim(*M14_XLIM)
    if M14_YLIM_SD is not None:
        axes[2, j].set_ylim(*M14_YLIM_SD)
    axes[2, j].grid(alpha=0.3); axes[2, j].legend(fontsize=7.6)

fig.suptitle("模块 14　同一信号强度下，peak 分布随 bg 的变化", fontsize=12.5)
fig.tight_layout(rect=[0, 0, 1, 0.95])
fig.savefig("pod_v30_m14_peak_vs_bg.png", dpi=120, bbox_inches="tight")
plt.show()

print("=" * 92)
print(f"加法假设检验（boost={_bo[_ib]:g}）：理想线 = μ(b,bg_min)+(bg−bg_min)，斜率应为 1")
print(f"{'N':>3}{'bg范围':>16}{'实测Δpeak':>11}{'Δbg':>8}{'实测斜率':>10}{'高bg残差':>11}")
for n in N_SHOTS_LIST:
    mu_b = SIG_RES[n]["mu"][_ib]
    dbg = _bg[-1] - _bg[_k0]
    dpk = mu_b[-1] - mu_b[_k0]
    print(f"{n:>3d}{f'{_bg[_k0]:.2f}→{_bg[-1]:.2f}':>16}{dpk:>11.2f}{dbg:>8.2f}"
          f"{dpk/max(dbg,1e-9):>10.3f}{mu_b[-1]-(mu_b[_k0]+dbg):>+11.2f}")
'''


# ============================================================ 模块 15：宏像元对比

M15_MD = r"""
## 模块 15 — 宏像元 3×9 vs 3×6 阈值对比

比较两种宏像元尺寸：3×9 = 27 个 SPAD、3×6 = 18 个 SPAD。

**口径（不写清会得到相反结论）**：
- 噪声在每个 SPAD 上**均匀**，只由 `p_eq`（单 SPAD、单发、单 bin 被点亮的平衡态概率）刻画，
  与宏像元多大无关。宏像元每 bin 底噪 `bg = n_tr · p_eq`，`n_tr = n_pix × N_shots`。
- 信号也按每 SPAD **均匀**处理（用户指定；若按像斑加权结论会反过来）。
  灵敏度判据取 `q_req = (T − bg) / n_tr`，越小越灵敏。

**出图**（三张独立图，不再并排）：
- ①a 同 bg 阈值精简版：只画 3×6@N=2、3×9@N=2、3×9@N=4
- ①b 同 bg 阈值全量版：全部配置
- ② `q_req`：全部配置（内容与原先右图相同）
"""

M15_PARAM = r'''
# ==================== 模块 15 绘图参数 ====================
# 三张独立图，不再并排连图：
#   ①a 阈值精简版（只画 3×6@N=2、3×9@N=2、3×9@N=4）
#   ①b 阈值全量版（全部配置）
#   ②  q_req（全部配置，内容与原先右图相同）
M15_XLIM    = None       # 【横轴】三张图共用的 bg 范围；None = 自动（线性轴）
M15_YLIM_T  = None       # 【①a / ①b 纵轴】阈值 T
M15_YLIM_Q  = None       # 【② 纵轴】q_req
M15_ANNOT   = True       # 【②】是否标注等 n_tr 的两个配置重合（打印表里仍会校验）
M15_FIGSIZE = (7.6, 4.8) # 【单张图】尺寸；三张各自独立输出
# 【①a 精简版】只画这些 (宏像元名, N_shots)。名字必须与 MACRO_CFGS 里的 name 完全一致。
M15_T_FOCUS = [
    ("3×6（18 SPAD）", 2),
    ("3×9（27 SPAD）", 2),
    ("3×9（27 SPAD）", 4),
]
'''

M15_PLOT = r'''
# ==================== 模块 15 绘图 ====================
_MARKER = {"3×9（27 SPAD）": "o", "3×6（18 SPAD）": "s"}


def _m15_iter(keys=None):
    # keys=None → 全部配置；否则只迭代给定的 (name, N) 列表。
    # n_tr 相同的第二条改画虚线 + 空心点，让底下那条透出来，才能看出"重合"。
    drawn_ntr = set()
    for cfg in MACRO_CFGS:
        mk = _MARKER.get(cfg["name"], "o")
        for n in cfg["shots"]:
            if keys is not None and (cfg["name"], n) not in keys:
                continue
            d = MACRO_RES.get((cfg["name"], n))
            if not d:
                continue
            ntr = cfg["n_pix"] * n
            dup = ntr in drawn_ntr
            drawn_ntr.add(ntr)
            col = _MACRO_COLOR[(cfg["name"], n)]
            sty = dict(ls="--" if dup else "-", marker=mk, ms=5.0 if dup else 3.2,
                       color=col, lw=1.4 if dup else 1.6,
                       mfc="none" if dup else col, mew=1.2 if dup else 0.0)
            tail = "（与上面同 n_tr，故重合）" if dup else ""
            yield d, sty, f"{cfg['name']} N={n}（n_tr={ntr}）{tail}", cfg["name"], n


def _m15_apply_axes(a, ylabel, title, ylim):
    a.set_xlabel("bg [计数/bin]")
    a.set_ylabel(ylabel)
    a.set_title(title, fontsize=11)
    if M15_XLIM is not None:
        a.set_xlim(*M15_XLIM)
    if ylim is not None:
        a.set_ylim(*ylim)
    a.grid(alpha=0.3)
    a.legend(fontsize=7.6, loc="best")


_focus = set(M15_T_FOCUS)

# ---- ①a：阈值精简版 ----
fig, ax = plt.subplots(figsize=M15_FIGSIZE)
for d, sty, lab, _, _ in _m15_iter(_focus):
    ax.plot(d["bg"], d["T"], **sty, label=lab)
_m15_apply_axes(ax, "检测阈值 T [计数/bin]",
                "①a 同 bg 下的检测阈值（精简：3×6@N=2、3×9@N=2/4）", M15_YLIM_T)
fig.suptitle(f"模块 15　①a　宏像元阈值精简版（FAR={FAR_LABEL[MACRO_FAR]}，"
             f"每档 {MACRO_N_MC:,} MC）", fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.94])
fig.savefig("pod_v30_m15_T_focus.png", dpi=120, bbox_inches="tight")
plt.show()

# ---- ①b：阈值全量版 ----
fig, ax = plt.subplots(figsize=M15_FIGSIZE)
for d, sty, lab, _, _ in _m15_iter(None):
    ax.plot(d["bg"], d["T"], **sty, label=lab)
_m15_apply_axes(ax, "检测阈值 T [计数/bin]",
                "①b 同 bg 下的检测阈值（全部配置）", M15_YLIM_T)
fig.suptitle(f"模块 15　①b　宏像元阈值全量版（FAR={FAR_LABEL[MACRO_FAR]}，"
             f"每档 {MACRO_N_MC:,} MC）", fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.94])
fig.savefig("pod_v30_m15_T_all.png", dpi=120, bbox_inches="tight")
plt.show()

# ---- ②：q_req（内容与原先右图相同，单独成图）----
fig, ax = plt.subplots(figsize=M15_FIGSIZE)
for d, sty, lab, _, _ in _m15_iter(None):
    ax.plot(d["bg"], d["q_req"], **sty, label=lab)
_m15_apply_axes(ax, r"$q_{req}=(T-\mathrm{bg})/n_{tr}$",
                "② 每 SPAD 每发所需额外点亮概率（越小越灵敏）", M15_YLIM_Q)
fig.suptitle(f"模块 15　②　q_req（FAR={FAR_LABEL[MACRO_FAR]}，"
             f"每档 {MACRO_N_MC:,} MC）", fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.94])
fig.savefig("pod_v30_m15_qreq.png", dpi=120, bbox_inches="tight")
plt.show()

print("=" * 96)
print("等 n_tr 校验：n_tr 相同的配置，纯噪声阈值应当完全一致（噪声只经 n_tr 进入）")
_by_ntr = {}
for cfg in MACRO_CFGS:
    for n in cfg["shots"]:
        d = MACRO_RES.get((cfg["name"], n))
        if d:
            _by_ntr.setdefault(cfg["n_pix"] * n, []).append((cfg["name"], n, d["T"]))
for ntr, items in sorted(_by_ntr.items()):
    if len(items) < 2:
        continue
    base = items[0][2]
    same = all(np.array_equal(base, it[2]) for it in items[1:])
    names = "、".join(f"{a} N={b}" for a, b, _ in items)
    print(f"  n_tr={ntr:>4d}：{names}　→ 阈值{'完全相同' if same else '不同（需排查）'}")
'''


# ============================================================================
#                            计算 / 载入缓存 cell
# ============================================================================

# 追加到模块 5 计算段末尾的小工具（分位数由 bincount 直接算，供绘图 cell 用）
M5_HELPERS = r'''

def _quantile_from_cnt(cnt, p):
    """由 peak 的 bincount 求 p 分位（返回整数计数值）。"""
    tot = cnt.sum()
    if tot <= 0:
        return 0.0
    return float(np.searchsorted(np.cumsum(cnt) / tot, p))


print(f"模块 5 就绪：NOISE_RES / THRESH（{len(BG_GRID)} 档 bg × N={N_SHOTS_LIST}，"
      f"每档 {N_MC_NOISE:,} MC）。模块 9 / 10 / 11 直接复用，不再重算。")
'''


M7_COMPUTE = r'''
# ==================== 模块 7 计算（只做换算，不跑 MC）====================
def equiv_distance(boost, D_ref=D_TARGET, p=PARAMS):
    """把 boost 折算成发射能量和反射率不变时的等效距离。"""
    if not np.isfinite(boost) or boost <= 0:
        return np.nan
    alpha = p["channel"]["alpha"]
    Ds = np.logspace(np.log10(0.3), np.log10(5000.0), 6000)
    vals = (D_ref**2 / Ds**2) * np.exp(-2*alpha*(Ds-D_ref))
    if boost > vals[0] or boost < vals[-1]:
        return np.nan
    return float(np.interp(-boost, -vals, Ds))


def _collect_critical(n_shots, far_tag, level):
    """把 POD_RES 里某 (N, FAR, PoD水平) 的临界点按 bg 升序排成表。

    列：0=bg  1=实测PoD  2=peak均值  3=临界能量[nJ]  4=等效距离[m]  5=阈值T  6=boost
    """
    rows = []
    for nt in NOISE_GRID[n_shots]:
        r = POD_RES.get((n_shots, float(nt)))
        if not r or "critical" not in r:
            continue
        rec = r["critical"].get(far_tag, {}).get(f"{level:.2f}")
        if not rec:
            continue
        b = float(rec["boost"])
        rows.append((float(r["noise"]), float(rec["pod"]), float(rec["peak_mean"]),
                     b * E_PULSE_BASE * 1e9, equiv_distance(b),
                     float(r["T_map"][far_tag]), b))
    a = np.asarray(rows, float)
    return a[np.argsort(a[:, 0])] if a.size else a


print(f"模块 7 就绪：POD_RES 共 {len(POD_RES)} 个 (N, bg) 档；"
      f"临界能量只解 FAR={[FAR_LABEL[f] for f in POD_FARS]}")
'''


M8_MD = r"""
## 模块 8 — 固定信号强度 × 全 bg 网格

对每个 `(N_shots, bg, boost)` 组合跑 MC，**保存完整的 peak 分布**（bincount）。
`boost` 固定不变、只让 bg 从 0.25 涨到 12，用来回答"同样的信号，底噪涨上来之后还剩多少"。

结果 `SIG_RES` 也是模块 14 的数据源，两个模块共用一份缓存。

- **boost**：回波能量倍率，`boost=0` 表示纯噪声（作为基线）。
- **Δpeak**：`peak_mean(boost) − peak_mean(0)`，信号带来的净抬升。
"""

M8_COMPUTE = r'''
# ==================== 模块 8 计算（多进程脚本 + 缓存）====================
BOOST_LIST_SIG = np.round(np.arange(0.0, 0.032 + 1e-12, 0.004), 6)
SIG_SCRIPT = "run_pod_v30_sig_scan.py"


def _try_load_sig(path):
    if not (USE_CACHE and os.path.exists(path)):
        return None
    try:
        z = np.load(path)
        if (int(z["n_mc"]) != N_MC_SIG
                or not np.allclose(z["grid_key"], BG_GRID)
                or not np.array_equal(z["n_shots_list"], np.asarray(N_SHOTS_LIST))
                or not np.allclose(z["boosts"], BOOST_LIST_SIG)):
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
    _zsig = _try_load_sig(_cand)
    if _zsig is not None:
        print(f"模块 8 命中缓存 {_cand}")
        break

if _zsig is None:
    import sys
    print("=" * 72)
    print("未找到完整 v30 信号缓存 → 自动调用多进程扫描（ProcessPool，吃满 CPU）")
    print(f"规模：{len(BG_GRID)} bg × {len(BOOST_LIST_SIG)} boost × N={N_SHOTS_LIST}"
          f" × {N_MC_SIG:,} MC")
    print("=" * 72)
    _rc = _run_cmd_stream([sys.executable, SIG_SCRIPT,
                           "--workers", str(int(N_WORKERS)),
                           "--n-mc", str(int(N_MC_SIG))])
    if _rc != 0:
        raise RuntimeError(f"{SIG_SCRIPT} 失败，请查看上方进度输出")
    _zsig = _try_load_sig(CACHE_SIG)
    if _zsig is None:
        raise RuntimeError(f"多进程信号扫描结束但仍无法载入完整缓存 {CACHE_SIG}")
    print(f"多进程信号扫描完成，已载入 {CACHE_SIG}")

SIG_RES = {}
for n in N_SHOTS_LIST:
    cnt = np.asarray(_zsig[f"peak_cnt_{n}"])
    mu = np.zeros((len(BOOST_LIST_SIG), len(BG_GRID)))
    sd = np.zeros_like(mu)
    for i in range(len(BOOST_LIST_SIG)):
        for k in range(len(BG_GRID)):
            s = peak_stats_from_cnt(cnt[i, k])
            mu[i, k] = s["mean"]; sd[i, k] = s["std"]
    SIG_RES[n] = dict(cnt=cnt, mu=mu, sd=sd)
print(f"模块 8 就绪：{len(BOOST_LIST_SIG)} boost × {len(BG_GRID)} bg × N={N_SHOTS_LIST}，"
      f"每档 {N_MC_SIG:,} MC（含完整 peak 分布）。模块 14 直接复用。")
'''

M8_PARAM = r'''
# ==================== 模块 8 绘图参数 ====================
M8_XLIM      = (0.0, 12.25)  # 【两排横轴】bg 范围
M8_YLIM_MU   = None          # 【上排纵轴】peak 均值；None = 自动
M8_YLIM_GAIN = None          # 【下排纵轴】归一化增益 Δpeak/boost；None = 自动
M8_CMAP      = "viridis"     # 【配色】boost 由小到大
M8_FARS      = [0.01, 0.05]  # 【上排叠加的阈值】只画 1% 与 5%（来自模块 5 的 THRESH）
M8_FIGSIZE   = (5.4 * len(N_SHOTS_LIST), 8.2)
'''

M8_PLOT = r'''
# ==================== 模块 8 绘图 ====================
_bg = np.asarray(BG_GRID, float)
_bo = np.asarray(BOOST_LIST_SIG, float)
_cb = plt.get_cmap(M8_CMAP)(np.linspace(0.08, 0.92, len(_bo)))
fig, axes = plt.subplots(2, len(N_SHOTS_LIST), figsize=M8_FIGSIZE, sharex=True)
axes = np.atleast_2d(axes)
for j, n in enumerate(N_SHOTS_LIST):
    mu = SIG_RES[n]["mu"]
    for i, b in enumerate(_bo):
        axes[0, j].plot(_bg, mu[i], "-", color=_cb[i], lw=1.4,
                        label=f"boost={b:g}" + ("（纯噪声）" if b == 0 else ""))
        if b > 0:   # 归一化增益：同样一份信号在不同 bg 下还剩多少
            axes[1, j].plot(_bg, (mu[i] - mu[0]) / b, "-", color=_cb[i], lw=1.4,
                            label=f"boost={b:g}")
    # 叠上 1% / 5% 检测阈值：peak 均值曲线跌到阈值以下，就是这份信号在该 bg 下失守的位置
    for far, ls in zip(M8_FARS, ["--", "-."]):
        axes[0, j].plot(_bg, THRESH[n]["T" + FAR_TAG[far]], ls, color="k", lw=1.5,
                        label=f"阈值 {FAR_LABEL[far]}")
    axes[0, j].set_title(f"N_shots={n}　peak 均值 vs 检测阈值", fontsize=11)
    axes[0, j].set_ylabel("peak 均值 [计数]")
    axes[1, j].set_title(f"N_shots={n}　归一化信号增益 Δpeak / boost", fontsize=11)
    axes[1, j].set_xlabel("bg [计数/bin]")
    axes[1, j].set_ylabel(r"$\Delta$peak / boost")
    for r in (0, 1):
        axes[r, j].set_xlim(*M8_XLIM); axes[r, j].grid(alpha=0.3)
    if M8_YLIM_MU is not None:
        axes[0, j].set_ylim(*M8_YLIM_MU)
    if M8_YLIM_GAIN is not None:
        axes[1, j].set_ylim(*M8_YLIM_GAIN)
axes[0, 0].legend(fontsize=6.8, ncol=2)
axes[1, 0].legend(fontsize=6.8, ncol=2)
fig.suptitle("模块 8　固定信号强度下，peak 与信号增益随 bg 的变化\n"
             "下排若随 bg 下滑，说明高底噪把同一份信号「吃掉」了（1-bit 抢占）",
             fontsize=12.2)
fig.tight_layout(rect=[0, 0, 1, 0.93])
fig.savefig("pod_v30_m8_sig_vs_bg.png", dpi=120, bbox_inches="tight")
plt.show()

print("=" * 84)
print("归一化增益随 bg 的衰减（越小说明信号被底噪吃得越狠）")
print(f"{'N':>3}{'boost':>9}{'bg=0.25':>10}{'bg=6':>9}{'bg=12':>9}{'12档/最低档':>13}")
_k_lo, _k_mid, _k_hi = 0, int(np.argmin(np.abs(_bg - 6.0))), len(_bg) - 1
for n in N_SHOTS_LIST:
    mu = SIG_RES[n]["mu"]
    for i, b in enumerate(_bo):
        if b == 0:
            continue
        g = (mu[i] - mu[0]) / b
        print(f"{n:>3d}{b:>9.3f}{g[_k_lo]:>10.2f}{g[_k_mid]:>9.2f}{g[_k_hi]:>9.2f}"
              f"{g[_k_hi]/max(g[_k_lo],1e-9):>13.3f}")

print("\n" + "=" * 84)
print("固定信号的「失守 bg」：peak 均值跌破阈值时的 bg（—— 表示全程都在阈值之上）")
print(f"{'N':>3}{'boost':>9}" + "".join(f"{'跌破'+FAR_LABEL[f]:>12}" for f in M8_FARS))
for n in N_SHOTS_LIST:
    mu = SIG_RES[n]["mu"]
    for i, b in enumerate(_bo):
        if b == 0:
            continue
        row = f"{n:>3d}{b:>9.3f}"
        for far in M8_FARS:
            d = mu[i] - THRESH[n]["T" + FAR_TAG[far]]
            k = np.where(d < 0)[0]
            row += f"{(f'{_bg[k[0]]:.2f}' if k.size else '——'):>12}"
        print(row)
'''


M15_COMPUTE = r'''
# ==================== 模块 15 计算（多进程脚本 + 缓存）====================
MACRO_SCRIPT = "compare_macro_v30.py"
MACRO_CACHE = "compare_macro_v30_cache.npz"
MACRO_FAR = 0.01                 # 主 FAR
MACRO_FAR_KEYS = [0.05, 0.01, 0.001, 100e-6]
MACRO_N_MC = 600_000             # ★ v30：由 v20 的 300k 提到 600k
MACRO_CFGS = [
    # shots 里给 3×6 多留一个 N=6：它的 n_tr = 18×6 = 108，与 3×9 N=4 相同，
    # 是「纯噪声只经 n_tr 进入」这条理论的等 n_tr 校验点。
    {"name": "3×9（27 SPAD）", "nx": 9, "ny": 3, "n_pix": 27, "shots": [1, 2, 4]},
    {"name": "3×6（18 SPAD）", "nx": 6, "ny": 3, "n_pix": 18, "shots": [1, 2, 4, 6]},
]
_MACRO_COLOR = {
    ("3×9（27 SPAD）", 1): "#9ecae1", ("3×9（27 SPAD）", 2): "#4292c6",
    ("3×9（27 SPAD）", 4): "#08519c",
    ("3×6（18 SPAD）", 1): "#fcae91", ("3×6（18 SPAD）", 2): "#fb6a4a",
    ("3×6（18 SPAD）", 4): "#cb181d", ("3×6（18 SPAD）", 6): "#67000d",
}

if not os.path.exists(MACRO_CACHE):
    import sys
    print("=" * 72)
    print(f"未找到 {MACRO_CACHE} → 调用 {MACRO_SCRIPT}（{MACRO_N_MC:,} MC/档，耗时较长）")
    print("=" * 72)
    _rc = _run_cmd_stream([sys.executable, MACRO_SCRIPT,
                           "--workers", str(int(N_WORKERS)),
                           "--n-mc", str(int(MACRO_N_MC))])
    if _rc != 0:
        raise RuntimeError(f"{MACRO_SCRIPT} 失败，请查看上方进度输出")

_zm = np.load(MACRO_CACHE)
_cfgs_raw = [tuple(int(v) for v in c) for c in np.asarray(_zm["cfgs"])]
_i_far = MACRO_FAR_KEYS.index(MACRO_FAR)
MACRO_RES = {}
for cfg in MACRO_CFGS:
    for n in cfg["shots"]:
        key = (cfg["nx"], cfg["ny"], n)
        if key not in _cfgs_raw:
            continue
        ci = _cfgs_raw.index(key)
        ok = np.asarray(_zm["done"])[ci]
        bg = np.asarray(_zm["bg_mc"])[ci][ok]
        T = np.asarray(_zm[f"thr_{_i_far}"])[ci][ok]
        o = np.argsort(bg)
        bg, T = bg[o], T[o]
        n_tr = cfg["n_pix"] * n
        MACRO_RES[(cfg["name"], n)] = dict(
            bg=bg, T=T, n_tr=n_tr, q_req=(T - bg) / n_tr,
            peak_mean=np.asarray(_zm["peak_mean"])[ci][ok][o],
        )
print(f"模块 15 就绪：{len(MACRO_RES)} 条曲线，FAR={MACRO_FAR:.0%}，"
      f"每档 {MACRO_N_MC:,} MC，来自 {MACRO_CACHE}")
'''
