# -*- coding: utf-8 -*-
"""增强 PoD_esti_v04.ipynb 的信号扫描与汇总模块。

只修改 notebook 源码，不执行任何计算。
"""

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
NOTEBOOK = ROOT / "PoD_esti_v04.ipynb"


MODULE7_MD = r"""## 模块 7 —— 第 3 步：逐噪声档精确求 PoD 50% / 90% 临界能量

v02 的 24 个全局对数能量点只适合展示趋势，交点主要依赖两个相距较远的点插值；
本模块改为对 **NOISE_GRID 中每个 0.25-noise 档**自适应求交点：

1. 先用少量粗扫描点包住 PoD 从 0 到 1 的过渡区；
2. 在过渡区补局部能量点，并用 probit 拟合求 PoD 50% / 90% 的候选 boost；
3. 在每个候选临界能量上重新做独立 Monte Carlo 验证；
4. 若验证 PoD 与目标偏差较大，根据 probit 斜率修正一次能量并重新验证。

每个临界点保存：实际 PoD、peak 均值、标准差和完整整数分布 `peak_cnt`。
因此模块 8 可以画出每个 noise 档的临界能量、临界 peak 均值、距离；本模块也会直接展示代表性
noise 档在临界能量处的 peak 分布，并用阈值竖线说明为什么能达到相应 PoD。

PoD 判据仍为：**信号窗内 peak ≥ T**。物理参数、阈值定义和信号窗均未修改。
"""


MODULE7_CODE = r'''# ---- PoD 专用子窗（只计算信号附近，前方保留暖机）----
from concurrent.futures import ThreadPoolExecutor, as_completed

POD_T_LO = T0_SIG - POD_WARM_NS * 1e-9
POD_T_HI = T0_SIG + SIG_POST_NS * 1e-9
TF_POD = np.arange(POD_T_LO, POD_T_HI, DT_FINE)
_sigmask = (TC_NS >= T0_SIG_NS - SIG_PRE_NS) & (TC_NS <= T0_SIG_NS + SIG_POST_NS)
IDX_SIG = np.where(_sigmask)[0]
CENTERS_SIG = CENTERS[IDX_SIG]
R_SIG_UNIT_POD = signal_photon_rate_fine(ECHO0, 1.0, TF_POD)
_NPH_BASE = np.trapezoid(R_SIG_UNIT_POD, TF_POD) * F_VALS.sum()

print(f"PoD 子窗：{POD_T_LO*1e9:.1f}–{POD_T_HI*1e9:.1f} ns，{TF_POD.size} 个细网格步")
print(f"每种 N_shots 对自己的完整 NOISE_GRID 求解："
      f"{[len(NOISE_GRID[n]) for n in N_SHOTS_LIST]} 档，noise 步长 0.25")
print(f"并行线程数 POD_WORKERS={POD_WORKERS}；临界点独立验证 {N_MC_POD_VERIFY:,} 次")


def sig_peaks(boost, n_shots, r_amb, n_real, seed):
    """返回给定 boost 下 n_real 次实现的信号窗 peak。"""
    f_arr = np.tile(F_VALS, n_shots)
    h = binary_macro_stepping(
        n_real, f_arr, R_SIG_UNIT_POD, TF_POD, r_amb, CENTERS_SIG,
        np.random.default_rng(seed), boost=boost,
    )
    return h.max(axis=1)


def _eval_boost_grid(boosts, n_shots, r_amb, n_real, seed0):
    """并行评估若干独立能量点；返回每点的 peak 样本。"""
    boosts = np.asarray(boosts, float)
    out = [None] * len(boosts)
    if POD_WORKERS <= 1:
        for i, boost in enumerate(boosts):
            out[i] = sig_peaks(boost, n_shots, r_amb, n_real, seed0 + 1009*i)
        return out
    with ThreadPoolExecutor(max_workers=POD_WORKERS) as pool:
        jobs = {
            pool.submit(sig_peaks, boost, n_shots, r_amb, n_real, seed0 + 1009*i): i
            for i, boost in enumerate(boosts)
        }
        for job in as_completed(jobs):
            out[jobs[job]] = job.result()
    return out


def _isotonic(p):
    """简单保序：消除有限 MC 导致的局部 PoD 下降。"""
    return np.maximum.accumulate(np.asarray(p, float))


def _crossing_logboost(boosts, pod, level):
    """在 log10(boost) 上找首次跨越；越界返回 NaN。"""
    order = np.argsort(boosts)
    x = np.log10(np.asarray(boosts)[order])
    p = _isotonic(np.asarray(pod)[order])
    if p[0] >= level:
        return float(x[0])
    if p[-1] < level:
        return np.nan
    i = int(np.searchsorted(p, level))
    dp = p[i] - p[i-1]
    w = 0.5 if dp <= 0 else (level - p[i-1]) / dp
    return float(x[i-1] + w*(x[i] - x[i-1]))


def _probit_fit(boosts, pod, n_real):
    """拟合 Phi^-1(PoD) = slope*log10(boost) + intercept。"""
    boosts = np.asarray(boosts, float)
    # Jeffreys 平滑，避免 0/1 在 norm.ppf 后成为无穷。
    success = np.rint(np.asarray(pod, float) * n_real)
    p = (success + 0.5) / (n_real + 1.0)
    transition = (p > 0.01) & (p < 0.99)
    if transition.sum() < 3:
        transition = np.argsort(np.abs(p - 0.5))[:min(5, len(p))]
    x = np.log10(boosts[transition])
    z = _norm.ppf(p[transition])
    slope, intercept = np.polyfit(x, z, 1)
    return float(slope), float(intercept)


def _critical_record(boost, n_shots, r_amb, T, level, slope, seed):
    """在候选临界能量独立验证；必要时按 probit 斜率修正一次。"""
    final = None
    for correction in range(2):
        pk = sig_peaks(boost, n_shots, r_amb, N_MC_POD_VERIFY, seed + 7919*correction)
        pod_actual = float((pk >= T).mean())
        final = {
            "boost": float(boost),
            "pod": pod_actual,
            "peak_mean": float(pk.mean()),
            "peak_std": float(pk.std()),
            "peak_cnt": np.bincount(pk, minlength=N_PIX_MACRO*n_shots + 1),
            "n_verify": int(pk.size),
        }
        if correction == 1 or abs(pod_actual - level) <= POD_VERIFY_TOL or slope <= 0:
            break
        p_smooth = (int((pk >= T).sum()) + 0.5) / (pk.size + 1.0)
        dx = (_norm.ppf(level) - _norm.ppf(p_smooth)) / slope
        # 单次修正最多 0.25 decade，避免有限 MC 偶然波动造成过冲。
        boost *= 10.0 ** float(np.clip(dx, -0.25, 0.25))
    return final


def solve_pod_noise(n_shots, k, seed0):
    """求一个 noise 档、两条 FAR 阈值下的 PoD50/90 临界点。"""
    R, Tr = NOISE_RES[n_shots], THRESH[n_shots]
    nt = float(R["noise_target"][k])
    n_tr = int(R["n_tr"])
    r_amb = float(R["r_det"][k] / PDE)
    T_map = {f"{far*1e6:.0f}": int(Tr[f"T{far*1e6:.0f}"][k]) for far in TARGET_FARS}
    if max(T_map.values()) > n_tr:
        return (n_shots, nt), {
            "noise": float(R["noise_mc"][k]), "e_lambda": float(R["e_lambda"][k]),
            "n_tr": n_tr, "T_map": T_map, "critical": {}, "invalid": "阈值超过二值硬上限",
        }

    # 全局粗扫只负责可靠包住过渡区，不直接作为最终交点。
    coarse_boost = np.logspace(POD_LOG_BOOST_MIN, POD_LOG_BOOST_MAX, N_POD_COARSE)
    coarse_pk = _eval_boost_grid(
        coarse_boost, n_shots, r_amb, N_MC_POD_COARSE, seed0,
    )
    coarse_pod = {
        tag: np.array([(pk >= T).mean() for pk in coarse_pk])
        for tag, T in T_map.items()
    }

    # 围绕四个粗交点补局部点；所有阈值共享同一批 peak 样本。
    roots0 = []
    for tag in T_map:
        for level in POD_LEVELS:
            x0 = _crossing_logboost(coarse_boost, coarse_pod[tag], level)
            if np.isfinite(x0):
                roots0.append(x0)
    if roots0:
        local_x = np.unique(np.concatenate([
            np.linspace(x0 - POD_LOCAL_HALF_DECADE, x0 + POD_LOCAL_HALF_DECADE,
                        N_POD_LOCAL_PER_ROOT)
            for x0 in roots0
        ]))
        local_boost = 10.0**local_x
        local_pk = _eval_boost_grid(
            local_boost, n_shots, r_amb, N_MC_POD_LOCAL, seed0 + 500_000,
        )
    else:
        local_boost = np.array([], float)
        local_pk = []

    critical = {}
    curve = {}
    for tag, T in T_map.items():
        boosts_fit = np.concatenate([coarse_boost, local_boost])
        pod_fit = np.concatenate([
            coarse_pod[tag],
            np.array([(pk >= T).mean() for pk in local_pk]),
        ])
        order = np.argsort(boosts_fit)
        boosts_fit, pod_fit = boosts_fit[order], pod_fit[order]
        curve[tag] = {"boost": boosts_fit, "pod": pod_fit}
        slope, intercept = _probit_fit(boosts_fit, pod_fit, N_MC_POD_LOCAL)
        critical[tag] = {}
        for level in POD_LEVELS:
            x_root = (_norm.ppf(level) - intercept) / slope if slope > 0 else np.nan
            if not np.isfinite(x_root):
                critical[tag][f"{level:.2f}"] = None
                continue
            boost = float(10.0**x_root)
            critical[tag][f"{level:.2f}"] = _critical_record(
                boost, n_shots, r_amb, T, level, slope,
                seed0 + int(level*10000) + int(tag)*101,
            )

    return (n_shots, nt), {
        "noise": float(R["noise_mc"][k]),
        "noise_target": nt,
        "e_lambda": float(R["e_lambda"][k]),
        "n_tr": n_tr,
        "T_map": T_map,
        "curve": curve,
        "critical": critical,
    }


# ---- 对完整 0.25-noise 网格求解；缓存键包含全部精度参数 ----
_pod_grid_key = np.concatenate([NOISE_GRID[n] for n in N_SHOTS_LIST])
_need = True
if USE_CACHE and os.path.exists(CACHE_POD):
    _z = np.load(CACHE_POD, allow_pickle=True)
    if (
        np.array_equal(_z["grid_key"], _pod_grid_key)
        and int(_z["n_coarse"]) == N_MC_POD_COARSE
        and int(_z["n_local"]) == N_MC_POD_LOCAL
        and int(_z["n_verify"]) == N_MC_POD_VERIFY
    ):
        POD_RES = _z["res"].item()
        _need = False
        print(f"已从缓存 {CACHE_POD} 载入逐 noise PoD 临界点")

if _need:
    POD_RES = {}
    _tall = time.time()
    n_total = sum(len(NOISE_GRID[n]) for n in N_SHOTS_LIST)
    print(f"PoD 临界点扫描：共 {n_total} 个 noise 档，步长 0.25；"
          f"每档同时求 100/10 ppm × PoD50/90")
    for n_shots in N_SHOTS_LIST:
        jobs = [(k, 7000 + n_shots*1_000_000 + k*20_000)
                for k in range(len(NOISE_GRID[n_shots]))]
        # 外层保持顺序，能量点在 solve_pod_noise 内部并行，避免线程嵌套。
        for done, (k, seed) in enumerate(jobs, 1):
            key, value = solve_pod_noise(n_shots, k, seed)
            POD_RES[key] = value
            if done == 1 or done % 5 == 0 or done == len(jobs):
                c100 = value.get("critical", {}).get("100", {})
                p90 = c100.get("0.90")
                msg = "无有效交点" if not p90 else (
                    f"E90={p90['boost']*E_PULSE_BASE*1e9:.3g} nJ，"
                    f"验证PoD={p90['pod']:.3f}，peak均值={p90['peak_mean']:.2f}"
                )
                elapsed = time.time() - _tall
                print(f"  [N_shots={n_shots} {done}/{len(jobs)}] "
                      f"noise={key[1]:.2f}：{msg}；累计 {elapsed/60:.1f} min")
    np.savez_compressed(
        CACHE_POD, res=np.array(POD_RES, dtype=object),
        grid_key=_pod_grid_key,
        n_coarse=N_MC_POD_COARSE, n_local=N_MC_POD_LOCAL,
        n_verify=N_MC_POD_VERIFY,
    )
    print(f"逐 noise PoD 扫描总用时 {(time.time()-_tall)/60:.1f} min；已写入 {CACHE_POD}")
'''


MODULE7_PLOT = r'''# ---- 模块 7 验证图：代表性噪声档的 PoD 曲线与临界 peak 分布 ----
fig, axes = plt.subplots(2, len(N_SHOTS_LIST), figsize=(8.2*len(N_SHOTS_LIST), 10.5))
if len(N_SHOTS_LIST) == 1:
    axes = axes.reshape(2, 1)

for col, n_shots in enumerate(N_SHOTS_LIST):
    ax_pod, ax_dist = axes[0, col], axes[1, col]
    grid = NOISE_GRID[n_shots]
    selected = grid[np.unique(np.linspace(0, len(grid)-1, 5).astype(int))]
    colors = plt.cm.viridis(np.linspace(0.08, 0.92, len(selected)))

    for color, nt in zip(colors, selected):
        r = POD_RES.get((n_shots, float(nt)))
        if not r or "curve" not in r:
            continue
        tag = "100"
        curve = r["curve"][tag]
        energy = curve["boost"] * E_PULSE_BASE * 1e9
        ax_pod.semilogx(energy, curve["pod"], "o-", ms=3, lw=1.2, color=color,
                        label=f"noise={nt:.2f}, T={r['T_map'][tag]}")
        for level, marker in [(0.50, "s"), (0.90, "*")]:
            rec = r["critical"][tag][f"{level:.2f}"]
            if rec:
                ax_pod.plot(rec["boost"]*E_PULSE_BASE*1e9, rec["pod"], marker,
                            color=color, ms=11 if level == 0.90 else 7,
                            mec="k", mew=0.4)

    ax_pod.axhline(0.50, color="0.35", ls="--", lw=1, label="目标 PoD 50%")
    ax_pod.axhline(0.90, color="0.35", ls=":", lw=1, label="目标 PoD 90%")
    ax_pod.set_xlabel("等效单脉冲发射能量 [nJ]（对数轴）")
    ax_pod.set_ylabel("验证/扫描 PoD")
    ax_pod.set_ylim(-0.03, 1.03)
    ax_pod.set_title(f"N_shots={n_shots}：局部加密后的交点（100 ppm 阈值）")
    ax_pod.legend(fontsize=7.5, ncol=2); ax_pod.grid(alpha=0.3, which="both")

    # 只展示代表性 noise 的 PoD90 临界分布；所有 noise 的 peak_cnt 均已存入缓存。
    for color, nt in zip(colors, selected):
        r = POD_RES.get((n_shots, float(nt)))
        if not r or "critical" not in r:
            continue
        rec = r["critical"].get("100", {}).get("0.90")
        if not rec:
            continue
        cnt = np.asarray(rec["peak_cnt"])
        x = np.arange(cnt.size)
        ax_dist.step(x, cnt/cnt.sum(), where="mid", color=color, lw=1.5,
                     label=(f"noise={nt:.2f}: mean={rec['peak_mean']:.2f}, "
                            f"PoD={rec['pod']:.3f}"))
        ax_dist.axvline(r["T_map"]["100"], color=color, ls=":", lw=1)
    ax_dist.set_xlabel("临界能量下的信号窗 peak [计数 / 1 ns bin]")
    ax_dist.set_ylabel("概率")
    ax_dist.set_title("PoD90 临界能量处的 peak 分布（同色虚线=T）")
    ax_dist.legend(fontsize=7.3); ax_dist.grid(alpha=0.3)

plt.suptitle(
    f"模块 7　临界能量交点及独立验证（每个临界点 {N_MC_POD_VERIFY:,} 次 MC，滤前）",
    fontsize=12.5,
)
plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.savefig("pod_v04_pod_critical_validation.png", dpi=120, bbox_inches="tight")
plt.show()
'''


MODULE8_MD = r"""## 模块 8 — 全部 0.25-noise 档的临界 peak、能量与距离

横轴直接使用模块 5 的完整 `NOISE_GRID`：
- N_shots=1：0.25–12，步长 0.25；
- N_shots=4：0.25–40，步长 0.25。

每个点都来自模块 7 在对应 noise 下的独立临界点验证，不再只计算 `[1, 5, 10]` 三档。
"""


MODULE8_CODE = r'''def equiv_distance(boost, D_ref=D_TARGET, p=PARAMS):
    """把 boost 折算成发射能量和反射率不变时的等效距离。"""
    if not np.isfinite(boost) or boost <= 0:
        return np.nan
    alpha = p["channel"]["alpha"]
    Ds = np.logspace(np.log10(0.3), np.log10(5000.0), 6000)
    vals = (D_ref**2 / Ds**2) * np.exp(-2*alpha*(Ds-D_ref))
    if boost > vals[0] or boost < vals[-1]:
        return np.nan
    return float(np.interp(-boost, -vals, Ds))


def collect_critical(n_shots, far_tag, level):
    rows = []
    for nt in NOISE_GRID[n_shots]:
        r = POD_RES.get((n_shots, float(nt)))
        if not r or "critical" not in r:
            continue
        rec = r["critical"].get(far_tag, {}).get(f"{level:.2f}")
        if not rec:
            continue
        rows.append((
            r["noise"], rec["boost"], rec["pod"], rec["peak_mean"],
            rec["peak_std"], r["T_map"][far_tag], equiv_distance(rec["boost"]),
        ))
    return np.asarray(rows, float)


fig, ax = plt.subplots(1, 3, figsize=(19, 5.7))
colors = {1: "tab:blue", 4: "tab:red"}
styles = {"100": "-", "10": "--"}

# ① 临界发射能量
for n_shots in N_SHOTS_LIST:
    for tag in ["100", "10"]:
        for level, alpha, width in [(0.50, 0.65, 1.4), (0.90, 1.0, 2.1)]:
            a = collect_critical(n_shots, tag, level)
            if a.size:
                ax[0].semilogy(
                    a[:, 0], a[:, 1]*E_PULSE_BASE*1e9,
                    styles[tag], color=colors[n_shots], alpha=alpha, lw=width,
                    label=f"N={n_shots}, {tag}ppm, PoD{level:.0%}",
                )
ax[0].set_xlabel("噪声均值 noise [计数 / 1 ns bin]（步长 0.25）")
ax[0].set_ylabel("临界等效单脉冲发射能量 [nJ]（对数轴）")
ax[0].set_title("① 每个 noise 的 PoD50/90 临界能量")
ax[0].legend(fontsize=7.2, ncol=2); ax[0].grid(alpha=0.3, which="both")

# ② 临界能量处的 peak 均值，并与阈值对照
for n_shots in N_SHOTS_LIST:
    for tag in ["100", "10"]:
        for level, alpha, width in [(0.50, 0.65, 1.4), (0.90, 1.0, 2.1)]:
            a = collect_critical(n_shots, tag, level)
            if a.size:
                ax[1].plot(
                    a[:, 0], a[:, 3], styles[tag], color=colors[n_shots],
                    alpha=alpha, lw=width,
                    label=f"N={n_shots}, {tag}ppm, PoD{level:.0%}",
                )
        a = collect_critical(n_shots, tag, 0.90)
        if a.size:
            ax[1].plot(a[:, 0], a[:, 5], ":", color=colors[n_shots], alpha=0.45,
                       lw=1.0, label=f"N={n_shots}, {tag}ppm 阈值 T")
ax[1].set_xlabel("噪声均值 noise [计数 / 1 ns bin]（步长 0.25）")
ax[1].set_ylabel("临界能量处 peak 均值 / 阈值 [计数]")
ax[1].set_title("② PoD50/90 临界 peak 均值及 threshold")
ax[1].legend(fontsize=6.8, ncol=2); ax[1].grid(alpha=0.3)

# ③ 等效距离
for n_shots in N_SHOTS_LIST:
    for tag in ["100", "10"]:
        for level, alpha, width in [(0.50, 0.65, 1.4), (0.90, 1.0, 2.1)]:
            a = collect_critical(n_shots, tag, level)
            if a.size:
                ax[2].plot(
                    a[:, 0], a[:, 6], styles[tag], color=colors[n_shots],
                    alpha=alpha, lw=width,
                    label=f"N={n_shots}, {tag}ppm, PoD{level:.0%}",
                )
ax[2].set_xlabel("噪声均值 noise [计数 / 1 ns bin]（步长 0.25）")
ax[2].set_ylabel("等效探测距离 [m]（发射能量、ρ=0.10 不变）")
ax[2].set_title("③ 每个 noise 的 PoD50/90 等效距离")
ax[2].legend(fontsize=7.2, ncol=2); ax[2].grid(alpha=0.3)

plt.suptitle(
    "模块 8　完整 0.25-noise 网格：临界能量、peak 均值与距离（滤前）",
    fontsize=12.5,
)
plt.tight_layout(rect=[0, 0, 1, 0.94])
plt.savefig("pod_v04_summary_dense_noise.png", dpi=120, bbox_inches="tight")
plt.show()

# 只打印稀疏检查表；完整 0.25 步长结果保存在 POD_RES / 缓存中。
print("="*116)
print(f"{'N':>3}{'noise':>8}{'FAR':>8}{'PoD目标':>8}{'PoD验证':>9}{'T':>5}"
      f"{'peak均值':>10}{'peak标准差':>11}{'能量[nJ]':>12}{'距离[m]':>10}")
for n_shots in N_SHOTS_LIST:
    stride = max(1, len(NOISE_GRID[n_shots]) // 12)
    for nt in NOISE_GRID[n_shots][::stride]:
        r = POD_RES.get((n_shots, float(nt)))
        if not r or "critical" not in r:
            continue
        for tag in ["100", "10"]:
            for level in POD_LEVELS:
                rec = r["critical"].get(tag, {}).get(f"{level:.2f}")
                if rec:
                    print(f"{n_shots:>3d}{r['noise']:>8.2f}{tag+'ppm':>8}"
                          f"{level:>8.0%}{rec['pod']:>9.3f}{r['T_map'][tag]:>5d}"
                          f"{rec['peak_mean']:>10.2f}{rec['peak_std']:>11.2f}"
                          f"{rec['boost']*E_PULSE_BASE*1e9:>12.3g}"
                          f"{equiv_distance(rec['boost']):>10.1f}")
'''


with NOTEBOOK.open("r", encoding="utf-8") as handle:
    notebook = json.load(handle)

for cell in notebook["cells"]:
    source = "".join(cell.get("source", []))

    if "import json, os, time" in source:
        source = source.replace(
            "import json, os, time",
            "import json, os, time\nfrom concurrent.futures import ThreadPoolExecutor, as_completed",
        )
        source = source.replace(
            "from scipy.stats import binom as _binom",
            "from scipy.stats import binom as _binom, norm as _norm",
        )
        source = source.replace(
            "MC_CHUNK    = 25_000",
            "MC_CHUNK    = 12_500     # 8 线程并行时减小分块，控制峰值内存\n"
            "NOISE_WORKERS = 8       # 纯噪声 MC 分块并行",
        )
        source = source.replace(
            "NOISE_POD  = [1.0, 5.0, 10.0]\n"
            "N_MC_POD   = 3000        # 每个能量档的 MC 条数（PoD≈0.5 时精度 ±0.9%）\n"
            "N_E_GRID   = 24          # 能量档数（对数网格）",
            "# ---- v04：完整 noise 网格上的自适应 PoD 交点 ----\n"
            "NOISE_POD = {n: NOISE_GRID[n].copy() for n in N_SHOTS_LIST}\n"
            "POD_WORKERS = 8          # 信号能量点并行评估\n"
            "N_POD_COARSE = 11        # 全局粗扫描能量点数，只负责包住过渡区\n"
            "N_MC_POD_COARSE = 300    # 每个粗扫描点的 MC 次数\n"
            "N_POD_LOCAL_PER_ROOT = 5 # 每个粗交点附近的局部能量点数\n"
            "N_MC_POD_LOCAL = 800     # 每个局部能量点的 MC 次数\n"
            "N_MC_POD_VERIFY = 5000   # 每个最终临界点的独立验证次数\n"
            "POD_LEVELS = [0.50, 0.90]\n"
            "POD_LOG_BOOST_MIN = -6.0\n"
            "POD_LOG_BOOST_MAX = 2.0\n"
            "POD_LOCAL_HALF_DECADE = 0.22\n"
            "POD_VERIFY_TOL = 0.02",
        )
        cell["source"] = source.splitlines(keepends=True)
    elif source.startswith("def run_noise_scan"):
        source = source.replace(
            "def run_noise_scan(n_shots, noise_grid, n_mc, chunk, seed0=2000, verbose_every=5):",
            "def _noise_chunk_stats(m, n_tr, r_det, inv_tab, seed):\n"
            "    \"\"\"单个纯噪声分块；返回可直接归并的充分统计量。\"\"\"\n"
            "    h = noise_macro_hist_fast(\n"
            "        m, n_tr, r_det, np.random.default_rng(seed), inv_tab=inv_tab,\n"
            "    )\n"
            "    a = h[:, I_STAT0:I_STAT1]\n"
            "    nz = a.mean(axis=1)\n"
            "    return (float(nz.sum()), float((nz*nz).sum()),\n"
            "            np.bincount(a.max(axis=1), minlength=n_tr + 2))\n"
            "\n"
            "\n"
            "def run_noise_scan(n_shots, noise_grid, n_mc, chunk, seed0=2000, verbose_every=5):",
        )
        source = source.replace(
            "        s1 = s2 = 0.0\n"
            "        for s in range(0, n_mc, chunk):\n"
            "            m = min(chunk, n_mc - s)\n"
            "            h = noise_macro_hist_fast(m, n_tr, r_det,\n"
            "                                      np.random.default_rng(seed0 + 1000*k + s),\n"
            "                                      inv_tab=inv_tab)\n"
            "            a = h[:, I_STAT0:I_STAT1]\n"
            "            nz = a.mean(axis=1)\n"
            "            s1 += nz.sum(); s2 += (nz*nz).sum()\n"
            "            res[\"peak_cnt\"][k] += np.bincount(a.max(axis=1), minlength=n_tr + 2)",
            "        s1 = s2 = 0.0\n"
            "        specs = [\n"
            "            (min(chunk, n_mc-s), n_tr, r_det, inv_tab, seed0 + 1000*k + s)\n"
            "            for s in range(0, n_mc, chunk)\n"
            "        ]\n"
            "        if NOISE_WORKERS <= 1:\n"
            "            parts = [_noise_chunk_stats(*spec) for spec in specs]\n"
            "        else:\n"
            "            with ThreadPoolExecutor(max_workers=NOISE_WORKERS) as pool:\n"
            "                parts = list(pool.map(lambda x: _noise_chunk_stats(*x), specs))\n"
            "        for p1, p2, pcnt in parts:\n"
            "            s1 += p1; s2 += p2\n"
            "            res[\"peak_cnt\"][k] += pcnt",
        )
        source = source.replace(
            "    print(f\"  预计耗时约 {_est/60:.0f} 分钟（按实测 ≈ 8 + 4.2×noise 秒/1e6 条估算）\")",
            "    _parallel_eff = 1.0 + 0.65*max(NOISE_WORKERS-1, 0)\n"
            "    print(f\"  单线程基准约 {_est/60:.0f} 分钟；{NOISE_WORKERS} 线程预计 \"\n"
            "          f\"{_est/_parallel_eff/60:.0f} 分钟（受内存带宽限制，不按线程数线性加速）\")",
        )
        cell["source"] = source.splitlines(keepends=True)
    elif source.startswith("## 模块 7 ——"):
        cell["source"] = MODULE7_MD.splitlines(keepends=True)
    elif source.startswith("# ---- PoD 专用的子窗"):
        cell["source"] = MODULE7_CODE.splitlines(keepends=True)
    elif source.startswith("# ---- PoD 曲线图 ----"):
        cell["source"] = MODULE7_PLOT.splitlines(keepends=True)
    elif source.startswith("## 模块 8 —"):
        cell["source"] = MODULE8_MD.splitlines(keepends=True)
    elif source.startswith("_NPH_BASE ="):
        cell["source"] = MODULE8_CODE.splitlines(keepends=True)

    if cell.get("cell_type") == "code":
        cell["execution_count"] = None
        cell["outputs"] = []

with NOTEBOOK.open("w", encoding="utf-8") as handle:
    json.dump(notebook, handle, ensure_ascii=False, indent=1)

print("已增强 PoD_esti_v04.ipynb 的模块 7/8；未执行任何 notebook cell。")
