# -*- coding: utf-8 -*-
"""升级 PoD_esti_v05：扩展 FAR 阈值到 ppm + 百分数，统一 v05 缓存命名。"""
import json
from pathlib import Path

NB = Path(__file__).resolve().parent / "PoD_esti_v05.ipynb"

# ---------- 新 cell 内容 ----------

CELL_PARAM_REPLACEMENTS = [
    (
        "TARGET_FARS = [100e-6, 10e-6]   # ★ v02：同时给 100 ppm 与 10 ppm 两条阈值\n",
        # 百分数 FAR 与 ppm 并存；tag 用纯 ASCII，避免字典键里出现 %
        """# ---- ★ v05：噪点率 / FAR 目标（ppm + 百分数）----
# tag 用于字典键与文件字段；label 用于图例显示
FAR_SPECS = [
    (10e-6,  \"10ppm\",  \"10 ppm\"),
    (100e-6, \"100ppm\", \"100 ppm\"),
    (0.001,  \"0p1pct\", \"0.1%\"),
    (0.005,  \"0p5pct\", \"0.5%\"),
    (0.01,   \"1pct\",   \"1%\"),
    (0.05,   \"5pct\",   \"5%\"),
]
TARGET_FARS = [v for v, _, _ in FAR_SPECS]
FAR_TAG     = {v: t for v, t, _ in FAR_SPECS}
FAR_LABEL   = {v: lab for v, _, lab in FAR_SPECS}
FAR_TAGS    = [t for _, t, _ in FAR_SPECS]
FAR_TAG_TO_LABEL = {t: lab for _, t, lab in FAR_SPECS}
"""
    ),
    (
        'CACHE_NOISE_FALLBACK = ["pod_esti_v04_cache_noise.npz"]\n'
        'CACHE_POD_FALLBACK   = ["pod_esti_v04_cache_pod.npz"]\n',
        # 噪声 MC 与 FAR 无关，仍允许一次性从 v04 迁入 v05；
        # PoD 依赖阈值集合，禁止再读 v04（FAR 集合已变）。
        'CACHE_NOISE_FALLBACK = ["pod_esti_v04_cache_noise.npz"]  # 仅噪声可迁入\n'
        'CACHE_POD_FALLBACK   = []  # ★ FAR 已扩展，禁止误用旧 PoD 缓存\n',
    ),
    (
        'print(f"  ★ 噪点率目标：{[f\'{f*1e6:.0f} ppm\' for f in TARGET_FARS]}，每档 {N_MC_NOISE:,} 次 MC")\n'
        'print(f"  ★ 并行：N_WORKERS={N_WORKERS}，MC_CHUNK={MC_CHUNK:,}，每 {CHECKPOINT_EVERY} 档增量落盘")\n'
        'print(f"  ★ 缓存主文件：{CACHE_NOISE} / {CACHE_POD}（可读 fallback：v04）")\n',
        'print(f"  ★ FAR 目标：{[FAR_LABEL[f] for f in TARGET_FARS]}，每档 {N_MC_NOISE:,} 次 MC")\n'
        'print(f"  ★ 并行：N_WORKERS={N_WORKERS}，MC_CHUNK={MC_CHUNK:,}，每 {CHECKPOINT_EVERY} 档增量落盘")\n'
        'print(f"  ★ 缓存主文件：{CACHE_NOISE} / {CACHE_POD}")\n'
        'print(f"  ★ 噪声可从 v04 fallback 迁入；PoD 因 FAR 变更必须用 v05 新缓存重跑")\n',
    ),
]

CELL_THRESH = r'''def far_threshold_from_cnt(cnt, target_far):
    """由 peak 的 bincount 求满足 P(peak ≥ T) < target_far 的最小整数 T。

    全程用【整数计数】比较（n_ge < target_far·n），避免浮点边界误判。
    返回 (T, 该 T 处实测 FAR, 该 T 处越阈次数, 生存函数数组)。
    """
    n = int(cnt.sum())
    n_ge = np.concatenate([[n], n - np.cumsum(cnt)[:-1]])
    lim = target_far * n
    ok = np.where(n_ge < lim)[0]
    sf = n_ge / n
    if ok.size == 0:
        return int(cnt.size), 0.0, 0, sf
    T = int(ok[0])
    return T, float(sf[T]), int(n_ge[T]), sf


def far_threshold_binom_indep(n_tr, p_bin, n_bins, target_far):
    """独立 Binomial 近似阈值（保守对照）。"""
    a_bin = 1.0 - (1.0 - target_far) ** (1.0 / n_bins)
    T = 0
    while T <= n_tr and _binom.sf(T - 1, n_tr, p_bin) > a_bin:
        T += 1
    return T


THRESH = {}
for n_shots in N_SHOTS_LIST:
    R = NOISE_RES[n_shots]
    ng = len(R["noise_target"])
    rec = {"noise": R["noise_mc"], "sigma_bin": np.zeros(ng)}
    for far in TARGET_FARS:
        tag = FAR_TAG[far]
        rec["T"+tag] = np.zeros(ng, dtype=int)
        rec["far"+tag] = np.zeros(ng)
        rec["nev"+tag] = np.zeros(ng, dtype=int)
        rec["Ti"+tag] = np.zeros(ng, dtype=int)
    for k in range(ng):
        rec["sigma_bin"][k] = np.sqrt(R["n_tr"] * R["p_eq"][k] * (1 - R["p_eq"][k]))
        for far in TARGET_FARS:
            tag = FAR_TAG[far]
            T, f_, nev, _ = far_threshold_from_cnt(R["peak_cnt"][k], far)
            rec["T"+tag][k] = T
            rec["far"+tag][k] = f_
            rec["nev"+tag][k] = nev
            rec["Ti"+tag][k] = far_threshold_binom_indep(
                R["n_tr"], R["p_eq"][k], N_STAT, far)
    THRESH[n_shots] = rec

print("="*140)
print(f"检测阈值汇总（每档 {N_MC_NOISE:,} 条 MC；FAR = { [FAR_LABEL[f] for f in TARGET_FARS] }）")
_hdr = f"{'N_shots':>8}{'noise':>8}"
for far in TARGET_FARS:
    _hdr += f"{'T@'+FAR_LABEL[far]:>10}"
_hdr += f"{'上限':>6}"
print(_hdr)
for n_shots in N_SHOTS_LIST:
    R, Tr = NOISE_RES[n_shots], THRESH[n_shots]
    step = max(1, len(Tr["noise"]) // 12)
    for k in range(0, len(Tr["noise"]), step):
        row = f"{n_shots:>8d}{Tr['noise'][k]:>8.3f}"
        for far in TARGET_FARS:
            row += f"{Tr['T'+FAR_TAG[far]][k]:>10d}"
        row += f"{R['n_tr']:>6d}"
        print(row)
'''

CELL_THRESH_PLOT = r'''# ---- 阈值图：6 条 FAR（2 条 ppm + 4 条百分数）----
fig, ax = plt.subplots(1, 3, figsize=(18.5, 5.4))
_cns = {1: "tab:blue", 4: "tab:red"}
# 线型按 FAR 从严到松区分
_ls_by_tag = {
    "100ppm": "-", "10ppm": "--",
    "5pct": "-.", "1pct": ":",
    "0p5pct": (0, (3, 1, 1, 1)), "0p1pct": (0, (1, 1)),
}

# ① noise–threshold
for n_shots in N_SHOTS_LIST:
    R, Tr = NOISE_RES[n_shots], THRESH[n_shots]
    c = _cns.get(n_shots, "k")
    for far in TARGET_FARS:
        tag = FAR_TAG[far]
        ax[0].plot(
            Tr["noise"], Tr["T"+tag], ls=_ls_by_tag.get(tag, "-"),
            color=c, lw=1.8,
            label=f"N={n_shots}, {FAR_LABEL[far]}",
        )
    ax[0].axhline(R["n_tr"], color=c, ls="-.", lw=1.0, alpha=0.4,
                  label=f"N={n_shots} 硬上限 {R['n_tr']}")
ax[0].set_xlabel("噪声均值 noise [计数 / 1 ns bin]")
ax[0].set_ylabel("检测阈值 T [计数 / 1 ns bin]")
ax[0].set_title("① noise–threshold（色=N_shots，线型=FAR）", fontsize=11)
ax[0].legend(fontsize=6.8, ncol=2, loc="upper left")
ax[0].grid(alpha=0.3)

# ② peak 生存函数 + 全部 FAR 水平线
_ns_sf = N_SHOTS_LIST[-1]
R = NOISE_RES[_ns_sf]
_sel = np.linspace(0, len(R["noise_target"])-1, 7).astype(int)
_cols = plt.cm.viridis(np.linspace(0.08, 0.92, len(_sel)))
for c, k in zip(_cols, _sel):
    cnt = R["peak_cnt"][k]; n = cnt.sum()
    n_ge = np.concatenate([[n], n - np.cumsum(cnt)[:-1]])
    Ts = np.arange(n_ge.size); m = n_ge > 0
    ax[1].semilogy(Ts[m], n_ge[m]/n, "-", color=c, lw=1.5,
                   label=f"noise={R['noise_mc'][k]:.2f}")
_far_colors = plt.cm.Reds(np.linspace(0.35, 0.95, len(TARGET_FARS)))
for far, fc in zip(TARGET_FARS, _far_colors):
    ax[1].axhline(far, color=fc, ls="--", lw=1.4,
                  label=f"目标 {FAR_LABEL[far]}")
ax[1].axhline(1.0/N_MC_NOISE, color="0.5", ls="-.", lw=1.1,
              label=f"MC 分辨极限 1/{N_MC_NOISE:,}")
ax[1].set_ylim(0.5/N_MC_NOISE, 1.5)
ax[1].set_xlabel("检测阈值 T [计数 / 1 ns bin]")
ax[1].set_ylabel("窗口级噪点率 P(peak ≥ T)")
ax[1].set_title(f"② peak 生存函数（N_shots={_ns_sf}）", fontsize=11)
ax[1].legend(fontsize=6.8, ncol=2)
ax[1].grid(alpha=0.3, which="both")

# ③ k_th = T/noise
for n_shots in N_SHOTS_LIST:
    Tr = THRESH[n_shots]
    for far in TARGET_FARS:
        tag = FAR_TAG[far]
        ax[2].plot(
            Tr["noise"], Tr["T"+tag]/np.maximum(Tr["noise"], 1e-9),
            ls=_ls_by_tag.get(tag, "-"), color=_cns.get(n_shots, "k"), lw=1.7,
            label=f"N={n_shots}, {FAR_LABEL[far]}",
        )
ax[2].axhline(5.0, color="0.4", ls=":", lw=1.4, label="v45 固定 k_th=5")
ax[2].set_xlabel("噪声均值 noise [计数 / 1 ns bin]")
ax[2].set_ylabel("阈值倍数 k_th = T / noise")
ax[2].set_ylim(0, 30)
ax[2].set_title("③ 阈值相对底噪的倍数", fontsize=11)
ax[2].legend(fontsize=6.8, ncol=2)
ax[2].grid(alpha=0.3)

plt.suptitle(
    f"模块 6　FAR 阈值曲线：{[FAR_LABEL[f] for f in TARGET_FARS]}（统计窗 {N_STAT} bins）",
    fontsize=12,
)
plt.tight_layout(rect=[0, 0, 1, 0.93])
plt.savefig("pod_v05_threshold.png", dpi=110, bbox_inches="tight")
plt.show()

# 相邻 FAR 之间的 ΔT
print("="*90)
print("相邻 FAR（从严到松排列的反序：从松到严）阈值抬升量 ΔT：")
_order = list(TARGET_FARS)  # FAR_SPECS 顺序：10ppm → 100ppm → 0.1% → 0.5% → 1% → 5%
# 按 FAR 数值从松到严排序打印差分
_sorted = sorted(TARGET_FARS, reverse=True)  # 大 FAR → 小 FAR
for n_shots in N_SHOTS_LIST:
    Tr = THRESH[n_shots]
    print(f"  N_shots={n_shots}:")
    for a, b in zip(_sorted[:-1], _sorted[1:]):
        d = Tr["T"+FAR_TAG[b]] - Tr["T"+FAR_TAG[a]]
        print(f"    {FAR_LABEL[a]} → {FAR_LABEL[b]}: 平均 ΔT=+{d.mean():.2f} "
              f"（范围 {d.min()}~{d.max()}）")
'''


def patch_5b(src: str) -> str:
    old = '''    threshold_100 = np.array([
        _threshold_from_cnt_5b(cnt, 100e-6) for cnt in R["peak_cnt"]
    ])
    threshold_10 = np.array([
        _threshold_from_cnt_5b(cnt, 10e-6) for cnt in R["peak_cnt"]
    ])
    ax.plot(
        R["noise_mc"], threshold_100, color="#e63946", lw=1.9, zorder=5,
        label="100 ppm 阈值",
    )
    ax.plot(
        R["noise_mc"], threshold_10, color="#9b2226", lw=1.9, ls="--",
        zorder=5, label="10 ppm 阈值",
    )'''
    new = '''    _far_line = {
        "100ppm": ("-",  "#e63946"),
        "10ppm":  ("--", "#9b2226"),
        "5pct":   ("-.", "#f4a261"),
        "1pct":   (":",  "#e9c46a"),
        "0p5pct": ((0, (3, 1, 1, 1)), "#2a9d8f"),
        "0p1pct": ((0, (1, 1)), "#264653"),
    }
    for far in TARGET_FARS:
        tag = FAR_TAG[far]
        ls, col = _far_line.get(tag, ("-", "k"))
        thr = np.array([_threshold_from_cnt_5b(cnt, far) for cnt in R["peak_cnt"]])
        ax.plot(R["noise_mc"], thr, color=col, ls=ls, lw=1.6, zorder=5,
                label=f"{FAR_LABEL[far]} 阈值")'''
    if old not in src:
        raise RuntimeError("5b threshold block not found")
    return src.replace(old, new)


def patch_pod(src: str) -> str:
    src = src.replace(
        '    """求一个 noise 档、两条 FAR 阈值下的 PoD50/90 临界点。"""',
        '    """求一个 noise 档、全部 FAR 阈值下的 PoD50/90 临界点。"""',
    )
    src = src.replace(
        '    T_map = {f"{far*1e6:.0f}": int(Tr[f"T{far*1e6:.0f}"][k]) for far in TARGET_FARS}',
        '    T_map = {FAR_TAG[far]: int(Tr["T"+FAR_TAG[far]][k]) for far in TARGET_FARS}',
    )
    src = src.replace(
        "    # 围绕四个粗交点补局部点；所有阈值共享同一批 peak 样本。",
        "    # 围绕各 FAR×PoD 粗交点补局部点；所有阈值共享同一批 peak 样本。",
    )
    # cache load/save with far_tags
    old_try = '''def _try_load_pod_cache(path, grid_key):
    if not os.path.exists(path):
        return None
    z = np.load(path, allow_pickle=True)
    if (
        np.array_equal(z["grid_key"], grid_key)
        and int(z["n_coarse"]) == N_MC_POD_COARSE
        and int(z["n_local"]) == N_MC_POD_LOCAL
        and int(z["n_verify"]) == N_MC_POD_VERIFY
    ):
        return z["res"].item()
    return None


def _save_pod_cache(path, res, grid_key):
    # 元组键无法直接进 npz；统一转成可 pickle 的 dict
    _atomic_savez(
        path, res=np.array(res, dtype=object),
        grid_key=grid_key,
        n_coarse=N_MC_POD_COARSE, n_local=N_MC_POD_LOCAL,
        n_verify=N_MC_POD_VERIFY,
    )'''
    new_try = '''def _try_load_pod_cache(path, grid_key):
    """仅当噪声网格、MC 精度与 FAR 标签集合全部一致时才接受缓存。"""
    if not os.path.exists(path):
        return None
    z = np.load(path, allow_pickle=True)
    if "far_tags" not in z.files:
        return None  # 旧缓存无 FAR 标签 → 一律作废
    if (
        np.array_equal(z["grid_key"], grid_key)
        and list(z["far_tags"]) == list(FAR_TAGS)
        and int(z["n_coarse"]) == N_MC_POD_COARSE
        and int(z["n_local"]) == N_MC_POD_LOCAL
        and int(z["n_verify"]) == N_MC_POD_VERIFY
    ):
        return z["res"].item()
    return None


def _save_pod_cache(path, res, grid_key):
    _atomic_savez(
        path, res=np.array(res, dtype=object),
        grid_key=grid_key,
        far_tags=np.array(FAR_TAGS),
        n_coarse=N_MC_POD_COARSE, n_local=N_MC_POD_LOCAL,
        n_verify=N_MC_POD_VERIFY,
    )'''
    if old_try not in src:
        raise RuntimeError("pod cache helpers not found")
    src = src.replace(old_try, new_try)

    src = src.replace(
        '                c100 = value.get("critical", {}).get("100", {})\n'
        '                p90 = c100.get("0.90")\n',
        '                c100 = value.get("critical", {}).get("100ppm", {})\n'
        '                p90 = c100.get("0.90")\n',
    )

    # completeness: also require each entry has all FAR tags
    old_complete = '''_expected_keys = {(ns, float(nt)) for ns in N_SHOTS_LIST for nt in NOISE_GRID[ns]}
_have_keys = set(POD_RES.keys()) if POD_RES else set()
_complete = POD_RES is not None and _expected_keys.issubset(_have_keys)'''
    new_complete = '''_expected_keys = {(ns, float(nt)) for ns in N_SHOTS_LIST for nt in NOISE_GRID[ns]}
_have_keys = set(POD_RES.keys()) if POD_RES else set()

def _pod_entry_has_all_fars(ent):
    crit = ent.get("critical", {})
    return all(tag in crit for tag in FAR_TAGS)

_complete = (
    POD_RES is not None
    and _expected_keys.issubset(_have_keys)
    and all(_pod_entry_has_all_fars(POD_RES[k]) for k in _expected_keys)
)'''
    if old_complete not in src:
        raise RuntimeError("pod complete check not found")
    src = src.replace(old_complete, new_complete)
    return src


CELL_MOD8 = r'''def equiv_distance(boost, D_ref=D_TARGET, p=PARAMS):
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


_ls_by_tag = {
    "100ppm": "-", "10ppm": "--",
    "5pct": "-.", "1pct": ":",
    "0p5pct": (0, (3, 1, 1, 1)), "0p1pct": (0, (1, 1)),
}
_cns = {1: "tab:blue", 4: "tab:red"}

# 主交付：只画 PoD90，避免 6×2×2 条线过载；PoD50 另存于表
fig, ax = plt.subplots(1, 3, figsize=(19, 5.8))
for n_shots in N_SHOTS_LIST:
    for far in TARGET_FARS:
        tag = FAR_TAG[far]
        a = collect_critical(n_shots, tag, 0.90)
        if not a.size:
            continue
        ls = _ls_by_tag.get(tag, "-")
        c = _cns[n_shots]
        ax[0].semilogy(a[:, 0], a[:, 1]*E_PULSE_BASE*1e9, ls=ls, color=c, lw=1.7,
                       label=f"N={n_shots}, {FAR_LABEL[far]}")
        ax[1].plot(a[:, 0], a[:, 3], ls=ls, color=c, lw=1.7,
                   label=f"N={n_shots}, {FAR_LABEL[far]} peak均值")
        ax[1].plot(a[:, 0], a[:, 5], ls=":", color=c, alpha=0.35, lw=1.0)
        ax[2].plot(a[:, 0], a[:, 6], ls=ls, color=c, lw=1.7,
                   label=f"N={n_shots}, {FAR_LABEL[far]}")

ax[0].set_xlabel("noise"); ax[0].set_ylabel("PoD90 临界能量 [nJ]")
ax[0].set_title("① PoD90 临界发射能量 vs noise")
ax[0].legend(fontsize=6.5, ncol=2); ax[0].grid(alpha=0.3, which="both")

ax[1].set_xlabel("noise"); ax[1].set_ylabel("peak 均值 / T [计数]")
ax[1].set_title("② PoD90 临界 peak 均值（点线≈对应 T）")
ax[1].legend(fontsize=6.2, ncol=2); ax[1].grid(alpha=0.3)

ax[2].set_xlabel("noise"); ax[2].set_ylabel("等效距离 [m]")
ax[2].set_title("③ PoD90 等效探测距离 vs noise")
ax[2].legend(fontsize=6.5, ncol=2); ax[2].grid(alpha=0.3)

plt.suptitle(
    f"模块 8　完整 0.25-noise × FAR{[FAR_LABEL[f] for f in TARGET_FARS]}（PoD90，滤前）",
    fontsize=12,
)
plt.tight_layout(rect=[0, 0, 1, 0.93])
plt.savefig("pod_v05_summary_dense_noise.png", dpi=120, bbox_inches="tight")
plt.show()

print("="*120)
print(f"{'N':>3}{'noise':>8}{'FAR':>10}{'PoD目标':>8}{'PoD验证':>9}{'T':>5}"
      f"{'peak均值':>10}{'能量[nJ]':>12}{'距离[m]':>10}")
for n_shots in N_SHOTS_LIST:
    stride = max(1, len(NOISE_GRID[n_shots]) // 10)
    for nt in NOISE_GRID[n_shots][::stride]:
        r = POD_RES.get((n_shots, float(nt)))
        if not r or "critical" not in r:
            continue
        for far in TARGET_FARS:
            tag = FAR_TAG[far]
            for level in POD_LEVELS:
                rec = r["critical"].get(tag, {}).get(f"{level:.2f}")
                if not rec:
                    continue
                print(f"{n_shots:>3d}{r['noise']:>8.2f}{FAR_LABEL[far]:>10}"
                      f"{level:>8.0%}{rec['pod']:>9.3f}{r['T_map'][tag]:>5d}"
                      f"{rec['peak_mean']:>10.2f}"
                      f"{rec['boost']*E_PULSE_BASE*1e9:>12.3g}"
                      f"{equiv_distance(rec['boost']):>10.1f}")
'''


def main():
    nb = json.loads(NB.read_text(encoding="utf-8"))
    changed = []
    for i, cell in enumerate(nb["cells"]):
        if cell.get("cell_type") != "code":
            continue
        src = "".join(cell.get("source", []))
        old = src

        if "TARGET_FARS = [100e-6, 10e-6]" in src and "CACHE_NOISE" in src:
            for a, b in CELL_PARAM_REPLACEMENTS:
                if a not in src:
                    raise RuntimeError(f"param replace miss: {a[:40]!r}")
                src = src.replace(a, b)
            changed.append(("PARAM", i))

        if src.startswith("def far_threshold_from_cnt"):
            src = CELL_THRESH
            changed.append(("THRESH", i))

        if src.startswith("# ---- 阈值图"):
            src = CELL_THRESH_PLOT
            changed.append(("THRESH_PLOT", i))

        if "模块 5b：noise–peak 密度条带" in src or src.startswith("# ---- 模块 5b"):
            src = patch_5b(src)
            changed.append(("5B", i))

        if "def solve_pod_noise" in src:
            src = patch_pod(src)
            changed.append(("POD", i))

        if src.startswith("def equiv_distance") and "collect_critical" in src:
            src = CELL_MOD8
            changed.append(("MOD8", i))

        if src != old:
            cell["source"] = src.splitlines(keepends=True)
        cell["execution_count"] = None
        cell["outputs"] = []

    # clear all outputs
    for cell in nb["cells"]:
        if cell.get("cell_type") == "code":
            cell["execution_count"] = None
            cell["outputs"] = []

    NB.write_text(json.dumps(nb, ensure_ascii=False, indent=1), encoding="utf-8")
    print("patched:", changed)


if __name__ == "__main__":
    main()
