# -*- coding: utf-8 -*-
"""
build_crosstalk_v24.py —— 生成 crosstalk_sim_v24.ipynb
========================================================
继承：v23 全部 cell（含 v13/v20/v21 + 对射/tcode/FPGA 波形）。
新增：一字滤波（连续三角度空间一致性）。

规则（用户定义）：
  连续三次仿真 = 连续三个角度（左 / 中 / 右）。
  对中间角度的某个峰，若它相对左右两侧最近峰的距离差
  均大于阈值 thr，则判为串扰/鬼影并丢掉。

缩写：
  XM（XtalkMark，串扰标记）
  FPGA（Field-Programmable Gate Array，现场可编程门阵列）
  TOF（Time of Flight，飞行时间）
"""
import json
import os

SRC_NB = "crosstalk_sim_v23.ipynb"
OUT_NB = "crosstalk_sim_v24.ipynb"


def code_cell(cid, source):
    return {"cell_type": "code", "id": cid, "metadata": {},
            "execution_count": None, "outputs": [],
            "source": source.splitlines(keepends=True)}


def md_cell(cid, source):
    return {"cell_type": "markdown", "id": cid, "metadata": {},
            "source": source.splitlines(keepends=True)}


OVERVIEW = """# 串扰仿真 v24 —— 一字滤波（连续三角度）

> **本 notebook = v23 原封不动 + 末尾追加 v24。**
> v23 证明了：tcode 管模组内鬼影，FPGA 管对射。
> v24 补上第三层：**角度上孤立的假点**，用一字滤波清掉。

## 一字滤波规则

连续三次仿真 = 左 / 中 / 右三个角度。对**中间**角度每个（XM 后仍保留的）峰：

```
dL = min |d_mid − d_left峰|
dR = min |d_mid − d_right峰|
若 dL > thr 且 dR > thr  →  判孤立假点，丢掉
否则                       →  保留
```

按**同一激光器通道**跨角度比较（三次仿真里的同一 `lid`）。

## 本版要证明的分工

| 手段 | 管什么 | 管不了什么 |
|---|---|---|
| tcode | 模组内稳定串扰 | 对射、角度孤立噪点 |
| FPGA 累计抖动 | 持续同频对射 | 偶然单角度闪烁 |
| **一字滤波** | 角度上孤立的假点 | 三角度都稳定出现的假峰（需靠前两层） |
"""


SECTION = """---
# v24 新增部分 —— 一字滤波

先看规则与函数，再跑三个对照场景。
"""


CELL_PARAMS = r"""# ============================================================================
# V24-1  一字滤波参数
# ============================================================================
LINE_THR_M      = 3.0     # 左右距离差阈值 [m]；两差都 > 此值 → 丢掉
LINE_DEMO_LASER = 5       # 解剖用激光器
LINE_FAKE_DIST  = 80.0    # 人为塞进「仅中间角度」的稳定假峰 [m]

# 三角度真目标距离（可轻微起伏；默认贴近平面墙）
LINE_D_LEFT   = 150.0
LINE_D_MID    = 150.0
LINE_D_RIGHT  = 150.0

print("V24-1 一字滤波参数")
print(f"  thr = {LINE_THR_M:.1f} m")
print(f"  三角度真距离 = [{LINE_D_LEFT:.0f}, {LINE_D_MID:.0f}, {LINE_D_RIGHT:.0f}] m")
print(f"  孤立假峰演示距离 = {LINE_FAKE_DIST:.0f} m（仅中间角度注入）")
"""


CELL_CORE = r"""# ============================================================================
# V24-2  一字滤波核心 + 三角度仿真辅助
# ============================================================================
def xm_kept_peaks(rs, lid):
    # XM 后仍保留的峰
    return [dict(q) for q in rs[lid]["peaks"] if not q["is_xtalk"]]


def nearest_abs_dist(d, peaks):
    if not peaks:
        return float("inf")
    return float(min(abs(d - p["dist"]) for p in peaks))


def line_filter_one(peaks_L, peaks_M, peaks_R, thr_m=None):
    # 对中间角度峰做一字判决；返回带 is_line_ghost / dL / dR 的新列表
    if thr_m is None:
        thr_m = LINE_THR_M
    out = []
    for q in peaks_M:
        dL = nearest_abs_dist(q["dist"], peaks_L)
        dR = nearest_abs_dist(q["dist"], peaks_R)
        iso = (dL > thr_m) and (dR > thr_m)
        qq = dict(q)
        qq["dL"] = dL
        qq["dR"] = dR
        qq["is_line_ghost"] = bool(iso)
        out.append(qq)
    return out


def apply_line_filter_wave(rr, line_peaks):
    # 在 XM after 波形上再清掉一字滤波判掉的峰
    after2 = rr["after"].copy()
    for q in line_peaks:
        if q.get("is_line_ghost"):
            after2[q["s"]:q["e"] + 1] = 0.0
    return after2


def line_filter_all(rs_L, rs_M, rs_R, thr_m=None):
    # 16 激光器全部跑一字滤波
    out = {}
    for lid in laser_ids:
        lp = line_filter_one(xm_kept_peaks(rs_L, lid),
                             xm_kept_peaks(rs_M, lid),
                             xm_kept_peaks(rs_R, lid),
                             thr_m=thr_m)
        out[lid] = {
            "peaks": lp,
            "after_xm": rs_M[lid]["after"],
            "after_line": apply_line_filter_wave(rs_M[lid], lp),
        }
    return out


def run_angle(fr, ratio, D, radar=False, use_fpga=False,
              n_levels=None, seed=None, phase_ns=None):
    # 单角度仿真；radar/use_fpga 可按角度开关
    if n_levels is None:
        n_levels = FPGA_N_LEVELS
    if seed is None:
        seed = FPGA_SEED
    cum = None
    inc = None
    if use_fpga:
        inc, cum = make_fpga_tables(n_levels=n_levels, seed=seed)
    with use_firings(fr):
        recs = detect_echoes(D)
    if radar:
        recs = inject_radar(recs, D, cum_table=cum, phase_ns=phase_ns)
    hs = build_hists(recs)
    rs = crosstalk_mark_all(hs, ratio)
    return {"recs": recs, "hs": hs, "rs": rs, "inc": inc, "cum": cum,
            "ratio": ratio, "D": D, "radar": radar}


def inject_stable_fake_peak(hs, lid, dist_m, amp_per_shot=1.0):
    # 往某激光器 4 次 shot 同一 bin 塞稳定假峰 → XM 会当成「真」留下
    b = int(round(dist_m / HIST_BIN_M))
    b = max(0, min(N_BINS - 1, b))
    for i in range(hs[lid]["shots"].shape[0]):
        hs[lid]["shots"][i, b] += amp_per_shot
    hs[lid]["add"] = hs[lid]["shots"].sum(axis=0)
    hs[lid]["max"] = hs[lid]["shots"].max(axis=0)
    # 源记录：标记为假峰，便于统计
    fake_rec = {
        "emit_laser": -2, "emit_kick": -1,
        "recv_laser": lid, "recv_kick": -1,
        "target_D": dist_m, "true_tof": 0.0, "t_echo": 0.0,
        "rec_tof": (dist_m * 2.0 / C_LIGHT),
        "rec_dist": dist_m,
        "is_true": False, "negligible": False,
        "is_radar": False, "is_fake_iso": True,
    }
    hs[lid]["src"].setdefault(b, []).append(fake_rec)
    return b


def recount_with_fake(hs, ratio):
    return crosstalk_mark_all(hs, ratio)


def summarize_line(tag, line_res, hs_M=None):
    n_keep = n_drop = 0
    n_true_keep = n_true_drop = 0
    n_fake_keep = n_fake_drop = 0
    n_ghost_keep = n_ghost_drop = 0
    for lid, rr in line_res.items():
        for q in rr["peaks"]:
            kind_fake = kind_true = kind_ghost = False
            if hs_M is not None:
                recs = [r for b in range(q["s"], q["e"] + 1)
                        for r in hs_M[lid]["src"].get(b, [])]
                kind_fake = any(r.get("is_fake_iso") or r.get("is_radar") for r in recs)
                kind_true = any(r.get("is_true") for r in recs)
                kind_ghost = (not kind_true) and (not kind_fake) and len(recs) > 0
            if q["is_line_ghost"]:
                n_drop += 1
                if kind_true and not kind_fake:
                    n_true_drop += 1
                if kind_fake:
                    n_fake_drop += 1
                if kind_ghost:
                    n_ghost_drop += 1
            else:
                n_keep += 1
                if kind_true and not kind_fake:
                    n_true_keep += 1
                if kind_fake:
                    n_fake_keep += 1
                if kind_ghost:
                    n_ghost_keep += 1
    print(f"\n【{tag}】一字滤波（thr={LINE_THR_M:.1f}m）")
    print(f"  中间角度 XM 后保留峰：丢掉 {n_drop}，留下 {n_keep}")
    print(f"  真峰     留下 {n_true_keep} / 误杀 {n_true_drop}")
    print(f"  假/对射  留下 {n_fake_keep} / 丢掉 {n_fake_drop}")
    print(f"  模组鬼影 留下 {n_ghost_keep} / 丢掉 {n_ghost_drop}")
    return {"keep": n_keep, "drop": n_drop,
            "true_keep": n_true_keep, "true_drop": n_true_drop,
            "fake_keep": n_fake_keep, "fake_drop": n_fake_drop,
            "ghost_keep": n_ghost_keep, "ghost_drop": n_ghost_drop}


print("V24-2 一字滤波函数已加载")
"""


CELL_PLOT = r"""# ============================================================================
# V24-3  一字滤波绘图
# ============================================================================
def plot_three_angle_strip(ang_L, ang_M, ang_R, line_res, lid, title,
                           thr_m=None, highlight_fake=None):
    # 三角度 × 距离散点：XM 后峰；中间角度标出一字判决
    if thr_m is None:
        thr_m = LINE_THR_M
    fig, ax = plt.subplots(figsize=(13, 4.2))
    rows = [("左", ang_L, 2), ("中", ang_M, 1), ("右", ang_R, 0)]
    for name, ang, y in rows:
        for q in xm_kept_peaks(ang["rs"], lid):
            ax.scatter([q["dist"]], [y], s=70, color=laser_color(lid),
                       zorder=3, edgecolors="k", lw=0.6)
        ax.axhline(y, color="0.85", lw=0.8)
        ax.text(-0.02, y, name, transform=ax.get_yaxis_transform(),
                ha="right", va="center", fontsize=11)

    # 中间角度：一字丢掉的画叉
    for q in line_res[lid]["peaks"]:
        if q["is_line_ghost"]:
            ax.plot(q["dist"], 1, "kx", ms=12, mew=2.2, zorder=5)
        else:
            ax.scatter([q["dist"]], [1], s=120, facecolors="none",
                       edgecolors="green", lw=2.0, zorder=4)

    ax.axvline(ang_M["D"], color="k", ls=":", lw=1.2, label=f"真目标≈{ang_M['D']:.0f}m")
    if highlight_fake is not None:
        ax.axvline(highlight_fake, color="#c0392b", ls="--", lw=1.2,
                   label=f"孤立假峰 {highlight_fake:.0f}m")
    ax.set_yticks([0, 1, 2])
    ax.set_yticklabels(["右", "中", "左"])
    ax.set_xlabel("记录距离 [m]")
    ax.set_title(title + f"  （绿圈=一字保留，黑叉=一字丢掉，thr={thr_m:.1f}m）")
    ax.legend(fontsize=8); ax.grid(alpha=0.25, axis="x")
    plt.tight_layout(); plt.show()


def plot_mid_before_after_line(ang_M, line_res, lid, title, D):
    # 中间角度：XM 后 vs 一字滤波后
    h = ang_M["hs"][lid]
    rr = line_res[lid]
    nz = np.flatnonzero((h["add"] > 0) | (rr["after_xm"] > 0) | (rr["after_line"] > 0))
    if len(nz) == 0:
        b_lo, b_hi = 0, N_BINS - 1
    else:
        pad = 40
        b_lo = max(0, int(nz.min()) - pad)
        b_hi = min(N_BINS - 1, int(nz.max()) + pad)
    x = np.arange(b_lo, b_hi + 1) * HIST_BIN_M

    fig, axes = plt.subplots(2, 1, figsize=(13, 6.5), sharex=True)
    ax = axes[0]
    ax.fill_between(x, 0, h["add"][b_lo:b_hi + 1], step="mid",
                    color="0.78", label="XM 前 hist_add")
    ax.plot(x, rr["after_xm"][b_lo:b_hi + 1], drawstyle="steps-mid",
            color="#2874a6", lw=2.0, label="XM 后")
    ax.axvline(D, color="k", ls=":", lw=1.2)
    ax.set_ylabel("计数"); ax.legend(fontsize=8); ax.grid(alpha=0.2)
    ax.set_title("(a) XM 之后（一字滤波之前）")

    ax = axes[1]
    ax.plot(x, rr["after_xm"][b_lo:b_hi + 1], drawstyle="steps-mid",
            color="0.7", lw=1.4, label="XM 后（参考）")
    ax.plot(x, rr["after_line"][b_lo:b_hi + 1], drawstyle="steps-mid",
            color=laser_color(lid), lw=2.2, label="一字滤波后")
    for q in rr["peaks"]:
        if q["is_line_ghost"]:
            ax.plot(q["dist"], max(q["add"], 0.3), "kx", ms=10, mew=2)
    ax.axvline(D, color="k", ls=":", lw=1.2, label=f"真目标 {D:.0f}m")
    ax.set_xlabel("记录距离 [m]"); ax.set_ylabel("计数")
    ax.set_title("(b) 一字滤波后（黑叉 = 角度孤立被丢掉）")
    ax.legend(fontsize=8); ax.grid(alpha=0.2)
    fig.suptitle(title, fontsize=13)
    plt.tight_layout(); plt.show()


print("V24-3 绘图函数已加载")
"""


# ---- 场景 A：孤立假峰 ----
CELL_A_DOC = """# V24-A 【证明一字滤波有用】中间角度塞孤立假峰

配置：

- 三角度都用 **tcode**，真目标都在 ~150 m（平面墙）→ 真峰跨角度对齐
- **不**开对射 / FPGA
- 只在**中间**角度、L5 上人为塞一个稳定假峰（4 次 shot 同位 → XM 滤不掉）

期望：一字滤波丢掉假峰，真峰留下。
"""

CELL_A = r"""# ============================================================================
# V24-A  孤立假峰演示
# ============================================================================
fr = ACTIVE_FR
ratio = ACTIVE_RATIO

A_L = run_angle(fr, ratio, LINE_D_LEFT,  radar=False)
A_M = run_angle(fr, ratio, LINE_D_MID,   radar=False)
A_R = run_angle(fr, ratio, LINE_D_RIGHT, radar=False)

# 只在中间角度 L5 塞假峰，并重跑 XM
inject_stable_fake_peak(A_M["hs"], LINE_DEMO_LASER, LINE_FAKE_DIST)
A_M["rs"] = recount_with_fake(A_M["hs"], ratio)

line_A = line_filter_all(A_L["rs"], A_M["rs"], A_R["rs"])
stat_A = summarize_line("A 孤立假峰", line_A, hs_M=A_M["hs"])

# 单通道明细
print(f"\n  L{LINE_DEMO_LASER} 中间角度一字判决明细：")
print(f"  {'距离[m]':>8} {'dL':>8} {'dR':>8} {'判决':>8}")
for q in line_A[LINE_DEMO_LASER]["peaks"]:
    dL = q["dL"] if np.isfinite(q["dL"]) else float("nan")
    dR = q["dR"] if np.isfinite(q["dR"]) else float("nan")
    judge = "丢掉" if q["is_line_ghost"] else "保留"
    print(f"  {q['dist']:8.1f} {dL:8.1f} {dR:8.1f} {judge:>8}")

plot_three_angle_strip(
    A_L, A_M, A_R, line_A, LINE_DEMO_LASER,
    title=f"图 V24-A1  L{LINE_DEMO_LASER} 三角度散点",
    highlight_fake=LINE_FAKE_DIST)
plot_mid_before_after_line(
    A_M, line_A, LINE_DEMO_LASER,
    title=f"图 V24-A2  L{LINE_DEMO_LASER} 中间角度：XM → 一字滤波",
    D=LINE_D_MID)
"""


# ---- 场景 B：稳定鬼影一字滤不掉 ----
CELL_B_DOC = """# V24-B 【证明一字滤波替不了 tcode】Excel 稳定鬼影

配置：

- 三角度都用 **Excel 码**（无 tcode），同一真距离
- 模组内鬼影在三个角度上位置几乎一样 → 空间上「对齐」

期望：一字滤波**几乎不丢鬼影**。稳定串扰必须靠 tcode，不是一字滤波。
"""

CELL_B = r"""# ============================================================================
# V24-B  Excel 稳定鬼影：一字滤波无效
# ============================================================================
B_L = run_angle(firings_excel, ACTIVE_RATIO, LINE_D_LEFT,  radar=False)
B_M = run_angle(firings_excel, ACTIVE_RATIO, LINE_D_MID,   radar=False)
B_R = run_angle(firings_excel, ACTIVE_RATIO, LINE_D_RIGHT, radar=False)

line_B = line_filter_all(B_L["rs"], B_M["rs"], B_R["rs"])
stat_B = summarize_line("B Excel 稳定鬼影", line_B, hs_M=B_M["hs"])

# 对照：中间角度 XM 后还有多少鬼峰
c_xm = count_peak_types(B_M["hs"], B_M["rs"])
print(f"  参考：中间角度 XM 后 纯鬼峰保留 = {c_xm[('纯鬼峰','保留')]}")

plot_three_angle_strip(
    B_L, B_M, B_R, line_B, LINE_DEMO_LASER,
    title=f"图 V24-B1  L{LINE_DEMO_LASER} Excel：鬼影三角度对齐 → 一字滤不掉")
plot_mid_before_after_line(
    B_M, line_B, LINE_DEMO_LASER,
    title=f"图 V24-B2  L{LINE_DEMO_LASER} Excel + 一字：波形几乎不变",
    D=LINE_D_MID)
"""


# ---- 场景 C：对射只出现在中间角度 ----
CELL_C_DOC = """# V24-C 【对射闪一下】只在中间角度出现的对射

配置：

- 三角度都用 tcode；左右**无对射**，中间**有对射**（无 FPGA）
- 对射峰在中间被 XM 留下（4 次同位），但左右没有对应峰

期望：一字滤波清掉对射；真峰保留。
（持续同频对射三角度都在时，仍需 FPGA；一字滤波管的是「闪一下」。）
"""

CELL_C = r"""# ============================================================================
# V24-C  中间角度对射闪现
# ============================================================================
C_L = run_angle(ACTIVE_FR, ACTIVE_RATIO, LINE_D_LEFT,  radar=False)
C_M = run_angle(ACTIVE_FR, ACTIVE_RATIO, LINE_D_MID,   radar=True, use_fpga=False)
C_R = run_angle(ACTIVE_FR, ACTIVE_RATIO, LINE_D_RIGHT, radar=False)

line_C = line_filter_all(C_L["rs"], C_M["rs"], C_R["rs"])
stat_C = summarize_line("C 中间对射闪现", line_C, hs_M=C_M["hs"])

plot_three_angle_strip(
    C_L, C_M, C_R, line_C, LINE_DEMO_LASER,
    title=f"图 V24-C1  L{LINE_DEMO_LASER} 对射只在中间 → 一字应丢掉",
    highlight_fake=RADAR_PHASE_NS * NS * C_LIGHT / 2.0)
plot_mid_before_after_line(
    C_M, line_C, LINE_DEMO_LASER,
    title=f"图 V24-C2  L{LINE_DEMO_LASER} 中间角度对射被一字清掉",
    D=LINE_D_MID)
"""


CELL_COMPARE = r"""# ============================================================================
# V24-4  三场景对比
# ============================================================================
rows = [
    ("A 孤立假峰(tcode+假峰)", stat_A),
    ("B Excel稳定鬼影",       stat_B),
    ("C 中间对射闪现",         stat_C),
]
print("=" * 88)
print(f"一字滤波对比（thr={LINE_THR_M:.1f}m，D≈{LINE_D_MID:.0f}m）")
print("=" * 88)
print(f"  {'场景':<28} {'真留':>6} {'真杀':>6} {'假留':>6} {'假丢':>6} "
      f"{'鬼留':>6} {'鬼丢':>6}")
print("  " + "-" * 78)
for name, st in rows:
    print(f"  {name:<28} {st['true_keep']:>6d} {st['true_drop']:>6d} "
          f"{st['fake_keep']:>6d} {st['fake_drop']:>6d} "
          f"{st['ghost_keep']:>6d} {st['ghost_drop']:>6d}")

fig, axes = plt.subplots(1, 2, figsize=(12, 4.2))
x = np.arange(len(rows))
labels = [n.split()[0] for n, _ in rows]
axes[0].bar(x - 0.2, [st["fake_keep"] for _, st in rows], 0.4,
            color="#c0392b", edgecolor="k", label="假/对射 · 留下")
axes[0].bar(x + 0.2, [st["fake_drop"] for _, st in rows], 0.4,
            color="#27ae60", edgecolor="k", label="假/对射 · 丢掉")
axes[0].set_xticks(x); axes[0].set_xticklabels(labels)
axes[0].set_title("假点 / 对射闪现")
axes[0].legend(fontsize=8); axes[0].grid(alpha=0.25, axis="y")

axes[1].bar(x - 0.2, [st["ghost_keep"] for _, st in rows], 0.4,
            color="#7f8c8d", edgecolor="k", label="模组鬼影 · 留下")
axes[1].bar(x + 0.2, [st["ghost_drop"] for _, st in rows], 0.4,
            color="#f39c12", edgecolor="k", label="模组鬼影 · 丢掉")
axes[1].set_xticks(x); axes[1].set_xticklabels(labels)
axes[1].set_title("模组内鬼影（B 应几乎全留）")
axes[1].legend(fontsize=8); axes[1].grid(alpha=0.25, axis="y")
fig.suptitle("图 V24-4  一字滤波效果对比")
plt.tight_layout(); plt.show()
"""


CELL_SUMMARY = r'''# ============================================================================
# V24 总结
# ============================================================================
print("=" * 82)
print("crosstalk_sim_v24 总结")
print("=" * 82)
print(f"""
【一字滤波】
  连续三角度（三次仿真）。中间峰相对左右最近峰的距离差均 > thr → 丢掉。
  本版 thr = {LINE_THR_M:.1f} m，按同一激光器通道跨角度比较。

【三场景】
  A  中间塞孤立假峰     → 假峰被一字丢掉，真峰留
  B  Excel 稳定鬼影     → 三角度对齐，一字几乎无效（必须靠 tcode）
  C  对射只在中间闪现   → 一字可清；若三角度持续对射仍需 FPGA

【分层滤噪（到 v24）】
  L1 tcode     : 模组内可构造串扰（散开 + 避真峰）
  L2 FPGA 累计 : 持续同频外来对射
  L3 一字滤波  : 角度上孤立的残余假点
  XM           : 单 sync 内 add/max 收割器
""")
'''


# ============================================================================
if not os.path.exists(SRC_NB):
    raise SystemExit(f"找不到 {SRC_NB}，请先生成 v23")

with open(SRC_NB, encoding="utf-8") as f:
    nb23 = json.load(f)

base = []
for c in nb23["cells"]:
    c = dict(c)
    if c["cell_type"] == "code":
        c["outputs"] = []
        c["execution_count"] = None
    base.append(c)

new = [
    md_cell("v24_section", SECTION),
    code_cell("v24_01_params", CELL_PARAMS),
    code_cell("v24_02_core", CELL_CORE),
    code_cell("v24_03_plot", CELL_PLOT),
    md_cell("v24_a_doc", CELL_A_DOC),
    code_cell("v24_a_run", CELL_A),
    md_cell("v24_b_doc", CELL_B_DOC),
    code_cell("v24_b_run", CELL_B),
    md_cell("v24_c_doc", CELL_C_DOC),
    code_cell("v24_c_run", CELL_C),
    code_cell("v24_compare", CELL_COMPARE),
    code_cell("v24_summary", CELL_SUMMARY),
]

nb = {
    "cells": [md_cell("v24_overview", OVERVIEW)] + base + new,
    "metadata": nb23.get("metadata", {}),
    "nbformat": 4,
    "nbformat_minor": 5,
}

with open(OUT_NB, "w", encoding="utf-8") as f:
    json.dump(nb, f, ensure_ascii=False, indent=1)

print(f"已生成 {OUT_NB}")
print(f"  继承 v23：{len(base)} cell")
print(f"  v24 新增：{1 + len(new)} cell")
