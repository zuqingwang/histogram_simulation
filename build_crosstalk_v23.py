# -*- coding: utf-8 -*-
"""
build_crosstalk_v23.py —— 生成 crosstalk_sim_v23.ipynb
========================================================
继承：v21 全部 cell（含 v13 + v20 XM + v21 tcode）。
新增：雷达对射波形链路，按教学顺序展开：

  1) 先加入雷达对射（最坏：同频固定相位）
  2) 只加 tcode：证明能滤模组内鬼影，但滤不掉对射
  3) 再加 FPGA 累计抖动：证明能把对射打散后被 XM 丢掉

每一步都画 hist_add / hist_max×ratio / 滤前滤后 波形图。

缩写：
  XM（XtalkMark，串扰标记）
  FPGA（Field-Programmable Gate Array，现场可编程门阵列）
  TOF（Time of Flight，飞行时间）
"""
import json
import os

SRC_NB = "crosstalk_sim_v21.ipynb"
OUT_NB = "crosstalk_sim_v23.ipynb"
TBL_R15 = os.path.join("docs", "tcode", "tcode_table_v22_r1.5_56ns.py")
TBL_R25 = os.path.join("docs", "tcode", "tcode_table_v22_r2.5_24ns.py")


def load_table(path):
    ns = {}
    with open(path, encoding="utf-8") as f:
        exec(compile(f.read(), path, "exec"), ns)
    return ns["TCODE_TABLE"], ns["TCODE_SEP_NS"], ns["TCODE_BUDGET_NS"]


T15, SEP_NS, B15 = load_table(TBL_R15)
T25, _, B25 = load_table(TBL_R25)
LASERS = sorted({l for l, _ in T15})


def table_literal(table):
    lines = []
    for lid in LASERS:
        kicks = sorted(k for l, k in table if l == lid)
        items = ", ".join(f"({lid},{k}): {table[(lid, k)]:>2d}" for k in kicks)
        lines.append(f"    {items},")
    return "\n".join(lines)


def code_cell(cid, source):
    return {"cell_type": "code", "id": cid, "metadata": {},
            "execution_count": None, "outputs": [],
            "source": source.splitlines(keepends=True)}


def md_cell(cid, source):
    return {"cell_type": "markdown", "id": cid, "metadata": {},
            "source": source.splitlines(keepends=True)}


# ============================================================================
OVERVIEW = f"""# 串扰仿真 v23 —— 雷达对射波形 + tcode + FPGA 累计抖动

> **本 notebook = v21 原封不动 + 末尾追加 v23。**
> v220 只有对射的「落点/聚类」示意，**没有** `hist_add` 滤前/滤后波形。
> v23 专门补上波形，并按下面顺序证明两件事。

## 教学顺序（本版核心）

| 步骤 | 配置 | 要证明什么 | 期望波形 |
|---|---|---|---|
| **S0** | Excel 码 + 雷达对射，无 FPGA | 基线：模组内鬼影 + 对射峰都稳 | 两者 XM 都滤不掉 |
| **S1** | **tcode** + 雷达对射，无 FPGA | tcode 只管模组内 | 鬼影消失，**对射仍在** |
| **S2** | tcode + **FPGA 累计抖动** + 对射 | FPGA 打散对射 | 鬼影和对射都消失 |

## FPGA 累计模型（全模组共用 g[k]）

```
ξ[k]            ∈ {{0, 8, ..., 8(N-1)}} ns     # 本 kick 新增；16 激光器相同
actual[k]       = global_delay + Σ_{{j=0..k}} ξ[j]
```

例：增量 `16ns` 后再 `8ns` → 累计 `16ns, 24ns`。
同 kick 内各激光器共用同一 `actual` → **不破坏 tcode 码差**；外来对射相对各 kick 错位 → 被 XM 滤掉。

## 两套 tcode（都保留）

- `XM_RATIO=1.5` → 预算 **{B15} ns**
- `XM_RATIO=2.5` → 预算 **{B25} ns**（本版演示默认用 1.5）

## 波形模型边界（重要）

本版每条回波仍是落入单个 1ns bin 的理想 δ 脉冲，幅度默认记 1。
`hist_add` 与 `hist_max` 是真实执行的离散直方图运算，但尚未卷积脉冲宽度、
IRF（Instrument Response Function，仪器响应函数）、探测器抖动或噪声。
因此当前 `SEP=12ns` 是用约束留出的峰宽裕度，不是由真实展宽波形仿真出来的。
"""


SECTION = """---
# v23 新增部分

按 **S0 → S1 → S2** 顺序往下跑。每一步都有波形图。
"""


CELL_PARAMS = f"""# ============================================================================
# 1. 参数命名
# ============================================================================
V23_RATIO_MODE = "1.5"          # "1.5" 或 "2.5"
V23_SEP_NS     = {SEP_NS}
V23_DEMO_D     = XM_DEMO_D      # 沿用 v13/v20 演示距离，便于对照
V23_DEMO_LASER = 5              # 单通道解剖用哪支激光器

# 雷达对射（最坏：同频、相位固定）
# 名义落点故意错开真目标（150m），避免叠成「混合峰」糊掉教学对比
RADAR_ENABLE      = True
RADAR_PHASE_NS    = 700.0       # 约 105m 表观距离（≠ V23_DEMO_D）
RADAR_TAG_EMIT    = -1          # 标记外来发射源

# FPGA 累计抖动（全模组共用 g[k]：不破坏 tcode 码差）
FPGA_ENABLE           = False   # S0/S1 关；S2 再开
FPGA_GLOBAL_DELAY_NS  = 8
FPGA_STEP_NS          = 8
FPGA_N_LEVELS         = 8       # 增量 ∈ 8×{{0..7}}
FPGA_SEED             = 7       # 选过：使每激光器 4 个 kick 的累计延时互异 ≥ SEP

print("1. 参数")
print(f"  演示距离 D={{V23_DEMO_D:.0f}}m，解剖激光器 L{{V23_DEMO_LASER}}")
print(f"  XM_RATIO 模式 = {{V23_RATIO_MODE}}")
print(f"  对射相位 = {{RADAR_PHASE_NS:.0f}} ns（表观约 {{RADAR_PHASE_NS*NS*C_LIGHT/2:.1f}} m）")
print(f"  FPGA 累计抖动默认 = {{'开' if FPGA_ENABLE else '关'}}（S2 再打开）")
"""


CELL_TCODE = f"""# ============================================================================
# 2. 导入双 tcode + 构建发光表
# ============================================================================
TCODE_R15 = {{
{table_literal(T15)}
}}
TCODE_R25 = {{
{table_literal(T25)}
}}

def tcode_r15(lid, kick, default=0):
    return TCODE_R15.get((lid, kick), default)

def tcode_r25(lid, kick, default=0):
    return TCODE_R25.get((lid, kick), default)

firings_excel = build_firings_tcode(make_tcode_fn("excel"))
firings_r15   = build_firings_tcode(tcode_r15)
firings_r25   = build_firings_tcode(tcode_r25)

if V23_RATIO_MODE == "2.5":
    ACTIVE_FN, ACTIVE_FR, ACTIVE_RATIO, ACTIVE_BUDGET = tcode_r25, firings_r25, 2.5, {B25}
    ACTIVE_TABLE = TCODE_R25
else:
    ACTIVE_FN, ACTIVE_FR, ACTIVE_RATIO, ACTIVE_BUDGET = tcode_r15, firings_r15, 1.5, {B15}
    ACTIVE_TABLE = TCODE_R15

print(f"2. 启用 tcode：ratio={{ACTIVE_RATIO}}，预算 {{ACTIVE_BUDGET}} ns，max(tx)={{max(ACTIVE_TABLE.values())}} ns")
"""


CELL_INJECT = r"""# ============================================================================
# 3. 雷达对射注入 + FPGA 累计延时
# ============================================================================
KICK_PERIOD_NS = KICK_SPACING / NS

def make_fpga_tables(n_levels=FPGA_N_LEVELS, seed=FPGA_SEED,
                     global_ns=FPGA_GLOBAL_DELAY_NS, step_ns=FPGA_STEP_NS):
    # 全模组共用 g[k]：每 kick 一个增量，再广播到 16 激光器
    # 返回 (increment[laser_idx, kick], cumulative[laser_idx, kick])，单位 ns
    rng = np.random.default_rng(seed)
    inc_k = rng.integers(0, n_levels, size=16) * step_ns
    cum_k = global_ns + np.cumsum(inc_k)
    inc = np.broadcast_to(inc_k, (N_LASERS, 16)).copy()
    cum = np.broadcast_to(cum_k, (N_LASERS, 16)).copy()
    return inc, cum


def fpga_delay_ns(lid, kick, cum_table=None):
    if cum_table is None:
        return 0.0
    return float(cum_table[lid - 1, kick])


def fpga_shot_spreads_ok(cum_table, sep_ns=None):
    # 每个激光器 4 个 shot kick 的累计延时两两 ≥ sep（否则对射双碰，ratio=1.5 留不住）
    if sep_ns is None:
        sep_ns = V23_SEP_NS
    for lid in laser_ids:
        vals = [cum_table[lid - 1, k] for k in SHOT_KICKS[lid]]
        for i in range(len(vals)):
            for j in range(i + 1, len(vals)):
                if abs(vals[i] - vals[j]) < sep_ns:
                    return False
    return True


def inject_radar(recs, D, cum_table=None, phase_ns=None):
    # 向回波列表注入外来雷达
    # 同频对射每个 kick 周期来一次：
    # rec_tof_ns = (phase_ns - actual_delay[l,k]) mod KICK_PERIOD_NS
    # 只有折回后落进接收窗 [0, TOF_WINDOW] 才会被记录
    # 无 FPGA 时 delay=0，4 次 shot 落同一点 → XM 滤不掉
    if phase_ns is None:
        phase_ns = RADAR_PHASE_NS
    extra = []
    for lid in laser_ids:
        for k in SHOT_KICKS[lid]:
            dly = fpga_delay_ns(lid, k, cum_table)
            tof_ns = float(np.mod(phase_ns - dly, KICK_PERIOD_NS))
            if 0.0 <= tof_ns <= (TOF_WINDOW / NS):
                extra.append({
                    "emit_laser": RADAR_TAG_EMIT, "emit_kick": -1,
                    "recv_laser": lid, "recv_kick": k,
                    "target_D": D, "true_tof": 0.0, "t_echo": 0.0,
                    "rec_tof": tof_ns * NS,
                    "rec_dist": tof_ns * NS * C_LIGHT / 2.0,
                    "is_true": False,
                    "negligible": False,
                    "is_radar": True,
                })
    return list(recs) + extra


def run_scenario(fr, ratio, D, use_fpga=False, n_levels=FPGA_N_LEVELS, seed=FPGA_SEED,
                 phase_ns=None):
    # 跑完整链路：发光表 → 模组回波 → 注入对射 → 直方图 → XM
    if phase_ns is None:
        phase_ns = RADAR_PHASE_NS
    cum = None
    inc = None
    if use_fpga:
        inc, cum = make_fpga_tables(n_levels=n_levels, seed=seed)
    with use_firings(fr):
        recs = detect_echoes(D)
    if RADAR_ENABLE:
        recs = inject_radar(recs, D, cum_table=cum, phase_ns=phase_ns)
    hs = build_hists(recs)
    rs = crosstalk_mark_all(hs, ratio)
    return {"recs": recs, "hs": hs, "rs": rs, "inc": inc, "cum": cum,
            "ratio": ratio, "D": D}


def count_peak_types(hs, rs):
    # 按峰内回波来源分类统计
    from collections import Counter
    c = Counter()
    for lid, rr in rs.items():
        for q in rr["peaks"]:
            recs = [r for b in range(q["s"], q["e"] + 1)
                    for r in hs[lid]["src"].get(b, [])]
            n_true = sum(1 for r in recs if r.get("is_true"))
            n_radar = sum(1 for r in recs if r.get("is_radar"))
            n_ghost = len(recs) - n_true - n_radar
            if n_true and not n_radar and not n_ghost:
                kind = "纯真峰"
            elif n_radar and not n_true and not n_ghost:
                kind = "纯对射峰"
            elif n_ghost and not n_true and not n_radar:
                kind = "纯鬼峰"
            else:
                kind = "混合峰"
            act = "丢弃" if q["is_xtalk"] else "保留"
            c[(kind, act)] += 1
    return c


def print_scene_stats(tag, hs, rs):
    c = count_peak_types(hs, rs)
    kinds = ["纯真峰", "纯鬼峰", "纯对射峰", "混合峰"]
    print(f"\n【{tag}】峰级统计（16 激光器合计）")
    print(f"  {'类型':>8} {'保留':>6} {'丢弃':>6}")
    for k in kinds:
        print(f"  {k:>8} {c[(k,'保留')]:>6d} {c[(k,'丢弃')]:>6d}")


# 若默认种子不满足「每激光器 4 shot 累计互异 ≥ SEP」，自动换种子
FPGA_INC, FPGA_CUM = make_fpga_tables()
if not fpga_shot_spreads_ok(FPGA_CUM):
    for _s in range(0, 500):
        FPGA_INC, FPGA_CUM = make_fpga_tables(seed=_s)
        if fpga_shot_spreads_ok(FPGA_CUM):
            FPGA_SEED = _s
            break
print("3. 对射注入与 FPGA 表已就绪（全模组共用 g[k]）")
print(f"  对射名义落点 {RADAR_PHASE_NS:.0f} ns ≈ {RADAR_PHASE_NS*NS*C_LIGHT/2:.1f} m"
      f"（真目标 {V23_DEMO_D:.0f} m，已错开）")
print(f"  FPGA_SEED={FPGA_SEED}，4-shot 打散检查 = {fpga_shot_spreads_ok(FPGA_CUM)}")
print("  共用累计延时前 8 kick: "
      + ", ".join(f"K{k}={FPGA_CUM[0,k]:.0f}" for k in range(8)))
"""


CELL_PLOT_FN = r"""# ============================================================================
# 4. 画波形的公用函数
# ============================================================================
from matplotlib.patches import Patch

# 全项目统一的回波类型视觉编码
V23_TRUE_COLOR  = "#27ae60"     # 真回波：绿色、粗黑边
V23_GHOST_COLOR = "#f5b041"     # 模组内鬼影：橙色、细黑边
V23_RADAR_COLOR = "#c0392b"     # 雷达对射：红色、斜线纹理


def v23_rec_kind(rec):
    if rec.get("is_true"):
        return "true"
    if rec.get("is_radar"):
        return "radar"
    return "ghost"


def v23_rec_style(rec, alpha=0.92):
    kind = v23_rec_kind(rec)
    if kind == "true":
        return dict(facecolor=V23_TRUE_COLOR, edgecolor="black",
                    linewidth=1.8, hatch=None, alpha=alpha)
    if kind == "radar":
        return dict(facecolor=V23_RADAR_COLOR, edgecolor="#7b241c",
                    linewidth=1.1, hatch="///", alpha=alpha)
    return dict(facecolor=V23_GHOST_COLOR, edgecolor="black",
                linewidth=0.8, hatch=None, alpha=alpha)


def v23_rec_label(rec):
    kind = v23_rec_kind(rec)
    if kind == "true":
        return f"真 L{rec['emit_laser']}K{rec['emit_kick']}"
    if kind == "radar":
        return f"对射 K{rec['recv_kick']}"
    return f"鬼 L{rec['emit_laser']}K{rec['emit_kick']}"


def v23_peak_for_bin(rr, b):
    return next((q for q in rr["peaks"] if q["s"] <= b <= q["e"]), None)


def v23_type_legend():
    return [
        Patch(facecolor=V23_TRUE_COLOR, edgecolor="black", linewidth=1.8,
              label="真实信号（绿、粗黑边）"),
        Patch(facecolor=V23_GHOST_COLOR, edgecolor="black", linewidth=0.8,
              label="模组内串扰鬼影（橙）"),
        Patch(facecolor=V23_RADAR_COLOR, edgecolor="#7b241c", hatch="///",
              label="雷达对射（红、斜线）"),
        plt.Line2D([0], [0], marker="x", color="black", ls="none",
                   markersize=9, markeredgewidth=2, label="被 XM 丢掉"),
    ]


def plot_waveform_chain(hs, rs, lid, ratio, title, D, radar_phase_ns=RADAR_PHASE_NS):
    # 单激光器全链路：逐回波堆叠 → add/max×ratio → 滤前/滤后
    h, rr = hs[lid], rs[lid]
    kicks = h["kicks"]
    nz = np.flatnonzero(h["add"] > 0)
    if len(nz) == 0:
        b_lo, b_hi = 0, N_BINS - 1
    else:
        pad = 40
        b_lo = max(0, int(nz.min()) - pad)
        b_hi = min(N_BINS - 1, int(nz.max()) + pad)
    x = np.arange(b_lo, b_hi + 1) * HIST_BIN_M
    true_m = D
    radar_m = radar_phase_ns * NS * C_LIGHT / 2.0

    fig, axes = plt.subplots(3, 1, figsize=(15, 11), sharex=True)

    # (a) 4 次 shot：每条回波单独成柱；同 bin 上下堆叠（对齐 cell12）
    ax = axes[0]
    for i, k in enumerate(kicks):
        y0 = i * 2.2
        ax.axhline(y0, color="0.82", lw=0.7)
        for b in range(b_lo, b_hi + 1):
            recs = [r for r in h["src"].get(b, []) if r["recv_kick"] == k]
            for level, rec in enumerate(recs):
                xc = b * HIST_BIN_M
                st = v23_rec_style(rec)
                ax.add_patch(plt.Rectangle(
                    (xc - 0.55, y0 + level * 0.82), 1.1, 0.76,
                    zorder=3, **st))
                ax.text(xc, y0 + level * 0.82 + 0.38, v23_rec_label(rec),
                        ha="center", va="center", fontsize=5.2, zorder=4)
        ax.text(-0.01, y0 + 0.38, f"shot{i}=K{k}",
                transform=ax.get_yaxis_transform(), ha="right",
                va="center", fontsize=8)
    ax.set_ylim(-0.1, len(kicks) * 2.2 + 0.4)
    ax.set_yticks([])
    ax.set_ylabel("4 次 shot（错开）")
    ax.set_title("(a) 逐回波堆叠柱：柱内写来源；绿=真，橙=模组鬼影，红斜线=雷达对射")
    ax.grid(alpha=0.2, axis="x")

    # (b) hist_add 按来源类型堆叠，再叠加 max 与门限
    ax = axes[1]
    type_count = {k: np.zeros(b_hi - b_lo + 1)
                  for k in ("true", "ghost", "radar")}
    for b in range(b_lo, b_hi + 1):
        for rec in h["src"].get(b, []):
            type_count[v23_rec_kind(rec)][b - b_lo] += echo_amp(rec)
    bottom = np.zeros_like(x, dtype=float)
    for kind, color, label, hatch in [
        ("true", V23_TRUE_COLOR, "真实信号", None),
        ("ghost", V23_GHOST_COLOR, "模组内鬼影", None),
        ("radar", V23_RADAR_COLOR, "雷达对射", "///"),
    ]:
        ax.bar(x, type_count[kind], bottom=bottom, width=HIST_BIN_M * 4,
               color=color, edgecolor="black", linewidth=0.5,
               hatch=hatch, alpha=0.88, label=label)
        bottom += type_count[kind]
    ax.plot(x, h["max"][b_lo:b_hi + 1], drawstyle="steps-mid",
            color="#2874a6", lw=1.6, label="hist_max")
    ax.plot(x, rr["thresh"][b_lo:b_hi + 1], drawstyle="steps-mid",
            color="purple", lw=1.4, ls="--", label=f"max × {ratio}")
    ax.set_ylabel("计数"); ax.grid(alpha=0.2)
    ax.set_title("(b) 滤前 hist_add 按来源着色；紫色门限高于峰顶 → XM 判串扰")
    ax.legend(fontsize=8, ncol=5)

    # (c) XM 后：保留峰仍按逐回波堆叠画；被丢弃峰淡化并打黑叉
    ax = axes[2]
    ymax = 0
    for b in range(b_lo, b_hi + 1):
        q = v23_peak_for_bin(rr, b)
        dropped = bool(q and q["is_xtalk"])
        for level, rec in enumerate(h["src"].get(b, [])):
            xc = b * HIST_BIN_M
            st = v23_rec_style(rec, alpha=0.20 if dropped else 0.92)
            ax.add_patch(plt.Rectangle(
                (xc - 0.55, level), 1.1, 0.86, zorder=3, **st))
            if not dropped:
                ax.text(xc, level + 0.43, v23_rec_label(rec),
                        ha="center", va="center", fontsize=5.2, zorder=4)
            ymax = max(ymax, level + 1)
    for q in rr["peaks"]:
        if q["is_xtalk"]:
            ax.plot(q["dist"], q["add"] + 0.20, "kx", ms=10, mew=2.0)
    ax.axvline(true_m, color="k", ls=":", lw=1.2, label=f"真目标 {true_m:.0f}m")
    ax.axvline(radar_m, color="#c0392b", ls="--", lw=1.2,
               alpha=0.45, label=f"对射名义 {radar_m:.0f}m")
    ax.set_xlabel("记录距离 rec_dist [m]")
    ax.set_ylabel("同距离回波堆叠层数")
    ax.set_ylim(0, ymax + 0.8)
    ax.set_title("(c) XM 后逐回波堆叠柱：被滤峰淡化，黑叉=滤除")
    ax.legend(handles=v23_type_legend(), fontsize=8, ncol=4)
    ax.grid(alpha=0.2, axis="x")

    fig.suptitle(title, fontsize=14, y=0.995)
    plt.tight_layout()
    plt.show()


def plot_grid_before_after(hs, rs, ratio, title, D):
    # 16 宫格：对齐 cell12 的逐回波堆叠柱风格，并保留 cell28 的黑叉滤除标记
    b_lo, b_hi = occupied_range(hs, pad=20)
    nrow = int(np.ceil(N_LASERS / 4))
    fig, axes = plt.subplots(nrow, 4, figsize=(22, 3.2 * nrow), sharex=True)
    axes = np.atleast_2d(axes)
    radar_m = RADAR_PHASE_NS * NS * C_LIGHT / 2.0
    ymax = 1

    for i, lid in enumerate(laser_ids):
        ax = axes[i // 4][i % 4]
        h, rr = hs[lid], rs[lid]
        for b in range(b_lo, b_hi + 1):
            q = v23_peak_for_bin(rr, b)
            dropped = bool(q and q["is_xtalk"])
            for level, rec in enumerate(h["src"].get(b, [])):
                xc = b * HIST_BIN_M
                st = v23_rec_style(rec, alpha=0.20 if dropped else 0.92)
                ax.add_patch(plt.Rectangle(
                    (xc - 0.55, level), 1.1, 0.86, zorder=3, **st))
                if not dropped:
                    ax.text(xc, level + 0.43, v23_rec_label(rec),
                            ha="center", va="center", fontsize=4.7, zorder=4)
                ymax = max(ymax, level + 1)
        for q in rr["peaks"]:
            if q["is_xtalk"]:
                ax.plot(q["dist"], q["add"] + 0.18, "kx",
                        ms=8, mew=1.7, zorder=6)
        ax.axvline(D, color="k", ls=":", lw=0.9, alpha=0.7)
        ax.axvline(radar_m, color="#c0392b", ls="--", lw=0.8, alpha=0.30)
        n_drop = sum(1 for q in rr["peaks"] if q["is_xtalk"])
        n_keep = sum(1 for q in rr["peaks"] if not q["is_xtalk"])
        ax.set_title(f"L{lid}: 保留{n_keep}/丢掉{n_drop}", fontsize=9)
        ax.tick_params(labelsize=7); ax.grid(alpha=0.16, axis="x")
    for j in range(N_LASERS, nrow * 4):
        axes[j // 4][j % 4].axis("off")
    for i in range(N_LASERS):
        axes[i // 4][i % 4].set_ylim(0, ymax + 0.8)
    fig.text(0.5, 0.01,
             f"逐回波堆叠柱；绿=真信号，橙=模组鬼影，红斜线=雷达对射；"
             f"黑叉=被 XM 滤除；黑点线=真目标 {D:.0f}m；XM_RATIO={ratio}",
             ha="center", fontsize=11)
    fig.legend(handles=v23_type_legend(), loc="upper center",
               ncol=4, fontsize=9)
    plt.suptitle(title, fontsize=14, y=0.995)
    plt.tight_layout(rect=[0, 0.03, 1, 0.93])
    plt.show()

print("4. 绘图函数已加载")
"""


# ---- S0 ----
CELL_S0_DOC = """# 5. 【S0】基线：Excel 码 + 雷达对射（无 FPGA）

期望：

- 模组内鬼影：每个激光器 4 次 kick 码相同 → 鬼影 4 次落同一 bin → XM 滤不掉
- 雷达对射：无 FPGA 抖动 → 4 次也落同一点 → XM 同样滤不掉
"""

CELL_S0 = r"""# ============================================================================
# 5. S0 运行 + 波形
# ============================================================================
S0 = run_scenario(firings_excel, ACTIVE_RATIO, V23_DEMO_D, use_fpga=False)
print_scene_stats("S0 Excel + 对射，无 FPGA", S0["hs"], S0["rs"])

plot_waveform_chain(
    S0["hs"], S0["rs"], V23_DEMO_LASER, ACTIVE_RATIO,
    title=f"图 V23-S0a  基线波形 —— L{V23_DEMO_LASER}，Excel 码 + 对射，无 FPGA",
    D=V23_DEMO_D)

plot_grid_before_after(
    S0["hs"], S0["rs"], ACTIVE_RATIO,
    title=f"图 V23-S0b  16 宫格基线 —— 鬼影和对射都应还在（XM_RATIO={ACTIVE_RATIO}）",
    D=V23_DEMO_D)
"""


# ---- S1 ----
CELL_S1_DOC = """# 6. 【S1】只加 tcode（仍无 FPGA）

期望：

- 模组内鬼影：被 tcode 打散 → `add/max≈1` → **被 XM 丢掉**
- 雷达对射：对方不共享我们的 tcode，相对落点仍固定 → **仍然留下**

> 这一步专门证明：**tcode 不能替代 FPGA 去打对射。**
"""

CELL_S1 = r"""# ============================================================================
# 6. S1 运行 + 波形
# ============================================================================
S1 = run_scenario(ACTIVE_FR, ACTIVE_RATIO, V23_DEMO_D, use_fpga=False)
print_scene_stats(f"S1 tcode({ACTIVE_RATIO}) + 对射，无 FPGA", S1["hs"], S1["rs"])

plot_waveform_chain(
    S1["hs"], S1["rs"], V23_DEMO_LASER, ACTIVE_RATIO,
    title=f"图 V23-S1a  只加 tcode —— L{V23_DEMO_LASER}：鬼影应散开被丢，对射峰应还在",
    D=V23_DEMO_D)

plot_grid_before_after(
    S1["hs"], S1["rs"], ACTIVE_RATIO,
    title=f"图 V23-S1b  16 宫格：tcode 清掉模组内鬼影后，红虚线处对射峰仍在",
    D=V23_DEMO_D)
"""


# ---- S2 ----
CELL_S2_DOC = """# 7. 【S2】tcode + FPGA 累计抖动

期望：

- 模组内鬼影：仍由 tcode 管住（**FPGA 全模组共用 g[k]**，同 kick 内码差不变）
- 雷达对射：累计延时让同一激光器 4 次落点两两错开 ≥ SEP → `add/max≈1` → **被 XM 丢掉**

FPGA 模型：`actual[k] = 8ns + Σ_{j≤k} ξ[j]`，ξ 对 16 激光器相同。
"""

CELL_S2 = r"""# ============================================================================
# 7. S2 运行 + FPGA 轨迹 + 波形
# ============================================================================
S2 = run_scenario(ACTIVE_FR, ACTIVE_RATIO, V23_DEMO_D,
                  use_fpga=True, n_levels=FPGA_N_LEVELS, seed=FPGA_SEED)
print_scene_stats(f"S2 tcode({ACTIVE_RATIO}) + FPGA 累计 + 对射", S2["hs"], S2["rs"])

# FPGA 增量 / 累计轨迹（全模组共用 → 各激光器相同）
fig, axes = plt.subplots(2, 1, figsize=(14, 7), sharex=True)
axes[0].step(range(16), S2["inc"][0], where="mid", lw=2.0, color="#1a5276")
axes[1].step(range(16), S2["cum"][0], where="mid", lw=2.0, color="#1a5276")
axes[0].set_ylabel("本 kick 新增 ξ [ns]")
axes[0].set_title(f"图 V23-S2a  FPGA 随机增量（seed={FPGA_SEED}，全模组共用）")
axes[0].grid(alpha=0.25)
axes[1].set_xlabel("全局 kick 序号"); axes[1].set_ylabel("累计延时 d [ns]")
axes[1].set_title("图 V23-S2b  实际延时 = 8ns 整体补偿 + 前面所有增量之和")
axes[1].grid(alpha=0.25)
plt.tight_layout(); plt.show()

# 对射 4 次落点被推开
lid = V23_DEMO_LASER
kicks = SHOT_KICKS[lid]
pos = np.mod(
    np.array([RADAR_PHASE_NS - S2["cum"][lid - 1, k] for k in kicks]),
    KICK_PERIOD_NS)
true_ns = 2.0 * V23_DEMO_D / C_LIGHT / NS
fig, ax = plt.subplots(figsize=(12, 3.2))
for i, (k, p) in enumerate(zip(kicks, pos)):
    ax.vlines(p, 0, 1, lw=3, color=V23_RADAR_COLOR,
              label=f"对射 K{k}: {p:.0f}ns")
    ax.scatter([p], [1], s=70, color=V23_RADAR_COLOR,
               edgecolor="#7b241c", hatch="///")
ax.axvline(true_ns, color="black", ls=":", lw=2.0,
           label=f"真实信号 {true_ns:.1f}ns")
ax.axvspan(true_ns - V23_SEP_NS, true_ns + V23_SEP_NS,
           color="black", alpha=0.10, label=f"真峰 ±{V23_SEP_NS}ns")
ax.set_ylim(0, 1.4); ax.set_xlabel("对射在直方图中的记录位置 [ns]")
ax.set_title(f"图 V23-S2c  L{lid}：红色=对射，黑点线=真实信号；同时检查是否撞真峰")
ax.legend(ncol=4); ax.grid(alpha=0.2)
plt.tight_layout(); plt.show()

plot_waveform_chain(
    S2["hs"], S2["rs"], V23_DEMO_LASER, ACTIVE_RATIO,
    title=f"图 V23-S2d  tcode + FPGA —— L{V23_DEMO_LASER}：鬼影和对射都应被 XM 清掉",
    D=V23_DEMO_D)

plot_grid_before_after(
    S2["hs"], S2["rs"], ACTIVE_RATIO,
    title=f"图 V23-S2e  16 宫格：双手段叠加后，只剩真目标（黑点线）",
    D=V23_DEMO_D)
"""


CELL_COMPARE = r"""# ============================================================================
# 8. 三步对比总表 + 残留柱状图
# ============================================================================
SCENES = [
    ("S0 Excel+对射", S0),
    (f"S1 tcode({ACTIVE_RATIO})+对射", S1),
    (f"S2 tcode+FPGA+对射", S2),
]

rows = []
for name, sc in SCENES:
    c = count_peak_types(sc["hs"], sc["rs"])
    rows.append({
        "name": name,
        "true_keep": c[("纯真峰", "保留")],
        "true_kill": c[("纯真峰", "丢弃")] + c[("混合峰", "丢弃")],
        "ghost_keep": c[("纯鬼峰", "保留")],
        "ghost_drop": c[("纯鬼峰", "丢弃")],
        "radar_keep": c[("纯对射峰", "保留")],
        "radar_drop": c[("纯对射峰", "丢弃")],
    })

print("=" * 88)
print(f"三步对比（D={V23_DEMO_D:.0f}m，XM_RATIO={ACTIVE_RATIO}）")
print("=" * 88)
print(f"  {'场景':<28} {'真保留':>6} {'真误杀':>6} {'鬼保留':>6} {'鬼丢掉':>6} "
      f"{'对射保留':>8} {'对射丢掉':>8}")
print("  " + "-" * 78)
for r in rows:
    print(f"  {r['name']:<28} {r['true_keep']:>6d} {r['true_kill']:>6d} "
          f"{r['ghost_keep']:>6d} {r['ghost_drop']:>6d} "
          f"{r['radar_keep']:>8d} {r['radar_drop']:>8d}")

fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
labels = [r["name"].split()[0] for r in rows]
axes[0].bar(labels, [r["ghost_keep"] for r in rows], color="#7f8c8d", edgecolor="k")
axes[0].set_title("模组内纯鬼峰 · 残留（保留）")
axes[0].set_ylabel("峰数")
axes[1].bar(labels, [r["radar_keep"] for r in rows], color="#c0392b", edgecolor="k")
axes[1].set_title("雷达对射峰 · 残留（保留）")
axes[2].bar(labels, [r["true_keep"] for r in rows], color="#1e8449", edgecolor="k")
axes[2].set_title("真目标峰 · 存活")
for ax in axes:
    ax.grid(alpha=0.25, axis="y")
    ax.tick_params(axis="x", rotation=15)
plt.suptitle("图 V23-对比  S0→S1→S2：tcode 清鬼影，FPGA 清对射，真目标始终在",
             fontsize=13)
plt.tight_layout()
plt.show()
"""


CELL_SUMMARY = f"""# ============================================================================
# 9. v23 总结
# ============================================================================
print("=" * 82)
print("crosstalk_sim_v23 总结")
print("=" * 82)
print(f'''
【回答你的问题】
  v220 有对射「落点/聚类」图，但没有 hist_add 滤前/滤后波形。
  v23 从一开始就加入雷达对射，并画出完整波形链路。

【证明链】
  S0  Excel + 对射          → 鬼影和对射都稳，XM 无效
  S1  + tcode               → 鬼影被滤，对射仍在（tcode 管不了对射）
  S2  + FPGA 累计抖动       → 对射落点被推开，XM 一并清掉

【手段分工】
  tcode  : 模组内串扰（同 kick / 跨 kick / 自身混叠），且避真峰
  FPGA   : 外来雷达对射（累计随机增量）
  XM     : 收割器（看 add/max）

【参数】
  当前演示 XM_RATIO = {{ACTIVE_RATIO}}，tcode 预算 {{ACTIVE_BUDGET}} ns
  另一套（切换 V23_RATIO_MODE）：1.5→{B15}ns / 2.5→{B25}ns

【下一步可做】
  - 一字滤波（连续三角度）
  - 距离扫描下的对射残留曲线
  - 把回波展宽（脉宽）加进仿真，复核 SEP
''')
"""


# v23 同时把继承来的 v21 cell28 改成项目统一的“逐回波堆叠柱 + 黑叉滤除”风格。
# 只替换该 cell 的绘图区；仿真、统计和用户参数全部保持不变。
V21_GRID_STYLE = r"""b_lo, b_hi = occupied_range(hs21, pad=15)
nrow = int(np.ceil(N_LASERS / 4))
fig, axes = plt.subplots(nrow, 4, figsize=(22, 3.2 * nrow), sharex=True)
axes = np.atleast_2d(axes)
ymax = 1

for i, lid in enumerate(laser_ids):
    ax = axes[i // 4][i % 4]
    h, rr = hs21[lid], rs21[lid]
    n_hit = n_miss = n_kill = 0
    for b in range(b_lo, b_hi + 1):
        q = next((qq for qq in rr["peaks"] if qq["s"] <= b <= qq["e"]), None)
        dropped = bool(q and q["is_xtalk"])
        for level, rec in enumerate(h["src"].get(b, [])):
            xc = b * HIST_BIN_M
            is_true = rec.get("is_true")
            color = "#27ae60" if is_true else "#f5b041"
            lw = 1.8 if is_true else 0.8
            alpha = 0.20 if dropped else 0.92
            ax.add_patch(plt.Rectangle(
                (xc - 0.55, level), 1.1, 0.86,
                facecolor=color, edgecolor="black", linewidth=lw,
                alpha=alpha, zorder=3))
            if not dropped:
                label = (f"真 L{rec['emit_laser']}K{rec['emit_kick']}"
                         if is_true else
                         f"鬼 L{rec['emit_laser']}K{rec['emit_kick']}")
                ax.text(xc, level + 0.43, label, ha="center", va="center",
                        fontsize=4.7, zorder=4)
            ymax = max(ymax, level + 1)
    for q in rr["peaks"]:
        kind, _ = peak_truth(h, q)
        if q["is_xtalk"]:
            ax.plot(q["dist"], q["add"] + 0.18, "kx",
                    ms=8, mew=1.7, zorder=6)
            if kind == "纯鬼峰":
                n_hit += 1
            else:
                n_kill += 1
                ax.plot(q["dist"], q["add"] + 0.18, "o",
                        mfc="none", mec="red", ms=13, mew=2, zorder=7)
        elif kind == "纯鬼峰":
            n_miss += 1
    ax.axvline(D_V21, color="k", ls=":", lw=1.0, alpha=0.6)
    ax.set_title(f"L{lid}：滤掉鬼影{n_hit}、残留{n_miss}"
                 + (f"、!误杀{n_kill}" if n_kill else ""), fontsize=9)
    ax.tick_params(labelsize=7); ax.grid(alpha=0.16, axis="x")

for j in range(N_LASERS, nrow * 4):
    axes[j // 4][j % 4].axis("off")
for i in range(N_LASERS):
    axes[i // 4][i % 4].set_ylim(0, ymax + 0.8)

from matplotlib.patches import Patch
legend_items = [
    Patch(facecolor="#27ae60", edgecolor="black", linewidth=1.8,
          label="真实信号（绿、粗黑边）"),
    Patch(facecolor="#f5b041", edgecolor="black", linewidth=0.8,
          label="模组内串扰鬼影（橙）"),
    plt.Line2D([0], [0], marker="x", color="black", ls="none",
               markersize=9, markeredgewidth=2, label="被 XM 滤除"),
    plt.Line2D([0], [0], marker="o", markerfacecolor="none",
               markeredgecolor="red", ls="none", markersize=10,
               markeredgewidth=2, label="误杀"),
]
fig.legend(handles=legend_items, loc="upper center", ncol=4, fontsize=9)
fig.text(0.5, 0.015,
         f"记录距离 rec_dist [m]；逐回波堆叠柱；黑点线=物体真实距离 {D_V21:.0f}m；"
         f"黑叉=被 XM 丢弃；红圈=误杀", ha="center", fontsize=11)
fig.text(0.008, 0.5, "同一记录距离处的回波堆叠层数",
         va="center", rotation="vertical", fontsize=11)
plt.suptitle(f"图G  v21 编码下的 XM 滤除 —— D={D_V21:.0f}m，XM_RATIO={XM_RATIO}",
             fontsize=14, y=0.995)
plt.tight_layout(rect=[0.015, 0.035, 1, 0.93])
plt.show()

"""


def restyle_v21_grid_cell(source):
    start = "b_lo, b_hi = occupied_range(hs21, pad=15)"
    end = "print_eval(st21,"
    if start not in source or end not in source:
        raise RuntimeError("找不到 v21 cell28 绘图区，无法套用统一回波风格")
    prefix = source.split(start, 1)[0]
    suffix = end + source.split(end, 1)[1]
    return prefix + V21_GRID_STYLE + suffix


# ============================================================================
if not os.path.exists(SRC_NB):
    raise SystemExit(f"找不到 {SRC_NB}")

with open(SRC_NB, encoding="utf-8") as f:
    nb21 = json.load(f)

base = []
for c in nb21["cells"]:
    c = dict(c)
    if c["cell_type"] == "code":
        c["outputs"] = []
        c["execution_count"] = None
    if c.get("id") == "v21_tc_figgrid":
        c["source"] = restyle_v21_grid_cell("".join(c["source"])).splitlines(keepends=True)
    base.append(c)

new = [
    md_cell("v23_section", SECTION),
    code_cell("v23_01_params", CELL_PARAMS),
    code_cell("v23_02_tcode", CELL_TCODE),
    code_cell("v23_03_inject", CELL_INJECT),
    code_cell("v23_04_plotfn", CELL_PLOT_FN),
    md_cell("v23_s0_doc", CELL_S0_DOC),
    code_cell("v23_s0_run", CELL_S0),
    md_cell("v23_s1_doc", CELL_S1_DOC),
    code_cell("v23_s1_run", CELL_S1),
    md_cell("v23_s2_doc", CELL_S2_DOC),
    code_cell("v23_s2_run", CELL_S2),
    code_cell("v23_compare", CELL_COMPARE),
    code_cell("v23_summary", CELL_SUMMARY),
]

nb = {
    "cells": [md_cell("v23_overview", OVERVIEW)] + base + new,
    "metadata": nb21.get("metadata", {}),
    "nbformat": 4,
    "nbformat_minor": 5,
}

with open(OUT_NB, "w", encoding="utf-8") as f:
    json.dump(nb, f, ensure_ascii=False, indent=1)

print(f"已生成 {OUT_NB}")
print(f"  继承 v21：{len(base)} cell")
print(f"  v23 新增：{1 + len(new)} cell")
