# -*- coding: utf-8 -*-
"""
生成 crosstalk_sim_v220.ipynb。

结构：
  v220 = v21 全部内容 + 重新整理的 v220 新增部分。
  不继承 v22 新增部分，避免重复和逻辑混乱。

缩写：
  XM（XtalkMark，串扰标记）
  FPGA（Field-Programmable Gate Array，现场可编程门阵列）
  TOF（Time of Flight，飞行时间）
"""
import json
import os

SRC_NB = "crosstalk_sim_v21.ipynb"
OUT_NB = "crosstalk_sim_v220.ipynb"
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
    return {
        "cell_type": "code",
        "id": cid,
        "metadata": {},
        "execution_count": None,
        "outputs": [],
        "source": source.splitlines(keepends=True),
    }


def md_cell(cid, source):
    return {
        "cell_type": "markdown",
        "id": cid,
        "metadata": {},
        "source": source.splitlines(keepends=True),
    }


OVERVIEW = f"""# 串扰仿真 v220 —— 双 tcode + 累计 FPGA 延时 + 雷达对射

> v220 从 v21 重新向后扩展，不继承 v22 的新增 cell，避免重复与逻辑交叉。
> v13、v20、v21 的既有内容完整保留。

## v220 新增部分的阅读顺序

1. 参数与变量命名
2. 导入发光时序并画图
3. 导入两套 tcode 并画码矩阵
4. 验证“散开 + 避真峰”并画约束余量
5. 距离扫描并画 tcode 滤除效果
6. 定义 FPGA 累计随机延时模型
7. 画随机增量与累计延时轨迹
8. 画一次雷达对射：真实信号 + 对射落点，看两者是否重叠
9. 解析 + 蒙特卡洛：FPGA 随机空间 vs 对射残留概率
10. 总结

## 两套保留方案

- `XM_RATIO=1.5`：tcode 预算 {B15}ns，模组内鬼影扫描残留为 0
- `XM_RATIO=2.5`：tcode 预算 {B25}ns，模组内鬼影扫描残留为 0

## FPGA 延时模型（本版纠正）

随机数不是每个 kick 的绝对延时，而是**增量**：

`actual_delay[l,k] = global_delay + sum(random_increment[l,0:k+1])`

例如增量依次为 `16ns, 8ns`，实际累计延时就是 `16ns, 24ns`。
"""


SECTION = """---
# v220 新增部分

下面严格按顺序执行；每个分析块后面都有对应图。
"""


CELL_1 = f"""# ============================================================================
# 1. 参数与变量命名
# ============================================================================
V220_RATIO_MODE = "1.5"       # "1.5" 或 "2.5"
V220_SEP_NS = {SEP_NS}         # 鬼影间隔、鬼影离真峰的最小距离

# FPGA 累计随机延时
FPGA_GLOBAL_DELAY_NS = 8      # 所有 kick 开始前的整体固定补偿
FPGA_STEP_NS = 8              # 每个随机增量的步长
FPGA_N_LEVELS = 3             # 增量可取 0, 8, ..., 8*(N-1)
FPGA_MC_TRIALS = 20000        # 雷达对射蒙卡次数
FPGA_RANDOM_SEED = 220

# 对射模型
RADAR_PHASE_NS = 1000         # 外来雷达名义落点，放在 2000ns 窗中间
RADAR_PULSE_SEP_NS = V220_SEP_NS

print("1. 参数与变量命名")
print(f"  V220_RATIO_MODE       = {{V220_RATIO_MODE}}")
print(f"  V220_SEP_NS           = {{V220_SEP_NS}} ns")
print(f"  FPGA_GLOBAL_DELAY_NS  = {{FPGA_GLOBAL_DELAY_NS}} ns")
print(f"  FPGA 随机增量集合      = 8×[0, ..., {{FPGA_N_LEVELS-1}}] ns")
print("  actual_delay[l,k] = global_delay + Σ increment[l,j], j=0..k")
"""


CELL_2 = r"""# ============================================================================
# 2. 导入并可视化当前发光时序（栅格风格，与 v13 图 1 统一）
# ============================================================================
# x 轴 = kick 序号(0~15)，y 轴 = 激光器，格内 = 该 kick 的 tx_trig_dly [ns]。
# 后面第 3 节的 tcode 码矩阵用同一套 (kick × 激光器) 坐标，两张图可以直接对照。
def plot_kick_grid(code_of, title, ax=None, unit="ns", show_text=True):
    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 9))
    for lid in laser_ids:
        for k in SHOT_KICKS[lid]:
            ax.add_patch(plt.Rectangle((k - 0.46, lid - 0.40), 0.92, 0.80,
                                       facecolor=laser_color(lid), edgecolor="k",
                                       linewidth=0.4, alpha=0.95, zorder=2))
            if show_text:
                ax.text(k, lid, f"{code_of(lid, k)}", fontsize=7,
                        ha="center", va="center", color="k", zorder=3)
    ax.set_xticks(range(16)); ax.set_yticks(laser_ids)
    ax.set_xlim(-0.6, 15.6); ax.set_ylim(0.4, N_LASERS + 0.6)
    ax.set_xlabel("kick 序号（第几个 kick 发光；不看时间轴）")
    ax.set_ylabel("激光器编号")
    ax.set_title(title)
    ax.grid(alpha=0.25, zorder=0); ax.set_axisbelow(True)
    return ax


plot_kick_grid(lambda lid, k: make_tcode_fn("excel")(lid, k, 0),
               f"图 V220-1  当前发光时序栅格（{N_LASERS} 激光器 × 16 kick，"
               f"格内 = Excel tx_trig_dly [ns]，色 = 激光器）")
plt.tight_layout()
plt.show()

print("2. 时序导入完成：")
for lid in laser_ids:
    print(f"  L{lid:<2d}: {SHOT_KICKS[lid]}")
"""


CELL_3 = f"""# ============================================================================
# 3. 导入双 tcode，并画码矩阵
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

firings_r15 = build_firings_tcode(tcode_r15)
firings_r25 = build_firings_tcode(tcode_r25)
firings_excel = build_firings_tcode(make_tcode_fn("excel"))

if V220_RATIO_MODE == "2.5":
    ACTIVE_TCODE, ACTIVE_CODE_FN = TCODE_R25, tcode_r25
    ACTIVE_FIRINGS, ACTIVE_RATIO, ACTIVE_BUDGET_NS = firings_r25, 2.5, {B25}
else:
    ACTIVE_TCODE, ACTIVE_CODE_FN = TCODE_R15, tcode_r15
    ACTIVE_FIRINGS, ACTIVE_RATIO, ACTIVE_BUDGET_NS = firings_r15, 1.5, {B15}

fig, axes = plt.subplots(1, 2, figsize=(17, 6.5))
for ax, (name, table, budget) in zip(
        axes,
        [("ratio=1.5", TCODE_R15, {B15}), ("ratio=2.5", TCODE_R25, {B25})]):
    mat = np.full((N_LASERS, 16), np.nan)
    for i, lid in enumerate(laser_ids):
        for k in SHOT_KICKS[lid]:
            mat[i, k] = table[(lid, k)]
    im = ax.imshow(mat, origin="lower", aspect="auto", cmap="viridis",
                   extent=[-0.5, 15.5, 0.5, 16.5], vmin=0, vmax={B15})
    for i, lid in enumerate(laser_ids):
        for k in SHOT_KICKS[lid]:
            ax.text(k, lid, str(table[(lid, k)]), ha="center", va="center",
                    fontsize=6.5, color="white" if table[(lid, k)] < {B15}*0.55 else "black")
    ax.set_xticks(range(16)); ax.set_yticks(laser_ids)
    ax.set_xlabel("kick"); ax.set_ylabel("激光器")
    ax.set_title(f"{{name}} 专用码表（预算 {{budget}}ns）")
    fig.colorbar(im, ax=ax, fraction=0.045, label="tx_trig_dly [ns]")
plt.suptitle("图 V220-2  两套 tcode 码矩阵", fontsize=14)
plt.tight_layout()
plt.show()

print(f"3. 当前启用：ratio={{ACTIVE_RATIO}}, tcode 预算={{ACTIVE_BUDGET_NS}}ns")
"""


CELL_4 = r"""# ============================================================================
# 4. tcode 约束验证：散开 + 避真峰，并画余量
# ============================================================================
def v220_ghost_classes(gap=CROSSTALK_MAX_GAP, dks=(0, 1)):
    rows = []
    for recv in laser_ids:
        for emit in laser_ids:
            if abs(emit - recv) > gap:
                continue
            for dk in dks:
                if emit == recv and dk == 0:
                    continue
                pairs = [(kr - dk, kr) for kr in SHOT_KICKS[recv]
                         if kr - dk in SHOT_KICKS[emit]]
                if len(pairs) >= 2:
                    kind = ("自身混叠" if emit == recv else
                            ("同kick串扰" if dk == 0 else "跨kick串扰"))
                    rows.append((emit, recv, dk, kind, pairs))
    return rows

def constraint_margins(code_fn):
    spread, avoid = [], []
    for emit, recv, dk, kind, pairs in v220_ghost_classes():
        diffs = sorted(code_fn(emit, ke, 0) - code_fn(recv, kr, 0)
                       for ke, kr in pairs)
        spread.append((kind, min(np.diff(diffs))))
        if dk == 0 and emit != recv:
            avoid.extend(abs(code_fn(emit, ke, 0) - code_fn(recv, kr, 0))
                         for ke, kr in pairs)
    return spread, avoid

fig, axes = plt.subplots(2, 2, figsize=(15, 10))
for col, (name, fn) in enumerate([("ratio=1.5 码表", tcode_r15),
                                  ("ratio=2.5 码表", tcode_r25)]):
    spread, avoid = constraint_margins(fn)
    spread_values = np.array([x[1] for x in spread])
    avoid_values = np.array(avoid)

    ax = axes[0, col]
    ax.hist(spread_values, bins=np.arange(0, max(spread_values)+3, 2),
            color="#2874a6", edgecolor="white")
    ax.axvline(V220_SEP_NS, color="red", ls="--", lw=2, label=f"要求 ≥{V220_SEP_NS}ns")
    ax.set_title(f"{name}：同类鬼影落点的最小间隔\n最小值={spread_values.min():.0f}ns")
    ax.set_xlabel("同一鬼影类内的最小间隔 [ns]"); ax.set_ylabel("鬼影类数量")
    ax.legend(); ax.grid(alpha=0.2)

    ax = axes[1, col]
    ax.hist(avoid_values, bins=np.arange(0, max(avoid_values)+3, 2),
            color="#1e8449", edgecolor="white")
    ax.axvline(V220_SEP_NS, color="red", ls="--", lw=2, label=f"要求 ≥{V220_SEP_NS}ns")
    ax.set_title(f"{name}：鬼影离真峰(0点)的距离\n最小值={avoid_values.min():.0f}ns")
    ax.set_xlabel("|鬼影码差| [ns]"); ax.set_ylabel("shot 数量")
    ax.legend(); ax.grid(alpha=0.2)

plt.suptitle("图 V220-3  tcode 两条硬约束：鬼影彼此散开，并且不撞真峰", fontsize=14)
plt.tight_layout()
plt.show()
"""


CELL_5 = r"""# ============================================================================
# 5. 加 tcode 后的效果：距离扫描
# ============================================================================
V220_DISTS = np.arange(5.0, 600.1, 5.0)
V220_CASES = [
    ("Excel + ratio1.5", firings_excel, 1.5, "#7f8c8d"),
    ("24ns码 + ratio2.5", firings_r25, 2.5, "#c0392b"),
    ("56ns码 + ratio1.5", firings_r15, 1.5, "#1e8449"),
]

def v220_sweep(firing_table, ratio):
    out = {k: [] for k in ("ghost_before", "ghost_after", "true_killed")}
    for dist in V220_DISTS:
        with use_firings(firing_table):
            recs = detect_echoes(dist)
        hs = build_hists(recs)
        rs = crosstalk_mark_all(hs, ratio)
        stat = evaluate(hs, rs)["peak"]
        out["ghost_before"].append(sum(stat[("纯鬼峰", a)] for a in ("保留", "丢弃")))
        out["ghost_after"].append(stat[("纯鬼峰", "保留")])
        out["true_killed"].append(sum(stat[(k, "丢弃")] for k in ("纯真峰", "混合峰")))
    return {k: np.asarray(v) for k, v in out.items()}

V220_SWEEP = {name: v220_sweep(fr, ratio) for name, fr, ratio, _ in V220_CASES}

fig, axes = plt.subplots(2, 1, figsize=(15, 9), sharex=True)
for name, _, ratio, color in V220_CASES:
    x = V220_SWEEP[name]
    residual = x["ghost_after"] / np.maximum(x["ghost_before"], 1) * 100
    axes[0].plot(V220_DISTS, residual, lw=2.2, color=color,
                 label=f"{name}（总残留 {x['ghost_after'].sum()}）")
    axes[1].plot(V220_DISTS, x["true_killed"], lw=2.0, color=color, label=name)
axes[0].set_ylabel("鬼影残留率 [%]"); axes[0].set_ylim(-2, 105)
axes[0].set_title("图 V220-4a  tcode + XM：鬼影残留率")
axes[0].legend(); axes[0].grid(alpha=0.25)
axes[1].set_xlabel("真实目标距离 [m]"); axes[1].set_ylabel("真目标误杀数")
axes[1].set_title("图 V220-4b  避真约束验收：误杀应恒为 0")
axes[1].legend(); axes[1].grid(alpha=0.25)
plt.tight_layout()
plt.show()
"""


FPGA_DOC = r"""# 6. FPGA 累计随机延时模型（纠正 v22）

对激光器 $l$ 和全局 kick $k$：

$$\xi_{l,k}\in\{0,8,16,\ldots,8(N-1)\}\text{ ns}$$

$$d_{l,k}=d_{\rm global}+\sum_{j=0}^{k}\xi_{l,j}$$

- $\xi$ 是每个 kick 新产生的**随机增量**
- $d$ 才是该 kick 实际使用的**累计延时**
- 例如增量为 `16ns, 8ns`，累计延时为 `16ns, 24ns`
- 固定的整体补偿 `d_global=8ns` 只平移所有落点，不改变多 shot 之间的间隔

外来雷达与我们同频、固定相位（最坏对射）时，它的脉冲每 `KICK_SPACING` 到一次，
所以在接收直方图中的记录位置要对 kick 周期取模：

$$t_{\rm radar,recorded}(l,k)=\left(\phi_{\rm radar}-d_{l,k}\right)\bmod T_{\rm kick}$$

只有当该值落在 $[0,\,T_{\rm window}]$ 内才会被记录。

因此累计随机 walk 会把原本固定的对射峰逐 kick 推向不同位置。

## 关键量：不是随机空间大小，而是相邻 shot 的延时差

XM 比的是同一个激光器 4 次 shot 之间的落点。落点差 = **累计延时之差**：

$$\Delta_{i}=d_{l,k_{i+1}}-d_{l,k_i}=\sum_{j=k_i+1}^{k_{i+1}}\xi_{l,j}$$

对**相邻 kick**（例如 L5 的 K0→K1），这个和只有**一项** $\xi$。
所以「累计 16 个 kick 的随机空间有多大」并不决定滤除效果，
决定性的是单个增量能不能一步跨过峰宽 `SEP`。
"""


CELL_7 = r"""# ============================================================================
# 7. 画 FPGA 随机增量与累计延时轨迹
# ============================================================================
def generate_fpga_cumulative(n_levels, rng, n_lasers=N_LASERS, n_kicks=16,
                             step_ns=None):
    step_ns = FPGA_STEP_NS if step_ns is None else step_ns
    increments = rng.integers(0, n_levels, size=(n_lasers, n_kicks)) * step_ns
    cumulative = FPGA_GLOBAL_DELAY_NS + np.cumsum(increments, axis=1)
    return increments, cumulative

rng_demo = np.random.default_rng(FPGA_RANDOM_SEED)
inc_demo, cum_demo = generate_fpga_cumulative(FPGA_N_LEVELS, rng_demo)

fig, axes = plt.subplots(2, 1, figsize=(15, 9), sharex=True)
for lid in [1, 5, 9, 13]:
    axes[0].step(range(16), inc_demo[lid-1], where="mid", lw=1.8,
                 label=f"L{lid}")
    axes[1].step(range(16), cum_demo[lid-1], where="mid", lw=2.0,
                 label=f"L{lid}")
axes[0].set_ylabel("本 kick 新增随机量 ξ [ns]")
axes[0].set_title("图 V220-5a  每个 kick 的随机增量（不是绝对延时）")
axes[0].legend(ncol=4); axes[0].grid(alpha=0.25)
axes[1].set_xlabel("全局 kick 序号"); axes[1].set_ylabel("实际累计延时 d [ns]")
axes[1].set_title("图 V220-5b  实际延时 = 8ns 整体补偿 + 前面所有随机增量之和")
axes[1].legend(ncol=4); axes[1].grid(alpha=0.25)
plt.tight_layout()
plt.show()

lid = 5
print(f"7. L{lid} 示例：")
for k in range(6):
    print(f"  K{k}: increment={inc_demo[lid-1,k]:>2d}ns, "
          f"cumulative={cum_demo[lid-1,k]:>3d}ns")
"""


CELL_8 = r"""# ============================================================================
# 8. 一次雷达对射：真实信号 + 对射落点 + XM 聚类
# ============================================================================
KICK_PERIOD_NS = KICK_SPACING / NS      # 外来雷达同频，落点对它取模
V220_DEMO_D = XM_DEMO_D                 # 沿用 v13/v20 演示距离，便于对照

def cluster_positions(positions, sep_ns):
    # 按峰宽把相邻落点聚成峰；返回 [[position...], ...]。
    if len(positions) == 0:
        return []
    pos = sorted(float(x) for x in positions)
    clusters = [[pos[0]]]
    for x in pos[1:]:
        if x - clusters[-1][-1] < sep_ns:
            clusters[-1].append(x)
        else:
            clusters.append([x])
    return clusters

def radar_positions_for_laser(lid, cumulative, phase_ns=RADAR_PHASE_NS,
                              wrap=True):
    # 同频对射：脉冲每 KICK_PERIOD_NS 到一次，记录位置对周期取模；
    # 只有落进 [0, TOF_WINDOW] 才会被记录。
    kicks = SHOT_KICKS[lid]
    positions = np.array([phase_ns - cumulative[lid-1, k] for k in kicks],
                         dtype=float)
    if wrap:
        positions = np.mod(positions, KICK_PERIOD_NS)
    visible = positions[(positions >= 0) & (positions <= TOF_WINDOW / NS)]
    return kicks, positions, visible


def true_signal_position_ns(D):
    # 真回波 rec_tof = 2D/c，与编码、与 FPGA 延时都无关（发射与参考零点一起平移）
    return 2.0 * D / C_LIGHT / NS


demo_lid = 5
demo_kicks, demo_pos, demo_visible = radar_positions_for_laser(demo_lid, cum_demo)
demo_clusters = cluster_positions(demo_visible, RADAR_PULSE_SEP_NS)
true_ns = true_signal_position_ns(V220_DEMO_D)

fig, axes = plt.subplots(3, 1, figsize=(15, 10), sharex=True)
colors = ["#2874a6", "#1e8449", "#d35400", "#7d3c98"]

# (a) 逐 shot：真实信号 + 对射，纵向错开，看它们是否叠在一起
ax = axes[0]
for i, (kick, pos) in enumerate(zip(demo_kicks, demo_pos)):
    y = i * 1.2
    ax.hlines(y, 0, TOF_WINDOW / NS, color="0.9", lw=0.8, zorder=0)
    ax.vlines(true_ns, y, y + 1, color="black", lw=3.5, zorder=3)
    ax.vlines(pos, y, y + 1, color=colors[i], lw=3.5, zorder=3)
    ax.scatter([pos], [y + 1], color=colors[i], s=60, zorder=4)
    overlap = abs(pos - true_ns) < RADAR_PULSE_SEP_NS
    ax.text(TOF_WINDOW / NS * 1.005, y + 0.4,
            f"K{kick}: 对射 {pos:.0f}ns，离真峰 {pos-true_ns:+.0f}ns"
            + ("  ← 叠在真峰上" if overlap else ""),
            fontsize=8, va="center",
            color="#c0392b" if overlap else "0.3")
ax.set_ylim(-0.2, len(demo_kicks) * 1.2 + 0.2)
ax.set_yticks([i * 1.2 + 0.5 for i in range(len(demo_kicks))])
ax.set_yticklabels([f"shot{i}" for i in range(len(demo_kicks))])
ax.set_title(f"图 V220-6a  L{demo_lid} 逐 shot：黑色 = 真实信号（D={V220_DEMO_D:.0f}m，"
             f"每 shot 都在同一位置），彩色 = 外来对射")
ax.grid(alpha=0.2, axis="x")

# (b) 4 次 shot 累加后的直方图：真信号累到 4，对射被推散
ax = axes[1]
bin_ns = HIST_BIN_NS
n_bin = int(TOF_WINDOW / NS / bin_ns)
edges = np.arange(n_bin + 1) * bin_ns
add_true = np.zeros(n_bin)
add_radar = np.zeros(n_bin)
add_true[min(int(true_ns / bin_ns), n_bin - 1)] += len(demo_kicks)
for pos in demo_visible:
    add_radar[min(int(pos / bin_ns), n_bin - 1)] += 1
centers = edges[:-1] + bin_ns / 2
ax.bar(centers, add_true, width=bin_ns * 3, color="black", label="真实信号 hist_add")
ax.bar(centers, add_radar, width=bin_ns * 3, color="#c0392b", label="对射 hist_add")
ax.axhline(RADAR_PULSE_SEP_NS * 0 + 1.5, color="#1e8449", ls="--", lw=1.2,
           label="ratio=1.5")
ax.axhline(2.5, color="#c0392b", ls="--", lw=1.2, label="ratio=2.5")
ax.set_ylabel("hist_add 计数")
ax.set_title("图 V220-6b  4 次 shot 累加：真信号 add=4；对射被累计延时打散后 add 变小")
ax.legend(ncol=4, fontsize=8); ax.grid(alpha=0.2)

# (c) 聚类后每个对射峰的重复次数
ax = axes[2]
for cluster in demo_clusters:
    center = np.mean(cluster)
    ax.bar(center, len(cluster), width=max(4, RADAR_PULSE_SEP_NS * 0.7),
           color="#c0392b" if len(cluster) >= 2 else "#1e8449",
           edgecolor="black")
    ax.text(center, len(cluster) + 0.08, f"add={len(cluster)}",
            ha="center", fontsize=9)
ax.axvline(true_ns, color="black", ls=":", lw=1.5,
           label=f"真实信号 {true_ns:.0f}ns")
ax.axvspan(true_ns - RADAR_PULSE_SEP_NS, true_ns + RADAR_PULSE_SEP_NS,
           color="black", alpha=0.10, label=f"真峰 ±{RADAR_PULSE_SEP_NS}ns")
ax.axhline(1.5, color="#1e8449", ls="--", label="ratio=1.5")
ax.axhline(2.5, color="#c0392b", ls="--", label="ratio=2.5")
ax.set_xlabel("在本机直方图中的记录位置 [ns]")
ax.set_ylabel("该峰累积的 shot 数")
ax.set_title("图 V220-6c  对射聚类峰的重复次数：add < ratio 才会被 XM 丢掉")
ax.legend(ncol=4, fontsize=8); ax.grid(alpha=0.2)
plt.tight_layout()
plt.show()

n_overlap = int(np.sum(np.abs(demo_pos - true_ns) < RADAR_PULSE_SEP_NS))
print(f"8. L{demo_lid}：真实信号在 {true_ns:.1f} ns（D={V220_DEMO_D:.0f}m）")
print(f"  对射落点 = {[round(float(p),1) for p in demo_pos]} ns")
print(f"  其中与真峰重叠（|Δ| < {RADAR_PULSE_SEP_NS}ns）的有 {n_overlap} 次")
if n_overlap:
    print("  -> 对射叠到真峰上：这一部分靠 XM 是分不开的，只能靠改相位/一字滤波")
"""


CELL_9 = r"""# ============================================================================
# 9. 为什么 N 小的时候滤不掉对射：解析 + 蒙特卡洛
# ============================================================================
# 结论先写在前面：
#   累计延时单调递增 => 同一激光器 4 次 shot 的对射落点按 kick 顺序单调排列。
#   所以只要【相邻两次 shot】的落点差 >= SEP，就不会聚成一个峰。
#   相邻两次 shot 的落点差 = 这两个 kick 之间所有增量之和 Δ。
#   对相邻 kick（如 L5 的 K0->K1），Δ 就是【单个增量 ξ】。
#   ξ ∈ 8×{0..N-1}；要 ξ >= SEP=12ns 必须 ξ >= 16ns，即 N >= 3。
#   => N=1,2 时相邻 kick 的对射【必然】重合，add>=2，ratio=1.5/2.5 都滤不掉。
#   这与「累计 16 个 kick 的随机空间有 N^16 种」无关。
FPGA_STEP_CONTRAST_NS = 16   # 仅作对照：把步长提到一步就能跨过峰宽

def survives_xm_unit_pulses(positions, ratio, sep_ns):
    # 每个 shot 一个单位脉冲；任一聚类峰的 shot 数 >= ratio 即有残留。
    clusters = cluster_positions(positions, sep_ns)
    return any(len(c) >= ratio for c in clusters)


def sum_increment_pmf(m, n_levels, step_ns):
    # m 个独立增量之和的分布；返回 (取值[ns], 概率)
    pmf = np.array([1.0])
    single = np.ones(n_levels) / n_levels
    for _ in range(m):
        nxt = np.zeros(len(pmf) + n_levels - 1)
        for i, p in enumerate(pmf):
            if p:
                nxt[i:i + n_levels] += p * single
        pmf = nxt
    return np.arange(len(pmf)) * step_ns, pmf


def analytic_survive_ratio15(n_levels, step_ns, sep_ns):
    # ratio=1.5 => 只要任意相邻两 shot 落点差 < sep 就残留
    per_laser = []
    for lid in laser_ids:
        ks = SHOT_KICKS[lid]
        p_all_separated = 1.0
        for ka, kb in zip(ks[:-1], ks[1:]):
            vals, pmf = sum_increment_pmf(kb - ka, n_levels, step_ns)
            p_all_separated *= float(pmf[vals >= sep_ns].sum())
        per_laser.append(1.0 - p_all_separated)
    return float(np.mean(per_laser))


def monte_carlo_radar_cumulative(n_levels, ratio, n_trials, rng, step_ns=None):
    survive_module = survive_channels = hit_channels = survive_given_hit = 0
    max_delay_sum = 0.0
    for _ in range(n_trials):
        increments, cumulative = generate_fpga_cumulative(
            n_levels, rng, step_ns=step_ns)
        max_delay_sum += cumulative.max()
        trial_survive = False
        # 外来雷达与我们同频：相对相位在【整个 kick 周期】内均匀
        phase = rng.uniform(0, KICK_PERIOD_NS)
        for lid in laser_ids:
            _, _, visible = radar_positions_for_laser(lid, cumulative, phase)
            hit = len(visible) > 0
            hit_channels += int(hit)
            channel_survives = survives_xm_unit_pulses(
                visible, ratio, RADAR_PULSE_SEP_NS)
            survive_channels += int(channel_survives)
            if hit:
                survive_given_hit += int(channel_survives)
            if channel_survives:
                trial_survive = True
        survive_module += int(trial_survive)
    n_ch = n_trials * N_LASERS
    return {
        "p_hit_channel": hit_channels / n_ch,
        "p_survive_channel": survive_channels / n_ch,
        "p_survive_given_hit": survive_given_hit / max(hit_channels, 1),
        "p_survive_module": survive_module / n_trials,
        "mean_max_delay": max_delay_sum / n_trials,
    }


N_LEVEL_LIST = [1, 2, 3, 4, 5, 6, 8]
MC_ROWS, MC_CONTRAST = [], []
rng_mc = np.random.default_rng(FPGA_RANDOM_SEED + 1)

print("9. 累计 FPGA 对射蒙卡")
print(f"   峰宽 SEP={RADAR_PULSE_SEP_NS}ns，步长={FPGA_STEP_NS}ns "
      f"=> 单个增量要跨过峰宽必须 ≥ {int(np.ceil(RADAR_PULSE_SEP_NS/FPGA_STEP_NS))*FPGA_STEP_NS}ns，"
      f"即 N ≥ {int(np.ceil(RADAR_PULSE_SEP_NS/FPGA_STEP_NS))+1}")
print(f"  {'N':>3} {'P(单增量≥SEP)':>13} {'MC残留|命中 r1.5':>16} "
      f"{'解析 r1.5':>10} {'MC残留|命中 r2.5':>16} {'命中率':>8} {'平均最大累计':>12}")
print("  " + "-" * 92)
for n_level in N_LEVEL_LIST:
    res = {}
    for ratio in (1.5, 2.5):
        res[ratio] = monte_carlo_radar_cumulative(
            n_level, ratio, FPGA_MC_TRIALS, rng_mc)
        MC_ROWS.append((n_level, ratio, res[ratio]))
    p_step = sum(1 for i in range(n_level)
                 if i * FPGA_STEP_NS >= RADAR_PULSE_SEP_NS) / n_level
    ana = analytic_survive_ratio15(n_level, FPGA_STEP_NS, RADAR_PULSE_SEP_NS)
    print(f"  {n_level:>3d} {p_step:>13.2f} {res[1.5]['p_survive_given_hit']:>16.2%} "
          f"{ana:>10.2%} {res[2.5]['p_survive_given_hit']:>16.2%} "
          f"{res[1.5]['p_hit_channel']:>8.1%} "
          f"{res[1.5]['mean_max_delay']:>10.0f}ns")

# 对照：把步长从 8ns 提到 16ns（一步即可跨过 12ns 峰宽），其余不变
rng_c = np.random.default_rng(FPGA_RANDOM_SEED + 2)
for n_level in N_LEVEL_LIST:
    MC_CONTRAST.append((n_level, monte_carlo_radar_cumulative(
        n_level, 1.5, max(FPGA_MC_TRIALS // 4, 1000), rng_c,
        step_ns=FPGA_STEP_CONTRAST_NS)))

fig, axes = plt.subplots(1, 3, figsize=(20, 5.8))

# (a) 残留概率
ax = axes[0]
for ratio, color in [(1.5, "#1e8449"), (2.5, "#c0392b")]:
    rows = [(n, r) for n, rr, r in MC_ROWS if rr == ratio]
    ax.plot([n for n, _ in rows], [r["p_survive_given_hit"] * 100 for _, r in rows],
            "o-", lw=2.2, color=color, label=f"MC 单通道 ratio={ratio}")
ax.plot(N_LEVEL_LIST,
        [analytic_survive_ratio15(n, FPGA_STEP_NS, RADAR_PULSE_SEP_NS) * 100
         for n in N_LEVEL_LIST],
        "k--", lw=1.6, label="解析 ratio=1.5（不计超窗）")
ax.plot([n for n, _ in MC_CONTRAST],
        [r["p_survive_given_hit"] * 100 for _, r in MC_CONTRAST],
        "s-.", lw=1.8, color="#2874a6",
        label=f"对照：步长={FPGA_STEP_CONTRAST_NS}ns, ratio=1.5")
ax.set_xlabel(f"每 kick 随机增量的档数 N（增量 = 步长×[0, N-1]）")
ax.set_ylabel("对射经 XM 后仍残留的概率 [%]")
ax.set_title("图 V220-7a  对射残留率（口径：该通道确实被对射打到）\n"
             "8ns 步长在 N≤2 时恒 100%：单个增量跨不过 12ns 峰宽")
ax.set_ylim(-3, 105)
ax.legend(fontsize=8); ax.grid(alpha=0.25)

# (b) 相邻 shot 落点差的分布（诊断图）
ax = axes[1]
demo_pairs = [(0, 1), (1, 2), (2, 6)]     # L5 的三段相邻 shot 间隔
for n_level, color in [(2, "#c0392b"), (3, "#e67e22"), (8, "#1e8449")]:
    vals_all, pmf_all = [], []
    for ka, kb in demo_pairs:
        vals, pmf = sum_increment_pmf(kb - ka, n_level, FPGA_STEP_NS)
        vals_all.append(vals); pmf_all.append(pmf)
    m = max(len(p) for p in pmf_all)
    agg = np.zeros(m)
    for p in pmf_all:
        agg[:len(p)] += p / len(pmf_all)
    ax.plot(np.arange(m) * FPGA_STEP_NS, agg * 100, "o-", color=color,
            lw=1.8, ms=4, label=f"N={n_level}")
ax.axvline(RADAR_PULSE_SEP_NS, color="black", ls="--", lw=2,
           label=f"峰宽 SEP={RADAR_PULSE_SEP_NS}ns")
ax.axvspan(0, RADAR_PULSE_SEP_NS, color="red", alpha=0.10)
ax.set_xlabel("相邻两次 shot 的对射落点差 Δ [ns]")
ax.set_ylabel("概率 [%]")
ax.set_title("图 V220-7b  为什么小 N 无效：Δ 落进红区就聚成一个峰\n"
             "（以 L5 的 K0→K1、K1→K2、K2→K6 三段平均）")
ax.legend(fontsize=8); ax.grid(alpha=0.25)

# (c) 延时预算
ax = axes[2]
rows15 = [(n, r) for n, rr, r in MC_ROWS if rr == 1.5]
ax.plot([n for n, _ in rows15], [r["mean_max_delay"] for _, r in rows15],
        "s-", lw=2.2, color="#7d3c98", label=f"步长 {FPGA_STEP_NS}ns")
ax.plot([n for n, _ in MC_CONTRAST], [r["mean_max_delay"] for _, r in MC_CONTRAST],
        "s-.", lw=1.8, color="#2874a6", label=f"步长 {FPGA_STEP_CONTRAST_NS}ns")
ax.axhline((KICK_SPACING - TOF_WINDOW) / NS, color="red", ls="--",
           label="单 kick 空闲预算 200ns")
ax.set_xlabel("随机增量档数 N")
ax.set_ylabel("16 kick 内平均最大累计延时 [ns]")
ax.set_title("图 V220-7c  代价：累计延时单调增长，必须守住 200ns 空闲预算")
ax.legend(fontsize=8); ax.grid(alpha=0.25)
plt.tight_layout()
plt.show()

print()
print("解读：")
print(f"  1) 8ns 步长 + {RADAR_PULSE_SEP_NS}ns 峰宽：单个增量只有取到 16ns 以上才算「分开」，")
print(f"     所以 N=1,2 残留 100%，与随机空间有 2^16 种无关 —— 卡的是【一步的大小】。")
print("  2) 大多数激光器有 3 个连续 kick（如 L5 的 K0,K1,K2），这两段间隔各只有 1 个增量，")
print("     是整条链上最难分开的地方，决定了整体残留率。")
print(f"  3) 把步长提到 {FPGA_STEP_CONTRAST_NS}ns，N=2 就能大幅压低残留，但累计延时也翻倍，")
print("     要同时看图 V220-7c 的 200ns 预算线。")
"""


SUMMARY = f"""# ============================================================================
# 10. v220 总结
# ============================================================================
print("=" * 82)
print("v220 总结")
print("=" * 82)
print(f'''
1. 模组内串扰：
   - ratio=1.5 使用 {B15}ns tcode
   - ratio=2.5 使用 {B25}ns tcode
   - 两套都满足：同类鬼影间隔 ≥ {SEP_NS}ns，且鬼影离真峰 ≥ {SEP_NS}ns

2. FPGA 模型已纠正：
   - random_increment 是本 kick 新增随机量
   - actual_delay = 8ns 整体补偿 + 从 K0 到当前 kick 的全部增量之和
   - 例：16ns 后再加 8ns，实际累计为 24ns

3. 雷达对射：
   - tcode 不控制外来雷达
   - 累计 FPGA 延时让对射峰在多个 shot 间移动，再由 XM 判断
   - ratio=1.5 对双碰敏感；ratio=2.5 允许双碰但不允许三碰
   - 对射有可能直接叠在真峰上（见图 V220-6），那一部分 XM 分不开

4. 决定滤除效果的不是随机空间大小，而是【一步的大小】：
   - 累计延时单调递增 -> 4 次 shot 落点单调排列 -> 只需相邻两次分开
   - 相邻 kick 之间只有一个增量 ξ，要跨过 {SEP_NS}ns 峰宽必须 ξ ≥ 16ns
   - 8ns 步长下 N=1,2 的 ξ 只能取 0 或 8ns，恒 < {SEP_NS}ns -> 残留 100%
   - 解析式与蒙卡在图 V220-7a 中重合，说明这是模型必然，不是统计噪声

5. 工程约束：
   - 增量档数 N 增大，随机性更强
   - 但累计延时随 kick 单调增长，必须同时检查 200ns 空闲预算
   - 提高步长（如 16ns）能用很小的 N 解决问题，代价是累计延时翻倍

6. 下一步：
   - tcode + 累计 FPGA 的联合模组内/外来雷达仿真（见 v23）
   - 连续三个角度的一字滤波（见 v24）
''')
"""


if not os.path.exists(SRC_NB):
    raise SystemExit(f"找不到 {SRC_NB}")

with open(SRC_NB, encoding="utf-8") as f:
    nb21 = json.load(f)

base_cells = []
for cell in nb21["cells"]:
    copied = dict(cell)
    if copied["cell_type"] == "code":
        copied["outputs"] = []
        copied["execution_count"] = None
    base_cells.append(copied)

new_cells = [
    md_cell("v220_section", SECTION),
    code_cell("v220_01_vars", CELL_1),
    code_cell("v220_02_timing", CELL_2),
    code_cell("v220_03_tcode", CELL_3),
    code_cell("v220_04_constraints", CELL_4),
    code_cell("v220_05_effect", CELL_5),
    md_cell("v220_06_fpga_doc", FPGA_DOC),
    code_cell("v220_07_fpga_walk", CELL_7),
    code_cell("v220_08_radar_example", CELL_8),
    code_cell("v220_09_radar_mc", CELL_9),
    code_cell("v220_10_summary", SUMMARY),
]

notebook = {
    "cells": [md_cell("v220_overview", OVERVIEW)] + base_cells + new_cells,
    "metadata": nb21.get("metadata", {}),
    "nbformat": 4,
    "nbformat_minor": 5,
}

with open(OUT_NB, "w", encoding="utf-8") as f:
    json.dump(notebook, f, ensure_ascii=False, indent=1)

print(f"已生成 {OUT_NB}")
print(f"  继承 v21：{len(base_cells)} cell")
print(f"  v220 新增：{1 + len(new_cells)} cell")
print(f"  新增绘图 cell：7 个")
