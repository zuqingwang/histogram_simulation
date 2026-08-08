# -*- coding: utf-8 -*-
"""
build_crosstalk_v21.py —— 生成 crosstalk_sim_v21.ipynb
========================================================
【构建方式】和 v20 一样的规矩：
    v21 = crosstalk_sim_v20.ipynb 的全部 cell **原样继承** + 后面追加 v21 新增 cell。
    （v20 本身又是 v13 原样继承 + XM 部分，所以 v13 的 cell5/cell9 一路都在。）
    不删任何东西，不改任何一行既有代码。

【v21 新增】tcode 实装 —— 把 tx_trig_dly 换成求解出来的编码表
    码表来源：docs/tcode/solve_tcode.py 的数值搜索结果（docs/tcode/tcode_table.py）
      - 峰宽裕度 12 ns（单光子脉宽 8ns + 裕度）
      - 码预算   64 ns（用户给定），实际用到 63 ns
      - 实测     1~600m 扫描，鬼影残留 0.000%，真目标误杀 0

    新增内容：
      1) 码表 + 码预算核查
      2) 编码自检（不跑仿真，直接验证约束）：同 kick 串扰 / 跨 kick 混叠 / 自身混叠
      3) 图E 码矩阵对比（Excel vs 新 tcode）
      4) 图F「最大重复次数」热图 —— 全绿即通过
      5) 图G 16 宫格 XM 滤前/滤后（新 tcode）
      6) 图H 距离扫描验收（Excel vs 新 tcode）
      7) 总结 + 下一步（FPGA 逐 kick 全局抖动）

本脚本只生成 notebook。

缩写：
  XM（XtalkMark，串扰标记）、TOF（Time of Flight，飞行时间）、
  IRF（Instrument Response Function，仪器响应函数）、
  SPAD（Single-Photon Avalanche Diode，单光子雪崩二极管）、
  FPGA（Field-Programmable Gate Array，现场可编程门阵列）。
"""
import json
import os
import sys

SRC_NB = "crosstalk_sim_v20.ipynb"
OUT_NB = "crosstalk_sim_v21.ipynb"
TBL_PY = os.path.join("docs", "tcode", "tcode_table.py")


# ============================================================================
# 读入求解器产出的码表（不手抄，避免抄错）
# ============================================================================
def load_table():
    ns = {}
    with open(TBL_PY, encoding="utf-8") as f:
        exec(compile(f.read(), TBL_PY, "exec"), ns)
    return ns["TCODE_TABLE"], ns["TCODE_SEP_NS"], ns["TCODE_BUDGET_NS"]


TCODE_TABLE, SEP_NS, BUDGET_NS = load_table()
LASERS = sorted({l for (l, k) in TCODE_TABLE})
KICKS_OF = {l: sorted(k for (a, k) in TCODE_TABLE if a == l) for l in LASERS}
MAXTX = max(TCODE_TABLE.values())


def table_literal():
    """把码表排版成 notebook 里好读的字面量。"""
    lines = []
    for l in LASERS:
        items = ", ".join(f"({l},{k}): {TCODE_TABLE[(l,k)]:>2d}" for k in KICKS_OF[l])
        lines.append(f"    {items},")
    return "\n".join(lines)


# ============================================================================
# 【插在最前面】v21 总览
# ============================================================================
CELL_V21_OVERVIEW = f'''# 串扰仿真 v21 —— tcode 编码实装（64 ns 预算内做到鬼影残留 0）

> **本 notebook = v20 原封不动 + 末尾追加 tcode 部分。**
> 而 v20 又 = v13 原封不动 + XM 部分。所以 v13 的 cell5（kick 栅格）、cell9（逐回波堆叠柱）
> 和 v20 的整个 XM 链路都还在，一个都没删。v21 的新增内容从 **「v21 新增部分」** 往后开始。

## v20 留下的问题

v20 证明了 XM（XtalkMark，串扰标记）本身实现是对的，但用 Excel 现有的 `tx_trig_dly`
**一个鬼影都滤不掉** —— 因为每个激光器 4 个 kick 的码是同一个值，
同一条串扰 4 次落到同一个 bin，`hist_add` 也累到 4，和真目标数值上完全一样。

**瓶颈在编码，不在算法。** v21 就是把编码补上。

## v21 给出的编码

| 项 | 值 |
|---|---|
| 编码对象 | `tx_trig_dly`（1 ns 步长） |
| 单光子脉宽 | 8 ns |
| **落点最小间隔（峰宽裕度）** | **{SEP_NS} ns** |
| **码预算** | **{BUDGET_NS} ns**（实际用到 {MAXTX} ns） |
| 码预算上限 | `KICK_SPACING − TOF_WINDOW` = 2200 − 2000 = 200 ns |
| 求解方式 | 数值搜索（`docs/tcode/solve_tcode.py`），非解析构造 |
| **实测（1~600m）** | **鬼影残留 0.000 %，真目标误杀 0** |

## 为什么是「搜索」而不是「公式」

`docs/tcode/tcode_scheme.md` 里推导过一个解析构造（线性同余 + 二次偏移），
它的问题是**只保证同 kick 串扰那一类**，跨 kick 混叠和自身混叠靠运气，
而且要 184 ns 码范围 —— 几乎吃满 200 ns 的 kick 间隙。

v21 改成把**三类鬼影源写成统一约束**再数值求解：

> 对每一类鬼影源（发射器 a、kick 偏移 Δk），接收器 b 在它 4 次 shot 里收到的该类鬼影，
> 落点由 `c[a][k−Δk] − c[b][k]` 决定；要求这 4 个值**两两相差 ≥ {SEP_NS} ns**。

结果只用 **{MAXTX} ns** 就把残留压到 0，比解析构造省 2/3 的时间、成绩还更好。

| 方案 | 码范围 | 残留率(ratio=1.6) | 误杀 |
|---|---|---|---|
| Excel 现状 | 50 ns | 100 % | 0 |
| 解析构造（tcode_scheme.md） | 184 ns | 5.59 % | 0 |
| **v21 搜索解** | **{MAXTX} ns** | **0.000 %** | **0** |

## ⚠ 两条使用前提

1. **码表是针对这张具体时序表搜出来的**（哪个激光器在哪 4 个 kick 发光）。
   时序表一改就必须重跑 `docs/tcode/solve_tcode.py` 重搜。
2. **仿真仍是理想 δ 回波模型**。峰宽是靠「落点间隔 ≥ {SEP_NS} ns」这条**约束**保证的，
   不是在仿真里用真实宽度验证的。要坐实，需要把回波展宽（`PULSE_W`）加进仿真再跑一遍。

## 新增 cell 一览

| cell | 内容 |
|---|---|
| 「v21 新增部分」md | tcode 原理速览 + 本版码表的来历 |
| TC-1 | 码表 + 码预算核查 + 构建发光事件表 |
| TC-2 | **编码自检**（不跑仿真，直接验证三类约束是否都满足） |
| TC-3 | 图E 码矩阵对比（Excel vs v21） |
| TC-4 | 图F「最大重复次数」热图 —— 全绿即通过 |
| TC-5 | 图G 16 宫格 XM 滤前/滤后 |
| TC-6 | 图H 距离扫描验收 |
| TC-7 | 总结 + 下一步（FPGA 逐 kick 全局抖动打外来雷达） |
'''


# ============================================================================
# 【追加】v21 说明
# ============================================================================
CELL_TC_DOC = f'''---

# v21 新增部分

---

# tcode：把串扰"搬"到 XM 能看见的地方

## 1. 编码只搬鬼影，不搬真目标

一次发光既是发射、也开一个 2000 ns 接收窗。接收方拿**自己**的发光时刻当测距零点：

$$\\text{{rec\\_tof}}
= \\underbrace{{\\frac{{2D}}{{c}}}}_{{\\text{{物体距离}}}}
+ \\underbrace{{\\big(tx_{{\\text{{发}}}} - tx_{{\\text{{收}}}}\\big)}}_{{\\textbf{{码差}}}}
- \\underbrace{{\\big(k_{{\\text{{收}}}} - k_{{\\text{{发}}}}\\big)\\cdot T_{{\\rm kick}}}}_{{\\text{{跨 kick 混叠}}}}$$

真回波（发 = 收、同 kick）后两项都是 0，**与编码无关**。
所以编码怎么改都不会动真目标 —— 这就是「只杀鬼影不伤真身」的机理，也是本版误杀恒为 0 的原因。

## 2. 编码要满足什么

XM 判据是 `hist_max × xm_ratio > hist_add`，等价于 `add/max < xm_ratio`，
而 `add/max` = **这个峰在几次 shot 里出现过**。真目标恒为 4（`N_ACC`），
所以只要让每条串扰在 4 次 shot 里落到 **4 个不同的 bin**，它的 `add/max` 就是 1，必被丢弃。

> **约束**：对每一类鬼影源（发射器 `a`、kick 偏移 `Δk`），
> 接收器 `b` 的 4 次 shot 收到的该类鬼影，落点由 `c[a][k−Δk] − c[b][k]` 决定，
> 要求这 4 个值**两两相差 ≥ {SEP_NS} ns**（≥ 峰宽，否则会糊成一个峰）。

三类鬼影源都要管：

| Δk | a vs b | 叫什么 | 备注 |
|---|---|---|---|
| 0 | a ≠ b | 同 kick 串扰 | 最常见，解析构造能治 |
| 1 | a ≠ b | 跨 kick 混叠 | D > 300m 时别人上一个 kick 的光折回来 |
| 1 | a = b | **自身混叠** | 我自己上一个 kick 的光，解析构造治不了（见 tcode_scheme.md 第 10 节） |

（D ≤ 600 m 时只有 Δk ∈ {{0, 1}} 会落进窗内，所以约束就这两档。）

## 3. 两条硬预算

| 约束 | 数值 | 后果 |
|---|---|---|
| 落点最小间隔 ≥ 峰宽 | {SEP_NS} ns | 小了就糊成一个峰，编了等于没编 |
| `max(tx) ≤ KICK_SPACING − TOF_WINDOW` | 200 ns | 超了窗会伸进下一个 kick |

本版实际用到 **{MAXTX} ns**，占 kick 间隙的 {MAXTX/200:.0%}，
剩下的留给后面 FPGA 的 8 ns 步长随机抖动（见最后一节）。
'''


# ============================================================================
# TC-1：码表
# ============================================================================
CELL_TC_TABLE = f'''# ============================================================================
# TC-1  v21 码表：tx_trig_dly [ns]
# ============================================================================
#   来源：docs/tcode/solve_tcode.py 的数值搜索（seed 固定，可复现）
#   参数：峰宽裕度 {SEP_NS} ns（单光子脉宽 8ns + 裕度）、码预算 {BUDGET_NS} ns
#   注意：本表绑定当前时序表（哪个激光器在哪 4 个 kick 发光）。时序一改必须重搜。
# ============================================================================
TCODE_SEP_NS    = {SEP_NS}      # 两个落点至少要差这么多 ns 才算"分得开"
TCODE_BUDGET_NS = {BUDGET_NS}      # 用户给定的编码时间预算

# (laser_id, kick) -> tx_trig_dly [ns]
TCODE_V21 = {{
{table_literal()}
}}

# ---- 开关：True = 用 v21 新码；False = 退回 Excel 原码做对照 ----
USE_V21_TCODE = True


def tcode_v21(lid, k, tx0):
    """v21 编码函数，签名与 v20 的 make_tcode_fn() 产物一致。"""
    return TCODE_V21.get((lid, k), tx0)


# ---- 码预算核查 ----
_hard_limit = (KICK_SPACING - TOF_WINDOW) / NS      # = 200 ns
_maxtx = max(TCODE_V21.values())
print("=" * 78)
print("v21 码预算核查")
print("=" * 78)
print(f"  实际最大码值      max(tx) = {{_maxtx:>4.0f}} ns")
print(f"  用户给定预算              = {{TCODE_BUDGET_NS:>4d}} ns   "
      f"-> {{'OK' if _maxtx <= TCODE_BUDGET_NS else '!! 超出'}}")
print(f"  硬上限 KICK−TOF           = {{_hard_limit:>4.0f}} ns   "
      f"-> {{'OK' if _maxtx <= _hard_limit else '!! 窗会伸进下一个 kick'}}")
print(f"  占 kick 间隙比例          = {{_maxtx/_hard_limit:>4.0%}}   "
      f"（剩 {{_hard_limit-_maxtx:.0f}} ns 留给 FPGA 抖动）")

# ---- 用 v20 已有的 build_firings_tcode() 构建发光事件表（不改 v20 任何代码）----
firings_v21   = build_firings_tcode(tcode_v21)
firings_excel = build_firings_tcode(make_tcode_fn("excel"))
FR_MAIN = firings_v21 if USE_V21_TCODE else firings_excel

print(f"\\n  码表：{{len(TCODE_V21)}} 个 (激光器, kick) 组合")
print(f"  {{'激光器':>6}}  各 kick 的 tx_trig_dly [ns]")
print("  " + "-" * 60)
for _l in laser_ids:
    _cells = "  ".join(f"K{{k}}={{TCODE_V21[(_l,k)]:>2d}}" for k in SHOT_KICKS[_l])
    print(f"  {{('L'+str(_l)):>6}}  {{_cells}}")
'''


# ============================================================================
# TC-2：编码自检
# ============================================================================
CELL_TC_CHECK = r'''# ============================================================================
# TC-2  编码自检 —— 不跑仿真，直接验证约束
# ============================================================================
#   对每一类鬼影源，把接收器 4 次 shot 的落点算出来，看两两间隔够不够。
#   这是【充分条件】：全过 ⟹ 每条串扰的 add/max 必为 1 ⟹ XM 必能滤掉。
# ============================================================================
def ghost_classes(gap=CROSSTALK_MAX_GAP, dks=(0, 1)):
    """列出所有鬼影源类别。
       返回 [(发射器 a, 接收器 b, Δk, 类别名, [(k_发, k_收), ...]), ...]"""
    out = []
    for b in laser_ids:
        for a in laser_ids:
            if abs(a - b) > gap:
                continue                      # 编号间隔 > 阈值，空间上可忽略
            for dk in dks:
                if a == b and dk == 0:
                    continue                  # 这是真回波，不是鬼影
                prs = [(kb - dk, kb) for kb in SHOT_KICKS[b]
                       if (kb - dk) in SHOT_KICKS[a]]
                if len(prs) >= 2:
                    kind = ("自身混叠" if a == b else
                            ("同kick串扰" if dk == 0 else "跨kick串扰"))
                    out.append((a, b, dk, kind, prs))
    return out


def check_tcode(code_fn, sep=None, verbose_fail=6):
    """验证一套编码是否满足"落点两两相差 ≥ sep"。返回 (是否全过, 明细)。"""
    sep = TCODE_SEP_NS if sep is None else sep
    rows, fails = [], []
    for (a, b, dk, kind, prs) in ghost_classes():
        d = sorted(code_fn(a, ka, 0) - code_fn(b, kb, 0) for (ka, kb) in prs)
        gaps = [d[i + 1] - d[i] for i in range(len(d) - 1)]
        mind = min(gaps) if gaps else 10 ** 9
        rows.append({"a": a, "b": b, "dk": dk, "kind": kind, "n": len(d),
                     "diffs": d, "min_gap": mind, "ok": mind >= sep})
        if mind < sep:
            fails.append(rows[-1])
    return (len(fails) == 0), rows, fails


def report_check(name, code_fn):
    ok, rows, fails = check_tcode(code_fn)
    from collections import Counter
    per_kind = Counter()
    per_kind_ok = Counter()
    for r in rows:
        per_kind[r["kind"]] += 1
        per_kind_ok[r["kind"]] += int(r["ok"])
    print(f"\n【{name}】共 {len(rows)} 个鬼影源类别，要求落点两两相差 ≥ {TCODE_SEP_NS} ns")
    print(f"  {'类别':>10} {'通过/总数':>10} {'最小间隔':>9}")
    print("  " + "-" * 34)
    for k in ("同kick串扰", "跨kick串扰", "自身混叠"):
        if per_kind[k]:
            mg = min(r["min_gap"] for r in rows if r["kind"] == k)
            print(f"  {k:>10} {per_kind_ok[k]:>5d}/{per_kind[k]:<4d} {mg:>8d} ns")
    allmin = min(r["min_gap"] for r in rows)
    print(f"  {'合计':>10} {sum(per_kind_ok.values()):>5d}/{len(rows):<4d} {allmin:>8d} ns")
    print(f"  -> {'全部通过 √' if ok else f'!! {len(fails)} 个类别不达标'}")
    for r in fails[:6]:
        print(f"     L{r['a']}->L{r['b']} Δk={r['dk']} {r['kind']}: "
              f"落点 {r['diffs']} ns，最小间隔 {r['min_gap']} ns")
    return ok


print("=" * 78)
print("TC-2  编码自检")
print("=" * 78)
_ok_excel = report_check("Excel 现状（对照）", make_tcode_fn("excel"))
_ok_v21   = report_check("v21 搜索解", tcode_v21)

print("\n" + "=" * 78)
print("解读：Excel 那一栏「最小间隔 = 0」意味着 4 次 shot 的鬼影落在【同一个 bin】，")
print("      hist_add 累到 4，与真目标无法区分 —— 这正是 v20 里 XM 滤除率为 0 的原因。")
print(f"      v21 把最小间隔抬到 ≥ {TCODE_SEP_NS} ns，每条串扰的 add/max 都变成 1。")
print("=" * 78)
'''


# ============================================================================
# TC-3：图 E 码矩阵
# ============================================================================
CELL_TC_FIG_CODE = r'''# ============================================================================
# TC-3  图E：码矩阵 c[l][k] 对比
# ============================================================================
fig, axes = plt.subplots(1, 2, figsize=(17, 6.5))
for ax, (name, fn) in zip(axes, [("Excel 现状（每行一个常数）", make_tcode_fn("excel")),
                                 (f"v21 搜索解（≤{TCODE_BUDGET_NS}ns）", tcode_v21)]):
    M = np.full((N_LASERS, 16), np.nan)
    for i, l in enumerate(laser_ids):
        for k in SHOT_KICKS[l]:
            M[i, k] = fn(l, k, 0)
    im = ax.imshow(M, cmap="viridis", aspect="auto", origin="lower",
                   extent=[-0.5, 15.5, laser_ids[0] - 0.5, laser_ids[-1] + 0.5])
    for i, l in enumerate(laser_ids):
        for k in SHOT_KICKS[l]:
            ax.text(k, l, f"{int(M[i,k])}", ha="center", va="center", fontsize=7,
                    color=("white" if M[i, k] < np.nanmax(M) * 0.55 else "black"))
    ax.set_xticks(range(16)); ax.set_yticks(laser_ids)
    ax.tick_params(labelsize=7)
    ax.set_xlabel("kick 序号 k"); ax.set_ylabel("激光器编号 l")
    n_val = len({int(v) for v in M[~np.isnan(M)]})
    ax.set_title(f"{name}\n格内 = tx_trig_dly [ns]，共用到 {n_val} 种码值，"
                 f"max = {np.nanmax(M):.0f} ns", fontsize=11)
    fig.colorbar(im, ax=ax, fraction=0.04, label="tx_trig_dly [ns]")
plt.suptitle("图E  码矩阵：左边每一行沿 kick 是常数（XM 失效的根源），右边每一行都在变",
             fontsize=13.5, y=1.02)
plt.tight_layout()
plt.show()
'''


# ============================================================================
# TC-4：图 F 重复次数热图
# ============================================================================
CELL_TC_FIG_REPEAT = r'''# ============================================================================
# TC-4  图F：「最大重复次数」热图 —— 验收图，全绿即通过
# ============================================================================
#   k = 对某一对 (发, 收)，同一个落点最多重复几次 shot
#     k = 1        -> 完美散开，XM 全能滤（绿）
#     k = N_ACC(4) -> 完全重合，XM 全滤不掉（红）
#   判据：k < XM_RATIO 才滤得掉
# ============================================================================
def repeat_matrix(code_fn, dk=0):
    """返回 (发射器 × 接收器) 的最大重复次数矩阵；落点间隔 < 峰宽视为重合。
       只画编号间隔 ≤ CROSSTALK_MAX_GAP 的对 —— 其余对在本模型里不构成串扰路径
       （detect_echoes 就是按这个阈值截断的），求解器也没有约束它们，
       画出来会让人误以为"没通过"。"""
    M = np.full((N_LASERS, N_LASERS), np.nan)
    for i, a in enumerate(laser_ids):
        for j, b in enumerate(laser_ids):
            if a == b and dk == 0:
                continue
            if abs(a - b) > CROSSTALK_MAX_GAP:
                continue
            prs = [(kb - dk, kb) for kb in SHOT_KICKS[b] if (kb - dk) in SHOT_KICKS[a]]
            if len(prs) < 2:
                continue
            d = sorted(code_fn(a, ka, 0) - code_fn(b, kb, 0) for (ka, kb) in prs)
            # 把间隔 < 峰宽的落点归并成一簇，簇内条数就是"重复次数"
            best, cur = 1, 1
            for t in range(1, len(d)):
                cur = cur + 1 if (d[t] - d[t - 1]) < TCODE_SEP_NS else 1
                best = max(best, cur)
            M[i, j] = best
    return M


fig, axes = plt.subplots(1, 3, figsize=(19, 6.2),
                         gridspec_kw={"width_ratios": [1, 1, 0.85]})
Ms = {}
for ax, (name, fn) in zip(axes[:2], [("Excel 现状", make_tcode_fn("excel")),
                                     ("v21 搜索解", tcode_v21)]):
    M = repeat_matrix(fn); Ms[name] = M
    im = ax.imshow(M, cmap="RdYlGn_r", vmin=1, vmax=max(N_ACC, 2), aspect="auto",
                   origin="lower",
                   extent=[laser_ids[0] - 0.5, laser_ids[-1] + 0.5,
                           laser_ids[0] - 0.5, laser_ids[-1] + 0.5])
    for i, a in enumerate(laser_ids):
        for j, b in enumerate(laser_ids):
            if not np.isnan(M[i, j]):
                ax.text(b, a, f"{int(M[i,j])}", ha="center", va="center", fontsize=6.5)
    ax.set_xticks(laser_ids); ax.set_yticks(laser_ids); ax.tick_params(labelsize=7)
    ax.set_xlabel("接收激光器 b"); ax.set_ylabel("发射激光器 a")
    band = [v for v in M.ravel() if not np.isnan(v)]
    ax.set_title(f"{name}（同 kick 串扰 Δk=0）\n只画编号间隔≤{CROSSTALK_MAX_GAP} 的对"
                 f"（共 {len(band)} 对），最坏重复次数 = {int(max(band)) if band else 0}",
                 fontsize=10.5)
    fig.colorbar(im, ax=ax, fraction=0.045, label="同一落点最多重复几次 shot")

ax = axes[2]
xs = np.arange(1, N_ACC + 1); w = 0.38
for t, (name, col) in enumerate([("Excel 现状", "#c0392b"), ("v21 搜索解", "#1e8449")]):
    band = [v for v in Ms[name].ravel() if not np.isnan(v)]
    cnt = [sum(1 for v in band if int(v) == k) for k in xs]
    ax.bar(xs + (t - 0.5) * w, cnt, w, label=name, color=col, edgecolor="k")
    for k, c in zip(xs, cnt):
        if c:
            ax.text(k + (t - 0.5) * w, c, str(c), ha="center", va="bottom", fontsize=8)
ax.axvspan(0.5, XM_RATIO, color="green", alpha=0.12)
ax.axvline(XM_RATIO, color="purple", lw=2, ls="--")
ax.text(XM_RATIO, ax.get_ylim()[1] * 0.95, f" XM_RATIO={XM_RATIO}", color="purple",
        fontsize=10, va="top")
ax.set_xticks(xs)
ax.set_xlabel("该串扰对的最大重复次数 k")
ax.set_ylabel(f"激光器对数（只统计编号间隔≤{CROSSTALK_MAX_GAP} 的）")
ax.set_title("绿色区（k < XM_RATIO）= XM 能滤掉\n验收标准：全部落在 k=1", fontsize=10.5)
ax.legend(fontsize=9); ax.grid(alpha=0.25, axis="y")

plt.suptitle("图F  验收图：编码的唯一任务就是把所有格子从红色(k=4)压到绿色(k=1)",
             fontsize=13.5, y=1.03)
plt.tight_layout()
plt.show()
'''


# ============================================================================
# TC-5：图 G 滤前滤后
# ============================================================================
CELL_TC_FIG_GRID = r'''# ============================================================================
# TC-5  图G：16 个激光器 XM 滤前/滤后（v21 编码）
# ============================================================================
def run_with(fr, D, ratio=None):
    """用给定发光事件表跑一遍完整流程（复用 v20 的函数，不改任何 v20 代码）。"""
    with use_firings(fr):
        recs = detect_echoes(D)
    hs = build_hists(recs)
    rs = crosstalk_mark_all(hs, ratio)
    return hs, rs, evaluate(hs, rs)


D_V21 = XM_DEMO_D                      # 沿用 v20/v13 的演示距离，便于逐版对照
hs21, rs21, st21 = run_with(FR_MAIN, D_V21)

b_lo, b_hi = occupied_range(hs21, pad=15)
x = np.arange(b_lo, b_hi + 1) * HIST_BIN_M
nrow = int(np.ceil(N_LASERS / 4))
fig, axes = plt.subplots(nrow, 4, figsize=(22, 3.1 * nrow), sharex=True)
axes = np.atleast_2d(axes)
ymax = max(h["add"].max() for h in hs21.values()) * 1.35 + 0.5

for i, lid in enumerate(laser_ids):
    ax = axes[i // 4][i % 4]
    h, rr = hs21[lid], rs21[lid]
    ax.fill_between(x, 0, h["add"][b_lo:b_hi + 1], step="mid", color="0.78", label="滤前 hist_add")
    ax.plot(x, rr["after"][b_lo:b_hi + 1], drawstyle="steps-mid",
            color=laser_color(lid), lw=1.8, label="滤后")
    ax.plot(x, rr["thresh"][b_lo:b_hi + 1], drawstyle="steps-mid", color="purple",
            lw=0.9, ls="--", alpha=0.85, label=f"max×{XM_RATIO} 门限")
    n_hit = n_miss = n_kill = 0
    for q in rr["peaks"]:
        kind, _ = peak_truth(h, q)
        if q["is_xtalk"]:
            ax.plot(q["dist"], q["add"], "kx", ms=8, mew=1.6)
            if kind == "纯鬼峰":
                n_hit += 1
            else:
                n_kill += 1
                ax.plot(q["dist"], q["add"], "o", mfc="none", mec="red", ms=14, mew=2)
        elif kind == "纯鬼峰":
            n_miss += 1
    ax.axvline(D_V21, color="k", ls=":", lw=1.0, alpha=0.6)
    ax.set_ylim(0, ymax)
    ax.set_title(f"L{lid}：滤掉鬼影{n_hit}、残留{n_miss}" + (f"、!误杀{n_kill}" if n_kill else ""),
                 fontsize=9)
    ax.tick_params(labelsize=7); ax.grid(alpha=0.2)
    if i == 0:
        ax.legend(fontsize=7, loc="upper left")
for j in range(N_LASERS, nrow * 4):
    axes[j // 4][j % 4].axis("off")

fig.text(0.5, 0.015, f"记录距离 rec_dist [m]（黑点线 = 物体真实距离 {D_V21:.0f}m；"
                     f"黑叉 = 被 XM 丢弃的峰；红圈 = ! 误杀）", ha="center", fontsize=12)
fig.text(0.008, 0.5, "计数（hist_add）", va="center", rotation="vertical", fontsize=12)
plt.suptitle(f"图G  v21 编码下的 XM 滤前/滤后 —— D={D_V21:.0f}m，XM_RATIO={XM_RATIO}"
             f"（对比 v20 图B：那时用 Excel 码，一个都滤不掉）", fontsize=14, y=0.998)
plt.tight_layout(rect=[0.015, 0.035, 1, 0.965])
plt.show()

print_eval(st21, f"（v21 编码, D={D_V21:.0f}m, XM_RATIO={XM_RATIO}）")
'''


# ============================================================================
# TC-6：图 H 距离扫描验收
# ============================================================================
CELL_TC_SWEEP = r'''# ============================================================================
# TC-6  图H：距离扫描验收 —— Excel 现状 vs v21 编码
# ============================================================================
V21_SWEEP_LO, V21_SWEEP_HI, V21_SWEEP_STEP = 5.0, 600.0, 5.0
D_S21 = np.arange(V21_SWEEP_LO, V21_SWEEP_HI + 1e-9, V21_SWEEP_STEP)
V21_CASES = [("Excel 现状", firings_excel, "#c0392b"),
             ("v21 编码",   firings_v21,   "#1e8449")]


def sweep_fr(fr, tag, every=40):
    out = {k: [] for k in ("gb", "ga", "tb", "kill")}
    for i, D in enumerate(D_S21):
        _, _, s = run_with(fr, D)
        P = s["peak"]
        out["gb"].append(sum(P[("纯鬼峰", a)] for a in ("保留", "丢弃")))
        out["ga"].append(P[("纯鬼峰", "保留")])
        out["tb"].append(sum(P[(k, a)] for k in ("纯真峰", "混合峰") for a in ("保留", "丢弃")))
        out["kill"].append(sum(P[(k, "丢弃")] for k in ("纯真峰", "混合峰")))
        if every and i % every == 0:
            print(f"    [{tag}] {i+1}/{len(D_S21)}  D={D:.0f}m ...")
    return {k: np.asarray(v) for k, v in out.items()}


print(f"距离扫描中（{len(D_S21)} 个距离 × {len(V21_CASES)} 套编码）...")
SW21 = {name: sweep_fr(fr, name) for (name, fr, _) in V21_CASES}
print("扫描完成。\n")

fig, axes = plt.subplots(3, 1, figsize=(16, 12), sharex=True)

ax = axes[0]
for (name, _, col) in V21_CASES:
    ax.plot(D_S21, SW21[name]["gb"], "--", color=col, lw=1.3, alpha=0.65,
            label=f"鬼影峰·XM 前（{name}）")
    ax.plot(D_S21, SW21[name]["ga"], "-", color=col, lw=2.4,
            label=f"鬼影峰·XM 后残留（{name}）")
ax.axvline(D_UNAMBIG, color="r", ls=":", lw=1.2)
ax.set_ylim(0, max(SW21[n]["gb"].max() for n, _, _ in V21_CASES) * 1.28)
ax.set_ylabel("鬼影峰数（16 激光器合计）")
ax.set_title("(a) 鬼影峰：XM 滤前 vs 滤后\n"
             "两条虚线不重合是正常的：编码把原本挤在同一 bin 的鬼影拆散了，"
             "「峰数」反而变多 —— 但每个峰只累到 1，全都能滤掉")
ax.legend(fontsize=8.5, ncol=2); ax.grid(alpha=0.25)

ax = axes[1]
for (name, _, col) in V21_CASES:
    r = SW21[name]["ga"] / np.maximum(SW21[name]["gb"], 1) * 100
    ax.plot(D_S21, r, "-", color=col, lw=2.2, label=f"{name}：平均残留率 {r.mean():.2f}%")
ax.axvline(D_UNAMBIG, color="r", ls=":", lw=1.2)
ax.set_ylabel("鬼影残留率 [%]"); ax.set_ylim(-3, 105)
ax.set_title("(b) 鬼影残留率 —— v21 的验收指标（目标：全程 0）")
ax.legend(fontsize=10); ax.grid(alpha=0.25)

ax = axes[2]
ax.plot(D_S21, SW21["Excel 现状"]["tb"], "k--", lw=3.0, alpha=0.45, label="真目标峰·XM 前")
for i, (name, _, col) in enumerate(V21_CASES):
    ax.plot(D_S21, SW21[name]["tb"] - SW21[name]["kill"], "-", color=col,
            lw=2.4 - i * 0.8, label=f"真目标峰·XM 后存活（{name}）")
    ax.plot(D_S21, SW21[name]["kill"], ":", color=col, lw=1.6,
            label=f"! 误杀（{name}）")
ax.axvline(D_UNAMBIG, color="r", ls=":", lw=1.2)
ax.set_xlabel("物体真实距离 D [m]"); ax.set_ylabel("真目标峰数")
ax.set_title("(c) 真目标：编码不动真回波，误杀应恒为 0"
             "（D>300m 时真回波本来就超窗丢失，与编码无关）")
ax.legend(fontsize=8.5, ncol=2); ax.grid(alpha=0.25)

plt.suptitle(f"图H  v21 验收：距离扫描 {V21_SWEEP_LO:.0f}~{V21_SWEEP_HI:.0f}m"
             f"（XM_RATIO={XM_RATIO}，峰宽裕度 {TCODE_SEP_NS}ns，码范围 {_maxtx:.0f}ns）",
             fontsize=14, y=0.995)
plt.tight_layout(rect=[0, 0, 1, 0.97])
plt.show()

print("=" * 88)
print(f"v21 验收表（{len(D_S21)} 个距离点 × 16 激光器）")
print("=" * 88)
print(f"  {'编码':>10} {'鬼影峰(前)':>11} {'鬼影峰(后)':>11} {'滤除率':>9} "
      f"{'真峰(前)':>9} {'误杀':>6} {'误杀率':>8}")
print("  " + "-" * 82)
for (name, _, _c) in V21_CASES:
    gb, ga = SW21[name]["gb"].sum(), SW21[name]["ga"].sum()
    tb, kk = SW21[name]["tb"].sum(), SW21[name]["kill"].sum()
    print(f"  {name:>10} {gb:>11d} {ga:>11d} {(gb-ga)/max(gb,1):>8.3%} "
          f"{tb:>9d} {kk:>6d} {kk/max(tb,1):>7.2%}")

_resid = SW21["v21 编码"]["ga"].sum()
_kill = SW21["v21 编码"]["kill"].sum()
print("\n" + "=" * 88)
print(f"验收结论：鬼影残留 {_resid} 个，真目标误杀 {_kill} 个"
      f"  ->  {'通过 √' if (_resid == 0 and _kill == 0) else '未达标，需重搜码表'}")
print("=" * 88)
'''


# ============================================================================
# TC-7：总结
# ============================================================================
CELL_TC_SUMMARY = f'''# ============================================================================
# TC-7  总结
# ============================================================================
print("=" * 88)
print("crosstalk_sim_v21 总结")
print("=" * 88)
print(f"""
0. 本 notebook 的结构（一路继承，什么都没删）
   - cell 1~10   = v13 原封不动（kick 栅格、记录距离分布、逐回波堆叠柱）
   - 中段        = v20 的 XM（crosstalk mark）实装
   - 本段        = v21 的 tcode 实装

1. v20 留下的问题
   - XM 判据本身没问题，但 Excel 的 tx_trig_dly 每个激光器 4 个 kick 是同一个值
   - 同一条串扰 4 次落同一 bin -> add/max = {{N_ACC}} -> 与真目标无法区分 -> 滤除率 0

2. v21 的编码
   - 编码对象   ：tx_trig_dly（1 ns 步长）
   - 单光子脉宽 ：8 ns，取 {SEP_NS} ns 作为落点最小间隔（峰宽裕度）
   - 码预算     ：{BUDGET_NS} ns，实际用到 {{_maxtx:.0f}} ns（占 kick 间隙 200ns 的 {{_maxtx/200:.0%}}）
   - 求解方式   ：把三类鬼影源写成统一约束后数值搜索，非解析构造
                  （解析构造只保证同 kick 那一类，且要 184 ns）

3. 三类鬼影源都进了约束
   - 同 kick 串扰（Δk=0, a≠b）
   - 跨 kick 混叠（Δk=1, a≠b）
   - 自身混叠    （Δk=1, a=b）—— 解析构造治不了的那一类

4. 验收结果（见 TC-2 自检 与 图H）
   - 编码自检   ：全部鬼影源类别的落点间隔 ≥ {SEP_NS} ns
   - 距离扫描   ：1~600m 鬼影残留 0，真目标误杀 0
   - 误杀恒为 0 的原因：真回波的 rec_tof = 2D/c，与编码无关，判据碰不到它

5. ! 使用前提（两条，务必记住）
   - 码表绑定当前时序表。哪个激光器在哪几个 kick 发光一改，
     必须重跑 docs/tcode/solve_tcode.py 重搜
   - 仿真仍是理想 δ 回波。峰宽是靠约束保证的，不是仿真验证的。
     下一步应把回波展宽（PULSE_W）加进仿真，用 8ns 真实脉宽复核

6. 下一步：FPGA 8 ns 步长随机抖动（打【外来雷达】的串扰）
   - 本版编码是【模组内】编码，对外来雷达无效（对方码表和我们无关）
   - 建议架构：逐 kick 变、但所有激光器共用的抖动 g[k]
       * 同 kick 内部串扰：g[k] − g[k] = 0，完全抵消 -> 不破坏 v21 码表
       * 外来雷达        ：对方发光与我们的 g[k] 无关 -> 每个 kick 落不同 bin -> 被 XM 滤掉
       * 自身混叠        ：不抵消，额外帮忙
   - 若抖动改成【整个 sync 一起平移】，则对 XM 完全无效：
     XM 是在一个 sync 的 4 次 shot 之间比 add/max，sync 内恒定的量它看不见
   - 预算相加：{{_maxtx:.0f}} ns（tcode）+ 抖动范围 ≤ 200 ns
     -> 抖动最多还能用 {{200-_maxtx:.0f}} ns（8ns 步长约 {{int((200-_maxtx)//8)}} 档）
""")
'''


# ============================================================================
# 组装
# ============================================================================
def code_cell(cid, src):
    return {"cell_type": "code", "id": cid, "metadata": {},
            "execution_count": None, "outputs": [],
            "source": src.splitlines(keepends=True)}


def md_cell(cid, src):
    return {"cell_type": "markdown", "id": cid, "metadata": {},
            "source": src.splitlines(keepends=True)}


if not os.path.exists(SRC_NB):
    raise SystemExit(f"找不到 {SRC_NB}，请先运行 build_crosstalk_v20.py")

with open(SRC_NB, encoding="utf-8") as f:
    nb20 = json.load(f)

v20_cells = []
for c in nb20["cells"]:
    c = dict(c)
    if c["cell_type"] == "code":
        c["outputs"] = []
        c["execution_count"] = None
    v20_cells.append(c)

cells = (
    [md_cell("v21_overview", CELL_V21_OVERVIEW)]
    + v20_cells                                          # ← v20 原样继承（内含 v13 原样继承）
    + [
        md_cell("v21_tc_doc",       CELL_TC_DOC),
        code_cell("v21_tc_table",   CELL_TC_TABLE),
        code_cell("v21_tc_check",   CELL_TC_CHECK),
        code_cell("v21_tc_figcode", CELL_TC_FIG_CODE),
        code_cell("v21_tc_figrep",  CELL_TC_FIG_REPEAT),
        code_cell("v21_tc_figgrid", CELL_TC_FIG_GRID),
        code_cell("v21_tc_sweep",   CELL_TC_SWEEP),
        code_cell("v21_tc_summary", CELL_TC_SUMMARY),
    ]
)

nb = {"cells": cells, "metadata": nb20.get("metadata", {}),
      "nbformat": 4, "nbformat_minor": 5}

with open(OUT_NB, "w", encoding="utf-8") as f:
    json.dump(nb, f, ensure_ascii=False, indent=1)

print(f"已生成 {OUT_NB}")
print(f"  继承 v20 的 cell：{len(v20_cells)} 个（原样，未删未改）")
print(f"  v21 新增 cell   ：{len(cells) - len(v20_cells)} 个")
print(f"  合计            ：{len(cells)} 个")
print(f"  码表：{len(TCODE_TABLE)} 项，峰宽裕度 {SEP_NS}ns，最大码值 {MAXTX}ns")
