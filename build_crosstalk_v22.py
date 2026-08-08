# -*- coding: utf-8 -*-
"""
build_crosstalk_v22.py —— 生成 crosstalk_sim_v22.ipynb
========================================================
继承：v21 全部 cell 原样保留（内含 v13 + v20 + v21）。
新增：避真峰 tcode 双门槛码表 + FPGA 对射分析，带清晰分块标题。

缩写：
  XM（XtalkMark，串扰标记）、TOF（Time of Flight，飞行时间）、
  FPGA（Field-Programmable Gate Array，现场可编程门阵列）、
  IRF（Instrument Response Function，仪器响应函数）。
"""
import json
import os

SRC_NB = "crosstalk_sim_v21.ipynb"
OUT_NB = "crosstalk_sim_v22.ipynb"
TBL_R15 = os.path.join("docs", "tcode", "tcode_table_v22_r1.5_56ns.py")
TBL_R25 = os.path.join("docs", "tcode", "tcode_table_v22_r2.5_24ns.py")


def load_tbl(path):
    ns = {}
    with open(path, encoding="utf-8") as f:
        exec(compile(f.read(), path, "exec"), ns)
    return ns["TCODE_TABLE"], ns["TCODE_SEP_NS"], ns["TCODE_BUDGET_NS"]


T15, SEP, B15 = load_tbl(TBL_R15)
T25, _, B25 = load_tbl(TBL_R25)
LASERS = sorted({l for (l, k) in T15})


def lit(tbl):
    lines = []
    kicks_of = {l: sorted(k for (a, k) in tbl if a == l) for l in LASERS}
    for l in LASERS:
        items = ", ".join(f"({l},{k}): {tbl[(l,k)]:>2d}" for k in kicks_of[l])
        lines.append(f"    {items},")
    return "\n".join(lines)


def code_cell(cid, src):
    return {"cell_type": "code", "id": cid, "metadata": {},
            "execution_count": None, "outputs": [],
            "source": src.splitlines(keepends=True)}


def md_cell(cid, src):
    return {"cell_type": "markdown", "id": cid, "metadata": {},
            "source": src.splitlines(keepends=True)}


# ===========================================================================
CELL_OVERVIEW = f'''# 串扰仿真 v22 —— 避真峰 tcode（双门槛）+ FPGA 对射分析

> **本 notebook = v21 原封不动 + 末尾追加 v22。**
> v21 含 v13（时序/回波）+ v20（XM）+ v21（第一版 tcode）。
> **从下面「v22 目录」开始才是本版新内容。**

## v22 相对 v21 的升级

1. **硬约束加严**：鬼影必须离真峰 ≥ {SEP} ns（避真），不能只散开
2. **两套门槛两套码表**（都要留着）：
   - `XM_RATIO = 1.5` → 码预算 **{B15} ns**（滤净实证下界）
   - `XM_RATIO = 2.5` → 码预算 **{B25} ns**（约束硬下界即可滤净）
3. **新增**：不加 tcode、只上 FPGA 随机时，**雷达对射**的残留概率（解析 + 蒙卡）

## v22 目录（清晰分块）

| 块 | 标题 | 做什么 |
|---|---|---|
| A | 路线图与门槛说明 | 1.5 vs 1.6、下界定义 |
| B | 变量命名 | SEP / 双码表 / 开关 |
| C | 导入码表 | 两张 `tx_trig_dly` 表 |
| D | tcode 约束自检 | 散开 + 避真 |
| E | 加了 tcode 的效果 | 距离扫描，1.5 表 vs 2.5 表 |
| F | FPGA 对射（无 tcode） | 解析公式 + 蒙卡曲线 |
| G | 总结与下一步 | 一字滤波预告 |
'''

CELL_A = f'''---
# v22 新增部分

---

# 【A】路线图与门槛说明

## A.1 编码滤噪分几层（回顾）

| 层 | 手段 | 管什么 |
|---|---|---|
| L1 | **tcode**（确定性） | 模组内串扰：散开 + **不撞真峰** |
| L2 | **FPGA 随机 delay** | 外来雷达对射 |
| L3 | 一字滤波（未做） | 角度上孤立的假点 |

## A.2 `XM_RATIO = 1.5` 和 `1.6` 几乎等价——但不是数学恒等

判据：`hist_max × ratio > hist_add` → 丢掉。

计数是**整数**时，常见情况：

| add | max | add/max | ratio=1.5 | ratio=1.6 |
|---|---|---|---|---|
| 1 | 1 | 1.0 | 丢掉 | 丢掉 |
| 2 | 1 | 2.0 | **保留** | **保留** |
| 3 | 2 | 1.5 | **保留** | 丢掉 |

你的直觉对：**大部分被拆散的串扰是单 shot 峰（add=max=1）**，此时 1.00001 和 1.99999 对判决没有区别——只要 `1 < ratio ≤ 2`，行为相同。

差别只出现在「半整数」边界（如 3/2=1.5）：1.6 会多杀一点碰撞峰，1.5 稍松。  
本版按你的要求用 **1.5**（与 1.6 同档），并和 **2.5** 对照。

## A.3 两种「下界」不要混

| 名称 | 含义 | 数值（SEP={SEP}ns） |
|---|---|---|
| **约束硬下界** | 散开+避真在数学上可行 | **24 ns**（再小放不下 4 个合法码差） |
| **ratio=2.5 滤净下界** | 距离扫描残留=0 | **{B25} ns**（本版取约束下界即可） |
| **ratio=1.5 滤净下界** | 距离扫描残留=0 | **{B15} ns**（双碰滤不掉，要更大码空间压生日碰撞） |

> ratio=1.5/1.6 **滤不掉** `add=2,max=1` 的双碰；要残留真正为 0，必须让这种碰撞几乎不出现 → 码预算显著大于 24ns。
'''

CELL_B = f'''# ============================================================================
# 【B】变量命名
# ============================================================================
#   SEP_NS          落点最小间隔 = 避真间隔 [ns]（峰宽裕度）
#   TCODE_R15       专为 XM_RATIO=1.5 滤净的码表（预算 {B15} ns）
#   TCODE_R25       专为 XM_RATIO=2.5 滤净的码表（预算 {B25} ns）
#   V22_RATIO_MODE  "1.5" | "2.5"  —— 当前启用哪张表 / 哪个门槛
#   FPGA_N_LEVELS   对射实验：每档 8ns，可选 0..N-1 共 N 档
# ============================================================================
SEP_NS = {SEP}                 # 单光子脉宽 8ns，取 12ns 裕度
V22_RATIO_MODE = "1.5"       # 改成 "2.5" 就切换到短码表 + 松门槛

# 对射 / FPGA 随机（本块 F 使用）
FPGA_STEP_NS   = 8           # FPGA 步长 [ns]
FPGA_N_LEVELS  = 8           # 默认可选档数；下面会扫 N=2..64
FPGA_N_TRIALS  = 5000        # 蒙卡次数（每个 N）

print("【B】变量命名")
print(f"  SEP_NS          = {{SEP_NS}} ns")
print(f"  V22_RATIO_MODE  = {{V22_RATIO_MODE!r}}  (可选 '1.5' / '2.5')")
print(f"  FPGA_STEP_NS    = {{FPGA_STEP_NS}} ns")
print(f"  FPGA_N_LEVELS   = {{FPGA_N_LEVELS}}  -> 随机范围 0..{{FPGA_STEP_NS*(FPGA_N_LEVELS-1)}} ns")
'''

CELL_C = f'''# ============================================================================
# 【C】导入码表 —— 两张 tx_trig_dly 都留着
# ============================================================================
TCODE_R15 = {{   # XM_RATIO=1.5 滤净，预算 {B15} ns
{lit(T15)}
}}
TCODE_R25 = {{   # XM_RATIO=2.5 滤净，预算 {B25} ns
{lit(T25)}
}}

def tcode_r15(lid, k, tx0=0):
    return TCODE_R15.get((lid, k), tx0)

def tcode_r25(lid, k, tx0=0):
    return TCODE_R25.get((lid, k), tx0)

# 按开关选用
if V22_RATIO_MODE == "2.5":
    TCODE_V22, tcode_v22, V22_XM_RATIO, V22_BUDGET = TCODE_R25, tcode_r25, 2.5, {B25}
else:
    TCODE_V22, tcode_v22, V22_XM_RATIO, V22_BUDGET = TCODE_R15, tcode_r15, 1.5, {B15}

firings_r15   = build_firings_tcode(tcode_r15)
firings_r25   = build_firings_tcode(tcode_r25)
firings_excel = build_firings_tcode(make_tcode_fn("excel"))
FR_V22 = firings_r15 if V22_RATIO_MODE == "1.5" else firings_r25

_hard = (KICK_SPACING - TOF_WINDOW) / NS
print("【C】码表导入与预算核查")
print(f"  当前模式 V22_RATIO_MODE={{V22_RATIO_MODE}}  ->  XM_RATIO={{V22_XM_RATIO}}, 预算 {{V22_BUDGET}} ns")
for name, tbl, br in [("R1.5", TCODE_R15, {B15}), ("R2.5", TCODE_R25, {B25})]:
    mx = max(tbl.values())
    print(f"  {{name}}: max(tx)={{mx}} ns / 预算 {{br}} ns / 硬上限 {{_hard:.0f}} ns  "
          f"-> {{'OK' if mx <= br and mx <= _hard else '!!'}}")
print(f"\\n  当前启用表（前 4 个激光器）：")
for _l in laser_ids[:4]:
    print(f"    L{{_l}}: " + "  ".join(f"K{{k}}={{TCODE_V22[(_l,k)]}}" for k in SHOT_KICKS[_l]))
'''

CELL_D = r'''# ============================================================================
# 【D】tcode 约束自检 —— 散开 + 避真
# ============================================================================
def ghost_classes_v22(gap=CROSSTALK_MAX_GAP, dks=(0, 1)):
    out = []
    for b in laser_ids:
        for a in laser_ids:
            if abs(a - b) > gap:
                continue
            for dk in dks:
                if a == b and dk == 0:
                    continue
                prs = [(kb - dk, kb) for kb in SHOT_KICKS[b]
                       if (kb - dk) in SHOT_KICKS[a]]
                if len(prs) >= 2:
                    kind = ("自身混叠" if a == b else
                            ("同kick串扰" if dk == 0 else "跨kick串扰"))
                    out.append((a, b, dk, kind, prs))
    return out


def check_spread_and_avoid(code_fn, sep=SEP_NS):
    """返回 (散开是否全过, 避真是否全过, 明细)。"""
    spread_fail, avoid_fail = [], []
    for (a, b, dk, kind, prs) in ghost_classes_v22():
        d = sorted(code_fn(a, ka, 0) - code_fn(b, kb, 0) for (ka, kb) in prs)
        gaps = [d[i+1] - d[i] for i in range(len(d)-1)]
        mind = min(gaps) if gaps else 10**9
        if mind < sep:
            spread_fail.append((a, b, dk, kind, d, mind))
        # 避真：仅同 kick（跨 kick 相对真峰还要减 T_kick，码预算内撞不到）
        if dk == 0 and a != b:
            for ka, kb in prs:
                diff = abs(code_fn(a, ka, 0) - code_fn(b, kb, 0))
                if diff < sep:
                    avoid_fail.append((a, b, ka, diff))
    return len(spread_fail) == 0, len(avoid_fail) == 0, spread_fail, avoid_fail


print("【D】约束自检（要求间隔 ≥ %d ns）" % SEP_NS)
for name, fn in [("Excel 现状", make_tcode_fn("excel")),
                 ("R1.5 码表", tcode_r15),
                 ("R2.5 码表", tcode_r25)]:
    ok_s, ok_a, fs, fa = check_spread_and_avoid(fn)
    print(f"  {name}: 散开={'通过' if ok_s else f'失败×{len(fs)}'}  "
          f"避真={'通过' if ok_a else f'失败×{len(fa)}'}")
'''

CELL_E = r'''# ============================================================================
# 【E】加了 tcode 的效果 —— 两套门槛对照扫描
# ============================================================================
V22_SWEEP_LO, V22_SWEEP_HI, V22_SWEEP_STEP = 5.0, 600.0, 5.0
D_S22 = np.arange(V22_SWEEP_LO, V22_SWEEP_HI + 1e-9, V22_SWEEP_STEP)

CASES_E = [
    ("Excel + ratio1.5", firings_excel, 1.5, "#7f8c8d"),
    ("R2.5表 + ratio2.5", firings_r25,   2.5, "#c0392b"),
    ("R1.5表 + ratio1.5", firings_r15,   1.5, "#1e8449"),
]


def sweep_case(fr, ratio, tag, every=40):
    out = {k: [] for k in ("gb", "ga", "tb", "kill")}
    for i, D in enumerate(D_S22):
        with use_firings(fr):
            recs = detect_echoes(D)
        hs = build_hists(recs)
        rs = crosstalk_mark_all(hs, ratio)
        s = evaluate(hs, rs)
        P = s["peak"]
        out["gb"].append(sum(P[("纯鬼峰", a)] for a in ("保留", "丢弃")))
        out["ga"].append(P[("纯鬼峰", "保留")])
        out["tb"].append(sum(P[(k, a)] for k in ("纯真峰", "混合峰") for a in ("保留", "丢弃")))
        out["kill"].append(sum(P[(k, "丢弃")] for k in ("纯真峰", "混合峰")))
        if every and i % every == 0:
            print(f"    [{tag}] {i+1}/{len(D_S22)} D={D:.0f}m")
    return {k: np.asarray(v) for k, v in out.items()}


print("【E】距离扫描中...")
SW22 = {}
for name, fr, ratio, col in CASES_E:
    SW22[name] = sweep_case(fr, ratio, name)
print("扫描完成。\n")

fig, axes = plt.subplots(2, 1, figsize=(15, 9), sharex=True)
ax = axes[0]
for name, fr, ratio, col in CASES_E:
    r = SW22[name]["ga"] / np.maximum(SW22[name]["gb"], 1) * 100
    ax.plot(D_S22, r, "-", color=col, lw=2.2,
            label=f"{name}  平均残留 {r.mean():.2f}%")
ax.axvline(D_UNAMBIG, color="r", ls=":", lw=1.2)
ax.set_ylabel("鬼影残留率 [%]"); ax.set_ylim(-3, 105)
ax.set_title("【E】模组内串扰：两套 tcode × 两种 XM_RATIO")
ax.legend(fontsize=9); ax.grid(alpha=0.25)

ax = axes[1]
for name, fr, ratio, col in CASES_E:
    ax.plot(D_S22, SW22[name]["kill"], "-", color=col, lw=2.0, label=f"{name} 误杀")
ax.axvline(D_UNAMBIG, color="r", ls=":", lw=1.2)
ax.set_xlabel("物体真实距离 D [m]"); ax.set_ylabel("真目标误杀数")
ax.set_title("误杀应恒为 0（避真约束 + 真回波与编码无关）")
ax.legend(fontsize=9); ax.grid(alpha=0.25)
plt.tight_layout(); plt.show()

print("=" * 78)
print("【E】验收汇总")
print(f"  {'方案':<22} {'鬼影前':>8} {'残留':>8} {'滤除率':>9} {'误杀':>6}")
print("  " + "-" * 60)
for name, fr, ratio, col in CASES_E:
    gb, ga = SW22[name]["gb"].sum(), SW22[name]["ga"].sum()
    kk = SW22[name]["kill"].sum()
    print(f"  {name:<22} {gb:>8d} {ga:>8d} {(gb-ga)/max(gb,1):>8.3%} {kk:>6d}")
'''

CELL_F_DOC = r'''# 【F】FPGA 对射分析（不加 tcode）

## F.1 场景

- **不考虑模组内串扰**（假设 tcode 已管住，或本实验故意关掉）
- 只看**外来雷达对射**：对方发光与我们无共享时序
- 我们的 FPGA 给每个 kick 一个随机延时（步长 8ns，共 N 档）
  - 简化模型：外来干扰在我们直方图里的落点 = 在 `M = N` 个离散档上均匀独立抽取
  - （累计随机 walk 会让相邻 kick 相关，但边缘分布仍覆盖这些档；独立模型是偏保守/清晰的一阶近似）

## F.2 XM 何时滤得掉

`N_ACC = 4` 次 shot，外来峰在各档的计数为多重集。

| 碰撞模式（4 次） | max 频数 | ratio=1.5 | ratio=2.5 |
|---|---|---|---|
| 1+1+1+1（全散开） | 1 | 丢掉 | 丢掉 |
| 2+1+1（恰一对碰） | 2 | **保留** | 丢掉 |
| 2+2 / 3+1 / 4 | ≥2 或 ≥3 | 保留 | 视情况 |

所以：

- **ratio=1.5**：只要出现任何双碰就会漏 → `P_survive = 1 - P(全散开)`
- **ratio=2.5**：要三碰及以上才漏 → 比 1.5 安全得多

## F.3 解析式

设 `M = N`（档数），`n = N_ACC = 4`：

$$P_{\text{全散开}} = \frac{M(M-1)(M-2)(M-3)}{M^4}\quad (M\ge 4)$$

$$P_{\text{survive}}(1.5) = 1 - P_{\text{全散开}}$$

$$P_{\text{survive}}(2.5) = P(\text{某档计数}\ge 3)$$

（下面代码里用蒙卡交叉验证解析式。）
'''

CELL_F = r'''# ============================================================================
# 【F】FPGA 对射：解析 + 蒙卡（无模组内串扰）
# ============================================================================
from math import comb

N_ACC_F = N_ACC   # =4


def p_all_distinct(M, n=N_ACC_F):
    if M < n:
        return 0.0
    p = 1.0
    for i in range(n):
        p *= (M - i) / M
    return p


def p_survive_15_analytic(M, n=N_ACC_F):
    """ratio=1.5：任何双碰即残留。"""
    return 1.0 - p_all_distinct(M, n)


def p_max_freq_ge(M, n, thr, n_mc=20000, rng=None):
    """蒙卡：n 次均匀抽 M 档，最大频数 ≥ thr 的概率。"""
    rng = np.random.default_rng(0) if rng is None else rng
    hits = 0
    for _ in range(n_mc):
        bins = rng.integers(0, M, size=n)
        # 频数
        _, cnt = np.unique(bins, return_counts=True)
        if cnt.max() >= thr:
            hits += 1
    return hits / n_mc


def p_survive_25_analytic_approx(M, n=N_ACC_F, n_mc=50000):
    """ratio=2.5：最大频数≥3 才残留。用高精度蒙卡当「解析参考」。"""
    return p_max_freq_ge(M, n, thr=3, n_mc=n_mc)


# ---- 扫 N = 档数 ----
N_LIST = [2, 3, 4, 5, 6, 8, 10, 12, 16, 24, 32, 48, 64]
print("【F】FPGA 对射残留概率（独立均匀模型，N_ACC=%d）" % N_ACC_F)
print(f"  {'N档':>6} {'范围ns':>8} {'P残留@1.5':>12} {'P残留@2.5':>12} {'蒙卡@1.5':>10}")
print("  " + "-" * 56)

rows = []
rng = np.random.default_rng(1)
for N in N_LIST:
    span = FPGA_STEP_NS * (N - 1)
    p15 = p_survive_15_analytic(N)
    p25 = p_survive_25_analytic_approx(N, n_mc=30000)
    p15_mc = p_max_freq_ge(N, N_ACC_F, thr=2, n_mc=FPGA_N_TRIALS, rng=rng)
    rows.append((N, span, p15, p25, p15_mc))
    print(f"  {N:>5d} {span:>7d}ns {p15:>11.3%} {p25:>11.3%} {p15_mc:>9.3%}")

Ns = np.array([r[0] for r in rows], dtype=float)
P15 = np.array([r[2] for r in rows])
P25 = np.array([r[3] for r in rows])
P15mc = np.array([r[4] for r in rows])

fig, ax = plt.subplots(figsize=(12, 5.5))
ax.plot(Ns, P15 * 100, "o-", color="#1e8449", lw=2.2, label="解析 P_survive (ratio=1.5)")
ax.plot(Ns, P15mc * 100, "x", color="#1e8449", ms=8, label="蒙卡 (ratio=1.5, 频数≥2)")
ax.plot(Ns, P25 * 100, "s-", color="#c0392b", lw=2.2, label="蒙卡/参考 P_survive (ratio=2.5, 频数≥3)")
ax.set_xlabel("FPGA 可选档数 N（步长 8ns，范围 = 8(N-1) ns）")
ax.set_ylabel("对射干扰在一次 sync 后仍被 XM 当成「峰」的概率 [%]")
ax.set_title("【F】不加 tcode、只靠 FPGA 随机打对射：档数越多，残留概率越低\n"
             "注意：ratio=1.5 对「双碰」零容忍，同样 N 下比 2.5 难打干净")
ax.set_xscale("log", base=2)
ax.grid(alpha=0.3)
ax.legend(fontsize=9)
plt.tight_layout(); plt.show()

print("""
解读：
  - N=8（范围 56ns）时，ratio=1.5 仍有相当概率因双碰漏网；
  - ratio=2.5 在同样 N 下残留低一个数量级以上（要三碰才漏）；
  - 这与模组内 tcode 的结论一致：1.5 档更吃码空间 / 随机空间。
  - tcode 管不了对射（对方无共享码表）；对射必须靠 FPGA 随机（或一字滤波）。
""")
'''

CELL_G = f'''# ============================================================================
# 【G】总结
# ============================================================================
print("=" * 78)
print("crosstalk_sim_v22 总结")
print("=" * 78)
print(f"""
【结构】
  v13 时序/回波 → v20 XM → v21 第一版 tcode → v22 避真 + 双门槛 + 对射

【A/B/C】变量与码表
  SEP = {SEP} ns（避真 = 散开间隔）
  R1.5 码表：预算 {B15} ns，配合 XM_RATIO=1.5，模组内扫描残留 → 0
  R2.5 码表：预算 {B25} ns，配合 XM_RATIO=2.5，模组内扫描残留 → 0
  两张表都留着，用 V22_RATIO_MODE 切换

【门槛】
  1.5 与 1.6：对「单 shot 鬼影」等价；对 add/max=1.5 的边界 1.6 略严
  1.5 滤不掉 add=2,max=1 的双碰 → 滤净需要更大码空间（本版 {B15} ns）
  2.5 能吃掉双碰 → 滤净可贴近约束硬下界 24 ns（本版 {B25} ns）

【F】对射
  不加 tcode，只靠 FPGA 的 N 档×8ns 随机：
  P_survive(1.5) = 1 - P(4 次全落不同档)
  P_survive(2.5) = P(某档至少 3 次)
  tcode 对对射无效；必须 FPGA 随机（逐 kick 变、建议全模组共用抖动以免拆掉模组内码差）

【未做 / 下一步】
  一字滤波：连续三角度，中间目标相对左右距离差都 > thr → 判鬼影丢掉
  tcode + FPGA 联合、累计随机 walk 精细模型
""")
'''


# ===========================================================================
if not os.path.exists(SRC_NB):
    raise SystemExit(f"找不到 {SRC_NB}")

with open(SRC_NB, encoding="utf-8") as f:
    nb21 = json.load(f)

v21_cells = []
for c in nb21["cells"]:
    c = dict(c)
    if c["cell_type"] == "code":
        c["outputs"] = []
        c["execution_count"] = None
    v21_cells.append(c)

cells = (
    [md_cell("v22_overview", CELL_OVERVIEW)]
    + v21_cells
    + [
        md_cell("v22_A", CELL_A),
        code_cell("v22_B", CELL_B),
        code_cell("v22_C", CELL_C),
        code_cell("v22_D", CELL_D),
        code_cell("v22_E", CELL_E),
        md_cell("v22_F_doc", CELL_F_DOC),
        code_cell("v22_F", CELL_F),
        code_cell("v22_G", CELL_G),
    ]
)

nb = {"cells": cells, "metadata": nb21.get("metadata", {}),
      "nbformat": 4, "nbformat_minor": 5}

with open(OUT_NB, "w", encoding="utf-8") as f:
    json.dump(nb, f, ensure_ascii=False, indent=1)

print(f"已生成 {OUT_NB}")
print(f"  继承 v21 cell：{len(v21_cells)}")
print(f"  v22 新增：{len(cells) - len(v21_cells)}")
print(f"  R1.5 预算 {B15} ns / R2.5 预算 {B25} ns / SEP {SEP} ns")
