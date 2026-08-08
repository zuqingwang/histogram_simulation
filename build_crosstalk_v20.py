# -*- coding: utf-8 -*-
"""
build_crosstalk_v20.py —— 生成 crosstalk_sim_v20.ipynb
========================================================
【构建方式】v20 = v13 的全部 cell **原样继承** + 后面追加 v20 新增 cell。
  本脚本直接读取 crosstalk_sim_v13.ipynb，把它的 10 个 cell 一字不改地搬过来
  （包括用户精心写的 cell5 kick 栅格图、cell9 逐回波堆叠柱），
  只在最前面插入一个 v20 总览 markdown，在最后面追加 v20 新增内容。
  ——— 不删 v13 的任何东西，也不改 v13 任何一行代码。

【v20 新增】crosstalk mark（XM，串扰标记）—— 只实装用户要求的【简化版】：
    hist_add = 一个激光器 N 次 shot 的累加波形
    hist_max = 逐 bin 取 N 次 shot 里的最大值（等价于"某一次的峰"）
    判据：  hist_max × xm_ratio > hist_add   →  该峰是串扰，丢弃
  硬件文档 docs/xm/xm_crosstalk_mark.md 里的分档 / flag0~3 / diff_flag / FIR
  等复杂分支 **本版一律不实装**，只在 markdown 里说明它们是后续可选的加强项。

  另外附带：
    - 滤前 / 滤后 波形都画出来（单激光器全链路解剖 + 16 宫格对比）
    - 统计：多少鬼影、多少串扰、XM 之后变成多少（峰级混淆 + 记录级分类）
    - TCODE_MODE 预览开关（默认 'excel' = 用户现状，不改任何用户参数）
    - 距离扫描 1~600m 的鬼影残留率曲线

本脚本只生成 notebook，不执行、不 debug（按用户要求）。

缩写：
  XM（XtalkMark，串扰标记）
  TOF / ToF（Time of Flight，飞行时间）
  SPAD（Single-Photon Avalanche Diode，单光子雪崩二极管）
  IRF（Instrument Response Function，仪器响应函数）
  FPGA（Field-Programmable Gate Array，现场可编程门阵列）
  LCG（Linear Congruential Generator，线性同余生成器）
"""
import json
import os

SRC_NB = "crosstalk_sim_v13.ipynb"
OUT_NB = "crosstalk_sim_v20.ipynb"


# ============================================================================
# 【插在最前面】v20 总览
# ============================================================================
CELL_V20_OVERVIEW = r'''# 串扰仿真 v20 —— 在 v13 之上加 crosstalk mark（XM）串扰滤除

> **本 notebook = v13 原封不动 + 末尾追加 XM 部分。**
> 下面 cell 1~10 全部是 v13 的原始内容（kick 栅格图、记录距离分布图、逐回波堆叠柱……
> 一个都没删、一行都没改）。v20 的新增内容从 **「v20 新增部分」** 那个标题往后开始。

## v20 做了什么

| | v13 | v20 新增 |
|---|---|---|
| 输出形式 | 回波**记录列表**（谁发的、被谁收到、算成多远） | 在记录之上再落成**真正的直方图波形** |
| 每个激光器 | 4 次发光的回波混在一起看 | 4 次 shot **分别留存** + 累加波形 `hist_add` |
| 算法 | 无（只忠实记录） | **crosstalk mark（XM）**：把串扰峰标出来并丢弃 |
| 结果 | 鬼影有多少 | 鬼影 **XM 前有多少 / XM 后剩多少 / 误杀了几个真目标** |

## 新增 cell 一览

| cell | 内容 |
|---|---|
| 「v20 新增部分」md | **XM 算法逐步讲解**（含具体数字例子），必读 |
| XM-1 | 直方图层：逐 shot 波形 → `hist_add` / `hist_max` |
| XM-2 | XM 核心（只有十几行，就是用户要求的那一条判据） |
| XM-3 | 统计：峰级混淆矩阵 + 记录级分类 |
| XM-4 | 图A 单激光器 XM 全链路解剖（4 次 shot → 判据 → 滤前/滤后） |
| XM-5 | 图B 16 个激光器滤前/滤后对比 |
| XM-6 | 图C 统计柱状图 + tcode 对比表 |
| XM-7 | 图D 距离扫描 1~600m 的鬼影残留率 |
| XM-8 | 总结 + v21（tcode 编码）预告 |

## ⚠ 先说结论

用 Excel 现有的 `tx_trig_dly`（每个激光器 4 个 kick 的码是**同一个值**）跑，
**XM 几乎滤不掉任何鬼影**。这不是 bug，是必然结果，原因在下面的 md 里推导。
要让 XM 生效，必须先做 **tcode（各激光器之间的发光时刻编码）**——那是 v21 的事，
原理和图解见 `docs/tcode/tcode_scheme.md`。
本版提供 `TCODE_MODE = "lcg2"` 预览开关（**默认仍是 `"excel"`，不改用户参数**），
它就是 tcode 文档里推荐的两层编码，用来证明 XM 本身实现是对的。
'''


# ============================================================================
# 【追加】v20 新增部分总说明（XM 算法详解）
# ============================================================================
CELL_XM_DOC = r'''---

# v20 新增部分

---

# crosstalk mark（XM）算法详解

## 0. 一句话

**串扰只出现在个别某一次发光里，真目标每一次发光都在。**

---

## 1. 我们手上有什么波形

一个激光器（比如 L5）在一个 sync 内发 **4 次光**（4 个 kick）。
v13 已经算出「每一次发光收到了哪些回波、各自算成多少米」，
v20 把它们**按 1 ns 一个 bin 落成波形**，于是每个激光器手上有：

| 波形 | 怎么来的 | 硬件里叫什么 |
|---|---|---|
| `shot[0..3]` | 第 0~3 次发光**各自**的 TOF 波形，**分别留存** | 单次测量结果 |
| `hist_add` | 4 条 shot **逐 bin 相加** | `hist_add` |
| `hist_max` | 4 条 shot **逐 bin 取最大值** | `hist_max` |

---

## 2. 判据（就是你要求的那一条）

$$\boxed{\;\text{hist\_max}[b] \times \texttt{xm\_ratio} \;>\; \text{hist\_add}[b]
\;\Longrightarrow\; \text{这个峰是串扰，丢弃}\;}$$

`xm_ratio` 默认取 **1.6**。

### 为什么只要存 `hist_max`，不用真的比较 4 条 shot

你的原话是「**某次**的某个峰 × ratio > hist_add 就判串扰」。
对同一个 bin，4 次里**只要有一次**满足就成立；而 4 个数里最容易满足的显然是**最大的那个**。
所以

$$\exists s:\; \text{shot}_s[b]\times r > \text{hist\_add}[b]
\quad\Longleftrightarrow\quad
\max_s \text{shot}_s[b]\times r > \text{hist\_add}[b]$$

**逐 bin 取最大值就够了 —— 这就是「maxhold」这个名字的由来**，
也是硬件只额外存一份 `hist_max` 而不是存 4 份波形的原因。
（本 notebook 里 4 条 shot 也全都留着，方便画图和核对。）

---

## 3. 一个具体数字例子

L5 的 4 次 shot 里出现了两个东西：

| bin | 是什么 | shot0 | shot1 | shot2 | shot3 | `hist_add` | `hist_max` |
|---|---|---|---|---|---|---|---|
| 1000（=150m） | 真目标回波 | 1 | 1 | 1 | 1 | **4** | **1** |
| 950（=142.5m） | 串扰（只有第 2 次撞上） | 0 | 0 | 1 | 0 | **1** | **1** |

代入判据（`xm_ratio = 1.6`）：

| bin | `hist_max × 1.6` | `hist_add` | 谁大 | 判决 |
|---|---|---|---|---|
| 1000 | 1 × 1.6 = **1.6** | **4** | add 大 | 不是串扰 → **保留** ✔ |
| 950 | 1 × 1.6 = **1.6** | **1** | max×r 大 | 是串扰 → **丢弃** ✔ |

---

## 4. `xm_ratio` 到底该取多少（**最重要的一节**）

把判据两边同除 `hist_max`：

$$\text{判串扰} \iff \frac{\text{hist\_add}[b]}{\text{hist\_max}[b]} < \texttt{xm\_ratio}$$

而 `hist_add / hist_max` 有非常清楚的物理含义 —— **这个峰在几次 shot 里出现过**：

| 情况 | `hist_add / hist_max` |
|---|---|
| 真目标（4 次全在） | **4**（= `N_ACC`） |
| 串扰只在 1 次里出现 | **1** |
| 串扰在 2 次里落到**同一个 bin** | **2** |
| 串扰在 k 次里落到同一个 bin | **k** |

于是 `xm_ratio` 的取值必须夹在中间：

$$\boxed{\;\underbrace{k_{\max}}_{\text{串扰最多重复几次}} \;<\; \texttt{xm\_ratio} \;<\; \underbrace{N_{\rm ACC}}_{\text{累加次数}}\;}$$

| `xm_ratio` | 能滤掉 | 真目标 | 说明 |
|---|---|---|---|
| 1.6（本版默认） | 只出现 1 次的串扰 | 安全（4 > 1.6，裕度 2.5×） | 你给的值 |
| 2.5 | 出现 ≤2 次的串扰 | 安全（4 > 2.5，裕度 1.6×） | 更狠，但保护裕度变小 |
| 3.5 | 出现 ≤3 次的串扰 | 危险 | 真目标少亮一次就被杀 |
| ≥ 4 | — | **真目标也被滤** | 禁止 |

**所以 XM 的全部能力上限，就是这条不等式的宽度，而它由累加次数 N_ACC 决定。**

---

## 5. ⚠ 由此推出：光有 XM 是不够的

串扰的落点（v13 已推导）：

$$\text{rec\_tof}
= \underbrace{\frac{2D}{c}}_{\text{物体距离}}
+ \underbrace{\big(tx_{\text{发}} - tx_{\text{收}}\big)}_{\textbf{码差}}
+ \underbrace{(k_{\text{发}} - k_{\text{收}})\cdot T_{\rm kick}}_{\text{跨 kick 混叠}}$$

而**真回波**：`rec_tof = 2D/c`，**与编码无关**（发射时刻和测距参考零点一起平移，差值抵消）。

现在看 Excel 的现状：**每个激光器在它的 4 个 kick 里 `tx_trig_dly` 是同一个值**（0 或 50 恒定）。
于是对固定的一对 (发, 收)：

- 码差 `tx_发 − tx_收` 每一次都**一模一样**
- ⟹ 同一条串扰 4 次都落到**同一个 bin**
- ⟹ `hist_add = 4`、`hist_max = 1`、`add/max = 4`
- ⟹ **和真目标数值上完全相同，XM 无法区分 → 滤不掉**

**结论：XM 是收割器，tcode 才是播种机。**
tcode 的任务就是让同一条串扰在 4 次 shot 里落到 **4 个不同的 bin**，
把 `add/max` 从 4 打到 1，XM 才有东西可收。

> 详细的 tcode 原理、构造方法和图解见 **`docs/tcode/tcode_scheme.md`**（v21 实装）。
> 本版可以把 `TCODE_MODE` 改成 `"lcg2"`（文档里推荐的两层编码）预览这个效果。

---

## 6. 本版**没有**实装的东西（都在硬件文档里，属于后续加强项）

`docs/xm/xm_crosstalk_mark.md` 里的完整硬件 XM 还有一堆分支，
本版**一律不实装**，避免把简单问题复杂化：

| 硬件分支 | 作用 | 本版为什么不做 |
|---|---|---|
| `XM_baseline_ratio` | 补偿单次基线与累加基线的差 | 理想 δ 模型本底恒为 0，没有基线要补 |
| 按 peak 值分 3 个档位 | 强弱信号用不同阈值 | 只有一条判据时没有分档的必要 |
| `flag0` 的三个计数分支 | 数片段里被标记点的个数 | δ 回波峰宽只有 1 bin，无「个数」可数 |
| `flag2` 峰位比较 | 防止把前后排的两个真目标误判 | 峰宽 1 bin 时峰位恒重合，恒成立 |
| `flag3` / `hist_diff_flag` | hdc 送来的强制通道 | 本模型没有 hdc |
| `hist_max` 每 4 bin 抽 1 | 省存储 | 会把单点 max 摊开到 4 个 bin，δ 模型下反而制造假标记 |
| FIR 平滑 | 真实波形去噪 | 理想波形无噪声 |

需要时再逐条加回来即可，判据主体（第 2 节那一条）不会变。
'''


# ============================================================================
# XM-1：直方图层
# ============================================================================
CELL_XM_HIST = r'''# ============================================================================
# XM-1  直方图层：把 v13 的"回波记录"落成真正的波形
# ============================================================================
#   每个接收激光器 L 有 N_ACC 次 shot（= 它的 4 个 kick），每次 shot 一条波形；
#   hist_add = 逐 bin 相加     hist_max = 逐 bin 取 4 次里的最大值
# ============================================================================
from collections import Counter
from contextlib import contextmanager

# ---- 直方图参数 ----
HIST_BIN_NS = 1.0                                       # bin 宽 [ns]（1ns 往返 = 15cm）
N_BINS      = int(round(TOF_WINDOW / NS / HIST_BIN_NS)) # = 2000 个 bin，覆盖 0~300m
HIST_BIN_M  = HIST_BIN_NS * NS * C_LIGHT / 2.0          # 1 bin 对应距离 [m]（≈0.1499）

# ---- 回波幅度：本版所有回波都记 1 个计数（最坏情况：串扰和真信号一样强）----
#   留一个按激光器编号间隔衰减的接口，默认不启用。
AMP_MODE   = "unit"                                     # "unit" | "gap"
AMP_BY_GAP = {0: 1.0, 1: 0.50, 2: 0.25}

# ---- 直方图里放不放"空间可忽略"的串扰（编号间隔 > CROSSTALK_MAX_GAP）----
#   False（默认）= 认为这些光根本进不了探测器，与 v13 的 PLOT_NEGLIGIBLE=False 一致。
XTALK_USE_NEGLIGIBLE = False

# ---- 每个激光器的 shot 索引（哪一次 kick 是第几次 shot）；与编码无关，只算一次 ----
SHOT_KICKS = {lid: sorted(k for (l, k, tx) in laser_fires_raw if l == lid) for lid in laser_ids}
SHOT_IDX   = {lid: {k: i for i, k in enumerate(ks)} for lid, ks in SHOT_KICKS.items()}
N_ACC      = int(np.median([len(ks) for ks in SHOT_KICKS.values()]))   # 累加次数


def echo_amp(rec):
    """一条回波在直方图里贡献多少计数。"""
    if AMP_MODE == "unit":
        return 1.0
    return AMP_BY_GAP.get(abs(rec["emit_laser"] - rec["recv_laser"]), 0.0)


def rec_bin(rec):
    """回波落在第几个 bin。"""
    return int(np.clip(np.floor(rec["rec_tof"] / NS / HIST_BIN_NS), 0, N_BINS - 1))


def build_hists(recs):
    """把 detect_echoes() 的记录列表落成波形。

    返回 {laser: {
        shots : (N_ACC, N_BINS)  每次 shot 的波形（全部留存）
        add   : (N_BINS,)        hist_add = 逐 bin 相加
        max   : (N_BINS,)        hist_max = 逐 bin 取最大
        src   : {bin: [rec,...]} 每个 bin 里都有哪些回波（做统计/溯源用）
        kicks : [kick,...]       shot 序号 -> kick 号
    }}
    """
    H   = {lid: np.zeros((len(SHOT_KICKS[lid]), N_BINS)) for lid in laser_ids}
    src = {lid: defaultdict(list) for lid in laser_ids}
    for r in recs:
        if (not r["is_true"]) and r["negligible"] and (not XTALK_USE_NEGLIGIBLE):
            continue                                    # 空间可忽略的光进不了探测器
        lid = r["recv_laser"]
        b   = rec_bin(r)
        H[lid][SHOT_IDX[lid][r["recv_kick"]], b] += echo_amp(r)
        src[lid][b].append(r)
    return {lid: {"shots": H[lid], "add": H[lid].sum(axis=0), "max": H[lid].max(axis=0),
                  "src": src[lid], "kicks": SHOT_KICKS[lid]} for lid in laser_ids}


# ---- tcode 预览接口 ----
#   v21 才正式做编码设计。这里提供几个开关，用来验证"XM 本身实现是对的、瓶颈在编码"。
#   原理与完整图解见 docs/tcode/tcode_scheme.md
#
#   "excel"  = 用 Excel 原值 0/50ns —— **默认，不改用户任何参数**
#   "none"   = 全 0（完全不编码，对照组）
#   "random" = 每个 (laser, kick) 独立随机码
#   "lcg"    = 第一层 线性同余码   c = ((l·(k+1)) mod P)·step
#              → 保证任意一对激光器的【同 kick】码差在各 kick 上互不相同
#   "lcg2"   = 第一层 + 第二层（docs/tcode 推荐方案）
#              c = ((l·(k+1)) mod P)·step + ((k²) mod Pg)·gstep
#              第二项对所有激光器相同：同 kick 码差里抵消（不破坏第一层），
#              自身混叠码差里不抵消（治第一层治不了的那一类）
TCODE_MODE   = "excel"
TCODE_SEED   = 20
TCODE_STEP   = 8        # 码步长 [ns]，必须 ≥ 回波峰宽，否则码不同也落同一 bin
TCODE_LEVELS = 16       # random 模式码级数 → 码范围 0 ~ (LEVELS-1)×STEP ns
TCODE_P      = 17       # 第一层素数模（需 > 激光器数，且与编号差互素）
TCODE_PG     = 9        # 第二层模
TCODE_GSTEP  = 8        # 第二层步长 [ns]
#   码预算：max(tx) 必须 ≤ KICK_SPACING − TOF_WINDOW = 200ns，否则窗会伸进下一个 kick
TCODE_BUDGET_NS = (KICK_SPACING - TOF_WINDOW) / NS


def make_tcode_fn(mode, seed=TCODE_SEED):
    """返回 f(laser_id, kick, tx_excel) -> tx_trig_dly [ns]。"""
    if mode == "excel":
        return lambda lid, k, tx0: tx0
    if mode == "none":
        return lambda lid, k, tx0: 0
    if mode == "random":
        rng = np.random.default_rng(seed)
        tbl = {(lid, k): int(rng.integers(0, TCODE_LEVELS)) * TCODE_STEP
               for (lid, k, tx0) in laser_fires_raw}
        return lambda lid, k, tx0: tbl[(lid, k)]
    if mode == "lcg":
        return lambda lid, k, tx0: ((lid * (k + 1)) % TCODE_P) * TCODE_STEP
    if mode == "lcg2":
        return lambda lid, k, tx0: (((lid * (k + 1)) % TCODE_P) * TCODE_STEP
                                    + ((k * k) % TCODE_PG) * TCODE_GSTEP)
    raise ValueError(f"未知 TCODE_MODE: {mode}")


def check_tcode_budget(mode):
    """码预算检查：max(tx) 超过 kick 间隙就会挤掉后一个 kick 的 TOF 窗。"""
    fn = make_tcode_fn(mode)
    mx = max(fn(lid, k, tx0) for (lid, k, tx0) in laser_fires_raw)
    ok = mx <= TCODE_BUDGET_NS
    print(f"  tcode '{mode:>6}'：max(tx) = {mx:>4d} ns  /  预算 {TCODE_BUDGET_NS:.0f} ns"
          f"  → {'OK' if ok else '!! 超预算，会挤掉下一个 kick 的窗'}")
    return ok


def build_firings_tcode(tcode_fn):
    """用给定 tcode 重建发光事件表（复用 v13 的 fire_time / ref_time，不改 v13 代码）。"""
    fs = []
    for (lid, k, tx0) in laser_fires_raw:
        tx = tcode_fn(lid, k, tx0)
        fs.append({"laser": lid, "kick": k, "tx": tx,
                   "t_fire": fire_time(lid, k, tx), "t_ref": ref_time(lid, k, tx)})
    fs.sort(key=lambda e: e["t_fire"])
    return fs


@contextmanager
def use_firings(fr):
    """临时把全局 firings 换成 fr，让 v13 的 detect_echoes() 直接可用。用完自动还原。"""
    global firings
    old = firings
    firings = fr
    try:
        yield
    finally:
        firings = old


# ---- 自检：看 D_CUMUL（v13 里用的 150m）处 L1 的波形 ----
XM_DEMO_D = D_CUMUL                       # 沿用 v13 cell9 的距离，便于对照
hists = build_hists(detect_echoes(XM_DEMO_D))

print("tcode 码预算检查（max(tx) 必须 ≤ kick 间隙，否则窗会伸进下一个 kick）：")
for _m in ["excel", "none", "random", "lcg", "lcg2"]:
    check_tcode_budget(_m)
print()

_lid = laser_ids[0]
_h   = hists[_lid]
_nz  = np.flatnonzero(_h["add"] > 1e-9)
print(f"直方图：{N_BINS} 个 bin × {HIST_BIN_NS:.0f}ns，每激光器 N_ACC = {N_ACC} 次 shot")
print(f"\nD={XM_DEMO_D:.0f}m 时 L{_lid} 的非零 bin：")
print(f"  {'bin':>5} {'距离[m]':>8} | " + " ".join(f"shot{i}" for i in range(len(_h['kicks'])))
      + f" | {'add':>5} {'max':>5} {'add/max':>8}  来源")
for b in _nz:
    col = " ".join(f"{_h['shots'][i, b]:>5.0f}" for i in range(_h["shots"].shape[0]))
    a, m = _h["add"][b], _h["max"][b]
    who = ",".join(sorted(set(("自" if r["is_true"] else "鬼")
                              + f"L{r['emit_laser']}K{r['emit_kick']}" for r in _h["src"][b])))
    print(f"  {b:>5d} {b*HIST_BIN_M:>8.2f} | {col} | {a:>5.0f} {m:>5.0f} {a/max(m,1e-9):>8.2f}  {who}")
print(f"\n  → add/max = 这个峰在几次 shot 里出现过。真目标应为 {N_ACC}，串扰应为 1。")
'''


# ============================================================================
# XM-2：XM 核心
# ============================================================================
CELL_XM_CORE = r'''# ============================================================================
# XM-2  crosstalk mark 核心 —— 就是那一条判据，十几行
# ============================================================================
#     hist_max[b] × XM_RATIO > hist_add[b]   →   b 处这个峰是串扰，丢弃
#
#   等价形式： hist_add[b] / hist_max[b] < XM_RATIO
#   物理含义： add/max = 该峰在几次 shot 里出现过；真目标 = N_ACC，串扰 = 1
#   取值约束： 串扰最大重复次数 < XM_RATIO < N_ACC
# ============================================================================
XM_RATIO = 1.6            # 用户设定值（xm_ratio）。改这个就能看到判据松紧的变化


def find_peaks(y, th=0.5):
    """把波形切成峰：连续 > th 的 bin 段。返回 [(start, end), ...]（闭区间）。
       理想 δ 模型下每个非零 bin 自成一个峰；把回波展宽后这里会自动合并成宽峰。"""
    above = np.asarray(y) > th
    if not above.any():
        return []
    idx = np.flatnonzero(above)
    brk = np.flatnonzero(np.diff(idx) > 1)
    return list(zip(np.concatenate(([idx[0]], idx[brk + 1])).tolist(),
                    np.concatenate((idx[brk], [idx[-1]])).tolist()))


def crosstalk_mark(h, ratio=None):
    """对一个激光器的直方图跑 XM。

    返回 dict：
      thresh  : hist_max × ratio          —— 画图用的判据门限曲线
      peaks   : 每个峰的判决明细
      after   : XM 之后的波形（被判串扰的峰整段清零）
      dropped : 被丢掉的峰
    """
    ratio = XM_RATIO if ratio is None else ratio
    add, mx = h["add"], h["max"]
    thresh = mx * ratio                              # ← 判据左边

    after, peaks = add.copy(), []
    for (s, e) in find_peaks(add):
        p        = s + int(np.argmax(add[s:e + 1]))  # 峰顶所在 bin
        add_p    = float(add[p])
        max_p    = float(mx[p])
        is_xtalk = bool(max_p * ratio > add_p)       # ★ 判据 ★
        if is_xtalk:
            after[s:e + 1] = 0.0                     # 丢弃整个峰
        peaks.append({"s": s, "e": e, "peak_bin": p, "dist": p * HIST_BIN_M,
                      "add": add_p, "max": max_p, "ratio": add_p / max(max_p, 1e-9),
                      "is_xtalk": is_xtalk})
    return {"thresh": thresh, "peaks": peaks, "after": after,
            "dropped": [q for q in peaks if q["is_xtalk"]]}


def crosstalk_mark_all(hs, ratio=None):
    return {lid: crosstalk_mark(h, ratio) for lid, h in hs.items()}


xm_res = crosstalk_mark_all(hists)

# ---- 自检：逐峰判决表 ----
print("=" * 92)
print(f"XM 判决表（D={XM_DEMO_D:.0f}m，XM_RATIO={XM_RATIO}，N_ACC={N_ACC}，"
      f"TCODE_MODE='{TCODE_MODE}'）")
print("=" * 92)
print(f"  {'接收器':>5} {'峰距离[m]':>10} {'hist_add':>9} {'hist_max':>9} "
      f"{'max×ratio':>10} {'add/max':>8} {'判决':>6}  该峰里有什么")
print("  " + "-" * 88)
for lid in laser_ids:
    for q in xm_res[lid]["peaks"]:
        who = Counter()
        for b in range(q["s"], q["e"] + 1):
            for r in hists[lid]["src"].get(b, []):
                who["真回波" if r["is_true"] else "鬼影"] += 1
        who_s = " + ".join(f"{k}×{v}" for k, v in sorted(who.items()))
        print(f"  L{lid:>4d} {q['dist']:>10.2f} {q['add']:>9.1f} {q['max']:>9.1f} "
              f"{q['max']*XM_RATIO:>10.2f} {q['ratio']:>8.2f} "
              f"{('串扰→丢' if q['is_xtalk'] else '保留'):>6}  {who_s}")
'''


# ============================================================================
# XM-3：统计
# ============================================================================
CELL_XM_STAT = r'''# ============================================================================
# XM-3  统计：多少鬼影、多少串扰，XM 之后变成多少
# ============================================================================
#   峰级真值（把峰里所有回波翻出来看）：
#     纯真峰 —— 只含真回波           → XM 若丢掉它 = 【误杀】，最严重的错误
#     混合峰 —— 真回波和鬼影叠一起   → XM 若丢掉它 也是误杀（真目标一起没了）
#     纯鬼峰 —— 只含鬼影             → XM 丢掉它 = 【命中】；保留 = 【漏滤】
#   记录级分类（沿用 v13 的判据）：
#     真回波 / 同kick串扰 / 跨kick串扰 / 自身混叠
# ============================================================================
def rec_kind(r):
    """一条回波属于哪一类。"""
    if r["is_true"]:
        return "真回波"
    if r["emit_laser"] == r["recv_laser"]:
        return "自身混叠"                       # 我自己别的 kick 的光超窗折回来
    return "同kick串扰" if r["emit_kick"] == r["recv_kick"] else "跨kick串扰"

GHOST_KINDS = ["同kick串扰", "跨kick串扰", "自身混叠"]
ALL_KINDS   = ["真回波"] + GHOST_KINDS


def peak_truth(h, q):
    """一个峰的真值类别 + 它包含的所有回波。"""
    recs = [r for b in range(q["s"], q["e"] + 1) for r in h["src"].get(b, [])]
    n_true = sum(1 for r in recs if r["is_true"])
    if n_true and len(recs) > n_true:
        return "混合峰", recs
    return ("纯真峰" if n_true else "纯鬼峰"), recs


def evaluate(hs, res):
    """汇总统计。"""
    st = {"peak": Counter(), "rec_before": Counter(), "rec_after": Counter(),
          "per_laser": defaultdict(Counter), "residual_dist": [], "killed_dist": []}
    for lid, r in res.items():
        for q in r["peaks"]:
            kind, recs = peak_truth(hs[lid], q)
            act = "丢弃" if q["is_xtalk"] else "保留"
            st["peak"][(kind, act)] += 1
            st["per_laser"][lid][(kind, act)] += 1
            for rec in recs:
                st["rec_before"][rec_kind(rec)] += 1
                if act == "保留":
                    st["rec_after"][rec_kind(rec)] += 1
            if kind == "纯鬼峰" and act == "保留":
                st["residual_dist"].append(q["dist"])
            if kind in ("纯真峰", "混合峰") and act == "丢弃":
                st["killed_dist"].append(q["dist"])
    return st


def print_eval(st, title=""):
    P, KINDS, ACTS = st["peak"], ["纯真峰", "混合峰", "纯鬼峰"], ["保留", "丢弃"]
    n_ghost = sum(P[("纯鬼峰", a)] for a in ACTS)
    n_true  = sum(P[(k, a)] for k in ("纯真峰", "混合峰") for a in ACTS)
    hit     = P[("纯鬼峰", "丢弃")]
    miss    = P[("纯鬼峰", "保留")]
    killed  = sum(P[(k, "丢弃")] for k in ("纯真峰", "混合峰"))

    print("=" * 84)
    print(f"XM 统计 {title}")
    print("=" * 84)
    print("\n【峰级混淆矩阵】")
    print(f"  {'真值':>8} {'保留':>7} {'丢弃':>7} {'小计':>7}")
    print(f"  {'-'*8} {'-'*7} {'-'*7} {'-'*7}")
    for k in KINDS:
        print(f"  {k:>8} {P[(k,'保留')]:>7d} {P[(k,'丢弃')]:>7d} "
              f"{P[(k,'保留')]+P[(k,'丢弃')]:>7d}")

    print("\n【关键指标】")
    print(f"  鬼影峰   ：XM 前 {n_ghost:>4d}  →  XM 后残留 {miss:>4d}   "
          f"（正确滤除 {hit}，滤除率 {hit/max(n_ghost,1):.1%}）")
    print(f"  真目标峰 ：XM 前 {n_true:>4d}  →  XM 后存活 {n_true-killed:>4d}   "
          f"（⚠ 误杀 {killed}，误杀率 {killed/max(n_true,1):.1%}）")

    print("\n【记录级分类】（一条记录 = 一次「某激光器某 kick 的光被我收到」）")
    print(f"  {'类别':>10} {'XM前':>7} {'XM后':>7} {'消除':>7} {'消除率':>8}")
    print(f"  {'-'*10} {'-'*7} {'-'*7} {'-'*7} {'-'*8}")
    for k in ALL_KINDS:
        b, a = st["rec_before"][k], st["rec_after"][k]
        print(f"  {k:>10} {b:>7d} {a:>7d} {b-a:>7d} {(b-a)/max(b,1):>7.1%}")
    gb = sum(st["rec_before"][k] for k in GHOST_KINDS)
    ga = sum(st["rec_after"][k] for k in GHOST_KINDS)
    print(f"  {'-'*10} {'-'*7} {'-'*7} {'-'*7} {'-'*8}")
    print(f"  {'鬼影合计':>10} {gb:>7d} {ga:>7d} {gb-ga:>7d} {(gb-ga)/max(gb,1):>7.1%}")

    if st["residual_dist"]:
        print(f"\n  残留鬼影峰所在距离：{sorted(set(round(d,1) for d in st['residual_dist']))} m")
    if st["killed_dist"]:
        print(f"  ⚠ 被误杀的真目标峰：{sorted(set(round(d,1) for d in st['killed_dist']))} m")


xm_stat = evaluate(hists, xm_res)
print_eval(xm_stat, f"（D={XM_DEMO_D:.0f}m, TCODE_MODE='{TCODE_MODE}', XM_RATIO={XM_RATIO}）")

_ng = sum(xm_stat["peak"][("纯鬼峰", a)] for a in ("保留", "丢弃"))
if _ng and xm_stat["peak"][("纯鬼峰", "丢弃")] == 0:
    print("\n" + "!" * 84)
    print("!! 一个鬼影都没滤掉 —— 这【不是 bug】，是必然结果。")
    print("!! Excel 里每个激光器 4 个 kick 的 tx_trig_dly 是同一个值 ⟹ 码差恒定")
    print("!! ⟹ 同一条串扰 4 次都落到同一个 bin ⟹ hist_add 也累到 4")
    print(f"!! ⟹ add/max = {N_ACC} > XM_RATIO = {XM_RATIO} ⟹ XM 认为它是真目标。")
    print("!! 验证：把 TCODE_MODE 改成 'lcg2' 重跑本 cell 之后的内容，鬼影会被打散。")
    print("!! 正式的编码设计见 docs/tcode/tcode_scheme.md（v21 实装）。")
    print("!" * 84)
'''


# ============================================================================
# XM-4：图 A 单激光器全链路
# ============================================================================
CELL_XM_FIG_CHAIN = r'''# ============================================================================
# XM-4  图A：单个激光器的 XM 全链路解剖（4 次 shot -> 判据 -> 滤前/滤后）
# ============================================================================
def occupied_range(hs, pad=20):
    lo, hi = N_BINS, 0
    for h in hs.values():
        nz = np.flatnonzero(h["add"] > 1e-9)
        if nz.size:
            lo, hi = min(lo, int(nz[0])), max(hi, int(nz[-1]))
    return (0, N_BINS - 1) if lo > hi else (max(0, lo - pad), min(N_BINS - 1, hi + pad))


# 选鬼影最多的激光器来看
XM_DEMO_LASER = max(laser_ids, key=lambda l: sum(
    1 for b, rs in hists[l]["src"].items() for r in rs if not r["is_true"]))
h, rr = hists[XM_DEMO_LASER], xm_res[XM_DEMO_LASER]
b_lo, b_hi = occupied_range({XM_DEMO_LASER: h}, pad=25)
x = np.arange(b_lo, b_hi + 1) * HIST_BIN_M

fig, axes = plt.subplots(3, 1, figsize=(15, 10.5), sharex=True)

# --- (a) 4 次 shot 各自的波形（全部留存）---
#   δ 回波只有 1 个 bin 宽，4 条曲线叠在一起根本看不出区别，
#   所以每条 shot 单独占一"泳道"（纵向平移 LANE），并用竖线+圆点画出每个回波。
ax = axes[0]
LANE = float(h["shots"].max()) * 1.6 + 0.5           # 泳道间距
for si, kick in enumerate(h["kicks"]):
    y0 = si * LANE
    col = plt.cm.viridis(si / max(len(h["kicks"]) - 1, 1))
    ax.axhline(y0, color="0.85", lw=0.8, zorder=0)
    row = h["shots"][si]
    for b in np.flatnonzero(row[b_lo:b_hi + 1] > 0) + b_lo:
        ax.vlines(b * HIST_BIN_M, y0, y0 + row[b], color=col, lw=2.2, zorder=3)
        ax.plot(b * HIST_BIN_M, y0 + row[b], "o", color=col, ms=6, zorder=4)
        ax.annotate(f"{row[b]:.0f}", (b * HIST_BIN_M, y0 + row[b]), fontsize=7,
                    ha="center", va="bottom", xytext=(0, 2), textcoords="offset points")
ax.set_yticks([si * LANE for si in range(len(h["kicks"]))])
ax.set_yticklabels([f"shot{si}\n(kick {k})" for si, k in enumerate(h["kicks"])], fontsize=8)
ax.set_ylim(-0.3 * LANE, len(h["kicks"]) * LANE)
ax.set_ylabel("每次 shot 各占一条泳道")
ax.set_title(f"(a) L{XM_DEMO_LASER} 的 {len(h['kicks'])} 次发光【分别】的波形"
             f"（每次一条泳道；柱高 = 该 bin 的计数）")
ax.grid(alpha=0.25, axis="x")

# --- (b) 判据：hist_max × ratio  vs  hist_add ---
ax = axes[1]
ax.fill_between(x, 0, h["add"][b_lo:b_hi + 1], step="mid", color="steelblue", alpha=0.35)
ax.plot(x, h["add"][b_lo:b_hi + 1], drawstyle="steps-mid", color="steelblue", lw=1.8,
        label=f"hist_add（{len(h['kicks'])} 次累加）")
ax.plot(x, h["max"][b_lo:b_hi + 1], drawstyle="steps-mid", color="green", lw=1.1, ls="-.",
        label="hist_max（逐 bin 取单次最大）")
ax.plot(x, rr["thresh"][b_lo:b_hi + 1], drawstyle="steps-mid", color="purple", lw=1.8, ls="--",
        label=f"判据门限 = hist_max × {XM_RATIO}")
for q in rr["peaks"]:
    if not (b_lo <= q["peak_bin"] <= b_hi):
        continue
    ax.annotate(f"add/max\n={q['ratio']:.1f}", (q["dist"], q["add"]), fontsize=7,
                ha="center", va="bottom", xytext=(0, 4), textcoords="offset points",
                color=("darkred" if q["is_xtalk"] else "darkgreen"))
    if q["is_xtalk"]:
        ax.plot(q["dist"], q["add"], "rv", ms=10)
ax.axhline(N_ACC, color="k", ls=":", lw=0.9)
ax.text(x[0], N_ACC, f" 真目标应有的高度 = N_ACC = {N_ACC}", fontsize=8, va="bottom")
ax.set_ylim(0, float(h["add"].max()) * 1.45 + 1.0)     # 给峰顶标注留出空间
ax.set_ylabel("计数")
ax.set_title(f"(b) XM 判据：紫色虚线（max×{XM_RATIO}）盖过蓝色 hist_add 的地方 -> 判为串扰（红三角）")
ax.legend(fontsize=8, ncol=3); ax.grid(alpha=0.25)

# --- (c) 滤前 / 滤后 ---
ax = axes[2]
ax.fill_between(x, 0, h["add"][b_lo:b_hi + 1], step="mid", color="0.75", label="XM【前】hist_add")
ax.plot(x, rr["after"][b_lo:b_hi + 1], drawstyle="steps-mid", color="crimson", lw=2.0,
        label="XM【后】")
for q in rr["dropped"]:
    if b_lo <= q["peak_bin"] <= b_hi:
        ax.plot(q["dist"], q["add"], "kx", ms=12, mew=2.2)
for q in rr["peaks"]:
    kind, _ = peak_truth(h, q)
    if not (b_lo <= q["peak_bin"] <= b_hi):
        continue
    col = {"纯真峰": "green", "混合峰": "orange", "纯鬼峰": "red"}[kind]
    ax.axvspan((q["s"] - 0.5) * HIST_BIN_M, (q["e"] + 0.5) * HIST_BIN_M, color=col, alpha=0.13)
ax.axvline(XM_DEMO_D, color="k", ls=":", lw=1.2)
ax.set_ylim(0, float(h["add"].max()) * 1.25 + 0.5)
ax.set_xlabel("记录距离 rec_dist [m]"); ax.set_ylabel("计数")
_nd = len(rr["dropped"])
ax.set_title(f"(c) 结果：灰=滤前，红=滤后，黑叉=被丢弃的峰；"
             f"底色 绿=纯真峰 橙=真假混合峰 红=纯鬼峰"
             + ("　—— 本例一个峰都没丢（红线完全盖住灰色）" if _nd == 0
                else f"　—— 本例丢弃了 {_nd} 个峰"))
ax.legend(fontsize=9); ax.grid(alpha=0.25)

plt.suptitle(f"图A  XM 全链路解剖 —— L{XM_DEMO_LASER}，物体 D={XM_DEMO_D:.0f}m，"
             f"XM_RATIO={XM_RATIO}，N_ACC={N_ACC}，TCODE_MODE='{TCODE_MODE}'",
             fontsize=13, y=0.995)
plt.tight_layout(rect=[0, 0, 1, 0.975])
plt.show()
'''


# ============================================================================
# XM-5：图 B 16 宫格
# ============================================================================
CELL_XM_FIG_GRID = r'''# ============================================================================
# XM-5  图B：全部 16 个激光器的【滤前 / 滤后】对比
# ============================================================================
b_lo, b_hi = occupied_range(hists, pad=15)
x = np.arange(b_lo, b_hi + 1) * HIST_BIN_M

nrow = int(np.ceil(N_LASERS / 4))
fig, axes = plt.subplots(nrow, 4, figsize=(22, 3.1 * nrow), sharex=True)
axes = np.atleast_2d(axes)
ymax = max(h["add"].max() for h in hists.values()) * 1.35 + 0.5

for i, lid in enumerate(laser_ids):
    ax = axes[i // 4][i % 4]
    h, rr = hists[lid], xm_res[lid]
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

    ax.axvline(XM_DEMO_D, color="k", ls=":", lw=1.0, alpha=0.6)
    ax.set_ylim(0, ymax)
    ax.set_title(f"L{lid}：滤掉鬼影{n_hit}、残留鬼影{n_miss}" + (f"、!误杀{n_kill}" if n_kill else ""),
                 fontsize=9)
    ax.tick_params(labelsize=7); ax.grid(alpha=0.2)
    if i == 0:
        ax.legend(fontsize=7, loc="upper left")

for j in range(N_LASERS, nrow * 4):
    axes[j // 4][j % 4].axis("off")

fig.text(0.5, 0.015, f"记录距离 rec_dist [m]（黑点线 = 物体真实距离 {XM_DEMO_D:.0f}m；"
                     f"黑叉 = 被 XM 丢弃的峰；红圈 = ! 误杀的真目标）", ha="center", fontsize=12)
fig.text(0.008, 0.5, "计数（hist_add）", va="center", rotation="vertical", fontsize=12)
plt.suptitle(f"图B  16 个激光器 XM 滤前/滤后对比 —— D={XM_DEMO_D:.0f}m，"
             f"XM_RATIO={XM_RATIO}，TCODE_MODE='{TCODE_MODE}'", fontsize=14, y=0.998)
plt.tight_layout(rect=[0.015, 0.035, 1, 0.965])
plt.show()
'''


# ============================================================================
# XM-6：图 C 统计
# ============================================================================
CELL_XM_FIG_STAT = r'''# ============================================================================
# XM-6  图C：统计可视化 + 不同 tcode / 不同距离的对比表
# ============================================================================
def run_pipeline(D, tcode_mode=None, ratio=None):
    """跑一遍完整流程（换 tcode 时临时替换全局 firings，不改 v13 代码）。"""
    if tcode_mode is None or tcode_mode == TCODE_MODE:
        recs = detect_echoes(D)
    else:
        with use_firings(build_firings_tcode(make_tcode_fn(tcode_mode))):
            recs = detect_echoes(D)
    hs = build_hists(recs)
    rs = crosstalk_mark_all(hs, ratio)
    return hs, rs, evaluate(hs, rs)


#   excel = 用户现状（默认）；lcg2 = docs/tcode 推荐的两层编码，用来预览编码生效的样子
TCODE_COMPARE = ["excel", "lcg2"]
ACTS = ["保留", "丢弃"]

fig = plt.figure(figsize=(18, 10))
gs = fig.add_gridspec(2, 3, hspace=0.40, wspace=0.28)

# --- (a) 记录级：XM 前后 ---
ax = fig.add_subplot(gs[0, 0])
before = [xm_stat["rec_before"][k] for k in ALL_KINDS]
after  = [xm_stat["rec_after"][k]  for k in ALL_KINDS]
xx, w = np.arange(len(ALL_KINDS)), 0.38
ax.bar(xx - w/2, before, w, color="0.7", edgecolor="k", label="XM 前")
ax.bar(xx + w/2, after,  w, color="crimson", edgecolor="k", alpha=0.85, label="XM 后")
for i, (b, a) in enumerate(zip(before, after)):
    ax.text(i - w/2, b, str(b), ha="center", va="bottom", fontsize=8)
    ax.text(i + w/2, a, str(a), ha="center", va="bottom", fontsize=8)
ax.set_xticks(xx); ax.set_xticklabels(ALL_KINDS, fontsize=9)
ax.set_ylabel("回波记录条数")
ax.set_title(f"(a) 记录级 XM 前后（D={XM_DEMO_D:.0f}m, tcode='{TCODE_MODE}'）", fontsize=10)
ax.legend(fontsize=8); ax.grid(alpha=0.25, axis="y")

# --- (b) 峰级混淆矩阵 ---
ax = fig.add_subplot(gs[0, 1])
KINDS = ["纯真峰", "混合峰", "纯鬼峰"]
M = np.array([[xm_stat["peak"][(k, a)] for a in ACTS] for k in KINDS], dtype=float)
im = ax.imshow(M, cmap="YlOrRd", aspect="auto")
for i in range(M.shape[0]):
    for j in range(M.shape[1]):
        ax.text(j, i, f"{int(M[i,j])}", ha="center", va="center", fontsize=14,
                color=("white" if M[i, j] > M.max() * 0.6 else "black"))
ax.set_xticks([0, 1]); ax.set_xticklabels(ACTS)
ax.set_yticks(range(3)); ax.set_yticklabels(KINDS)
ax.set_title("(b) 峰级混淆矩阵\n(纯鬼峰,丢弃)=命中   (纯鬼峰,保留)=漏滤\n"
             "(纯真峰/混合峰,丢弃)=! 误杀", fontsize=9)
fig.colorbar(im, ax=ax, fraction=0.045)

# --- (c) 逐激光器 ---
ax = fig.add_subplot(gs[0, 2])
hit  = [xm_stat["per_laser"][l][("纯鬼峰", "丢弃")] for l in laser_ids]
miss = [xm_stat["per_laser"][l][("纯鬼峰", "保留")] for l in laser_ids]
kill = [sum(xm_stat["per_laser"][l][(k, "丢弃")] for k in ("纯真峰", "混合峰")) for l in laser_ids]
ax.bar(laser_ids, hit, color="seagreen", edgecolor="k", label="鬼影峰·正确滤除")
ax.bar(laser_ids, miss, bottom=hit, color="salmon", edgecolor="k", label="鬼影峰·残留")
ax.bar(laser_ids, kill, bottom=np.array(hit) + np.array(miss), color="red", edgecolor="k",
       hatch="//", label="! 真目标被误杀")
ax.set_xlabel("接收激光器"); ax.set_ylabel("峰数")
ax.set_xticks(laser_ids); ax.tick_params(labelsize=7)
ax.set_title("(c) 逐激光器鬼影峰去留", fontsize=10)
ax.legend(fontsize=7.5); ax.grid(alpha=0.25, axis="y")

# --- (d)(e) 两种 tcode × 四个距离 ---
for col, mode in enumerate(TCODE_COMPARE):
    ax = fig.add_subplot(gs[1, col])
    gb, ga, tk = [], [], []
    for D in D_TEST:                      # D_TEST 来自 v13 的验证 cell
        _, _, s = run_pipeline(D, mode)
        gb.append(sum(s["peak"][("纯鬼峰", a)] for a in ACTS))
        ga.append(s["peak"][("纯鬼峰", "保留")])
        tk.append(sum(s["peak"][(k, "丢弃")] for k in ("纯真峰", "混合峰")))
    xx = np.arange(len(D_TEST))
    ax.bar(xx - 0.21, gb, 0.42, color="0.7", edgecolor="k", label="鬼影峰·XM 前")
    ax.bar(xx + 0.21, ga, 0.42, color="crimson", edgecolor="k", label="鬼影峰·XM 后残留")
    ax.plot(xx + 0.21, tk, "k^--", ms=8, label="! 误杀真目标峰")
    for i in range(len(D_TEST)):
        ax.text(xx[i] - 0.21, gb[i], str(gb[i]), ha="center", va="bottom", fontsize=8)
        ax.text(xx[i] + 0.21, ga[i], str(ga[i]), ha="center", va="bottom", fontsize=8)
    ax.set_xticks(xx); ax.set_xticklabels([f"{d:.0f}m" for d in D_TEST])
    ax.set_ylabel("鬼影峰数（16 激光器合计）")
    ax.set_title(f"({'de'[col]}) TCODE_MODE='{mode}'"
                 + ("  ← 用户现状：码不随 kick 变 -> XM 无效"
                    if mode == "excel" else "  ← 预览：码随 kick 变 -> XM 生效"), fontsize=9.5)
    ax.legend(fontsize=8); ax.grid(alpha=0.25, axis="y")

# --- (f) 文字说明 ---
ax = fig.add_subplot(gs[1, 2]); ax.axis("off")
ax.text(0.0, 1.0,
        "为什么 'excel' 这一栏几乎滤不掉？\n"
        "──────────────────────────\n"
        "串扰落点：\n"
        "  rec_tof = 2D/c + (tx_发 − tx_收)\n"
        "                 + (k_发 − k_收)·T_kick\n"
        "真回波：\n"
        "  rec_tof = 2D/c      （与编码无关）\n\n"
        "Excel 里每个激光器 4 个 kick 的 tx 是同一个值\n"
        "  => (tx_发 − tx_收) 每次都一样\n"
        "  => 同一条串扰 4 次落到同一个 bin\n"
        f"  => hist_add 也累到 {N_ACC}，add/max = {N_ACC}\n"
        f"  => {N_ACC} > XM_RATIO = {XM_RATIO}，XM 判它是真目标\n\n"
        "tcode 的任务因此非常明确：\n"
        f"  让每一对(发,收)的码差在 {N_ACC} 个 kick 里互不相同，\n"
        "  且差值间隔 ≥ 峰宽，真正落到不同 bin。\n"
        "  -> 原理与图解见 docs/tcode/tcode_scheme.md",
        fontsize=9.5, va="top", linespacing=1.5)

plt.suptitle(f"图C  XM 统计总览（XM_RATIO={XM_RATIO}，N_ACC={N_ACC}；"
             f"裕度区间 = (串扰重复次数, {N_ACC})）", fontsize=14, y=0.985)
plt.show()

# ---- 文字统计表 ----
print("\n" + "=" * 96)
print("多距离 × 多 tcode 统计表")
print("=" * 96)
print(f"  {'tcode':>8} {'D[m]':>6} {'真峰':>5} {'混峰':>5} {'鬼峰':>5} | "
      f"{'鬼峰丢弃':>8} {'鬼峰残留':>8} {'滤除率':>7} | {'误杀':>5} {'误杀率':>7}")
print("  " + "-" * 92)
for mode in TCODE_COMPARE:
    for D in D_TEST:
        _, _, s = run_pipeline(D, mode)
        P = s["peak"]
        n_pure = sum(P[("纯真峰", a)] for a in ACTS)
        n_mix  = sum(P[("混合峰", a)] for a in ACTS)
        n_gh   = sum(P[("纯鬼峰", a)] for a in ACTS)
        hit, miss = P[("纯鬼峰", "丢弃")], P[("纯鬼峰", "保留")]
        kill = sum(P[(k, "丢弃")] for k in ("纯真峰", "混合峰"))
        print(f"  {mode:>8} {D:>6.0f} {n_pure:>5d} {n_mix:>5d} {n_gh:>5d} | "
              f"{hit:>8d} {miss:>8d} {hit/max(n_gh,1):>6.1%} | "
              f"{kill:>5d} {kill/max(n_pure+n_mix,1):>6.1%}")
'''


# ============================================================================
# XM-7：图 D 距离扫描
# ============================================================================
CELL_XM_SWEEP = r'''# ============================================================================
# XM-7  图D：距离扫描 —— 物体从近到远走一遍，看 XM 前后鬼影残留 & 真目标存活
#   这是最终目标（噪点率 10ppm）的雏形指标。
# ============================================================================
SWEEP_LO, SWEEP_HI, SWEEP_STEP = 5.0, 600.0, 5.0     # 600m = 2 × D_UNAMBIG
D_SWEEP = np.arange(SWEEP_LO, SWEEP_HI + 1e-9, SWEEP_STEP)


def sweep(tcode_mode, progress_every=30):
    out = {k: [] for k in ("ghost_before", "ghost_after", "true_before", "true_killed")}
    for i, D in enumerate(D_SWEEP):
        _, _, s = run_pipeline(D, tcode_mode)
        P = s["peak"]
        out["ghost_before"].append(sum(P[("纯鬼峰", a)] for a in ACTS))
        out["ghost_after"].append(P[("纯鬼峰", "保留")])
        out["true_before"].append(sum(P[(k, a)] for k in ("纯真峰", "混合峰") for a in ACTS))
        out["true_killed"].append(sum(P[(k, "丢弃")] for k in ("纯真峰", "混合峰")))
        if progress_every and i % progress_every == 0:
            print(f"    [{tcode_mode}] {i+1}/{len(D_SWEEP)}  D={D:.0f}m ...")
    return {k: np.asarray(v) for k, v in out.items()}


print(f"距离扫描中（{len(D_SWEEP)} 个距离 × {len(TCODE_COMPARE)} 种 tcode，稍等）...")
SW = {m: sweep(m) for m in TCODE_COMPARE}
print("扫描完成。")

fig, axes = plt.subplots(3, 1, figsize=(16, 12), sharex=True)
COLS = {"excel": "tab:blue", "random": "tab:orange", "none": "tab:green"}

ax = axes[0]
ax.plot(D_SWEEP, SW[TCODE_COMPARE[0]]["ghost_before"], "k--", lw=1.4,
        label="鬼影峰·XM 前")
for m in TCODE_COMPARE:
    ax.plot(D_SWEEP, SW[m]["ghost_after"], "-", color=COLS.get(m), lw=1.8,
            label=f"鬼影峰·XM 后残留（tcode='{m}'）")
ax.axvline(D_UNAMBIG, color="r", ls=":", lw=1.2)
ax.text(D_UNAMBIG, ax.get_ylim()[1]*0.9, f" D_UNAMBIG={D_UNAMBIG:.0f}m", color="r", fontsize=9)
ax.set_ylabel("鬼影峰数（16 激光器合计）")
ax.set_title("(a) 鬼影峰：XM 滤除前 vs 滤除后")
ax.legend(fontsize=9); ax.grid(alpha=0.25)

ax = axes[1]
for m in TCODE_COMPARE:
    resid = SW[m]["ghost_after"] / np.maximum(SW[m]["ghost_before"], 1)
    ax.plot(D_SWEEP, resid * 100, "-", color=COLS.get(m), lw=1.8, label=f"残留率 tcode='{m}'")
ax.axhline(0, color="k", lw=0.8)
ax.set_ylabel("鬼影残留率 [%]")
ax.set_title("(b) 鬼影残留率（越低越好；最终目标是压到噪点率 10ppm 量级）")
ax.legend(fontsize=9); ax.grid(alpha=0.25)

ax = axes[2]
ax.plot(D_SWEEP, SW[TCODE_COMPARE[0]]["true_before"], "k--", lw=1.4, label="真目标峰·XM 前")
for m in TCODE_COMPARE:
    ax.plot(D_SWEEP, SW[m]["true_before"] - SW[m]["true_killed"], "-", color=COLS.get(m),
            lw=1.8, label=f"真目标峰·XM 后存活（tcode='{m}'）")
    ax.plot(D_SWEEP, SW[m]["true_killed"], "--", color=COLS.get(m), lw=1.1, alpha=0.7,
            marker="x", ms=4, label=f"! 误杀（tcode='{m}'）")
ax.axvline(D_UNAMBIG, color="r", ls=":", lw=1.2)
ax.set_xlabel("物体真实距离 D [m]"); ax.set_ylabel("真目标峰数")
ax.set_title("(c) 真目标：XM 后还剩几个（D>300m 时真回波本来就超窗丢失，与 XM 无关）")
ax.legend(fontsize=8, ncol=2); ax.grid(alpha=0.25)

plt.suptitle(f"图D  距离扫描 {SWEEP_LO:.0f}~{SWEEP_HI:.0f}m（步长 {SWEEP_STEP:.0f}m）"
             f" —— XM 滤前/滤后，XM_RATIO={XM_RATIO}", fontsize=14, y=0.995)
plt.tight_layout(rect=[0, 0, 1, 0.97])
plt.show()

print("\n" + "=" * 88)
print(f"距离扫描汇总（{len(D_SWEEP)} 个距离点 × 16 个激光器）")
print("=" * 88)
print(f"  {'tcode':>8} {'鬼影峰(前)':>11} {'鬼影峰(后)':>11} {'滤除率':>8} "
      f"{'真峰(前)':>9} {'误杀':>6} {'误杀率':>8}")
print("  " + "-" * 84)
for m in TCODE_COMPARE:
    gb, ga = SW[m]["ghost_before"].sum(), SW[m]["ghost_after"].sum()
    tb, tk = SW[m]["true_before"].sum(), SW[m]["true_killed"].sum()
    print(f"  {m:>8} {gb:>11d} {ga:>11d} {(gb-ga)/max(gb,1):>7.2%} "
          f"{tb:>9d} {tk:>6d} {tk/max(tb,1):>7.2%}")
'''


# ============================================================================
# XM-8：总结
# ============================================================================
CELL_XM_SUMMARY = r'''# ============================================================================
# XM-8  总结
# ============================================================================
print("=" * 88)
print("crosstalk_sim_v20 总结")
print("=" * 88)
print(f"""
0. 本 notebook 的结构
   - cell 1~10 = v13 原封不动（kick 栅格、记录距离分布、逐回波堆叠柱，一个没删）
   - 之后 = v20 新增的 XM 部分

1. 实装的 XM（只有用户要求的那一条判据）
   - 每个激光器 {N_ACC} 次 shot 全部留存 → hist_add（累加）、hist_max（逐 bin 取最大）
   - 判据：hist_max × {XM_RATIO} > hist_add  →  该峰是串扰，丢弃
   - 等价：add/max < {XM_RATIO}；而 add/max = 该峰在几次 shot 里出现过
   - 硬件文档里的分档 / flag0~3 / diff_flag / FIR / max 抽取，本版一律未实装

2. 阈值裕度
   - 必须满足： 串扰最大重复次数 < XM_RATIO({XM_RATIO}) < N_ACC({N_ACC})
   - 现在的 {XM_RATIO} 只能滤掉"{N_ACC} 次里只出现 1 次"的串扰
   - 想连"出现 2 次"的也滤掉 → XM_RATIO 抬到 2.5；但真目标保护裕度从 2.5× 降到 1.6×

3. 本版结论
   - TCODE_MODE='excel'（用户现状）：码不随 kick 变 → 串扰每次落同一 bin
     → add/max = {N_ACC} → XM 滤除率 ≈ 0（图C(d)、图D 已证实）
   - TCODE_MODE='lcg2'（预览，docs/tcode 推荐的两层编码）：码随 kick 变
     → 串扰被打散 → XM 立刻生效（图C(e)、图D）
   - 说明 XM 实现无误，**瓶颈在编码不在算法**

4. 下一版 v21：tcode（各激光器之间的 tx_trig_dly 编码）
   - 目标：对任意一对激光器 (a,b)，码差 c[a][k] − c[b][k] 在 {N_ACC} 个 kick 上互不相同，
     且两两间隔 ≥ 峰宽 → 同一条串扰散到 {N_ACC} 个不同 bin → add/max = 1 → 被 XM 全滤
   - 推荐构造（docs/tcode 已验证）：
       c[l][k] = ((l*(k+1)) mod 17)*8ns  +  ((k*k) mod 9)*8ns
       第一项治同 kick 串扰；第二项对所有激光器相同，治自身混叠
   - 已验证成绩（1~600m 扫描，误杀 0）：
       鬼影残留率  XM_RATIO=1.6 → 5.6%    XM_RATIO=2.5 → 0.53%
   - 构造推导、约束、图解、验收指标：见 docs/tcode/tcode_scheme.md

5. 以后可选的加强项（本版只提不做）
   - 次大值方案 add − max − max2：裕度从 (1,N) 扩到 (2,N)，能吃掉重复 2 次的串扰
   - 软减法 out = max(add − α·max, 0)：真目标峰上叠了串扰时不会连真目标一起丢
   - n_hit 逐 bin 命中计数：N=4 只要 2bit/bin，比 hist_max 更省，但丢幅度信息
""")
'''


# ============================================================================
# 组装：v13 全部 cell 原样继承 + 前后追加
# ============================================================================
def code_cell(cid, src):
    return {"cell_type": "code", "id": cid, "metadata": {},
            "execution_count": None, "outputs": [],
            "source": src.splitlines(keepends=True)}


def md_cell(cid, src):
    return {"cell_type": "markdown", "id": cid, "metadata": {},
            "source": src.splitlines(keepends=True)}


if not os.path.exists(SRC_NB):
    raise SystemExit(f"找不到 {SRC_NB}，无法继承 v13 的 cell。")

with open(SRC_NB, encoding="utf-8") as f:
    nb13 = json.load(f)

# v13 的 cell 一字不改地搬过来（连 outputs 一起清掉，保持干净的未执行状态）
v13_cells = []
for c in nb13["cells"]:
    c = dict(c)
    if c["cell_type"] == "code":
        c["outputs"] = []
        c["execution_count"] = None
    v13_cells.append(c)

cells = (
    [md_cell("v20_overview", CELL_V20_OVERVIEW)]
    + v13_cells                                          # ← v13 原样继承，不删不改
    + [
        md_cell("v20_xm_doc",       CELL_XM_DOC),
        code_cell("v20_xm_hist",    CELL_XM_HIST),
        code_cell("v20_xm_core",    CELL_XM_CORE),
        code_cell("v20_xm_stat",    CELL_XM_STAT),
        code_cell("v20_xm_chain",   CELL_XM_FIG_CHAIN),
        code_cell("v20_xm_grid",    CELL_XM_FIG_GRID),
        code_cell("v20_xm_figstat", CELL_XM_FIG_STAT),
        code_cell("v20_xm_sweep",   CELL_XM_SWEEP),
        code_cell("v20_xm_summary", CELL_XM_SUMMARY),
    ]
)

nb = {
    "cells": cells,
    "metadata": nb13.get("metadata", {}),
    "nbformat": 4,
    "nbformat_minor": 5,
}

with open(OUT_NB, "w", encoding="utf-8") as f:
    json.dump(nb, f, ensure_ascii=False, indent=1)

print(f"已生成 {OUT_NB}")
print(f"  继承 v13 的 cell：{len(v13_cells)} 个（原样，未删未改）")
print(f"  v20 新增 cell   ：{len(cells) - len(v13_cells)} 个")
print(f"  合计            ：{len(cells)} 个")
