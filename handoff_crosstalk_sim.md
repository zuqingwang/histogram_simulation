# 交接文档 —— 当前工作：`crosstalk_sim`（串扰仿真 + tcode 编码）

> 本文档文件名：`handoff_crosstalk_sim.md`（工作名 = `crosstalk_sim`；禁止 `handoff_现在工作.md`）。
> **本线包含 tcode（发光时刻编码）工具链**（`tcode_calculator*.ipynb`、`docs/tcode/`），不另拆工作名，除非用户明确要求。
> 写给**完全没有上下文的新会话**。文件名、模块号、变量名、参数值均写全。
> 最后更新：2026-08-04。
> 配套流水日志：`worklog_crosstalk_sim.md`。
> 兄弟工作线：`handoff_lidar_histogram_sim.md`、`handoff_PoD_esti.md`。

---

## 0. 缩写表

| 缩写 | 英文全称 | 含义 |
|---|---|---|
| XM | XtalkMark | 串扰标记滤除：挑「只在少数发光里出现」的峰 |
| tcode | timing code / 发光时刻编码 | 让各 kick 的发射时刻码变化，打散鬼影落点 |
| ToF | Time of Flight | 飞行时间 |
| SPAD | Single-Photon Avalanche Diode | 单光子雪崩二极管 |
| IRF | Instrument Response Function | 仪器响应函数 |
| FPGA | Field-Programmable Gate Array | 现场可编程门阵列（本线用其模拟累计抖动） |
| LCG | Linear Congruential | 线性同余（早期 tcode 文档构造；现生产多用离散字母表） |
| SEP | separation | 最小码差（本项目默认 12 ns） |
| MC | Monte Carlo | 蒙特卡洛 |

---

## 1. 我们在做什么任务

仓库：`E:\claude temp\Histogram-simulation`。

工作线 **B：模组内串扰 + 雷达对射 + 编码滤噪**，主产物 **`crosstalk_sim_v42.ipynb`**。

回答：
1. Excel 固定码为什么滤不掉鬼影；
2. tcode（离散字母表）如何把鬼影 `hist_add/hist_max` 打散到可被 XM 滤掉；
3. 对射同型号同 tcode 时，如何靠双方 FPGA / 角度抖动 / 一字滤波清掉；
4. 在给定字母表与 (N,M,L) 下，编码可搜索性与最优配置。

**边界**：本线用理想 δ 回波 + XM；能量展宽/SPAD 物理内核在 `lidar_histogram_sim`；FAR/PoD 在 `PoD_esti`。**禁止**为 PoD 任务改本线参数，也禁止本线去改 lidar/PoD 物理值。

---

## 2. 已经完成了什么

### 2.1 版本地图（摘要）

| 阶段 | 文件 | 要点 |
|---|---|---|
| 早期 | `crosstalk_sim_v01`–`v03` | 串扰起步 |
| 长焦 A 组 | v10–v13 | 16 激光器 kick 栅格、δ 回波、鬼影检测 |
| XM + tcode | v20–v21 | XM 滤噪；tcode 实装 |
| 对射 / FPGA | v220、v22–v24 | 双 tcode、FPGA、对射 MC、一字滤波教学 |
| 干净重写 | v30 | 参数字典 + 四步滤噪（i→iv） |
| 离散字母表 | v40 | tcode 从连续整数 → 有限字母表（100 ns 预算） |
| 可变 kick 间隔 | v41 | kick 仅触发；间隔 = max(enc)+2.2 μs |
| **当前** | **`crosstalk_sim_v42.ipynb`** | 画图 100 ns 压缩、角度抖动、consider_gap、(N,M,L) 蒙卡 |

构造脚本：`build_crosstalk_v*.py` 到 **v41** 止；**v42 无 builder（手工演进）**。

### 2.2 tcode 工具链（本线附属）

| 文件 | 作用 |
|---|---|
| `tcode_calculator.ipynb` | v1：连续整数搜索；已知零残留表 |
| `tcode_calculator_v2.ipynb` | v2：**默认 v40 离散字母表**；min-conflicts |
| `build_tcode_calculator.py` / `build_tcode_calculator_v2.py` | 生成器 |
| `docs/tcode/tcode_scheme.md` | 原理（偏 v20/v21 LCG；与现生产有代差） |
| `docs/tcode/solve_tcode.py` | 码表求解主脚本 |
| `docs/tcode/tcode_table*.py` | 各预算/ratio 零残留表 |

**当前仿真采用的码表**（v42 cell 05 内嵌加载）：
- `docs/tcode/tcode_table_v40_r2.5_L5_100ns.py` → `{0,25,50,75,100}`
- `docs/tcode/tcode_table_v40_r1.5_L17_100ns.py` → 17 档 linspace 0…100

验收口径（`tcode_calculator_v2`）：**1~600 m 扫描，XM 后鬼影残留=0 且真峰误杀=0**。

### 2.3 v42 模块结构（25 cell）

| Cell | 内容 |
|---|---|
| 01 | 全局参数字典 `TIMING`/`SCENE`/`XM`/`TCODE`/`FPGA`/`ANGLE_JITTER`/`RADAR`/`LINE` |
| 03–04 | Excel 发光时序；**kick 栅格** `plot_kick_grid` |
| 05 | tcode 码表 + kick×激光器矩阵 |
| 06 | 核心仿真与绘图（可变 kick、XM、对射、堆叠柱） |
| 07 | 一字滤波 `apply_line_filter` |
| 08–18 | 四步演示 i→iv + 对比图 |
| 19 | 总结（**文案仍写 v41，待改**） |
| 20–21 | consider_gap True/False 可搜索性（已跑） |
| 22–23 | **(N,M,L) 蒙卡优化**（**cell 23 未执行**） |

### 2.4 四步滤噪逻辑（已跑通演示，D=150 m）

```
i   Excel 固定码     → 鬼影+对射均难滤
ii  +tcode           → 模组鬼影应打散（单距演示仍可能有残留峰）
iii +双方 FPGA 抖动  → 对射可被 XM 丢掉
iv  +一字滤波        → 三角度清对射孤点
```

### 2.5 当前关键参数（v42 cell 01）

```
TIMING.tof_fixed_us = 2.2
TIMING.tof_window_ns = 2000.0      # 计算窗
TIMING.plot_tof_ns = 100.0         # 仅画图压缩；计算仍 2000 ns
SCENE.demo_D_m = 150.0；demo_laser = 5
XM.ratio = 2.5；hist_bin_ns = 1.0；use_pulse_width = False
TCODE.ratio_mode = "2.5"；consider_gap = True；max_gap = 2
TCODE.sep_ns = 12；budget_ns = 100
TCODE.alphabet_r25 = [0,25,50,75,100]
FPGA.step_ns = 8，n_levels = 8，global_delay_ns = 8
ANGLE_JITTER.enable = True，step_ns = 8，n_levels = 8
RADAR.enable = True，phase_ns = 700.0，同型号同 tcode
LINE.thr_m = 3.0（演示 iv 启用；字典默认 enable=False）
N_MC = 10000（cell 23）；N∈{3,5,7,9}，M∈{1,4,8}，L∈{1,4,8}
```

### 2.6 绘图规范锚点（全项目统一）

| 图类 | 参考 / v42 落点 |
|---|---|
| kick 栅格 | 风格承自 v23；v42 **cell 04** `plot_kick_grid` |
| tcode 矩阵 | 承自 v220；v42 **cell 05** imshow |
| 回波堆叠柱 | 承自 v23；v42 cell 06：绿真 / 橙鬼 / **红+斜线对射** |
| 滤除标记 | 黑叉=滤除；规范还要求红空心圈=误杀（**v42 尚未画红圈**） |

---

## 3. 目前面临的问题、卡在哪里

1. **`crosstalk_sim_v42.ipynb` cell 23 (N,M,L) 蒙卡未跑** → Pareto / 加权最优尚无结果。
2. **无 `build_crosstalk_v42.py`** → v42 相对 v41 的差异未脚本化，难复现生成。
3. **cell 19 总结仍打印 v41**。
4. **consider_gap 随机搜索**（cell 21）：30 次试验零噪点均为 0/30；最优噪点率约 0.15–0.24，**不是** tcode 计算器「1~600 m 零残留」口径。
5. **单距 D=150 m 演示** 与 **全距离扫描零残留验收** 不要混谈。
6. **`tcode_scheme.md` 仍以 LCG 为主**；生产以 v40 离散字母表 + `solve_tcode.py` 为准（文档代差）。
7. v42 **未实现红空心圈误杀标记**（与仓库绘图规范尚有缺口）。

---

## 4. 下一步计划

1. 跑通 cell 23（可先降 `N_MC` smoke test）；结果写入 worklog。
2. 补 `build_crosstalk_v42.py`，固化 v41→v42 差异。
3. 修正 cell 19 总结为 v42；补红空心圈误杀标记。
4. 将 consider_gap 随机搜索与 `tcode_calculator_v2` / `solve_tcode.py` 正式求解对接。
5. 可选：更新 `docs/tcode/tcode_scheme.md`，增补 v40 离散字母表章节。
6. 长 MC 按规则三加缓存 + 多线程（目前串扰线尚未完全按 PoD v05 方式固化）。

---

## 5. 踩过的坑（不要再踩）

1. **Excel 固定码** → 鬼影 add/max=4 → XM 滤不掉；必须逐 kick 变码。
2. **码步长 < 峰宽** → 码不同仍落同一宽峰；SEP≥12 ns；离散字母表拉大档距（v40 动机）。
3. **对射只靠 tcode** → 同型号同 tcode 对射稳留；需双方 FPGA + 可选角度抖动/一字滤波。
4. **对射标识** → 不能只靠竖虚线；必须 **红色 + 斜线纹理**，并与真回波同图画、标是否同峰宽。
5. **真/鬼/对射** → 绿 / 橙 / 红斜线；柱内标激光器+kick。
6. **滤除 vs 误杀** → 黑叉=滤除；红空心圈=误杀真峰。
7. **画图 100 ns** → 仅示意压缩；**计算窗仍 2000 ns**；勿把 `plot_tof_ns` 当物理窗。
8. **单距演示 ≠ 全扫描验收**。
9. **工作线耦合** → 禁止改 PoD/lidar 时顺手改串扰，反之亦然。
10. **v220 不继承 v22**：对射/FPGA 主线从 v21→v220；v22–v24 是另一支教学链。

---

## 6. 给新会话的最短指令

1. 打开 `crosstalk_sim_v42.ipynb` + 本文件 + `worklog_crosstalk_sim.md`。
2. 码表求解看 `tcode_calculator_v2.ipynb` 与 `docs/tcode/solve_tcode.py`。
3. 先确认 cell 23 是否已跑；未跑则按用户要求启动或降采样试跑。
4. 改 SPAD/能量物理 → `lidar_histogram_sim`；改 FAR/PoD → `PoD_esti`。
5. 会话结束前更新本文件与 worklog。
