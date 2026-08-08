# worklog —— `crosstalk_sim`（串扰仿真 + tcode 编码）

> 顶部 = 现状与当前任务（可覆盖更新）；下部 = 历史记录（只追加，不删改）。
> 交接文档见 `handoff_crosstalk_sim.md`。
> **tcode_calculator / `docs/tcode/` 归入本工作名**，不另开 `worklog_tcode`，除非用户另行指定。
> 本日志于 2026-08-04 补建；早期版本条目为根据 notebook / build 脚本回溯整理。

---

# 一、现状 / 当前任务

**当前版本：v42（2026-08-04）**
**主产物：`crosstalk_sim_v42.ipynb`（25 cell）**
**tcode 主工具：`tcode_calculator_v2.ipynb` + `docs/tcode/solve_tcode.py` + v40 离散字母表**

## 1. 目标一句话

在长焦 A 组多激光器时序下，用 **tcode（发光时刻编码）+ XM（XtalkMark，串扰标记）+ FPGA/角度抖动 + 一字滤波**，滤除模组内串扰鬼影与雷达对射，并搜索可行的离散字母表与 (N,M,L) 配置。

## 2. 当前可运行状态

- **`crosstalk_sim_v42.ipynb`**：cell 0–22 已有输出；**cell 23 (N,M,L) 蒙卡未执行**。
- 四步滤噪演示（D=150 m）已跑通；consider_gap 可搜索性试验（cell 21）已跑。
- **无 `build_crosstalk_v42.py`**（v42 手工演进）。
- tcode：生产码表为 v40 离散字母表（ratio 2.5 → L5；ratio 1.5 → L17；预算 100 ns）。
- 兄弟线：`lidar_histogram_sim_v45`（物理内核）、`PoD_esti_v05`（FAR/PoD）——**勿交叉改参**。

## 3. 当前生效关键参数（v42）

```
tof_window_ns = 2000（计算）；plot_tof_ns = 100（仅画图）
demo_D_m = 150；XM.ratio = 2.5；SEP = 12 ns；budget = 100 ns
TCODE.ratio_mode = "2.5"；consider_gap = True；max_gap = 2
alphabet_r25 = [0,25,50,75,100]
ANGLE_JITTER.enable = True；RADAR.enable = True
N_MC = 10000（cell 23 待跑）
```

## 4. 待办

- [ ] 执行 cell 23 (N,M,L) 蒙卡（或先降 N_MC smoke）
- [ ] 补 `build_crosstalk_v42.py`
- [ ] cell 19 总结改成 v42；补红空心圈误杀标记
- [ ] consider_gap 搜索对接 `solve_tcode.py` / `tcode_calculator_v2`
- [ ] 可选：更新 `docs/tcode/tcode_scheme.md` 补离散字母表章节

---

# 二、历史记录（只追加，不删改）

## 文档补建 —— 2026-08-04

### 新增
- `handoff_crosstalk_sim.md`、`worklog_crosstalk_sim.md`（本文件）。
- 明确 tcode 工具链归属本工作名。

---

## v42 —— 画图压缩 / 角度抖动 / consider_gap / (N,M,L)

### 新增
- 时序链按 **100 ns** 压缩排布（计算仍 2000 ns）。
- **角度抖动** `ANGLE_JITTER`（默认 enable）。
- consider_gap True/False 编码可搜索性试验（cell 20–21，已跑）。
- (N,M,L) 蒙卡优化框架（cell 22–23；**23 未跑**）。

### 删减 / 待修
- 无独立 builder；cell 19 总结仍写 v41。
- 红空心圈误杀标记尚未实现。

### 运行结论（已有）
- 四步演示 D=150 m、XM.ratio=2.5、consider_gap=True：iii 步对射可清、iv 步一字滤波后真留完整、鬼留 0。
- consider_gap 随机 30 次：零噪点 0/30；最优噪点率约 0.15–0.24（非全扫描零残留口径）。

### 关键参数（相对 v41）
| 项 | 说明 |
|---|---|
| `plot_tof_ns` | 新增画图压缩窗 100 ns |
| `ANGLE_JITTER` | 新增 |
| `TCODE.opt_step_ns` | NML 搜索用 16 ns 档 |

---

## v41 —— 可变 kick 间隔

### 新增
- kick 仅作触发；间隔 = max(enc) + 2.2 μs（`tof_fixed_us`）。
- `build_crosstalk_v41.py`。

---

## v40 —— 离散字母表 tcode

### 新增
- tcode 从连续整数改为有限字母表；预算 100 ns。
- 码表：`tcode_table_v40_r2.5_L5_100ns.py`、`tcode_table_v40_r1.5_L17_100ns.py`。
- `tcode_calculator_v2.ipynb`：默认离散字母表 + min-conflicts。
- 动机：码步长必须 ≥ 峰宽（SEP），连续细码会被 IRF 展宽「粘」回同一峰。

### 运行结论
- ratio=2.5、L5：1~600 m 扫描零残留（计算器验收口径）。
- ratio=1.5、L17：零残留，档数更多。

---

## v30 —— 干净重写四步滤噪

### 新增
- 参数字典化；四步：i Excel → ii tcode → iii 双方 FPGA → iv 一字滤波。
- `build_crosstalk_v30.py`。

---

## v220 / v22–v24 —— 对射与 FPGA 支线

### 新增
- **v220**（自 v21，不继承 v22）：双 tcode、累计 FPGA、对射 MC；tcode 矩阵图风格锚点。
- **v23**：对射波形教学链；kick 栅格 / 堆叠柱 / 滤除黑叉 的绘图风格锚点。
- **v24**：一字滤波（三角度）。

### 踩坑
- 对射同型号同 tcode → 仅靠编码不够。
- 对射必须红+斜线，不能只靠虚线。

---

## v20–v21 —— XM + tcode 实装

### 新增
- v20：XM（XtalkMark）。
- v21：tcode 实装；早期码表 `docs/tcode/tcode_table.py`；`tcode_scheme.md`（LCG 两层公式）。
- `tcode_calculator.ipynb`（连续整数搜索）。

### 关键结论
- Excel 固定码 → 鬼影与真目标 add/max 同为 4 → 滤不掉。
- 变码后鬼影残留可大幅下降（方案文档记载量级：97.5% → 数 % 或更低，视 xm_ratio）。

---

## v10–v13 —— 长焦 A 组基础仿真

### 新增
- 16 激光器、Excel 时序、kick 栅格、δ 回波、鬼影检测框架。

---

## v01–v03 —— 起步

### 新增
- 串扰仿真最初版本（含 `*_debugged` 副本）。
