# worklog —— `lidar_histogram_sim`（激光雷达直方图物理仿真）

> 顶部 = 现状与当前任务（可覆盖更新）；下部 = 历史记录（只追加，不删改）。
> 交接文档见 `handoff_lidar_histogram_sim.md`。
> 本日志于 2026-08-04 补建；早期版本条目为根据 notebook / build 脚本回溯整理，非当时逐日流水。

---

# 一、现状 / 当前任务

**当前版本：v45（2026-08-04）**
**主产物：`lidar_histogram_sim_v45.ipynb`（约 60 cell）**

## 1. 目标一句话

完整光链路 + SPAD（Single-Photon Avalanche Diode，单光子雪崩二极管）RC 恢复 + **1 ns HDC（Hardware Data Channel，硬件数据通道）二值采样**，产出宏像元直方图、前沿测距、SNR（Signal-to-Noise Ratio，信噪比）与能量/距离扫描分析；并为 `PoD_esti` 提供可提取的物理内核。

## 2. 当前可运行状态

- **`lidar_histogram_sim_v45.ipynb`**：当前主产物；模块 19/20 已有输出；体积约 3.6 MB（含大量图）。
- 构造链：`build_v*_from_*.py`、`patch_v44.py`；v45 在 v44 上手工增强模块 19。
- **停用模块（v43 起）**：8、8b、9、17、15、16（100 ppm）。
- 兄弟线：`PoD_esti_v05.ipynb`（FAR/PoD）、`crosstalk_sim_v42.ipynb`（串扰+tcode）——**勿交叉改本文件物理参数**。

## 3. 当前生效关键参数

```
λ=905 nm；P_peak=235 W；τ_r=0.7 ns，τ_f=1.9 ns
D=30 m，ρ=0.10；E_lambda=0.68 W/m²/nm
PDE=0.30，DCR=0，IRF σ=100 ps，τ_RC=8.7315 ns，Vth_frac=0.60
bin=1 ns，N_shots=4；宏像元 9×3，macro_cap=108
仿真窗：ToF−50~+100 ns；护带 ≈8.5 ns
dt_fine 实际：模块 0b = 200 ps（PARAMS 字典仍写 10 ps）
```

## 4. 待办

- [ ] 若重跑模块 19/20 长 MC：按规则三加缓存 + 多线程
- [ ] 物理参数若变更：同步评估是否影响 `PoD_esti`
- [ ] 可选：抽离模块 19/20 为轻量脚本，减小打开成本

---

# 二、历史记录（只追加，不删改）

## 文档补建 —— 2026-08-04

### 新增
- `handoff_lidar_histogram_sim.md`、`worklog_lidar_histogram_sim.md`（本文件）。
- 与 `PoD_esti` / `crosstalk_sim` 交接文档交叉引用。

### 说明
- 以下版本条目为回溯整理，细节以各版 notebook 内 markdown 为准。

---

## v45 —— 模块 19 增强（能量扫描诊断）

### 新增 / 改动
- 模块 19：环境底噪诊断；指数分段 boost 网格；扣背景重心。
- 基于 `lidar_histogram_sim_v44.ipynb` 复制手改（无独立 `build_v45`）。

### 运行结论 / 现象
- 线性 boost 加密 ≠ peak/area 均匀采样。
- 二值 8 ns 过阈窗导致 walk 形态与模拟前端 RC 曲线不一致。

---

## v44 —— `patch_v44.py`（自 v43）

### 新增
- 模块 12b 整数采样可视化增强。
- 模块 19 双 y 轴。
- 模块 20：`c/D²` 拟合 + 平方反比起始距离。

---

## v40–v43 —— 分析线 fork（自 v32）

### 新增
- 停用一批叶子演示模块；追加能量扫描 / SNR–距离分析（后称模块 19/20）。
- **首次停用模块 16（100 ppm）**。
- v41：模块 **0b** 独立时间窗（实际 `dt_fine=200 ps`）。
- v41+：展示降强 `DEMO11B_SCALE` / `MOD12_SCALE`。

### 脚本
- `build_v40_from_v32.py`、`build_v41_from_v32.py`、`build_v42_from_v32.py`、`build_v43_from_v32.py`。

---

## v30–v32 —— 二值采样分水岭

### 新增
- **v30**：统计范式从 timestamp 改为 **1 ns HDC 二值**；`spad_binary_trace`；过阈窗 ≈8 ns；`macro_cap=27×N_shots`。
- **v31**：`return_attrib`；左扩 guard band 修边界假象。
- **v32**：信号/环境按比例公平归因；模块 17/18（后被 19 取代）。

### 踩坑
- 无护带 → 窗左端背景系统性偏低。
- 二值硬上限导致强信号削顶；加 shot 不缓解饱和占比。

---

## v20–v21 —— RC 恢复与检测链

### 新增
- v20：`Vth_frac` 10%→60%；`τ_RC=8.7315 ns`；`g(Vov)`；模块 13 SNR。
- v21：模块 14 前沿定时、15 SNR 分布、16 **100 ppm** 噪点率（后停用）。

### 脚本
- `build_v20_from_v15.py`、`build_v21_from_v20.py`。

---

## v10–v15 —— 单 SPAD 蒙卡扩展

### 新增
- 在解析光链路上加逐光子 MC、IRF、阵列、RC 接口等（细节见各版 notebook）。

---

## v0–v4 —— 解析光链路

### 新增
- 激光源 → TX → 大气/光斑/ToF → 朗伯目标 → RX 链路预算的解析框架。
