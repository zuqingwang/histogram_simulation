# 交接文档 —— 当前工作：`lidar_histogram_sim`（激光雷达直方图物理仿真）

> 本文档文件名：`handoff_lidar_histogram_sim.md`（工作名 = `lidar_histogram_sim`；禁止 `handoff_现在工作.md`）。
> 写给**完全没有上下文的新会话**。文件名、模块号、变量名、参数值均写全。
> 最后更新：2026-08-04。
> 配套流水日志：`worklog_lidar_histogram_sim.md`。
> 兄弟工作线：`handoff_PoD_esti.md`（探测概率）、`handoff_crosstalk_sim.md`（串扰 + tcode）。

---

## 0. 缩写表

| 缩写 | 英文全称 | 含义 |
|---|---|---|
| SPAD | Single-Photon Avalanche Diode | 单光子雪崩二极管 |
| PDE | Photon Detection Efficiency | 光子探测效率（本项目 `PDE_max = 0.30`） |
| ToF | Time of Flight | 飞行时间 |
| IRF | Instrument Response Function | 仪器响应函数（高斯 σ = 100 ps） |
| HDC | Hardware Data Channel | 硬件数据通道：1 ns 时钟采样 → 比 60% 阈值 → 出 0/1 |
| Vov | Over Voltage | 过电压 = V_bias − V_br（`Vov_max = 3.3 V`） |
| FAR | False Alarm Rate | 虚警率 / 噪点率 |
| SNR | Signal-to-Noise Ratio | 信噪比 |
| BRDF | Bidirectional Reflectance Distribution Function | 双向反射分布函数 |
| DCR | Dark Count Rate | 暗计数率（本项目 = 0） |

---

## 1. 我们在做什么任务

仓库：`E:\claude temp\Histogram-simulation`。

工作线 **A：激光雷达直方图物理仿真**，主产物 **`lidar_histogram_sim_v45.ipynb`**。

目标：从激光源 → TX → 大气/目标 → RX → SPAD 阵列 → **1 ns 二值采样直方图**，给出可解释的测距/SNR/能量扫描分析。  
这是整仓的**物理内核来源**；`PoD_esti` 从本线提取 SPAD + 二值引擎；`crosstalk_sim` 为独立并行线（理想 δ 回波 + 编码滤噪），不要交叉改参数。

---

## 2. 已经完成了什么

### 2.1 版本地图（摘要）

| 阶段 | 文件 | 要点 |
|---|---|---|
| 解析光链路 | v0–v4 | 激光→TX→信道→目标→RX |
| 单 SPAD 蒙卡 | v10–v15 | 逐光子 MC、IRF、阵列、RC |
| RC + SNR | v20–v21 | `Vth_frac=0.60`、`τ_RC=8.7315 ns`、模块 13–16 |
| **二值分水岭** | **v30** | timestamp → **HDC 二值**；`spad_binary_trace`；`macro_cap` |
| 护带 + 归因 | v31–v32 | 左扩 guard band；信号/环境公平归因 |
| 分析线 fork | v40–v45 | 停用叶子模块；追加模块 19/20 |
| **当前** | **`lidar_histogram_sim_v45.ipynb`** | 模块 19 增强（底噪诊断、指数 boost、扣背景重心） |

构造脚本：`build_v*_from_*.py`、`patch_v44.py`（v45 为 v44 复制后手改）。

### 2.2 v45 关键模块（60 cell；勿整文件硬读，约 3.6 MB）

| 模块 | Cell（约） | 内容 |
|---|---|---|
| 0 / 0b | 2–3 / 25 | `PARAMS`；全局仿真/绘图窗（`dt_fine` 实际 200 ps） |
| 1–5 | 4–13 | 激光 / TX / 信道 / 目标 / RX |
| 6 / 6b | 14–17 | 单 SPAD 光子率；阵列俯视图 |
| 7 / 7b / 7c | 18–23 | 硬死时间 / RC timestamp 引擎 / g(Vov)（7c 仍执行） |
| **9b** | **31–32** | **核心：`spad_binary_trace`（二值引擎）** |
| 9c | 33–34 | 单 SPAD 二值折线 |
| 10 / 11 / 11b | 35–40 | 宏像元 9×3；全宏直方图+归因；展示折线 |
| 12 / 12b | 45–48 | 热图；亚 ns 前沿测距 |
| 13 / 14 | 49–52 | SNR；检测阈值 + 前沿定时 |
| **19** | **58** | 能量扫描 → 前沿/重心 → peak/area |
| **20** | **59** | SNR & 信号率 vs 距离（`c/D²` 拟合） |

**v43 起停用（转 markdown，不执行）**：模块 8、8b、9、17、15、**16（100 ppm）**。

### 2.3 当前关键参数（不得擅自改）

```
λ = 905 nm；P_peak = 235 W；τ_r=0.7 ns，τ_f=1.9 ns
主目标：D = 30 m，ρ = 0.10（ToF ≈ 200.1 ns）
环境光：E_lambda = 0.68 W/m²/nm，enable=True
SPAD：PDE=0.30，DCR=0，IRF σ=100 ps，τ_RC=8.7315 ns，Vth_frac=0.60
直方图：bin=1 ns，N_shots=4；仿真窗 ToF−50~+100 ns
宏像元：9×3=27 SPAD；macro_cap = 27×4 = 108
护带：guard_band ≈ T_OVER + 5×jitter ≈ 8.5 ns（tf_gen 左扩，centers 不动）
```

注意：`PARAMS["hist"]["dt_fine"]=10 ps` 仍写在字典里，**实际 MC 用模块 0b 的 200 ps**。

### 2.4 与兄弟线关系

| 线 | 关系 |
|---|---|
| `PoD_esti` | 提取 cell 3 物理参数 + 光链路 + cell 32 `spad_binary_trace`；FAR/PoD 由 PoD 线承担 |
| `crosstalk_sim` | 独立；共享绘图风格约定，无代码 import |

---

## 3. 目前面临的问题、卡在哪里

1. **模块 16（100 ppm）已停用**：Poisson 窗口级 FAR 假设不严谨；正式 FAR/PoD 走 `PoD_esti`。
2. **模块 19 自检**：二值 8 ns 过阈窗把脉冲涂成方波 → walk 形态与模拟前端 RC 曲线不一致；线性 boost 加密 ≠ peak/area 均匀。
3. **notebook 体积大**：禁止整文件 Read；用 `json.load` 按 cell 导出 source。
4. **尚无本线专用缓存/多线程规则落地到 v45 本身**（规则三目前主要在 `PoD_esti` v05 固化）；若重跑模块 19/20 长 MC，应补缓存。
5. 本 handoff/worklog **本次会话才补建**；历史版本细节以 notebook markdown 与 `worklog_lidar_histogram_sim.md` 为准。

---

## 4. 下一步计划

1. 若继续物理分析：优先模块 **19/20**（能量扫描、SNR–距离）；长 MC 加缓存 + 并行。
2. 任何物理参数变更：先说明是否同步 `PoD_esti`（PoD 线默认「未擅自改 v45 参数」）。
3. 不要为了 PoD/串扰任务去改本 notebook 的核心物理值，除非用户明确要求。
4. 可选：为模块 19/20 抽独立轻量脚本，避免反复打开 3.6 MB notebook。

---

## 5. 踩过的坑（不要再踩）

1. **边界假象（「最开始没噪声」）**：二值引擎在 `[t_av, t_av+T_OVER]` 置 1；若光子生成窗从 `t_lo` 起，前 ~8 ns 缺上游雪崩尾巴 → 背景偏低。**修法**：`tf_gen` 左扩 `guard_band ≥ T_OVER`，采样 `centers` 不动。
2. **二值硬上限**：单 SPAD 每 bin 每 shot 至多 1；宏像元 `macro_cap=108`；加 shot 不能缓解饱和占比（分子分母同乘）。
3. **`dt_fine` 双轨**：PARAMS 写 10 ps，0b 用 200 ps——读参不要只看 PARAMS。
4. **展示降强**：`DEMO11B_SCALE=0.3`、`MOD12_SCALE=0.5` 仅展示，不影响 11/13/14/19/20 分析。
5. **100 ppm 模块不可直接当 PoD 真值**：已停用；FAR 用 `PoD_esti` MC。
6. **大 notebook**：先 `json` 读 cell source，勿整文件硬读。
7. **float32 坑在 PoD 线**（非本 notebook）：`rng.random(float32)` 可精确为 0；PoD 快速引擎已改 `standard_exponential`。

---

## 6. 给新会话的最短指令

1. 打开 `lidar_histogram_sim_v45.ipynb` + 本文件 + `worklog_lidar_histogram_sim.md`。
2. 物理内核看 **模块 9b / cell 32 `spad_binary_trace`**；分析看模块 19/20。
3. 改 FAR/PoD → 去 `PoD_esti`；改编码/滤串扰 → 去 `crosstalk_sim`。
4. 会话结束前更新本文件与对应 worklog。
