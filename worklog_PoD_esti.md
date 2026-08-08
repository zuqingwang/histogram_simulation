# worklog —— `PoD_esti`（探测概率估计）

> 顶部 = 现状与当前任务（可覆盖更新）；下部 = 历史记录（只追加，不删改）。
> 交接文档见 `handoff_PoD_esti.md`（工作名 `PoD_esti`，禁止再用 `handoff_现在工作.md`）。

---

# 一、现状 / 当前任务

**当前版本：v10（2026-08-08，已跑通）—— 基于 v05，不基于 v06**
**主产物：`PoD_esti_v10.ipynb` + `pod_esti_v10_core.py` + `run_pod_v10_scan.py`**

## 1. 目标一句话（v10）

在 **per-shot `hist_i`** 架构下重算：每次实现记录最多 4 发的 `hist_i`，
`hist_add = sum(hist_1…hist_N)`；**noise** = 单次 `hist_i` 底噪，**bg** = `hist_add` 底噪，
**peak** 一律在 `hist_add` 上统计。三张新图已出。**不复用 v05/v06 缓存。**

## 1b. 口径（v10）

| 符号 | 含义 |
|---|---|
| `hist_i` | 第 i 发、宏像元 27 SPAD 的二值累加直方图（计数 ≤27） |
| `hist_add` | 前 N 发之和；N∈{1,2,4} 由 **同一次 N=4 仿真的前缀和** 得到 |
| **noise** | 单次 `hist_i` 统计窗均值（环境标准；与 N 无关） |
| **bg** | `hist_add` 统计窗均值 |
| **peak** | `hist_add` 上统计窗（纯噪声）或信号窗（含信号）的最大 bin |

## 2. 当前可运行状态

- **全量 MC 已完成**：noise 311 s（48×200k），signal 270 s（48×5×8k），20 进程。
- **`PoD_esti_v10.ipynb` 已 `nbconvert --execute`**，图 `pod_v10_fig1/2/3*.png` 已写出。
- 中位 bg/noise：**N=1→1.00，N=2→2.00，N=4→4.00**（口径自洽）。

## 3. 当前生效关键参数（v10）

```
物理参数      与 v05 / v45 一致（未擅自改）
宏像元        9×3 = 27 SPAD；N_SHOTS_MAX=4；分析 N∈{1,2,4}
噪声网格      按单次 noise：0.25→12 步长 0.25（48 档）；N=2/4 由前缀和派生
N_MC_NOISE    200_000 / 档
N_MC_SIG      8_000 / (noise, boost)
BOOST_LIST    [0, 0.004, 0.008, 0.016, 0.032]
并行          ProcessPoolExecutor，默认 20 进程
缓存          pod_esti_v10_cache_noise.npz / pod_esti_v10_cache_signal.npz
```

## 4. 待办

- [x] `pod_esti_v10_core.py`：`hist_i` + 前缀和 API
- [x] `run_pod_v10_scan.py`：纯噪声 + 固定信号扫描（新缓存）
- [x] `PoD_esti_v10.ipynb`：三组新图 + execute
- [ ] （可选）v10 架构上接回六档 FAR / PoD50·90
- [ ] （后续）固件阈值查找表；DCR 等

## 5. v10 关键结论（简）

1. **peak–bg 形状**：按 N 归一后 N=2/4 相对 N=1 的 RMS 残差约 **1.52 / 2.59** 计数，**不完全重合**（有系统弯曲差异）。
2. **bg+5σ vs T@1%**：T 系统性 **高于** bg+5σ（mean Δ ≈ +3.1/+4.3/+5.8），经验规则偏松。
3. **固定信号**：`peak_mean` vs noise 的 R²≈0.997–1.000（很直）；但斜率随信号增强下降（N=1：1.20→0.68），**不是可加平移**；std 的线性更差。

## 5. 仍有效的 v02 定量结论（抽样）

| N_shots | noise | T@100ppm | T@10ppm | E@PoD90(100ppm) [nJ] | PoD90 等效距离 [m] |
|---|---|---|---|---|---|
| 1 | 1 | 9 | 10 | 23.6 | 86.6 |
| 1 | 5 | 17 | 19 | 89.9 | 44.6 |
| 1 | 10 | 23 | 24 | 634 | 16.8 |
| 4 | 1 | 9 | 10 | 4.10 | 205.6 |
| 4 | 5 | 19 | 21 | 7.63 | 151.4 |
| 4 | 10 | 28 | 30 | 10.5 | 129.4 |

完整表与生存函数见 `PoD_esti_v02.ipynb` / `worklog` 历史 v02 段。

---

# 二、历史记录（只追加，不删改）

## v10 —— 2026-08-08　per-shot hist_i / hist_add 全量重算（基于 v05）

**新增**
- `pod_esti_v10_core.py`：`noise_hists_per_shot`、`binary_macro_stepping_per_shot`、前缀和统计；import 不跑 MC。
- `run_pod_v10_scan.py`：noise/signal/all；禁止旧缓存；稀疏检查点。
- `build_pod_esti_v10.py` → `PoD_esti_v10.ipynb`（图1 peak–bg；图2 bg+5σ vs T@1%；图3 固定信号）。
- 缓存：`pod_esti_v10_cache_noise.npz`（48×200k，311 s）、`pod_esti_v10_cache_signal.npz`（48×5×8k，270 s）。

**未改**
- 物理参数（D、ρ、脉冲、bin 宽等）与 v05 一致。
- 不基于 v06；不读 v05 缓存。

**运行结论**
- bg/noise 中位 = 1/2/4（N=1/2/4）。
- 归一 peak–bg：N=2/4 相对 N=1 RMS = 1.52 / 2.59。
- T@1% − (bg+5σ) 均值 = +3.13 / +4.29 / +5.78（N=1/2/4）→ 5σ 规则偏松。
- 固定信号：mean–noise R²≥0.997，但斜率随 E 下降（死时间/饱和）；非纯平移。

**踩过的坑**
- `r_det_for_noise` 必须传 `n_tr=27`。
- 检查点 `_save_noise` 不能假设 `rows` 含全部 k（KeyError: 4）。
- 信号路径用 `binary_macro_stepping_per_shot` + 前缀和，勿对 N=1/2 各扫一遍。

## v06 —— 2026-08-08　noise / bg 双轴出图

**新增**
- `PoD_esti_v06.ipynb`（自 v05 复制）。
- 口径说明 markdown + 派生 `noise_ambient` / `bg` / `_axis_x` / `_pod_x_from_rec`。
- 模块 5（noise–peak）、5b（密度条带）、6（阈值）、8（临界能量汇总）经 `AXIS_KINDS` 各出 noise 轴与 bg 轴一套图；savefig 带 `_noise` / `_bg` 后缀（`pod_v06_*.png`）。
- 模块 7 横轴是能量，**不做**双循环；图例改为标 `bg=…/noise≈…`。

**未改**
- 物理参数、FAR、MC 条数、缓存文件名（仍用 v05 缓存）。

**踩过的坑**
- 误把模块 7（能量轴）也包进 `AXIS_KINDS` 双循环 → 两套图完全相同；已撤销。
- 模块 8 仅改 `set_xlabel` 不够：`collect_critical` 必须按轴从 `e_lambda` 折合 ambient noise，否则 noise 轴仍是 bg 数值。
- `fig.savefig` 与 `plt.savefig` 都要改，否则轴后缀漏掉（5b 曾漏改 `fig.savefig`）。

## v05-mp —— 2026-08-04　改用多进程，彻底解决 CPU 吃不满

### 现象
- 线程版调到 `POD_BIN_WORKERS=20` 后 CPU 反而只剩 **16%**（约 3 核）。

### 根因
- `binary_macro_stepping` 是「Python 层逐细网格步 `for` 循环 + 小数组 NumPy」结构
  （PoD 子窗 361 步，单块数组仅 250×27 或 250×108）。
- 这种形态绝大部分时间持有 GIL（Global Interpreter Lock，全局解释器锁），
  **线程池加再多路也换不出并行**，只能吃到十几到三十几个百分点。
- 结论更正：v05-speed 条目里「外层×内层 4×5 可得 3–5×」的估计过于乐观，
  实际线程方案对本内核基本无效。

### 新增
- **`build_pod_core_v05.py`**：从 notebook 提取计算 cell（2/4/6/8/9/17/22/25）
  生成 **`pod_esti_v05_core.py`**（可被子进程 import 的纯计算内核，不含绘图）。
- **`run_pod_scan_v05.py`**：`ProcessPoolExecutor` 多进程全量扫描，
  支持 `--workers` / `--limit` / `--checkpoint-every`，断点续跑，实时 flush 进度。
- notebook cell 25 重算分支加提示：全量重算请改用 `run_pod_scan_v05.py`。

### 关键实现约束
- Windows 用 spawn 启动子进程，worker 代码必须在**可 import 的模块级**，
  不能是 notebook 内的闭包，所以才需要 `pod_esti_v05_core.py`。
- 子进程 import 内核会载 noise 缓存（约 5 s）；
  `POD_CORE_QUIET=1` 静音，避免 20 份启动日志刷屏。
- 主进程必须先确保 `pod_esti_v05_cache_noise.npz` 存在，
  否则 20 个子进程会各自触发一遍全量噪声 MC。

### 实测
- 冒烟：6 档 / 6 进程 ≈ 0.5 min（含约 5 s 内核 import）。
- **全量完成：193 档 / 20 进程 = 24.9 min，CPU 稳定 100%**；
  `pod_esti_v05_cache_pod.npz` 已含 **208/208 档 × 六档 FAR**。
- 对照：线程版跑同样内容 CPU 只有 16–30%，估计需数小时。

### 同批修掉的残留 bug（模块 7 验证图 cell 26）
- `tag = "100"`、`r["T_map"]["100"]` 仍是 FAR 扩档前的旧键 → 必然 `KeyError: '100'`。
- 已改为 `REF_FAR_TAG = "100ppm"`，并在标题注明「以 100 ppm 为例，六档汇总见模块 8」。

### 踩坑：脚本改 .ipynb 后编辑器不会自动 reload
- 现象：磁盘上 cell 26 已是 `REF_FAR_TAG`，但运行仍报 `KeyError: '100'`（traceback 里还是 `tag = "100"`）。
- 根因：Cursor / VS Code 的 notebook 编辑器持有自己的内存副本，用脚本改磁盘文件后**不会自动刷新**。
- 风险：此时在编辑器里按 Ctrl+S，会把旧内存副本写回磁盘，**覆盖掉所有脚本修改**。
- 正确做法：脚本改完 `.ipynb` 后，关闭该 notebook 标签页（不保存）再重新打开，或用命令面板的
  `File: Revert File`，然后 Restart Kernel 重跑。

### 踩坑：噪声缓存与 PoD 缓存的 `res` 结构不同
- `pod_esti_v05_cache_noise.npz` 的 `res` 是 `{n_shots: {...数组...}}`，
  顶层只有 2 个键（1 和 4），**真实档数看 `grid_key` 长度（208）**。
- `pod_esti_v05_cache_pod.npz` 的 `res` 是 `{(n_shots, noise): {...}}`，顶层就是 208 个键。
- 用 `len(res)` 去数噪声缓存会误判成「只剩 2 档」，别据此重跑噪声 MC。

---

## v05-cpu —— 2026-08-04　修复 CPU 占用只有 ~30%

### 原因
- 嵌套线程池 `4×5`：GIL 下有效并行常塌成约 5 路 → 总占用 ≈25–35%。
- `CHECKPOINT_EVERY=1`：每档压缩写整个 PoD 缓存，主线程长时间占锁/写盘，worker 空转。
- OpenBLAS/MKL 默认再开多线程，与自管线程池互相掐。

### 改动
- `POD_BIN_WORKERS = N_WORKERS`（20），`POD_WORKERS = 1`（单层外层并行）。
- `CHECKPOINT_EVERY = 8`；快照在锁外写盘。
- `import numpy` **之前**设置 `OMP/MKL/OPENBLAS_NUM_THREADS=1`。

### 操作
- 必须 **Restart Kernel** 后再跑（环境变量对已加载的 numpy 无效）。
- 可用 `pod_esti_v05_cache_pod.partial.npz` 断点续跑。

---

## v05-cache-win —— 2026-08-04　修复 Windows 原子写缓存 WinError 2

### 踩坑 / 修法
- 错误做法：`tmp = path + ".tmp"` + `np.savez_compressed(tmp)` → 实际写出 `*.tmp.npz`，`os.replace(tmp, path)` 报 `FileNotFoundError: WinError 2`。
- 正确做法：`tmp = path + ".tmp.npz"`。
- 已改：`PoD_esti_v05.ipynb` cell 17 `_atomic_savez`；`upgrade_pod_esti_v05_cache_parallel.py` 同步。

---

## docs-sibling —— 2026-08-04　补建兄弟线 handoff / worklog

### 新增
- `handoff_lidar_histogram_sim.md` / `worklog_lidar_histogram_sim.md`
- `handoff_crosstalk_sim.md` / `worklog_crosstalk_sim.md`（含 tcode_calculator 与 `docs/tcode/`）
- `handoff_PoD_esti.md` 增加兄弟线交叉引用。

---

## v05-speed —— 2026-08-04　PoD 加速（不改 MC 精度参数）

### 新增 / 改动
- `POD_BIN_WORKERS=4`：外层同时跑 4 档 noise；`POD_WORKERS=N_WORKERS//4=5` 给档内。
- `POD_MC_CHUNK=250`：大 MC（尤其 `N_MC_POD_VERIFY=5000`）切块后统一进线程池。
- `_verify_critical_batch`：六个 FAR × PoD50/90 的临界验证改成批量并行（最多两轮），不再逐点串行。
- `_eval_mc_jobs`：粗扫/局部/验证共用同一扁平任务队列，避免嵌套线程池。
- 脚本：`upgrade_pod_esti_v05_speed.py`

### 未改
- `N_MC_POD_COARSE/LOCAL/VERIFY`、噪声网格、物理参数、FAR 六档集合均未动。

### 预期
- 相对加速前 PoD：约 **3–5×**；全量 PoD 粗估从 5–6 h 降到约 **1–2 h**（视 CPU/内存带宽而定）。
- 噪声有缓存时仍接近秒级。

---

## v05-FAR-order —— 2026-08-04　FAR_SPECS 顺序更正（从严到松）

### 改动
- `FAR_SPECS` 顺序改为：`10 ppm → 100 ppm → 0.1% → 0.5% → 1% → 5%`
- 标签集合未变，仅排列与图/表遍历顺序变化；PoD 缓存键 `far_tags` 顺序随之变化，旧顺序缓存会自动失效。

### 更正
- 前条「v05-FAR」所写顺序 `100ppm, 10ppm, 5%, …` 作废，以本条为准。

---

## v05-FAR —— 2026-08-04　扩展虚警率阈值至六档；PoD 缓存失效

### 新增
- `FAR_SPECS`：`(值, ASCII标签, 图例标签)`  
  `100ppm / 10ppm / 5pct(5%) / 1pct(1%) / 0p5pct(0.5%) / 0p1pct(0.1%)`
- `FAR_TAG` / `FAR_LABEL` / `FAR_TAGS`：统一键名，禁止再用 `f"{far*1e6:.0f}"`
- 升级脚本：`upgrade_pod_esti_v05_far.py`
- 模块 5b / 6 / 7 / 8：阈值曲线、表、PoD 临界点均按六档展开
- PoD 缓存键增加 `far_tags` 校验；不匹配则整份作废重算

### 删减 / 停用
- **`CACHE_POD_FALLBACK` 清空**：不再读 `pod_esti_v04_cache_pod.npz`（旧缓存只有 100/10 ppm）
- 仍保留 **noise** 的 v04 fallback（纯噪声 MC 与 FAR 无关）

### 参数改动
| 项 | 由 | 到 |
|---|---|---|
| TARGET_FARS / FAR_SPECS | `[100e-6, 10e-6]` | 六档：100ppm, 10ppm, 5%, 1%, 0.5%, 0.1% |
| CACHE_POD_FALLBACK | `[pod_esti_v04_cache_pod.npz]` | `[]` |

### 踩过的坑
- 百分比 FAR 若用 `far*1e6` 做标签会得到荒谬整数；必须用显式 ASCII tag。
- 模块 7 验证图曾硬编码 `critical.get("100")`，已改为 `"100ppm"`。
- `critical` 种子曾写 `int(tag)`，标签改成 `100ppm` 后会立刻崩溃；已改为 `FAR_TAGS.index(tag)`。
- **`np.savez_compressed(path+".tmp")` 在 Windows 会再追加 `.npz`** → `os.replace` 报 WinError 2；临时文件必须写成 `path+".tmp.npz"`。

### 运行结论
- 代码语法 OK；**尚未启动全量 MC**。
- 用户要求后续模拟重跑：noise 可迁移 v04→v05；**PoD 必须重跑**。
- 并行说明：`N_WORKERS=20` 只加速「同档内的分块/能量点」；noise 档与 PoD 档外层仍串行；临界验证 `N_MC_POD_VERIFY` 不走线程池。
- 粗估（加速后）：noise 有缓存 ≈ 秒级；PoD 全量六档约 **1–2 小时**。

---

## v01 —— 2026-08-03　`PoD_esti.ipynb` 首个可运行版本

### 新增
- **`build_pod_esti.py`**：生成 `PoD_esti.ipynb` 的构造脚本（25 cell）。
- **`PoD_esti.ipynb`**：模块 0–8 全流程。
  | 模块 | 内容 |
  |---|---|
  | 0 | 参数与常数（v45 物理参数 + `PoD_esti` 专用参数区） |
  | 1 | 光链路移植（脉冲、发射/接收光学、大气、朗伯目标、像元收集比例、环境光速率） |
  | 2 | 时间窗、bin、掐头去尾、宏像元、`T_OVER` |
  | 3a | 精确引擎 `spad_binary_trace`（v45 cell 32 原样） |
  | 3b | 快速引擎 A `noise_macro_hist_fast` + 快速引擎 B `binary_macro_stepping` |
  | 3c | **三引擎一致性验证**（每 bin 均值 / peak 分布 / bin 间自相关） |
  | 4 | 纯噪声单条波形演示（标注掐头去尾、noise、peak） |
  | 5 | **第 1 步**：噪声强度扫描 → noise–peak 曲线 |
  | 6 | **第 2 步**：100 ppm 阈值求解 → noise–threshold 曲线（含独立 Binomial 保守对照） |
  | 7 | **第 3 步**：能量扫描 → PoD 曲线 → E@PoD50 / E@PoD90 |
  | 8 | 汇总表 + 能量门槛/等效距离总图 |
- **`run_notebook.py`**：基于 `nbclient` 的 notebook 执行器（本机 `jupyter nbconvert` CLI 不可用）。
- 输出：`pod_esti_engine_check.png`、`pod_esti_noise_waveform.png`、`pod_esti_noise_peak.png`、
  `pod_esti_threshold.png`、`pod_esti_pod_curves.png`、`pod_esti_summary.png`；
  缓存 `pod_esti_cache_noise.npz`、`pod_esti_cache_pod.npz`。

### 设计决策（本次会话与用户确认）
| 编号 | 决策 |
|---|---|
| D1 | 波形取 **N_shots = 1 与 4 两种都做**，画对比 |
| D2 | 保持 0–200 ns 采集窗，**目标距离改为 15 m**（ToF ≈ 100 ns，居中）；仅本文件生效，不动 v45 |
| D3 | 噪声强度按 **`E_lambda` 倍数** 对数扫描 1e-2 … 1e1 |
| D4 | PoD 判据 = **峰值 ≥ T 且峰位落在信号窗内**（近似为"信号窗内最大 bin ≥ T"） |
| D5 | 能量横轴用 **等效单脉冲发射能量 [nJ]**，副列给 boost / ρ_eff / 入射光子数 |
| D6 | 100 ppm 求解用 **N_MC = 1e6**，结果缓存到 `.npz` |

### 关键技术方案：两个快速引擎
精确引擎是逐光子 Python 循环（RC 恢复有时序依赖），约 5–8 ms/条，1e6 条需 1.5 h，不可接受。解法：

- **快速 A（仅纯环境光）**：环境速率恒定 ⇒ 雪崩序列是**更新过程（renewal process）**，
  瞬时强度 `h(Δ) = r_amb·PDE·g(1−e^{−Δ/τ})`。预先数值积分 `H(Δ)=∫h`，
  再用 `Δ = H⁻¹(−ln U)`（`np.interp` 反函数插值）一次性并行出样。
  **首个间隔特殊处理**：起始时刻完全恢复 ⇒ 服从 `Exp(r_amb·PDE)`，
  与精确引擎 `last = -1e30` 的初值严格对应。
  **覆盖技巧**：bin 被点亮 ⟺ 它之前最近一次雪崩距它 < `T_OVER`；
  把每次雪崩的窗口**在下一次雪崩处截断**得到互不重叠的区间，
  再用「差分数组 + `np.bincount` + `cumsum`」直接累进 `(实现数, bin 数)` 矩阵，无需逐轨迹展开。
- **快速 B（含信号）**：按细网格时间步同步推进，在「实现 × SPAD」两维向量化。
  一步内到达 `n ~ Poisson(μ)` 个光子、各自以 `φ` 触发 ⇒ 该步至少触发一次的概率 `= 1 − e^{−μφ}`，
  且步内至多一次雪崩（首个触发后 Vov=0 ⇒ g(0)=0）。因 `t − t_last` 恒为步长整数倍，
  `φ` 做成查表 `phi[age]`。**与精确引擎逐步严格等价，不是近似。**

### 一致性验证结果（模块 3c，纯环境光基准档，N_shots=4）
| 引擎 | 条数 | 每 bin 均值 | peak 均值 | peak 标准差 | 相邻 bin 相关 |
|---|---|---|---|---|---|
| 精确 | 1500 | 5.1117 | 10.784 | 1.465 | 0.8704 |
| 快速 A | 150000 | 5.0971 | 10.787 | 1.493 | 0.8715 |
| 快速 B | 1500 | 5.0661 | 10.795 | 1.484 | 0.8718 |

带信号（boost=3e-3）：精确峰值 12.783 ± 2.986 vs 快速 B 12.842 ± 3.081，差 0.4σ。
自相关函数从 lag=0 线性降到 **lag = 8 ns（= T_OVER）处精确归零**，三引擎完全重合——
这正是 8 ns 矩形过阈窗应有的三角形自相关。

### 踩过的坑（本版新记录，全部已修）

1. **等效距离算反了方向**（严重）。
   - 错误：`tgt = 1.0/boost`，把"回波弱 160 倍仍可探测"解成"目标要挪到 1.2 m"。
   - 正确：应求 `link(D)/link(D_ref) = boost`，**boost 越小 ⇒ 距离越远**。修后 4 发累加
     100 klux 的 PoD50 等效距离由 1.2 m 变成 188 m（物理上才合理）。

2. **100 ppm 阈值判定被浮点边界误判**。
   - 错误：用 `sf < 1e-4` 比较浮点生存函数。`1.0 - 19998/20000` 算出
     `9.999999999998899e-05`，比 `1e-4` 略小，于是恰好 2/20000 = 1e-4 的阈值被错误接受。
   - 正确：全程用**整数计数**比较 `n_ge < target_far * n`。修后 N_shots=4/E_λ×1 的阈值由 19→20（缩水测试口径）。

3. **`np.interp` 反解 PoD 时按 PoD 排序会乱序**。
   - 错误：`np.interp(level, pod[argsort(pod)], log10(boost)[...])`。MC 涨落让 PoD 非单调时，
     排序会把不相邻的能量点配到一起。
   - 正确：先按 boost 排序、对 PoD 做 `np.maximum.accumulate` 保序（isotonic）修正，
     再在首个跨越点上插值。

4. **误以为"近距离像面光斑更大"**。
   - 错误猜想：把目标从 30 m 挪到 15 m 会让像斑变大、收集比例下降。
   - 实测：15 m 与 30 m 的像斑长轴**都是 798.9 µm，Σf_pix 都是 0.0597，完全相同**。
   - 根因：系统远早于远场（瑞利距离 z_R ≈ 0.40 m ≪ 目标距离），
     光斑 `y_D ≈ 2θ_y·D` 与距离成正比，故像面光斑 `s = (y_D/D)·f_RX = 2θ_y·f_RX` **与距离无关**。
   - 结论：30 m → 15 m 的唯一实质变化是回波强 4 倍（1/D²），结论完全可移植。

5. **理论底噪公式只是低速率近似**。
   - `noise ≈ n_tr·[1−exp(−r_det·T_OVER)]` 假设雪崩是速率 `r_det` 的 Poisson 过程。
   - 但 RC 恢复带来 8 ns 计数死区，雪崩不可能挨太近。`E_λ×10` 档实测底噪 39.23
     vs 理论 41.74（**低 6%**），这是真实的死区抑制，不是 bug。
   - **阈值一律以 MC 实测为准，理论线只做量级核对。**

6. **`MC_CHUNK = 100_000` 撞内存墙**。
   - 现象：快速 A 峰值内存约 1.2 GB，CPU 利用率只有 43%，1e6 条明显掉速。
   - 根因：`chunk × n_tr = 1.08e7` 的 float64 中间数组有约 10 个（每个 86 MB），
     加上 `chunk × (nbins+1)` 的两个 bincount 输出（各 160 MB）。
   - 修法：`MC_CHUNK` 改为 **25_000**（峰值约 300 MB）。改这个值**不会让缓存失效**
     （缓存键只含 `amb_mults` / `n_mc` / `n_shots_list`）。

7. **本机 `jupyter nbconvert` CLI 不可用**。
   - 现象：`jupyter nbconvert` 报 `Jupyter command 'jupyter-nbconvert' not found`，
     但 `import nbconvert, nbclient` 正常。
   - 修法：用 `run_notebook.py`（`nbclient.NotebookClient`）执行，并设 `MPLBACKEND=Agg`。

8. **PowerShell 默认 GBK 编码**（承自 v00）。
   - 含中文/`µ`/emoji 的 `print()` 输出到 stdout 会抛 `UnicodeEncodeError`。
   - 修法：跑任何脚本前先 `$env:PYTHONIOENCODING="utf-8"`。

---

## v00 —— 2026-08-03　立项、读码、写规则与文档

### 新增
- `.cursor/rules/session-handoff.mdc`：全局规则一，每次会话结束必须写/更新 `handoff_现在工作.md`。
- `.cursor/rules/worklog.mdc`：全局规则二，每项工作从开工起维护 `worklog_<工作名>.md`。
- `CLAUDE.md` 追加「文档维护」小节，把上述两条规则同步给 Claude Code。
- `handoff_现在工作.md`：完整交接文档。
- `worklog_PoD_esti.md`：本文件。

### 改动
- 无代码改动。未触碰任何既有 `.ipynb`。

### 读码结论（`lidar_histogram_sim_v45.ipynb`，共 60 cell）
- **唯一核心引擎**是 cell 32 的 `spad_binary_trace()`：光子 Poisson 到达 → 逐光子按
  `p_fire = PDE_max·g(1−e^{−Δt/τ_RC})` 判雪崩（雪崩即复位 Vov）→ 每次雪崩铺 8 ns 过阈窗 →
  1 ns 时钟采样落在窗内即记 1。单 SPAD 每 bin 每 shot ∈ {0,1}。
- 第 1、2 步（纯噪声）**只依赖**：`r_amb_ph`（cell 15）、SPAD 参数（cell 21）、
  `spad_binary_trace`（cell 32）、宏像元 SPAD 数与 `N_shots`。
  **不需要**整条发射/接收光学链路（环境光速率与像元收集比例 `f_pix` 无关）。
- 第 3 步（信号）才需要 cell 5/7/9/11/13/15 的完整链路和 `signal_photon_rate_fine`。
- v45 cell 56（模块 16）里有 100 ppm 阈值的 Poisson 解析推导，但**已被停用且假设不严谨**，只能当对照。

### 踩过的坑
1. **大 notebook 不能整文件读**：v45 有 3.6 MB，绝大部分是 base64 图片。
   正确做法是写小脚本 `json.load` 后只导出 `cell['source']` 到临时 txt 再读。
2. **PowerShell 默认 GBK 编码**：把含中文/emoji 的 cell 源码 `print()` 到 stdout 会抛
   `UnicodeEncodeError`。必须写文件时显式 `encoding='utf-8'`，或设 `PYTHONIOENCODING=utf-8`。
3. **0–200 ns 窗与 30 m 目标冲突**：30 m 的 ToF = 200.14 ns 落在窗外。
   立项时就得定死窗口与目标距离的关系，否则模块 6/7 写到一半才发现。

### 关键数值（承自 v45）
- `T_OVER = 8.001 ns`，`macro_cap = 27 × N_shots`。
- 默认环境档 `r_amb = 2.036e7 ph/s`，`r_det = 6.107e6 cps`，`nc_base ≈ 5.1` 计数/bin（N_shots=4）。
- 底噪必须用带展宽的式子 `μ_win = 1 − exp(−r_det·T_OVER)`，
  不能用 `r_det·bin_width`（会低估约 `T_OVER/bin_width = 8` 倍）。

---

## v04（2026-08-03）—— 从 v02 基线扩噪声网格 + 自适应 PoD 交点

### 新增 / 改动
- `build_pod_esti_v04.py`：以 `PoD_esti_v02.ipynb` 为基线生成 v04。
- `enhance_pod_esti_v04.py`：改写模块 7/8。
- 噪声网格：
  - N_shots=1：noise 0.25→12，步长 0.25（48 档）
  - N_shots=4：noise 0.25→40，步长 0.25（160 档）
- 模块 7：不再用稀疏 24 点能量网格硬插值；改为每档自适应局部加密 + probit，
  并对 PoD50/90 临界能量做独立验证，保存 peak 均值与 `peak_cnt`。
- 模块 8：横轴改为完整 0.25-noise 网格，画出临界能量 / 临界 peak 均值 / 等效距离。
- 缓存：`pod_esti_v04_cache_noise.npz`、`pod_esti_v04_cache_pod.npz`（已跑出，grid_key 长度 208）。
- 并行曾设为 2→8 线程；后续由 v05 升为 20。

### 明确放弃
- v03 的“平衡态起步省暖机”路线：实测不比 50 ns 暖机更快，且曾引入统计偏差；已停用。


---

## v05（2026-08-03）—— 模块 5b noise–peak 密度条带

### 新增
- `build_pod_esti_v05.py`：从 `PoD_esti_v04.ipynb` 生成 v05，不执行任何 cell。
- `PoD_esti_v05.ipynb` 模块 5b：
  - 参考 `C:\Users\wangzuqing\Desktop\noise-peak\np_short.ipynb` cell 14；
  - 每个 0.25-noise 档单独成列；
  - 方块宽度、颜色深度表示该列内 peak 相对概率；
  - 叠加 100 ppm 和 10 ppm 阈值；
  - N_shots=1 与 N_shots=4 分面显示。

### 保留 / 缓存
- 模块 5 原 noise–peak 三联图未删除。
- v05 保持使用 v04 缓存文件名：
  `pod_esti_v04_cache_noise.npz`、`pod_esti_v04_cache_pod.npz`。
- 生成时未运行仿真；v05 代码语法检查通过，所有输出为空。

### 踩过的坑
- 用户要求直接复用缓存，但当前仓库扫描不到任何 `*cache_noise*.npz`。
  因此代码已配置为复用 v04 缓存，但若执行时缓存仍不存在，原模块 5 会自动重算；
  不能声称当前已经具备可复用缓存。

### 更正
- 上一条“缓存不存在”的判断错误：专用 Glob 搜索未返回二进制文件，但
  `Test-Path` 已确认两个文件都存在；`np.load` 也成功读出。
- `pod_esti_v04_cache_noise.npz`：`grid_key.shape=(208,)`，包含完整 v04 noise 网格。
- `pod_esti_v04_cache_pod.npz`：`grid_key.shape=(208,)`，包含完整 v04 PoD 网格。
- v05 保持相同缓存名和缓存键，因此可以直接复用，不需要重跑 MC。

---

## v05b（2026-08-03）—— 20 线程 + 增量缓存规则固化

### 新增 / 改动
- `upgrade_pod_esti_v05_cache_parallel.py`：升级 `PoD_esti_v05.ipynb`。
- `.cursor/rules/pod-cache-parallel.mdc`：全局规则三（耗时仿真必须缓存 + 多线程）。
- `CLAUDE.md` 同步写入规则三。
- `N_WORKERS = 20`，`MC_CHUNK = 5000`。
- 主缓存改名为 `pod_esti_v05_cache_*.npz`；仍可读取 v04 缓存并同步到 v05。
- 增量检查点 `*.partial.npz`：每档落盘，支持断点续跑；写盘使用临时文件 + `os.replace`。

### 未做
- 未重新跑完整 MC；预期首次运行若命中 v04 缓存会直接载入并写出 v05 主缓存。

---

## v02 —— 2026-08-03　噪声线性加密网格 + 双 FAR 阈值（已全量跑通）

### 新增
- uild_pod_esti_v02.py → PoD_esti_v02.ipynb（可 Restart & Run All）。
- 噪声档改为目标 noise 线性网格；用更新过程平衡态反解 E_lambda。
- 同时给出 100 ppm 与 10 ppm 阈值。
- 快速引擎 A：O(1) 反函数直查表 + float32（约 295×）。
- 缓存只存 peak 的 bincount。

### 当时参数（已被 v04/v05 扩展，此处仅作历史）
- N_shots=1：noise 0.25→10，步长 0.25（40 档）
- N_shots=4：noise 0.50→20，步长 0.50（40 档）
- PoD 环境档：NOISE_POD=[1.0,5.0,10.0]
- 全量约 55 分钟；缓存 pod_esti_v02_cache_*.npz

### 关键结论（摘要）
- 100→10 ppm 阈值通常只再抬 1–3 计数。
- N_shots=4 显著降低 PoD90 所需能量（同 noise=10：约 10.5 nJ vs N=1 的 634 nJ）。
- 独立 Binomial 阈值偏保守 0–3 计数。

### 说明
- 本条为 2026-08-04 补记：v02 详细结论原先写在顶部现状区，刷新顶部时改为摘要引用；按规则另起历史条目保留，不删改旧文。

---

## 更正 —— 交接文档命名（2026-08-04）

- **错误做法**：写成 `handoff_现在工作.md`。
- **正确做法 / 现行约定**：`handoff_<工作名>.md`，本工作为 `handoff_PoD_esti.md`。
- 全局规则已改：`.cursor/rules/session-handoff.mdc`、`CLAUDE.md` 规则一。
- 历史条目里若仍出现旧文件名，以本条更正为准，不再回溯篡改旧记录。
