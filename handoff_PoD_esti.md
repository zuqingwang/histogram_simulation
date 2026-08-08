# 交接文档 —— 当前工作：`PoD_esti`（探测概率估计）

> 本文档文件名：`handoff_PoD_esti.md`（工作名 = `PoD_esti`，与主产物 `PoD_esti_v*.ipynb` 一致；禁止 `handoff_现在工作.md`）。
> 本文档写给**完全没有上下文的新会话**。所有文件名、模块号、变量名、参数值均写全，不使用会话内指代。
> 最后更新：2026-08-08（`PoD_esti` **v10**：per-shot `hist_i` / `hist_add` 重算；**基于 v05，不基于 v06**）。
> 配套流水日志：`worklog_PoD_esti.md`。
> 兄弟工作线：`handoff_lidar_histogram_sim.md`、`handoff_crosstalk_sim.md`、`handoff_peak_vs_noise.md`。

---

## 0. 缩写表

| 缩写 | 英文全称 | 含义 |
|---|---|---|
| SPAD | Single-Photon Avalanche Diode | 单光子雪崩二极管 |
| PDE | Photon Detection Efficiency | 光子探测效率（`PDE_max = 0.30`） |
| ToF | Time of Flight | 飞行时间 |
| FAR | False Alarm Rate | 虚警率 / 噪点率 |
| PoD | Probability of Detection | 探测概率 |
| MC | Monte Carlo | 蒙特卡洛 |

---

## 1. 我们在做什么任务

仓库：`E:\claude temp\Histogram-simulation`。

当前工作线 **A′：探测概率估计**，主产物 **`PoD_esti_v10.ipynb`**（基于 **v05** 分叉，**不要**基于 v06）。

### v10 目标
在 **per-shot 直方图**架构下全量重算（**不复用** v05/v06 缓存），并回答三张新图相关的问题：

1. 纯噪声 **peak vs bg**（N=1/2/4）曲线形状是否一致；
2. **bg + 5·std(peak)** 与 **1% FAR 阈值** 的关系；
3. 固定信号、noise 线性增长时，`hist_add` 上 peak 分布是否仅平移，均值/std 是否线性。

### 口径（v10）

| 符号 | 含义 |
|---|---|
| `hist_i` | 第 i 发、宏像元 27 SPAD 的二值累加直方图（单发计数 ≤27） |
| `hist_add(N)` | `sum(hist_1…hist_N)`；N∈{1,2,4} 由 **同一次 N=4 仿真的前缀和**得到 |
| **noise** | 单次 `hist_i` 统计窗均值（环境标准；与 N 无关） |
| **bg** | `hist_add` 统计窗均值 |
| **peak** | 一律在 `hist_add` 上统计（纯噪声→统计窗；含信号→信号窗） |

冒烟已验证：目标 noise=0.25/0.50/0.75 时，实测 bg(N=1/2/4) ≈ 1×/2×/4× noise。

---

## 2. 已经完成了什么

### 2.1 版本地图

| 版本 | 文件 | 状态 |
|---|---|---|
| v05 | `PoD_esti_v05.ipynb` | 计算基线；6 档 FAR；PoD 全量已跑完 |
| v06 | `PoD_esti_v06.ipynb` | 仅作图口径（noise/bg 双轴）；缓存仍用 v05；**v10 不基于它** |
| **v10** | **`PoD_esti_v10.ipynb`** | **当前主产物**：`hist_i` 架构 + 三张新图；**新缓存重算** |

### 2.2 关键文件（v10）

| 文件 | 作用 |
|---|---|
| `pod_esti_v10_core.py` | 自 v05 core 分叉；`noise_hists_per_shot` / `binary_macro_stepping_per_shot` / `hist_add_from_prefix` / `stats_from_hist_i`；**import 不跑 MC** |
| `run_pod_v10_scan.py` | 多进程扫描：`noise` / `signal` / `all`；20 进程；增量 `.partial.npz` |
| `build_pod_esti_v10.py` | 生成 `PoD_esti_v10.ipynb`（改图改此文件） |
| `PoD_esti_v10.ipynb` | 读缓存出图（图1/2/3） |
| `pod_esti_v10_cache_noise.npz` | 纯噪声 48 档 × N_MC=200000（新算） |
| `pod_esti_v10_cache_signal.npz` | 固定 boost × 扫 noise（新算） |

### 2.3 参数（不得擅自改物理量）

```
物理参数      与 v05 / v45 一致
宏像元        9×3 = 27 SPAD
N_SHOTS_MAX   4；分析 N_SHOTS_LIST = [1, 2, 4]（前缀和）
NOISE_GRID_AMB  0.25→12 步长 0.25（48 档，按单次 noise）
N_MC_NOISE    200_000 / 档
N_MC_SIG      8_000 / (noise, boost)
BOOST_LIST    [0, 0.004, 0.008, 0.016, 0.032]
并行          ProcessPoolExecutor 默认 20
缓存          pod_esti_v10_cache_*.npz；CACHE_*_FALLBACK = []（禁止旧缓存）
```

### 2.4 运行方式

```powershell
$env:PYTHONIOENCODING="utf-8"
python run_pod_v10_scan.py noise --workers 20          # 纯噪声
python run_pod_v10_scan.py signal --workers 20         # 固定信号
python run_pod_v10_scan.py all --workers 20            # 两者
python run_pod_v10_scan.py noise --limit 3 --n-mc 8000 # 冒烟
python build_pod_esti_v10.py
python -m nbconvert --to notebook --execute --inplace PoD_esti_v10.ipynb
```

---

## 3. 目前面临的问题、卡在哪里

**主线已跑通。** 全量缓存齐，`PoD_esti_v10.ipynb` 已 execute。

定量结论摘要：
- bg/noise 中位 = **1.00 / 2.00 / 4.00**（N=1/2/4）。
- 图1：按 N 归一后 peak–bg **不完全重合**（相对 N=1 的 RMS：N=2→1.52，N=4→2.59 计数）。
- 图2：T@1% 系统性高于 bg+5σ（mean Δ ≈ +3.1/+4.3/+5.8）→ **5σ 经验规则偏松**。
- 图3：`peak_mean`–noise R²≈0.997–1.000，但斜率随信号增强下降（N=1：1.20→0.68）→ **非可加纯平移**。

尚未做：v10 架构上的完整六档 FAR / PoD50·90 能量扫描。

---

## 4. 下一步计划

1. 用户审阅三张图与上述结论，确认是否接受「归一后形状不完全一致 / 5σ 偏松 / 均值近线性但斜率随信号降」。
2. （可选）在 v10 的 `hist_i` 架构上接回 FAR 六档与 PoD 临界能量。
3. （可选）把图1归一残差做成对 bg 的函数曲线，定位差异主要出现在高噪声端还是低噪声端。

---

## 5. 踩过的坑（不要再踩）

1. **v10 必须基于 v05，不要从 v06 复制**：v06 只改作图轴，没有 `hist_i`。
2. **不要复用 `pod_esti_v05_cache_*.npz`**：口径与派生方式已变；`CACHE_*_FALLBACK=[]`。
3. **N=1/2 不要各自独立扫**：始终仿 N=4，前缀和得 N=1/2，省算力且共享随机实现。
4. **`r_det_for_noise(noise, n_tr)` 必须传 `n_tr=27`**（单次宏像元），不是 `27×N`。
5. **检查点保存必须容忍稀疏 `rows`**（完成一部分就落盘时会缺 key）。
6. **import `pod_esti_v10_core` 不应触发百万条 MC**：扫描只走 `run_pod_v10_scan.py`。
7. 编辑器改完 `.ipynb` 后必须 Revert/重开，避免旧内存副本覆盖磁盘。
