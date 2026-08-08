# 交接文档 —— 当前工作：`PoD_esti`（探测概率估计）

> 本文档文件名：`handoff_PoD_esti.md`（工作名 = `PoD_esti`，与主产物 `PoD_esti_v*.ipynb` 一致；禁止 `handoff_现在工作.md`）。
> 本文档写给**完全没有上下文的新会话**。所有文件名、模块号、变量名、参数值均写全，不使用会话内指代。
> 最后更新：2026-08-08（`PoD_esti` **v06**：noise/bg 双轴出图；计算缓存仍为 v05）。
> 配套流水日志：`worklog_PoD_esti.md`。
> 兄弟工作线：`handoff_lidar_histogram_sim.md`、`handoff_crosstalk_sim.md`、`handoff_peak_vs_noise.md`。

---

## 0. 缩写表

| 缩写 | 英文全称 | 含义 |
|---|---|---|
| SPAD | Single-Photon Avalanche Diode | 单光子雪崩二极管 |
| PDE | Photon Detection Efficiency | 光子探测效率（`PDE_max = 0.30`） |
| ToF | Time of Flight | 飞行时间 |
| IRF | Instrument Response Function | 仪器响应函数（高斯 σ = 100 ps） |
| FAR | False Alarm Rate | 虚警率 / 噪点率 |
| PoD | Probability of Detection | 探测概率 |
| MC | Monte Carlo | 蒙特卡洛 |
| ppm | parts per million | 100 ppm = 1e-4；10 ppm = 1e-5 |

---

## 1. 我们在做什么任务

仓库：`E:\claude temp\Histogram-simulation`。

当前工作线 **A′：探测概率估计**，主产物 **`PoD_esti_v06.ipynb`**（作图口径）；计算/缓存基线仍是 **`PoD_esti_v05`**。  
从 `lidar_histogram_sim_v45.ipynb` 提取 SPAD + 二值采样内核，回答：

1. 纯噪声下，不同噪声强度的 peak 分布与 peak 曲线（**同时按 noise 与 bg 两套横轴出图**）；
2. 使 FAR 分别低于 **10 ppm / 100 ppm / 0.1% / 0.5% / 1% / 5%** 的整数阈值 T；
3. 在各 T 下，PoD = 50% / 90% 所需能量、临界 peak 均值、等效距离。

### 口径（v06，与 `peak_vs_noise` v02 一致）
| 符号 | 含义 |
|---|---|
| **noise** | 环境标准：折合 **N_shots=1、27 SPAD、每 1 ns bin** 的平衡态底；与发数无关 |
| **bg** | 当前 N_shots 下统计窗实测 baseline（=`noise_mc`）；N=1 时 ≈noise，N=4 时 ≈4·noise |

旧字段名 `noise` / `noise_mc` / `NOISE_GRID` 多数对应 **bg（累加后的底）**。

---

## 2. 已经完成了什么

### 2.1 版本地图

| 版本 | 文件 | 状态 |
|---|---|---|
| v01 | `PoD_esti.ipynb` | 保留 |
| v02 | `PoD_esti_v02.ipynb` | **已跑通**；稀疏 noise/PoD 档；55 min 全量 |
| v03 | `PoD_esti_v03.ipynb` | 探索失败，勿继续用 |
| v04 | `PoD_esti_v04.ipynb` | 扩到完整 0.25-noise 网格 + 自适应 PoD；**缓存已产出（仅 100/10 ppm）** |
| v05 | `PoD_esti_v05.ipynb` | 6 档 FAR；20 线程；增量缓存；**PoD 全量已跑完** |
| **v06** | **`PoD_esti_v06.ipynb`** | **当前主产物（作图）**：在 v05 上拆分 noise/bg 双轴；**缓存仍用 v05** |

### 2.2 关键文件

| 文件 | 作用 |
|---|---|
| **`PoD_esti_v06.ipynb`** | **当前作图主产物**（noise + bg 双轴） |
| `PoD_esti_v05.ipynb` | 计算/缓存基线；FAR 六档 |
| `upgrade_pod_esti_v06_bg.py` / `fix_pod_esti_v06_axes.py` | v06 升级与轴修复 |
| `build_pod_esti_v05.py` / `upgrade_pod_esti_v05_*.py` | v05 生成与 FAR/并行升级 |
| `run_pod_scan_v05.py` / `pod_esti_v05_core.py` | 多进程 PoD 全量扫描 |
| `pod_esti_v05_cache_noise.npz` / `pod_esti_v05_cache_pod.npz` | **v05/v06 共用**主缓存 |
| `worklog_PoD_esti.md` | 工作日志 |

### 2.3 `PoD_esti_v06.ipynb` 模块要点（相对 v05）

| 模块 | 内容 |
|---|---|
| 0–4 | 与 v05 相同（参数/链路/引擎）；缓存名仍为 `pod_esti_v05_cache_*.npz` |
| ★口径 | markdown + 派生 `noise_ambient`/`bg`/`AXIS_KINDS`/`_axis_x` |
| 5 / 5b / 6 | 原图各再出一套 **bg 横轴**；savefig `pod_v06_*_{noise\|bg}.png` |
| 7 | 能量轴验证图，**不**双循环；图例标 bg 与折合 noise |
| 8 | `collect_critical` 按轴取 `_pod_x_from_rec`（bg=`r["noise"]`；noise 由 `e_lambda` 折合） |

### 2.4 当前参数（不得擅自改物理量）

```
NOISE_GRID[1] = 0.25→12 / 0.25
NOISE_GRID[4] = 0.25→40 / 0.25
N_MC_NOISE = 1e6
FAR_SPECS = [
  (10e-6,  "10ppm",  "10 ppm"),
  (100e-6, "100ppm", "100 ppm"),
  (0.001,  "0p1pct", "0.1%"),
  (0.005,  "0p5pct", "0.5%"),
  (0.01,   "1pct",   "1%"),
  (0.05,   "5pct",   "5%"),
]
D_TARGET = 15.0 m, RHO_TARGET = 0.10
CACHE_NOISE = pod_esti_v05_cache_noise.npz
CACHE_POD   = pod_esti_v05_cache_pod.npz
CACHE_NOISE_FALLBACK = [pod_esti_v04_cache_noise.npz]   # 允许，noise 与 FAR 无关
CACHE_POD_FALLBACK   = []                               # ★ 禁止读 v04 PoD
```

### 2.5 缓存策略（重要）

| 缓存 | 可否复用 v04 | 说明 |
|---|---|---|
| 纯噪声 `*_cache_noise.npz` | **可以** | peak 分布与 FAR 无关；v05 可读 v04 并写出 `pod_esti_v05_cache_noise.npz` |
| PoD `*_cache_pod.npz` | **不可以** | FAR 已从 2 档扩到 6 档；`CACHE_POD_FALLBACK=[]`；缓存键含 `far_tags` |

此前 v05 若曾用 v04 PoD fallback，那是旧行为；**当前代码已切断**。

### 2.6 多进程扫描工具链（PoD 全量重算的唯一正确姿势）

| 文件 | 作用 |
|---|---|
| `build_pod_core_v05.py` | 从 notebook 计算 cell（2/4/6/8/9/17/22/25）生成纯内核 |
| `pod_esti_v05_core.py` | 自动生成，勿手改；可被子进程 import，不含绘图 |
| `run_pod_scan_v05.py` | `ProcessPoolExecutor` 全量扫描，断点续跑 |

```powershell
$env:PYTHONIOENCODING="utf-8"
python run_pod_scan_v05.py --workers 20          # 全量
python run_pod_scan_v05.py --limit 8 --workers 8 # 冒烟
```

**状态：已跑完，208/208 档 × 六档 FAR，耗时 24.9 min，CPU 100%。**

### 2.7 推荐运行

```powershell
Remove-Item Env:MPLBACKEND -ErrorAction SilentlyContinue
$env:PYTHONIOENCODING="utf-8"
python run_notebook.py PoD_esti_v05.ipynb PoD_esti_v05_out.ipynb 20000
```

预期：模块 5 从 v05 noise 缓存秒级命中；模块 6 出六档阈值；模块 7/8 **应直接命中 `pod_esti_v05_cache_pod.npz`（208 档已跑完）**，不再重算。若该缓存被删或参数变更，改用 §2.6 的多进程脚本重跑，**不要**在 notebook 里干等线程版。

---

## 3. 目前面临的问题、卡在哪里

1. **`PoD_esti_v06.ipynb` 代码已就绪、语法通过，但尚未整本 execute 出图**。
   应在确认 `pod_esti_v05_cache_*.npz` 命中后，只重跑口径 cell + 作图 cell（5/5b/6/7/8），勿触发全量 MC。
2. 编辑器若仍开着旧 v05/v06 标签页，改完磁盘后必须 **Revert / 重开**，否则跑的是内存旧副本。
3. **v03 路线已废弃**。
4. 物理局限仍在：未加 DCR/串扰；`boost` 揉合能量/反射率/距离；PoD 只判探测到。

---

## 4. 下一步计划

1. 打开 `PoD_esti_v06.ipynb`：Restart Kernel → 跑到模块 5 命中 noise 缓存 → 跑口径 cell → 跑 5b/6/7/8 出双轴图。
2. 核对 N=1 两轴几乎重合、N=4 上 noise 轴约为 bg 轴的 1/4。
3. 固件阈值查找表：拟合 `T(noise)` / 或按需拟合 `T(bg)`（须写清口径）。
4. 可选：直接扫距离；加 DCR；扫描宏像元尺寸。

---

## 5. 踩过的坑（不要再踩）

### 工程 / 工具
1. **大 notebook 不能整文件硬读**：先 `json.load` 只导出 `cell['source']`。
2. **PowerShell 默认 GBK**：中文输出要 `PYTHONIOENCODING=utf-8`。
3. **`MPLBACKEND=Agg` 会阻止 notebook 内嵌图**：跑完要展示图时先清掉该环境变量。
4. **Glob 可能漏掉二进制 `.npz`**：不要仅凭 Glob 说“缓存不存在”；用 `Test-Path` / `np.load` 确认。
5. **长 MC 不能只在全部结束后存盘**：必须增量检查点 + 断点续跑；写盘用临时文件再 `os.replace`。
6. **线程数拉高必须减小 `MC_CHUNK`**：20 线程时用 5000；不要假设线性加速。
7. **FAR 标签禁止再用 `f"{far*1e6:.0f}"`**：对 5% 等会变成荒谬数字；一律用 `FAR_TAG` / `FAR_SPECS` 的 ASCII 标签（`100ppm`/`5pct`/…）。
8. **FAR 扩展后旧 PoD 缓存一律作废**：不要把 `pod_esti_v04_cache_pod.npz` 加回 fallback。
9. **`np.savez_compressed` 临时路径必须带 `.npz`**：写成 `path+".tmp"` 时 NumPy 会再追加成 `path+".tmp.npz"`，`os.replace` 在 Windows 上报 WinError 2；正确写法是 `path+".tmp.npz"`。
10. **PoD 内核受 GIL 限制，线程池无效**：`binary_macro_stepping` 是 Python 逐步循环 + 小数组 NumPy，线程再多 CPU 也只有 16–30%。**全量重算必须用多进程脚本 `run_pod_scan_v05.py --workers 20`**（实测 24.9 min，CPU 100%）；notebook 内的线程版只当缓存命中用。
11. **notebook 与多进程脚本不要同时跑 PoD**：两者写同一个 `pod_esti_v05_cache_pod.npz`，会互相覆盖并抢 CPU。
12. **改了 notebook 计算 cell（2/4/6/8/9/17/22/25）后必须重跑 `build_pod_core_v05.py`**，否则多进程脚本用的还是旧内核。
13. **用脚本改完 `.ipynb` 后，编辑器不会自动 reload**：Cursor/VS Code 的 notebook 编辑器有独立内存副本，此时运行的仍是旧代码（典型症状：磁盘已改成 `REF_FAR_TAG`，报错却还是 `tag = "100"`）。必须关闭标签页重开或 `File: Revert File`，再 Restart Kernel。**切勿在未 reload 时按 Ctrl+S，会把旧副本覆盖回磁盘。**
14. **噪声缓存 `res` 是 `{n_shots: ...}` 结构，顶层只有 2 个键**；真实档数看 `grid_key` 长度（208）。不要用 `len(res)` 判断它是否完整，PoD 缓存才是 `{(n_shots, noise): ...}`。

### 物理 / 统计
9. **边界假象**：生成网格左扩护带（`WARM_NS=50`），统计再掐头 24 ns。
10. **底噪用 `1-exp(-r_det·T_OVER)` / 更新过程平衡态**，不要用 `r_det·bin_width`。
11. **FAR 不能假设 bin 独立**：一律 MC；阈值比较用整数 `n_ge < target_far*n`。
12. **PoD 插值先保序**：`np.maximum.accumulate`，再找跨越点。
13. **等效距离方向**：`boost` 越小 ⇒ 距离越远。
14. **v03 平衡态起步**：实现不当会偏统计且更慢；v04/v05 继续用暖机快速引擎。

---

## 6. 给新会话的最短指令

1. 打开 **`PoD_esti_v06.ipynb`** 与 `worklog_PoD_esti.md`。
2. 确认 `pod_esti_v05_cache_noise.npz` / `pod_esti_v05_cache_pod.npz` 存在；命中后只出图。
3. 看图时分清 **noise（环境标准）** vs **bg（实测 baseline）**；不要把旧 `noise_mc` 当环境标准。
4. 全量重算 PoD 仍用 `run_pod_scan_v05.py --workers 20`（不要在 notebook 线程版里干等）。
5. 会话结束前更新本文件与 `worklog_PoD_esti.md`。

---

## 文档命名更正（2026-08-04）

- **错误做法**：交接文档写成 `handoff_现在工作.md`（把「现在工作」四个汉字当文件名）。
- **正确做法**：`handoff_<工作名>.md`，工作名与 notebook/工作标识一致。
  例：工作 `PoD_esti` → `handoff_PoD_esti.md`；日志仍是 `worklog_PoD_esti.md`。
- 用户已将文件重命名为 `handoff_PoD_esti.md`；全局规则 `session-handoff.mdc` 与 `CLAUDE.md` 已同步修改。
