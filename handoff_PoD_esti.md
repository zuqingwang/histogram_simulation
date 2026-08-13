# 交接文档 —— 当前工作：`PoD_esti`（探测概率估计）

> 文件名：`handoff_PoD_esti.md`（禁止 `handoff_现在工作.md`）。
> 最后更新：2026-08-11（模块 15 拆成三张独立图；此前 2026-08-10 **v30** 文件结构大改 + 全量重算）。
> 流水日志：`worklog_PoD_esti.md`；理论与公式推导：`theory_PoD_esti_v30.md`。
>
> **当前主产物是 `PoD_esti_v30.ipynb`（57 cell）。物理一条没改，只重排结构。**
> 模块 15 现出三张独立图：`pod_v30_m15_T_focus.png`（精简）、`pod_v30_m15_T_all.png`（全量）、
> `pod_v30_m15_qreq.png`（原右图）；精简版曲线由 `M15_T_FOCUS` 控制。
---

## ★ 新会话必读的四件事

1. **不要手改 `PoD_esti_v30.ipynb`。** 它由 `build_pod_esti_v30.py` 生成，手改会在下次
   build 时被覆盖。要改分析层就改 **`v30_cells.py`**；要改装配顺序或对 v20 内核的补丁，
   就改 **`build_pod_esti_v30.py`**。改完跑：
   `python build_pod_esti_v30.py; python build_pod_core_v30.py`

2. **横轴只用 `bg`。** `noise`（单发 `hist_i` 每 bin 均值 `= bg/N`）从 v30 起
   **不再作为任何图的横轴**，只用来描述单次直方图的底噪。

3. **阈值 `T` 是整数计数，而且必须是。** `hist_add` 是 27×N 条二值轨迹求和，
   判定 `peak >= T`，所以 `T=10.25` 与 `T=11` 完全等价。阈值曲线像阶梯是因为
   bg 连续步进 0.25 而 T 只能跳 1，**不是 bug**。v20 用插值造"连续阈值"的模块 12 已作废。

4. **这台机器的 PowerShell 不支持 `&&`**，串联命令用 `;`；带引号的 `python -c` 单行脚本
   基本转义不过去，需要多语句就写成 `.py` 文件。

5. **PoD 临界点求根器在 v30 被重写过**（v20 版本会把失败结果静默写进缓存，
   造成模块 7/12/13 曲线出现 3–5 倍毛刺）。动这块代码前先读「踩过的坑」第 27 条。
   改完必须跑 `python inspect_v30_cache.py` 看第 ③ 节的求根质量表，
   **「超容差」那一列必须是 0**。

---

## 0. 缩写表

| 缩写 | 英文全称 | 含义 |
|---|---|---|
| SPAD | Single-Photon Avalanche Diode | 单光子雪崩二极管 |
| FAR | False Alarm Rate | 虚警率 / 噪点率 |
| PoD | Probability of Detection | 探测概率 |
| MC | Monte Carlo | 蒙特卡洛 |
| bg | background（本项目口径） | `hist_add` 统计窗均值 |
| noise | 单次底 | 单次 `hist_i` 统计窗均值 |

---

## 1. 我们在做什么任务

仓库：`E:\claude temp\Histogram-simulation`。

主产物：**`PoD_esti_v30.ipynb`**（57 cell：code 37 / markdown 20）

在 1 bit（二值）SPAD 宏像元 + 多发累加的条件下，把"底噪多高 → 阈值定多少 →
要多强回波才能以给定概率检出 → 因此能测多远"这条链算通，并给出可复现的 MC 证据。

### v30 的文件角色分工（**这张表是新会话最需要的**）

| 文件 | 角色 |
|---|---|
| `v30_cells.py` | **分析层 cell 源码**（模块 5–15 的计算/绘图参数/绘图三段式）。要改图先改这里 |
| `build_pod_esti_v30.py` | **装配器**：从 v20 逐字抽物理内核 + 打 5 个带断言的补丁 + 拼上 `v30_cells` → notebook |
| `PoD_esti_v30.ipynb` | **生成物，不要手改** |
| `build_pod_core_v30.py` | 从 notebook cell 2/4/6/8/9/17/18/22 抽出 → `pod_esti_v30_core.py` |
| `pod_esti_v30_core.py` | **生成物**，供所有多进程脚本 import |
| `run_pod_v30_noise_scan.py` | 噪声扫描（模块 5 的数据源，也是模块 9/10/11 的数据源） |
| `run_pod_v30_pod_scan.py` | PoD 临界能量（模块 6 的数据源，也是模块 7/12/13 的） |
| `run_pod_v30_sig_scan.py` | 固定信号扫描（模块 8 的数据源，也是模块 14 的） |
| `compare_macro_v30.py` | 宏像元 3×9 vs 3×6（模块 15 的数据源） |
| `theory_PoD_esti_v30.md` | **全部理论与公式推导**（notebook 里不再有） |
| `inspect_v30_cache.py` | 四份缓存体检 + 阈值抽样表 |
| `check_v30_modules.py` | 按序无头实跑所有 code cell + matplotlib 缺字检查 |
| `show_nb_cells.py` | 打 cell 地图 |

### v30 模块表

每个分析模块都是 `[计算/载入缓存] → [绘图参数] → [绘图]` 三段式；绘图 cell 只读缓存。

| 模块 | 内容 | 数据来源 |
|---|---|---|
| 0–4 | 参数 / 光链路 / 时间窗 / SPAD 引擎 / 波形示例 | 现算 |
| 5 | 纯噪声 bg 扫描 → peak 分布 + 阈值（只叠 1%/5%）+ 密度条带 | **唯一大重算** |
| 6 | PoD 临界能量（只解 0.5/1/5/10% FAR） | PoD 扫描 |
| 7 | 全 bg 汇总：临界 peak / 能量 / 等效距离 | 复用模块 6 |
| 8 | 固定信号 × 全 bg：peak 均值 + 归一化增益 Δpeak/boost | 信号扫描 |
| 9 | 阈值倍数 ρ_N = T_N/T_1 | **复用模块 5** |
| 10 | hist 内 std / peak 均值 / peak std vs bg（叠二项与 Gumbel 解析） | **复用模块 5** |
| 11 | 有效 z 值 (T−μ_peak)/σ_peak | **复用模块 5** |
| 12 | PoD50/90 所需信号 | 复用模块 6 |
| 13 | 平方反比测远 | 复用模块 6 |
| 14 | 同信号不同 bg 的 peak 分布 + "peak+bg"加法理想线 | **复用模块 8** |
| 15 | 宏像元 3×9 vs 3×6（线性轴；含 3×6@N=6 等 n_tr 校验；**三张独立图** ①a 精简 / ①b 全量 / ② q_req） | 宏像元缓存 |
| 回收站 | v20 模块 6 图、v20 模块 12（连续阈值） | 整段 python 注释 |

### v20 → v30 模块编号对照（v30 重排过，找旧结论时看这张表）

| v20 模块 | v30 去向 |
|---|---|
| 5（noise 横轴扫描）+ 9（bg 横轴扫描） | **合并成 v30 模块 5**（统一 bg 横轴、0.25 步长、密度条带同 cell、只叠 1%/5% 阈值） |
| 6（六条 FAR 阈值大图） | **作废** → 回收站，python 注释保留 |
| 7（PoD 临界能量） | v30 模块 6（只解 4 条 FAR，交点图放大到过渡区，peak 分布改实心） |
| 8（全 bg 汇总） | v30 模块 7 |
| 10（固定信号 × bg 网格） | v30 模块 8（叠 1%/5% 阈值 + 新增「失守 bg」表） |
| 11（阈值倍数 ρ） | v30 模块 9（删掉 noise 横轴那张） |
| 13（有效 z 值） | v30 模块 11（删掉 noise 横轴那张，新增均值参考线） |
| 12（插值造连续阈值） | **作废**（口径错误，见坑 21） |
| 12C 那张图 | 保留为 v30 模块 12 |
| 14（测距） | v30 模块 13 |
| 15（同信号不同 bg 的 peak 分布） | v30 模块 14（含「peak+bg」理想线） |
| 16（宏像元对比） | v30 模块 15（改线性轴、600k MC、补 3×6@N=6） |
| 新增 | v30 模块 10（hist std / peak 均值 / peak std vs bg，复用模块 5 的 1e6 MC） |

### 口径（v30 收紧）

| 符号 | 含义 |
|---|---|
| `hist_i` | 第 i 发宏像元直方图 |
| `hist_add(N)` | 前 N 发之和；N∈{1,2,4} |
| **bg** | `hist_add` 统计窗每 bin 均值；`noise_target`/`noise_mc` 字段存的就是目标/实测 **bg** |
| noise | 单次 `hist_i` 统计窗均值 `= bg/N`。**不再作横轴** |
| **peak** | 在 `hist_add` 上统计的窗内最大计数 |
| **T** | 检测阈值，**整数计数** |
| `n_tr` | 轨迹数 `= n_pix × N_shots` |

### 与 v10 的关键差异

| | v10 | v11 |
|---|---|---|
| 扫轴 | `NOISE_GRID_AMB`，各 N 目标 bg=`AMB×N` | **`BG_GRID`** 对所有 N 相同 |
| bg 步长 | N=1→0.25；N=2→0.5；N=4→1.0 | **一律 0.25**（0.25→12，48 档） |
| 仿真 | 一次 AMB 仿 4 发前缀和 | 每档 `(N, bg)`：`noise_amb=bg/N` |
| 新分析 | 模块 9 | + **模块 10** 阈值倍数 ρ |

### 模块 10（新增）——要回答的问题

同 bg 下 $\rho_{N/1}=T_N/T_1$ 是否近似常数？倍数从哪来？

**已有解析答案**（完整推导见 `theory_peak_bg_multishot.md`）：

- 同 bg 下单 bin 分布是**精确二项** $\mathrm{Binomial}(27N,\ \mathrm{bg}/(27N))$：
  均值恒为 bg，方差 $\mathrm{bg}(1-\mathrm{bg}/(27N))$。
- bg 只锁均值不锁形状；N=1 单发速率 4 倍 → 二值饱和 + 死时间更强 → 欠离散 → 尾轻 → **阈值更低**。
- 大偏差修正项 $(T-\mathrm{bg})^2/(2\cdot 27N)$，$\propto 1/N$，是倍数的唯一来源。
- **ρ 不是常数**：@FAR=1% 从 1.033（bg=1）单调升到 1.195（bg=12）；@10 ppm 升到 1.286。
- peak 均值差 +2.3%（bg=1）→ +13.6%（bg=12）；peak std 比 1.05→1.52。
- N=1 在 bg=12、FAR=10 ppm 时 T=26，距二值硬上限 27 仅余 1。

---

## 2. 已经完成了什么

| 版本 | 文件 | 状态 |
|---|---|---|
| v05 | `PoD_esti_v05.ipynb` | 完整基线；保留 |
| v10 | `PoD_esti_v10.ipynb` | 统一 AMB 前缀和版；保留 |
| v11 | `PoD_esti_v11.ipynb` | 统一 bg 网格 + 模块 10；保留 |
| **v20** | **`PoD_esti_v20.ipynb`** | **当前主产物** |

### 关键文件

| 文件 | 作用 |
|---|---|
| `PoD_esti_v20.ipynb` | **主 notebook**（51 cell，模块 0–17） |
| `upgrade_pod_esti_v20_from_v11.py` | v11→v20 升级脚本（新模块源码都在这里，改 notebook 优先改它再重跑） |
| `build_pod_core_v20.py` / `pod_esti_v20_core.py` | 多进程内核（改过 notebook 计算 cell 后必须重新导出） |
| `run_pod_v20_noise_scan.py` | 噪声 ProcessPool；任务键 `(N,bg)` |
| `run_pod_v20_pod_scan.py` | PoD ProcessPool |
| `run_pod_v20_sig_scan.py` | **v20 新增**：模块 9.3/15 的信号扫描，1296 任务 ProcessPool，11.7 min |
| `check_v20_modules.py` | **v20 新增**：不用 Jupyter 做语法自检（`--syntax`）与无头真跑（`--all`） |
| `pod_esti_v20_cache_*.npz` | 主缓存；noise/pod 有 **v11 fallback**，读到即同步写回，不会重算 |
| `PoD_esti_v11.ipynb` 及 `*_v11_*` 一整套 | 上一版，保留可跑 |
| `theory_peak_bg_multishot.md` | **阈值倍数解析模型**（推导 + 定量表 + 对照清单） |
| `theory_peak_bg_multishot.py` | 上文的数值脚本，产出 `theory_peak_bg_multishot_fig.png` |
| `check_same_bg_two_ways.py` | 定向 MC：同 bg=4 的 A(noise=1×4) vs B(noise=4×1) |
| `check_bin_correlation.py` | bin 间 ACF 与有效独立 bin 数 M_eff 诊断 |
| `scan_hist_std_peak.py` | 三联图扫描（hist内std / peak均值 / peak std），带缓存+多进程 |
| `compare_macro_3x9_vs_3x6.py` | **宏像元 3×9 vs 3×6 阈值对比**，带缓存+多进程；已全量跑完 |
| `compare_macro_3x9_vs_3x6_cache.npz` | 上文缓存（24 档 p_eq × 7 配置 × 200,000 MC）；键数组 `amb` = `27·p_eq`，**改命名时不要动它的数值** |
| `compare_macro_3x9_vs_3x6.png` / `_log.txt` | 6 联图 / 完整数值表（日志是**无 BOM 的 UTF-8**，用 UTF-8 打开） |
| `check_engine_vs_v45.py` | **引擎一致性核对**：与 `lidar_histogram_sim_v45.ipynb` 逐行 + 比特级 + 统计级 |
| `theory_engine_equivalence.md` | **引擎一致性的理论模型**：泊松稀释→更新过程→逆变换→并集恒等式→更新-回报→n_tr 折叠→O(dt²) 离散化 |
| `theory_engine_equivalence.py` | 上文的逐环节数值检验（T1–T5b），日志 `theory_engine_equivalence_log.txt` |

### 参数（物理量未改）

```
N_SHOTS_LIST=[1,2,4]
BG_GRID = 0.25→12 / 0.25（48）；NOISE_GRID[n]=BG_GRID
仿真：noise_amb = bg / N
N_MC_NOISE=1e6；N_WORKERS=20；MC_CHUNK=5000
★ FAR_SPECS = 7 条：10ppm/100ppm/0.1%/0.5%/1%/5%/10%（10% 为 v30 新增）
★ POD_FARS  = [0.5%, 1%, 5%, 10%]  只对这四条求 PoD 临界能量，阈值仍七条全算
★ 模块 8/14：BOOST=[0,0.004,…,0.032]（9 档），N_MC_SIG=20_000（v20 为 8000）
★ 模块 15：MACRO_N_MC=600_000（v20 为 300k）；3×6 额外跑 N=6 做等 n_tr 校验
★ RUN_ENGINE_CHECK=False（模块 3c 默认跳过）
CACHE_* = pod_esti_v30_cache_*.npz + compare_macro_v30_cache.npz
★ CACHE_*_FALLBACK 一律为空：FAR 列表与 res 结构都变了，旧缓存不可复用
```

### 缓存现状（v30 全量重算，2026-08-10）

| 缓存 | 规模 | 状态 |
|---|---|---|
| `pod_esti_v30_cache_noise.npz` | 48 bg × N∈{1,2,4} × 1,000,000 MC | 144/144 done，**18 min** |
| `pod_esti_v30_cache_pod.npz` | 144 档 PoD 临界点（4 条 FAR × PoD50/90） | done，**11.4 min**（用修好的求根器重跑） |
| `pod_esti_v30_cache_signal.npz` | 48 bg × 9 boost × 3 N × 20,000 MC | 1296/1296 done，24.9 min |
| `compare_macro_v30_cache.npz` | 24 档 p_eq × 7 配置 × 600,000 MC | 168/168 done，24.4 min |

总机时约 **79 min**。`scan_hist_std_peak_cache.npz` **已不再使用**：`hist_std` 并入噪声扫描的
充分统计量（`stats_from_hist_i` 新增 `hist_std_sum`），模块 10 直接复用那批 1e6 MC。

### 验证现状

- `python build_pod_esti_v30.py` 自带 AST 自检：37 个 code cell **0 语法错误**。
- `python check_v30_modules.py` → 按序无头实跑全部 code cell：**37/37 通过，0.3 min，无缺字警告**。
- `python inspect_v30_cache.py` → 四份缓存完成度、阈值抽样表，以及
  **PoD 求根质量体检**（逐 N × FAR × PoD 等级统计 `|验证 PoD − 目标|`）。
  当前结果：**1152 个临界点全部收敛，超容差 0 个**，最大误差 0.0198、中位约 0.005。
- （**还没在 Jupyter 里 Restart & Run All**，只做过无头实跑。）

### 已完成的旁支分析：宏像元 3×9 vs 3×6（`compare_macro_3x9_vs_3x6.py`）

口径（**必须先说清，否则结论会反过来**）：环境光与信号**都按每 SPAD 均匀**处理。
于是 `bg = n_tr·p_eq`，宏像元收到的信号 ∝ `n_tr = n_pix·N`，
灵敏度判据 `q_req = (T − bg)/n_tr`（每 SPAD 每发需额外贡献的点亮概率，越小越灵敏）。

**环境光的唯一自变量是 `p_eq`** —— 单个 SPAD、单发、单个 1 ns bin 被点亮的平衡态概率，
与宏像元多大、累加多少发都无关。脚本里扫 24 档 `p_eq = 0.0185 → 0.4444`，
对应照度 38 → 1326 klux。（代码用 `N_PIX_REF·p_eq = 0.5…12` 打整齐网格，`N_PIX_REF=27`，
这个参考量只是"给档位起个整数名字的标尺"，不是任何配置的物理量。
早期版本把它叫 `noise27`，已弃用；老日志里出现 `noise27=x` 就是 `p_eq = x/27`。）

**两种横轴口径，回答的是不同问题，引用比值时必须写清是哪一种：**

| 口径 | 图 | 含义 |
|---|---|---|
| **同一片天光**（横轴 = 照度 klux） | 图① | 各配置 bg 不同，bg ∝ `n_tr` |
| **同一个 bg**（横轴 = 各配置自身 bg，log 轴） | 图②–⑥ | 同一 bg 意味着各配置处在**不同**天光下 |

1. **纯噪声阈值只取决于 `n_tr = n_pix × N_shots`**。引擎把「SPAD 数」和「shot 数」
   折进同一个轨迹数维度（见 `noise_macro_hist_fast` docstring）。
   MC 实证：3×6@N=6 与 3×9@N=4（都是 n_tr=108）在 24 档上 **T@1% 最大差 0 计数**；
   图②–⑤ 两条线完全重合、图⑥ 紫线恒为 1.000。
2. **同一片天光**下 3×6@N=2（n_tr=36）vs 3×9@N=4（n_tr=108），FAR=1%：绝对阈值前者低一半以上
   （T=15 vs 32 @ p_eq=0.148，因为底噪只有 1/3），但 `q_req` 前者要差 **1.8–2.0 倍**。
3. **同一个 bg** 下同样这一对：差距放大到 **2.3–2.7 倍**
   （bg=4 时 2.70，bg=16 时 2.25）。因为把 3×6@N=2 拉到同样的 bg，等于让它处在 3 倍强的天光里。
4. 同 N 下 3×6 相对 3×9（同天光）：bg 恒为 2/3，T 低 1–19 计数，`q_req` 高 **1.19–1.37 倍**
   （围绕 √(27/18)=1.22，抖动来自整数阈值量化）。
5. 图② 的新发现：**`T` 对自身 `bg` 作图时，7 种配置几乎重合成一条曲线** ——
   给定底噪，阈值基本只由 bg 决定，对 `n_tr` 只有弱依赖。
   图③ 的 `T − bg` 在高 bg 处饱和（n_tr=108 → ≈20，n_tr=36 → ≈12），
   因为 `T−bg ≈ z·√(n_tr·p_eq(1−p_eq))` 而 `p_eq→0.44` 时 `p(1−p)` 已接近极大值 0.25。
6. 设计取舍：**宏像元缩小 1.5 倍，用 1.5 倍发数可精确换回同样噪声性能**，代价是帧率。

### 已完成的核对：引擎与 `lidar_histogram_sim_v45.ipynb` 一致（`check_engine_vs_v45.py`）

结论：**PoD 没有把 SPAD 简化成 8 ns 硬死时间，与 v45 是同一套引擎。**

- `pod_esti_v11_core.spad_binary_trace` 是 v45 **cell 32（模块 9b）**同名函数默认路径的
  逐行移植；归一化后仅存的差异是内联临时变量 `p_fire` 与两句合并成一行。
- 同一 `default_rng(seed)` 下 60/60 条轨迹逐 bin 完全相同。
- 扫描实际用的快速引擎 `noise_macro_hist_fast` 是精确逐光子引擎的连续时间极限，
  p_bin 与解析值一致到 **≤0.4%**（noise27=0.5/2/6/12 四档）。
- RC 恢复是显式建模的：`vov_frac = 1 − exp(−Δt/τ_RC)`，
  触发概率 `PDE_max·g(vov_frac)`（`g` 为 exp 型凹函数、k=3），每次雪崩把 Vov 打回 0。

**v45 内部并存两套引擎**，是读出方式不同，不是新旧版本：

| | v45 模块 7b `simulate_spad_shot_rc` | v45 模块 9b / PoD `spad_binary_trace` |
|---|---|---|
| 读出 | timestamp 计数（多 bit） | 每 1 ns 采样点 0/1（1 bit） |
| Vth_frac 用法 | 雪崩时判 `vov_frac ≥ 0.60` 才计数 | 只以 `T_OVER = −τ·ln(1−Vth) = 8.00 ns` 进入 |
| 复位策略 | `reset_mode='count'`：亚阈雪崩**不复位** | 每次雪崩都复位，过阈窗**顺延堆积** |

PoD_esti 研究的就是 1 bit 读出，用模块 9b 是对的。

### 为什么三个引擎**必然**一致（理论，见 `theory_engine_equivalence.md`）

四份代码只对应一个数学对象：`spad_binary_trace`（逐光子）与 `binary_macro_stepping`（步进）
是同一个离散模型的两种写法，`noise_macro_hist_fast`（快速）是连续模型的精确采样器，
`p_bin_equilibrium` 是同一连续模型的平衡态闭式解。链条：

1. 环境光子是齐次泊松、每个以 `PDE·g(vov)` 独立触发 ⟹ **泊松稀释**给出条件强度
   `h(Δ) = r_det·g(1−e^{−Δ/τ})`，**只依赖年龄 Δ**（雪崩把 Vov 打回固定的 0）。
2. 年龄依赖强度 ⟹ **更新过程**，`S(Δ)=e^{−H(Δ)}`。
3. `Δ = H⁻¹(E), E~Exp(1)` 的分布精确等于 S ⟹ 快速引擎采的就是这个律。
4. `∪[aₖ,aₖ+T) = ⊔[aₖ,min(aₖ+T,aₖ₊₁))` ⟹ 差分数组 + cumsum 是**恒等变形**。
5. **更新-回报定理**：`p_bin = (1/μ)∫₀^{T_OVER}S = E[min(X,T_OVER)]/E[X]`。
6. 纯噪声下 n_pix×N 条轨迹 i.i.d. ⟹ 单 bin 边缘精确为 `Binomial(n_tr, p_bin)`，
   联合律只依赖 `n_tr` ⟹ **这就是「3×6@N=6 ≡ 3×9@N=4」的证明**。

**离散化误差是 O(dt²)**（把离散模型也精确解出来算的，不是 MC）：dt 每减半误差降 4.00 倍；
生产的 `DT_FINE=200 ps` 处 p_bin 偏差仅 1e-4 相对量级，低于解析式自身的数值积分精度。
原因是 `h(0)=r_det·g(0)=0`——刚雪崩完触发率恰为 0，把"一步内两次雪崩"这个主误差源压掉了。

抖动的影响解析上约 `f_X(T_OVER)·σ²`，相对量级 **3.3e-5**，可忽略。

**推导覆盖不到的地方**：有信号时条件强度依赖绝对时间 `t`，不再是更新过程，
上述 2–5 全部失效、没有闭式，只能用步进引擎——这不是偷懒，是数学上不允许。

---

### v20 新模块的定量结论（都已实测，不是预测）

**模块 11（要求 1 后半）**：单条 hist 内的 std 贴合二项解析 √(bg(1−bg/27N))，
明显低于纯泊松 √bg —— 1 bit SPAD 的亚泊松压缩。bg=10.75 时 N=1 实测 2.480（解析 2.544）、
N=4 实测 3.018（解析 3.111）；实测系统性略低 2–3%，来自 bin 间正相关使样本 std 下偏。
peak 标准差随 bg 先升后平（N=1 在 bg≈9 后掉头，p_eq 越过 0.5 附近方差被压）。
图 ③ 另叠 **Gumbel 极值解析** σ_peak≈(π/√6)σ_bin/z_M（z_M=(peakμ−bg)/σ_bin，无拟合参数）：
大 bg（≳6）实测与解析贴合到几个百分点（bg=10.75 N=1 实测 1.315 / 解析 1.346），
小 bg 计数少、分布离散强右偏、未进极值渐近区，解析系统性偏低（bg=0.25 实测 0.662 / 解析 0.227）。

**模块 12（要求 1「不要阶梯」）**：用户澄清"不要阶梯"是**画法** —— 阈值采样点要用**直线相连成折线**，
不要用 `step` 阶梯函数（带竖直立边）。图 A 整数阈值已改成 `plot(bg, T, "-", marker=".")` 折线，
连续阈值 $T_c$ 退成细虚线趋势线叠加。另：`far_threshold_continuous(cnt, far)` 在
log(生存函数) 上线性插值求实数阈值 $T_c$，满足 $T_{整数}=\lceil T_c\rceil$，与原整数阈值逐档一致。
1e6 MC 下 10 ppm 那条在小 bg 处有 6/48（N=1）、3/48（N=2）档落到 MC 分辨极限，记 NaN 不外推。

**模块 13（要求 2）**：FAR=1%、PoD90、bg=10.75 时临界发射能量
N=1 需 **175.6 nJ**、N=2 需 22.57 nJ、N=4 只需 **8.42 nJ**（相差 20 倍）。
同条件净峰高 `S_net` 分别 11.71 / 16.55 / 19.49 计数。
**注意 `S_net = peak_mean − bg` 会高估信号贡献**：即使没有信号，
在 15 个 bin 的信号窗里取最大值本身就已高于 bg 约 2σ。
干净的口径是 `ΔS = peak_mean − peak_mean(boost=0)`，模块 13 会自动引用模块 9.3 的基线算出来。

**模块 14（要求 3）**：α=0.1 /km 时大气衰减修正 **< 3%**，本仓库距离段内平方反比就是主要规律。
FAR=1%、PoD90、bg=6.25 时纯 1/D² 测距 N=1 为 56.4 m、N=2 为 112.3 m、N=4 为 168.0 m。
曲线逐点有几米抖动，来自整数阈值阶梯 + PoD 临界点求解的 MC 噪声（`POD_VERIFY_TOL=0.02`），
**看趋势不要抠单点**。

**模块 15（要求 4）—— 三条明确结论**：

1. **peak 均值不是简单地加上 bg。** boost=0.016 时 N=1 的 Δμ 从 bg=0.25 的 6.77
   掉到 bg=10.25 的 **2.69**（掉 60%）；N=4 从 27.82 掉到 22.45（掉 19%）。
   主因用"**抢占**"模型解释：1 bit SPAD 一个 bin 只能亮一次，环境光先点亮就轮不到信号，
   净增量 ∝ (1−p_eq)，于是 Δμ(bg)/Δμ(bg_min) ≈ (1−p_eq(bg))/(1−p_eq(bg_min))
   （图 15B 下排的黑虚线，无拟合参数）。实测比它掉得更快
   （N=1 在 bg=10.25 处实测 0.397 vs 模型 0.626），多出来的来自极值竞争。
2. **peak 标准差加信号后变大，不是变小。** σ(b)/σ(0) 从低 bg 的 2.8–5.3 降到高 bg 的
   1.04–1.64，但**始终 ≥ 1**。
3. **分布形状确实在变。** 纯噪声 peak 右偏（偏度 0.2–1.0，极值分布特征），
   加信号后偏度掉到 **≈0.05**、接近对称的二项形状。
   所以只报"均值 + std"不足以描述 peak，必须看完整分布（图 15A）。

---

## 3. 目前面临的问题、卡在哪里

**没有卡住的地方。** v30 已生成，四份缓存已全量重算。剩下的都是可选收尾：

1. **还没在 Jupyter 里 Restart & Run All**。缓存都在，只花在绘图渲染上，不会重算 MC。
   编辑器须先 **Revert/重开** `PoD_esti_v30.ipynb` 再 Restart → Run All。
2. ~~10% FAR 档的 PoD 交点偏噪~~ **已解决**。曾以为是「整数阈值在低计数区的固有现象」，
   实际是 PoD 求根器的 bug（初值偏 + 验证只做一次固定修正 + 无条件接受），
   详见下面「踩过的坑」第 6 条。修好后 1152 个临界点全部收敛，模块 7/12/13 的毛刺消失。
3. **模块 10 面板③的 Gumbel 极值解析在高 bg 段系统性偏低**（实测 σ_peak > 解析）。
   推测是 `T_OVER=8 ns` 造成的 bin 间正相关让有效独立 bin 数 `M_eff` 明显小于 152，
   而解析式仍按 152 个独立 bin 算。这是**模型不够，不是数据错**，待补 `M_eff` 修正项。
4. **（待查，v11 起就挂着）** `theory_engine_equivalence.py` 的 T5b 里，步进引擎
   `binary_macro_stepping` 在 dt=800/200 ps 两档都比解析值**偏低**（2.2σ / 0.7σ，
   未达显著但方向一致），怀疑是它「先出 bin、再处理本步雪崩」的半步对齐。
   纯噪声扫描用快速引擎不受影响，但 **PoD 信号支路
   （`binary_macro_stepping_per_shot`）用的正是它**，建议加大样本复核。

---

## 4. 下一步计划

1. 打开 `PoD_esti_v30.ipynb`，**Revert/重开 → Restart → Run All**，确认渲染无误。
2. 逐图核对绘图参数是否合适：每个模块的「绘图参数」cell 里都注明了每个参数管哪根轴，
   看着不顺眼直接改那个 cell 再跑绘图 cell，**不会触发重算**。
3. 用模块 9/10/11 的 MC 曲线逐条对照 `theory_PoD_esti_v30.md` 与
   `theory_peak_bg_multishot.md` 第 10 节的检查清单，把 ρ̄、残差写回 `worklog_PoD_esti.md`。
4. 复核 `binary_macro_stepping` 的半步对齐（见上第 4 条）。
5. 给模块 10 面板③的 Gumbel 解析补 `M_eff` 修正（见上第 3 条）：
   先用噪声扫描里的 peak 分布反解「等效独立 bin 数」，看它随 bg 怎么走。
6. （可选）把「亚阈雪崩不复位」的非顺延变体做成开关，验证同 bg 下阈值曲线的二阶差异；
   解析上 p_bin 差异从 p_eq=0.0185 的 −0.4% 增大到 p_eq=0.4444 的 −13.7%，
   但因 `r_det` 是由目标 bg 反解的，这个差异主要被吸收进 klux 换算里。

### 全部重跑命令（正常不需要，缓存都在）

```powershell
$env:PYTHONIOENCODING="utf-8"
# 改过 v30_cells.py / build_pod_esti_v30.py 之后必须先重建
python build_pod_esti_v30.py                                   # → PoD_esti_v30.ipynb（含 AST 自检）
python build_pod_core_v30.py                                   # → pod_esti_v30_core.py

# 四个扫描必须按顺序（PoD 依赖噪声扫描的阈值）
python run_pod_v30_noise_scan.py --workers 20                  # 144 档 × 1e6 MC，约 18 min
python run_pod_v30_pod_scan.py   --workers 20                  # 144 档，约 8 min
python run_pod_v30_sig_scan.py   --workers 20 --n-mc 20000     # 1296 档
python compare_macro_v30.py      --workers 20 --n-mc 600000    # 7 配置 × 24 档

python inspect_v30_cache.py                                    # 缓存体检
python check_v30_modules.py                                    # 无头实跑 + 缺字检查
```

**注意**：这台机器的 PowerShell **不支持 `&&`**，串联用 `;`。

---

## 5. 踩过的坑（不要再踩）

1. **「新增三部分」≠删掉 v05**（v10 曾做错）。  
2. **v10 的 `AMB×N` 使 N=2/4 的 bg 步长变粗** → 同 bg 对比不公平；v11 改为统一 `BG_GRID`。  
3. **禁止读 v10/v05 缓存**（网格口径已变）。  
4. `r_det_for_noise(noise, n_tr)` 反解单次 noise 时 **`n_tr=27`**。  
5. 噪声/PoD 全量必须 **ProcessPool**（GIL）。  
6. 冒烟 `--n-mc 2000` 若写入主缓存会污染正式 1e6 结果 → 冒烟后须删 `pod_esti_v11_cache_noise*.npz`。  
7. `subprocess.run` 缓冲 → 用 `_run_cmd_stream`。  
8. **不要直接画整数阈值的比值 `T_N/T_1`**：T≈10 时 1 个计数就是 10%，锯齿会淹没真实趋势。
   看 `T₄−T₁`（整数差），或用连续阈值插值（`theory_peak_bg_multishot.py` 的 `thr_continuous()`）。  
9. 高斯闭式 `T≈bg+zσ` 在低 bg 严重偏低（bg=1、FAR=1%：4.75 vs 精确 7），只能看趋势不能取数。  
10. **不要用「各发峰不对齐」推断多发 peak 更小**。运算是 `max_j Σ_i h_i[j]`，
    不是 `Σ_i max_j h_i[j]`。同 bg 下胜负只看求和后的方差：本仿真 SPAD 是 1 bit（欠离散），
    N 大方差大 → peak 大、阈值高。40 万次 MC 实测：bg=4 时 noise=1×4 的 peak=9.11 >
    noise=4×1 的 8.71。详见 `theory_peak_bg_multishot.md` 3.4 节。  
11. **bin 之间是【正】相关，不是负相关**（曾写反过）。引擎判据「点亮 ⟺ 最近雪崩距今 <
    T_OVER=8 ns」使一次雪崩点亮约 8 个连续 bin，ACF≈max(0,1−L/8)。
    后果：单 bin 边缘分布仍精确为 Binomial(27N, p_eq)，但**有效独立 bin 数 M_eff≈46–76**，
    远小于名义 152。凡是用「152 bin 独立」估 peak 的地方都会偏高约 0.6–1.2 计数。
12. **比较不同宏像元尺寸前，先讲清"信号是均匀还是按像斑加权"**。
    按像斑 `FX=[0.0014,0.0152,0.084,0.234,0.330,…]` 加权时，9 列砍到 6 列只丢 1.8% 信号
    → 结论是 3×6 更灵敏；按均匀处理时 3×6 只收到 2/3 信号 → 结论反转成 3×6 差 1.2 倍。
    **两个结论都对，只是口径不同**，不写清口径的对比图是误导。当前采用**均匀**口径。
13. **不要用 `python xxx.py | Tee-Object` 跑长扫描**。中途中断 shell 会杀掉管道消费端，
    python 父进程卡死在 stdout 写入上：worker 继续烧满 CPU，但检查点和日志都停住不动，
    看上去像"跑得很慢"，其实永远不会结束。要留日志就用 `*> 日志.txt` 直接重定向。
13b. **PowerShell 落中文日志会写坏编码**。`python … *> log.txt` 按 ANSI 落盘，
    `… | Out-File -Encoding utf8` 在本机实际写成 UTF-16，两种都会让日志变乱码。
    **正确做法**：`cmd /c "set PYTHONIOENCODING=utf-8&& python … > log.txt 2>&1"`，
    由 cmd 原样透传字节。（`*>` 只在纯 ASCII 输出时才安全。）
14. **不要以为 PoD 把 SPAD 简化成了 8 ns 硬死时间**。已逐行 + 比特级核对过，
    与 v45 模块 9b 完全一致，RC 恢复与恢复期部分灵敏都在模型里（见第 2 节末）。
    v45 模块 7b 那套 `reset_mode='count'` 的计数引擎是**另一种读出方式**，不是"正确版本"。
15. **新版换缓存名时，内核里的加载顺序也要同步改。** `build_pod_core_v20.py` 的 FOOTER
    一开始只试主缓存名 `pod_esti_v20_cache_noise.npz`，结果 `pod_esti_v20_core.py` 里
    `NOISE_RES={}`、`THRESH={}`，多进程脚本会**从零重算 144 档 × 1e6 MC**。
    已改成按 `[CACHE_NOISE, *CACHE_NOISE_FALLBACK, CACHE_NOISE_CKPT]` 顺序找，命中旧版后同步写回。
16. **PowerShell 不支持 heredoc（`<<EOF`），也不支持 `&&` / `&` 串命令。**
    写多行 git commit message 要先把正文 `Write` 成临时文件再 `git commit -F 文件`；
    串命令用 `cmd /c "a && b"` 或 PowerShell 的 `;`。
17. **nbformat 4.5 要求每个 cell 有唯一 `id`。** 用脚本造新 cell 时忘了加会让 Jupyter 报
    schema 警告。`upgrade_pod_esti_v20_from_v11.py` 里用 `_next_id()` 统一发号。
18. **`S_net = peak_mean − bg` 会高估信号贡献。** 即使 boost=0，在 15 个 bin 的信号窗里
    取最大值本身也已高于 bg 约 2σ。要看"信号净带来多少"，用
    `ΔS = peak_mean − peak_mean(boost=0, 同 bg 同 N 同窗口)`（模块 13 已自动这么算）。
19. **不要凭直觉断言"信号会把 peak 钉死、σ 变小"。** 我在写模块 15 说明时先这么猜过，
    8,000 MC × 1296 档的实测否掉了它：σ(b)/σ(0) **始终 ≥ 1**。信号自己带来的二项涨落
    叠在噪声之上，并没有减小 peak 的离散度。
19b. **v11 的模块 9.3 会把「跑了一半的 `.partial.npz`」当成完整缓存收下**。
    `_try_load_sig_m9` 只校验 `n_mc` / `grid_key` / `boosts`，**没校验完成度**，
    于是断点续跑留下的 `pod_esti_v11_cache_signal.partial.npz`（实测只有 72/432 档有数据）
    会被静默采纳，剩下的档位全是零。已删掉那个残留文件；
    **v20 的 `_try_load_sig_m9` 已加 `done_<N>` 完成度校验**，不会再犯。
    教训：缓存的校验键里必须包含「完成度」，不能只包含「参数是否匹配」。
20. **notebook cell 之间有隐式依赖**：模块 13/14 用模块 11 定义的 `_COLORS_N`，
    模块 14 用 cell 28 的 `equiv_distance`，模块 15 用 cell 35 的 `SIG_M9`。
    Run All 顺序没问题，但单独跑某个 cell 会 NameError。
    `check_v20_modules.py` 的默认 cell 列表已按依赖排好序，验证时直接用 `--all`。
    **v30 已改善**：`_COLORS_N` 移到模块 0 全局参数里定义，不再依赖某个分析模块先跑过。

### v30 新增的坑

21. **【最重要】"阈值应该是 0.25 的小数"是误解，不要照做。**
    `hist_add` 是 27×N 条**二值**轨迹求和，取值只能是整数 0…n_tr，判定是 `peak >= T`，
    所以 `T=10.25` 与 `T=11` 是**同一个判定**，0.25 粒度的阈值在计数轴上没有物理意义。
    阈值曲线像阶梯的真正原因：**bg 是均值、可以 0.25 连续步进，而 T 只能按 1 计数跳**。
    **v20 模块 12 用插值造"连续实数阈值"是错误方向，v30 已整块作废进回收站。**
    唯一能出现 0.25 粒度的真实口径是除以 `N_shots` 看"每发平均计数"（粒度 1/N，N=4 → 0.25），
    但它与 bg 不同轴，混画会串口径。

22. **进度打印里不要写死 FAR tag。** `run_pod_v30_pod_scan.py` 的 `_progress_msg()`
    原样从 v20 抄过来，写死了 `critical["100ppm"]`。v30 把 `T_map` 限制到 `POD_FARS` 之后，
    `critical` 里根本没有 `100ppm` 键 → **每一档都误报「无有效交点」**，看着像全盘失败，
    其实数据完全正常。已改成从 `POD_FARS` 自动挑第一个存在的 tag。
    **教训：日志里引用配置项要跟着配置走，否则会制造假警报，浪费排查时间。**

23. **`PoD_esti_v30.ipynb` 是生成物，手改会丢。** 改分析层去 `v30_cells.py`，
    改装配/补丁去 `build_pod_esti_v30.py`，然后重新 build。
    build 脚本里的补丁都用 `must_replace()` 带命中数断言，v20 结构一变就直接报错，
    不会悄悄产出错的 v30。

24. **模块 5 的计算段必须拆成两个 cell**（noise 扫描 / 阈值）。合成一个 cell 的话，
    `build_pod_core_v30.py` 没法在"自动开跑"处截断——截了会把后面的阈值函数一起切掉。

25. **`σ_peak` 的方向别写反**（这是坑 19 的再强调，我在写 v30 理论文档时又写错了一次）。
    有信号时 `σ(b)/σ(0)` **始终 ≥ 1**，不是变小。信号自身的二项涨落叠在噪声之上，
    盖过了"位置被钉住"省下的极值离散度。

26. ~~**10% FAR 是 v30 新加的档，天生噪。** T 只有 3–6 个计数时 PoD 随 boost 是粗台阶，
    probit 拟合吃力，个别档验证 PoD 落到 0.27 / 0.999。不是 bug，是整数阈值在低计数区的固有现象。~~
    **【这条判断是错的，见第 27 条】** 把真 bug 归因成"固有现象"，差点让它留在成品里。
    教训：当"固有现象"这个解释能让你不用动手修时，先怀疑这个解释本身。

27. **【最重要】PoD 临界点求根器会把失败结果静默写进缓存 —— 已修，但要知道它长什么样。**

    **现象**：模块 7 的临界 peak / 临界发射能量 / 等效测距三张图上，5% 与 10% FAR 曲线
    出现大量 3–5 倍的尖刺，相邻 bg 档反复横跳。例：`N=4, bg=8.75, 5%` 解出 `boost=0.0287`，
    而邻档 `bg=9.00` 只要 `0.00779`。

    **怎么发现的**：不是看图猜，是去缓存里把每个临界点的 **验证 PoD** 调出来对比目标值。
    发现失败点的验证 PoD 是 **1.000 或 0.68**，根本不是 0.90 ——
    **错误一直明写在数据里，只是从来没人检查过。**

    **根因三层，缺一不可**：
    - 粗扫描太稀：`N_POD_COARSE=11` 个点铺满 4 个数量级（0.4 decade 间距），每点仅 300 MC，
      粗交点本身就能偏半个数量级；
    - 局部窗太窄：`POD_LOCAL_HALF_DECADE=0.22`，粗交点一偏，加密网格**整段错过真根**；
    - 验证不设防：`_verify_critical_batch` 只做 **一次** Newton 步、步长还夹在 ±0.25 decade
      （初值偏 0.6 decade 时数学上就追不回来），然后**无条件接受**。

    **修法**（`build_pod_esti_v30.py` 的 `patch_pod_scan` + `_SOLVER_NEW`）：
    粗扫 11→15 点、300→600 MC；局部窗 0.22→0.35 decade；probit 改**逐 level 局部拟合**
    （只用交点 ±0.6 decade 内的点，不让远处 PoD≈0/≈1 的饱和点拽偏斜率）；
    验证改成**带括号的多轮迭代**（有括号走 probit 割线、退化则二分；无括号用斜率做 Newton 步
    并主动外扩），最多 6 轮，每轮所有活跃候选一起并行评估，最后取最接近目标的点，
    并把 `pod_err` / `verify_rounds` 存进记录。

    **结果**：1152 个临界点全部收敛，超容差 0 个，最大误差 0.0198。模块 7/12/13 毛刺消失。

    **通用教训（下次写任何迭代求根都适用）**：
    ① 解完一定要回头验一遍，不合格就继续迭代或报错，**绝不能静默接受**；
    ② 固定次数 + 夹死步长的"修正"不是求根，是许愿；
    ③ 把解的质量指标（这里是 `pod_err`）**存进缓存**，这样体检脚本能直接查，
       不用靠肉眼看曲线毛刺来发现问题——`inspect_v30_cache.py` 第 ③ 节现在就是干这个的。
