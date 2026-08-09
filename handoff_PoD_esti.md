# 交接文档 —— 当前工作：`PoD_esti`（探测概率估计）

> 文件名：`handoff_PoD_esti.md`（禁止 `handoff_现在工作.md`）。
> 最后更新：2026-08-09（**v20**：在 v11 上追加模块 11–17，把用户的 4 条要求补齐；
> 三份全量扫描缓存已全部跑完；无头验证 0 错误）。
> 流水日志：`worklog_PoD_esti.md`。
>
> **当前主产物是 `PoD_esti_v20.ipynb`。** v11 的内容一条都没删，物理参数一个都没改。
> 剩下的唯一主线动作是在 Jupyter 里 Restart & Run All 跑一遍确认渲染（缓存都在，不会重算）。

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

主产物：**`PoD_esti_v20.ipynb`**（51 cell）
- 由 `PoD_esti_v11.ipynb` 经 `upgrade_pod_esti_v20_from_v11.py` 升级。
- v11（38 cell）= v10 全模块 + 统一 bg 网格 + 模块 10；**v20 在其后追加模块 11–17，一条不删。**

### v20 追加的模块（对应用户的 4 条要求）

| 模块 | cell | 内容 | 要求 |
|---|---|---|---|
| 11 | 38/39 | 每条 hist 内 152 bin 的 std 均值 / peak 均值 / peak 标准差 随 bg | 1 后半 |
| 12 | 40/41 | **连续（实数）阈值曲线**，折线不再是整数阶梯 | 1「不要阶梯」 |
| 13 | 42/43 | FAR=5% / 1% / 100 ppm 下 PoD50 与 PoD90 所需信号的均值 | 2 |
| 14 | 44/45 | 平方反比测远（纯 1/D² 与含大气衰减两种口径） | 3 |
| 15 | 46/47 | 同信号强度、不同 bg 时 peak 分布怎么变 | 4 |
| 16 | 48/49 | 宏像元 3×9 vs 3×6 阈值对比 | 旁支并入 |
| 17 | 50 | 理论汇总（阈值倍数 + 引擎一致性 + 口径说明） | 旁支并入 |

### 口径（v11 起沿用）

| 符号 | 含义 |
|---|---|
| `hist_i` | 第 i 发宏像元直方图 |
| `hist_add(N)` | 前 N 发之和；N∈{1,2,4} |
| **noise** | 单次 `hist_i` 统计窗均值 |
| **bg** | `hist_add` 统计窗均值；`noise_target`/`noise_mc` 字段表示目标/实测 **bg** |
| **peak** | 在 `hist_add` 上统计 |

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
模块 9.3/15：BOOST=[0,0.004,…,0.032]（9 档），N_MC_SIG_M9=8000
模块 11：N_MC_HSP=100_000
模块 13/14：FAR_M13=[5%, 1%, 100ppm]，LEVELS_M13=[0.50, 0.90]
CACHE_* = pod_esti_v20_cache_*.npz（noise/pod fallback 到 pod_esti_v11_cache_*.npz）
```

### 缓存现状（全部已跑完，2026-08-09）

| 缓存 | 规模 | 状态 |
|---|---|---|
| `pod_esti_v11_cache_noise.npz` | 48 bg × N∈{1,2,4} × 1,000,000 MC | 144/144 done |
| `pod_esti_v11_cache_pod.npz` | 144 档 PoD 临界点 | done（本轮补跑 32 档，13.8 min） |
| `pod_esti_v20_cache_signal.npz` | 48 bg × 9 boost × 3 N × 8,000 MC | 1296/1296 done（11.7 min） |
| `scan_hist_std_peak_cache.npz` | 48 bg × 3 N × 100,000 MC | done（模块 11 用） |
| `compare_macro_3x9_vs_3x6_cache.npz` | 24 档 p_eq × 7 配置 × 200,000 MC | done（模块 16 用） |

### 验证现状

`python check_v20_modules.py --syntax` → 51 cell，**0 语法错误**。
`python check_v20_modules.py --all` → 无头真跑 cell 28/31/33/35/37/39/41/43/45/47/49，
**0 运行期错误**，全部 `pod_v20_*.png` 已落盘。
（这只是无头验证，**还没在 Jupyter 里 Restart & Run All**。）

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

**模块 12（要求 1「不要阶梯」）**：新函数 `far_threshold_continuous(cnt, far)` 在
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

**没有卡住的地方。** 所有扫描已跑完，所有新模块已无头验证通过。剩下的都是可选收尾：

1. **还没在 Jupyter 里 Restart & Run All**。缓存都在，预计只花在绘图和渲染上，不会重算 MC。
   编辑器须先 **Revert/重开** `PoD_esti_v20.ipynb` 再 Restart → Run All。
2. 模块 14 的测距曲线逐点抖动几米（原因见上），要更平滑得加大 `N_MC_POD_VERIFY`
   并收紧 `POD_VERIFY_TOL`，代价是 PoD 扫描时间线性增长。
3. **（待查）** `theory_engine_equivalence.py` 的 T5b 里，步进引擎 `binary_macro_stepping`
   在 dt=800/200 ps 两档都比解析值**偏低**（2.2σ / 0.7σ，未达显著但方向一致），
   怀疑是它「先出 bin、再处理本步雪崩」的半步对齐。纯噪声扫描用快速引擎不受影响，
   但 **PoD 信号支路（`binary_macro_stepping_per_shot`）用的正是它**，建议加大样本复核。

---

## 4. 下一步计划

1. 打开 `PoD_esti_v20.ipynb`，**Revert/重开 → Restart → Run All**，确认渲染无误。
   若中途要重建内核：`python build_pod_core_v20.py`。
2. 用模块 10/12 的 MC 曲线逐条对照 `theory_peak_bg_multishot.md` 第 10 节的 7 条检查清单，
   把 ρ̄、残差、是否常数写回 `worklog_PoD_esti.md`。
   （模块 12 的连续阈值正是为这件事准备的 —— 别再用整数阈值比值。）
3. 复核 `binary_macro_stepping` 的半步对齐（见上第 3 条）。
4. （可选）加大 `N_MC_POD_VERIFY`、收紧 `POD_VERIFY_TOL`，抹平模块 14 测距曲线的抖动。
5. （可选）`compare_macro_3x9_vs_3x6.py` 若要 100 ppm 阈值，需 `--n-mc 2000000` 重跑
   （缓存键含 n_mc，会自动判失效重算）。
6. （可选）把「亚阈雪崩不复位」的非顺延变体做成开关，验证同 bg 下阈值曲线的二阶差异；
   解析上 p_bin 差异从 p_eq=0.0185 的 −0.4% 增大到 p_eq=0.4444 的 −13.7%，
   但因 `r_det` 是由目标 bg 反解的，这个差异主要被吸收进 klux 换算里。

### 全部重跑命令（正常不需要，缓存都在）

```powershell
$env:PYTHONIOENCODING="utf-8"
python build_pod_core_v20.py
python run_pod_v20_noise_scan.py --workers 20                 # 144 档 × 1e6 MC
python run_pod_v20_pod_scan.py   --workers 20                 # 约 14 min
python run_pod_v20_sig_scan.py   --workers 20                 # 约 12 min
python scan_hist_std_peak.py     --workers 20 --n-mc 100000   # 模块 11
python check_v20_modules.py --syntax
python check_v20_modules.py --all
```

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
