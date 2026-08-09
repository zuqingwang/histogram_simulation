# 交接文档 —— 当前工作：`PoD_esti`（探测概率估计）

> 文件名：`handoff_PoD_esti.md`（禁止 `handoff_现在工作.md`）。
> 最后更新：2026-08-09（**v11**：统一 bg 步长 0.25 + 模块 10 阈值倍数分析；
> 新增宏像元 3×9 vs 3×6 对比；完成与 `lidar_histogram_sim_v45.ipynb` 的引擎一致性核对）。
> 流水日志：`worklog_PoD_esti.md`。
>
> **未完成的主任务**：`PoD_esti_v11.ipynb` 尚未整本全量执行。用户会主动提醒继续跑，
> 具体命令见本文第 4 节第 1 条，以及 `worklog_PoD_esti.md` 顶部「6. 待续跑任务」。

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

主产物：**`PoD_esti_v11.ipynb`**  
- 基于完整 **`PoD_esti_v10.ipynb`**（v05 全模块 + hist_i + 模块 9）。  
- **v11 核心**：各 N 的目标 **bg 网格统一**为步长 0.25；新增 **模块 10** 分析同 bg 下阈值倍数。

### 口径（v11）

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
| **v11** | **`PoD_esti_v11.ipynb`** | **当前主产物** |

### 关键文件

| 文件 | 作用 |
|---|---|
| `PoD_esti_v11.ipynb` | 主 notebook（38 cell，含模块 10） |
| `upgrade_pod_esti_v11_from_v10.py` | v10→v11 升级脚本 |
| `build_pod_core_v11.py` / `pod_esti_v11_core.py` | 多进程内核 |
| `run_pod_v11_noise_scan.py` | 噪声 ProcessPool；任务键 `(N,bg)` |
| `run_pod_v11_pod_scan.py` | PoD ProcessPool |
| `pod_esti_v11_cache_*.npz` | 新缓存；`FALLBACK=[]`，**禁止**读 v10 |
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
CACHE_* = pod_esti_v11_cache_*.npz
```

### 冒烟验证（已通过）

`python run_pod_v11_noise_scan.py --workers 4 --limit 2 --n-mc 2000`  
同 bg=0.25 时 N=1/2/4 的 peakμ≈1.65–1.67（接近，符合低噪近似 ρ→1）。  
**冒烟缓存已删除**，正式跑须全量重算。

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

## 3. 目前面临的问题、卡在哪里

1. **尚未整本 Run All**：噪声 48×3=144 档 ×1e6 + PoD 同规模，耗时长。
   用户表示扫描由他自己择机跑，会主动来提醒继续。命令见第 4 节第 1 条。
2. 模块 10 的定量结论（ρ̄、残差、是否常数）要等噪声缓存 + THRESH 算完才有。  
3. 编辑器须 **Revert/重开** `PoD_esti_v11.ipynb` 后再 Restart → Run All。
4. 旁支分析（宏像元对比、引擎核对）**都已跑完**，不阻塞主线。

---

## 4. 下一步计划

1. **【待续跑，用户会提醒】整本执行 `PoD_esti_v11.ipynb`**。推荐先在命令行把两个耗时扫描跑完，
   再回 notebook 里 Restart & Run All（那时缓存已在，模块 9/10 直接出图）：

   ```powershell
   $env:PYTHONIOENCODING="utf-8"
   python build_pod_core_v11.py            # 只有改过 notebook 才需要重新导出内核
   python run_pod_v11_noise_scan.py --workers 20    # 48 bg × N∈{1,2,4} = 144 档 × 1e6 MC
   python run_pod_v11_pod_scan.py  --workers 20     # 依赖上一步的阈值
   ```

   然后打开 `PoD_esti_v11.ipynb`，**Revert/重开 → Restart → Run All**。
2. 用模块 10 的 MC 曲线逐条对照 `theory_peak_bg_multishot.md` 第 10 节的 7 条检查清单，
   把 ρ̄、残差、是否常数写回 `worklog_PoD_esti.md`。
3. （可选）`compare_macro_3x9_vs_3x6.py` 若要 100 ppm 阈值，需 `--n-mc 2000000` 重跑
   （缓存键含 n_mc，会自动判失效重算）。
4. （可选）把「亚阈雪崩不复位」的非顺延变体做成开关，验证同 bg 下阈值曲线的二阶差异；
   解析上 p_bin 差异从 p_eq=0.0185 的 −0.4% 增大到 p_eq=0.4444 的 −13.7%，
   但因 `r_det` 是由目标 bg 反解的，这个差异主要被吸收进 klux 换算里。
5. （待查）`theory_engine_equivalence.py` 的 T5b 里，步进引擎 `binary_macro_stepping`
   在 dt=800/200 ps 两档都比解析值**偏低**（2.2σ / 0.7σ，未达显著但方向一致），
   怀疑是它「先出 bin、再处理本步雪崩」的半步对齐。纯噪声扫描用快速引擎不受影响，
   但 **PoD 信号支路（`binary_macro_stepping_per_shot`）用的正是它**，建议加大样本复核。

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
