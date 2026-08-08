# handoff_peak_vs_noise

> 面向**完全没有上下文的新会话**。本文只写事实与结论，过程流水账见 `worklog_peak_vs_noise.md`。

---

## 1. 我们在做什么任务

### 项目背景
本仓库（`Histogram-simulation`）模拟一台 **TCSPC**（Time-Correlated Single Photon Counting，
时间相关单光子计数）直方图式激光雷达的接收链路：
**SPAD**（Single-Photon Avalanche Diode，单光子雪崩二极管）阵列 + 二值采样 + 1 ns 直方图。

仓库里有若干条平行工作线，各自有独立的 handoff / worklog：

| 工作名 | 主产物 | 职责 |
|---|---|---|
| `lidar_histogram_sim` | `lidar_histogram_sim_v45.ipynb` | 物理内核源头（光链路、SPAD、二值采样） |
| `PoD_esti` | `PoD_esti_v05.ipynb` | 检测概率（Probability of Detection）与 FAR 阈值 |
| `crosstalk_sim` | `crosstalk_sim_v23.ipynb` 等 | 串扰 / 雷达对射与 tcode 编码 |
| **`peak_vs_noise`（本文）** | **`peak_vs_noise_v02.ipynb`**（v01 保留） | **peak 分布随环境噪声的演化规律** |

### 当前工作名与目标
工作名 **`peak_vs_noise`**，要回答一个具体问题：

> **给定固定的信号强度，环境噪声逐渐增强时，
> 信号峰值 peak 的分布如何变化？是线性变化吗？**

### 目标产物
1. `peak_vs_noise_scan.py` —— 多进程扫描脚本，产出缓存。全量 210 档 / 20 进程 / **实测 21.6 min**。
2. `peak_vs_noise_v01_cache.npz` —— 扫描结果缓存（**65 KB**，存的是 bincount）；**v02 复用，不重扫**。
3. `build_peak_vs_noise_v01.py` → `peak_vs_noise_v01.ipynb`（旧版，仅 noise 轴）。
4. **`build_peak_vs_noise_v02.py` → `peak_vs_noise_v02.ipynb`（当前主产物）**：
   凡横轴为 noise 的图都再画一张 bg 轴图；已 `nbconvert --execute` 通过（约 47 s）。
   **由 builder 生成，不要直接手编 notebook**，否则重生成会覆盖。

### 已经得到的答案（结论摘要）

**不是线性的。但必须分两个层次说，否则会得出错误印象：**

| 层次 | 现象 | 数值 |
|---|---|---|
| ① 单条 `<peak>`–noise 曲线的形状 | **看起来几乎是直线**，严格说是次线性（饱和型） | 全程一次拟合 R² = **0.9966…1.0000**；偏离低噪声切线最多 17.1% n_tr（纯噪声档）/ 2.2% n_tr（信号档 N=4） |
| ② 信号与噪声能否叠加 | **明确不能**，这是强非线性所在 | `∂<peak>/∂noise` 随信号强度**下降 55%**（N=1 与 N=4 都是 55%） |
| ③ 信号的净贡献 | 被持续压缩 | Δ 衰减到无噪声时的 **25%–50%** |
| ④ 检测层面的后果 | 衰减最剧烈 | 可分辨度 d′ 最差只剩低噪声时的 **16%** |

**层次 ② 是最干净的判据**：若 `<peak> = f(信号) + g(noise)` 可线性叠加，
则 `∂<peak>/∂noise` 必与信号强度**无关**。实测该斜率随信号从 0 加到 25.58 nJ 单调降 55%
（N=1：1.624 → 0.736；N=4：1.385 → 0.619），直接否定可分离性。

**物理机制（两条，叠加）**
- **二值采样硬上限**：每个 SPAD 每 bin 至多 1 个计数，`peak ≤ n_tr = 27 × N_shots`，
  越接近上限增速越慢 —— 解释层次 ①的弯曲。
- **死时间抢占**：噪声光子先触发某个 SPAD 后，该 SPAD 在 RC 恢复期
  （`τ_RC = 8.73 ns`，过阈窗 `T_OVER ≈ 8 ns`）内无法再响应信号光子。
  信号与噪声**竞争同一批 SPAD**，不是各自累加 —— 解释层次 ②③。

**工程含义**：`<peak> ≈ a + b·noise` 作为估算是可用的近似；
但**不能**由此推论"信号和噪声可以分别算再相加"。
这也是 `PoD_esti` 中"临界能量随 noise 快速上升"的微观原因。

### 关键名词口径（全项目一致，务必沿用；v02 起拆分 noise / bg）
- **noise（环境标准）**：折合到 **N_shots=1、宏像元 27 SPAD、每 1 ns bin** 的平衡态底计数。
  **与发数无关**。由 `E_lambda → r_det → 27·p_bin_equilibrium` 得到。
- **bg（实测 baseline）**：当前 `N_shots` 下，统计窗内累加直方图的实测平均底
  （缓存字段 `noise_mc` / `noisemc_{n}`）。**N_shots=1 时 bg≈noise；N_shots=4 时 bg≈4·noise**。
- **旧代码里叫 `noise` 的量**（如目标网格、`noise_mc`）多数其实是 **bg（累加后的底）**。
- **peak**：单次测量中，信号窗内宏像元直方图的**最大 bin 计数**。
- **二值采样硬上限** `n_tr = 27 × N_shots`：每个 SPAD 每 bin 最多贡献 1 个计数，
  所以 **N_shots=1 时 peak ≤ 27，N_shots=4 时 peak ≤ 108**。
- **boost**：信号强度倍率。单脉冲能量 = `boost × E_PULSE_BASE`，
  **`E_PULSE_BASE` = 799.4 nJ**（由 `P_peak=235 W` 与双指数 τ_r=0.7 / τ_f=1.9 ns 积分得出）。
- **FAR**（False Alarm Rate，虚警率）：检测阈值 T 的选取依据。
  阈值曲线取自 `PoD_esti_v05` 的百万条纯噪声 MC，共六档：
  `10ppm / 100ppm / 0p1pct / 0p5pct / 1pct / 5pct`。

---

## 2. 已经完成了什么

### 2.1 计算内核：直接复用，不复制
`peak_vs_noise_scan.py` 顶部 `import pod_esti_v05_core as core`，
**没有复制任何物理参数**。所以本工作的物理与 `PoD_esti_v05` 逐项一致。
用到的 core 接口只有这几个：

| core 接口 | 用途 |
|---|---|
| `core._peaks_chunk(boost, n_shots, r_amb, n_real, seed)` | 跑一块 MC，返回 peak 数组（**主力**） |
| `core.NOISE_RES[n_shots]["r_det"][k]` | 第 k 档 noise 对应的探测率；`r_amb = r_det / core.PDE` |
| `core.NOISE_GRID` | noise 网格（N=1：48 档；N=4：160 档） |
| `core.THRESH[n_shots]["T"+tag]` | 六档 FAR 阈值曲线 |
| `core.E_PULSE_BASE`、`core.PDE`、`core.FAR_TAGS` | 换算与标签 |
| `core._atomic_savez` | 原子写缓存 |

**注意**：`import pod_esti_v05_core` 会载入 `pod_esti_v05_cache_noise.npz`
并计算六档 FAR 阈值，约耗时 5–10 s。若该噪声缓存缺失会触发百万条噪声 MC，
所以**必须先保证 `pod_esti_v05_cache_noise.npz` 存在**。

### 2.2 扫描设计（落在 `peak_vs_noise_scan.py`）

| 项 | 取值 | 位置 |
|---|---|---|
| `BOOST_LIST` | `[0.0, 0.004, 0.008, 0.016, 0.032]` | 模块顶部常量 |
| 对应脉冲能量 | `0 / 3.20 / 6.40 / 12.79 / 25.58` nJ | 第 0 档（boost=0）**就是纯噪声基线** |
| noise 网格 | 直接取 `core.NOISE_GRID`（N=1 共 48 档，N=4 共 160 档） | `_job()` |
| 额外参考档 | `k = K_NOFLOOR = -1`，表示 **noise = 0（`r_amb = 0`）** | `_job()` |
| `N_MC_DEFAULT` | 8000 条 / (n_shots, noise, boost) | 可用 `--n-mc` 覆盖 |
| `MC_CHUNK` | 2000（控制单进程峰值内存） | `_peak_bincount()` |
| 并行 | `ProcessPoolExecutor`，20 进程；任务粒度 = 1 个 noise 档（5 个 boost 一起算） | `main()` |
| 缓存 | `peak_vs_noise_v01_cache.npz` / 检查点 `...partial.npz`（每 10 档） | `_save()` / `_load()` |

**存储格式的关键决定**：缓存里存的不是原始 peak 样本，而是 **peak 的 `bincount`**。
因为 peak 取值域只有 `0…n_tr`（≤108），`bincount` 就是**完整分布的充分统计量**，
既能事后算任意分位数与分布图，整个扫描又只占约 1 MB。

缓存内的字段（`n` 取 1 和 4）：

| 字段 | 形状 | 含义 |
|---|---|---|
| `grid_key` / `boosts` / `n_mc` | — | **缓存键**，任一项变化则旧缓存作废 |
| `noise_{n}` | (n_noise,) | 目标 noise（横轴） |
| `noisemc_{n}` | (n_noise,) | MC 实测 noise |
| `cnt_{n}` | (n_boost, n_noise, n_tr+2) | peak 的 bincount |
| `cnt0_{n}` | (n_boost, n_tr+2) | **noise = 0** 的无噪声纯信号参考档 |
| `done_{n}` / `done0_{n}` | (n_noise,) / 标量 | 完成标记，支持断点续跑 |
| `T_{n}` | (6, n_noise) | 六档 FAR 阈值，顺序同 `far_tags` |

### 2.3 分析 notebook（当前主产物 `peak_vs_noise_v02.ipynb`）
由 `build_peak_vs_noise_v02.py` 生成。相对 v01：**每个原 noise 轴图都再画一张 bg 轴图**。
载入时派生 `noise`（环境标准）与 `bg`（=`noise_mc`）；`AXIS_KINDS` 双循环；`TANGENT_BY_AXIS` 按轴分存。

v01 模块结构仍适用（模块编号相同）；v02 每个作图模块输出 ×2（noise 轴 + bg 轴）。

### 2.4 运行方式
```powershell
$env:PYTHONIOENCODING="utf-8"
python peak_vs_noise_scan.py          # 全量 210 档；已有 v01 缓存则秒退
python build_peak_vs_noise_v02.py     # 重新生成 v02 notebook（改图/改分析只改 builder）
python -m nbconvert --to notebook --execute --inplace peak_vs_noise_v02.ipynb
```
或者打开 `peak_vs_noise_v02.ipynb`，Restart Kernel → Run All（只读缓存）。

注意：`jupyter nbconvert` 在本机 PowerShell 里**不在 PATH**，必须用 `python -m nbconvert`。

---

## 3. 目前面临的问题、卡在哪里

**没有卡住的问题。** v01 扫描缓存可用；v02 双轴出图已 execute 通过。
实测：N=1 中位 bg/noise=1.00；N=4 中位 bg/noise=4.00。

需要注意的现存局限（不是 bug，是范围限定）：
1. 只扫了**一个目标距离 D = 15.0 m、一种反射率 ρ = 0.10**（沿用 `PoD_esti_v05` 场景）。
   距离/反射率变化时结论的**定性**形状应不变（机制是硬上限 + 死时间抢占），但**定量**数值会变。
2. `DCR`（Dark Count Rate，暗计数率）设为 0，噪声全部来自环境光。
3. 每档 8000 条 MC，均值的统计误差约 σ/√8000（σ ≈ 1.5–4.5 计数 ⇒ 约 0.02–0.05 计数），
   足够看均值趋势；但分布的极远尾部（如 1e-4 分位）不够，需要时要加大 `--n-mc`。
   这也是模块 2 局部斜率曲线看起来有抖动的原因（尤其 N=4），趋势本身是可靠的。
4. 结论的**定量**数值绑定在本扫描范围（N=1 的 noise ≤ 12，N=4 的 noise ≤ 40）。
   `<peak>` 曲线"接近直线"这一点**依赖于扫描范围**：范围若一直外推到硬上限附近，
   饱和弯曲会显著得多。引用 R² 时务必带上 noise 范围。

---

## 4. 下一步计划

按优先级：

1. **（可选）扩展到多个距离**：把 `D_TARGET` 做成扫描维度，
   验证"次线性 + 次可加"的结论在近距离（强信号、易饱和）与远距离（弱信号）两端都成立。
   改动位置：`peak_vs_noise_scan.py` 的 `_job()`，需要 core 支持按距离重建 `ECHO0` / `R_SIG_UNIT_POD`。
2. **（可选）解析近似**：给出 `<peak>(noise, boost)` 的经验闭式，供固件快速估算。
3. **（可选）加大尾部精度**：`python peak_vs_noise_scan.py --n-mc 200000`。
4. **若改了 `PoD_esti_v05` / v06 的计算 cell**：先同步 `pod_esti_v05_core.py`，再决定是否重扫本工作。
5. 兄弟线：`PoD_esti_v06.ipynb` 已按同一口径加 bg 轴图；出图前需在 notebook 中命中 v05 缓存后重跑作图 cell。

---

## 5. 踩过的坑（不要再踩）

### 坑 1：沿用别的模块的 boost 范围 → 信号全在饱和区
- **错误做法**：按 `PoD_esti` 的 `POD_LOG_BOOST_MIN/MAX` 直觉取 `boost ≈ O(1)`。
- **现象**：N=1 时 peak 均值 13→25，几乎顶满硬上限 27；噪声的影响被硬上限压掉，看不出任何规律。
- **根因**：`E_PULSE_BASE = 799.4 nJ` 是很强的脉冲，15 m / ρ=0.10 目标下
  `boost = 0.05` 就足以点亮 27 个 SPAD 的一半。
- **正确做法**：任何新的信号强度扫描都要**先探测再定档**。
  本工作最终取 `boost = 0.004…0.032`，使低噪声下 peak 占硬上限的 10%–40%。

### 坑 2：`E_PULSE_BASE` 不是参数表里的 800 nJ
- **错误做法**：看到 `PARAMS["laser"]["E_pulse"] = 800e-9` 就以为单脉冲能量是 800 nJ 并手推 boost。
- **现象**：第一版 `BOOST_LIST` 把探测表的「能量列 ÷ 1000」误当成 boost，
  写成 `0.0032…0.0256`（实际 E = 2.56…20.47 nJ，与预期的 3.2…25.6 nJ 不符）。
- **根因**：`amp_mode = "peak"`，实际能量由 `P_peak = 235 W` 与双指数波形积分得出 = **799.4 nJ**；
  与参数表里的 800e-9 数值接近纯属巧合，该参数在 `amp_mode="peak"` 下**不生效**。
- **正确做法**：boost ↔ 能量一律写成 `boost * core.E_PULSE_BASE * 1e9` 现算，绝不手推。

### 坑 3：用 `ThreadPoolExecutor` 跑 MC → CPU 吃不满
- **错误做法**：沿用早期 `PoD_esti` 的多线程方案。
- **现象**：CPU 总占用只有 16%–30%。
- **根因**：`core.binary_macro_stepping` 是「Python 层逐细网格步循环 + 小数组 NumPy」，
  几乎全程持有 **GIL**（Global Interpreter Lock，全局解释器锁）。
- **正确做法**：用 `ProcessPoolExecutor`（本工作与 `run_pod_scan_v05.py` 都是这样）。
  同时在 core 顶部把 `OMP_NUM_THREADS` 等 BLAS 变量设为 1（**必须在 `import numpy` 之前**），
  否则 20 进程 × N 线程互相抢核反而更慢。

### 坑 4：`np.savez_compressed` + `os.replace` 的 Windows 陷阱
- **错误做法**：临时文件名不带 `.npz` 后缀。
- **现象**：`FileNotFoundError: [WinError 2]`。
- **根因**：`np.savez_compressed` 会给不以 `.npz` 结尾的路径**自动追加 `.npz`**，
  于是 `os.replace(tmp, path)` 找不到 `tmp`。
- **正确做法**：临时文件名显式写成 `xxx.tmp.npz`。本工作直接复用已修好的 `core._atomic_savez`。

### 坑 5：改了 notebook 后编辑器仍跑旧代码
- **现象**：磁盘上的 notebook 已修好，但执行时仍报早已修掉的错。
- **根因**：Jupyter / Cursor 的 notebook 编辑器缓存了内存中的旧版本。
- **正确做法**：**不保存**地关闭 notebook 标签页 → 重新打开（或 `File: Revert File`）→ Restart Kernel。
  尤其在跑完 `build_peak_vs_noise_v01.py` 重新生成 notebook 之后，这一步是必需的。

### 坑 6（生成器相关）：在 `r"""..."""` 里嵌套三引号 / 用 `\"` 转义
- **错误做法**：在 `build_peak_vs_noise_v01.py` 的 `code(r"""...""")` 内部写
  `\"\"\"docstring\"\"\"` 或 `print("""...""")`。
- **现象**：raw string 里 `\"` 会**保留反斜杠**，生成出语法错误的 cell；
  嵌套的 `"""` 会**提前终止**外层字符串。
- **正确做法**：生成的代码内部一律用**三单引号** `'''...'''`。
  改完 builder 后跑一遍逐 cell `compile()` 校验（本会话用过临时脚本 `_check_nb.py`，已删）。

### 坑 7：matplotlib 图上用 mathtext 中文下标 / 数学角括号
- **错误做法**：`ax.set_ylabel("⟨peak⟩$_{纯信号}$")`。
- **现象**：两个毛病。① mathtext 对中文支持很差，会渲染失败或大量警告；
  ② `⟨` `⟩`（U+27E8/27E9）**在 Microsoft YaHei 里缺字形**，
  报 `UserWarning: Glyph 10216 missing from font(s)`，图上显示为方框。
- **正确做法**：数学模式 `$...$` 内只放拉丁字母/符号（如 `$n_{tr}$`），中文写在外面；
  角括号用 ASCII 写成 `<peak>`。

### 坑 8：builder 里 notebook 的 `source` 用 `split("\n")`
- **错误做法**：`"source": src.split("\n")`。
- **现象**：Jupyter 把整个 cell 的**所有行拼成一行** → `SyntaxError: invalid syntax`。
  更坑的是，如果自己写的逐 cell 校验脚本自行补了换行符，会**误判为通过**。
- **根因**：nbformat 的 `source` 是「**每行自带换行符**」的字符串列表（末行不带）。
- **正确做法**：见 `build_peak_vs_noise_v01.py` 的 `_lines()`：
  `[ln + "\n" for ln in body.split("\n")[:-1]] + [body.split("\n")[-1]]`。
  顺带给每个 cell 加 `id` 字段，消除 `MissingIDFieldWarning`。

### 坑 9：在 `r"""..."""` 里写 `\\n` 想给图标题换行
- **现象**：图标题里出现**字面的 "\n" 两个字符**。
- **根因**：r-string 中 `\\n` 是 3 个字符，经 JSON 落到 notebook 源码仍是 `\\n`，
  Python 解析后值是 `\` + `n`，不是换行符。
- **正确做法**：在 r-string 里就写单个 `\n`。

### 坑 10：`jupyter nbconvert` 不在 PATH
- **现象**：`jupyter : 无法将"jupyter"项识别为 cmdlet`；
  改用 `python -m jupyter nbconvert` 则报 `Jupyter command jupyter-nbconvert not found`。
- **正确做法**：用 **`python -m nbconvert`**。

### 坑 11（最重要，方法论）：先写死结论文字，再拿数据去"验证"
- **错误做法**：在 notebook 结论 cell 里预先写死"⟨peak⟩ 随 noise 明显次线性、
  高噪声端斜率只剩低噪声端的一小部分"。
- **现象**：实测出来 **全程线性拟合 R² = 0.9966…1.0000**，信号档偏离线性外推最多只有
  2.2% n_tr（N=4），预设措辞与数据明显不符，等于把结论的强度写反了。
- **根因**：把"机制上必然非线性"直接等同于"曲线形状上明显弯曲"。
  实际上二值硬上限造成的曲率在本扫描范围内很弱；真正的强非线性表现在
  **可分离性**（斜率随信号强度变）与**信号净增量衰减**上，不在单条曲线的曲率上。
- **正确做法**：结论里的**判定措辞与全部数值都由实测自动生成**（现在 cell 15 就是这样）；
  非线性要分层论证（曲线曲率 / 可分离性 / 净增量 / d′），不能只看一条曲线直不直。
