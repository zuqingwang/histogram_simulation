# handoff_peak_distribution

> 面向**完全没有上下文的新会话**。过程流水账见 `worklog_peak_distribution.md`。

---

## 1. 我们在做什么任务

工作名 **`peak_distribution`**（文件名用正确拼写；用户口令里的 `distrubution` 是笔误）。

要回答：

> **若信号强度 ×2，peak 分布的众数 / 中位数 / 均值 / p90 是否也 ×2？**
> **PoD50 / PoD90 的位置附近，上述量又如何缩放？**

与姊妹工作 `peak_vs_noise`（固定信号、扫噪声）正交：本工作**固定噪声、扫信号**，专门检验比例缩放。

物理内核：`pod_esti_v05_core`（不复制参数）。  
PoD（Probability of Detection，检测概率）= `P(peak ≥ T)`，`T` 取自 `PoD_esti_v05` 的 FAR 阈值。

### 目标产物
| 文件 | 作用 |
|---|---|
| `peak_distribution_scan.py` | 多进程扫描 |
| `peak_distribution_v01_cache.npz` | 缓存（bincount） |
| `build_peak_distribution.py` | notebook 生成器 |
| `peak_distribution.ipynb` | 分析出图（已 execute） |

---

## 2. 已经完成了什么

### 扫描设计
- `N_shots ∈ {1, 4}`
- noise 代表档：N=1 → `0.5/2/5/8`；N=4 → `1/5/15/30`
- boost：`0` + `logspace(0.001→0.08, 25)` 共 26 档；相邻「加倍」对实际倍率 ≈ **2.076**（对数步长）
- 每点 **8000** 条 MC；`ProcessPoolExecutor` 20 进程；全量 8 档约 **5–6 min**
- 缓存键：`boosts` + `noise_*` + `n_mc`

### notebook 模块（14 cell）
1. 载入 + 统计量（mode/median/mean/p90）+ `(b,≈2b)` 对索引  
2. 分布对照图（含「横轴×2」线性假设虚线）  
3. 归一化缩放比 `norm = [stat(hi)/stat(lo)]/(E_hi/E_lo)`（线性 ⇔ =1）  
4. 净增量 Δ 及其 norm  
5. PoD–能量曲线 + PoD50/PoD90 临界能量  
6. 在 PoD50(FAR=100ppm) 工作点把信号×2 的缩放  
7. 自动结论

### 关键实测结论（v01）

判定量 = `[stat(hi)/stat(lo)] / (E_hi/E_lo)`，线性 ⇔ **1**。

| 检验对象 | 归一化缩放比（汇总） | 是否「也×2」 |
|---|---|---|
| mode / median / mean / p90 | 中位 **0.60**，范围 0.48–0.96 | **否** |
| 净增量 Δ（扣噪声本底） | 中位 **0.93**，弱信号≈1，强信号→0.6 | **弱信号近似是，强信号否** |
| PoD50 工作点再×2 | 中位 **0.69** | **否** |

PoD50/PoD90 示例（FAR=100ppm）：
- N=1, noise=0.5：E50≈10.1 nJ，E90≈17.2 nJ，E90/E50≈1.70
- N=4, noise=5：E50≈5.0 nJ，E90≈7.5 nJ，E90/E50≈1.50
- N=1, noise≥5：部分 E90 越出 boost 网格上限（nan）

**一句话**：peak 分布的众数/中位数等**不会**随信号×2 而×2；扣掉噪声后的净增量只在小信号区近似线性；PoD50/90 是能量位置，在该工作点加倍信号，peak 统计仍不会×2。

---

## 3. 目前面临的问题、卡在哪里

无阻塞。局限：
1. 只扫了代表性 noise，不是全网格。
2. N=1 高噪声下 PoD90 可能超出 boost 上限（0.08 ≈ 64 nJ）。
3. `(b,2b)` 对实际倍率是 2.076 不是精确 2 —— 已用归一化处理。

---

## 4. 下一步计划

1. （可选）把 boost 上限提到 0.2，补齐 N=1 高噪声的 PoD90。  
2. （可选）与 `peak_vs_noise` 联立，做 `(noise, boost)` 二维插值表供固件。  
3. 若改了 `PoD_esti_v05` 计算 cell → 先 `build_pod_core_v05.py` 再重跑本扫描。

运行：
```powershell
$env:PYTHONIOENCODING="utf-8"
python peak_distribution_scan.py
python build_peak_distribution.py
python -m nbconvert --to notebook --execute --inplace peak_distribution.ipynb
```

---

## 5. 踩过的坑（不要再踩）

1. **对数网格的「×2」不是精确 2** → 必须用 `norm = raw_ratio / (E_hi/E_lo)`，期望值才是 1。  
2. **直接对 mode/median 做 ratio 会严重低估线性** → 噪声本底拉高了分母；要同时看净增量 Δ。  
3. **PoD50/90「位置」不是 peak 轴上的点** → 是能量轴上的临界值；「信号×2 → 位置×2」这句话本身不成立，应改问工作点处统计量缩放。  
4. notebook `source` 必须每行自带 `\n`（见 `build_peak_distribution.py` 的 `_lines()`）。  
5. 用 `python -m nbconvert`，不要指望 PATH 里的 `jupyter`。  
6. 结论文字必须由实测数自动生成；弱信号 Δ≈线性、强信号次线性，不能一刀切写「全部非线性」。
