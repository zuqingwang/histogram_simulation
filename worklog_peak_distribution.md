# worklog_peak_distribution

工作名：`peak_distribution`
主产物：`peak_distribution.ipynb` + `peak_distribution_scan.py`

> 文件名用正确拼写 `distribution`（用户口令里的 `distrubution` 是笔误）。

---

## 顶部：现状 / 当前任务

### 当前版本
v01（已跑通并写出结论）。

### 要回答的问题
**若信号强度 ×2，peak 分布的众数 / 中位数 / 均值 / p90 是否也 ×2？**
**PoD50 / PoD90 位置附近如何缩放？**

### 当前可运行状态
全部跑通。缓存 `peak_distribution_v01_cache.npz` 已齐；notebook 已 `nbconvert --execute`。

### 关键参数
| 参数 | 取值 |
|---|---|
| N_shots | 1, 4 |
| noise | N=1: 0.5/2/5/8；N=4: 1/5/15/30 |
| boost | 0 + logspace(0.001→0.08, 25)，共 26 档；×2 对实际倍率 ≈2.076 |
| N_MC | 8000 |
| FAR（PoD） | 100ppm, 10ppm |
| 并行 | ProcessPoolExecutor 20 |

### 结论摘要
- mode/median/mean/p90：**不会 ×2**（归一化中位 0.60）
- 净增量 Δ：弱信号近似线性，强信号次线性（归一化中位 0.93）
- PoD50 工作点再×2：仍不会（归一化中位 0.69）

### 待办
- [x] 扫描脚本
- [x] 全量扫描
- [x] notebook 出图与缩放比判定
- [x] handoff
- [ ]（可选）提高 boost 上限以覆盖 N=1 高噪声 PoD90

---

## 下部：历史记录（只追加）

### v01 — 本会话

**新增**
- `peak_distribution_scan.py`、`peak_distribution_v01_cache.npz`
- `build_peak_distribution.py` → `peak_distribution.ipynb`（14 cell，已执行）
- `handoff_peak_distribution.md`、本 worklog

**关键数值**
- 判定量 `norm = [stat(hi)/stat(lo)]/(E_hi/E_lo)`，线性 ⇔ 1
- 直接统计量 norm 中位 0.60（范围 0.48–0.96）
- 净增量 Δ norm 中位 0.93（弱信号≈1，强信号→0.6）
- PoD50 工作点 norm 中位 0.69
- 例：N=1 noise=0.5，FAR=100ppm：E50≈10.1 nJ，E90≈17.2 nJ（比≈1.70）

**踩过的坑**
- 对数网格「×2」实为 ×2.076 → 必须归一化
- 直接对含噪声本底的 mode/median 做比会严重偏离；要同时看 Δ
- PoD50/90 是能量位置，不是 peak 轴位置
- 弱信号 Δ 接近线性，不能把结论写成「全部非线性」
- **图中中文全是方框**：`peak_distribution.ipynb` 只读缓存、不 import `pod_esti_v05_core`，
  因此没带上 core 里的 YaHei 字体设置。修法：在首个 code cell 显式设
  `matplotlib.rcParams["font.sans-serif"] = ["Microsoft YaHei", ...]`
  并 `axes.unicode_minus = False`；顺带把 `×→『』Δ≥` 等易缺字形符号改成 ASCII。
