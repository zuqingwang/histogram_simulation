# -*- coding: utf-8 -*-
# 检查模块 15a：同 bg 下阈值比 T_i/T_j 是否等于单 bin 标准差比 σ_i/σ_j，
# 其中 σ = sqrt(bg*(1 - bg/n_tr)) = sqrt(bg - bg^2/(N_shots*N_SPAD))。
import sys
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

import numpy as np

Z = np.load("compare_macro_v30_cache.npz")
FAR_KEYS = [0.05, 0.01, 0.001, 100e-6]
I_FAR = FAR_KEYS.index(0.01)   # 模块 15 主 FAR = 1%

# cfgs 存的是 (nx, ny, N_shots)
CFGS = [tuple(int(v) for v in c) for c in np.asarray(Z["cfgs"])]
DONE = np.asarray(Z["done"])
BG = np.asarray(Z["bg_mc"])
T = np.asarray(Z[f"thr_{I_FAR}"])

# 15a 三条：(name, nx, ny, N, n_pix)
FOCUS = [
    ("3×6@N=2", 6, 3, 2, 18),
    ("3×9@N=2", 9, 3, 2, 27),
    ("3×9@N=4", 9, 3, 4, 27),
]

def load(nx, ny, n, n_pix):
    key = (nx, ny, n)
    ci = CFGS.index(key)
    ok = DONE[ci]
    bg = BG[ci][ok]
    thr = T[ci][ok]
    o = np.argsort(bg)
    bg, thr = bg[o].astype(float), thr[o].astype(float)
    n_tr = n_pix * n
    # 用户写的公式：σ = sqrt(bg - bg^2 / (#shot * #spad)) = sqrt(bg(1-bg/n_tr))
    sig = np.sqrt(np.maximum(bg * (1.0 - bg / n_tr), 0.0))
    return dict(bg=bg, T=thr, n_tr=n_tr, sig=sig)

D = {name: load(nx, ny, n, npx) for name, nx, ny, n, npx in FOCUS}

# 三条曲线的 bg 网格不完全相同（各配置扫描的 p_eq 相同，但 bg=n_tr·p_eq 不同）
# 取公共 bg 区间，对参考曲线插值到另一条的 bg 上再比。
PAIRS = [
    ("3×9@N=4", "3×9@N=2"),   # n_tr 108 / 54 = 2
    ("3×9@N=4", "3×6@N=2"),   # n_tr 108 / 36 = 3
    ("3×9@N=2", "3×6@N=2"),   # n_tr 54 / 36 = 1.5
]

print("=" * 100)
print("模块 15a：阈值比 vs 单 bin 标准差比")
print("σ = √[bg(1 − bg/n_tr)] = √(bg − bg²/(#shot·#spad))")
print("FAR = 1%；在公共 bg 区间上，把分子曲线插到分母的 bg 网格再取比值")
print("=" * 100)

for num, den in PAIRS:
    A, B = D[num], D[den]
    # 公共 bg：落在两条都覆盖的区间内，用分母的采样点
    lo = max(A["bg"].min(), B["bg"].min())
    hi = min(A["bg"].max(), B["bg"].max())
    m = (B["bg"] >= lo) & (B["bg"] <= hi) & (B["bg"] > 0.5)  # 跳过极低 bg（T 很小、整数量化狠）
    bg = B["bg"][m]
    T_b, s_b = B["T"][m], B["sig"][m]
    T_a = np.interp(bg, A["bg"], A["T"])
    s_a = np.interp(bg, A["bg"], A["sig"])
    rT = T_a / T_b
    rS = s_a / s_b
    # 也看 (T−bg) 比 —— 若阈值 ≈ bg + zσ，则超额比更该贴近 σ 比
    rE = (T_a - bg) / np.maximum(T_b - bg, 1e-9)

    print(f"\n--- {num} / {den} ---")
    print(f"  n_tr: {A['n_tr']} / {B['n_tr']} = {A['n_tr']/B['n_tr']:.3f}")
    print(f"  公共 bg 档数 {bg.size}，范围 {bg.min():.2f} → {bg.max():.2f}")
    print(f"  {'bg':>7} {'T_num':>7} {'T_den':>7} {'T比':>7} "
          f"{'σ_num':>7} {'σ_den':>7} {'σ比':>7} {'(T−bg)比':>9} {'T比/σ比':>9}")
    # 抽若干代表性档打印
    pick = np.unique(np.linspace(0, bg.size - 1, min(12, bg.size)).astype(int))
    for i in pick:
        print(f"  {bg[i]:7.2f} {T_a[i]:7.1f} {T_b[i]:7.1f} {rT[i]:7.3f} "
              f"{s_a[i]:7.3f} {s_b[i]:7.3f} {rS[i]:7.3f} {rE[i]:9.3f} {rT[i]/rS[i]:9.3f}")
    print(f"  汇总（全公共档）：")
    print(f"    T比     中位 {np.median(rT):.3f}  均值 {rT.mean():.3f}  "
          f"范围 {rT.min():.3f}–{rT.max():.3f}")
    print(f"    σ比     中位 {np.median(rS):.3f}  均值 {rS.mean():.3f}  "
          f"范围 {rS.min():.3f}–{rS.max():.3f}")
    print(f"    (T−bg)比 中位 {np.median(rE):.3f}  均值 {rE.mean():.3f}  "
          f"范围 {rE.min():.3f}–{rE.max():.3f}")
    print(f"    |T比 − σ比|     中位 {np.median(np.abs(rT-rS)):.3f}  均值 {np.mean(np.abs(rT-rS)):.3f}")
    print(f"    |(T−bg)比 − σ比| 中位 {np.median(np.abs(rE-rS)):.3f}  均值 {np.mean(np.abs(rE-rS)):.3f}")

# 低 bg 极限：bg ≪ n_tr ⇒ σ≈√bg，同 bg 下 σ比 → 1；T比也应接近 1
print("\n" + "=" * 100)
print("低 bg 极限（bg ≪ n_tr）：σ≈√bg，同 bg 下 σ比→1，若阈值≈bg+zσ 则 T比也→1")
print("=" * 100)
for num, den in PAIRS:
    A, B = D[num], D[den]
    lo = max(A["bg"].min(), B["bg"].min())
    hi = min(A["bg"].max(), B["bg"].max(), 3.0)  # 只看 bg≤3
    m = (B["bg"] >= lo) & (B["bg"] <= hi) & (B["bg"] > 0.5)
    if m.sum() < 2:
        print(f"  {num}/{den}: 低 bg 公共档不足")
        continue
    bg = B["bg"][m]
    rT = np.interp(bg, A["bg"], A["T"]) / B["T"][m]
    rS = np.interp(bg, A["bg"], A["sig"]) / B["sig"][m]
    print(f"  {num}/{den}（bg≤3，{m.sum()} 档）：T比中位={np.median(rT):.3f}，"
          f"σ比中位={np.median(rS):.3f}")
