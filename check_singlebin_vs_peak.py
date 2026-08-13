"""分离验证三个不同的统计量，别再把它们混为一谈。

  A. 单 shot、单 bin 的计数：27 个 SPAD 各自「亮/不亮」之和。
     若 27 个 SPAD 的点亮概率都相同，它是二项分布 Binomial(27, p)，方差 27p(1-p)，p=0.5 见顶。
     但 f_pix（像斑在 27 个 SPAD 上的收集比例）差异极大 ⇒ 各 SPAD 的 p_t 不同
     ⇒ 实际是【泊松二项分布】，方差 Σ p_t(1-p_t) = 27·p̄(1-p̄) − 27·Var_t(p_t)，**恒低于**二项。
  B. N=4 累加、单 bin 的计数：4 个独立 shot 相加，方差应为 A 的 4 倍。
  C. peak：152 个 bin 上的【最大值】。与 A、B 是不同的统计量。

每个能量档的「最亮 bin」按【该档自己的】平均波形定，不能用最强档的（那时已是 30 ns 平台，峰位漂移）。
"""
import os
import sys

os.environ.setdefault("POD_CORE_QUIET", "1")
sys.stdout.reconfigure(encoding="utf-8")

import numpy as np
import pod_esti_v30_core as core
import run_peak_energy_scan as scan

N_MC = 6000
N_SHOTS = 4
N_PIX = core.N_PIX_MACRO
N_TR = N_PIX * N_SHOTS

# 本脚本记录的结论是在真实像斑轮廓下测的，显式锁定 F_PIX_REAL，
# 不要跟着 scan.F_PIX_MODE（现默认 uniform）漂移，否则和已存的输出对不上。
f = scan.F_PIX_REAL
print(f"宏像元 {N_PIX} 个 SPAD 的收集比例 f_pix：")
print(f"  最大 {f.max():.4g}，最小 {f.min():.4g}，均值 {f.mean():.4g}，"
      f"相对离散度 std/mean = {f.std() / f.mean():.2f}")
print(f"  最大的 5 个：{np.sort(f)[::-1][:5].round(4)}")
print(f"n_tr = {N_TR}（{N_PIX} SPAD × {N_SHOTS} shots），每档 {N_MC:,} 次 MC，bg = 0\n")

grid = [0.005, 0.01, 0.0175, 0.03, 0.05, 0.09, 0.15, 0.25, 0.5, 1.0, 2.6, 5.0, 12.93, 33.4]

print("A/B = 单 bin（该档自己的最亮 bin）　C = peak（152 bin 取最大）")
print(f"{'boost':>8} | {'A: 单shot单bin':>30} | {'B: N=4单bin':>16} | {'C: peak':>16}")
print(f"{'':>8} | {'p̄':>6}{'实测σ':>7}{'二项σ':>7}{'实测/二项':>9} | {'占比':>7}{'σ':>8} | {'占比':>7}{'σ':>8}")
print("-" * 84)

rows = []
for b in grid:
    rng = np.random.default_rng(4242)
    h1, h4, pk = [], [], []
    done = 0
    while done < N_MC:
        m = min(2000, N_MC - done)
        hi = core.binary_macro_stepping_per_shot(
            m, f, N_SHOTS, core.R_SIG_UNIT_GEN, core.TF_GEN,
            0.0, core.CENTERS, rng, boost=float(b))
        add = hi.sum(axis=1)[:, core.I_STAT0:core.I_STAT1]
        h1.append(hi[:, 0, core.I_STAT0:core.I_STAT1])
        h4.append(add)
        pk.append(add.max(axis=1))
        done += m
    h1 = np.concatenate(h1).astype(float)
    h4 = np.concatenate(h4).astype(float)
    pk = np.concatenate(pk).astype(float)

    j = int(h4.mean(axis=0).argmax())          # 该档自己的最亮 bin
    a1, a4 = h1[:, j], h4[:, j]

    pbar = a1.mean() / N_PIX
    sd_binom1 = np.sqrt(N_PIX * pbar * (1 - pbar))
    rows.append((b, pbar, a1.std(), sd_binom1, a4.mean() / N_TR, a4.std(),
                 pk.mean() / N_TR, pk.std(), a1.var()))
    print(f"{b:>8.4g} | {pbar:>6.3f}{a1.std():>7.3f}{sd_binom1:>7.3f}"
          f"{a1.std() / max(sd_binom1, 1e-9):>9.3f} | "
          f"{a4.mean() / N_TR:>7.3f}{a4.std():>8.3f} | {pk.mean() / N_TR:>7.3f}{pk.std():>8.3f}")

r = np.array(rows)
print("\n" + "=" * 84)
i1 = int(r[:, 2].argmax())
i4 = int(r[:, 5].argmax())
ic = int(r[:, 7].argmax())
print(f"A 单shot单bin  σ 最大 {r[i1, 2]:.3f} @ p̄ = {r[i1, 1]:.3f}"
      f"（纯二项理论应 @ p=0.5，极大 {np.sqrt(N_PIX * 0.25):.3f}）")
print(f"B N=4 单bin   σ 最大 {r[i4, 5]:.3f} @ 占比 {r[i4, 4]:.3f}")
print(f"C peak       σ 最大 {r[ic, 7]:.3f} @ 占比 {r[ic, 6]:.3f}")

print(f"\nB 是否等于 A 的 2 倍（4 个独立 shot 相加 ⇒ σ 应 ×2）：")
print(f"  B_σ / A_σ 各档：{np.round(r[:, 5] / np.maximum(r[:, 2], 1e-9), 3)}")

print(f"\nA 实测方差低于二项的部分 ⇒ 归因于 f_pix 不均匀，隐含 Var_t(p_t)：")
print(f"  {'boost':>8}{'p̄':>8}{'二项方差':>10}{'实测方差':>10}{'隐含Var_t(p)':>13}{'√Var_t/p̄':>10}")
for k in range(len(r)):
    vb = N_PIX * r[k, 1] * (1 - r[k, 1])
    vt = max(vb - r[k, 8], 0.0) / N_PIX
    print(f"  {r[k, 0]:>8.4g}{r[k, 1]:>8.3f}{vb:>10.3f}{r[k, 8]:>10.3f}{vt:>13.4f}"
          f"{np.sqrt(vt) / max(r[k, 1], 1e-9):>10.2f}")
