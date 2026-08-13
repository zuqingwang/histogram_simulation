"""判定实验：peak 标准差拱顶偏离 0.5，到底是 f_pix 不均匀造成的，还是「取最大值」造成的？

做法：同一套物理、同一个总收集比例 Σf_pix，只改分配方式，跑两组对照
  REAL    —— 真实像斑轮廓（27 个 SPAD 的 f_pix 相差 240 倍）
  UNIFORM —— 27 个 SPAD 平分（每个 = Σf_pix / 27），此时各轨迹点亮概率严格相同

每组各测两个统计量
  SB   单个 bin 的计数（取该能量档自己的最亮 bin）。若各轨迹 p 相同，它是严格二项分布，
       方差 n_tr·p(1-p)，必定在占比 0.5 见顶——这是数学必然，不是经验问题。
  PEAK 152 个 bin 上的最大值。

四种组合的拱顶位置就能把两个机制拆开：
  UNIFORM 的 SB 顶在 0.5   → 校验代码与理论自洽
  UNIFORM 的 PEAK 顶在哪   → 「取最大值」单独贡献多少左移
  REAL 的 SB 顶在哪        → 「f_pix 不均匀」单独贡献多少左移
"""
from __future__ import annotations

import os
import sys
from concurrent.futures import ProcessPoolExecutor

os.environ.setdefault("POD_CORE_QUIET", "1")
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

import numpy as np

import pod_esti_v30_core as core
import run_peak_energy_scan as scan

N_MC = 6000
N_SHOTS = 4
N_PIX = core.N_PIX_MACRO
N_TR = N_PIX * N_SHOTS

# 必须显式取 F_PIX_REAL：scan.F_PIX 跟着 F_PIX_MODE 走，默认已是均匀，
# 用它会把本脚本变成「均匀 vs 均匀」的空对照。
F_REAL = scan.F_PIX_REAL
F_UNIF = np.full(N_PIX, F_REAL.sum() / N_PIX)      # 总量相同，只改分配

# 拱顶落在 boost 0.01–0.05，这一段要密
GRID = np.unique(np.concatenate([
    np.geomspace(0.002, 0.008, 5),
    np.geomspace(0.008, 0.06, 13),
    np.geomspace(0.06, 0.5, 7),
    np.geomspace(0.5, 20.0, 6),
]).round(8))


def _one(args):
    tag, f_arr, b = args
    rng = np.random.default_rng(20260810)
    h1, h4, pk = [], [], []
    done = 0
    while done < N_MC:
        m = min(2000, N_MC - done)
        hi = core.binary_macro_stepping_per_shot(
            m, f_arr, N_SHOTS, core.R_SIG_UNIT_GEN, core.TF_GEN,
            0.0, core.CENTERS, rng, boost=float(b))
        add = hi.sum(axis=1)[:, core.I_STAT0:core.I_STAT1]
        h1.append(hi[:, 0, core.I_STAT0:core.I_STAT1])
        h4.append(add)
        pk.append(add.max(axis=1))
        done += m
    h1 = np.concatenate(h1).astype(float)
    h4 = np.concatenate(h4).astype(float)
    pk = np.concatenate(pk).astype(float)

    j = int(h4.mean(axis=0).argmax())              # 该档自己的最亮 bin
    sb1, sb4 = h1[:, j], h4[:, j]
    pbar = sb1.mean() / N_PIX                      # 单 shot 单 bin 的占比
    return dict(tag=tag, boost=b,
                sb_p=sb4.mean() / N_TR, sb_sd=sb4.std(),
                sb_binom=np.sqrt(N_TR * (sb4.mean() / N_TR) * (1 - sb4.mean() / N_TR)),
                pk_p=pk.mean() / N_TR, pk_sd=pk.std(),
                pbar1=pbar)


def _arch(rows, key_sd, key_p):
    i = int(np.argmax([r[key_sd] for r in rows]))
    return rows[i][key_sd], rows[i][key_p], rows[i]["boost"]


if __name__ == "__main__":
    print("=" * 100)
    print("判定实验：peak 标准差拱顶左移，是 f_pix 不均匀、还是「取最大值」造成的？")
    print("=" * 100)
    print(f"  REAL    f_pix: 最大 {F_REAL.max():.4g} 最小 {F_REAL.min():.4g} "
          f"离散度 std/mean = {F_REAL.std() / F_REAL.mean():.3f}")
    print(f"  UNIFORM f_pix: 全部 = {F_UNIF[0]:.4g}          离散度 std/mean = 0.000")
    print(f"  两组 Σf_pix 相同（{F_REAL.sum():.4f} vs {F_UNIF.sum():.4f}），"
          f"n_tr = {N_TR}，每档 {N_MC:,} 次 MC，bg = 0")
    print(f"  能量网格 {GRID.size} 档：boost {GRID[0]:.4g} → {GRID[-1]:.4g}\n")

    jobs = ([("REAL", F_REAL, b) for b in GRID]
            + [("UNIFORM", F_UNIF, b) for b in GRID])
    with ProcessPoolExecutor(max_workers=20) as ex:
        res = list(ex.map(_one, jobs))

    out = {}
    for tag in ("REAL", "UNIFORM"):
        out[tag] = sorted([r for r in res if r["tag"] == tag], key=lambda r: r["boost"])

    for tag in ("REAL", "UNIFORM"):
        print("=" * 100)
        print(f"【{tag}】  SB = N=4 单 bin 计数　PEAK = 152 bin 取最大")
        print(f"  {'boost':>9} | {'SB占比':>7}{'SB实测σ':>9}{'SB二项σ':>9}{'比值':>7} | "
              f"{'PEAK占比':>9}{'PEAK σ':>8} | {'PEAK/SB σ':>10}")
        for r in out[tag]:
            print(f"  {r['boost']:>9.4g} | {r['sb_p']:>7.3f}{r['sb_sd']:>9.3f}"
                  f"{r['sb_binom']:>9.3f}{r['sb_sd'] / max(r['sb_binom'], 1e-9):>7.3f} | "
                  f"{r['pk_p']:>9.3f}{r['pk_sd']:>8.3f} | "
                  f"{r['pk_sd'] / max(r['sb_sd'], 1e-9):>10.3f}")

    print("\n" + "=" * 100)
    print("拱顶位置汇总（关键结论）")
    print("=" * 100)
    print(f"  {'配置':>9} {'统计量':>6} {'拱顶σ':>9} {'拱顶处占比':>11} {'此处boost':>11}")
    for tag in ("REAL", "UNIFORM"):
        for lbl, ks, kp in (("SB", "sb_sd", "sb_p"), ("PEAK", "pk_sd", "pk_p")):
            sd, p, b = _arch(out[tag], ks, kp)
            print(f"  {tag:>9} {lbl:>6} {sd:>9.3f} {p:>11.3f} {b:>11.4g}")
    print(f"\n  纯二项理论：拱顶应在占比 0.500，峰值 √(n_tr/4) = {np.sqrt(N_TR * 0.25):.3f}")
