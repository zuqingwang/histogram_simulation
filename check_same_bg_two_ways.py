# -*- coding: utf-8 -*-
"""同一个 bg，两种造法，peak 到底谁大？—— 用项目自己的引擎直接 MC 验证。

对照组：
  A: noise=1 × N=4 shots  → bg = 4
  B: noise=4 × N=1 shot   → bg = 4

顺带回答「4 发的峰不会落在同一个 bin」这个直觉：
  - 打印每发 hist_i 各自的 peak、以及它们的 bin 位置是否重合
  - 打印「假如 4 发的峰完美对齐」的上界，与真实 hist_add peak 对比

用法：
    $env:PYTHONIOENCODING="utf-8"
    python check_same_bg_two_ways.py --n-mc 20000
"""
from __future__ import annotations

import argparse
import os
import time

os.environ.setdefault("POD_CORE_QUIET", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np

import pod_esti_v11_core as core

BG = 4.0
I0, I1 = core.I_STAT0, core.I_STAT1


def run_case(noise_per_shot, n_shots, n_mc, chunk, seed0, tag):
    """返回 hist_add 的 peak 样本 + 每发 hist_i 的 peak / 峰位。"""
    r_det = float(core.r_det_for_noise(float(noise_per_shot), core.N_PIX_MACRO))
    inv_tab = core.build_inv_table(r_det)
    pk_add, pk_shot, arg_shot, bg_mc = [], [], [], []
    done, part = 0, 0
    t0 = time.time()
    while done < n_mc:
        m = min(chunk, n_mc - done)
        rng = np.random.default_rng(seed0 + 7919 * part)
        hi = core.noise_hists_per_shot(m, n_shots, r_det, rng, inv_tab=inv_tab)
        win = hi[:, :, I0:I1]                 # (m, n_shots, N_STAT)
        add = win.sum(axis=1)                 # hist_add
        pk_add.append(add.max(axis=1))
        bg_mc.append(add.mean(axis=1))
        pk_shot.append(win.max(axis=2))       # (m, n_shots) 每发自己的 peak
        arg_shot.append(win.argmax(axis=2))   # (m, n_shots) 每发峰位
        done += m
        part += 1
        print(f"    [{tag}] {done}/{n_mc}  ({time.time()-t0:.0f}s)", flush=True)
    return (np.concatenate(pk_add), np.concatenate(pk_shot),
            np.concatenate(arg_shot), np.concatenate(bg_mc))


def summarize(pk, label):
    pk = np.asarray(pk, float)
    return (f"{label:<26} n={pk.size:>7d}  均值={pk.mean():7.4f}  "
            f"std={pk.std():6.4f}  中位={np.median(pk):5.1f}  "
            f"p99={np.percentile(pk,99):5.1f}  "
            f"p99.9={np.percentile(pk,99.9):5.1f}  max={pk.max():4.0f}")


def far_thr(pk, far):
    """满足 P(peak>=T)<=far 的最小整数 T。"""
    pk = np.asarray(pk, int)
    cnt = np.bincount(pk)
    n = cnt.sum()
    sf = np.concatenate([[1.0], 1.0 - np.cumsum(cnt) / n])  # sf[t]=P(X>=t)
    ok = np.where(sf <= far)[0]
    return int(ok[0]) if ok.size else -1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-mc", type=int, default=20000)
    ap.add_argument("--chunk", type=int, default=5000)
    args = ap.parse_args()

    print("=" * 78)
    print(f"同 bg = {BG:g} 的两种造法，统计窗 {core.N_STAT} bins，"
          f"宏像元 {core.N_PIX_MACRO} SPAD")
    print("=" * 78)

    print("\n[A] noise=1 × 4 shots")
    pkA, pkA_shot, argA, bgA = run_case(BG / 4, 4, args.n_mc, args.chunk, 101, "A")
    print("\n[B] noise=4 × 1 shot")
    pkB, pkB_shot, argB, bgB = run_case(BG, 1, args.n_mc, args.chunk, 202, "B")

    print("\n" + "-" * 78)
    print("① 实测 bg 校验（两边必须一致，否则对比无效）")
    print(f"  A: bg = {bgA.mean():.4f}   B: bg = {bgB.mean():.4f}")

    print("\n② hist_add 的 peak 分布")
    print("  " + summarize(pkA, "A: noise=1 ×4 shots"))
    print("  " + summarize(pkB, "B: noise=4 ×1 shot"))
    d = pkA.mean() - pkB.mean()
    print(f"  → A − B = {d:+.4f} 计数（{100*d/pkB.mean():+.2f}%）"
          f"　{'A 更大' if d > 0 else 'B 更大'}")

    print("\n③ 单 bin 计数的均值/方差（同 bg 下形状是否相同）")
    for tag, ns, npershot in [("A", 4, BG / 4), ("B", 1, BG)]:
        n_tr = core.N_PIX_MACRO * ns
        p = npershot / core.N_PIX_MACRO
        v = n_tr * p * (1 - p)
        print(f"  {tag}: Binomial(n_tr={n_tr:>3d}, p={p:.4f}) → "
              f"均值={n_tr*p:.3f}  方差={v:.4f}  σ={np.sqrt(v):.4f}  "
              f"Fano={1-p:.4f}")

    print("\n④ FAR 阈值（由本次 MC 直接反解）")
    for far, lab in [(0.05, "5%"), (0.01, "1%"), (0.001, "0.1%")]:
        ta, tb = far_thr(pkA, far), far_thr(pkB, far)
        print(f"  FAR={lab:>5}:  T(A, 4shot)={ta:>3d}   T(B, 1shot)={tb:>3d}   "
              f"差={ta-tb:+d}")

    print("\n" + "-" * 78)
    print("⑤ 回应「4 发的峰不会落在同一个 bin」")
    print(f"  A 组每发 hist_i 自己的 peak：均值 {pkA_shot.mean():.3f}，"
          f"（noise=1 的单发峰）")
    aligned = pkA_shot.sum(axis=1)
    print(f"  若 4 发的峰【完美对齐】，peak 会是 {aligned.mean():.3f}（上界，不可能达到）")
    print(f"  实际 hist_add 的 peak = {pkA.mean():.3f}"
          f"　→ 只有对齐上界的 {100*pkA.mean()/aligned.mean():.1f}%")
    same = (argA == argA[:, [0]]).all(axis=1).mean()
    print(f"  4 发峰位完全重合的比例 = {100*same:.3f}%（确实几乎不重合）")
    print(f"  但 B 组（noise=4 单发）的 peak = {pkB.mean():.3f}，"
          f"依然比 A 组的 {pkA.mean():.3f} 小")
    print("  ⇒ 「峰不对齐」是对的，但它不是决定胜负的量：")
    print("     hist_add 的 peak 不是由各发的 peak 叠出来的，")
    print("     而是【先逐 bin 求和、再取最大】。逐 bin 求和后均值都是 bg，")
    print("     谁的方差大谁的 peak 大。")

    # ---------------- 作图 ----------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import rcParams
    rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
    rcParams["axes.unicode_minus"] = False

    fig, ax = plt.subplots(1, 2, figsize=(13.0, 4.8))

    hi_max = int(max(pkA.max(), pkB.max())) + 2
    bins = np.arange(-0.5, hi_max + 0.5)
    ax[0].hist(pkA, bins=bins, density=True, alpha=0.55, color="tab:red",
               label=f"A: noise=1 ×4 shots（μ={pkA.mean():.2f}, σ={pkA.std():.2f}）")
    ax[0].hist(pkB, bins=bins, density=True, alpha=0.55, color="tab:blue",
               label=f"B: noise=4 ×1 shot（μ={pkB.mean():.2f}, σ={pkB.std():.2f}）")
    ax[0].axvline(BG, color="k", ls=":", lw=1.4, label=f"bg = {BG:g}（两组相同）")
    ax[0].set_xlabel("hist_add 的 peak [计数]")
    ax[0].set_ylabel("概率密度")
    ax[0].set_title(f"同 bg={BG:g}：4 发低噪的 peak 更大、尾更重")
    ax[0].legend(fontsize=8.5)
    ax[0].grid(alpha=0.3)

    ax[1].semilogy(np.arange(hi_max + 1),
                   [np.mean(pkA >= t) for t in range(hi_max + 1)],
                   "-o", color="tab:red", ms=3.5, lw=1.6, label="A: noise=1 ×4")
    ax[1].semilogy(np.arange(hi_max + 1),
                   [np.mean(pkB >= t) for t in range(hi_max + 1)],
                   "-s", color="tab:blue", ms=3.5, lw=1.6, label="B: noise=4 ×1")
    for far, lab in [(0.05, "5%"), (0.01, "1%"), (0.001, "0.1%")]:
        ax[1].axhline(far, color="0.5", ls="--", lw=0.9)
        ax[1].text(0.3, far * 1.15, f"FAR={lab}", fontsize=7.5, color="0.35")
    ax[1].set_xlabel("阈值 T [计数]")
    ax[1].set_ylabel("P(peak ≥ T) = FAR")
    ax[1].set_title("生存函数：A 的尾整体在 B 右侧，故阈值更高")
    ax[1].set_ylim(max(1.0 / pkA.size, 1e-6), 1.5)
    ax[1].legend(fontsize=8.5)
    ax[1].grid(alpha=0.3, which="both")

    fig.suptitle(
        f"同一个 bg={BG:g} 的两种造法（{args.n_mc:,} 次 MC，统计窗 {core.N_STAT} bins）",
        fontsize=12)
    fig.tight_layout()
    fig.savefig("check_same_bg_two_ways.png", dpi=130, bbox_inches="tight")
    print("\n图已保存 → check_same_bg_two_ways.png")


if __name__ == "__main__":
    main()
