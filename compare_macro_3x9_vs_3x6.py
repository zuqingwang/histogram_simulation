# -*- coding: utf-8 -*-
"""宏像元 3×9（27 SPAD）与 3×6（18 SPAD）的阈值对比。

物理设定
--------
* **环境光在 SPAD 上均匀分布**：每个 SPAD 的探测率 r_det 相同，与宏像元多大无关。
  于是宏像元每 bin 的底噪 bg = (SPAD 数 × N_shots) × p_eq，其中 p_eq 是
  **单个 SPAD 单发**在一个 1 ns bin 内被点亮的平衡态概率。
* **信号同样按均匀处理**：每个 SPAD 收到的信号强度相同，宏像元收到的信号 ∝ n_pix。
  （实际像斑在 x 方向很窄，但本次分析按均匀算。）
* **纯噪声时引擎把「SPAD 数」和「shot 数」折进同一个维度** n_tr = n_pix × N_shots
  （见 core.noise_macro_hist_fast 的 docstring）。所以在同一环境光下，
  **阈值只是 n_tr 的函数**：3×6 跑 6 发 与 3×9 跑 4 发（都是 n_tr=108）完全等价。

在「噪声均匀 + 信号均匀」下，信号要跨过阈值所需的**每 SPAD 每发的额外点亮概率**为
    q_req = (T − bg) / (n_pix · N) = T/n_tr − p_eq
它同样只是 n_tr 与 p_eq 的函数，是本脚本的灵敏度判据（越小越灵敏）。

扫描自变量与横轴
----------------
本模型里环境光**完全由 p_eq 一个数刻画**（单个 SPAD、单发、单个 1 ns bin 被点亮的
平衡态概率），它与宏像元多大、累加多少发都无关。扫描就沿 p_eq 打网格。

作图时用两种横轴，各自回答不同问题：
* **图①：横轴 = 环境光照度 [klux]**（由 p_eq 反解 E_lambda 换算而来）。
  回答「同一片天光下，各配置的底噪分别是多少」——bg = n_tr · p_eq，正比于 n_tr。
* **图②–⑥：横轴 = 该配置自己的 bg**。
  回答「工作在同一个底噪水平上时，各配置的阈值/灵敏度如何」。
  注意同一个 bg 对不同配置意味着**不同的环境光**（n_tr 大的配置要更弱的天光才能到同样 bg），
  且各配置的 bg 覆盖范围不同（bg 上限 = n_tr · p_eq_max）。

缓存 + 多进程（项目规则三）。用法：
    $env:PYTHONIOENCODING="utf-8"
    python compare_macro_3x9_vs_3x6.py --workers 20 --n-mc 300000
    python compare_macro_3x9_vs_3x6.py --limit 3 --n-mc 5000    # 冒烟
"""
from __future__ import annotations

import argparse
import functools
import os
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

print = functools.partial(print, flush=True)  # noqa: A001

os.environ.setdefault("POD_CORE_QUIET", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import numpy as np

import pod_esti_v11_core as core

CACHE = "compare_macro_3x9_vs_3x6_cache.npz"
CACHE_CKPT = "compare_macro_3x9_vs_3x6_cache.partial.npz"

# (x 方向 SPAD 数, y 方向 SPAD 数, N_shots)
CONFIGS = [
    (9, 3, 1), (9, 3, 2), (9, 3, 4),
    (6, 3, 1), (6, 3, 2), (6, 3, 4), (6, 3, 6),
]
CFG_COLOR = {
    (9, 3, 1): "#9ecae1", (9, 3, 2): "#4292c6", (9, 3, 4): "#08519c",
    (6, 3, 1): "#fcae91", (6, 3, 2): "#fb6a4a", (6, 3, 4): "#cb181d",
    (6, 3, 6): "#67000d",
}
CFG_LS = {9: "-", 6: "--"}

FAR_MAIN = 0.01           # 主 FAR（1%）
FAR_KEYS = [0.05, 0.01, 0.001, 100e-6]

# ---- 环境光扫描档位 ----
# 环境光只由 p_eq（单 SPAD／单发／单 bin 的平衡态点亮概率）决定。为了让档位落在整齐的数上，
# 用一个【参考宏像元】（N_PIX_REF 个 SPAD、单发）的每 bin 底噪 ref = N_PIX_REF·p_eq 来打网格。
# 这个参考量只是"给档位起个整数名字"的标尺，不是任何配置的物理量。
N_PIX_REF = 27                                  # 参考宏像元 = 现行 3×9
AMB_REF_STEP = 0.5                              # 参考底噪步长
AMB_REF_MAX = 12.0                              # 参考底噪上限 ⇒ p_eq 上限 = 12/27 = 0.4444
E_LAMBDA_100KLUX = 0.68                         # W/m²/nm，PARAMS 基准（≈100 klux）


def label(cfg):
    nx, ny, n = cfg
    return f"{ny}×{nx}（{nx*ny} SPAD） N={n}"


# ---------------------------------------------------------------- 信号权重
def sig_weight(cfg):
    """信号按均匀处理：宏像元累加后收到的信号总量 ∝ n_pix × N_shots = n_tr。"""
    nx, ny, n = cfg
    return float(nx * ny * n)


# ---------------------------------------------------------------- MC 单档
def _job(a):
    ci, ai, nx, ny, n_shots, amb_ref, n_mc, chunk, seed0 = a
    n_pix = nx * ny
    n_tr = n_pix * n_shots
    # 环境光只由 p_eq = amb_ref / N_PIX_REF 决定；反解成单 SPAD 探测率 r_det，
    # 各配置共用同一个 r_det，保证它们看的是同一片天光
    r_det = float(core.r_det_for_noise(amb_ref, N_PIX_REF))
    inv_tab = core.build_inv_table(r_det)
    i0, i1 = core.I_STAT0, core.I_STAT1

    peak_cnt = np.zeros(n_tr + 2, dtype=np.int64)
    bg_sum = 0.0
    done, part = 0, 0
    while done < n_mc:
        m = min(chunk, n_mc - done)
        rng = np.random.default_rng(seed0 + 7919 * part)
        h = core.noise_macro_hist_fast(m, n_tr, r_det, rng, inv_tab=inv_tab)
        a_ = h[:, i0:i1]
        bg_sum += float(a_.mean(axis=1).sum())
        peak_cnt += np.bincount(a_.max(axis=1).astype(np.int64),
                                minlength=peak_cnt.size)
        done += m
        part += 1

    v = np.arange(peak_cnt.size, dtype=float)
    tot = peak_cnt.sum()
    mu = float((v * peak_cnt).sum() / tot)
    sd = float(np.sqrt(max((v * v * peak_cnt).sum() / tot - mu * mu, 0.0)))
    thr = {f: float(core.far_threshold_from_cnt(peak_cnt, f)[0]) for f in FAR_KEYS}

    return dict(ci=ci, ai=ai, n_tr=n_tr, amb_ref=float(amb_ref),
                p_eq=float(amb_ref) / N_PIX_REF,
                r_det=r_det, e_lambda=float(core.e_lambda_for_r_det(r_det)),
                bg_mc=bg_sum / tot, peak_mean=mu, peak_std=sd, thr=thr)


# ---------------------------------------------------------------- 缓存
_FIELDS = ["bg_mc", "peak_mean", "peak_std", "e_lambda"] + [f"thr_{i}" for i in range(len(FAR_KEYS))]


def _empty(nc, na):
    d = {k: np.zeros((nc, na)) for k in _FIELDS}
    d["done"] = np.zeros((nc, na), dtype=bool)
    return d


def _save(path, res, amb, n_mc):
    tmp = path + ".tmp.npz"
    np.savez_compressed(tmp, amb=amb, n_mc=n_mc,
                        cfgs=np.asarray(CONFIGS), fars=np.asarray(FAR_KEYS),
                        **res)
    os.replace(tmp, path)


def _load(path, amb, n_mc):
    if not os.path.exists(path):
        return None
    try:
        z = np.load(path)
        if (int(z["n_mc"]) != int(n_mc) or not np.array_equal(z["cfgs"], np.asarray(CONFIGS))
                or z["amb"].shape != amb.shape or not np.allclose(z["amb"], amb)):
            return None
        out = {k: np.array(z[k]) for k in _FIELDS}
        out["done"] = np.array(z["done"])
        return out
    except Exception:
        return None


# ---------------------------------------------------------------- 主流程
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=20)
    ap.add_argument("--n-mc", type=int, default=300_000)
    ap.add_argument("--chunk", type=int, default=5_000)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--checkpoint-every", type=int, default=10)
    args = ap.parse_args()

    # amb 数组 = 参考底噪 N_PIX_REF·p_eq，同时充当缓存键
    amb = np.round(np.arange(AMB_REF_STEP, AMB_REF_MAX + 1e-9, AMB_REF_STEP), 4)
    na, nc = len(amb), len(CONFIGS)

    print("=" * 88)
    print("环境光只由 p_eq 刻画 = 单 SPAD／单发／单个 1 ns bin 的平衡态点亮概率")
    print(f"  扫描 {na} 档：p_eq = {amb[0]/N_PIX_REF:.4f} → {amb[-1]/N_PIX_REF:.4f}"
          f"（用 {N_PIX_REF}·p_eq = {amb[0]:g}→{amb[-1]:g} 打整齐网格）")
    print("噪声每 SPAD 均匀 ⇒ bg = n_tr · p_eq；信号也按每 SPAD 均匀 ⇒ 收到的信号 ∝ n_tr")
    print("灵敏度判据 q_req = (T − bg)/n_tr = 每 SPAD 每发需额外贡献的点亮概率（越小越灵敏）")
    print("配置与轨迹数 n_tr = n_pix × N_shots：")
    for cfg in CONFIGS:
        nx, ny, n = cfg
        ntr = nx * ny * n
        print(f"  {label(cfg):>24}  n_tr = {nx*ny}×{n} = {ntr:>3}"
              f"   bg = {ntr}·p_eq，覆盖 bg ∈ [{ntr*amb[0]/N_PIX_REF:.2f},"
              f" {ntr*amb[-1]/N_PIX_REF:.1f}]")
    print("=" * 88)

    res = _load(CACHE, amb, args.n_mc) or _load(CACHE_CKPT, amb, args.n_mc)
    if res is None:
        res = _empty(nc, na)
        print("未找到缓存，从零开始")
    else:
        print(f"命中缓存，已完成 {int(res['done'].sum())}/{nc*na} 档")

    todo = [(ci, ai, cfg[0], cfg[1], cfg[2], float(amb[ai]),
             args.n_mc, args.chunk, 91000 + 1373 * ci + 17 * ai)
            for ci, cfg in enumerate(CONFIGS) for ai in range(na)
            if not res["done"][ci, ai] and not (args.limit and ai >= args.limit)]
    # 单档耗时 ≈ n_tr × 环境光强度，重活先派发可显著缩短尾部空转
    todo.sort(key=lambda j: -(j[2] * j[3] * j[4]) * j[5])
    print(f"环境光 {na} 档 × 配置 {nc} 种 × {args.n_mc:,} MC；待算 {len(todo)}；"
          f"workers={args.workers}")

    if todo:
        t0 = time.time()
        dn = 0
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            futs = [pool.submit(_job, j) for j in todo]
            for fut in as_completed(futs):
                r = fut.result()
                ci, ai = r["ci"], r["ai"]
                res["bg_mc"][ci, ai] = r["bg_mc"]
                res["peak_mean"][ci, ai] = r["peak_mean"]
                res["peak_std"][ci, ai] = r["peak_std"]
                res["e_lambda"][ci, ai] = r["e_lambda"]
                for i, f in enumerate(FAR_KEYS):
                    res[f"thr_{i}"][ci, ai] = r["thr"][f]
                res["done"][ci, ai] = True
                dn += 1
                el = time.time() - t0
                eta = el / dn * (len(todo) - dn)
                print(f"  [{dn}/{len(todo)} {100*dn/len(todo):5.1f}%] "
                      f"{label(CONFIGS[ci]):>24} p_eq={r['p_eq']:.4f} → "
                      f"bg={r['bg_mc']:6.3f} peakμ={r['peak_mean']:6.2f} "
                      f"T@1%={r['thr'][0.01]:5.0f}"
                      f"　已用 {el/60:.1f} min，剩 {eta/60:.1f} min")
                if dn % args.checkpoint_every == 0 or dn == len(todo):
                    _save(CACHE_CKPT, res, amb, args.n_mc)
        _save(CACHE, res, amb, args.n_mc)
        if os.path.exists(CACHE_CKPT):
            try:
                os.remove(CACHE_CKPT)
            except OSError:
                pass
        print(f"[扫描完成] → {CACHE}，{(time.time()-t0)/60:.1f} min")

    report(res, amb, args.n_mc)


# ---------------------------------------------------------------- 报告与作图
def _curve_vs_bg(cfg, T, bg, ok):
    """把某个配置的曲线改用【它自己的 bg】做自变量，返回 (bg, q_req, T)，bg 单调升序。

    同一环境光下 bg = n_tr·p_eq 随配置不同而不同，所以要比较「同一 bg 下」的表现，
    必须先各自换成以 bg 为自变量，再插值到公共网格。
    """
    ci = CONFIGS.index(cfg)
    m = ok[ci]
    x = bg[ci][m]
    q = (T[ci][m] - x) / sig_weight(cfg)
    o = np.argsort(x)
    return x[o], q[o], T[ci][m][o]


def report(res, amb, n_mc):
    im = FAR_KEYS.index(FAR_MAIN)
    T = res[f"thr_{im}"]
    bg = res["bg_mc"]
    ok = res["done"]

    # ---- 等价性核对：3×6 N=6 与 3×9 N=4（都是 n_tr=108） ----
    i_96 = CONFIGS.index((6, 3, 6))
    i_94 = CONFIGS.index((9, 3, 4))
    m = ok[i_96] & ok[i_94]
    if m.any():
        d = np.abs(T[i_96][m] - T[i_94][m])
        print(f"\n[等价性核对] n_tr=108 的两种实现（3×6@N=6 vs 3×9@N=4）："
              f"T@1% 最大差 {d.max():.0f} 计数，bg 最大差 "
              f"{np.abs(bg[i_96][m]-bg[i_94][m]).max():.4f}"
              f"（应为 0 / MC 噪声量级 → 证实阈值只取决于 n_tr）")

    # ---- 主对比表 ----
    ia = CONFIGS.index((6, 3, 2))
    ib = CONFIGS.index((9, 3, 4))
    print(f"\n{'='*100}")
    print("主对比 A（同一片天光）：3×6 @ N=2（n_tr=36） vs 3×9 @ N=4（n_tr=108），FAR=1%")
    print(f"{'p_eq':>7} {'klux':>6} | {'3×6 N=2':^26} | {'3×9 N=4':^26} | {'所需信号比':>10}")
    print(f"{'':>7} {'':>6} | {'bg':>7}{'T':>6}{'T-bg':>7}{'T/bg':>6} |"
          f" {'bg':>7}{'T':>6}{'T-bg':>7}{'T/bg':>6} | {'(6,2)/(9,4)':>11}")
    for k in range(len(amb)):
        if not (ok[ia, k] and ok[ib, k]):
            continue
        sa = (T[ia, k] - bg[ia, k]) / sig_weight(CONFIGS[ia])
        sb = (T[ib, k] - bg[ib, k]) / sig_weight(CONFIGS[ib])
        klux = res["e_lambda"][ib, k] / E_LAMBDA_100KLUX * 100
        print(f"{amb[k]/N_PIX_REF:7.4f} {klux:6.0f} | {bg[ia,k]:7.3f}{T[ia,k]:6.0f}"
              f"{T[ia,k]-bg[ia,k]:7.3f}{T[ia,k]/bg[ia,k]:6.2f} |"
              f" {bg[ib,k]:7.3f}{T[ib,k]:6.0f}{T[ib,k]-bg[ib,k]:7.3f}"
              f"{T[ib,k]/bg[ib,k]:6.2f} | {sa/sb:11.3f}")

    # ---- 同 N 下 3×6 / 3×9（同一片天光） ----
    print(f"\n{'='*100}")
    print("同一片天光、同 N_shots 下 3×6 相对 3×9（FAR=1%）：bg 恒为 2/3，看 T 与所需信号")
    print(f"{'p_eq':>7} {'klux':>6} | "
          + " | ".join(f"N={n}: T6  T9  ΔT  信号比" for n in (1, 2, 4)))
    for k in range(0, len(amb), 2):
        cells = []
        good = True
        for n in (1, 2, 4):
            i6, i9 = CONFIGS.index((6, 3, n)), CONFIGS.index((9, 3, n))
            if not (ok[i6, k] and ok[i9, k]):
                good = False
                break
            s6 = (T[i6, k] - bg[i6, k]) / sig_weight((6, 3, n))
            s9 = (T[i9, k] - bg[i9, k]) / sig_weight((9, 3, n))
            cells.append(f"{T[i6,k]:3.0f} {T[i9,k]:3.0f} {T[i6,k]-T[i9,k]:4.0f} {s6/s9:6.3f}")
        if good:
            klux = res["e_lambda"][i9, k] / E_LAMBDA_100KLUX * 100
            print(f"{amb[k]/N_PIX_REF:7.4f} {klux:6.0f} | " + " | ".join(cells))

    # ---- 主对比 B：同一个 bg（各配置对应不同天光）----
    print(f"\n{'='*100}")
    print("主对比 B（同一个 bg，各配置对应不同天光）：把各配置的曲线插到公共 bg 网格上，FAR=1%")
    print(f"{'bg':>6} | " + " | ".join(f"{label(c):>22}" for c in
                                       ((6, 3, 2), (9, 3, 4))) + " | 所需信号比")
    print(f"{'':>6} | " + " | ".join(f"{'T':>7}{'T-bg':>7}{'q_req':>8}"
                                     for _ in range(2)) + " | (6,2)/(9,4)")
    bga, qa, Ta = _curve_vs_bg((6, 3, 2), T, bg, ok)
    bgb, qb, Tb = _curve_vs_bg((9, 3, 4), T, bg, ok)
    for x in (2, 4, 6, 8, 10, 12, 14, 16):
        if not (bga.min() <= x <= bga.max() and bgb.min() <= x <= bgb.max()):
            continue
        ta, qa_ = np.interp(x, bga, Ta), np.interp(x, bga, qa)
        tb, qb_ = np.interp(x, bgb, Tb), np.interp(x, bgb, qb)
        print(f"{x:6.1f} | {ta:7.2f}{ta-x:7.2f}{qa_:8.4f} |"
              f" {tb:7.2f}{tb-x:7.2f}{qb_:8.4f} | {qa_/qb_:11.3f}")

    # ---- 作图 ----
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import rcParams
    rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
    rcParams["axes.unicode_minus"] = False

    klux = res["e_lambda"].max(axis=0) / E_LAMBDA_100KLUX * 100   # 各配置相同，取任一行
    XLAB_BG = "bg = 该配置 hist_add 统计窗每 bin 均值 [计数]"

    fig, axes = plt.subplots(2, 3, figsize=(17.5, 9.4))
    ax = axes.ravel()

    # ---- 图①：横轴 = 环境光照度，回答「同一片天光下各配置底噪多少」----
    a = ax[0]
    for cfg in CONFIGS:
        ci = CONFIGS.index(cfg)
        msk = ok[ci]
        if not msk.any():
            continue
        a.plot(klux[msk], bg[ci][msk], CFG_LS[cfg[0]], color=CFG_COLOR[cfg],
               marker="o", ms=2.4, lw=1.7, label=label(cfg))
    a.set_xlabel("环境光照度 [klux]（等价于 p_eq = 单 SPAD 单发每 bin 点亮概率）")
    a.set_ylabel("bg [计数/bin]")
    a.set_title("① 同一片天光下的底噪 bg = n_tr·p_eq（噪声每 SPAD 均匀，纯按 n_tr 缩放）",
                fontsize=10.5)
    a.grid(alpha=0.3)
    a.legend(fontsize=7.0)

    # ---- 图②–⑤：横轴 = 该配置自己的 bg ----
    def curves_vs_bg(a, yfun, ylab, title):
        for cfg in CONFIGS:
            ci = CONFIGS.index(cfg)
            msk = ok[ci]
            if not msk.any():
                continue
            a.plot(bg[ci][msk], yfun(ci, msk), CFG_LS[cfg[0]], color=CFG_COLOR[cfg],
                   marker="o", ms=2.4, lw=1.7, label=label(cfg))
        a.set_xscale("log")
        a.set_xlabel(XLAB_BG)
        a.set_ylabel(ylab)
        a.set_title(title, fontsize=10.5)
        a.grid(alpha=0.3, which="both")
        a.legend(fontsize=7.0)

    curves_vs_bg(ax[1], lambda ci, m: T[ci][m], "T@FAR=1% [计数]",
                 "② 阈值 T vs 自身 bg（n_tr 越大，同 bg 下阈值越高）")
    curves_vs_bg(ax[2], lambda ci, m: T[ci][m] - bg[ci][m], "T − bg [计数]",
                 "③ 阈值余量 T − bg（信号必须填上的部分）")
    curves_vs_bg(ax[3], lambda ci, m: T[ci][m] / np.maximum(bg[ci][m], 1e-9), "T / bg",
                 "④ 阈值/底噪比")
    curves_vs_bg(ax[4],
                 lambda ci, m: (T[ci][m] - bg[ci][m]) / sig_weight(CONFIGS[ci]),
                 "q_req = (T−bg)/n_tr   [每 SPAD 每发的点亮概率]",
                 "⑤ 所需信号 q_req（信号也均匀 → 收到的信号 ∝ n_tr）")

    # ---- 图⑥：成对比值，插值到公共 bg 网格 ----
    a = ax[5]
    pairs = [((6, 3, 2), (9, 3, 4)), ((6, 3, 4), (9, 3, 4)),
             ((6, 3, 2), (9, 3, 2)), ((6, 3, 6), (9, 3, 4))]
    for (ca, cb), col in zip(pairs, ["tab:red", "tab:orange", "tab:green", "tab:purple"]):
        bga, qa, _ = _curve_vs_bg(ca, T, bg, ok)
        bgb, qb, _ = _curve_vs_bg(cb, T, bg, ok)
        if bga.size < 2 or bgb.size < 2:
            continue
        lo = max(bga.min(), bgb.min())
        hi = min(bga.max(), bgb.max())
        if not (hi > lo):
            continue
        xs = np.geomspace(lo, hi, 80)
        a.plot(xs, np.interp(xs, bga, qa) / np.interp(xs, bgb, qb), "-", lw=1.9,
               color=col, label=f"{label(ca)}  ÷  {label(cb)}")
    a.axhline(1.0, color="k", lw=1.0, ls=":")
    a.set_xscale("log")
    a.set_xlabel("公共 bg（两个配置都工作在这个底噪上）[计数/bin]")
    a.set_ylabel("所需信号强度之比（<1 表示前者更灵敏）")
    a.set_title("⑥ 成对比较：同一 bg 下折算到相同探测能力所需的信号", fontsize=10.5)
    a.grid(alpha=0.3, which="both")
    a.legend(fontsize=7.2)

    fig.suptitle(
        f"宏像元 3×9（27 SPAD）vs 3×6（18 SPAD）阈值对比　FAR=1%，每档 {n_mc:,} 次纯噪声 MC\n"
        f"噪声与信号均按每 SPAD 均匀 → bg = n_tr·p_eq、收到的信号 ∝ n_tr；"
        f"阈值只取决于 n_tr = n_pix×N_shots（故 3×6@N=6 与 3×9@N=4 两条曲线完全重合）\n"
        f"图① 横轴 = 环境光照度（同一片天光）；图②–⑥ 横轴 = 各配置自身的 bg"
        f"（同一 bg 对应不同天光）",
        fontsize=12)
    fig.tight_layout()
    fig.savefig("compare_macro_3x9_vs_3x6.png", dpi=130, bbox_inches="tight")
    print("\n图已保存 → compare_macro_3x9_vs_3x6.png")


if __name__ == "__main__":
    main()
