# -*- coding: utf-8 -*-
"""快速体检 v30 的三份缓存：噪声 / PoD / 信号，以及宏像元缓存。

用法：python inspect_v30_cache.py
"""
import os
import sys

os.environ.setdefault("POD_CORE_QUIET", "1")
sys.stdout.reconfigure(encoding="utf-8")

import numpy as np

import pod_esti_v30_core as core


def sec(t):
    print("\n" + "=" * 92)
    print(t)
    print("=" * 92)


sec("① 噪声缓存 pod_esti_v30_cache_noise.npz")
if not core.NOISE_RES:
    print("  [缺] 未载入。先跑 run_pod_v30_noise_scan.py")
else:
    for n in core.N_SHOTS_LIST:
        R = core.NOISE_RES[n]
        d = R["done"]
        print(f"  N={n}: 完成 {int(d.sum())}/{len(d)}   n_tr={R['n_tr']}   "
              f"hist_std[bg=0.25/6/12] = {R['hist_std'][0]:.3f} / "
              f"{R['hist_std'][23]:.3f} / {R['hist_std'][47]:.3f}")
    print(f"  主缓存 {os.path.getsize(core.CACHE_NOISE)/1e6:.1f} MB；"
          f"检查点残留：{os.path.exists(core.CACHE_NOISE_CKPT)}")

sec("② 阈值 THRESH（整数计数）")
if not core.THRESH:
    print("  [缺] 噪声缓存不完整，THRESH 未构建")
else:
    tags = ["T5pct", "T1pct", "T0p5pct", "T0p1pct", "T100ppm", "T10ppm", "T10pct"]
    hdr = f"{'bg':>6} |" + "".join(f"{t[1:]:>8}" for t in tags)
    for n in core.N_SHOTS_LIST:
        print(f"\n  --- N_shots={n}（n_tr={core.N_PIX_MACRO*n}）---")
        print(hdr)
        for k in (0, 3, 11, 23, 35, 47):
            row = f"{core.BG_GRID[k]:6.2f} |"
            for t in tags:
                row += f"{int(core.THRESH[n][t][k]):>8d}"
            print(row)

sec("③ PoD 缓存 pod_esti_v30_cache_pod.npz")
if os.path.exists(core.CACHE_POD):
    z = np.load(core.CACHE_POD, allow_pickle=True)
    res = z["res"].item()
    print(f"  共 {len(res)} 个 (N, bg) 档；文件 "
          f"{os.path.getsize(core.CACHE_POD)/1e6:.1f} MB")
    ok = sum(1 for v in res.values() if v.get("critical"))
    print(f"  其中有临界解的：{ok}")
    k0 = sorted(res)[0]
    print(f"  样例 {k0}：T_map = {res[k0].get('T_map')}")

    # 求根质量：验证 PoD 与目标等级的偏差。v20 的求根器会把 PoD=1.000 / 0.68 的
    # 点当成 PoD90 存进来，导致模块 7 曲线出现 3–5 倍尖刺，所以这里必须逐点体检。
    print(f"\n  --- 求根质量（|验证PoD − 目标| 应 ≤ {core.POD_VERIFY_TOL}）---")
    print(f"  {'N':>2} {'FAR':>5} {'level':>6} {'档数':>5} {'中位|err|':>9} "
          f"{'最大|err|':>9} {'超容差':>7} {'平均轮数':>8}")
    worst = []
    for n in core.N_SHOTS_LIST:
        for tag in core.POD_FAR_TAGS:
            for lv in ("0.50", "0.90"):
                errs, rounds = [], []
                for (nn, bg), v in res.items():
                    if nn != n:
                        continue
                    c = (v.get("critical") or {}).get(tag, {}).get(lv)
                    if not c:
                        continue
                    e = abs(c.get("pod_err", c["pod"] - float(lv)))
                    errs.append(e)
                    rounds.append(c.get("verify_rounds", 1))
                    if e > core.POD_VERIFY_TOL:
                        worst.append((e, n, tag, lv, bg, c["pod"], c["boost"]))
                if not errs:
                    continue
                errs = np.array(errs)
                bad = int((errs > core.POD_VERIFY_TOL).sum())
                print(f"  {n:>2} {core.FAR_TAG_TO_LABEL[tag]:>5} {lv:>6} {len(errs):>5} "
                      f"{np.median(errs):>9.4f} {errs.max():>9.4f} {bad:>7} "
                      f"{np.mean(rounds):>8.2f}")
    if worst:
        print(f"\n  [注意] 共 {len(worst)} 个点超容差，最差 5 个：")
        for e, n, tag, lv, bg, pod, b in sorted(worst, reverse=True)[:5]:
            print(f"    N={n} bg={bg:5.2f} {core.FAR_TAG_TO_LABEL[tag]:>5} PoD{lv}: "
                  f"验证PoD={pod:.3f} boost={b:.4g}")
    else:
        print("\n  [好] 全部临界点都收敛在容差内")
else:
    print("  [缺] 先跑 run_pod_v30_pod_scan.py")

sec("④ 信号缓存 pod_esti_v30_cache_signal.npz")
if os.path.exists(core.CACHE_SIG):
    z = np.load(core.CACHE_SIG)
    print(f"  n_mc={int(z['n_mc']):,}  boosts={np.asarray(z['boosts'])}")
    for n in core.N_SHOTS_LIST:
        c_ = z[f"peak_cnt_{n}"]
        d_ = z[f"done_{n}"] if f"done_{n}" in z.files else (c_.sum(axis=2) > 0)
        print(f"  N={n}: shape={c_.shape}  完成 {int(np.sum(d_))}/{d_.size}")
else:
    print("  [缺] 先跑 run_pod_v30_sig_scan.py")

sec("⑤ 宏像元缓存 compare_macro_v30_cache.npz")
MC = "compare_macro_v30_cache.npz"
if os.path.exists(MC):
    z = np.load(MC)
    print(f"  n_mc={int(z['n_mc']):,}  配置={[tuple(int(v) for v in c) for c in z['cfgs']]}")
    print(f"  完成 {int(np.sum(z['done']))}/{z['done'].size}")
else:
    print("  [缺] 先跑 compare_macro_v30.py")
print()
