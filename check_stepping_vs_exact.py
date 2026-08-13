# -*- coding: utf-8 -*-
"""校验：快速步进引擎 `binary_macro_stepping_per_shot` 与精确引擎 `spad_binary_trace`
在**各档信号能量**下给出一致的累加波形，重点看饱和区能否压出削顶平台。

背景：v30 原实现只保留「最近一次雪崩时刻」`tcov`，判据 `0 <= c - tcov < T_OVER`。
饱和时 SPAD 在同一 bin 内连续雪崩，`tcov` 被改写成带抖动的更晚时刻 ⇒ `c - tcov < 0`
⇒ 本该被前一次雪崩 8 ns 窗口点亮的 bin 被判成灭的，平台压不出来。
修法见 `pod_esti_v30_core._cover_commit`。

精确引擎是逐光子 Python 循环，boost 越大越慢，所以只跑到 boost=50。
    python -u check_stepping_vs_exact.py
"""
from __future__ import annotations

import os
import sys
import time

os.environ.setdefault("POD_CORE_QUIET", "1")
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

import numpy as np

import pod_esti_v30_core as core

N_REAL = 60
BOOSTS = (0.002, 0.01, 0.05, 0.2, 0.5, 2.0, 10.0, 50.0)

F_REAL = core.FPIX[:, core.M_PEAK * core.MACRO_BY:(core.M_PEAK + 1) * core.MACRO_BY].ravel()
F_UNIF = np.full(core.N_PIX_MACRO, F_REAL.sum() / core.N_PIX_MACRO)

print("=" * 104)
print(f"快速步进引擎 vs 精确引擎　N_shots=1，n_tr={core.N_PIX_MACRO}，bg=0，"
      f"每档 {N_REAL} 次实现")
print(f"T_OVER = {core.T_OVER * 1e9:.3f} ns，抖动 = {core.JIT * 1e12:.0f} ps，"
      f"细网格步长 = {(core.TF_GEN[1] - core.TF_GEN[0]) * 1e12:.0f} ps")
print("=" * 104)

for tag, F in (("uniform", F_UNIF), ("real", F_REAL)):
    print(f"\n--- f_pix = {tag}（Σf = {F.sum():.5f}）---")
    print(f"  {'boost':>8} │ {'平均波形峰值':>12} {'峰后10ns平台':>12} {'打满27的bin数':>13} "
          f"{'peak均值':>9} │ 同左（精确） │ {'平台相对差':>10}")
    for b in BOOSTS:
        out = {}
        for name in ("fast", "exact"):
            rng = np.random.default_rng(20260810)
            t0 = time.time()
            if name == "fast":
                h = core.binary_macro_stepping_per_shot(
                    N_REAL, F, 1, core.R_SIG_UNIT_GEN, core.TF_GEN, 0.0,
                    core.CENTERS, rng, boost=b)[:, 0, :].astype(float)
            else:
                h = core.macro_hist_exact(
                    N_REAL, F, core.R_SIG_UNIT_GEN, 0.0, core.TF_GEN,
                    core.CENTERS, rng, boost=b)
            w = h.mean(axis=0)
            j = int(w.argmax())
            out[name] = dict(
                mx=w.max(),
                plat=w[j:j + 10].mean(),
                nfull=int((w > 0.99 * core.N_PIX_MACRO).sum()),
                pk=h[:, core.I_STAT0:core.I_STAT1].max(axis=1).mean(),
                sec=time.time() - t0)
        f_, e_ = out["fast"], out["exact"]
        rel = abs(f_["plat"] - e_["plat"]) / max(e_["plat"], 1e-9)
        print(f"  {b:>8.4g} │ {f_['mx']:>12.2f} {f_['plat']:>12.2f} {f_['nfull']:>13} "
              f"{f_['pk']:>9.2f} │ {e_['mx']:>6.2f} {e_['plat']:>6.2f} "
              f"{e_['nfull']:>3} {e_['pk']:>6.2f} │ {rel:>9.1%}")

print("\n" + "=" * 104)
print("判据：饱和档（boost >= 0.5）两套引擎的平台值应当都逼近 n_tr=27，相对差在百分之几以内。")
print("修复前快速引擎在 boost=10 只能给到 19.11（精确 27.00），相对差 29%。")
