# -*- coding: utf-8 -*-
"""核对 PoD_esti_v11 的 SPAD 引擎与 lidar_histogram_sim_v45.ipynb 是否一致。

要回答的问题：v45 里 SPAD 不是简单的 8 ns 硬死时间，而是
  · 雪崩后过电压 Vov 归零，按 RC 指数恢复  Vov(Δ)/Vov_max = 1 − e^{−Δ/τ_RC}
  · 恢复期**仍有部分响应能力**，触发概率 = PDE_max · g(Vov)，g 为 exp 型凹函数
  · 每次雪崩把 Vov 重新打回 0，重新开始恢复（窗口顺延堆积）
这套模型在 PoD_esti_v11 里有没有被简化？

三步核对：
  A. 源码级：把 v45 cell 32 的 spad_binary_trace 抽出来，与 PoD 版逐行 diff。
  B. 比特级：同一 rng 种子下跑两份实现，输出必须逐 bin 完全相同。
  C. 统计级：PoD 实际跑扫描用的是【快速引擎】（更新过程 + H⁻¹ 直查表），
     它是精确逐光子引擎的连续时间极限。比对三者的 p_bin（每 bin 点亮概率）：
       ① p_bin_equilibrium 解析值
       ② 精确逐光子引擎 MC
       ③ 快速引擎 MC
     另外附带说明 v45 里【另一套】计数引擎 simulate_spad_shot_rc 的差别。

用法：
    $env:PYTHONIOENCODING="utf-8"
    python check_engine_vs_v45.py
"""
from __future__ import annotations

import difflib
import functools
import json
import os
import re
import textwrap
import time

print = functools.partial(print, flush=True)  # noqa: A001

os.environ.setdefault("POD_CORE_QUIET", "1")

import numpy as np

import pod_esti_v11_core as core

NB_V45 = "lidar_histogram_sim_v45.ipynb"


# ------------------------------------------------------------------ 取 v45 源码
def extract_v45_binary_trace():
    nb = json.load(open(NB_V45, encoding="utf-8"))
    src = None
    for c in nb["cells"]:
        s = "".join(c.get("source", []))
        if c["cell_type"] == "code" and "def spad_binary_trace" in s:
            src = s
            break
    if src is None:
        raise RuntimeError("v45 里没找到 spad_binary_trace")
    i = src.index("def spad_binary_trace")
    j = src.find("\ndef over_waveform", i)
    return src[i: j if j > 0 else len(src)]


def normalize(fn_src):
    """去掉注释、docstring、空行，只留可执行语句，便于逐行 diff。"""
    lines = []
    in_doc = False
    for ln in fn_src.splitlines():
        s = ln.rstrip()
        if not s.strip():
            continue
        q = s.strip().count('"""')
        if in_doc:
            if q:
                in_doc = False
            continue
        if s.strip().startswith('"""'):
            if q == 1:
                in_doc = True
            continue
        s = re.sub(r"\s+#.*$", "", s)
        if not s.strip() or s.strip().startswith("#"):
            continue
        lines.append(s.rstrip())
    return lines


def main():
    v45_src = extract_v45_binary_trace()

    # ---------------------------------------------------------- A. 源码级 diff
    print("=" * 96)
    print("A. 源码级对比：v45 cell 32 的 spad_binary_trace（默认路径） vs PoD_esti_v11")
    print("=" * 96)

    # v45 的函数体里，默认路径包在 `if not return_attrib:` 下面；把归因分支切掉
    cut = v45_src.find("    # ===== v31/v32 归因路径")
    v45_default = v45_src[:cut] if cut > 0 else v45_src
    v45_lines = [l for l in normalize(v45_default)
                 if "return_attrib" not in l]
    v45_lines = [l[4:] if l.startswith("        ") else l for l in v45_lines]

    import inspect
    pod_lines = normalize(inspect.getsource(core.spad_binary_trace))

    # 只比较「物理判定」那几行，缩进差异忽略
    key = lambda ls: [re.sub(r"\s+", " ", l).strip() for l in ls]
    d = list(difflib.unified_diff(key(v45_lines), key(pod_lines),
                                  "v45", "PoD_v11", lineterm="", n=1))
    core_stmts = ["vov_frac = 1.0 - np.exp(-d) if d < 700 else 1.0",
                  "d = (t - last) * inv_tau",
                  "n_ph = rng.poisson(mu)",
                  "av.append(t)",
                  "last = t"]
    print("关键物理语句是否两边都在：")
    for st in core_stmts:
        a = any(st in l for l in key(v45_lines))
        b = any(st in l for l in key(pod_lines))
        print(f"  {'OK ' if (a and b) else '差异'}  {st:<52} v45={a} PoD={b}")
    only = [l for l in d if l.startswith(("+", "-")) and not l.startswith(("+++", "---"))]
    print(f"\n归一化后仅存的行级差异 {len(only)} 条：")
    for l in only:
        print("   ", l)

    # ---------------------------------------------------------- B. 比特级复现
    print("\n" + "=" * 96)
    print("B. 比特级：同一 rng 种子跑两份实现，输出应逐 bin 完全相同")
    print("=" * 96)
    ns = {"np": np, "spad_response_g": core.spad_response_g}
    exec(v45_default + "\n", ns)          # noqa: S102 — 就是要跑 v45 的原始定义
    v45_fn = ns["spad_binary_trace"]

    tf = core.TF_GEN
    centers = core.CENTERS
    r_sig = np.zeros_like(tf)
    r_det_demo = float(core.r_det_for_noise(4.0, core.N_PIX_MACRO))
    r_amb_demo = r_det_demo / core.PDE

    same = 0
    tot = 60
    n_lit_a = n_lit_b = 0
    for s in range(tot):
        a = v45_fn(r_sig, r_amb_demo, tf, centers, core.PDE, core.TAU_RC,
                   core.VTH_FRAC, core.JIT, np.random.default_rng(1234 + s),
                   core.T_OVER, 0.0, core.RESP_SHAPE, core.RESP_K)
        b = core.spad_binary_trace(r_sig, r_amb_demo, tf, centers, core.PDE, core.TAU_RC,
                                   core.VTH_FRAC, core.JIT, np.random.default_rng(1234 + s),
                                   core.T_OVER, 0.0, core.RESP_SHAPE, core.RESP_K)
        same += int(np.array_equal(a, b))
        n_lit_a += int(a.sum()); n_lit_b += int(b.sum())
    print(f"  r_det = {r_det_demo:.4e} cps（对应 27 SPAD 单发 bg=4）")
    print(f"  {tot} 条轨迹逐 bin 完全一致：{same}/{tot}"
          f"　点亮 bin 总数 v45={n_lit_a} PoD={n_lit_b}")

    # ---------------------------------------------------------- C. 快速引擎一致性
    print("\n" + "=" * 96)
    print("C. 统计级：解析 / 精确逐光子引擎 / 快速引擎（扫描实际使用）三者的 p_bin")
    print("=" * 96)
    print(f"{'noise27':>8} {'r_det[cps]':>12} {'解析 p_bin':>11} {'精确引擎':>10} "
          f"{'快速引擎':>10} {'精确/解析':>10} {'快速/解析':>10}")
    i0, i1 = core.I_STAT0, core.I_STAT1
    for noise27 in (0.5, 2.0, 6.0, 12.0):
        r_det = float(core.r_det_for_noise(noise27, core.N_PIX_MACRO))
        r_amb = r_det / core.PDE
        p_ana = core.p_bin_equilibrium(r_det)[0]

        n_exact = 400 if noise27 >= 6 else 1500
        rng = np.random.default_rng(20260809)
        acc = 0
        for _ in range(n_exact):
            tr = core.spad_binary_trace(r_sig, r_amb, tf, centers, core.PDE, core.TAU_RC,
                                        core.VTH_FRAC, core.JIT, rng, core.T_OVER, 0.0,
                                        core.RESP_SHAPE, core.RESP_K)
            acc += int(tr[i0:i1].sum())
        p_exact = acc / (n_exact * (i1 - i0))

        h = core.noise_macro_hist_fast(200_000, 1, r_det, np.random.default_rng(7))
        p_fast = float(h[:, i0:i1].mean())

        print(f"{noise27:8.2f} {r_det:12.4e} {p_ana:11.5f} {p_exact:10.5f} "
              f"{p_fast:10.5f} {p_exact/p_ana:10.4f} {p_fast/p_ana:10.4f}")

    # ---------------------------------------------------------- D. 另一套引擎
    print("\n" + "=" * 96)
    print("D. 提醒：v45 里还有【另一套】引擎 simulate_spad_shot_rc（模块 7b，timestamp 计数模型）")
    print("=" * 96)
    print(textwrap.dedent("""\
        两套引擎在 v45 内部就是并存的、面向不同读出方式：
          · 模块 7b  simulate_spad_shot_rc（reset_mode='count'）
              雪崩后先判 vov_frac >= Vth_frac(0.60) 才【计数并复位】；
              低于阈值的「亚阈雪崩」既不计数、也【不复位】RC。
              输出 = 光子时间戳列表 → 直方图（v21 谱系，多 bit 读出）。
          · 模块 9b  spad_binary_trace（v30 起，PoD_esti 用的就是这个）
              不做 Vth 判定；每次雪崩都把 Vov 打回 0、并把输出拉高 t_over=8 ns，
              窗口顺延堆积。输出 = 每 1 ns 采样点的 0/1（1 bit 读出）。
              Vth_frac 只以导出量 T_OVER = −τ_RC·ln(1−Vth_frac) = 8.00 ns 的形式进入。
        这两者自洽：1 bit 前端里「输出为 1」与「Vov 尚未恢复到 Vth」是同一段 8 ns，
        期间再来雪崩会再次把结电容放电 → 窗口顺延。PoD_esti 研究的正是 1 bit 读出，
        因此用模块 9b 是对的；它与模块 7b 的差别不是「简化」，而是读出方式不同。
        """))


if __name__ == "__main__":
    t0 = time.time()
    main()
    print(f"\n用时 {time.time()-t0:.1f} s")
