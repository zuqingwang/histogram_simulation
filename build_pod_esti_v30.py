# -*- coding: utf-8 -*-
"""从 PoD_esti_v20.ipynb 组装出 PoD_esti_v30.ipynb。

设计原则
--------
1. **物理内核逐字沿用**：模块 0–4 与引擎、扫描函数都从 v20 原样抽出，只打少量
   有明确理由的补丁（见 PATCH_* 函数），保证 v30 与 v20 物理完全一致。
2. **分析层重写**：所有分析模块改成三段式
   [计算/载入缓存] → [绘图参数] → [绘图]，绘图 cell 只读缓存。
3. **markdown 瘦身**：只留基本描述 + 缩写；理论与公式推导移到 theory_PoD_esti_v30.md。

用法：
    python build_pod_esti_v30.py
"""
from __future__ import annotations

import ast
import json
import sys
import textwrap

import v30_cells as C

SRC = "PoD_esti_v20.ipynb"
DST = "PoD_esti_v30.ipynb"

# v20 里要逐字搬过来的 cell 序号
CELL_PARAMS      = 2    # 模块 0 全局参数
CELL_OPTICS      = 4    # 模块 1 光链路
CELL_WINDOW      = 6    # 模块 2 时间窗 / 宏像元
CELL_ENGINE_A    = 8    # 模块 3 器件响应
CELL_ENGINE_B    = 9    # 模块 3 三个引擎 + stats_from_hist_i
CELL_ENGCHECK    = 11   # 3c 一致性验证
CELL_INVCHECK    = 13   # 3d noise → E_lambda 反解校验
CELL_WAVEFORM    = 15   # 模块 4 纯噪声波形
CELL_NOISESCAN   = 17   # 模块 5 噪声扫描
CELL_THRESH      = 22   # 阈值函数 + THRESH
CELL_PODSCAN     = 25   # 模块 6 PoD 扫描
CELL_TRASH_THR   = 23   # → 回收站：v20 的六条 FAR 阈值大图
CELL_TRASH_M12   = 41   # → 回收站：v20 模块 12（连续阈值）


# ---------------------------------------------------------------- 工具
def _lines(text: str) -> list[str]:
    """把源码切成 nbformat 要求的行列表（除最后一行外都带 \\n）。"""
    text = text.strip("\n") + "\n"
    out = text.splitlines(keepends=True)
    if out and out[-1].endswith("\n"):
        out[-1] = out[-1][:-1]
    return out


def md(text: str) -> dict:
    return {"cell_type": "markdown", "metadata": {}, "source": _lines(text)}


def code(text: str) -> dict:
    return {"cell_type": "code", "execution_count": None, "metadata": {},
            "outputs": [], "source": _lines(text)}


def must_replace(src: str, old: str, new: str, what: str, count: int = 1) -> str:
    """带断言的替换：v20 结构一变就立刻报错，而不是悄悄产出错误的 v30。"""
    if src.count(old) != count:
        raise SystemExit(f"[补丁失败] {what}：预期命中 {count} 次，实际 "
                         f"{src.count(old)} 次\n---- 期望片段 ----\n{old}")
    return src.replace(old, new, count)


def rename_v20_to_v30(src: str) -> str:
    for a, b in [("pod_esti_v20_cache_", "pod_esti_v30_cache_"),
                 ("run_pod_v20_", "run_pod_v30_"),
                 ("pod_v20_", "pod_v30_"),
                 ("pod_esti_v20_core", "pod_esti_v30_core"),
                 ("PoD_esti_v20.ipynb", "PoD_esti_v30.ipynb")]:
        src = src.replace(a, b)
    return src


# ---------------------------------------------------------------- 各 cell 的补丁
def patch_params(s: str) -> str:
    """模块 0：新增 10% FAR、PoD 只解 4 条、新缓存名、绘图配色、引擎校验开关。"""
    # ① FAR_SPECS 增加 10%
    s = must_replace(
        s,
        '    (0.05,   "5pct",   "5%"),\n]',
        '    (0.05,   "5pct",   "5%"),\n'
        '    (0.10,   "10pct",  "10%"),      # ★ v30 新增\n]',
        "FAR_SPECS 增加 10%",
    )
    # ② PoD 只对 4 条 FAR 求临界能量 + 信号扫描 MC 数
    s = must_replace(
        s,
        "FAR_TAG_TO_LABEL = {t: lab for _, t, lab in FAR_SPECS}\n",
        "FAR_TAG_TO_LABEL = {t: lab for _, t, lab in FAR_SPECS}\n"
        "# ★ v30：阈值七条都算（由 peak 分布直接得到，几乎不花钱），\n"
        "#        但 PoD 临界能量只对下面这四条求解，机时省一大半。\n"
        "POD_FARS = [0.005, 0.01, 0.05, 0.10]\n"
        "POD_FAR_TAGS = [FAR_TAG[f] for f in POD_FARS]\n",
        "新增 POD_FARS",
    )
    # ③ 信号扫描 MC 数（模块 8 / 14 共用）
    s = must_replace(
        s,
        "N_MC_NOISE  = 1_000_000",
        "N_MC_SIG    = 20_000            # ★ v30：模块 8/14 固定信号扫描（v20 为 8000）\n"
        "N_MC_NOISE  = 1_000_000",
        "新增 N_MC_SIG",
    )
    # ③b PoD 求根参数加密：v20 的粗网格 11 点铺满 4 个数量级（0.4 decade 间距）+
    #     每点仅 300 MC，粗交点本身就能偏半个数量级，而局部窗只有 ±0.22 decade，
    #     于是常常整段错过真根。加密粗扫、放宽局部窗，并给迭代验证留轮数。
    s = must_replace(
        s,
        "N_POD_COARSE = 11        # 全局粗扫描能量点数，只负责包住过渡区\n"
        "N_MC_POD_COARSE = 300    # 每个粗扫描点的 MC 次数\n"
        "N_POD_LOCAL_PER_ROOT = 5 # 每个粗交点附近的局部能量点数\n",
        "N_POD_COARSE = 15        # ★ v30：11→15，粗网格间距从 0.4 缩到 ~0.29 decade\n"
        "N_MC_POD_COARSE = 600    # ★ v30：300→600，粗交点定位误差减半\n"
        "N_POD_LOCAL_PER_ROOT = 7 # ★ v30：5→7\n",
        "加密 PoD 粗扫描",
    )
    s = must_replace(
        s,
        "POD_LOCAL_HALF_DECADE = 0.22\n"
        "POD_VERIFY_TOL = 0.02\n",
        "POD_LOCAL_HALF_DECADE = 0.35   # ★ v30：0.22→0.35，粗交点偏一点也还能罩住真根\n"
        "POD_VERIFY_TOL = 0.02\n"
        "POD_VERIFY_ROUNDS = 6          # ★ v30：临界点迭代验证的最大轮数（原实现只有 1 步）\n",
        "放宽局部窗 + 新增 POD_VERIFY_ROUNDS",
    )
    # ④ 缓存全部换新名，且不再回落到 v11/v20（本版要求全量重算）
    s = must_replace(
        s,
        '# ★ v20：物理内核与网格与 v11 逐字相同 → v11 缓存可直接复用（读到后同步写回 v20 主名）\n'
        'CACHE_NOISE_FALLBACK = ["pod_esti_v11_cache_noise.npz"]\n'
        'CACHE_POD_FALLBACK   = ["pod_esti_v11_cache_pod.npz"]\n',
        '# ★ v30：FAR 列表与 res 结构都变了（新增 10% 与 hist_std），旧缓存一律不复用。\n'
        'CACHE_NOISE_FALLBACK = []\n'
        'CACHE_POD_FALLBACK   = []\n',
        "清空缓存 fallback",
    )
    s = must_replace(
        s,
        'print(f"  ★ v11 缓存已登记为 fallback：读到即同步写回 v20 主名，不会重算")\n',
        'print(f"  ★ v30 缓存与 v11/v20 不互通（FAR 列表与 res 字段都变了），本版全量重算")\n',
        "修正过时的 fallback 提示",
    )
    # ⑤ 全局绘图配色 + 引擎一致性验证开关
    s += textwrap.dedent('''

        # ---- ★ v30：全局绘图约定与开关 ----
        # 三个 N_shots 在全项目所有图里用同一套颜色，方便跨图对照
        _COLORS_N = {1: "tab:blue", 2: "tab:green", 4: "tab:red"}
        # 模块 3c 的精确引擎 vs 快速引擎比对很费时，且已在 v20 验证为 bit 级一致。
        # 需要重新验证时改成 True。
        RUN_ENGINE_CHECK = False
        # 模块 3d 的 noise → E_lambda 反解校验单次约 2 分钟，闭合误差已验证 <0.2%。
        RUN_INVERSE_CHECK = False

        print(f"  ★ v30：FAR 共 {len(FAR_SPECS)} 条；PoD 临界能量只解 "
              f"{[FAR_LABEL[f] for f in POD_FARS]}")
        print(f"  ★ v30：信号扫描 {N_MC_SIG:,} MC/档；RUN_ENGINE_CHECK={RUN_ENGINE_CHECK}")
        ''')
    return s


def patch_engine_b(s: str) -> str:
    """cell 9：stats_from_hist_i 增加「单条 hist 内 152 bin 的 std」充分统计。

    这样模块 10 就能直接复用模块 5 的 1e6 MC，不必再单跑一次 10 万条的小扫描。
    """
    return must_replace(
        s,
        "            peak_cnt=np.bincount(pk, minlength=n_tr + 2).astype(np.int64),\n"
        "        )",
        "            peak_cnt=np.bincount(pk, minlength=n_tr + 2).astype(np.int64),\n"
        "            # ★ v30：单条 hist_add 在统计窗内 152 个 bin 上的样本 std，\n"
        "            #        累加后除以条数就是模块 10 的「hist 内 std 均值」\n"
        "            hist_std_sum=float(a.std(axis=1).sum()),\n"
        "        )",
        "stats_from_hist_i 增加 hist_std_sum",
    )


def patch_engine_check(s: str) -> str:
    """cell 11：整段挂到 RUN_ENGINE_CHECK 开关下，默认不跑。"""
    return (
        "# ★ v30：一致性验证默认跳过（已在 v20 验证通过）。开关在模块 0。\n"
        "if not RUN_ENGINE_CHECK:\n"
        '    print("模块 3c 引擎一致性验证：已在 v20 验证为 bit 级一致，本次跳过。")\n'
        '    print("  要重跑，把模块 0 的 RUN_ENGINE_CHECK 改成 True。")\n'
        "else:\n"
        + textwrap.indent(s.strip("\n"), "    ")
    )


def patch_inverse_check(s: str) -> str:
    """cell 13：反解校验单次约 2 分钟，同样挂开关，默认不跑。"""
    return (
        "# ★ v30：反解精度校验默认跳过（已在 v20 验证闭合误差 <0.2%）。开关在模块 0。\n"
        "if not RUN_INVERSE_CHECK:\n"
        '    print("模块 3d noise → E_lambda 反解校验：已验证闭合误差 <0.2%，本次跳过。")\n'
        '    print("  要重跑，把模块 0 的 RUN_INVERSE_CHECK 改成 True。")\n'
        "else:\n"
        + textwrap.indent(s.strip("\n"), "    ")
    )


def patch_noise_scan(s: str) -> str:
    """cell 17：串行路径同步支持 hist_std，并让缓存校验认得新字段。"""
    # ① 分块统计增加 hist std
    s = must_replace(
        s,
        "    nz = a.mean(axis=1)\n"
        "    return (float(nz.sum()), float((nz*nz).sum()),\n"
        "            np.bincount(a.max(axis=1), minlength=n_tr + 2))",
        "    nz = a.mean(axis=1)\n"
        "    return (float(nz.sum()), float((nz*nz).sum()),\n"
        "            np.bincount(a.max(axis=1), minlength=n_tr + 2),\n"
        "            float(a.std(axis=1).sum()))   # ★ v30：hist 内 std",
        "_noise_chunk_stats 增加 hist std",
    )
    # ② res 结构新增 hist_std
    s = must_replace(
        s,
        '               "p_eq": np.zeros(ng), "noise_mc": np.zeros(ng), "noise_std": np.zeros(ng),',
        '               "p_eq": np.zeros(ng), "noise_mc": np.zeros(ng), "noise_std": np.zeros(ng),\n'
        '               "hist_std": np.zeros(ng),   # ★ v30',
        "run_noise_scan res 增加 hist_std",
    )
    # ③ 累加与写入
    s = must_replace(s, "        s1 = s2 = 0.0", "        s1 = s2 = s3 = 0.0",
                     "累加器增加 s3")
    s = must_replace(
        s,
        "        for p1, p2, pcnt in parts:\n"
        "            s1 += p1; s2 += p2\n",
        "        for p1, p2, pcnt, p3 in parts:\n"
        "            s1 += p1; s2 += p2; s3 += p3\n",
        "解包增加 p3",
    )
    s = must_replace(
        s,
        '        res["noise_std"][k] = np.sqrt(max(s2/n_mc - (s1/n_mc)**2, 0.0))',
        '        res["noise_std"][k] = np.sqrt(max(s2/n_mc - (s1/n_mc)**2, 0.0))\n'
        '        res["hist_std"][k] = s3 / n_mc',
        "写入 hist_std",
    )
    # ④ 缓存校验：缺 hist_std 的旧缓存一律判为不可用
    s = must_replace(
        s,
        "            and np.allclose(z[\"grid_key\"], grid_key)):\n"
        "        return z[\"res\"].item()",
        "            and np.allclose(z[\"grid_key\"], grid_key)):\n"
        "        _r = z[\"res\"].item()\n"
        "        # ★ v30：没有 hist_std 字段的一律当作旧缓存丢弃\n"
        "        if any(\"hist_std\" not in _r.get(_n, {}) for _n in N_SHOTS_LIST):\n"
        "            return None\n"
        "        return _r",
        "缓存校验要求 hist_std",
    )
    return s


_SOLVER_NEW = '''
def _probit_fit_local(boosts, pod, n_real, level, half_decade=0.6):
    """只用经验交点附近的点做 probit 拟合。

    ★ v30 修复：全域拟合会被 4 个数量级上的饱和点（PoD≈0 与 PoD≈1）拽偏，
    5% FAR 档的初值经常偏半个数量级以上。
    """
    x0 = _crossing_logboost(boosts, pod, level)
    if not np.isfinite(x0):
        return _probit_fit(boosts, pod, n_real)
    x = np.log10(np.asarray(boosts, float))
    sel = np.abs(x - x0) <= half_decade
    if sel.sum() >= 3:
        return _probit_fit(np.asarray(boosts, float)[sel],
                           np.asarray(pod, float)[sel], n_real)
    return _probit_fit(boosts, pod, n_real)


def _next_root_guess(hist, level, slope, max_step=0.5):
    """由已验证的 (log10 boost, PoD) 历史给出下一个试探点。

    有括号（一点低于目标、一点高于目标）就在 probit 空间做割线，
    割线跑出括号则退回二分；没有括号就用拟合斜率做 Newton 步并主动向外扩。
    """
    n = float(N_MC_POD_VERIFY)
    _clip = lambda p: min(max(p, 0.5 / n), 1.0 - 0.5 / n)
    below = [h for h in hist if h[1] < level]
    above = [h for h in hist if h[1] > level]
    if below and above:
        lo = max(below, key=lambda h: h[0])
        hi = min(above, key=lambda h: h[0])
        if hi[0] > lo[0]:
            zl, zh, zt = _norm.ppf(_clip(lo[1])), _norm.ppf(_clip(hi[1])), _norm.ppf(level)
            x = (lo[0] + (zt - zl) / (zh - zl) * (hi[0] - lo[0])
                 if zh > zl else 0.5 * (lo[0] + hi[0]))
            if not (lo[0] < x < hi[0]):
                x = 0.5 * (lo[0] + hi[0])
            return float(x)
    x0, p0 = hist[-1][0], hist[-1][1]
    s = slope if (slope and slope > 0) else 2.0
    dx = float(np.clip((_norm.ppf(level) - _norm.ppf(_clip(p0))) / s, -max_step, max_step))
    if dx == 0.0:
        dx = max_step if p0 < level else -max_step
    return float(x0 + dx)


def _verify_critical_batch(cands, n_shots, r_amb, seed0):
    """多轮批量迭代求根，把每个 (FAR, PoD 等级) 临界点解到验证 PoD 落进容差。

    ★ v30 修复：v20 只做一次 Newton 步、步长夹在 ±0.25 decade，初值偏 0.5 decade
    以上时根本追不回来，却仍然无条件接受结果。表现是模块 7 的临界能量曲线出现
    3–5 倍的毛刺，验证 PoD 实测 0.68 或 1.000 而不是 0.90。
    现在每轮把所有活跃候选一起并行评估（保持吞吐），最多 POD_VERIFY_ROUNDS 轮，
    最终取历史上最接近目标的那个点，并把 pod_err 一并存进记录备查。
    """
    if not cands:
        return {}
    state = [{"c": c, "i": i, "x": float(np.log10(c["boost"])), "hist": [], "done": False}
             for i, c in enumerate(cands)]

    for rnd in range(POD_VERIFY_ROUNDS):
        act = [s for s in state if not s["done"]]
        if not act:
            break
        pks = _eval_mc_jobs(
            [(float(10.0 ** s["x"]), N_MC_POD_VERIFY,
              seed0 + 7919 * s["i"] + 1_000_003 * rnd) for s in act],
            n_shots, r_amb,
        )
        for s, pk in zip(act, pks):
            c = s["c"]
            rec = _pk_to_record(10.0 ** s["x"], pk, c["T"], n_shots)
            s["hist"].append((s["x"], rec["pod"], rec))
            if abs(rec["pod"] - c["level"]) <= POD_VERIFY_TOL:
                s["done"] = True
                continue
            nx = _next_root_guess(s["hist"], c["level"], c["slope"])
            if not np.isfinite(nx):
                s["done"] = True
            else:
                s["x"] = nx

    finals = {}
    for s in state:
        c = s["c"]
        best = min(s["hist"], key=lambda h: abs(h[1] - c["level"]))
        rec = dict(best[2])
        rec["verify_rounds"] = len(s["hist"])
        rec["pod_err"] = float(best[1] - c["level"])
        finals[(c["tag"], f"{c['level']:.2f}")] = rec
    return finals
'''


def patch_pod_scan(s: str) -> str:
    """cell 25：① 临界能量只对 POD_FARS 求解；② 换成带括号的迭代求根（修毛刺）。"""
    s = must_replace(
        s,
        "    T_map = {FAR_TAG[far]: int(Tr[\"T\" + FAR_TAG[far]][k]) for far in TARGET_FARS}",
        "    # ★ v30：只对 POD_FARS 求 PoD 交点（阈值本身七条都在 THRESH 里）\n"
        "    T_map = {FAR_TAG[far]: int(Tr[\"T\" + FAR_TAG[far]][k]) for far in POD_FARS}",
        "T_map 限制到 POD_FARS",
    )
    # ② 整体替换旧的一次性 Newton 验证器
    i0 = s.index("def _verify_critical_batch(")
    i1 = s.index("def solve_pod_noise(")
    s = s[:i0] + _SOLVER_NEW.strip("\n") + "\n\n\n" + s[i1:]
    # ③ 初值改用「交点附近的局部 probit 拟合」，并在拟合明显跑偏时退回经验交点
    s = must_replace(
        s,
        "        slope, intercept = _probit_fit(boosts_fit, pod_fit, N_MC_POD_LOCAL)\n"
        "        for level in POD_LEVELS:\n"
        "            x_root = (_norm.ppf(level) - intercept) / slope if slope > 0 else np.nan\n",
        "        for level in POD_LEVELS:\n"
        "            # ★ v30：逐 level 做局部 probit 拟合，别让远处的饱和点拽偏初值\n"
        "            slope, intercept = _probit_fit_local(\n"
        "                boosts_fit, pod_fit, N_MC_POD_LOCAL, level)\n"
        "            x_root = (_norm.ppf(level) - intercept) / slope if slope > 0 else np.nan\n"
        "            x_emp = _crossing_logboost(boosts_fit, pod_fit, level)\n"
        "            if np.isfinite(x_emp) and (not np.isfinite(x_root)\n"
        "                                       or abs(x_root - x_emp) > 0.5):\n"
        "                x_root = x_emp   # 拟合外推得离谱时，经验交点更可信\n",
        "初值改为局部 probit 拟合",
    )
    return s


# ---------------------------------------------------------------- 新写的 markdown
MD_OVERVIEW = """
# PoD_esti v30 —— 宏像元 SPAD 直方图的检测阈值与探测概率

**这个文件在做什么**：在 1 bit（二值）SPAD 宏像元 + 多发累加的条件下，
把「底噪有多高」「阈值该定多少」「要多强的回波才能以给定概率被检出」「因此能测多远」
这一整条链算通，并给出可复现的蒙特卡洛证据。

**贯穿全文的横轴只有一个：bg**。所有对比都在同一个 bg 上做。

## 文件结构

| 模块 | 内容 | 数据来源 |
|---|---|---|
| 0–4 | 参数、光链路、时间窗、SPAD 引擎、波形示例 | 现算（很快） |
| 5 | 纯噪声 bg 扫描 → peak 分布 + 检测阈值 | **本文件唯一的大重算** |
| 6–7 | PoD 临界能量与全 bg 汇总 | PoD 扫描 |
| 8 | 固定信号 × 全 bg 网格 | 信号扫描 |
| 9–11 | 阈值倍数 ρ、三个统计量、有效 z 值 | **复用模块 5** |
| 12–13 | 所需信号、平方反比测远 | **复用模块 6** |
| 14 | 同信号不同 bg 的 peak 分布 | **复用模块 8** |
| 15 | 宏像元 3×9 vs 3×6 | 独立脚本缓存 |
| 回收站 | 已废弃代码，注释保留 | — |

每个分析模块都是三段式：**计算/载入缓存 → 绘图参数 → 绘图**。
只想调图不想重算的话，改「绘图参数」那个 cell 再跑「绘图」cell 就行。

## 缩写

| 缩写 | 全称 | 含义 |
|---|---|---|
| SPAD | Single-Photon Avalanche Diode | 单光子雪崩二极管 |
| TCSPC | Time-Correlated Single Photon Counting | 时间相关单光子计数 |
| ToF | Time of Flight | 飞行时间 |
| PDE | Photon Detection Efficiency | 光子探测效率 |
| IRF | Instrument Response Function | 仪器响应函数 |
| BRDF | Bidirectional Reflectance Distribution Function | 双向反射分布函数 |
| FAR | False Alarm Rate | 虚警率，`P(peak ≥ T)` |
| PoD | Probability of Detection | 探测概率 |
| MC | Monte Carlo | 蒙特卡洛 |
| EVT | Extreme Value Theory | 极值理论 |
| RC | Resistor-Capacitor | 阻容（过偏压恢复模型）|

## 核心记号

| 记号 | 定义 |
|---|---|
| `hist_i` | 单发直方图，每 bin 取值 0…n_pix |
| `hist_add` | N 发累加直方图，每 bin 取值 0…n_tr |
| `n_tr` | 轨迹数 `= n_pix × N_shots`（3×9、N=4 → 108）|
| **`bg`** | `hist_add` 在统计窗内每 bin 的平均计数。**全项目统一横轴** |
| `noise` | 单发 `hist_i` 每 bin 的平均计数，`= bg / N_shots`。只用来描述单次直方图的底噪 |
| `peak` | 一条 `hist_add` 在统计窗 152 个 bin 内的最大计数 |
| `T` | 检测阈值，**整数计数**（见下） |
| `boost` | 回波能量倍率，`boost=0` 为纯噪声基线 |

> **为什么阈值 T 只能是整数**：`hist_add` 是 27×N 条二值轨迹求和，取值只能是整数 0…n_tr，
> 判定就是 `peak >= T`，所以 `T=10.25` 与 `T=11` 完全等价。阈值曲线看着像阶梯，
> 是因为 **bg 以 0.25 连续步进而 T 只能按 1 跳**，不是画法问题。

## 四条原始需求分别在哪个模块回答

| 需求 | 说明 | 去哪看 |
|---|---|---|
| ① 扫 N=1/2/4，bg 步长 0.25，给阈值曲线；每个 bg 给「单条 hist 的 std 均值」「peak 均值」「peak 标准差」 | 阈值曲线在 **模块 5**；三个统计量在 **模块 10**（①②③ 三个子图） | 模块 5、10 |
| ② 5% / 1% FAR 下，PoD50 与 PoD90 所需的信号均值，随 bg 与 N_shots | **模块 12**（上排净峰高 `S_net`，下排临界发射能量 `E_crit`） | 模块 12 |
| ③ 平方反比下的测距能力 | **模块 13**（①纯 `1/D²`，②再叠大气衰减 `e^{-2αD}`） | 模块 13 |
| ④ peak 分布怎么随 bg 变；peak 均值是不是「信号峰 + bg」直接相加；同信号下 peak 标准差怎么变 | **模块 14** 三排分别回答：①分布形状 ②加法假设检验（含理想「信号峰+bg」斜率 1 参考线）③有信号 vs 无信号的 σ_peak | 模块 14 |

④ 的结论先说在这里：**peak 均值不是简单相加**，实测始终低于理想线（斜率约 0.77–0.85，
N 越大压得越狠），因为 SPAD 是 1 bit——某个 bin 已被环境光点亮时，信号光子再来也不增加计数
（抢占效应）；**σ_peak 在有信号时反而略大于无信号**，信号自身带来的二项涨落盖过了
「峰位被钉住」带来的收窄。

理论推导、公式与历史踩坑不在本文件里，见 `theory_PoD_esti_v30.md`、
`worklog_PoD_esti.md`、`handoff_PoD_esti.md`。
"""

MD_M0 = """
## 模块 0 — 全局参数与常数

整个 notebook 唯一的参数入口：光链路、SPAD 器件、时间窗、宏像元、bg 网格、
FAR 目标、MC 规模、并行线程与缓存文件名。**改任何物理量都只改这个 cell。**

- **FAR_SPECS**：七条虚警率目标（10 ppm → 10%），阈值全算。
- **POD_FARS**：只对这四条（0.5% / 1% / 5% / 10%）求 PoD 临界能量。
- **BG_GRID**：0.25 → 12，步长 0.25，共 48 档；仿真时单发 `noise = bg / N_shots`。
- **RUN_ENGINE_CHECK**：模块 3c 引擎比对开关，默认 `False`。
"""

MD_M1 = """
## 模块 1 — 光链路：从发射脉冲到单个 SPAD 上的光子率

从 `lidar_histogram_sim_v45.ipynb` 逐行移植，函数与参数一致。依次给出激光脉冲波形、
目标回波、环境光背景，以及经光学系统与 **PDE**（Photon Detection Efficiency，光子探测效率）
后落到单个 SPAD 上的光子到达率。

- **BRDF**（Bidirectional Reflectance Distribution Function，双向反射分布函数）：目标散射模型。
"""

MD_M2 = """
## 模块 2 — 时间窗、宏像元与掐头去尾

定义采集窗、1 ns bin 网格、宏像元（3×9 = 27 个 SPAD），以及**统计窗**：
采集窗前后各掐掉 `TRIM_NS`，只在中间 152 个 bin 上取 peak 和算 bg，
避开暖机段与截断边缘。全项目所有统计量都在这个统计窗里算。
"""

MD_M3 = """
## 模块 3 — SPAD 二值采样引擎

SPAD 是 1 bit 器件：一个 1 ns bin 内无论来多少光子，最多只记 1 次。雪崩之后过偏压按
**RC**（Resistor-Capacitor，阻容）规律恢复，探测能力随之同步回升；期间若再次雪崩，
电压被打回 0 重新恢复。不是简单的固定死时间。

| 引擎 | 用途 |
|---|---|
| `spad_binary_trace` | 精确版（v45 原样），作基准 |
| `noise_macro_hist_fast` | 快速版 A：$H^{-1}$ 直查表，仅纯环境光 |
| `binary_macro_stepping` | 快速版 B：同步时间步进，含信号 |
| `stats_from_hist_i` | 由单发直方图派生各 N 的充分统计量 |

三个引擎为什么等价，见 `theory_engine_equivalence.md`。
"""

MD_M3C = """
### 3c 引擎一致性验证（默认跳过）

逐条比对精确引擎与快速引擎的输出。已在 v20 验证为 bit 级一致，因此默认
`RUN_ENGINE_CHECK = False` 直接跳过（这一步很费时）。要重新验证，把模块 0 的开关改成 `True`。
"""

MD_M3D = """
### 3d `noise → E_lambda` 反解精度校验

给定目标底噪，反解所需的环境光谱辐照度 `E_lambda`，再正向算回底噪，检查闭合误差。
单次约 2 分钟，已验证闭合误差 <0.2%，因此默认 `RUN_INVERSE_CHECK = False` 跳过。
要重新验证，把模块 0 的开关改成 `True`。
"""

MD_M4 = """
## 模块 4 — 纯噪声波形长什么样

画一条典型的纯环境光 `hist_add`，直观看清底噪起伏幅度与 peak 的位置，
为后面的阈值讨论先建立感觉。
"""

MD_TRASH = """
## 回收站（不运行）

以下是 v20 中已被替换或判定无用的代码，**整段注释保留只为可追溯**。
变量名与 v30 已经不一致，不要取消注释直接跑。

- v20 模块 6 的六条 FAR 阈值大图 → v30 并入模块 5，且只画 1% / 5%。
- v20 模块 12 的连续（实数）阈值 → 已确认阈值必须是整数计数，整块作废；
  其中 12C 的「有效 z 值」被保留，成为 v30 模块 11。
"""


# ---------------------------------------------------------------- 组装
def main() -> None:
    nb = json.load(open(SRC, encoding="utf-8"))
    cells_v20 = nb["cells"]

    def src(i: int) -> str:
        c = cells_v20[i]
        if c["cell_type"] != "code":
            raise SystemExit(f"[结构不符] v20 cell {i} 不是 code cell")
        return "".join(c["source"])

    def commented(i: int, title: str) -> str:
        body = textwrap.indent(rename_v20_to_v30(src(i)).strip("\n"), "# ",
                               lambda _l: True)
        return f"# ================= {title} =================\n{body}\n"

    P = lambda s: rename_v20_to_v30(s)  # noqa: E731

    out: list[dict] = [
        md(MD_OVERVIEW),

        md(MD_M0),
        code(patch_params(P(src(CELL_PARAMS)))),

        md(MD_M1),
        code(P(src(CELL_OPTICS))),

        md(MD_M2),
        code(P(src(CELL_WINDOW))),

        md(MD_M3),
        code(P(src(CELL_ENGINE_A))),
        code(patch_engine_b(P(src(CELL_ENGINE_B)))),

        md(MD_M3C),
        code(patch_engine_check(P(src(CELL_ENGCHECK)))),

        md(MD_M3D),
        code(patch_inverse_check(P(src(CELL_INVCHECK)))),

        md(MD_M4),
        code(P(src(CELL_WAVEFORM))),

        # ---- 模块 5：唯一的大重算 ----
        # 计算段拆成两个 cell：① 噪声 MC 扫描 ② 由 peak 分布反解阈值。
        # 拆开是为了让 build_pod_core_v30.py 能在「自动开跑」之前干净地截断 ①。
        md(C.M5_MD),
        code(patch_noise_scan(P(src(CELL_NOISESCAN)))),
        code(P(src(CELL_THRESH)).strip("\n") + "\n" + C.M5_HELPERS),
        code(C.M5_PARAM),
        code(C.M5_PLOT),

        # ---- 模块 6：PoD 临界能量 ----
        md(C.M6_MD),
        code(patch_pod_scan(P(src(CELL_PODSCAN)))),
        code(C.M6_PARAM),
        code(C.M6_PLOT),

        # ---- 模块 7：全 bg 汇总 ----
        md(C.M7_MD),
        code(C.M7_COMPUTE),
        code(C.M7_PARAM),
        code(C.M7_PLOT),

        # ---- 模块 8：固定信号 × bg ----
        md(C.M8_MD),
        code(C.M8_COMPUTE),
        code(C.M8_PARAM),
        code(C.M8_PLOT),

        # ---- 模块 9–11：全部复用模块 5 ----
        md(C.M9_MD),
        code(C.M9_PARAM),
        code(C.M9_PLOT),

        md(C.M10_MD),
        code(C.M10_PARAM),
        code(C.M10_PLOT),

        md(C.M11_MD),
        code(C.M11_PARAM),
        code(C.M11_PLOT),

        # ---- 模块 12–13：复用模块 6 ----
        md(C.M12_MD),
        code(C.M12_PARAM),
        code(C.M12_PLOT),

        md(C.M13_MD),
        code(C.M13_PARAM),
        code(C.M13_PLOT),

        # ---- 模块 14：复用模块 8 ----
        md(C.M14_MD),
        code(C.M14_PARAM),
        code(C.M14_PLOT),

        # ---- 模块 15：宏像元 ----
        md(C.M15_MD),
        code(C.M15_COMPUTE),
        code(C.M15_PARAM),
        code(C.M15_PLOT),

        # ---- 回收站 ----
        md(MD_TRASH),
        code(commented(CELL_TRASH_THR, "v20 模块 6：六条 FAR 阈值大图（作废）")
             + "\n"
             + commented(CELL_TRASH_M12, "v20 模块 12：连续阈值（作废）")),
    ]

    nb_out = {
        "cells": out,
        "metadata": nb.get("metadata", {}),
        "nbformat": nb.get("nbformat", 4),
        "nbformat_minor": nb.get("nbformat_minor", 5),
    }

    # ---- 自检：每个 code cell 必须能被 ast 解析 ----
    bad = []
    for i, c in enumerate(out):
        if c["cell_type"] != "code":
            continue
        text = "".join(c["source"])
        try:
            ast.parse(text)
        except SyntaxError as e:
            bad.append((i, e))
    if bad:
        for i, e in bad:
            print(f"[语法错误] cell {i}: {e}", file=sys.stderr)
        raise SystemExit(f"共 {len(bad)} 个 cell 语法不通过，未写出文件")

    with open(DST, "w", encoding="utf-8") as f:
        json.dump(nb_out, f, ensure_ascii=False, indent=1)

    n_code = sum(c["cell_type"] == "code" for c in out)
    n_md = sum(c["cell_type"] == "markdown" for c in out)
    print(f"[完成] {DST}：{len(out)} cells（code {n_code} / markdown {n_md}），"
          f"全部通过语法自检")


if __name__ == "__main__":
    main()
