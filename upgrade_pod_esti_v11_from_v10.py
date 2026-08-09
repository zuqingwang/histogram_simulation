# -*- coding: utf-8 -*-
"""从 PoD_esti_v10 升级为 v11。

核心变更：
1. 各 N 的目标 bg 网格统一为 BG_GRID（步长 0.25，0.25→12）
2. 每档仿真 noise_amb = bg_target / N（不再用 AMB×N 导致 N=2/4 步长变粗）
3. 新增模块 10：同 bg 下阈值倍数 ρ=T_N/T_1 是否近似常数，并给出物理解释
4. 新缓存名 pod_esti_v11_*；禁止复用 v10 缓存

用法：先 Copy-Item PoD_esti_v10.ipynb PoD_esti_v11.ipynb，再
  python upgrade_pod_esti_v11_from_v10.py
"""
from __future__ import annotations

import json
import re

PATH = "PoD_esti_v11.ipynb"
nb = json.load(open(PATH, encoding="utf-8"))


def src(i):
    return "".join(nb["cells"][i].get("source", []))


def set_src(i, text):
    body = text.strip("\n")
    lines = body.split("\n")
    nb["cells"][i]["source"] = [ln + "\n" for ln in lines[:-1]] + [lines[-1]]
    nb["cells"][i]["outputs"] = []
    nb["cells"][i]["execution_count"] = None


def clear_all_outputs():
    for c in nb["cells"]:
        if c["cell_type"] == "code":
            c["outputs"] = []
            c["execution_count"] = None


def md_cell(text):
    body = text.strip("\n")
    lines = body.split("\n")
    return {
        "cell_type": "markdown",
        "id": f"v11md{abs(hash(body)) % 10**10}",
        "metadata": {},
        "source": [ln + "\n" for ln in lines[:-1]] + [lines[-1]],
    }


def code_cell(text):
    body = text.strip("\n")
    lines = body.split("\n")
    return {
        "cell_type": "code",
        "id": f"v11cd{abs(hash(body)) % 10**10}",
        "metadata": {},
        "execution_count": None,
        "outputs": [],
        "source": [ln + "\n" for ln in lines[:-1]] + [lines[-1]],
    }


# =====================================================================
# 0) 标题
# =====================================================================
s0 = src(0)
s0 = s0.replace("PoD_esti v10", "PoD_esti v11")
s0 = s0.replace("**v10**", "**v11**")
if "v11 相对 v10" not in s0:
    s0 += """

---

## ★ v11 相对 v10 的增量

1. **统一 bg 步长 0.25**：`BG_GRID = 0.25→12 / 0.25`，N=1/2/4 共用同一目标 bg 网格（v10 中 N=2/4 因 `AMB×N` 步长变成 0.5/1.0）。
2. 每档仿真：`noise_amb = bg_target / N`（同 bg 公平对比不同 shot）。
3. **新增模块 10**：同 bg 下阈值倍数 ρ=T_N/T_1 是否近似常数，并解释倍数来源。
4. **不复用** v10 缓存；物理参数数值与 v10/v05 一致。
"""
set_src(0, s0)

# =====================================================================
# 1) 参数 cell
# =====================================================================
s2 = src(2)
s2 = s2.replace("★ v10：一次仿真发满 4 发", "★ v11：一次最多仿 4 发（按 N 分别扫 bg）")
s2 = s2.replace("★ v10：由 hist_i 前缀和派生", "★ v11：N∈{1,2,4}；同 bg 网格对比")

old_grid = """# ★ v10：按单次 noise（环境标准）扫；各 N 的目标 bg≈N·noise（键仍用 bg 刻度，兼容模块 5–8）
NOISE_GRID_AMB = np.round(np.arange(0.25, 12.0 + 1e-9, 0.25), 4)  # 48 档
NOISE_GRID = {n: np.round(NOISE_GRID_AMB * n, 4) for n in N_SHOTS_LIST}"""
new_grid = """# ★ v11：各 N 目标 bg 网格统一（步长 0.25）；仿真时 noise_amb = bg / N
BG_GRID = np.round(np.arange(0.25, 12.0 + 1e-9, 0.25), 4)  # 48 档，统一 bg
NOISE_GRID = {n: BG_GRID.copy() for n in N_SHOTS_LIST}       # 目标 = bg（兼容旧键 noise_target）
NOISE_GRID_AMB = BG_GRID.copy()  # 兼容旧变量名；v11 中表示统一 bg 网格，不再表示「单次 ambient 扫轴」"""
if old_grid not in s2:
    raise SystemExit("NOISE_GRID block not found in cell 2")
s2 = s2.replace(old_grid, new_grid)

s2 = s2.replace("pod_esti_v10_cache_", "pod_esti_v11_cache_")
s2 = s2.replace("★ v10 禁止复用旧缓存，全量重算", "★ v11 禁止复用 v10/旧缓存，全量重算")
s2 = s2.replace(
    'print(f"  ★ v10 单次 noise 网格：{NOISE_GRID_AMB[0]:g}→{NOISE_GRID_AMB[-1]:g}，共 {NOISE_GRID_AMB.size} 档；各 N 目标 bg≈N·noise")',
    'print(f"  ★ v11 统一 BG_GRID：{BG_GRID[0]:g}→{BG_GRID[-1]:g}，共 {BG_GRID.size} 档，步长 0.25；noise_amb=bg/N")',
)
# 打印各 N 网格时说明 amb
s2 = s2.replace(
    'print(f"  ★ N_shots={_ns} 的噪声网格：noise = {_g[0]:g} → {_g[-1]:g}，"',
    'print(f"  ★ N_shots={_ns} 目标 bg = {_g[0]:g} → {_g[-1]:g}（noise_amb=bg/{_ns}），"',
)
# 全局 v10→v11 字面（谨慎：只改注释/打印中的版本标记）
s2 = s2.replace("CACHE_SIG   = \"pod_esti_v11_cache_signal.npz\"  # 模块 9 固定信号",
                "CACHE_SIG   = \"pod_esti_v11_cache_signal.npz\"  # 模块 9 固定信号")
set_src(2, s2)

# =====================================================================
# 2) 模块 5 markdown：说明统一 bg
# =====================================================================
s16 = src(16)
s16 = s16.replace("v10", "v11")
if "统一 bg" not in s16:
    s16 += """

### ★ v11 网格口径

- 横轴对比量统一为 **bg**（`hist_add` 统计窗均值），步长 **0.25**，范围 0.25→12。
- N=1/2/4 **共用同一 `BG_GRID`**；仿真时单次底 `noise_amb = bg / N`。
- 字段 `noise_target` / `noise_mc` 仍表示 **目标/实测 bg**（兼容模块 6–9 旧作图）。
"""
set_src(16, s16)

# =====================================================================
# 3) 替换 run_noise_scan_v10_amb → v11 按 (N, bg) 扫
# =====================================================================
s17 = src(17)

NEW_SCAN = r'''
def run_noise_scan_v11_bg(bg_grid, n_mc, chunk, seed0=2000, verbose_every=1,
                          res_all=None, on_progress=None):
    """★ v11：按统一 bg 网格扫；对每个 N 单独设 noise_amb=bg/N。

    返回 {N: res_dict}；res_dict 字段兼容旧作图（noise_mc=实测 bg，noise_target=目标 bg）。
    """
    grid = np.asarray(bg_grid, float)
    ng = len(grid)
    if res_all is None:
        res_all = {}
    for n in N_SHOTS_LIST:
        n_tr = N_PIX_MACRO * n
        if n not in res_all:
            res_all[n] = {
                "n_shots": n, "n_tr": n_tr,
                "noise_target": grid.copy(),                 # 目标 bg
                "noise_amb_target": np.round(grid / n, 6),  # 对应单次 noise
                "r_det": np.zeros(ng), "e_lambda": np.zeros(ng), "p_eq": np.zeros(ng),
                "noise_mc": np.zeros(ng), "noise_std": np.zeros(ng),  # = bg
                "noise_amb_mc": np.zeros(ng), "noise_amb_std": np.zeros(ng),
                "peak_cnt": np.zeros((ng, n_tr + 2), dtype=np.int64),
                "done": np.zeros(ng, dtype=bool),
            }
        elif "done" not in res_all[n]:
            res_all[n]["done"] = np.array(
                [int(c.sum()) > 0 for c in res_all[n]["peak_cnt"]], dtype=bool)

    t_start = time.time()
    jobs = [(n, k) for n in N_SHOTS_LIST for k in range(ng)
            if not bool(res_all[n]["done"][k])]
    n_jobs = len(jobs)
    print(f"v11 噪声扫描：{ng} bg × N={list(N_SHOTS_LIST)} = {ng*len(N_SHOTS_LIST)} 档，"
          f"待算 {n_jobs}，每档 {n_mc:,} MC", flush=True)
    for ji, (n, k) in enumerate(jobs):
        bg_t = float(grid[k])
        nt_amb = bg_t / n
        r_det = float(r_det_for_noise(float(nt_amb), N_PIX_MACRO))
        e_lam = float(e_lambda_for_r_det(r_det))
        p_eq = float(p_bin_equilibrium(r_det)[0])
        inv_tab = build_inv_table(r_det)

        acc = dict(noise_sum=0.0, noise_sumsq=0.0, bg_sum=0.0, bg_sumsq=0.0,
                   peak_cnt=np.zeros(N_PIX_MACRO * n + 2, dtype=np.int64), nn=0)
        done_m, part = 0, 0
        while done_m < n_mc:
            m = min(chunk, n_mc - done_m)
            seeds = [seed0 + 10007 * (n * 1000 + k) + 104729 * part + 17 * t
                     for t in range(NOISE_WORKERS)]
            ms = [m // NOISE_WORKERS + (1 if t < m % NOISE_WORKERS else 0)
                  for t in range(NOISE_WORKERS)]

            def _one(args, _n=n, _rd=r_det, _it=inv_tab):
                mm, sd = args
                if mm <= 0:
                    return None
                rng = np.random.default_rng(sd)
                hi = noise_hists_per_shot(mm, _n, _rd, rng, inv_tab=_it)
                return stats_from_hist_i(hi, n_shots_list=[_n])

            with ThreadPoolExecutor(max_workers=NOISE_WORKERS) as pool:
                parts = list(pool.map(_one, zip(ms, seeds)))
            for st in parts:
                if st is None:
                    continue
                b = st[n]
                acc["noise_sum"] += b["noise_sum"]; acc["noise_sumsq"] += b["noise_sumsq"]
                acc["bg_sum"] += b["bg_sum"]; acc["bg_sumsq"] += b["bg_sumsq"]
                acc["peak_cnt"] += b["peak_cnt"]; acc["nn"] += b["n"]
            done_m += m; part += 1

        R = res_all[n]; nn = max(acc["nn"], 1)
        R["r_det"][k] = r_det; R["e_lambda"][k] = e_lam; R["p_eq"][k] = p_eq
        R["noise_amb_mc"][k] = acc["noise_sum"] / nn
        R["noise_amb_std"][k] = float(np.sqrt(max(
            acc["noise_sumsq"]/nn - (acc["noise_sum"]/nn)**2, 0.0)))
        R["noise_mc"][k] = acc["bg_sum"] / nn
        R["noise_std"][k] = float(np.sqrt(max(
            acc["bg_sumsq"]/nn - (acc["bg_sum"]/nn)**2, 0.0)))
        R["peak_cnt"][k] = acc["peak_cnt"]
        R["done"][k] = True
        if on_progress is not None:
            on_progress(res_all, n, k)
        if (ji % verbose_every) == 0 or ji == n_jobs - 1:
            el = time.time() - t_start
            eta = el / (ji + 1) * (n_jobs - ji - 1)
            pk = peak_stats_from_cnt(R["peak_cnt"][k])
            print(f"  [{ji+1}/{n_jobs}] N={n} bg={bg_t:.2f}（amb={nt_amb:.3f}）→ "
                  f"bg_mc={R['noise_mc'][k]:.3f} peakμ={pk['mean']:.2f}  "
                  f"已用 {el/60:.1f} min，预计剩余 {eta/60:.1f} min", flush=True)
    return res_all


# 兼容旧名
run_noise_scan_v10_amb = run_noise_scan_v11_bg
'''

# 删掉旧函数定义（从 def run_noise_scan_v10_amb 到 return res_all 后的空行）
m = re.search(
    r"\ndef run_noise_scan_v10_amb\(.*?\n    return res_all\n",
    s17, flags=re.S)
if not m:
    raise SystemExit("run_noise_scan_v10_amb not found")
s17 = s17[:m.start()] + "\n" + NEW_SCAN + s17[m.end():]

# 替换自动开跑段
old_boot = """# ---- ★ v10：按单次 noise 网格扫一次，前缀和填满 N=1/2/4（主缓存 + 检查点）----
_grid_key = np.asarray(NOISE_GRID_AMB, float)
NOISE_RES = None
_loaded_from = None
for _cand in [CACHE_NOISE, *CACHE_NOISE_FALLBACK, CACHE_NOISE_CKPT]:
    NOISE_RES = _try_load_noise_cache(_cand, _grid_key)
    if NOISE_RES is not None:
        _loaded_from = _cand
        break

def _noise_is_complete(res_all):
    for n in N_SHOTS_LIST:
        if n not in res_all:
            return False
        r = res_all[n]
        if "done" in r:
            if not np.all(r["done"]):
                return False
        else:
            if not all(int(c.sum()) > 0 for c in r["peak_cnt"]):
                return False
            r["done"] = np.ones(len(r["noise_target"]), dtype=bool)
        # v10：档数必须等于 AMB 网格
        if len(r["noise_target"]) != len(NOISE_GRID_AMB):
            return False
    return True


if NOISE_RES is not None and _noise_is_complete(NOISE_RES):
    print(f"已从缓存 {_loaded_from} 载入纯噪声 MC（每档 {N_MC_NOISE:,} 条，完整，v10 hist_i）")
    if _loaded_from != CACHE_NOISE:
        _save_noise_cache(CACHE_NOISE, NOISE_RES, _grid_key)
        print(f"已同步写入主缓存 {CACHE_NOISE}")
else:
    # ★ v10：缺缓存时自动拉起多进程脚本（Run All 可连续跑完，不报错中断）
    # 不用 notebook 内 ThreadPool：GIL 导致吃不满 CPU。
    import sys
    print("=" * 72)
    print("未找到完整 v10 噪声缓存 → 自动调用多进程扫描（ProcessPool，吃满 CPU）")
    print("下方实时打印子进程进度……")
    print("=" * 72)
    _rc = _run_cmd_stream([sys.executable, "build_pod_core_v10.py"])
    if _rc != 0:
        raise RuntimeError("build_pod_core_v10.py 失败，无法继续噪声扫描")
    _rc = _run_cmd_stream(
        [sys.executable, "run_pod_v10_noise_scan.py",
         "--workers", str(int(N_WORKERS))])
    if _rc != 0:
        raise RuntimeError("run_pod_v10_noise_scan.py 失败，请查看上方进度输出")
    NOISE_RES = _try_load_noise_cache(CACHE_NOISE, _grid_key)
    if NOISE_RES is None or not _noise_is_complete(NOISE_RES):
        raise RuntimeError(
            f"多进程扫描结束但仍无法载入完整缓存 {CACHE_NOISE}")
    print(f"多进程扫描完成，已载入 {CACHE_NOISE}")"""

new_boot = """# ---- ★ v11：统一 BG_GRID，按 (N, bg) 扫（主缓存 + 检查点）----
_grid_key = np.asarray(BG_GRID, float)
NOISE_RES = None
_loaded_from = None
for _cand in [CACHE_NOISE, *CACHE_NOISE_FALLBACK, CACHE_NOISE_CKPT]:
    NOISE_RES = _try_load_noise_cache(_cand, _grid_key)
    if NOISE_RES is not None:
        _loaded_from = _cand
        break

def _noise_is_complete(res_all):
    for n in N_SHOTS_LIST:
        if n not in res_all:
            return False
        r = res_all[n]
        if "done" in r:
            if not np.all(r["done"]):
                return False
        else:
            if not all(int(c.sum()) > 0 for c in r["peak_cnt"]):
                return False
            r["done"] = np.ones(len(r["noise_target"]), dtype=bool)
        if len(r["noise_target"]) != len(BG_GRID):
            return False
        if not np.allclose(r["noise_target"], BG_GRID, atol=1e-6):
            return False
    return True


if NOISE_RES is not None and _noise_is_complete(NOISE_RES):
    print(f"已从缓存 {_loaded_from} 载入纯噪声 MC（每档 {N_MC_NOISE:,} 条，完整，v11 统一 bg）")
    if _loaded_from != CACHE_NOISE:
        _save_noise_cache(CACHE_NOISE, NOISE_RES, _grid_key)
        print(f"已同步写入主缓存 {CACHE_NOISE}")
else:
    import sys
    print("=" * 72)
    print("未找到完整 v11 噪声缓存 → 自动调用多进程扫描（ProcessPool，吃满 CPU）")
    print("下方实时打印子进程进度……")
    print("=" * 72)
    _rc = _run_cmd_stream([sys.executable, "build_pod_core_v11.py"])
    if _rc != 0:
        raise RuntimeError("build_pod_core_v11.py 失败，无法继续噪声扫描")
    _rc = _run_cmd_stream(
        [sys.executable, "run_pod_v11_noise_scan.py",
         "--workers", str(int(N_WORKERS))])
    if _rc != 0:
        raise RuntimeError("run_pod_v11_noise_scan.py 失败，请查看上方进度输出")
    NOISE_RES = _try_load_noise_cache(CACHE_NOISE, _grid_key)
    if NOISE_RES is None or not _noise_is_complete(NOISE_RES):
        raise RuntimeError(
            f"多进程扫描结束但仍无法载入完整缓存 {CACHE_NOISE}")
    print(f"多进程扫描完成，已载入 {CACHE_NOISE}")"""

if old_boot not in s17:
    raise SystemExit("noise boot block not found")
s17 = s17.replace(old_boot, new_boot)
s17 = s17.replace("v10", "v11")  # residual labels in this cell
# undo accidental rename of function if double-applied — ensure v11 name
s17 = s17.replace("run_noise_scan_v11_bg = run_noise_scan_v11_bg",
                  "run_noise_scan_v10_amb = run_noise_scan_v11_bg")
set_src(17, s17)

# =====================================================================
# 4) 模块 7：_peaks_chunk 只仿 n_shots；脚本名 v11；打印 bg 步长
# =====================================================================
s25 = src(25)
s25 = s25.replace(
    'print(f"每种 N_shots 对自己的完整 NOISE_GRID 求解："\n'
    '      f"{[len(NOISE_GRID[n]) for n in N_SHOTS_LIST]} 档，noise 步长 0.25")',
    'print(f"每种 N_shots 对统一 BG_GRID 求解："\n'
    '      f"{[len(NOISE_GRID[n]) for n in N_SHOTS_LIST]} 档，bg 步长 0.25；noise_amb=bg/N")',
)
s25 = s25.replace(
    '''def _peaks_chunk(boost, n_shots, r_amb, n_real, seed):
    """★ v10：仿 N_SHOTS_MAX 发 hist_i，再取前 n_shots 前缀和的 peak。"""
    rng = np.random.default_rng(seed)
    hist_i = binary_macro_stepping_per_shot(
        n_real, F_VALS, N_SHOTS_MAX, R_SIG_UNIT_POD, TF_POD, r_amb, CENTERS_SIG,
        rng, boost=boost,
    )
    return hist_add_from_prefix(hist_i, n_shots).max(axis=1)''',
    '''def _peaks_chunk(boost, n_shots, r_amb, n_real, seed):
    """★ v11：按当前 N 仿 n_shots 发（r_amb 已对应 noise=bg/N）。"""
    rng = np.random.default_rng(seed)
    hist_i = binary_macro_stepping_per_shot(
        n_real, F_VALS, n_shots, R_SIG_UNIT_POD, TF_POD, r_amb, CENTERS_SIG,
        rng, boost=boost,
    )
    return hist_i.sum(axis=1).max(axis=1)''',
)
s25 = s25.replace("build_pod_core_v10.py", "build_pod_core_v11.py")
s25 = s25.replace("run_pod_v10_pod_scan.py", "run_pod_v11_pod_scan.py")
s25 = s25.replace("run_pod_v10_noise_scan.py", "run_pod_v11_noise_scan.py")
s25 = s25.replace("v10 噪声", "v11 噪声")
s25 = s25.replace("完整 v10", "完整 v11")
s25 = s25.replace("未找到完整 v10", "未找到完整 v11")
# pod grid key already uses NOISE_GRID which equals BG_GRID per N
s25 = s25.replace("0.25-noise", "0.25-bg")
set_src(25, s25)

# =====================================================================
# 5) 模块 8/9 文案与 9.3 网格
# =====================================================================
for i in [27, 29, 30, 32, 34]:
    si = src(i)
    si2 = si.replace("v10", "v11")
    si2 = si2.replace("0.25-noise", "0.25-bg")
    if i == 27:
        si2 = si2.replace(
            "- N_shots=1：0.25–12，步长 0.25；\n- N_shots=4：0.25–40，步长 0.25。",
            "- N=1/2/4：**统一** bg=0.25–12，步长 0.25（v11）。",
        )
    if si2 != si:
        set_src(i, si2)

s35 = src(35)
s35 = s35.replace("v10", "v11")
s35 = s35.replace("NOISE_GRID_AMB", "BG_GRID")
s35 = s35.replace("扫 noise（hist_i 前缀和）", "扫统一 bg（每 N：noise_amb=bg/N）")
# 9.3 内循环：按 bg 档，对每个 N 用对应 amb
# 原逻辑是一次 AMB 仿 N_SHOTS_MAX 前缀和；改为每 N 单独（与模块 5 一致）
if "binary_macro_stepping_per_shot" in s35 and "N_SHOTS_MAX" in s35:
    # 保守：保留结构但把网格改成 BG_GRID；若仍用前缀和则 amb 语义不对。
    # 重写 9.3 主循环关键关键
    pass
set_src(35, s35)

# 更彻底：替换 9.3 扫描主体中「按 AMB + 前缀和」为「按 bg × 每 N」
s35 = src(35)
# 若仍含 N_SHOTS_MAX 前缀和扫描，整段替换为 v11 版
if "N_SHOTS_MAX" in s35 and "BOOST_LIST_M9" in s35:
    # 找到开始扫描的 print 到保存缓存之间较难精确；用函数式重写整个 cell 后半
    # 简化：在 cell 开头加注释，并把仿真循环里的 amb 改为按 N
    s35_new = '''# ---- 模块 9.3：固定信号 × 统一 bg 网格（v11：每 N 用 noise_amb=bg/N）----
BOOST_LIST_M9 = np.round(np.arange(0.0, 0.032 + 1e-12, 0.004), 6)
N_MC_SIG_M9 = 8000

def _try_load_sig_m9(path):
    if not (USE_CACHE and os.path.exists(path)):
        return None
    try:
        z = np.load(path, allow_pickle=True)
        if (int(z["n_mc"]) != N_MC_SIG_M9
                or not np.allclose(z["grid_key"], BG_GRID)
                or not np.array_equal(z["n_shots_list"], np.asarray(N_SHOTS_LIST))
                or not np.allclose(z["boosts"], BOOST_LIST_M9)):
            return None
        return z
    except Exception:
        return None

_zsig = None
for _cand in [CACHE_SIG, CACHE_SIG_CKPT]:
    _zsig = _try_load_sig_m9(_cand)
    if _zsig is not None:
        print(f"模块 9.3 命中缓存 {_cand}")
        break

if _zsig is not None:
    SIG_M9 = {
        n: dict(
            peak_cnt=np.asarray(_zsig[f"peak_cnt_{n}"]),
            peak_mean=np.zeros((len(BOOST_LIST_M9), len(BG_GRID))),
            peak_std=np.zeros((len(BOOST_LIST_M9), len(BG_GRID))),
        ) for n in N_SHOTS_LIST}
    for n in N_SHOTS_LIST:
        for i in range(len(BOOST_LIST_M9)):
            for k in range(len(BG_GRID)):
                s = peak_stats_from_cnt(SIG_M9[n]["peak_cnt"][i, k])
                SIG_M9[n]["peak_mean"][i, k] = s["mean"]
                SIG_M9[n]["peak_std"][i, k] = s["std"]
else:
    print(f"开始模块 9.3：{len(BG_GRID)} bg × {len(BOOST_LIST_M9)} boost × "
          f"{N_MC_SIG_M9} MC × N={N_SHOTS_LIST}", flush=True)
    _cnt = {n: np.zeros((len(BOOST_LIST_M9), len(BG_GRID), N_PIX_MACRO * n + 2),
                        dtype=np.int64) for n in N_SHOTS_LIST}
    _done = np.zeros(len(BG_GRID), dtype=bool)
    _t0 = time.time()
    for k, bg in enumerate(BG_GRID):
        for n in N_SHOTS_LIST:
            nt_amb = float(bg) / n
            r_det = float(r_det_for_noise(nt_amb, N_PIX_MACRO))
            r_amb = r_det / PDE
            for ib, boost in enumerate(BOOST_LIST_M9):
                rng = np.random.default_rng(9000 + 1009 * k + 17 * n + 101 * ib)
                hi = binary_macro_stepping_per_shot(
                    N_MC_SIG_M9, F_VALS, n, R_SIG_UNIT_POD, TF_POD, r_amb, CENTERS_SIG,
                    rng, boost=float(boost),
                )
                pk = hi.sum(axis=1).max(axis=1).astype(int)
                _cnt[n][ib, k] = np.bincount(pk, minlength=N_PIX_MACRO * n + 2)
        _done[k] = True
        _done_n = int(_done.sum())
        _el = time.time() - _t0
        _eta = _el / _done_n * (len(BG_GRID) - _done_n)
        print(f"  sig [{_done_n}/{len(BG_GRID)}] bg={bg:.2f}  "
              f"已用 {_el/60:.1f} min，预计剩余 {_eta/60:.1f} min", flush=True)
        if (_done_n % 4 == 0) or _done_n == len(BG_GRID):
            np.savez_compressed(
                CACHE_SIG_CKPT,
                grid_key=BG_GRID, boosts=np.asarray(BOOST_LIST_M9, float),
                n_mc=N_MC_SIG_M9, n_shots_list=np.asarray(N_SHOTS_LIST),
                **{f"peak_cnt_{n}": _cnt[n] for n in N_SHOTS_LIST})
    np.savez_compressed(
        CACHE_SIG,
        grid_key=BG_GRID, boosts=np.asarray(BOOST_LIST_M9, float),
        n_mc=N_MC_SIG_M9, n_shots_list=np.asarray(N_SHOTS_LIST),
        **{f"peak_cnt_{n}": _cnt[n] for n in N_SHOTS_LIST})
    if os.path.exists(CACHE_SIG_CKPT):
        try:
            os.remove(CACHE_SIG_CKPT)
        except OSError:
            pass
    SIG_M9 = {
        n: dict(peak_cnt=_cnt[n],
                peak_mean=np.zeros((len(BOOST_LIST_M9), len(BG_GRID))),
                peak_std=np.zeros((len(BOOST_LIST_M9), len(BG_GRID))))
        for n in N_SHOTS_LIST}
    for n in N_SHOTS_LIST:
        for i in range(len(BOOST_LIST_M9)):
            for k in range(len(BG_GRID)):
                s = peak_stats_from_cnt(SIG_M9[n]["peak_cnt"][i, k])
                SIG_M9[n]["peak_mean"][i, k] = s["mean"]
                SIG_M9[n]["peak_std"][i, k] = s["std"]
    print(f"模块 9.3 完成 → {CACHE_SIG}", flush=True)

_bg_t = BG_GRID
fig, ax = plt.subplots(2, len(N_SHOTS_LIST), figsize=(5.2 * len(N_SHOTS_LIST), 8.0), sharex=True)
for j, n in enumerate(N_SHOTS_LIST):
    for i, b in enumerate(BOOST_LIST_M9):
        ax[0, j].plot(_bg_t, SIG_M9[n]["peak_mean"][i], lw=1.2, label=f"b={b:g}")
        ax[1, j].plot(_bg_t, SIG_M9[n]["peak_std"][i], lw=1.2)
    ax[0, j].set_title(f"N={n} peak均值"); ax[0, j].set_ylabel("peak mean")
    ax[1, j].set_title(f"N={n} peak std"); ax[1, j].set_xlabel("bg"); ax[1, j].set_ylabel("peak std")
    ax[0, j].legend(fontsize=7, ncol=2); ax[0, j].grid(True, alpha=0.3); ax[1, j].grid(True, alpha=0.3)
fig.suptitle("模块 9.3a　固定信号：peak 均值/std 随 bg（v11 统一步长）", fontsize=12)
fig.tight_layout(); fig.savefig("pod_v11_m9_sig_vs_bg.png", dpi=120, bbox_inches="tight"); plt.show()

# 平移检验：相对 boost=0 的均值增量
fig, ax = plt.subplots(1, len(N_SHOTS_LIST), figsize=(5.0 * len(N_SHOTS_LIST), 4.2))
if len(N_SHOTS_LIST) == 1:
    ax = [ax]
for j, n in enumerate(N_SHOTS_LIST):
    base = SIG_M9[n]["peak_mean"][0]
    for i, b in enumerate(BOOST_LIST_M9[1:], start=1):
        ax[j].plot(_bg_t, SIG_M9[n]["peak_mean"][i] - base, lw=1.3, label=f"b={b:g}")
    ax[j].set_xlabel("bg"); ax[j].set_ylabel("Δpeak_mean"); ax[j].set_title(f"N={n}")
    ax[j].legend(fontsize=7); ax[j].grid(True, alpha=0.3)
fig.suptitle("模块 9.3b　分布平移检验", fontsize=12)
fig.tight_layout(); fig.savefig("pod_v11_m9_shift.png", dpi=120, bbox_inches="tight"); plt.show()

print("模块 9.3 线性斜率摘要（Δpeak / boost，对 bg 平均）：")
for n in N_SHOTS_LIST:
    base = SIG_M9[n]["peak_mean"][0]
    for i, b in enumerate(BOOST_LIST_M9[1:], start=1):
        if b <= 0:
            continue
        slope = np.mean((SIG_M9[n]["peak_mean"][i] - base) / b)
        print(f"  N={n} boost={b:g} → 平均斜率 {slope:.2f}")
'''
    set_src(35, s35_new)

# 模块 9 图文件名
for i in [31, 33]:
    si = src(i)
    si = si.replace("pod_v10_", "pod_v11_")
    si = si.replace("v10", "v11")
    set_src(i, si)

# =====================================================================
# 6) 新增模块 10：阈值倍数对比 + 物理解释
# =====================================================================
md10 = md_cell(r'''
## 模块 10 — ★ v11：同 bg 下不同 N 的阈值倍数

### 问题
在 **统一 bg 网格（步长 0.25）** 上，比较各 N 的 FAR 阈值曲线 $T_N(\mathrm{bg})$：
比值
$$\rho_{N/1}(\mathrm{bg}) = \frac{T_N(\mathrm{bg})}{T_1(\mathrm{bg})}$$
是否近似与 bg 无关的常数？若是，这个倍数从哪里来？

### 理想对照：独立泊松叠加 → $\rho=1$
若每 bin 近似独立 Poisson，且无死时间/饱和，则
$$\mathrm{hist\_add}(N)\ \big|_{\mathrm{noise}=\mathrm{bg}/N}
\;\stackrel{d}{=}\;
\mathrm{hist\_add}(1)\ \big|_{\mathrm{noise}=\mathrm{bg}}.$$
同 bg 下 peak 同分布 ⇒ $T_N=T_1$ ⇒ **$\rho\equiv 1$**。  
因此：若仿真里 $\rho$ 明显偏离 1，只能来自 **非线性器件效应**（不是简单的「N 倍累加」）。

### 高斯极值近似：为何 $\rho$ 可能随 bg 缓慢变化
把 peak 近似为 $\mu(\mathrm{bg})+\sigma(\mathrm{bg})\cdot Z$，FAR 分位
$$T \approx \mu + z_{\mathrm{FAR}}\sigma.$$
若 Poisson 型涨落 $\mu\approx a\,\mathrm{bg}$、$\sigma\approx\sqrt{b\,\mathrm{bg}}$，则
$$T \approx a\,\mathrm{bg} + z\sqrt{b\,\mathrm{bg}}.$$
此时 $\rho_{N/1}(\mathrm{bg})$ **一般不是严格常数**（$\sqrt{\mathrm{bg}}$ 项权重随 bg 变）；  
bg 很大时 $\rho\to a_N/a_1$。可用「拟合常数 $\bar\rho$」看残差是否小。

### 本仿真偏离 $\rho=1$ 的物理来源
1. **SPAD 死时间 / 饱和**：低 rate 多发叠加 ≠ 高 rate 单发（同均值，尾部分布不同）。  
2. **宏像元有限 SPAD 数（27）与 binary stepping**：二值硬上限与相关计数。  
3. **峰提取窗内 bin 相关**：不是独立 Poisson 极值。

下面用模块 5/6 的 `THRESH`（同 bg 网格）直接算 $\rho$，并与「$\rho=1$」「拟合常数」对比。
''')

code10 = code_cell(r'''
# ---- 模块 10：同 bg 阈值倍数 ρ=T_N/T_1 ----
_FAR_M10 = [0.01, 100e-6, 10e-6]  # 1%、100ppm、10ppm
_COLORS_N = {1: "tab:blue", 2: "tab:green", 4: "tab:red"}

# 对齐到统一 BG_GRID（THRESH["noise"] 字段存的是实测 bg）
bg_ref = np.asarray(BG_GRID, float)
T_curves = {}  # (n, far) -> T(bg_ref) 线性插值到目标 bg（通常已对齐）
for n in N_SHOTS_LIST:
    Tr = THRESH[n]
    bg_m = np.asarray(Tr["noise"], float)  # 实测 bg
    for far in _FAR_M10:
        tag = FAR_TAG[far]
        T = np.asarray(Tr["T" + tag], float)
        # 按目标网格：优先用 noise_target 索引；否则插值
        tgt = np.asarray(NOISE_RES[n]["noise_target"], float)
        if len(tgt) == len(bg_ref) and np.allclose(tgt, bg_ref, atol=1e-6):
            T_curves[(n, far)] = T.copy()
        else:
            T_curves[(n, far)] = np.interp(bg_ref, bg_m, T)

# --- 图 A：T vs bg ---
fig, axes = plt.subplots(1, len(_FAR_M10), figsize=(5.2 * len(_FAR_M10), 4.6), sharey=False)
for ax, far in zip(axes, _FAR_M10):
    for n in N_SHOTS_LIST:
        ax.plot(bg_ref, T_curves[(n, far)], "-", color=_COLORS_N[n], lw=2.0, label=f"N={n}")
    ax.set_xlabel("bg（统一步长 0.25）")
    ax.set_ylabel("阈值 T")
    ax.set_title(f"T vs bg　FAR={FAR_LABEL[far]}")
    ax.grid(True, alpha=0.3); ax.legend()
fig.suptitle("模块 10A　同 bg 下 N=1/2/4 阈值曲线", fontsize=12)
fig.tight_layout()
fig.savefig("pod_v11_m10_T_vs_bg.png", dpi=120, bbox_inches="tight")
plt.show()

# --- 图 B：ρ vs bg ---
fig, axes = plt.subplots(1, len(_FAR_M10), figsize=(5.2 * len(_FAR_M10), 4.6), sharey=True)
rho_summary = {}
for ax, far in zip(axes, _FAR_M10):
    T1 = np.maximum(T_curves[(1, far)], 1e-9)
    ax.axhline(1.0, color="k", ls=":", lw=1.2, label="理想泊松 ρ=1")
    for n in N_SHOTS_LIST:
        if n == 1:
            continue
        rho = T_curves[(n, far)] / T1
        rho_summary[(n, far)] = rho
        rho_bar = float(np.nanmean(rho))
        ax.plot(bg_ref, rho, "-", color=_COLORS_N[n], lw=2.0,
                label=f"N={n}/1  (均值={rho_bar:.3f})")
        ax.axhline(rho_bar, color=_COLORS_N[n], ls="--", lw=1.0, alpha=0.7)
    ax.set_xlabel("bg")
    ax.set_ylabel(r"$\rho_{N/1}=T_N/T_1$")
    ax.set_title(f"倍数 vs bg　FAR={FAR_LABEL[far]}")
    ax.grid(True, alpha=0.3); ax.legend(fontsize=8)
fig.suptitle("模块 10B　阈值倍数随 bg 是否近似常数", fontsize=12)
fig.tight_layout()
fig.savefig("pod_v11_m10_rho_vs_bg.png", dpi=120, bbox_inches="tight")
plt.show()

# --- 图 C：相对拟合常数的残差 ---
fig, axes = plt.subplots(1, len(_FAR_M10), figsize=(5.2 * len(_FAR_M10), 4.2), sharey=True)
print("=" * 72)
print("模块 10　阈值倍数摘要（相对 N=1）")
print("=" * 72)
for ax, far in zip(axes, _FAR_M10):
    T1 = np.maximum(T_curves[(1, far)], 1e-9)
    for n in N_SHOTS_LIST:
        if n == 1:
            continue
        rho = T_curves[(n, far)] / T1
        rho_bar = float(np.nanmean(rho))
        resid = rho - rho_bar
        rel = resid / max(abs(rho_bar), 1e-9)
        ax.plot(bg_ref, rel * 100, "-", color=_COLORS_N[n], lw=1.8, label=f"N={n}")
        print(f"FAR={FAR_LABEL[far]:>7s}  N={n}/1："
              f"  ρ̄={rho_bar:.4f}  "
              f"ρ范围=[{np.nanmin(rho):.4f},{np.nanmax(rho):.4f}]  "
              f"相对残差 rms={100*np.sqrt(np.nanmean(rel**2)):.2f}%  "
              f"max|Δρ/ρ̄|={100*np.nanmax(np.abs(rel)):.2f}%")
    ax.axhline(0, color="k", ls=":", lw=1.0)
    ax.set_xlabel("bg"); ax.set_ylabel("相对残差 [%]")
    ax.set_title(f"(ρ−ρ̄)/ρ̄　FAR={FAR_LABEL[far]}")
    ax.grid(True, alpha=0.3); ax.legend()
fig.suptitle("模块 10C　相对「常数倍数」假设的残差", fontsize=12)
fig.tight_layout()
fig.savefig("pod_v11_m10_rho_resid.png", dpi=120, bbox_inches="tight")
plt.show()

# --- 额外：T≈bg+z·σ 形状检验（用 peak 的 mean/std）---
fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.4))
for n in N_SHOTS_LIST:
    R = NOISE_RES[n]
    bg = np.asarray(R["noise_mc"], float)
    st = [peak_stats_from_cnt(c) for c in R["peak_cnt"]]
    mu = np.array([s["mean"] for s in st])
    sd = np.array([s["std"] for s in st])
    ax[0].plot(bg, mu / np.maximum(bg, 1e-9), "-", color=_COLORS_N[n], lw=1.8, label=f"N={n}")
    ax[1].plot(bg, sd / np.sqrt(np.maximum(bg, 1e-9)), "-", color=_COLORS_N[n], lw=1.8, label=f"N={n}")
ax[0].set_xlabel("bg"); ax[0].set_ylabel("peak_mean / bg")
ax[0].set_title("均值系数 a≈μ/bg（泊松理想≈1）")
ax[0].grid(True, alpha=0.3); ax[0].legend()
ax[1].set_xlabel("bg"); ax[1].set_ylabel(r"peak_std / √bg")
ax[1].set_title("涨落系数 √b≈σ/√bg")
ax[1].grid(True, alpha=0.3); ax[1].legend()
fig.suptitle("模块 10D　解释倍数：同 bg 下 μ、σ 的尺度（非线性使 N 间不等）", fontsize=12)
fig.tight_layout()
fig.savefig("pod_v11_m10_mu_sigma_scale.png", dpi=120, bbox_inches="tight")
plt.show()

print("""
【解读提纲】
1. 若 ρ 全程贴近 1：同 bg 下多发叠加≈单发，器件近似线性泊松。
2. 若 ρ 近似水平但 ≠1：存在近似「常数倍数」；倍数来自饱和/死时间导致的
   同均值不同尾部 —— 低 rate×N 与高 rate×1 的 peak 分位不同。
3. 若 ρ 随 bg 明显倾斜：与高斯极值公式中 √bg 项一致，或饱和程度随 rate 变化。
4. 看 10D：若同 bg 下 μ/bg、σ/√bg 在不同 N 间分开，则 T 倍数可直接由分位差解释。
""")
''')

nb["cells"].extend([md10, code10])

# =====================================================================
# 7) 全文残留脚本名 / 图名
# =====================================================================
for i, c in enumerate(nb["cells"]):
    if c["cell_type"] != "code":
        continue
    s = src(i)
    s2 = s.replace("build_pod_core_v10.py", "build_pod_core_v11.py")
    s2 = s2.replace("run_pod_v10_noise_scan.py", "run_pod_v11_noise_scan.py")
    s2 = s2.replace("run_pod_v10_pod_scan.py", "run_pod_v11_pod_scan.py")
    s2 = s2.replace("pod_v10_", "pod_v11_")
    if s2 != s:
        set_src(i, s2)

clear_all_outputs()
json.dump(nb, open(PATH, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
print(f"OK → {PATH}，cells={len(nb['cells'])}")
