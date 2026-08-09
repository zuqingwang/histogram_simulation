# -*- coding: utf-8 -*-
"""从完整 PoD_esti_v05 升级为 v10（保留全部模块，新增 hist_i 与模块 9）。

用法：先 Copy-Item PoD_esti_v05.ipynb PoD_esti_v10.ipynb，再
  python upgrade_pod_esti_v10_from_v05.py
"""
from __future__ import annotations

import json
import re

PATH = "PoD_esti_v10.ipynb"
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


def insert_cells(at, cells):
    for j, c in enumerate(cells):
        nb["cells"].insert(at + j, c)


def md_cell(text):
    body = text.strip("\n")
    lines = body.split("\n")
    return {
        "cell_type": "markdown",
        "id": f"v10md{len(nb['cells'])}",
        "metadata": {},
        "source": [ln + "\n" for ln in lines[:-1]] + [lines[-1]],
    }


def code_cell(text):
    body = text.strip("\n")
    lines = body.split("\n")
    return {
        "cell_type": "code",
        "id": f"v10cd{len(nb['cells'])}",
        "metadata": {},
        "execution_count": None,
        "outputs": [],
        "source": [ln + "\n" for ln in lines[:-1]] + [lines[-1]],
    }


# =====================================================================
# 0) 标题
# =====================================================================
s0 = src(0)
s0 = s0.replace("PoD_esti v05", "PoD_esti v10")
s0 = s0.replace("**v05**", "**v10**")
if "hist_i" not in s0:
    s0 += """

---

## ★ v10 相对 v05 的增量（主体仍是 v05 全模块）

1. 每次实现记录 `hist_i`（最多 4 发）；N∈{1,2,4} 的 `hist_add` 取前缀和，避免重复仿真。
2. **noise** = 单次 `hist_i` 底；**bg** = `hist_add` 底；**peak** 在 `hist_add` 上统计。
3. **新增模块 9**：peak–bg 形状；1%阈值 / bg+5σ / peak均值 分三张图按 N 对比；固定信号扫 noise。
4. **不复用** v05 缓存；物理参数数值与 v05 一致。
"""
set_src(0, s0)

# =====================================================================
# 1) 参数 cell：N_SHOTS、网格、缓存
# =====================================================================
s2 = src(2)
s2 = s2.replace("N_SHOTS_LIST = [1, 4]    # 两种累加发数都做",
                "N_SHOTS_MAX  = 4                 # ★ v10：一次仿真发满 4 发\n"
                "N_SHOTS_LIST = [1, 2, 4]         # ★ v10：由 hist_i 前缀和派生")

old_grid = """NOISE_GRID = {
    1: np.round(np.arange(0.25, 12.0 + 1e-9, 0.25), 4),   # 48 档，上限 27
    4: np.round(np.arange(0.25, 40.0 + 1e-9, 0.25), 4),   # 160 档，上限 108
}"""
new_grid = """# ★ v10：按单次 noise（环境标准）扫；各 N 的目标 bg≈N·noise（键仍用 bg 刻度，兼容模块 5–8）
NOISE_GRID_AMB = np.round(np.arange(0.25, 12.0 + 1e-9, 0.25), 4)  # 48 档
NOISE_GRID = {n: np.round(NOISE_GRID_AMB * n, 4) for n in N_SHOTS_LIST}"""
if old_grid not in s2:
    raise SystemExit("NOISE_GRID block not found")
s2 = s2.replace(old_grid, new_grid)

s2 = s2.replace('CACHE_NOISE = "pod_esti_v05_cache_noise.npz"',
                'CACHE_NOISE = "pod_esti_v10_cache_noise.npz"')
s2 = s2.replace('CACHE_POD   = "pod_esti_v05_cache_pod.npz"',
                'CACHE_POD   = "pod_esti_v10_cache_pod.npz"\n'
                'CACHE_SIG   = "pod_esti_v10_cache_signal.npz"  # 模块 9 固定信号')
s2 = s2.replace(
    'CACHE_NOISE_FALLBACK = ["pod_esti_v04_cache_noise.npz"]  # 仅噪声可迁入\n'
    'CACHE_POD_FALLBACK   = []  # ★ FAR 已扩展，禁止误用旧 PoD 缓存\n'
    'CACHE_NOISE_CKPT = "pod_esti_v05_cache_noise.partial.npz"\n'
    'CACHE_POD_CKPT   = "pod_esti_v05_cache_pod.partial.npz"',
    'CACHE_NOISE_FALLBACK = []  # ★ v10 禁止复用旧缓存，全量重算\n'
    'CACHE_POD_FALLBACK   = []\n'
    'CACHE_NOISE_CKPT = "pod_esti_v10_cache_noise.partial.npz"\n'
    'CACHE_POD_CKPT   = "pod_esti_v10_cache_pod.partial.npz"\n'
    'CACHE_SIG_CKPT   = "pod_esti_v10_cache_signal.partial.npz"',
)
# 打印区补充
if "NOISE_GRID_AMB" not in s2.split("print")[-1]:
    s2 = s2.replace(
        'print(f"宏像元 {MACRO_BX}×{MACRO_BY} = {MACRO_BX*MACRO_BY} 个 SPAD；N_shots 取 {N_SHOTS_LIST}")',
        'print(f"宏像元 {MACRO_BX}×{MACRO_BY} = {MACRO_BX*MACRO_BY} 个 SPAD；'
        'N_SHOTS_MAX={N_SHOTS_MAX}；分析 N={N_SHOTS_LIST}")\n'
        'print(f"  ★ v10 单次 noise 网格：{NOISE_GRID_AMB[0]:g}→{NOISE_GRID_AMB[-1]:g}，'
        '共 {NOISE_GRID_AMB.size} 档；各 N 目标 bg≈N·noise")',
    )
set_src(2, s2)

# =====================================================================
# 2) 引擎 cell 9：追加 hist_i API
# =====================================================================
s9 = src(9)
if "noise_hists_per_shot" not in s9:
    insert = r'''

# ============================================================================
# ★ v10：per-shot hist_i 与前缀和
# ============================================================================
def noise_hists_per_shot(n_real, n_shots, r_det, rng, inv_tab=None):
    """纯噪声 hist_i：(n_real, n_shots, NBINS)，每 shot 计数 ∈[0,27]。"""
    if inv_tab is None:
        inv_tab = build_inv_table(r_det)
    out = np.zeros((n_real, n_shots, NBINS), dtype=np.int32)
    for s in range(n_shots):
        out[:, s, :] = noise_macro_hist_fast(
            n_real, N_PIX_MACRO, r_det, rng, inv_tab=inv_tab)
    return out


def binary_macro_stepping_per_shot(n_real, f_pix, n_shots, r_sig_unit, tgrid, r_amb,
                                   centers, rng, boost=1.0, tau_rc=TAU_RC, t_over=T_OVER,
                                   pde=PDE, jitter=JIT, resp_shape=RESP_SHAPE, resp_k=RESP_K):
    """信号+环境 hist_i：(n_real, n_shots, len(centers))。"""
    f_arr = np.tile(np.asarray(f_pix, float), int(n_shots))
    dt = tgrid[1] - tgrid[0]
    n_tr = f_arr.size
    n_pix = int(np.asarray(f_pix).size)
    nb = len(centers)
    k_max = int(np.ceil(20.0 * tau_rc / dt))
    phi = pde * spad_response_g(1.0 - np.exp(-np.arange(k_max + 1) * dt / tau_rc),
                                resp_shape, resp_k)
    age = np.full((n_real, n_tr), k_max, dtype=np.int32)
    tcov = np.full((n_real, n_tr), -1e30)
    hist_i = np.zeros((n_real, n_shots, nb), dtype=np.int32)
    mu_all = (r_sig_unit[:, None] * f_arr[None, :] * boost + r_amb) * dt
    ib = 0
    for i in range(tgrid.size):
        t = tgrid[i]
        while ib < nb and centers[ib] < t:
            d = centers[ib] - tcov
            lit = ((d >= 0) & (d < t_over)).reshape(n_real, n_shots, n_pix).sum(axis=2)
            hist_i[:, :, ib] = lit
            ib += 1
        p = -np.expm1(-mu_all[i][None, :] * phi[age])
        fire = rng.random((n_real, n_tr)) < p
        age = np.minimum(age + 1, k_max)
        if fire.any():
            age[fire] = 1
            nf = int(fire.sum())
            tcov[fire] = t + (rng.normal(0.0, jitter, nf) if jitter > 0 else 0.0)
    while ib < nb:
        d = centers[ib] - tcov
        lit = ((d >= 0) & (d < t_over)).reshape(n_real, n_shots, n_pix).sum(axis=2)
        hist_i[:, :, ib] = lit
        ib += 1
    return hist_i


def hist_add_from_prefix(hist_i, n_shots):
    return hist_i[:, :n_shots, :].sum(axis=1)


def stats_from_hist_i(hist_i, n_shots_list=None, i0=None, i1=None):
    """由 N_SHOTS_MAX 发 hist_i 派生各 N 的 noise/bg/peak 充分统计。"""
    if n_shots_list is None:
        n_shots_list = N_SHOTS_LIST
    if i0 is None:
        i0 = I_STAT0
    if i1 is None:
        i1 = I_STAT1
    n_real, n_max, _ = hist_i.shape
    shot_nz = hist_i[:, :, i0:i1].mean(axis=2)
    out = {}
    for n in n_shots_list:
        hadd = hist_add_from_prefix(hist_i, n)
        a = hadd[:, i0:i1]
        bg = a.mean(axis=1)
        pk = a.max(axis=1)
        nz = shot_nz[:, :n].mean(axis=1)
        n_tr = N_PIX_MACRO * n
        out[n] = dict(
            n=n_real,
            noise_sum=float(nz.sum()), noise_sumsq=float((nz*nz).sum()),
            bg_sum=float(bg.sum()), bg_sumsq=float((bg*bg).sum()),
            peak_cnt=np.bincount(pk, minlength=n_tr + 2).astype(np.int64),
        )
    return out

'''
    s9 = s9.replace(
        'print("引擎与更新过程工具就绪：")',
        insert + '\nprint("引擎与更新过程工具就绪：")',
    )
    s9 = s9.replace(
        'print("  · binary_macro_stepping —— 快速 B，同步时间步进，含信号")',
        'print("  · binary_macro_stepping —— 快速 B，同步时间步进，含信号")\n'
        'print("  · noise_hists_per_shot / binary_macro_stepping_per_shot —— ★v10 hist_i")\n'
        'print("  · hist_add_from_prefix / stats_from_hist_i —— ★v10 前缀和")',
    )
set_src(9, s9)

# =====================================================================
# 3) 模块 5 扫描：改为按 AMB 扫一次，前缀和填满 N=1/2/4
# =====================================================================
s17 = src(17)

# 替换 _noise_chunk_stats 为 hist_i 版，并重写 run_noise_scan 入口逻辑的主循环部分
# 采用：新增 run_noise_scan_v10_amb，替换开跑循环

if "run_noise_scan_v10_amb" not in s17:
    # 在 run_noise_scan 定义之后、开跑之前插入新函数，并替换开跑循环
    amb_fn = r'''

def run_noise_scan_v10_amb(noise_amb_grid, n_mc, chunk, seed0=2000, verbose_every=5,
                           res_all=None, on_progress=None):
    """★ v10：按单次 noise 扫；每档仿 N_SHOTS_MAX 发 hist_i，前缀和得到 N=1/2/4。

    返回 {N: res_dict}，res_dict 字段兼容原 run_noise_scan（noise_mc=实测 bg）。
    额外字段：noise_amb_mc（实测单次 noise）、noise_amb_target。
    """
    grid = np.asarray(noise_amb_grid, float)
    ng = len(grid)
    if res_all is None:
        res_all = {}
    for n in N_SHOTS_LIST:
        n_tr = N_PIX_MACRO * n
        if n not in res_all:
            res_all[n] = {
                "n_shots": n, "n_tr": n_tr,
                "noise_target": np.round(grid * n, 4),  # 兼容旧键：目标 bg≈N·noise
                "noise_amb_target": grid.copy(),
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
    for k, nt_amb in enumerate(grid):
        if all(bool(res_all[n]["done"][k]) for n in N_SHOTS_LIST):
            continue
        # r_det 按单次 27 SPAD 的 noise 反解
        r_det = float(r_det_for_noise(float(nt_amb), N_PIX_MACRO))
        e_lam = float(e_lambda_for_r_det(r_det))
        p_eq = float(p_bin_equilibrium(r_det)[0])
        inv_tab = build_inv_table(r_det)

        acc = {n: dict(noise_sum=0.0, noise_sumsq=0.0, bg_sum=0.0, bg_sumsq=0.0,
                       peak_cnt=np.zeros(N_PIX_MACRO * n + 2, dtype=np.int64), nn=0)
               for n in N_SHOTS_LIST}
        done_m, part = 0, 0
        while done_m < n_mc:
            m = min(chunk, n_mc - done_m)
            seeds = [seed0 + 10007 * k + 104729 * part + 17 * t
                     for t in range(NOISE_WORKERS)]
            ms = [m // NOISE_WORKERS + (1 if t < m % NOISE_WORKERS else 0)
                  for t in range(NOISE_WORKERS)]
            def _one(args):
                mm, sd = args
                if mm <= 0:
                    return None
                rng = np.random.default_rng(sd)
                hi = noise_hists_per_shot(mm, N_SHOTS_MAX, r_det, rng, inv_tab=inv_tab)
                return stats_from_hist_i(hi)
            with ThreadPoolExecutor(max_workers=NOISE_WORKERS) as pool:
                parts = list(pool.map(_one, zip(ms, seeds)))
            for st in parts:
                if st is None:
                    continue
                for n in N_SHOTS_LIST:
                    a, b = acc[n], st[n]
                    a["noise_sum"] += b["noise_sum"]; a["noise_sumsq"] += b["noise_sumsq"]
                    a["bg_sum"] += b["bg_sum"]; a["bg_sumsq"] += b["bg_sumsq"]
                    a["peak_cnt"] += b["peak_cnt"]; a["nn"] += b["n"]
            done_m += m; part += 1

        for n in N_SHOTS_LIST:
            R = res_all[n]; a = acc[n]; nn = max(a["nn"], 1)
            R["r_det"][k] = r_det; R["e_lambda"][k] = e_lam; R["p_eq"][k] = p_eq
            R["noise_amb_mc"][k] = a["noise_sum"] / nn
            R["noise_amb_std"][k] = float(np.sqrt(max(
                a["noise_sumsq"]/nn - (a["noise_sum"]/nn)**2, 0.0)))
            R["noise_mc"][k] = a["bg_sum"] / nn          # ★ 兼容旧图：noise_mc = bg
            R["noise_std"][k] = float(np.sqrt(max(
                a["bg_sumsq"]/nn - (a["bg_sum"]/nn)**2, 0.0)))
            R["peak_cnt"][k] = a["peak_cnt"]
            R["done"][k] = True
        if on_progress is not None:
            on_progress(res_all, k)
        if (k % verbose_every) == 0 or k == ng - 1:
            el = time.time() - t_start
            eta = el / (k + 1) * (ng - k - 1)
            pk4 = peak_stats_from_cnt(res_all[4]["peak_cnt"][k])
            print(f"  [amb {k+1:>3d}/{ng}] noise={nt_amb:.2f} → "
                  f"bg(N=1/2/4)="
                  f"{res_all[1]['noise_mc'][k]:.3f}/"
                  f"{res_all[2]['noise_mc'][k]:.3f}/"
                  f"{res_all[4]['noise_mc'][k]:.3f}  "
                  f"peakμ(N=4)={pk4['mean']:.2f}  [{el:.0f}s, 剩约{eta:.0f}s]")
    return res_all

'''
    # 插在 "# ---- 估算总耗时并开跑" 之前
    marker = "# ---- 估算总耗时并开跑"
    if marker not in s17:
        raise SystemExit("run marker not found in cell 17")
    s17 = s17.replace(marker, amb_fn + "\n" + marker)

    # 替换开跑主体：用 v10 amb 扫描
    # 找到 from for _cand in ... 到写主缓存结束，替换为新逻辑
    old_run_start = "_grid_key = np.concatenate([np.asarray(NOISE_GRID[n]) for n in N_SHOTS_LIST])"
    if old_run_start not in s17:
        raise SystemExit("grid_key line not found")

    # 从 _grid_key 到 cell 末尾（THRESH 不在本 cell）—— cell 17 以写完 CACHE_NOISE 结束
    idx = s17.find(old_run_start)
    head = s17[:idx]
    new_run = r'''
# ---- ★ v10：按单次 noise 网格扫一次，前缀和填满 N=1/2/4（主缓存 + 检查点）----
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
    if NOISE_RES is None:
        NOISE_RES = {}
        print("未找到匹配的 v10 噪声缓存，开始全新扫描（hist_i 前缀和）")
    else:
        print(f"从 {_loaded_from} 载入部分结果，断点续跑")
    print(f"纯噪声 MC：AMB {len(NOISE_GRID_AMB)} 档 × N={N_SHOTS_LIST}（前缀和）× {N_MC_NOISE:,} 条")
    _ckpt_counter = {"n": 0}
    def _on_progress(res_all, k):
        _ckpt_counter["n"] += 1
        if (_ckpt_counter["n"] % CHECKPOINT_EVERY) == 0:
            _save_noise_cache(CACHE_NOISE_CKPT, res_all, _grid_key)
    _tall = time.time()
    NOISE_RES = run_noise_scan_v10_amb(
        NOISE_GRID_AMB, N_MC_NOISE, MC_CHUNK, res_all=NOISE_RES, on_progress=_on_progress)
    print(f"总用时 {time.time()-_tall:.0f} s")
    _save_noise_cache(CACHE_NOISE, NOISE_RES, _grid_key)
    print(f"已写入主缓存 {CACHE_NOISE}")
    if os.path.exists(CACHE_NOISE_CKPT):
        try: os.remove(CACHE_NOISE_CKPT)
        except OSError: pass
'''
    s17 = head + new_run
set_src(17, s17)

# =====================================================================
# 4) 模块 5/5b/6 作图：适配 N=3 列；xlabel 标明 bg
# =====================================================================
for i in (18, 20, 23):
    s = src(i)
    s = s.replace("figsize=(17.5, 5.2)", "figsize=(5.5*len(N_SHOTS_LIST), 5.2)")
    s = s.replace("figsize=(18.5, 5.4)", "figsize=(5.8*len(N_SHOTS_LIST), 5.4)")
    # 1×3 固定子图 → 按 N 数量
    s = s.replace("fig, ax = plt.subplots(1, 3, figsize=(5.5*len(N_SHOTS_LIST), 5.2))",
                  "fig, ax = plt.subplots(1, len(N_SHOTS_LIST), figsize=(5.5*len(N_SHOTS_LIST), 5.2))")
    if "subplots(1, 3," in s and "N_SHOTS_LIST" in s:
        s = re.sub(r"plt\.subplots\(1, 3, figsize=\([^)]+\)\)",
                   "plt.subplots(1, len(N_SHOTS_LIST), "
                   "figsize=(5.5*len(N_SHOTS_LIST), 5.2))", s)
    s = s.replace('a.set_xlabel("noise")',
                  'a.set_xlabel("bg（hist_add 均值；≈N·noise）")')
    s = s.replace("ax.set_xlabel(\"noise\")",
                  "ax.set_xlabel(\"bg（hist_add 均值；≈N·noise）\")")
    s = s.replace('ax[0].set_xlabel("noise")',
                  'ax[0].set_xlabel("bg（≈N·noise）")')
    s = s.replace('ax[2].set_xlabel("noise")',
                  'ax[2].set_xlabel("bg（≈N·noise）")')
    s = s.replace('pod_v05_', 'pod_v10_')
    set_src(i, s)

for i in (26, 28):
    s = src(i).replace("pod_v05_", "pod_v10_")
    set_src(i, s)

# 模块 7 _peaks_chunk：改用 hist_i 前缀和
s25 = src(25)
old_pk = '''def _peaks_chunk(boost, n_shots, r_amb, n_real, seed):
    """单块 MC，无内部并行。"""
    f_arr = np.tile(F_VALS, n_shots)
    h = binary_macro_stepping(
        n_real, f_arr, R_SIG_UNIT_POD, TF_POD, r_amb, CENTERS_SIG,
        np.random.default_rng(seed), boost=boost,
    )
    return h.max(axis=1)'''
new_pk = '''def _peaks_chunk(boost, n_shots, r_amb, n_real, seed):
    """★ v10：仿 N_SHOTS_MAX 发 hist_i，再取前 n_shots 前缀和的 peak。"""
    rng = np.random.default_rng(seed)
    hist_i = binary_macro_stepping_per_shot(
        n_real, F_VALS, N_SHOTS_MAX, R_SIG_UNIT_POD, TF_POD, r_amb, CENTERS_SIG,
        rng, boost=boost,
    )
    return hist_add_from_prefix(hist_i, n_shots).max(axis=1)'''
if old_pk in s25:
    s25 = s25.replace(old_pk, new_pk)
    set_src(25, s25)
else:
    print("WARN: _peaks_chunk pattern not found; leave as-is")

# =====================================================================
# 5) 追加模块 9（三部分新图）
# =====================================================================
# 先清掉若已存在的旧模块 9
nb["cells"] = [c for c in nb["cells"]
               if "模块 9" not in "".join(c.get("source", []))[:40]]

mod9 = []
mod9.append(md_cell(r"""
## 模块 9 — ★ v10 新增：hist_i / hist_add 专项分析

本模块**追加**在 v05 全流程之后，不替代模块 5–8。

口径回顾：
- `hist_i`：第 i 发宏像元直方图；`hist_add(N)=sum(hist_1..hist_N)`（前缀和）
- **noise**：单次 `hist_i` 统计窗均值；**bg**：`hist_add` 统计窗均值；**peak** 在 `hist_add` 上

三部分：
1. 纯噪声 **peak vs bg**（N=1/2/4）形状是否一致  
2. **三张对比图**：① 1% FAR 阈值 ② bg+5·std(peak) ③ peak 均值 —— 各图画齐 N=1/2/4  
3. 固定信号、noise 线性增长：分布是否平移、均值/std 是否线性  
"""))

mod9.append(md_cell("### 9.1　纯噪声 peak–bg 曲线（N=1/2/4）"))
mod9.append(code_cell(r"""
# 直接复用模块 5 的 NOISE_RES（v10：noise_mc 字段 = 实测 bg）
_COLORS_N = {1: "tab:blue", 2: "tab:green", 4: "tab:red"}

fig, ax = plt.subplots(1, 3, figsize=(16.5, 4.8))
for n in N_SHOTS_LIST:
    R = NOISE_RES[n]
    st = [peak_stats_from_cnt(c) for c in R["peak_cnt"]]
    bg = R["noise_mc"]
    mu = np.array([s["mean"] for s in st])
    p50 = np.array([s["p50"] for s in st])
    ax[0].plot(bg, mu, "-", color=_COLORS_N[n], lw=1.8, label=f"N={n} mean")
    ax[0].plot(bg, p50, ":", color=_COLORS_N[n], lw=1.1, alpha=0.75)
ax[0].set_xlabel("bg"); ax[0].set_ylabel("peak")
ax[0].set_title("peak vs bg（实线均值，点线中位）")
ax[0].legend(fontsize=8); ax[0].grid(alpha=0.3)

for n in N_SHOTS_LIST:
    R = NOISE_RES[n]
    st = [peak_stats_from_cnt(c) for c in R["peak_cnt"]]
    bg = R["noise_mc"]
    mu = np.array([s["mean"] for s in st])
    ax[1].plot(bg / n, mu / n, "-", color=_COLORS_N[n], lw=1.8, label=f"N={n}")
ax[1].set_xlabel("bg/N ≈ noise"); ax[1].set_ylabel("peak_mean / N")
ax[1].set_title("按 N 归一（重合 ⇒ 形状一致）")
ax[1].legend(fontsize=8); ax[1].grid(alpha=0.3)

R1 = NOISE_RES[1]
st1 = [peak_stats_from_cnt(c) for c in R1["peak_cnt"]]
x1 = R1["noise_mc"] / 1.0
y1 = np.array([s["mean"] for s in st1]) / 1.0
for n in [2, 4]:
    R = NOISE_RES[n]
    st = [peak_stats_from_cnt(c) for c in R["peak_cnt"]]
    x = R["noise_mc"] / n
    y = np.array([s["mean"] for s in st]) / n
    ax[2].plot(x, y - np.interp(x, x1, y1), "-", color=_COLORS_N[n], lw=1.6,
               label=f"N={n} − N=1")
ax[2].axhline(0, color="0.4", lw=1)
ax[2].set_xlabel("bg/N ≈ noise"); ax[2].set_ylabel("Δ(peak_mean/N)")
ax[2].set_title("相对 N=1 归一曲线残差")
ax[2].legend(fontsize=8); ax[2].grid(alpha=0.3)
fig.suptitle("模块 9.1　纯噪声 peak–bg（v10 hist_i 前缀和）", fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.93])
fig.savefig("pod_v10_m9_peak_vs_bg.png", dpi=120, bbox_inches="tight")
plt.show()
"""))

mod9.append(md_cell(r"""
### 9.2　分三张图对比 N=1 / 2 / 4

横轴均为 **bg**。三张图分别画：

1. **1% FAR 阈值 T**  
2. **bg + 5·std(peak)**  
3. **peak 均值**  

便于直接看不同发数曲线差异（而不是叠在同一坐标系里挤在一起时难读）。
"""))
mod9.append(code_cell(r"""
_FAR1 = 0.01
_series = {}
for n in N_SHOTS_LIST:
    R = NOISE_RES[n]
    st = [peak_stats_from_cnt(c) for c in R["peak_cnt"]]
    bg = np.asarray(R["noise_mc"], float)
    mu = np.array([s["mean"] for s in st])
    sd = np.array([s["std"] for s in st])
    T1 = np.array([far_threshold_from_cnt(c, _FAR1)[0] for c in R["peak_cnt"]])
    _series[n] = dict(bg=bg, mu=mu, sd=sd, T1=T1, y5=bg + 5.0 * sd)

# --- 图 A：1% 阈值 ---
fig, ax = plt.subplots(figsize=(7.2, 4.8))
for n in N_SHOTS_LIST:
    s = _series[n]
    ax.plot(s["bg"], s["T1"], "-", color=_COLORS_N[n], lw=2.0, label=f"N={n}")
ax.set_xlabel("bg"); ax.set_ylabel("T @ FAR=1%")
ax.set_title("对比图 A　1% FAR 阈值 vs bg（N=1/2/4）")
ax.legend(); ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig("pod_v10_m9_compare_T1pct.png", dpi=120, bbox_inches="tight")
plt.show()

# --- 图 B：bg+5σ ---
fig, ax = plt.subplots(figsize=(7.2, 4.8))
for n in N_SHOTS_LIST:
    s = _series[n]
    ax.plot(s["bg"], s["y5"], "-", color=_COLORS_N[n], lw=2.0, label=f"N={n}")
ax.set_xlabel("bg"); ax.set_ylabel("bg + 5·std(peak)")
ax.set_title("对比图 B　bg+5·std(peak) vs bg（N=1/2/4）")
ax.legend(); ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig("pod_v10_m9_compare_bg5std.png", dpi=120, bbox_inches="tight")
plt.show()

# --- 图 C：peak 均值 ---
fig, ax = plt.subplots(figsize=(7.2, 4.8))
for n in N_SHOTS_LIST:
    s = _series[n]
    ax.plot(s["bg"], s["mu"], "-", color=_COLORS_N[n], lw=2.0, label=f"N={n}")
ax.set_xlabel("bg"); ax.set_ylabel("peak_mean")
ax.set_title("对比图 C　peak 均值 vs bg（N=1/2/4）")
ax.legend(); ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig("pod_v10_m9_compare_peak_mean.png", dpi=120, bbox_inches="tight")
plt.show()

print("T@1% − (bg+5σ) 摘要：")
for n in N_SHOTS_LIST:
    s = _series[n]
    d = s["T1"] - s["y5"]
    print(f"  N={n}: mean Δ={d.mean():+.2f}, max|Δ|={np.max(np.abs(d)):.2f}")
"""))

mod9.append(md_cell(r"""
### 9.3　固定信号强度，noise 线性增长

对若干固定 boost，扫单次 noise；peak 在信号窗的 `hist_add` 上统计。  
看：分布是否仅平移；均值/std 是否随 noise 线性。

数据由下方 cell 现场计算（或读 `CACHE_SIG`）；**不复用** v05/旧 peak_vs_noise 缓存。
"""))
mod9.append(code_cell(r"""
# ---- 模块 9.3：固定信号 × 扫 noise（hist_i 前缀和）----
BOOST_LIST_M9 = [0.0, 0.004, 0.008, 0.016, 0.032]
N_MC_SIG_M9 = 8000
MC_CHUNK_SIG_M9 = 1000
_SIG_SEED0 = 720_000

def _load_sig_cache(path):
    if not os.path.exists(path):
        return None
    z = np.load(path, allow_pickle=True)
    if (int(z["n_mc"]) != N_MC_SIG_M9
            or not np.allclose(z["boosts"], BOOST_LIST_M9)
            or not np.allclose(z["grid_key"], NOISE_GRID_AMB)
            or not np.array_equal(z["n_shots_list"], np.asarray(N_SHOTS_LIST))):
        print(f"信号缓存键不匹配，忽略 {path}")
        return None
    return z

_zsig = _load_sig_cache(CACHE_SIG)
if _zsig is None:
    _zsig = _load_sig_cache(CACHE_SIG_CKPT)

if _zsig is not None and np.all(_zsig["done"]):
    print(f"已载入信号缓存 {CACHE_SIG if os.path.exists(CACHE_SIG) else CACHE_SIG_CKPT}")
    SIG_M9 = {n: dict(
        peak_cnt=np.asarray(_zsig[f"peak_cnt_{n}"]),
        peak_mean=np.zeros((_zsig["boosts"].size, _zsig["grid_key"].size)),
        peak_std=np.zeros((_zsig["boosts"].size, _zsig["grid_key"].size)),
    ) for n in N_SHOTS_LIST}
    for n in N_SHOTS_LIST:
        for i in range(len(BOOST_LIST_M9)):
            for k in range(len(NOISE_GRID_AMB)):
                s = peak_stats_from_cnt(SIG_M9[n]["peak_cnt"][i, k])
                SIG_M9[n]["peak_mean"][i, k] = s["mean"]
                SIG_M9[n]["peak_std"][i, k] = s["std"]
else:
    print(f"开始模块 9.3 信号扫描：{len(NOISE_GRID_AMB)} noise × {len(BOOST_LIST_M9)} boost × "
          f"{N_MC_SIG_M9} MC（N_SHOTS_MAX={N_SHOTS_MAX} 前缀和）")
    _cnt = {n: np.zeros((len(BOOST_LIST_M9), len(NOISE_GRID_AMB), N_PIX_MACRO*n+2), dtype=np.int64)
            for n in N_SHOTS_LIST}
    _done = np.zeros(len(NOISE_GRID_AMB), dtype=bool)
    if _zsig is not None:
        _done = np.asarray(_zsig["done"], dtype=bool)
        for n in N_SHOTS_LIST:
            _cnt[n] = np.asarray(_zsig[f"peak_cnt_{n}"])
    _t0 = time.time()
    for k, nt in enumerate(NOISE_GRID_AMB):
        if _done[k]:
            continue
        r_det = float(r_det_for_noise(float(nt), N_PIX_MACRO))
        r_amb = r_det / PDE
        for ib, boost in enumerate(BOOST_LIST_M9):
            left, part = N_MC_SIG_M9, 0
            while left > 0:
                m = min(MC_CHUNK_SIG_M9, left)
                rng = np.random.default_rng(_SIG_SEED0 + 1009*k + 10007*ib + 104729*part)
                hi = binary_macro_stepping_per_shot(
                    m, F_VALS, N_SHOTS_MAX, R_SIG_UNIT_POD, TF_POD, r_amb, CENTERS_SIG,
                    rng, boost=float(boost))
                for n in N_SHOTS_LIST:
                    pk = hist_add_from_prefix(hi, n).max(axis=1)
                    _cnt[n][ib, k] += np.bincount(pk, minlength=_cnt[n].shape[-1])
                left -= m; part += 1
        _done[k] = True
        if (k % 4) == 0 or k == len(NOISE_GRID_AMB) - 1:
            _atomic_savez(CACHE_SIG_CKPT,
                          grid_key=NOISE_GRID_AMB, boosts=np.asarray(BOOST_LIST_M9, float),
                          n_mc=N_MC_SIG_M9, n_shots_list=np.asarray(N_SHOTS_LIST),
                          done=_done,
                          **{f"peak_cnt_{n}": _cnt[n] for n in N_SHOTS_LIST})
            print(f"  sig [{k+1}/{len(NOISE_GRID_AMB)}] noise={nt:.2f}  ({time.time()-_t0:.0f}s)")
    _atomic_savez(CACHE_SIG,
                  grid_key=NOISE_GRID_AMB, boosts=np.asarray(BOOST_LIST_M9, float),
                  n_mc=N_MC_SIG_M9, n_shots_list=np.asarray(N_SHOTS_LIST),
                  done=_done,
                  **{f"peak_cnt_{n}": _cnt[n] for n in N_SHOTS_LIST})
    if os.path.exists(CACHE_SIG_CKPT):
        try: os.remove(CACHE_SIG_CKPT)
        except OSError: pass
    print(f"已写入 {CACHE_SIG}")
    SIG_M9 = {n: dict(peak_cnt=_cnt[n],
                      peak_mean=np.zeros((len(BOOST_LIST_M9), len(NOISE_GRID_AMB))),
                      peak_std=np.zeros((len(BOOST_LIST_M9), len(NOISE_GRID_AMB))))
              for n in N_SHOTS_LIST}
    for n in N_SHOTS_LIST:
        for i in range(len(BOOST_LIST_M9)):
            for k in range(len(NOISE_GRID_AMB)):
                s = peak_stats_from_cnt(SIG_M9[n]["peak_cnt"][i, k])
                SIG_M9[n]["peak_mean"][i, k] = s["mean"]
                SIG_M9[n]["peak_std"][i, k] = s["std"]

_E_NJ = np.asarray(BOOST_LIST_M9) * E_PULSE_BASE * 1e9
_noise_t = NOISE_GRID_AMB

fig, ax = plt.subplots(2, len(N_SHOTS_LIST), figsize=(5.2*len(N_SHOTS_LIST), 8.0), sharex=True)
for j, n in enumerate(N_SHOTS_LIST):
    d = SIG_M9[n]
    for i, e in enumerate(_E_NJ):
        ls = "-" if BOOST_LIST_M9[i] > 0 else "--"
        ax[0, j].plot(_noise_t, d["peak_mean"][i], ls=ls, lw=1.4, label=f"E={e:.2f}nJ")
        ax[1, j].plot(_noise_t, d["peak_std"][i], ls=ls, lw=1.4, label=f"E={e:.2f}nJ")
    ax[0, j].set_title(f"N={n} peak均值 vs noise"); ax[1, j].set_title(f"N={n} peak std vs noise")
    ax[1, j].set_xlabel("noise（单次）")
    ax[0, j].legend(fontsize=7); ax[1, j].legend(fontsize=7)
    ax[0, j].grid(alpha=0.3); ax[1, j].grid(alpha=0.3)
fig.suptitle("模块 9.3a　固定信号：peak 均值/std 随 noise", fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.95])
fig.savefig("pod_v10_m9_sig_mean_std.png", dpi=120, bbox_inches="tight")
plt.show()

# 分布平移检验：N=4、中等信号
_n, _ib = 4, min(2, len(BOOST_LIST_M9)-1)
fig, ax = plt.subplots(1, 2, figsize=(13.0, 4.6))
_ks = np.unique(np.linspace(0, len(_noise_t)-1, 6).astype(int))
for k in _ks:
    cnt = SIG_M9[_n]["peak_cnt"][_ib, k]
    s = peak_stats_from_cnt(cnt)
    x = np.arange(cnt.size); pmf = cnt / max(cnt.sum(), 1)
    ax[0].step(x, pmf, where="mid", lw=1.2, label=f"noise={_noise_t[k]:.2f}")
    ax[1].step(x - s["mean"], pmf, where="mid", lw=1.2, label=f"noise={_noise_t[k]:.2f}")
ax[0].set_title(f"N={_n}, E={_E_NJ[_ib]:.2f}nJ 原始分布"); ax[0].legend(fontsize=7); ax[0].grid(alpha=0.3)
ax[1].set_title("中心化（重叠 => 近似纯平移）"); ax[1].legend(fontsize=7); ax[1].grid(alpha=0.3)
fig.suptitle("模块 9.3b　分布平移检验", fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.93])
fig.savefig("pod_v10_m9_sig_shift.png", dpi=120, bbox_inches="tight")
plt.show()

print(f"{'N':>3}{'E[nJ]':>10}{'slope':>10}{'R2':>8}")
for n in N_SHOTS_LIST:
    for i, e in enumerate(_E_NJ):
        y = SIG_M9[n]["peak_mean"][i]
        coef = np.polyfit(_noise_t, y, 1)
        yhat = np.polyval(coef, _noise_t)
        r2 = 1 - np.sum((y-yhat)**2) / max(np.sum((y-y.mean())**2), 1e-30)
        print(f"{n:>3}{e:>10.2f}{coef[0]:>10.3f}{r2:>8.4f}")
"""))

nb["cells"].extend(mod9)

# 全局替换残余 v05 缓存名提示
for i, c in enumerate(nb["cells"]):
    if c["cell_type"] != "code":
        continue
    s = src(i)
    s2 = s.replace("pod_esti_v05_cache", "pod_esti_v10_cache")
    s2 = s2.replace("pod_v05_", "pod_v10_")
    if s2 != s:
        set_src(i, s2)

clear_all_outputs()
json.dump(nb, open(PATH, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
print(f"wrote {PATH}, cells={len(nb['cells'])}")

# 语法检查关键 cell
import ast
for i in [2, 9, 17, 25] + list(range(len(nb["cells"]) - 7, len(nb["cells"]))):
    if nb["cells"][i]["cell_type"] != "code":
        continue
    try:
        ast.parse(src(i))
        print(f"cell {i}: OK")
    except SyntaxError as e:
        print(f"cell {i}: SYNTAX {e}")
        lines = src(i).splitlines()
        ln = e.lineno or 1
        for j in range(max(0, ln-2), min(len(lines), ln+2)):
            print(f"  {j+1}| {lines[j]}")
