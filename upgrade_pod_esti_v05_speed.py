# -*- coding: utf-8 -*-
"""加速 PoD_esti_v05：外层多档并行 + 验证批量化 + MC 分块喂满线程。"""
import json
import re
from pathlib import Path

nb_path = Path("PoD_esti_v05.ipynb")
nb = json.loads(nb_path.read_text(encoding="utf-8"))


def set_src(cell, text: str):
    cell["source"] = text.splitlines(keepends=True)


# ---------- cell 2 ----------
c2 = "".join(nb["cells"][2]["source"])
if "POD_BIN_WORKERS" not in c2:
    c2 = re.sub(
        r"N_WORKERS\s*=\s*20[^\n]*\nMC_CHUNK\s*=\s*5_000[^\n]*\nNOISE_WORKERS\s*=\s*N_WORKERS\n",
        "N_WORKERS   = 20         # 全机线程预算\n"
        "MC_CHUNK    = 5_000      # 纯噪声分块；20 线程时控制峰值内存\n"
        "NOISE_WORKERS = N_WORKERS\n"
        "POD_BIN_WORKERS = 4      # PoD 外层同时跑几档 noise（×内层≈20）\n"
        "POD_MC_CHUNK = 250       # 信号 MC 再切块，喂满内层线程\n",
        c2, count=1,
    )
    c2 = re.sub(
        r"POD_WORKERS\s*=\s*N_WORKERS[^\n]*\n",
        "POD_WORKERS = max(1, N_WORKERS // POD_BIN_WORKERS)  # 每档内 boost/分块并行\n",
        c2, count=1,
    )
else:
    print("cell2 already has POD_BIN_WORKERS")

if "POD_BIN_WORKERS={POD_BIN_WORKERS}" not in c2:
    c2 = re.sub(
        r'print\(f"  ★ 并行：N_WORKERS=\{N_WORKERS\}，MC_CHUNK=\{MC_CHUNK:,\}，每 \{CHECKPOINT_EVERY\} 档增量落盘"\)',
        'print(f"  ★ 并行：N_WORKERS={N_WORKERS}，噪声分块={MC_CHUNK:,}；'
        'PoD 外层×内层={POD_BIN_WORKERS}×{POD_WORKERS}，POD_MC_CHUNK={POD_MC_CHUNK}；'
        '每 {CHECKPOINT_EVERY} 档增量落盘")',
        c2, count=1,
    )
set_src(nb["cells"][2], c2)
assert "POD_BIN_WORKERS" in c2 and "POD_MC_CHUNK" in c2
print("patched cell 2")


# ---------- import threading ----------
for i, c in enumerate(nb["cells"]):
    if c.get("cell_type") != "code":
        continue
    s = "".join(c["source"])
    if "ThreadPoolExecutor" in s and "from concurrent.futures" in s:
        if "import threading" not in s:
            s = s.replace(
                "from concurrent.futures import ThreadPoolExecutor, as_completed",
                "from concurrent.futures import ThreadPoolExecutor, as_completed\nimport threading",
            )
            set_src(nb["cells"][i], s)
            print("added import threading in cell", i)
        break


# ---------- cell 25 functions ----------
NEW_POD_FUNCS = r'''print(f"PoD 子窗：{POD_T_LO*1e9:.1f}–{POD_T_HI*1e9:.1f} ns，{TF_POD.size} 个细网格步")
print(f"每种 N_shots 对自己的完整 NOISE_GRID 求解："
      f"{[len(NOISE_GRID[n]) for n in N_SHOTS_LIST]} 档，noise 步长 0.25")
print(f"并行：外层 POD_BIN_WORKERS={POD_BIN_WORKERS} × 内层 POD_WORKERS={POD_WORKERS}；"
      f"MC 分块 POD_MC_CHUNK={POD_MC_CHUNK}；临界验证 {N_MC_POD_VERIFY:,} 次/点")


def _peaks_chunk(boost, n_shots, r_amb, n_real, seed):
    """单块 MC，无内部并行。"""
    f_arr = np.tile(F_VALS, n_shots)
    h = binary_macro_stepping(
        n_real, f_arr, R_SIG_UNIT_POD, TF_POD, r_amb, CENTERS_SIG,
        np.random.default_rng(seed), boost=boost,
    )
    return h.max(axis=1)


def sig_peaks(boost, n_shots, r_amb, n_real, seed):
    """兼容入口：大 n_real 时按 POD_MC_CHUNK 切开，用 POD_WORKERS 并行。"""
    n_real = int(n_real)
    if n_real <= POD_MC_CHUNK or POD_WORKERS <= 1:
        return _peaks_chunk(boost, n_shots, r_amb, n_real, seed)
    specs = []
    for s in range(0, n_real, POD_MC_CHUNK):
        m = min(POD_MC_CHUNK, n_real - s)
        specs.append((boost, n_shots, r_amb, m, seed + 104729 * s))
    with ThreadPoolExecutor(max_workers=POD_WORKERS) as pool:
        parts = list(pool.map(lambda sp: _peaks_chunk(*sp), specs))
    return np.concatenate(parts)


def _eval_mc_jobs(job_specs, n_shots, r_amb):
    """统一并行入口。
    job_specs: [(boost, n_real, seed), ...] → 与输入等长的 peak 数组列表。
    每个 job 再按 POD_MC_CHUNK 切开，全部丢进同一个线程池，避免嵌套池。
    """
    flat, owners = [], []
    for j, (boost, n_real, seed) in enumerate(job_specs):
        n_real = int(n_real)
        if n_real <= 0:
            continue
        for s in range(0, n_real, POD_MC_CHUNK):
            m = min(POD_MC_CHUNK, n_real - s)
            flat.append((float(boost), n_shots, r_amb, m, int(seed) + 104729 * s))
            owners.append(j)
    out = [None] * len(job_specs)
    if not flat:
        return [np.zeros(0, dtype=int) for _ in job_specs]
    if POD_WORKERS <= 1:
        parts = [_peaks_chunk(*sp) for sp in flat]
    else:
        with ThreadPoolExecutor(max_workers=POD_WORKERS) as pool:
            parts = list(pool.map(lambda sp: _peaks_chunk(*sp), flat))
    buckets = [[] for _ in job_specs]
    for own, pk in zip(owners, parts):
        buckets[own].append(pk)
    for j, segs in enumerate(buckets):
        out[j] = np.concatenate(segs) if segs else np.zeros(0, dtype=int)
    return out


def _eval_boost_grid(boosts, n_shots, r_amb, n_real, seed0):
    """并行评估若干独立能量点；返回每点的 peak 样本。"""
    boosts = np.asarray(boosts, float)
    jobs = [(float(b), n_real, seed0 + 1009 * i) for i, b in enumerate(boosts)]
    return _eval_mc_jobs(jobs, n_shots, r_amb)


def _isotonic(p):
    """简单保序：消除有限 MC 导致的局部 PoD 下降。"""
    return np.maximum.accumulate(np.asarray(p, float))


def _crossing_logboost(boosts, pod, level):
    """在 log10(boost) 上找首次跨越；越界返回 NaN。"""
    order = np.argsort(boosts)
    x = np.log10(np.asarray(boosts)[order])
    p = _isotonic(np.asarray(pod)[order])
    if p[0] >= level:
        return float(x[0])
    if p[-1] < level:
        return np.nan
    i = int(np.searchsorted(p, level))
    dp = p[i] - p[i - 1]
    w = 0.5 if dp <= 0 else (level - p[i - 1]) / dp
    return float(x[i - 1] + w * (x[i] - x[i - 1]))


def _probit_fit(boosts, pod, n_real):
    """拟合 Phi^-1(PoD) = slope*log10(boost) + intercept。"""
    boosts = np.asarray(boosts, float)
    success = np.rint(np.asarray(pod, float) * n_real)
    p = (success + 0.5) / (n_real + 1.0)
    transition = (p > 0.01) & (p < 0.99)
    if transition.sum() < 3:
        transition = np.argsort(np.abs(p - 0.5))[:min(5, len(p))]
    x = np.log10(boosts[transition])
    z = _norm.ppf(p[transition])
    slope, intercept = np.polyfit(x, z, 1)
    return float(slope), float(intercept)


def _pk_to_record(boost, pk, T, n_shots):
    return {
        "boost": float(boost),
        "pod": float((pk >= T).mean()),
        "peak_mean": float(pk.mean()),
        "peak_std": float(pk.std()),
        "peak_cnt": np.bincount(pk, minlength=N_PIX_MACRO * n_shots + 1),
        "n_verify": int(pk.size),
    }


def _verify_critical_batch(cands, n_shots, r_amb, seed0):
    """批量验证全部 FAR×PoD 临界点（两轮内可并行），避免逐点串行 5000 次 MC。
    cands: list of dict(tag, level, boost, T, slope)
    返回 {(tag, level_key): record}
    """
    if not cands:
        return {}
    boosts = [c["boost"] for c in cands]
    pks = _eval_mc_jobs(
        [(b, N_MC_POD_VERIFY, seed0 + 7919 * i) for i, b in enumerate(boosts)],
        n_shots, r_amb,
    )
    active = []
    finals = {}
    for i, c in enumerate(cands):
        rec = _pk_to_record(boosts[i], pks[i], c["T"], n_shots)
        level, slope = c["level"], c["slope"]
        key = (c["tag"], f"{level:.2f}")
        if abs(rec["pod"] - level) <= POD_VERIFY_TOL or slope <= 0:
            finals[key] = rec
            continue
        p_smooth = (int((pks[i] >= c["T"]).sum()) + 0.5) / (pks[i].size + 1.0)
        dx = (_norm.ppf(level) - _norm.ppf(p_smooth)) / slope
        new_boost = float(boosts[i] * 10.0 ** float(np.clip(dx, -0.25, 0.25)))
        active.append({**c, "boost": new_boost, "_key": key, "_i": i})
    if active:
        boosts2 = [a["boost"] for a in active]
        pks2 = _eval_mc_jobs(
            [(b, N_MC_POD_VERIFY, seed0 + 1_000_003 + 7919 * a["_i"])
             for a, b in zip(active, boosts2)],
            n_shots, r_amb,
        )
        for a, b, pk in zip(active, boosts2, pks2):
            finals[a["_key"]] = _pk_to_record(b, pk, a["T"], n_shots)
    return finals


def solve_pod_noise(n_shots, k, seed0):
    """求一个 noise 档、全部 FAR 阈值下的 PoD50/90 临界点。"""
    R, Tr = NOISE_RES[n_shots], THRESH[n_shots]
    nt = float(R["noise_target"][k])
    n_tr = int(R["n_tr"])
    r_amb = float(R["r_det"][k] / PDE)
    T_map = {FAR_TAG[far]: int(Tr["T" + FAR_TAG[far]][k]) for far in TARGET_FARS}
    if max(T_map.values()) > n_tr:
        return (n_shots, nt), {
            "noise": float(R["noise_mc"][k]), "e_lambda": float(R["e_lambda"][k]),
            "n_tr": n_tr, "T_map": T_map, "critical": {}, "invalid": "阈值超过二值硬上限",
        }

    coarse_boost = np.logspace(POD_LOG_BOOST_MIN, POD_LOG_BOOST_MAX, N_POD_COARSE)
    coarse_pk = _eval_boost_grid(
        coarse_boost, n_shots, r_amb, N_MC_POD_COARSE, seed0,
    )
    coarse_pod = {
        tag: np.array([(pk >= T).mean() for pk in coarse_pk])
        for tag, T in T_map.items()
    }

    roots0 = []
    for tag in T_map:
        for level in POD_LEVELS:
            x0 = _crossing_logboost(coarse_boost, coarse_pod[tag], level)
            if np.isfinite(x0):
                roots0.append(x0)
    if roots0:
        local_x = np.unique(np.concatenate([
            np.linspace(x0 - POD_LOCAL_HALF_DECADE, x0 + POD_LOCAL_HALF_DECADE,
                        N_POD_LOCAL_PER_ROOT)
            for x0 in roots0
        ]))
        local_boost = 10.0 ** local_x
        local_pk = _eval_boost_grid(
            local_boost, n_shots, r_amb, N_MC_POD_LOCAL, seed0 + 500_000,
        )
    else:
        local_boost = np.array([], float)
        local_pk = []

    critical = {tag: {} for tag in T_map}
    curve = {}
    cands = []
    for tag, T in T_map.items():
        boosts_fit = np.concatenate([coarse_boost, local_boost])
        pod_fit = np.concatenate([
            coarse_pod[tag],
            np.array([(pk >= T).mean() for pk in local_pk]) if len(local_pk) else np.array([], float),
        ])
        order = np.argsort(boosts_fit)
        boosts_fit, pod_fit = boosts_fit[order], pod_fit[order]
        curve[tag] = {"boost": boosts_fit, "pod": pod_fit}
        slope, intercept = _probit_fit(boosts_fit, pod_fit, N_MC_POD_LOCAL)
        for level in POD_LEVELS:
            x_root = (_norm.ppf(level) - intercept) / slope if slope > 0 else np.nan
            if not np.isfinite(x_root):
                critical[tag][f"{level:.2f}"] = None
                continue
            cands.append({
                "tag": tag, "level": level, "T": T, "slope": slope,
                "boost": float(10.0 ** x_root),
            })

    verified = _verify_critical_batch(cands, n_shots, r_amb, seed0 + 700_000)
    for (tag, lk), rec in verified.items():
        critical[tag][lk] = rec
    for tag in T_map:
        for level in POD_LEVELS:
            critical[tag].setdefault(f"{level:.2f}", None)

    return (n_shots, nt), {
        "noise": float(R["noise_mc"][k]),
        "noise_target": nt,
        "e_lambda": float(R["e_lambda"][k]),
        "n_tr": n_tr,
        "T_map": T_map,
        "curve": curve,
        "critical": critical,
    }
'''

c25 = "".join(nb["cells"][25]["source"])
start = c25.find('print(f"PoD 子窗：')
mid = c25.find("# ---- 对完整 0.25-noise 网格求解")
if start < 0 or mid < 0:
    raise SystemExit(f"markers not found start={start} mid={mid}")
c25 = c25[:start] + NEW_POD_FUNCS + "\n\n" + c25[mid:]

# Replace main loop — match by unique anchors
loop_start = c25.find("    _tall = time.time()")
loop_end = c25.find("    if os.path.exists(CACHE_POD_CKPT):")
if loop_start < 0 or loop_end < 0:
    raise SystemExit(f"loop markers not found {loop_start} {loop_end}")
# include the ckpt removal block
loop_end2 = c25.find("\n\n", c25.find("except OSError:", loop_end))
if loop_end2 < 0:
    # take through pass of remove
    loop_end2 = c25.find("            pass\n", loop_end)
    if loop_end2 < 0:
        raise SystemExit("cannot find end of loop")
    loop_end2 += len("            pass\n")

NEW_LOOP = '''    _tall = time.time()
    n_total = len(_expected_keys)
    print(f"PoD 临界点扫描：共 {n_total} 个 noise 档，步长 0.25；"
          f"每档求 {len(FAR_TAGS)} 档 FAR × PoD50/90；"
          f"外层×内层 = {POD_BIN_WORKERS}×{POD_WORKERS}，MC 分块 {POD_MC_CHUNK}")
    pending = []
    for n_shots in N_SHOTS_LIST:
        for k in range(len(NOISE_GRID[n_shots])):
            key = (n_shots, float(NOISE_GRID[n_shots][k]))
            if key in POD_RES:
                continue
            seed = 7000 + n_shots * 1_000_000 + k * 20_000
            pending.append((n_shots, k, seed))
    _ckpt_lock = threading.Lock()
    _done_n = 0

    def _progress_msg(value):
        c100 = value.get("critical", {}).get("100ppm", {})
        p90 = c100.get("0.90")
        if not p90:
            return "无有效交点"
        return (f"E90={p90['boost']*E_PULSE_BASE*1e9:.3g} nJ，"
                f"验证PoD={p90['pod']:.3f}，peak均值={p90['peak_mean']:.2f}")

    def _consume(n_shots, k, seed):
        return solve_pod_noise(n_shots, k, seed)

    if not pending:
        print("无需计算：待跑列表为空")
    elif POD_BIN_WORKERS <= 1 or len(pending) <= 1:
        for n_shots, k, seed in pending:
            key, value = _consume(n_shots, k, seed)
            POD_RES[key] = value
            _done_n += 1
            if (_done_n % CHECKPOINT_EVERY) == 0:
                _save_pod_cache(CACHE_POD_CKPT, POD_RES, _pod_grid_key)
            if _done_n == 1 or _done_n % 5 == 0 or _done_n == len(pending):
                elapsed = time.time() - _tall
                print(f"  [N_shots={n_shots} noise={key[1]:.2f}] "
                      f"{_progress_msg(value)}；累计 {elapsed/60:.1f} min；"
                      f"已完成 {len(POD_RES)}/{n_total}")
    else:
        with ThreadPoolExecutor(max_workers=POD_BIN_WORKERS) as pool:
            futs = {
                pool.submit(_consume, ns, k, seed): (ns, k)
                for ns, k, seed in pending
            }
            for fut in as_completed(futs):
                key, value = fut.result()
                with _ckpt_lock:
                    POD_RES[key] = value
                    _done_n += 1
                    if (_done_n % CHECKPOINT_EVERY) == 0:
                        _save_pod_cache(CACHE_POD_CKPT, POD_RES, _pod_grid_key)
                    if _done_n == 1 or _done_n % 5 == 0 or _done_n == len(pending):
                        elapsed = time.time() - _tall
                        print(f"  [N_shots={key[0]} noise={key[1]:.2f}] "
                              f"{_progress_msg(value)}；累计 {elapsed/60:.1f} min；"
                              f"已完成 {len(POD_RES)}/{n_total}")
    _save_pod_cache(CACHE_POD, POD_RES, _pod_grid_key)
    print(f"逐 noise PoD 扫描总用时 {(time.time()-_tall)/60:.1f} min；已写入 {CACHE_POD}")
    if os.path.exists(CACHE_POD_CKPT):
        try:
            os.remove(CACHE_POD_CKPT)
        except OSError:
            pass
'''

c25 = c25[:loop_start] + NEW_LOOP + c25[loop_end2:]
set_src(nb["cells"][25], c25)
print("patched cell 25")

for i in [2, 25]:
    compile("".join(nb["cells"][i]["source"]), f"cell_{i}", "exec")
print("syntax OK")

nb_path.write_text(json.dumps(nb, ensure_ascii=False, indent=1), encoding="utf-8")
print("written", nb_path)
