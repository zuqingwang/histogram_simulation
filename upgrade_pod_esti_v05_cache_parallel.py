# -*- coding: utf-8 -*-
"""升级 PoD_esti_v05：20 线程并行 + 耗时 MC 增量落盘。不执行仿真。"""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
NB = ROOT / "PoD_esti_v05.ipynb"

# ---- 新参数块片段 ----
PARAM_OLD = '''MC_CHUNK    = 12_500     # 8 线程并行时减小分块，控制峰值内存
NOISE_WORKERS = 8       # 纯噪声 MC 分块并行

# ---- 第 3 步：能量扫描与 PoD ----
#   ★ v02：PoD 的环境档直接按【目标 noise】指定（两种 N_shots 用同一组 noise，便于公平对比）
# ---- v04：完整 noise 网格上的自适应 PoD 交点 ----
NOISE_POD = {n: NOISE_GRID[n].copy() for n in N_SHOTS_LIST}
POD_WORKERS = 8          # 信号能量点并行评估'''

PARAM_NEW = '''# ---- ★ v05：并行与缓存策略（CPU 支持 20 线程）----
N_WORKERS   = 20         # 纯噪声 / PoD 统一线程数
MC_CHUNK    = 5_000      # 20 线程并发时减小分块，控制峰值内存
NOISE_WORKERS = N_WORKERS
CHECKPOINT_EVERY = 1     # 每完成几档就落盘一次（防中断丢失）

# ---- 第 3 步：能量扫描与 PoD ----
# ---- v04：完整 noise 网格上的自适应 PoD 交点 ----
NOISE_POD = {n: NOISE_GRID[n].copy() for n in N_SHOTS_LIST}
POD_WORKERS = N_WORKERS  # 信号能量点并行评估'''

CACHE_OLD = '''# ---- 缓存 ----
USE_CACHE   = True
CACHE_NOISE = "pod_esti_v04_cache_noise.npz"
CACHE_POD   = "pod_esti_v04_cache_pod.npz"'''

CACHE_NEW = '''# ---- ★ v05 缓存：主文件 + 兼容读取旧版 + 增量检查点 ----
USE_CACHE = True
CACHE_NOISE = "pod_esti_v05_cache_noise.npz"
CACHE_POD   = "pod_esti_v05_cache_pod.npz"
CACHE_NOISE_FALLBACK = ["pod_esti_v04_cache_noise.npz"]
CACHE_POD_FALLBACK   = ["pod_esti_v04_cache_pod.npz"]
CACHE_NOISE_CKPT = "pod_esti_v05_cache_noise.partial.npz"
CACHE_POD_CKPT   = "pod_esti_v05_cache_pod.partial.npz"'''

PRINT_EXTRA = '''print(f"  ★ 噪点率目标：{[f'{f*1e6:.0f} ppm' for f in TARGET_FARS]}，每档 {N_MC_NOISE:,} 次 MC")'''
PRINT_EXTRA_NEW = '''print(f"  ★ 噪点率目标：{[f'{f*1e6:.0f} ppm' for f in TARGET_FARS]}，每档 {N_MC_NOISE:,} 次 MC")
print(f"  ★ 并行：N_WORKERS={N_WORKERS}，MC_CHUNK={MC_CHUNK:,}，每 {CHECKPOINT_EVERY} 档增量落盘")
print(f"  ★ 缓存主文件：{CACHE_NOISE} / {CACHE_POD}（可读 fallback：v04）")'''

NOISE_SCAN_CELL = r'''def _noise_chunk_stats(m, n_tr, r_det, inv_tab, seed):
    """单个纯噪声分块；返回可直接归并的充分统计量。"""
    h = noise_macro_hist_fast(
        m, n_tr, r_det, np.random.default_rng(seed), inv_tab=inv_tab,
    )
    a = h[:, I_STAT0:I_STAT1]
    nz = a.mean(axis=1)
    return (float(nz.sum()), float((nz*nz).sum()),
            np.bincount(a.max(axis=1), minlength=n_tr + 2))


def peak_stats_from_cnt(cnt):
    """由 peak 的 bincount 精确算出均值/标准差/各分位数。"""
    v = np.arange(cnt.size, dtype=float)
    n = cnt.sum()
    mean = (v * cnt).sum() / n
    var = (v*v * cnt).sum() / n - mean**2
    cum = np.cumsum(cnt) / n
    q = lambda p: float(np.searchsorted(cum, p))
    return dict(n=int(n), mean=mean, std=np.sqrt(max(var, 0.0)),
                p01=q(0.01), p50=q(0.50), p99=q(0.99),
                p999=q(0.999), p9999=q(0.9999),
                pmax=float(np.nonzero(cnt)[0].max()) if np.any(cnt) else 0.0)


def _atomic_savez(path, **kwargs):
    """先写临时文件再替换，避免中断写出半截缓存。
    注意：np.savez_compressed 若路径不以 .npz 结尾会自动追加 .npz，
    因此临时文件必须自带 .npz 后缀，否则 os.replace 会找不到文件（WinError 2）。
    """
    path = str(path)
    tmp = path + ".tmp.npz"
    np.savez_compressed(tmp, **kwargs)
    os.replace(tmp, path)


def _try_load_noise_cache(path, grid_key):
    if not os.path.exists(path):
        return None
    z = np.load(path, allow_pickle=True)
    if (int(z["n_mc"]) == N_MC_NOISE and list(z["n_shots_list"]) == list(N_SHOTS_LIST)
            and z["grid_key"].shape == grid_key.shape
            and np.allclose(z["grid_key"], grid_key)):
        return z["res"].item()
    return None


def _save_noise_cache(path, res, grid_key):
    _atomic_savez(path, res=np.array(res, dtype=object),
                  n_mc=N_MC_NOISE, n_shots_list=np.array(N_SHOTS_LIST),
                  grid_key=grid_key)


def run_noise_scan(n_shots, noise_grid, n_mc, chunk, seed0=2000, verbose_every=5,
                   res=None, start_k=0, on_progress=None):
    """对一组【目标 noise】跑纯噪声 MC；支持从 start_k 断点续跑。"""
    n_tr = N_PIX_MACRO * n_shots
    ng = len(noise_grid)
    if res is None:
        res = {"n_shots": n_shots, "n_tr": n_tr,
               "noise_target": np.asarray(noise_grid, float),
               "r_det": np.zeros(ng), "e_lambda": np.zeros(ng),
               "p_eq": np.zeros(ng), "noise_mc": np.zeros(ng), "noise_std": np.zeros(ng),
               "peak_cnt": np.zeros((ng, n_tr + 2), dtype=np.int64),
               "done": np.zeros(ng, dtype=bool)}
    elif "done" not in res:
        # 兼容旧缓存：peak_cnt 有样本即视为已完成
        res["done"] = np.array([int(c.sum()) > 0 for c in res["peak_cnt"]], dtype=bool)

    t_start = time.time()
    for k, nt in enumerate(noise_grid):
        if k < start_k or bool(res["done"][k]):
            continue
        r_det = r_det_for_noise(float(nt), n_tr)
        res["r_det"][k] = r_det
        res["e_lambda"][k] = e_lambda_for_r_det(r_det)
        res["p_eq"][k] = p_bin_equilibrium(r_det)[0]
        inv_tab = build_inv_table(r_det)
        s1 = s2 = 0.0
        res["peak_cnt"][k][:] = 0
        specs = [
            (min(chunk, n_mc - s), n_tr, r_det, inv_tab, seed0 + 1000*k + s)
            for s in range(0, n_mc, chunk)
        ]
        if NOISE_WORKERS <= 1:
            parts = [_noise_chunk_stats(*spec) for spec in specs]
        else:
            with ThreadPoolExecutor(max_workers=NOISE_WORKERS) as pool:
                parts = list(pool.map(lambda x: _noise_chunk_stats(*x), specs))
        for p1, p2, pcnt in parts:
            s1 += p1; s2 += p2
            res["peak_cnt"][k] += pcnt
        res["noise_mc"][k] = s1 / n_mc
        res["noise_std"][k] = np.sqrt(max(s2/n_mc - (s1/n_mc)**2, 0.0))
        res["done"][k] = True
        if on_progress is not None:
            on_progress(res, k)
        if k == 0 or k == ng-1 or (k+1) % verbose_every == 0:
            el = time.time() - t_start
            remain = max(int((~res["done"]).sum()), 0)
            done_n = int(res["done"].sum())
            eta = el / max(done_n - start_k, 1) * remain
            pk = peak_stats_from_cnt(res["peak_cnt"][k])
            print(f"  [N_shots={n_shots} {done_n:>3d}/{ng}] 目标 noise={nt:>6.2f} → "
                  f"实测 {res['noise_mc'][k]:>6.3f}（E_λ={res['e_lambda'][k]:.4f}，"
                  f"≈{res['e_lambda'][k]/0.68*100:>5.0f} klux）  "
                  f"peak 中位={pk['p50']:>5.1f} 99.99%={pk['p9999']:>5.1f}  "
                  f"[已用 {el:.0f}s, 剩约 {eta:.0f}s]")
    return res


# ---- 估算总耗时并开跑（主缓存 + fallback + 增量检查点）----
_grid_key = np.concatenate([np.asarray(NOISE_GRID[n]) for n in N_SHOTS_LIST])
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
            # 兼容 v04：有 peak 样本即视为完成
            if not all(int(c.sum()) > 0 for c in r["peak_cnt"]):
                return False
            r["done"] = np.ones(len(r["noise_target"]), dtype=bool)
    return True


if NOISE_RES is not None and _noise_is_complete(NOISE_RES):
    print(f"已从缓存 {_loaded_from} 载入纯噪声 MC 结果（每档 {N_MC_NOISE:,} 条，完整）")
    if _loaded_from != CACHE_NOISE:
        _save_noise_cache(CACHE_NOISE, NOISE_RES, _grid_key)
        print(f"已同步写入主缓存 {CACHE_NOISE}")
else:
    if NOISE_RES is None:
        NOISE_RES = {}
        print("未找到匹配的完整/部分噪声缓存，开始全新扫描")
    else:
        print(f"从 {_loaded_from} 载入部分结果，断点续跑")
    _est = sum((8 + 4.2*float(nt)) * (N_MC_NOISE/1e6)
               for n in N_SHOTS_LIST for nt in NOISE_GRID[n])
    print(f"纯噪声 MC：{len(N_SHOTS_LIST)} 种 N_shots × 各 "
          f"{[len(NOISE_GRID[n]) for n in N_SHOTS_LIST]} 档 × {N_MC_NOISE:,} 条")
    # 经验：NumPy 内存带宽受限，20 线程有效加速大约 3–5×，不是线性
    _parallel_eff = min(N_WORKERS, 1.0 + 0.22*max(N_WORKERS-1, 0))
    print(f"  单线程基准约 {_est/60:.0f} 分钟；{N_WORKERS} 线程预计 "
          f"{_est/_parallel_eff/60:.0f} 分钟（内存带宽受限，非线性加速）")

    def _ckpt_noise(res_all):
        _save_noise_cache(CACHE_NOISE_CKPT, res_all, _grid_key)

    _tall = time.time()
    _ckpt_counter = {"n": 0}

    def _on_progress(n_shots):
        def _cb(res, k):
            NOISE_RES[n_shots] = res
            _ckpt_counter["n"] += 1
            if (_ckpt_counter["n"] % CHECKPOINT_EVERY) == 0:
                _ckpt_noise(NOISE_RES)
        return _cb

    for _ns in N_SHOTS_LIST:
        _prev = NOISE_RES.get(_ns)
        NOISE_RES[_ns] = run_noise_scan(
            _ns, NOISE_GRID[_ns], N_MC_NOISE, MC_CHUNK,
            res=_prev, on_progress=_on_progress(_ns),
        )
        _ckpt_noise(NOISE_RES)
        print(f"  N_shots={_ns} 已检查点写入 {CACHE_NOISE_CKPT}")
    print(f"总用时 {time.time()-_tall:.0f} s")
    _save_noise_cache(CACHE_NOISE, NOISE_RES, _grid_key)
    print(f"已写入主缓存 {CACHE_NOISE}")
    if os.path.exists(CACHE_NOISE_CKPT):
        try:
            os.remove(CACHE_NOISE_CKPT)
        except OSError:
            pass
'''

POD_TAIL_OLD_START = "# ---- 对完整 0.25-noise 网格求解；缓存键包含全部精度参数 ----"

POD_TAIL_NEW = r'''# ---- 对完整 0.25-noise 网格求解；主缓存 + fallback + 增量检查点 ----
_pod_grid_key = np.concatenate([NOISE_GRID[n] for n in N_SHOTS_LIST])


def _try_load_pod_cache(path, grid_key):
    if not os.path.exists(path):
        return None
    z = np.load(path, allow_pickle=True)
    if (
        np.array_equal(z["grid_key"], grid_key)
        and int(z["n_coarse"]) == N_MC_POD_COARSE
        and int(z["n_local"]) == N_MC_POD_LOCAL
        and int(z["n_verify"]) == N_MC_POD_VERIFY
    ):
        return z["res"].item()
    return None


def _save_pod_cache(path, res, grid_key):
    # 元组键无法直接进 npz；统一转成可 pickle 的 dict
    _atomic_savez(
        path, res=np.array(res, dtype=object),
        grid_key=grid_key,
        n_coarse=N_MC_POD_COARSE, n_local=N_MC_POD_LOCAL,
        n_verify=N_MC_POD_VERIFY,
    )


POD_RES = None
_loaded_pod_from = None
for _cand in [CACHE_POD, *CACHE_POD_FALLBACK, CACHE_POD_CKPT]:
    POD_RES = _try_load_pod_cache(_cand, _pod_grid_key)
    if POD_RES is not None:
        _loaded_pod_from = _cand
        break

_expected_keys = {(ns, float(nt)) for ns in N_SHOTS_LIST for nt in NOISE_GRID[ns]}
_have_keys = set(POD_RES.keys()) if POD_RES else set()
_complete = POD_RES is not None and _expected_keys.issubset(_have_keys)

if _complete:
    print(f"已从缓存 {_loaded_pod_from} 载入逐 noise PoD 临界点（完整 "
          f"{len(_have_keys)} 档）")
    if _loaded_pod_from != CACHE_POD:
        _save_pod_cache(CACHE_POD, POD_RES, _pod_grid_key)
        print(f"已同步写入主缓存 {CACHE_POD}")
else:
    if POD_RES is None:
        POD_RES = {}
        print("未找到匹配的完整/部分 PoD 缓存，开始全新扫描")
    else:
        miss = len(_expected_keys - _have_keys)
        print(f"从 {_loaded_pod_from} 载入 {len(_have_keys)} 档，断点续跑剩余 {miss} 档")
    _tall = time.time()
    n_total = len(_expected_keys)
    print(f"PoD 临界点扫描：共 {n_total} 个 noise 档，步长 0.25；"
          f"每档同时求 100/10 ppm × PoD50/90；并行 POD_WORKERS={POD_WORKERS}")
    _ckpt_n = 0
    for n_shots in N_SHOTS_LIST:
        jobs = [(k, 7000 + n_shots*1_000_000 + k*20_000)
                for k in range(len(NOISE_GRID[n_shots]))]
        for done, (k, seed) in enumerate(jobs, 1):
            key = (n_shots, float(NOISE_GRID[n_shots][k]))
            if key in POD_RES:
                continue
            key, value = solve_pod_noise(n_shots, k, seed)
            POD_RES[key] = value
            _ckpt_n += 1
            if (_ckpt_n % CHECKPOINT_EVERY) == 0:
                _save_pod_cache(CACHE_POD_CKPT, POD_RES, _pod_grid_key)
            if done == 1 or done % 5 == 0 or done == len(jobs):
                c100 = value.get("critical", {}).get("100", {})
                p90 = c100.get("0.90")
                msg = "无有效交点" if not p90 else (
                    f"E90={p90['boost']*E_PULSE_BASE*1e9:.3g} nJ，"
                    f"验证PoD={p90['pod']:.3f}，peak均值={p90['peak_mean']:.2f}"
                )
                elapsed = time.time() - _tall
                print(f"  [N_shots={n_shots} {done}/{len(jobs)}] "
                      f"noise={key[1]:.2f}：{msg}；累计 {elapsed/60:.1f} min；"
                      f"已完成 {len(POD_RES)}/{n_total}")
        _save_pod_cache(CACHE_POD_CKPT, POD_RES, _pod_grid_key)
    _save_pod_cache(CACHE_POD, POD_RES, _pod_grid_key)
    print(f"逐 noise PoD 扫描总用时 {(time.time()-_tall)/60:.1f} min；已写入 {CACHE_POD}")
    if os.path.exists(CACHE_POD_CKPT):
        try:
            os.remove(CACHE_POD_CKPT)
        except OSError:
            pass
'''


def main():
    nb = json.loads(NB.read_text(encoding="utf-8"))
    changed = []
    for i, cell in enumerate(nb["cells"]):
        if cell.get("cell_type") != "code":
            continue
        src = "".join(cell.get("source", []))
        old = src

        if PARAM_OLD in src:
            src = src.replace(PARAM_OLD, PARAM_NEW)
        if CACHE_OLD in src:
            src = src.replace(CACHE_OLD, CACHE_NEW)
        if PRINT_EXTRA in src and "N_WORKERS=" not in src:
            src = src.replace(PRINT_EXTRA, PRINT_EXTRA_NEW)

        if src.startswith("def _noise_chunk_stats") or src.startswith("def run_noise_scan"):
            src = NOISE_SCAN_CELL

        if POD_TAIL_OLD_START in src:
            head, _ = src.split(POD_TAIL_OLD_START, 1)
            src = head + POD_TAIL_NEW

        # 去掉模块 7 里重复 import（已在模块 0）
        if src.startswith("# ---- PoD 专用子窗"):
            src = src.replace(
                "from concurrent.futures import ThreadPoolExecutor, as_completed\n\n",
                "",
            )

        if src != old:
            cell["source"] = src.splitlines(keepends=True)
            cell["execution_count"] = None
            cell["outputs"] = []
            changed.append(i)

    # 清空全部 code 输出，避免旧结果干扰
    for cell in nb["cells"]:
        if cell.get("cell_type") == "code":
            cell["execution_count"] = None
            cell["outputs"] = []

    NB.write_text(json.dumps(nb, ensure_ascii=False, indent=1), encoding="utf-8")
    print(f"已升级 {NB.name}，修改 cell 索引：{changed}")


if __name__ == "__main__":
    main()
