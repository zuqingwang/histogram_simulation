# -*- coding: utf-8 -*-
"""从 PoD_esti_v20.ipynb 提取计算内核 → pod_esti_v20_core.py（供多进程 import）。

只提取计算 cell；cell 17 截断到自动开跑之前，避免 import 时触发线程版噪声 MC。
"""
import json
from pathlib import Path

NB = Path("PoD_esti_v20.ipynb")
OUT = Path("pod_esti_v20_core.py")
CALC_CELLS = [2, 4, 6, 8, 9, 17, 22, 25]
CELL17_CUT = "# ---- ★ v20：统一 BG_GRID"
CELL25_CUT = "POD_RES = None"

HEADER = '''# -*- coding: utf-8 -*-
"""PoD_esti v20 计算内核（由 build_pod_core_v20.py 自动生成，请勿手改）。

用途：供 run_pod_v20_noise_scan.py / run_pod_v20_pod_scan.py import。
import 时不自动跑噪声 MC（噪声扫描请用多进程脚本）。
环境变量 POD_CORE_QUIET=1 时静音 print。
"""
import builtins as _builtins
import os as _os

_os.environ.setdefault("MPLBACKEND", "Agg")
for _k in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    _os.environ.setdefault(_k, "1")

_QUIET = _os.environ.get("POD_CORE_QUIET") == "1"
_REAL_PRINT = _builtins.print
if _QUIET:
    _builtins.print = lambda *a, **k: None

'''

FOOTER = '''

def _build_thresh_from_noise(noise_res):
    """由 NOISE_RES 构建 THRESH（六档 FAR）。"""
    thresh = {}
    for n_shots in N_SHOTS_LIST:
        R = noise_res[n_shots]
        ng = len(R["noise_target"])
        rec = {"noise": R["noise_mc"], "sigma_bin": np.zeros(ng)}
        for far in TARGET_FARS:
            tag = FAR_TAG[far]
            rec["T"+tag] = np.zeros(ng, dtype=int)
            rec["far"+tag] = np.zeros(ng)
            rec["nev"+tag] = np.zeros(ng, dtype=int)
            rec["Ti"+tag] = np.zeros(ng, dtype=int)
        for k in range(ng):
            rec["sigma_bin"][k] = np.sqrt(R["n_tr"] * R["p_eq"][k] * (1 - R["p_eq"][k]))
            for far in TARGET_FARS:
                tag = FAR_TAG[far]
                T, f_, nev, _ = far_threshold_from_cnt(R["peak_cnt"][k], far)
                rec["T"+tag][k] = T
                rec["far"+tag][k] = f_
                rec["nev"+tag][k] = nev
                rec["Ti"+tag][k] = far_threshold_binom_indep(
                    R["n_tr"], R["p_eq"][k], N_STAT, far)
        thresh[n_shots] = rec
    return thresh


def _noise_cache_complete(res_all):
    if not res_all:
        return False
    for n in N_SHOTS_LIST:
        if n not in res_all:
            return False
        r = res_all[n]
        if "done" in r:
            if not np.all(r["done"]):
                return False
        elif not all(int(c.sum()) > 0 for c in r["peak_cnt"]):
            return False
        if len(r["noise_target"]) != len(BG_GRID):
            return False
        if not np.allclose(r["noise_target"], BG_GRID, atol=1e-6):
            return False
    return True


_grid_key_noise = np.asarray(BG_GRID, float)
# ★ v20：主缓存名换成 pod_esti_v20_*，但物理内核与网格和 v11 逐字相同，
# 所以按 [主名, fallback…, 检查点] 的顺序找；命中旧版后同步写回 v20 主名（规则三）。
NOISE_RES = None
for _cand in [CACHE_NOISE, *CACHE_NOISE_FALLBACK, CACHE_NOISE_CKPT]:
    NOISE_RES = _try_load_noise_cache(_cand, _grid_key_noise)
    if NOISE_RES is not None:
        if _cand != CACHE_NOISE:
            try:
                _save_noise_cache(CACHE_NOISE, NOISE_RES, _grid_key_noise)
            except Exception:
                pass
        break
if NOISE_RES is None:
    NOISE_RES = {}
THRESH = _build_thresh_from_noise(NOISE_RES) if _noise_cache_complete(NOISE_RES) else {}

_builtins.print = _REAL_PRINT
'''


def main():
    nb = json.loads(NB.read_text(encoding="utf-8"))
    parts = [HEADER]
    for idx in CALC_CELLS:
        cell = nb["cells"][idx]
        if cell["cell_type"] != "code":
            raise SystemExit(f"cell {idx} 不是 code")
        src = "".join(cell["source"])
        if idx == 17:
            cut = src.find(CELL17_CUT)
            if cut < 0:
                raise SystemExit("cell 17 未找到截断标记")
            src = src[:cut].rstrip() + "\n"
        if idx == 25:
            cut = src.find(CELL25_CUT)
            if cut < 0:
                raise SystemExit("cell 25 未找到截断标记")
            src = src[:cut].rstrip() + "\n"
        if idx == 22:
            ti = src.find("\nTHRESH = {}")
            if ti >= 0:
                src = src[:ti].rstrip() + "\n"
        parts.append(f"\n# ===== 源自 PoD_esti_v20.ipynb cell {idx} =====\n")
        parts.append(src)
    parts.append(FOOTER)
    OUT.write_text("".join(parts), encoding="utf-8")
    print(f"wrote {OUT} ({OUT.stat().st_size} bytes)")


if __name__ == "__main__":
    main()
