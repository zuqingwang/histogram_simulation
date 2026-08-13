# -*- coding: utf-8 -*-
"""从 PoD_esti_v30.ipynb 提取计算内核 → pod_esti_v30_core.py（供多进程脚本 import）。

只提取计算 cell。两处必须截断，否则 import 时就会触发线程版 MC：
  · cell 17（噪声扫描）截到「自动开跑」之前
  · cell 22（PoD 扫描）截到 POD_RES 之前
cell 18 只取阈值函数，THRESH 由本文件的 FOOTER 统一构建。
"""
import json
from pathlib import Path

NB = Path("PoD_esti_v30.ipynb")
OUT = Path("pod_esti_v30_core.py")

# 见 show_nb_cells.py 输出：2 参数 / 4 光链路 / 6 时间窗 / 8,9 引擎 /
#                          17 噪声扫描 / 18 阈值 / 22 PoD 扫描
CALC_CELLS = [2, 4, 6, 8, 9, 17, 18, 22]
CUT = {
    17: "# ---- ★ v20：统一 BG_GRID",
    18: "\nTHRESH = {}",
    22: "POD_RES = None",
}

HEADER = '''# -*- coding: utf-8 -*-
"""PoD_esti v30 计算内核（由 build_pod_core_v30.py 自动生成，请勿手改）。

供 run_pod_v30_{noise,pod,sig}_scan.py 与 compare_macro_v30.py import。
import 时不自动跑 MC。环境变量 POD_CORE_QUIET=1 时静音 print。
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
    """由 NOISE_RES 构建 THRESH（★ v30：七档 FAR，含新增的 10%）。"""
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
        if "hist_std" not in r:          # ★ v30：旧结构缓存一律判为不完整
            return False
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
# ★ v30：全量重算，CACHE_NOISE_FALLBACK 已清空；只找主缓存与检查点。
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
            raise SystemExit(f"cell {idx} 不是 code cell")
        src = "".join(cell["source"])
        if idx in CUT:
            cut = src.find(CUT[idx])
            if cut < 0:
                raise SystemExit(f"cell {idx} 未找到截断标记：{CUT[idx]!r}")
            src = src[:cut].rstrip() + "\n"
        parts.append(f"\n# ===== 源自 PoD_esti_v30.ipynb cell {idx} =====\n")
        parts.append(src)
    parts.append(FOOTER)
    OUT.write_text("".join(parts), encoding="utf-8")
    print(f"wrote {OUT} ({OUT.stat().st_size} bytes)")


if __name__ == "__main__":
    main()
