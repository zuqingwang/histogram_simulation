# -*- coding: utf-8 -*-
"""从 PoD_esti_v05.ipynb 提取纯计算内核，生成 pod_esti_v05_core.py。

生成的模块可被子进程 import（Windows spawn 需要 module-level 可导入），
用于 run_pod_scan_v05.py 的多进程 PoD 扫描。

只提取计算 cell，不含任何绘图 cell：
  cell 2  参数与常数
  cell 4  光链路
  cell 6  SPAD 器件参数
  cell 8  spad_response_g
  cell 9  三个引擎 + noise/E_lambda 互换
  cell 17 纯噪声扫描（含缓存载入，构建 NOISE_RES）
  cell 22 FAR 阈值（构建 THRESH）
  cell 25 PoD 函数部分（含 PoD 缓存读写，截断到主扫描循环之前）
"""
import json
from pathlib import Path

NB = Path("PoD_esti_v05.ipynb")
OUT = Path("pod_esti_v05_core.py")

CALC_CELLS = [2, 4, 6, 8, 9, 17, 22, 25]
# cell 25 保留到 _save_pod_cache 定义结束；`POD_RES = None` 起是主扫描循环，
# 那部分由 run_pod_scan_v05.py 用多进程重新实现。
CELL25_CUT = "POD_RES = None"

HEADER = '''# -*- coding: utf-8 -*-
"""PoD_esti v05 计算内核（由 build_pod_core_v05.py 自动生成，请勿手改）。

来源：PoD_esti_v05.ipynb 的计算 cell，绘图 cell 一律不提取。
用途：供 run_pod_scan_v05.py 的 ProcessPoolExecutor 子进程 import。

import 本模块会：
  1. 重建全部物理参数与引擎；
  2. 载入 pod_esti_v05_cache_noise.npz（或 v04 fallback）构建 NOISE_RES；
  3. 计算 THRESH（六档 FAR 阈值）。
若噪声缓存缺失会触发全量噪声 MC —— 子进程绝不能撞上这种情况，
所以 run_pod_scan_v05.py 会先在主进程确保缓存存在。

环境变量 POD_CORE_QUIET=1 时静音 import 期间的 print。
"""
import builtins as _builtins
import os as _os

_os.environ.setdefault("MPLBACKEND", "Agg")
# 自管进程并行时禁止 BLAS 再开线程，否则 20 进程 × N 线程会互相抢核
for _k in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    _os.environ.setdefault(_k, "1")

_QUIET = _os.environ.get("POD_CORE_QUIET") == "1"
_REAL_PRINT = _builtins.print
if _QUIET:
    _builtins.print = lambda *a, **k: None

'''

FOOTER = '''

# 恢复 print，避免污染宿主进程
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
        if idx == 25:
            cut = src.find(CELL25_CUT)
            if cut < 0:
                raise SystemExit("cell 25 未找到主扫描循环的截断标记")
            src = src[:cut].rstrip() + "\n"
        parts.append(f"\n# {'='*70}\n# ===== 源自 PoD_esti_v05.ipynb cell {idx} =====\n# {'='*70}\n")
        parts.append(src)
        if not src.endswith("\n"):
            parts.append("\n")
    parts.append(FOOTER)

    text = "".join(parts)
    compile(text, str(OUT), "exec")
    OUT.write_text(text, encoding="utf-8")
    print(f"已生成 {OUT}（{len(text):,} 字符，来自 cell {CALC_CELLS}）")


if __name__ == "__main__":
    main()
