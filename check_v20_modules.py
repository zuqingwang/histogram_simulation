# -*- coding: utf-8 -*-
"""不用 Jupyter，直接验证 PoD_esti_v20.ipynb 的 code cell。

两种模式：

  --syntax          只对每个 code cell 做 compile()，列出 cell 清单和语法错误（秒级）。
  <cell 编号…>      无头执行指定的 code cell（matplotlib 用 Agg，plt.show 被短路），
                    图会正常落盘、表会正常打印，用来确认没有运行期错误。

执行模式下的上下文由 pod_esti_v20_core 提供（它已经跑过 cell 2/4/6/8/9/17/22/25 的
计算部分，并载入了 NOISE_RES / THRESH）；POD_RES 按 [v20 主缓存, v11 fallback, 检查点]
的顺序另行载入。

注意 cell 之间有依赖（例如模块 13/14 要先跑模块 11 拿到 _COLORS_N、
模块 14 要先跑 cell 28 拿到 equiv_distance、模块 15 要先跑 cell 35 拿到 SIG_M9），
所以按 notebook 顺序给编号，或直接用 --all。

用法：
    python check_v20_modules.py --syntax
    python check_v20_modules.py --all
    python check_v20_modules.py 39 41 49
"""
from __future__ import annotations

import json
import os
import sys

sys.stdout.reconfigure(encoding="utf-8")
os.environ.setdefault("POD_CORE_QUIET", "1")
os.environ.setdefault("MPLBACKEND", "Agg")

NB = "PoD_esti_v20.ipynb"
DEFAULT_CELLS = [28, 31, 33, 35, 37, 39, 41, 43, 45, 47, 49]


def syntax_check(nb) -> int:
    bad = 0
    for i, c in enumerate(nb["cells"]):
        s = "".join(c["source"])
        if c["cell_type"] == "code":
            try:
                compile(s, f"cell{i}", "exec")
            except SyntaxError as e:
                bad += 1
                print(f"!! SYNTAX ERROR cell {i} line {e.lineno}: {e.msg}")
                print("   >>", (e.text or "").rstrip())
        print(f"{i:02d} {c['cell_type'][:4]:<4} {len(s):5d}  "
              f"{s.strip().split(chr(10))[0][:88]}")
    print("syntax errors:", bad)
    return bad


def run_cells(nb, cells) -> None:
    import subprocess

    import numpy as np
    import pod_esti_v20_core as core

    g = dict(vars(core))
    g["__name__"] = "__main__"

    def _run_cmd_stream(cmd):
        env = os.environ.copy()
        env["PYTHONUNBUFFERED"] = "1"
        env["PYTHONIOENCODING"] = "utf-8"
        env.setdefault("POD_CORE_QUIET", "1")
        print(f"[run] {' '.join(str(c) for c in cmd)}", flush=True)
        p = subprocess.Popen([str(c) for c in cmd], stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT, text=True,
                             encoding="utf-8", errors="replace", bufsize=1, env=env)
        for line in p.stdout:
            print(line, end="", flush=True)
        return int(p.wait())

    g["_run_cmd_stream"] = _run_cmd_stream

    pod = None
    for cand in [core.CACHE_POD, *core.CACHE_POD_FALLBACK, core.CACHE_POD_CKPT]:
        pod = core._try_load_pod_cache(cand, core._pod_grid_key)
        if pod is not None:
            print(f"POD_RES ← {cand}（{len(pod)} 档）")
            break
    if pod is None:
        print("!! 未找到 PoD 缓存，模块 8/13/14 会跳过或报错")
    g["POD_RES"] = pod if pod is not None else {}
    print(f"NOISE_RES: {sorted(core.NOISE_RES.keys())}，"
          f"THRESH: {sorted(core.THRESH.keys())}")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    matplotlib.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
    matplotlib.rcParams["axes.unicode_minus"] = False
    plt.show = lambda *a, **k: plt.close("all")
    g["plt"] = plt
    g["np"] = np

    for i in cells:
        c = nb["cells"][i]
        if c["cell_type"] != "code":
            print(f"--- cell {i} 不是 code，跳过")
            continue
        src = "".join(c["source"])
        print("\n" + "=" * 90)
        print(f">>> CELL {i}: {src.strip().split(chr(10))[0][:80]}")
        print("=" * 90, flush=True)
        exec(compile(src, f"{NB}#cell{i}", "exec"), g)

    print("\n全部指定 cell 执行完毕，无异常。")


def main():
    nb = json.load(open(NB, encoding="utf-8"))
    args = sys.argv[1:]
    if args and args[0] == "--syntax":
        raise SystemExit(1 if syntax_check(nb) else 0)
    cells = DEFAULT_CELLS if (not args or args[0] == "--all") else [int(a) for a in args]
    run_cells(nb, cells)


if __name__ == "__main__":
    main()
