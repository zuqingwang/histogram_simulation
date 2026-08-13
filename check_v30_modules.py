# -*- coding: utf-8 -*-
"""无头校验 PoD_esti_v30.ipynb：按顺序实跑所有 code cell，报告耗时、缺字与异常。

用法：
    python check_v30_modules.py                 # 全部跑
    python check_v30_modules.py --from 16       # 从第 16 个 cell 开始（前面已验证过）
    python check_v30_modules.py --to 20         # 只跑到第 20 个
    python check_v30_modules.py --list          # 只列 cell 地图，不执行

注意：本脚本会真的调用扫描脚本（若缓存缺失），可能非常耗时。
建议先把 run_pod_v30_*.py 跑完再执行。
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
import traceback
import warnings

os.environ.setdefault("MPLBACKEND", "Agg")
sys.stdout.reconfigure(encoding="utf-8")

NB = "PoD_esti_v30.ipynb"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--from", dest="i0", type=int, default=0)
    ap.add_argument("--to", dest="i1", type=int, default=10**9)
    ap.add_argument("--list", action="store_true")
    args = ap.parse_args()

    nb = json.load(open(NB, encoding="utf-8"))
    cells = nb["cells"]

    if args.list:
        for i, c in enumerate(cells):
            src = "".join(c["source"])
            head = next((l.strip() for l in src.splitlines() if l.strip()), "")
            print(f"{i:3d} {c['cell_type'][:4].upper():4s} {head[:80]}")
        return 0

    # 收集 matplotlib 的缺字警告（字体里没有该字形时会 warn）
    glyph_warnings: list[str] = []

    def _showwarning(message, category, filename, lineno, file=None, line=None):
        msg = str(message)
        if "Glyph" in msg or "missing from" in msg or "font" in msg.lower():
            glyph_warnings.append(msg)

    warnings.showwarning = _showwarning

    ns: dict = {"__name__": "__main__"}
    t_all = time.time()
    n_run = 0
    for i, c in enumerate(cells):
        if c["cell_type"] != "code" or not (args.i0 <= i <= args.i1):
            continue
        src = "".join(c["source"])
        head = next((l.strip() for l in src.splitlines() if l.strip()), "")
        print(f"\n{'='*84}\n[cell {i}] {head[:76]}\n{'='*84}", flush=True)
        t0 = time.time()
        try:
            exec(compile(src, f"<{NB} cell {i}>", "exec"), ns)
        except Exception:
            traceback.print_exc()
            print(f"\n[失败] cell {i} 抛异常，已用时 {time.time()-t0:.1f}s")
            return 1
        n_run += 1
        dt = time.time() - t0
        flag = "  ← 慢" if dt > 60 else ""
        print(f"[cell {i} 通过] {dt:.1f}s{flag}", flush=True)

    print(f"\n{'='*84}")
    print(f"[全部通过] 执行 {n_run} 个 code cell，总用时 {(time.time()-t_all)/60:.1f} min")
    if glyph_warnings:
        print(f"\n[字体警告] {len(glyph_warnings)} 条缺字，需修：")
        for m in sorted(set(glyph_warnings))[:30]:
            print("   ", m)
        return 2
    print("[字体检查] 无缺字警告")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
