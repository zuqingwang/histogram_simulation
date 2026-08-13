# -*- coding: utf-8 -*-
# 一次性迁移：给 v01 的真实光斑缓存补上 f_pix_mode / 新版 grid_key 字段，
# 使其在加入 F_PIX_MODE 开关之后仍能被 run_peak_energy_scan.py 认成有效缓存。
import os
import sys

os.environ["PVE_FPIX_MODE"] = "real"
os.environ.setdefault("POD_CORE_QUIET", "1")
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

import numpy as np

import run_peak_energy_scan as S

SRC = "peak_vs_energy_cache_real.npz"
z = np.load(SRC, allow_pickle=False)
R = {k: z[k] for k in z.files}
z.close()

old_key = str(R["grid_key"])
old_boosts = R["boosts"]
assert old_boosts.size == S.NB and np.allclose(old_boosts, S.BOOSTS), \
    f"网格不一致：缓存 {old_boosts.size} 档 vs 当前 real 网格 {S.NB} 档"

R["grid_key"] = S.GRID_KEY
R["f_pix_mode"] = "real"
R["f_pix"] = S.F_PIX_REAL

S._save(R, SRC)
print(f"旧 grid_key: {old_key}")
print(f"新 grid_key: {S.GRID_KEY}")
print(f"完成档数: {int(R['done'].sum())}/{S.NB}  →  {SRC}")
