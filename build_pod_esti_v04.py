# -*- coding: utf-8 -*-
"""以 PoD_esti_v02.ipynb 为基线生成 PoD_esti_v04.ipynb。

本轮只修改纯噪声扫描网格：
- N_shots=1：noise 0.25 到 12，步长 0.25；
- N_shots=4：noise 0.25 到 40，步长 0.25。

不执行 notebook，只生成并清空所有代码单元输出。
"""

import json
import runpy
from pathlib import Path


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "PoD_esti_v02.ipynb"
TARGET = ROOT / "PoD_esti_v04.ipynb"


REPLACEMENTS = {
    "# PoD_esti v02 ——": "# PoD_esti v04 ——",
    "## v02 相对 v01（`PoD_esti.ipynb`）改了什么": (
        "## v04 基于 v02：本轮只修改噪声扫描范围"
    ),
    "末尾是 `PoD_esti` 专用参数区，v02 改动的项已标 ★。": (
        "末尾是 `PoD_esti` 专用参数区；v04 仅修改噪声扫描网格，其余参数沿用 v02。"
    ),
    "    1: np.round(np.arange(0.25, 10.0 + 1e-9, 0.25), 4),   # 40 档，上限 27": (
        "    1: np.round(np.arange(0.25, 12.0 + 1e-9, 0.25), 4),   # 48 档，上限 27"
    ),
    "    4: np.round(np.arange(0.50, 20.0 + 1e-9, 0.50), 4),   # 40 档，上限 108": (
        "    4: np.round(np.arange(0.25, 40.0 + 1e-9, 0.25), 4),   # 160 档，上限 108"
    ),
    "（N_shots=1 用 0.25→10 步长 0.25，N_shots=4 用 0.5→20 步长 0.5，各 40 档），": (
        "（N_shots=1 用 0.25→12、N_shots=4 用 0.25→40，步长均为 0.25），"
    ),
    "v02 要跑 40 档 × 2 种 N_shots × 1e6 条，非快速引擎不可。": (
        "v04 要跑 48 + 160 个噪声档、每档 1e6 条，非快速引擎不可。"
    ),
    "40 档 × 1e6 条的 peak 原始样本要 320 MB": (
        "最多 160 档 × 1e6 条的 peak 原始样本体量很大"
    ),
    'CACHE_NOISE = "pod_esti_v02_cache_noise.npz"': (
        'CACHE_NOISE = "pod_esti_v04_cache_noise.npz"'
    ),
    'CACHE_POD   = "pod_esti_v02_cache_pod.npz"': (
        'CACHE_POD   = "pod_esti_v04_cache_pod.npz"'
    ),
    'print(f"【PoD_esti v02 汇总】': 'print(f"【PoD_esti v04 汇总】',
    "pod_v02_engine_check.png": "pod_v04_engine_check.png",
    "pod_v02_noise_waveform.png": "pod_v04_noise_waveform.png",
    "pod_v02_noise_peak.png": "pod_v04_noise_peak.png",
    "pod_v02_threshold.png": "pod_v04_threshold.png",
    "pod_v02_pod_curves.png": "pod_v04_pod_curves.png",
    "pod_v02_summary.png": "pod_v04_summary.png",
}


def transform_source(source):
    text = "".join(source)
    for old, new in REPLACEMENTS.items():
        text = text.replace(old, new)
    return text.splitlines(keepends=True)


with SOURCE.open("r", encoding="utf-8") as handle:
    notebook = json.load(handle)

for cell in notebook["cells"]:
    cell["source"] = transform_source(cell.get("source", []))
    if cell.get("cell_type") == "code":
        cell["execution_count"] = None
        cell["outputs"] = []

with TARGET.open("w", encoding="utf-8") as handle:
    json.dump(notebook, handle, ensure_ascii=False, indent=1)

print(
    "已生成 PoD_esti_v04.ipynb（未执行）："
    "N_shots=1 为 0.25→12，N_shots=4 为 0.25→40，步长均为 0.25。"
)

# 继续写入 v04 的自适应 PoD 扫描与完整 0.25-noise 汇总模块。
# 该脚本只改 notebook 源码，不执行任何 cell。
runpy.run_path(str(ROOT / "enhance_pod_esti_v04.py"), run_name="__main__")
