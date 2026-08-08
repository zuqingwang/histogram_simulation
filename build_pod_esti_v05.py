# -*- coding: utf-8 -*-
"""以 PoD_esti_v04.ipynb 为基线生成 PoD_esti_v05.ipynb。

本轮只新增模块 5b 的 noise–peak 密度条带图；模块 5 原图保留。
v05 继续使用 v04 的缓存文件名，不执行 notebook。
"""

import json
import uuid
from pathlib import Path


ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "PoD_esti_v04.ipynb"
TARGET = ROOT / "PoD_esti_v05.ipynb"


MD_5B = r"""### 模块 5b —— noise–peak 密度条带

参考 `np_short.ipynb` cell 14 的表现方式：

- 横轴仍是 noise，每个 0.25-noise 档单独成列；
- 纵轴是该次实现的统计窗 peak；
- 每个小方块表示该 noise 档下出现过的一个整数 peak；
- **越宽、颜色越深表示该 peak 在本列中的相对概率越高**；
- 红色实线和深红色虚线分别是由同一批 `peak_cnt` 反解的 100 ppm、10 ppm 阈值。

密度按列归一化，只用于展示每个 noise 档内部的 peak 分布形状，不比较不同 noise 档的绝对样本数。
"""


CODE_5B = r'''# ---- 模块 5b：noise–peak 密度条带（直接复用模块 5 的 NOISE_RES）----
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection


def _threshold_from_cnt_5b(cnt, target_far):
    """由 peak bincount 求满足 P(peak >= T) < target_far 的最小整数 T。"""
    n = int(np.sum(cnt))
    n_ge = np.concatenate([[n], n - np.cumsum(cnt)[:-1]])
    ok = np.where(n_ge < target_far*n)[0]
    return int(ok[0]) if ok.size else int(len(cnt))


NOISE_STEP_5B = 0.25
MAX_HALF_W_5B = 0.38 * NOISE_STEP_5B
MIN_HALF_W_5B = 0.10 * NOISE_STEP_5B
DENSITY_POWER_5B = 0.5
HALF_H_5B = 0.52
CMAP_5B = plt.cm.Blues

fig, axes = plt.subplots(
    1, len(N_SHOTS_LIST), figsize=(8.2*len(N_SHOTS_LIST), 6.8), dpi=140,
    squeeze=False,
)
axes = axes[0]

for ax, n_shots in zip(axes, N_SHOTS_LIST):
    R = NOISE_RES[n_shots]
    rects, colors = [], []

    for k, noise in enumerate(R["noise_mc"]):
        counts = np.asarray(R["peak_cnt"][k], dtype=np.float64)
        peaks = np.flatnonzero(counts)
        if peaks.size == 0:
            continue
        cmax = counts[peaks].max()
        for peak in peaks:
            rel = float((counts[peak] / cmax) ** DENSITY_POWER_5B)
            half_w = MIN_HALF_W_5B + (MAX_HALF_W_5B - MIN_HALF_W_5B) * rel
            rgba = list(CMAP_5B(0.20 + 0.80*rel))
            rgba[3] = 0.25 + 0.75*rel
            rects.append(
                mpatches.Rectangle(
                    (noise-half_w, peak-HALF_H_5B),
                    2*half_w, 2*HALF_H_5B,
                )
            )
            colors.append(rgba)

    ax.add_collection(
        PatchCollection(
            rects, facecolors=colors, edgecolors="none", antialiased=False,
        )
    )

    threshold_100 = np.array([
        _threshold_from_cnt_5b(cnt, 100e-6) for cnt in R["peak_cnt"]
    ])
    threshold_10 = np.array([
        _threshold_from_cnt_5b(cnt, 10e-6) for cnt in R["peak_cnt"]
    ])
    ax.plot(
        R["noise_mc"], threshold_100, color="#e63946", lw=1.9, zorder=5,
        label="100 ppm 阈值",
    )
    ax.plot(
        R["noise_mc"], threshold_10, color="#9b2226", lw=1.9, ls="--",
        zorder=5, label="10 ppm 阈值",
    )
    ax.axhline(
        R["n_tr"], color="k", lw=1.1, ls="-.", alpha=0.65,
        label=f"二值硬上限 = {R['n_tr']}",
    )

    peak_max = max(
        int(np.flatnonzero(np.asarray(cnt))[−1])
        for cnt in R["peak_cnt"] if np.any(cnt)
    )
    ax.set_xlim(
        max(0.0, float(np.min(R["noise_mc"]))-0.3),
        float(np.max(R["noise_mc"]))+0.4,
    )
    ax.set_ylim(0, min(R["n_tr"]+2, peak_max+2))
    ax.set_xlabel("噪声均值 noise [计数 / 1 ns bin]")
    ax.set_ylabel("统计窗 peak [计数 / 1 ns bin]（越宽/越深 = 本列内越密）")
    ax.set_title(
        f"N_shots={n_shots}：noise–peak 密度条带\n"
        f"{len(R['noise_mc'])} 个 noise 档，每档 {N_MC_NOISE:,} 次 MC",
        fontsize=10.5,
    )
    ax.grid(True, alpha=0.25)
    ax.legend(loc="upper left", fontsize=8, framealpha=0.92)

sm = plt.cm.ScalarMappable(cmap=CMAP_5B, norm=plt.Normalize(0, 1))
sm.set_array([])
cbar = fig.colorbar(sm, ax=axes.tolist(), fraction=0.025, pad=0.02)
cbar.set_label("列内相对密度（0–1）")
fig.suptitle(
    f"模块 5b　纯环境光 noise–peak 密度条带（统计窗 {N_STAT} bins，滤前）",
    fontsize=12,
)
fig.subplots_adjust(left=0.07, right=0.93, bottom=0.11, top=0.87, wspace=0.22)
fig.savefig("pod_v05_noise_peak_density_strips.png", dpi=120, bbox_inches="tight")
plt.show()
'''

# JSON 中不能出现 Unicode 负号作为 Python 下标运算符。
CODE_5B = CODE_5B.replace("−1", "-1")


def make_cell(cell_type, source):
    cell = {
        "cell_type": cell_type,
        "id": uuid.uuid4().hex[:8],
        "metadata": {},
        "source": source.splitlines(keepends=True),
    }
    if cell_type == "code":
        cell["execution_count"] = None
        cell["outputs"] = []
    return cell


with SOURCE.open("r", encoding="utf-8") as handle:
    notebook = json.load(handle)

new_cells = []
inserted = False
for cell in notebook["cells"]:
    source = "".join(cell.get("source", []))
    source = source.replace("# PoD_esti v04 ——", "# PoD_esti v05 ——")
    source = source.replace("【PoD_esti v04 汇总】", "【PoD_esti v05 汇总】")
    # 图文件改用 v05 名称；缓存文件 pod_esti_v04_cache_*.npz 故意保持不变。
    source = source.replace("pod_v04_", "pod_v05_")
    cell["source"] = source.splitlines(keepends=True)
    if cell.get("cell_type") == "code":
        cell["execution_count"] = None
        cell["outputs"] = []
    new_cells.append(cell)

    if (
        cell.get("cell_type") == "code"
        and 'plt.savefig("pod_v05_noise_peak.png"' in source
    ):
        new_cells.append(make_cell("markdown", MD_5B))
        new_cells.append(make_cell("code", CODE_5B))
        inserted = True

if not inserted:
    raise RuntimeError("未找到模块 5 原 noise–peak 绘图 cell，未生成 v05")

notebook["cells"] = new_cells
with TARGET.open("w", encoding="utf-8") as handle:
    json.dump(notebook, handle, ensure_ascii=False, indent=1)

print(
    "已生成 PoD_esti_v05.ipynb：保留模块 5 原图并新增 5b 密度条带；"
    "继续使用 v04 缓存；未执行 notebook。"
)
