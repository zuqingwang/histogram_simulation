# -*- coding: utf-8 -*-
"""PoD_esti_v05 → v06：口径拆分 noise / bg，凡 noise 轴图都新增一张 bg 轴图。

约定
  · noise：环境标准 = 折合 N_shots=1、宏像元 27、每 1 ns bin 的平衡态底计数（与发数无关）
  · bg   ：当前 N_shots 累加波形统计窗实测 baseline（= NOISE_RES[*]['noise_mc']）
  · N_shots=1 时 bg ≈ noise

做法：复制 notebook，插入口径说明 + 派生数组；把关键作图 cell 包成
for _axis_key in ('noise','bg') 双画。不重跑 MC，复用 v05 缓存。
"""
import json
import pathlib
import re
import shutil

ROOT = pathlib.Path(__file__).resolve().parent
SRC = ROOT / "PoD_esti_v05.ipynb"
DST = ROOT / "PoD_esti_v06.ipynb"

shutil.copy2(SRC, DST)
nb = json.load(open(DST, encoding="utf-8"))


def src_of(i):
    return "".join(nb["cells"][i].get("source", []))


def set_src(i, text):
    lines = text.split("\n")
    nb["cells"][i]["source"] = [ln + "\n" for ln in lines[:-1]] + ([lines[-1]] if lines else [])
    nb["cells"][i]["outputs"] = []
    nb["cells"][i]["execution_count"] = None


# ---- 找插入点：NOISE_RES 载入完成、THRESH 汇总之后，模块 5 作图之前 ----
insert_at = None
for i, c in enumerate(nb["cells"]):
    s = src_of(i)
    if c["cell_type"] == "code" and "模块 5" in s and "noise_mc" in s and "set_xlabel" in s:
        insert_at = i
        break
if insert_at is None:
    raise SystemExit("找不到模块 5 作图 cell")

MD = r"""## ★ v06 口径：`noise`（环境标准）vs `bg`（波形实测 baseline）

| 符号 | 含义 | 与 N_shots / 宏像元的关系 |
|---|---|---|
| **noise** | 环境噪声标准：折合到 **N_shots=1、宏像元 27 SPAD、每 1 ns bin** 的平衡态底计数 | **无关**（同一 `E_lambda` / `r_det` 下固定） |
| **bg** | 当前配置下，统计窗内累加直方图的**实测 baseline 均值**（= `noise_mc`） | **有关**（N_shots=4 时 bg≈4·noise） |

约定：**N_shots=1 时 bg ≈ noise**。本版凡横/纵轴原为 noise（或 noise 相关）的图，都**再画一张以 bg 为轴**的图。
"""

CODE = r"""# ---- v06：由 NOISE_RES 派生环境标准 noise 与实测 bg ----
AXIS_KINDS = [
    ("noise", "环境标准 noise（折合 N_shots=1 的底计数 / 1 ns bin）"),
    ("bg",    "实测 baseline bg（当前 N_shots 累加波形统计窗均值 / 1 ns bin）"),
]

def _ambient_noise_from_r_det(r_det):
    r_det = np.asarray(r_det, float)
    out = np.zeros_like(r_det)
    for i, r in enumerate(r_det):
        if not np.isfinite(r) or r <= 0:
            out[i] = 0.0
        else:
            out[i] = 27.0 * p_bin_equilibrium(float(r))[0]
    return out

for _ns in N_SHOTS_LIST:
    R = NOISE_RES[_ns]
    R["bg"] = np.asarray(R["noise_mc"], float)
    if "r_det" in R:
        R["noise_ambient"] = _ambient_noise_from_r_det(R["r_det"])
    else:
        # 兜底：用目标累加底 / N_shots（平衡态下与上式一致）
        R["noise_ambient"] = np.asarray(R["noise_target"], float) / float(_ns)
    print(f"N_shots={_ns}: noise(环境标准) "
          f"{R['noise_ambient'][0]:.3f}->{R['noise_ambient'][-1]:.3f}; "
          f"bg {R['bg'][0]:.3f}->{R['bg'][-1]:.3f}; "
          f"中位 bg/noise="
          f"{np.nanmedian(R['bg']/np.maximum(R['noise_ambient'],1e-12)):.2f}")

def _axis_x(R, key):
    '''作图用横轴：noise -> noise_ambient；bg -> bg。'''
    return R["noise_ambient"] if key == "noise" else R["bg"]
"""

def make_cell(cell_type, text, idx):
    if cell_type == "markdown":
        return {"cell_type": "markdown", "id": f"v06_{idx}",
                "metadata": {}, "source": [ln + "\n" for ln in text.strip("\n").split("\n")[:-1]]
                + [text.strip("\n").split("\n")[-1]]}
    return {"cell_type": "code", "id": f"v06_{idx}",
            "execution_count": None, "metadata": {}, "outputs": [],
            "source": [ln + "\n" for ln in text.strip("\n").split("\n")[:-1]]
            + [text.strip("\n").split("\n")[-1]]}


nb["cells"].insert(insert_at, make_cell("code", CODE, "code"))
nb["cells"].insert(insert_at, make_cell("markdown", MD, "md"))

# 重新定位模块 5（插入后 index +2）
plot_cells = []
for i, c in enumerate(nb["cells"]):
    if c["cell_type"] != "code":
        continue
    s = src_of(i)
    # 含 noise 轴标签或用 noise_mc 作横轴的主要作图 cell
    if "set_xlabel" not in s and 'xlabel("noise' not in s and "xlabel('noise" not in s:
        # 也可能 ylabel 含 noise
        if not ("set_ylabel" in s and "noise" in s):
            if "noise_mc" not in s or "plot" not in s:
                continue
    # 收窄：真正画 noise 相关图的模块
    markers = [
        'set_xlabel("噪声均值 noise',
        'set_xlabel("noise")',
        "set_xlabel(\"noise\")",
        'ax.set_xlabel("noise")',
        'xlabel("noise")',
        '噪声均值 noise',
    ]
    if any(m in s for m in markers) or (
        "noise_mc" in s and "set_xlabel" in s and ("peak" in s.lower() or "阈值" in s or "PoD" in s or "klux" in s or "k_th" in s or "距离" in s)
    ):
        plot_cells.append(i)

print("将包装的作图 cell:", plot_cells)


def wrap_dual_axis(src: str) -> str:
    """把 cell 包进 for _axis_key, _axis_label in AXIS_KINDS，并替换横轴数据与标签。"""
    if "for _axis_key, _axis_label in AXIS_KINDS" in src:
        return src
    # 替换常见横轴取值
    body = src
    body = body.replace('x = R["noise_mc"]', 'x = _axis_x(R, _axis_key)')
    body = body.replace("x = R['noise_mc']", "x = _axis_x(R, _axis_key)")
    # plot(R["noise_mc"], ...) 形式
    body = re.sub(
        r'R\["noise_mc"\]',
        '_axis_x(R, _axis_key)',
        body,
    )
    body = re.sub(
        r"R\['noise_mc'\]",
        '_axis_x(R, _axis_key)',
        body,
    )
    # THRESH 里 Tr["noise"] 也是 noise_mc 拷贝——阈值图横轴
    body = re.sub(
        r'Tr\["noise"\]',
        '_axis_x(NOISE_RES[n_shots], _axis_key)',
        body,
    )
    # 模块 8 等用 nt / noise 字段作横轴的：保留 PoD 结果里的 noise 字段作 bg，
    # 同时提供 ambient。这里若出现 xs = ... noise ... 再个案处理。
    # xlabel
    body = body.replace(
        'set_xlabel("噪声均值 noise [计数 / 1 ns bin]（统计窗 152 个 bin 的平均）")',
        'set_xlabel(_axis_label)',
    )
    body = body.replace(
        'set_xlabel("噪声均值 noise [计数 / 1 ns bin]")',
        'set_xlabel(_axis_label)',
    )
    body = body.replace(
        'ax.set_xlabel("噪声均值 noise [计数 / 1 ns bin]")',
        'ax.set_xlabel(_axis_label)',
    )
    body = body.replace(
        'ax[0].set_xlabel("噪声均值 noise [计数 / 1 ns bin]")',
        'ax[0].set_xlabel(_axis_label)',
    )
    body = body.replace(
        'ax[2].set_xlabel("噪声均值 noise [计数 / 1 ns bin]")',
        'ax[2].set_xlabel(_axis_label)',
    )
    body = body.replace(
        'ax[0].set_xlabel("noise")',
        'ax[0].set_xlabel(_axis_key)',
    )
    body = body.replace(
        'ax[1].set_xlabel("noise")',
        'ax[1].set_xlabel(_axis_key)',
    )
    body = body.replace(
        'ax[2].set_xlabel("noise")',
        'ax[2].set_xlabel(_axis_key)',
    )
    # suptitle 加标签
    if "plt.suptitle" in body and "【" not in body:
        body = body.replace("plt.suptitle(", 'plt.suptitle(f"【{_axis_key}】" + ', 1)
        # 上面可能把 f"【{_axis_key}】" + f"..." 搞坏；改更稳的方式
    # 撤回笨重替换，改用前缀 print
    # 缩进整 cell
    lines = body.split("\n")
    indented = "\n".join(("    " + ln if ln.strip() else ln) for ln in lines)
    wrapped = (
        "for _axis_key, _axis_label in AXIS_KINDS:\n"
        "    print(f\"\\n===== 轴 = {_axis_key}：{_axis_label} =====\")\n"
        + indented
    )
    return wrapped


# 修正 wrap：不要破坏已有 f-string suptitle。用简单前缀即可。
def wrap_dual_axis_v2(src: str) -> str:
    if "for _axis_key, _axis_label in AXIS_KINDS" in src:
        return src
    body = src
    body = re.sub(r'R\["noise_mc"\]', '_axis_x(R, _axis_key)', body)
    body = re.sub(r"R\['noise_mc'\]", '_axis_x(R, _axis_key)', body)
    # 注意：THRESH 汇总打印等非作图不要误伤——只在本函数用于已筛选的 plot_cells
    body = body.replace(
        'set_xlabel("噪声均值 noise [计数 / 1 ns bin]（统计窗 152 个 bin 的平均）")',
        'set_xlabel(_axis_label)',
    )
    body = body.replace(
        'a.set_xlabel("噪声均值 noise [计数 / 1 ns bin]")',
        'a.set_xlabel(_axis_label)',
    )
    body = body.replace(
        'ax.set_xlabel("噪声均值 noise [计数 / 1 ns bin]")',
        'ax.set_xlabel(_axis_label)',
    )
    body = body.replace(
        'ax[0].set_xlabel("噪声均值 noise [计数 / 1 ns bin]")',
        'ax[0].set_xlabel(_axis_label)',
    )
    body = body.replace(
        'ax[2].set_xlabel("噪声均值 noise [计数 / 1 ns bin]")',
        'ax[2].set_xlabel(_axis_label)',
    )
    body = body.replace('ax[0].set_xlabel("noise")', 'ax[0].set_xlabel(_axis_key)')
    body = body.replace('ax[1].set_xlabel("noise")', 'ax[1].set_xlabel(_axis_key)')
    body = body.replace('ax[2].set_xlabel("noise")', 'ax[2].set_xlabel(_axis_key)')
    # Tr["noise"] 用于画阈值曲线横轴
    body = re.sub(
        r'(ax\d*\.plot\()\s*Tr\["noise"\]',
        r'\1_axis_x(NOISE_RES[n_shots], _axis_key)',
        body,
    )
    body = re.sub(
        r'plot\(Tr\["noise"\]',
        'plot(_axis_x(NOISE_RES[n_shots], _axis_key)',
        body,
    )
    lines = body.split("\n")
    indented = "\n".join(("    " + ln if ln.strip() else ln) for ln in lines)
    return (
        "for _axis_key, _axis_label in AXIS_KINDS:\n"
        "    print(f\"\\n===== 轴 = {_axis_key}：{_axis_label} =====\")\n"
        + indented
    )


for i in plot_cells:
    old = src_of(i)
    new = wrap_dual_axis_v2(old)
    set_src(i, new)
    print(f"  wrapped cell {i}, len {len(old)} → {len(new)}")

# 模块 8：横轴常为从 POD_RES 取的 noise 字段（即旧 noise_mc）。需要双画。
for i, c in enumerate(nb["cells"]):
    if c["cell_type"] != "code":
        continue
    s = src_of(i)
    if "POD_RES" in s and "set_xlabel" in s and "noise" in s and "for _axis_key" not in s:
        if "PoD90" in s or "临界能量" in s or "等效距离" in s or "peak 均值" in s:
            # 构造 xs_noise / xs_bg：用 e_lambda 反算 ambient，bg 用 r['noise']
            body = s
            # 在作图前注入：若循环 POD_RES，把横轴做成可选
            # 保守策略：整 cell 双循环，并把用作横轴的 noise 列表换成 ambient 或 bg
            # 常见写法：xs.append(r['noise']) 或 nt
            if "for _axis_key" in body:
                continue
            # 注入辅助：根据 e_lambda 算 ambient
            prefix = (
                "def _pod_axis_value(r, key):\n"
                "    bg = float(r.get('noise', np.nan))\n"
                "    if key == 'bg':\n"
                "        return bg\n"
                "    el = float(r.get('e_lambda', np.nan))\n"
                "    if not np.isfinite(el) or el <= 0:\n"
                "        # 无 e_lambda 时用 bg/N_shots 兜底（需外层提供 n_shots）\n"
                "        return bg  # 临时；下面按 N 再除\n"
                "    r_det = PDE * R_AMB_BASE * (el / PARAMS['ambient']['E_lambda'])\n"
                "    return 27.0 * p_bin_equilibrium(r_det)[0]\n\n"
            )
            body = prefix + body
            # 替换 xlabel noise
            body = body.replace('set_xlabel("noise")', 'set_xlabel(_axis_key)')
            body = body.replace("set_xlabel('noise')", "set_xlabel(_axis_key)")
            # 若有 xs.append( something noise ) — 个案：读源
            lines = body.split("\n")
            indented = "\n".join(("    " + ln if ln.strip() else ln) for ln in lines)
            body = (
                "for _axis_key, _axis_label in AXIS_KINDS:\n"
                "    print(f\"\\n===== 轴 = {_axis_key}：{_axis_label} =====\")\n"
                + indented
            )
            set_src(i, body)
            print(f"  wrapped PoD summary cell {i}")

# 更新标题 markdown 第一格加 v06 说明
for i, c in enumerate(nb["cells"]):
    if c["cell_type"] == "markdown" and "PoD_esti" in src_of(i)[:200]:
        s = src_of(i)
        if "v06" not in s[:300]:
            s = s.replace("v05", "v06", 1)
            s = (
                "**v06**：在 v05 基础上拆分 **noise（环境标准）** / **bg（实测 baseline）**，"
                "凡 noise 轴图都加画一张 bg 轴图。缓存仍用 `pod_esti_v05_cache_*.npz`。\n\n"
                + s
            )
            set_src(i, s)
            print(f"  retitled intro cell {i}")
        break

json.dump(nb, open(DST, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
print(f"已写入 {DST}")
