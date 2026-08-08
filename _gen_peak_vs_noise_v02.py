# -*- coding: utf-8 -*-
"""从 build_peak_vs_noise_v01.py 生成 v02：口径拆分 noise / bg，凡 noise 轴图都加一张 bg 轴图。

约定（与用户确认一致）
  · noise：环境噪声标准 = 折合到 N_shots=1、宏像元 27 SPAD、每 1 ns bin 的平衡态底计数。
           与发数无关。由 e_lambda → r_det → 27·p_bin 得到。
  · bg   ：当前 N_shots 下，统计窗内累加波形的实测 baseline 均值（= 缓存里的 noise_mc）。
  · N_shots=1 时 bg ≈ noise；N_shots=4 时 bg ≈ 4·noise（平衡态）。
"""
import pathlib
import re
import shutil

ROOT = pathlib.Path(__file__).resolve().parent
SRC = ROOT / "build_peak_vs_noise_v01.py"
DST = ROOT / "build_peak_vs_noise_v02.py"

text = SRC.read_text(encoding="utf-8")

# ---- 版本与文件名 ----
text = text.replace("peak_vs_noise_v01", "peak_vs_noise_v02")
text = text.replace("工作名 `peak_vs_noise`，v01", "工作名 `peak_vs_noise`，v02")
text = text.replace("生成 peak_vs_noise_v02.ipynb", "生成 peak_vs_noise_v02.ipynb")
# 缓存仍可读 v01（数据兼容，只是画图口径变了）
text = text.replace(
    'CACHE = "peak_vs_noise_v02_cache.npz"',
    'CACHE = "peak_vs_noise_v01_cache.npz"  # v02 复用 v01 扫描缓存，只改画图口径',
)
text = text.replace(
    "python peak_vs_noise_scan.py\n```\n\"\"\")",
    "python peak_vs_noise_scan.py   # 产物 peak_vs_noise_v01_cache.npz，v02 直接读\n```\n\"\"\")",
)

# ---- 名词口径 ----
old_gloss = '''- **noise**：宏像元（macro pixel，9×3 = 27 个 SPAD）在 **1 ns 直方图 bin** 上的**平均累加计数**。
  它是本工作的横轴，由环境光谱辐照度 `E_lambda` 经 `r_det_for_noise()` 精确反解得到。'''
new_gloss = '''- **noise（环境标准）**：折合到 **N_shots=1、宏像元 27 SPAD、每 1 ns bin** 的平衡态底计数。
  **与发数无关**。由 `E_lambda → r_det → 27·p_bin_equilibrium` 得到。
- **bg（波形实测 baseline）**：当前 `N_shots` 下，统计窗内累加直方图的实测平均底计数
  （= 缓存字段 `noise_mc`）。**N_shots=1 时 bg ≈ noise；N_shots=4 时 bg ≈ 4·noise**。
- 本版凡横/纵轴原为 noise 的图，都**并排再画一张按 bg 为轴**的图。'''
if old_gloss not in text:
    raise SystemExit("gloss block not found")
text = text.replace(old_gloss, new_gloss)

# ---- 载入 cell：计算 noise / bg ----
old_load = '''DATA = {}
for n in N_LIST:
    DATA[n] = {
        "noise":    np.asarray(_z[f"noise_{n}"], float),    # 目标 noise（横轴）
        "noise_mc": np.asarray(_z[f"noisemc_{n}"], float),  # MC 实测 noise
        "cnt":      np.asarray(_z[f"cnt_{n}"]),   # (n_boost, n_noise, n_tr+2)
        "cnt0":     np.asarray(_z[f"cnt0_{n}"]),  # (n_boost, n_tr+2) —— noise=0
        "done":     np.asarray(_z[f"done_{n}"], bool),
        "T":        np.asarray(_z[f"T_{n}"]),     # (n_far, n_noise) FAR 阈值
        "n_tr":     27 * n,
    }
    if not DATA[n]["done"].all():
        print(f"⚠ N_shots={n}：{(~DATA[n]['done']).sum()} / "
              f"{len(DATA[n]['done'])} 档尚未算完（缓存可能是检查点）")
'''

new_load = '''def _ambient_noise_from_e_lambda(e_lam):
    """环境标准 noise：N_shots=1 宏像元每 1 ns bin 的平衡态底计数。"""
    e_lam = np.asarray(e_lam, float)
    out = np.zeros_like(e_lam, dtype=float)
    e0 = float(core.PARAMS["ambient"]["E_lambda"])
    for i, e in enumerate(e_lam):
        if not np.isfinite(e) or e <= 0:
            out[i] = 0.0
            continue
        r_det = core.PDE * core.R_AMB_BASE * (e / e0)
        out[i] = 27.0 * core.p_bin_equilibrium(r_det)[0]
    return out


DATA = {}
for n in N_LIST:
    noise_target = np.asarray(_z[f"noise_{n}"], float)       # 扫描时的累加目标底（历史字段）
    noise_mc = np.asarray(_z[f"noisemc_{n}"], float)         # 实测累加 baseline = bg
    e_lam = np.asarray(_z[f"elam_{n}"], float)
    noise_amb = _ambient_noise_from_e_lambda(e_lam)          # 环境标准 noise
    DATA[n] = {
        "noise":    noise_amb,          # ★ v02 横轴口径 A：环境标准
        "bg":       noise_mc,           # ★ v02 横轴口径 B：实测 baseline
        "noise_target": noise_target,   # 旧累加目标，仅备查
        "noise_mc": noise_mc,
        "e_lambda": e_lam,
        "cnt":      np.asarray(_z[f"cnt_{n}"]),
        "cnt0":     np.asarray(_z[f"cnt0_{n}"]),
        "done":     np.asarray(_z[f"done_{n}"], bool),
        "T":        np.asarray(_z[f"T_{n}"]),
        "n_tr":     27 * n,
    }
    if not DATA[n]["done"].all():
        print(f"⚠ N_shots={n}：{(~DATA[n]['done']).sum()} / "
              f"{len(DATA[n]['done'])} 档尚未算完（缓存可能是检查点）")
    print(f"N={n}: 环境标准 noise 范围 {noise_amb[0]:.3f}→{noise_amb[-1]:.3f}；"
          f"bg 范围 {noise_mc[0]:.3f}→{noise_mc[-1]:.3f}；"
          f"中位 bg/noise={np.nanmedian(noise_mc/np.maximum(noise_amb,1e-12)):.2f}")
'''

if old_load not in text:
    raise SystemExit("load block not found")
text = text.replace(old_load, new_load)

# ---- 打印里的 N_shots 行 ----
text = text.replace(
    '''    print(f"N_shots={n}: noise {d['noise'][0]:.2f} → {d['noise'][-1]:.2f}"
          f"（{len(d['noise'])} 档），二值硬上限 n_tr = {d['n_tr']}")''',
    '''    print(f"N_shots={n}: noise(环境标准) {d['noise'][0]:.2f} → {d['noise'][-1]:.2f}；"
          f"bg {d['bg'][0]:.2f} → {d['bg'][-1]:.2f}"
          f"（{len(d['noise'])} 档），二值硬上限 n_tr = {d['n_tr']}")''',
)

# ---- 在第一个分析 code 前插入 AXIS_KINDS 与双图画辅助说明 ----
# 策略：把所有 `x = d["noise"]` 换成按 AXIS_KINDS 循环。
# 用较稳妥的方式：在 STAT 计算之后插入 AXIS_KINDS，再把作图模式改为双重循环。

marker = 'STAT[n] = {"mean": m, "std": s, "mean0": m0, "std0": s0}\n'
insert = '''STAT[n] = {"mean": m, "std": s, "mean0": m0, "std0": s0}

# v02：所有原 noise 轴图都要再画一张 bg 轴图
AXIS_KINDS = [
    ("noise", "环境标准 noise（折合 N_shots=1 的底计数 / 1 ns bin）"),
    ("bg",    "实测 baseline bg（当前 N_shots 累加波形统计窗均值 / 1 ns bin）"),
]
'''
if marker not in text:
    raise SystemExit("STAT marker not found")
text = text.replace(marker, insert, 1)

# 将 `x = d["noise"]` 替换为循环结构太脆弱。
# 改用：在每个主要作图 code 块开头，如果只有一次 `x = d["noise"]`，
# 包一层 for _axis_key, _axis_label in AXIS_KINDS。
# 更稳：全局替换 x 取值，并复制整段 fig 逻辑——对 builder 字符串做正则包装。

def wrap_plot_code_blocks(src: str) -> str:
    """找到含 x = d[\"noise\"] 的 code(r\"\"\"...\"\"\") 块，包成 AXIS_KINDS 双循环。"""
    # 拆出 code(r""" ... """) 块
    pattern = re.compile(r'code\(r"""(.*?)"""\)', re.S)
    out = []
    last = 0
    for m in pattern.finditer(src):
        out.append(src[last:m.start()])
        body = m.group(1)
        if 'x = d["noise"]' in body or "x = d['noise']" in body:
            # 避免重复包装
            if "AXIS_KINDS" in body and "for _axis_key" in body:
                out.append(m.group(0))
            else:
                body2 = body.replace('x = d["noise"]', 'x = d[_axis_key]')
                body2 = body2.replace("x = d['noise']", "x = d[_axis_key]")
                # xlabel 里写死的 noise 说明改为用 _axis_label
                body2 = re.sub(
                    r'ax\.set_xlabel\("环境噪声 noise（宏像元每 1 ns bin 的平均累加计数）"\)',
                    'ax.set_xlabel(_axis_label)',
                    body2,
                )
                body2 = re.sub(
                    r'ax\.set_xlabel\("环境噪声 noise（计数 / 1 ns bin）"\)',
                    'ax.set_xlabel(_axis_label)',
                    body2,
                )
                body2 = re.sub(
                    r'ax\.set_xlabel\("noise（计数 / 1 ns bin）"\)',
                    'ax.set_xlabel(_axis_label)',
                    body2,
                )
                # 标题里补充轴口径
                body2 = re.sub(
                    r'(fig\.suptitle\(")',
                    r'\1【"+_axis_key+"】',
                    body2,
                )
                body2 = re.sub(
                    r"(fig\.suptitle\(f\")",
                    r'\1【{_axis_key}】',
                    body2,
                )
                body2 = re.sub(
                    r"(fig\.suptitle\(\nf\")",
                    r'fig.suptitle(\nf"【{_axis_key}】',
                    body2,
                )
                # 子图 title 里的 noise= 对 bg 轴改名
                body2 = body2.replace(
                    'ax.set_title(f"noise = {x[k]:.2f}"',
                    'ax.set_title(f"{_axis_key} = {x[k]:.2f}"',
                )
                wrapped = (
                    'code(r"""\n'
                    'for _axis_key, _axis_label in AXIS_KINDS:\n'
                    '    print(f\"\\n===== 横轴 = {_axis_key}：{_axis_label} =====\")\n'
                    + _indent_block(body2, 4)
                    + '\n""")'
                )
                out.append(wrapped)
        else:
            out.append(m.group(0))
        last = m.end()
    out.append(src[last:])
    return "".join(out)


def _indent_block(block: str, n: int) -> str:
    pad = " " * n
    lines = block.split("\n")
    # 保留首行若为空
    return "\n".join(pad + ln if ln.strip() != "" else ln for ln in lines)


text2 = wrap_plot_code_blocks(text)

# 修复可能的双重 code(r""" 包装问题：wrap 函数已用 code(r""" 包裹
# 但原 text 里的 code(r""" 被替换时，pattern 匹配的是内层——检查 DST

DST.write_text(text2, encoding="utf-8")
print(f"已写入 {DST}")
print("请运行: python build_peak_vs_noise_v02.py")
