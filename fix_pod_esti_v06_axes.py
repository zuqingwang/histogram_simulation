# -*- coding: utf-8 -*-
"""修 PoD_esti_v06.ipynb：
1) 模块 7（能量轴图）取消无意义的 AXIS 双循环；
2) 模块 8 collect_critical 真正按 noise/bg 取横轴；
3) 各 noise 轴图的 savefig 按轴名区分，避免互相覆盖。
"""
import ast
import json
import re

SRC = "PoD_esti_v05.ipynb"
DST = "PoD_esti_v06.ipynb"

nb = json.load(open(DST, encoding="utf-8"))
v05 = json.load(open(SRC, encoding="utf-8"))


def src_of(i):
    return "".join(nb["cells"][i].get("source", []))


def set_src(i, text):
    body = text.strip("\n")
    lines = body.split("\n")
    nb["cells"][i]["source"] = [ln + "\n" for ln in lines[:-1]] + [lines[-1]]


def find_v05_code(substr):
    for c in v05["cells"]:
        if c["cell_type"] != "code":
            continue
        s = "".join(c.get("source", []))
        if substr in s:
            return s
    raise RuntimeError(f"v05 cell not found: {substr}")


# ---- 1) 模块 7：恢复 v05，仅更新图例口径（横轴不是 noise，不做双循环）----
mod7 = find_v05_code("模块 7 验证图")
mod7 = mod7.replace(
    'label=f"noise={nt:.2f}, T={r[\'T_map\'][tag]}"',
    'label=(f"bg={float(r[\'noise\']):.2f}/noise≈{_ambient_from_pod_rec(r):.2f}, '
    'T={r[\'T_map\'][tag]}")',
)
mod7 = mod7.replace(
    'label=(f"noise={nt:.2f}: mean={rec[\'peak_mean\']:.2f}, "',
    'label=(f"bg={float(r[\'noise\']):.2f}: mean={rec[\'peak_mean\']:.2f}, "',
)
mod7 = mod7.replace(
    'pod_v05_pod_critical_validation.png',
    'pod_v06_pod_critical_validation.png',
)
helper = (
    "def _ambient_from_pod_rec(r):\n"
    "    '''由 POD 记录的 e_lambda 折合环境标准 noise。'''\n"
    "    el = float(r.get('e_lambda', np.nan))\n"
    "    if np.isfinite(el) and el > 0:\n"
    "        r_det = PDE * R_AMB_BASE * (el / PARAMS['ambient']['E_lambda'])\n"
    "        return float(27.0 * p_bin_equilibrium(r_det)[0])\n"
    "    return float('nan')\n\n"
)
set_src(28, helper + mod7)

# ---- 2) 模块 8：完整重写（从文件读模板，避免嵌套三引号）----
mod8_path = "_pod_v06_mod8.py"
with open(mod8_path, "w", encoding="utf-8") as f:
    f.write(
        "def _pod_x_from_rec(r, n_shots, axis_key):\n"
        "    '''POD_RES 横轴：bg=r[\"noise\"](=noise_mc)；noise=环境标准。'''\n"
        "    bg = float(r['noise'])\n"
        "    if axis_key == 'bg':\n"
        "        return bg\n"
        "    el = float(r.get('e_lambda', np.nan))\n"
        "    if np.isfinite(el) and el > 0:\n"
        "        r_det = PDE * R_AMB_BASE * (el / PARAMS['ambient']['E_lambda'])\n"
        "        return float(27.0 * p_bin_equilibrium(r_det)[0])\n"
        "    return bg / float(n_shots)\n"
        "\n"
        "\n"
        "for _axis_key, _axis_label in AXIS_KINDS:\n"
        "    print(f'\\n===== 轴 = {_axis_key}：{_axis_label} =====')\n"
        "\n"
        "    def equiv_distance(boost, D_ref=D_TARGET, p=PARAMS):\n"
        "        '''把 boost 折算成发射能量和反射率不变时的等效距离。'''\n"
        "        if not np.isfinite(boost) or boost <= 0:\n"
        "            return np.nan\n"
        "        alpha = p['channel']['alpha']\n"
        "        Ds = np.logspace(np.log10(0.3), np.log10(5000.0), 6000)\n"
        "        vals = (D_ref**2 / Ds**2) * np.exp(-2*alpha*(Ds-D_ref))\n"
        "        if boost > vals[0] or boost < vals[-1]:\n"
        "            return np.nan\n"
        "        return float(np.interp(-boost, -vals, Ds))\n"
        "\n"
        "    def collect_critical(n_shots, far_tag, level):\n"
        "        rows = []\n"
        "        for nt in NOISE_GRID[n_shots]:\n"
        "            r = POD_RES.get((n_shots, float(nt)))\n"
        "            if not r or 'critical' not in r:\n"
        "                continue\n"
        "            rec = r['critical'].get(far_tag, {}).get(f'{level:.2f}')\n"
        "            if not rec:\n"
        "                continue\n"
        "            x = _pod_x_from_rec(r, n_shots, _axis_key)\n"
        "            rows.append((\n"
        "                x, rec['boost'], rec['pod'], rec['peak_mean'],\n"
        "                rec['peak_std'], r['T_map'][far_tag], equiv_distance(rec['boost']),\n"
        "            ))\n"
        "        return np.asarray(rows, float)\n"
        "\n"
        "    _ls_by_tag = {\n"
        "        '100ppm': '-', '10ppm': '--',\n"
        "        '5pct': '-.', '1pct': ':',\n"
        "        '0p5pct': (0, (3, 1, 1, 1)), '0p1pct': (0, (1, 1)),\n"
        "    }\n"
        "    _cns = {1: 'tab:blue', 4: 'tab:red'}\n"
        "\n"
        "    fig, ax = plt.subplots(1, 3, figsize=(19, 5.8))\n"
        "    for n_shots in N_SHOTS_LIST:\n"
        "        for far in TARGET_FARS:\n"
        "            tag = FAR_TAG[far]\n"
        "            a = collect_critical(n_shots, tag, 0.90)\n"
        "            if not a.size:\n"
        "                continue\n"
        "            ls = _ls_by_tag.get(tag, '-')\n"
        "            c = _cns[n_shots]\n"
        "            ax[0].semilogy(a[:, 0], a[:, 1]*E_PULSE_BASE*1e9, ls=ls, color=c, lw=1.7,\n"
        "                           label=f'N={n_shots}, {FAR_LABEL[far]}')\n"
        "            ax[1].plot(a[:, 0], a[:, 3], ls=ls, color=c, lw=1.7,\n"
        "                       label=f'N={n_shots}, {FAR_LABEL[far]} peak均值')\n"
        "            ax[1].plot(a[:, 0], a[:, 5], ls=':', color=c, alpha=0.35, lw=1.0)\n"
        "            ax[2].plot(a[:, 0], a[:, 6], ls=ls, color=c, lw=1.7,\n"
        "                       label=f'N={n_shots}, {FAR_LABEL[far]}')\n"
        "\n"
        "    ax[0].set_xlabel(_axis_label); ax[0].set_ylabel('PoD90 临界能量 [nJ]')\n"
        "    ax[0].set_title(f'★ PoD90 临界发射能量 vs {_axis_key}')\n"
        "    ax[0].legend(fontsize=6.5, ncol=2); ax[0].grid(alpha=0.3, which='both')\n"
        "\n"
        "    ax[1].set_xlabel(_axis_label); ax[1].set_ylabel('peak 均值 / T [计数]')\n"
        "    ax[1].set_title('★ PoD90 临界 peak 均值（点线≈对应 T）')\n"
        "    ax[1].legend(fontsize=6.2, ncol=2); ax[1].grid(alpha=0.3)\n"
        "\n"
        "    ax[2].set_xlabel(_axis_label); ax[2].set_ylabel('等效距离 [m]')\n"
        "    ax[2].set_title(f'★ PoD90 等效探测距离 vs {_axis_key}')\n"
        "    ax[2].legend(fontsize=6.5, ncol=2); ax[2].grid(alpha=0.3)\n"
        "\n"
        "    plt.suptitle(\n"
        "        f'模块 8　完整 0.25-noise × FAR{[FAR_LABEL[f] for f in TARGET_FARS]}'\n"
        "        f'（PoD90，滤前，横轴={_axis_key}）',\n"
        "        fontsize=12,\n"
        "    )\n"
        "    plt.tight_layout(rect=[0, 0, 1, 0.93])\n"
        "    plt.savefig(f'pod_v06_summary_dense_{_axis_key}.png', dpi=120, bbox_inches='tight')\n"
        "    plt.show()\n"
        "\n"
        "    print('='*120)\n"
        "    print(f\"{'N':>3}{_axis_key:>8}{'FAR':>10}{'PoD目标':>8}{'PoD验证':>9}{'T':>5}\"\n"
        "          f\"{'peak均值':>10}{'能量[nJ]':>12}{'距离[m]':>10}\")\n"
        "    for n_shots in N_SHOTS_LIST:\n"
        "        stride = max(1, len(NOISE_GRID[n_shots]) // 10)\n"
        "        for nt in NOISE_GRID[n_shots][::stride]:\n"
        "            r = POD_RES.get((n_shots, float(nt)))\n"
        "            if not r or 'critical' not in r:\n"
        "                continue\n"
        "            xval = _pod_x_from_rec(r, n_shots, _axis_key)\n"
        "            for far in TARGET_FARS:\n"
        "                tag = FAR_TAG[far]\n"
        "                for level in POD_LEVELS:\n"
        "                    rec = r['critical'].get(tag, {}).get(f'{level:.2f}')\n"
        "                    if not rec:\n"
        "                        continue\n"
        "                    print(f'{n_shots:>3d}{xval:>8.2f}{FAR_LABEL[far]:>10}'\n"
        "                          f'{level:>8.0%}{rec[\"pod\"]:>9.3f}{r[\"T_map\"][tag]:>5d}'\n"
        "                          f'{rec[\"peak_mean\"]:>10.2f}'\n"
        "                          f'{rec[\"boost\"]*E_PULSE_BASE*1e9:>12.3g}'\n"
        "                          f'{equiv_distance(rec[\"boost\"]):>10.1f}')\n"
    )
set_src(30, open(mod8_path, encoding="utf-8").read())

# ---- 3) cell 20/22/25：savefig 加轴名后缀 ----
def add_axis_to_savefig(s):
    def repl(m):
        path, rest = m.group(1), m.group(2)
        if "{_axis_key}" in path:
            return m.group(0)
        if not path.endswith(".png"):
            return m.group(0)
        base = path[:-4].replace("pod_v05_", "pod_v06_")
        return f'plt.savefig(f"{base}_{{_axis_key}}.png"{rest})'

    s = re.sub(r'plt\.savefig\(\s*"([^"]+)"([^)]*)\)', repl, s)
    s = s.replace('a.set_xlabel("noise")', "a.set_xlabel(_axis_label)")
    s = s.replace("a.set_xlabel('noise')", "a.set_xlabel(_axis_label)")
    s = s.replace('ax[j].set_xlabel("noise")', "ax[j].set_xlabel(_axis_label)")
    s = s.replace('ax.set_xlabel("noise")', "ax.set_xlabel(_axis_label)")
    # 常见中文标签
    for old in [
        'set_xlabel("noise（实测）")',
        'set_xlabel("noise (实测)")',
        'set_xlabel("实测 noise")',
        'set_xlabel("noise")',
    ]:
        s = s.replace(old, "set_xlabel(_axis_label)")
    return s


for i in (20, 22, 25):
    set_src(i, add_axis_to_savefig(src_of(i)))

json.dump(nb, open(DST, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
print("updated", DST)

for i in (19, 20, 22, 25, 28, 30):
    s = src_of(i)
    try:
        ast.parse(s)
        print(f"cell {i}: OK syntax, len={len(s)}")
    except SyntaxError as e:
        print(f"cell {i}: SYNTAX ERROR {e}")
        lines = s.splitlines()
        ln = e.lineno or 1
        for j in range(max(0, ln - 3), min(len(lines), ln + 2)):
            print(f"  {j+1}: {lines[j]}")
