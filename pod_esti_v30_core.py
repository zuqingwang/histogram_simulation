# -*- coding: utf-8 -*-
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


# ===== 源自 PoD_esti_v30.ipynb cell 2 =====
import json, os, time
# ★ 必须在 import numpy/scipy 之前：多线程自管并行时禁止 BLAS 再开线程
for _k in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_k, "1")
from concurrent.futures import ThreadPoolExecutor, as_completed
import threading
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.special import erf
from scipy.stats import binom as _binom, norm as _norm

for _f in ["Microsoft YaHei", "SimHei", "Noto Sans CJK SC", "Source Han Sans SC"]:
    try:
        matplotlib.rcParams["font.sans-serif"] = [_f]; break
    except Exception:
        continue
matplotlib.rcParams["axes.unicode_minus"] = False

C_LIGHT  = 2.99792458e8      # 光速 [m/s]
H_PLANCK = 6.62607015e-34    # 普朗克常数 [J·s]

# ============================ 物理参数（与 v45 / v01 逐项一致，未改） ============================
PARAMS = {
    "laser": {
        "wavelength": 905e-9,   # 波长 [m]
        "amp_mode": "peak",     # "peak"=用 P_peak 定幅 / "energy"=用 E_pulse 定幅
        "P_peak": 235.0,        # 峰值光功率 [W]
        "E_pulse": 800e-9,      # 目标单脉冲能量 [J]（amp_mode="energy" 时生效）
        "tau_r": 0.7e-9,        # 双指数上升时间常数 [s]
        "tau_f": 1.9e-9,        # 双指数下降时间常数 [s]
        "x_L": 60e-6,           # 发光面 X 尺寸（慢轴, 1/e² 全宽）[m]
        "y_L": 1150e-6,         # 发光面 Y 尺寸（快轴, 1/e² 全宽）[m]
    },
    "tx": {"f_TX": 36e-3, "D_TX": 12.7e-3, "eta_coupling": 0.87, "T_TX": 0.92},
    "channel": {"alpha": 0.1e-3},          # 消光系数 @905 nm [1/m]（0.1/km ≈ 晴朗）
    "rx": {"f_RX": 25e-3, "D_RX": 13e-3, "eta_RX": 0.82,
           "T_RX": 0.90, "T_filter": 0.90, "filter_bw": 12e-9},
    "spad": {
        "PDE": 0.30,            # PDE_max，满过电压时的峰值光子探测效率
        "DCR": 0.0e3,           # 单像元暗计数率 [cps]
        "jitter_sigma": 100e-12,# IRF 高斯 σ [s]
        "tau_rc": 8.7315e-9,    # RC 恢复时间常数 τ = R·C_J [s]（使计数死区 = 8 ns）
        "Vov_max": 3.3,         # 满过电压 [V]
        "Vth_frac": 0.60,       # 计数所需最小过电压占比
        "reset_mode": "count",
        "resp_shape": "exp",    # 响应函数 g(Vov) 形状
        "resp_k": 3.0,
    },
    "spad_array": {"pitch": 10e-6, "Nx": 9, "Ny": 120, "fill_factor": 1.00},
    "ambient": {
        "enable": True,
        "E_lambda": 0.68,       # 905 nm 处光谱辐照度 [W/m²/nm]，≈100 klux 白天（基准值）
        "surface_rho": 0.10,    # 被观测面反射率（用于环境光漫反射）
    },
    "hist": {"bin_width": 1e-9, "seed": 0},
}
E_PHOTON = H_PLANCK * C_LIGHT / PARAMS["laser"]["wavelength"]

# ============================ PoD_esti 专用参数区（可调） ============================
# ---- 场景 ----
D_TARGET   = 15.0        # 目标距离 [m] → ToF ≈ 100.07 ns，落在 0–200 ns 窗正中
RHO_TARGET = 0.10        # 目标反射率（与 v45 主回波一致）

# ---- 采集窗与统计窗 ----
WIN_LO_NS  = 0.0         # 采集窗起点 [ns]
WIN_HI_NS  = 200.0       # 采集窗终点 [ns]
TRIM_NS    = 24.0        # 掐头去尾各多少 ns（= 3 × 过阈窗宽 8 ns），统计时剔除
WARM_NS    = 50.0        # 光子生成的左侧护带（暖机）[ns]，采样点不动
DT_FINE    = 200e-12     # 逐光子引擎的细网格步长 [s]（同 v45）

# ---- 宏像元（macro pixel）----
MACRO_BX   = 9           # 沿短边 x 的 SPAD 数
MACRO_BY   = 3           # 沿长边 y 的 SPAD 数   ⇒ 27 个 SPAD
N_SHOTS_MAX  = 4                 # ★ v11：一次最多仿 4 发（按 N 分别扫 bg）
N_SHOTS_LIST = [1, 2, 4]         # ★ v11：N∈{1,2,4}；同 bg 网格对比

# ---- ★ v02：第 1、2 步的噪声档，改为【目标 noise 线性等间距】 ----
#   键 = N_shots；值 = 目标底噪（宏像元每 1 ns bin 的平均计数）。
#   两条网格都覆盖到各自二值硬上限（27 / 108）的约 20%，并都包含 noise = 9 与 10。
# ★ v11：各 N 目标 bg 网格统一（步长 0.25）；仿真时 noise_amb = bg / N
BG_GRID = np.round(np.arange(0.25, 12.0 + 1e-9, 0.25), 4)  # 48 档，统一 bg
NOISE_GRID = {n: BG_GRID.copy() for n in N_SHOTS_LIST}       # 目标 = bg（兼容旧键 noise_target）
NOISE_GRID_AMB = BG_GRID.copy()  # 兼容旧变量名；v11 中表示统一 bg 网格，不再表示「单次 ambient 扫轴」
# ---- ★ v05：噪点率 / FAR 目标（ppm + 百分数）----
# tag 用于字典键与文件字段；label 用于图例显示
FAR_SPECS = [
    (10e-6,  "10ppm",  "10 ppm"),
    (100e-6, "100ppm", "100 ppm"),
    (0.001,  "0p1pct", "0.1%"),
    (0.005,  "0p5pct", "0.5%"),
    (0.01,   "1pct",   "1%"),
    (0.05,   "5pct",   "5%"),
    (0.10,   "10pct",  "10%"),      # ★ v30 新增
]
TARGET_FARS = [v for v, _, _ in FAR_SPECS]
FAR_TAG     = {v: t for v, t, _ in FAR_SPECS}
FAR_LABEL   = {v: lab for v, _, lab in FAR_SPECS}
FAR_TAGS    = [t for _, t, _ in FAR_SPECS]
FAR_TAG_TO_LABEL = {t: lab for _, t, lab in FAR_SPECS}
# ★ v30：阈值七条都算（由 peak 分布直接得到，几乎不花钱），
#        但 PoD 临界能量只对下面这四条求解，机时省一大半。
POD_FARS = [0.005, 0.01, 0.05, 0.10]
POD_FAR_TAGS = [FAR_TAG[f] for f in POD_FARS]
N_MC_SIG    = 20_000            # ★ v30：模块 8/14 固定信号扫描（v20 为 8000）
N_MC_NOISE  = 1_000_000         # 每档纯噪声 MC 条数（10 ppm 需 ≥1e6 才有约 10 个越阈事件）
# 分块大小：快速引擎 A 峰值内存 ≈ chunk × n_tr × 4 B × 约 8 个中间数组。
# chunk=100_000 会到 GB 级并因内存压力显著掉速；25_000 约 150 MB，实测最快。
# ---- ★ v05：并行与缓存策略（CPU 支持 20 线程）----
N_WORKERS   = 20         # 全机线程预算
MC_CHUNK    = 5_000      # 纯噪声分块；20 线程时控制峰值内存
NOISE_WORKERS = N_WORKERS
POD_BIN_WORKERS = N_WORKERS  # ★ 单层：外层 20 档并行；内层不再嵌套线程池
POD_MC_CHUNK = 250       # 信号 MC 再切块，喂满内层线程
CHECKPOINT_EVERY = 8     # ★ 每 8 档落盘一次（原先每档压缩写盘会饿死 CPU）

# ---- 第 3 步：能量扫描与 PoD ----
# ---- v04：完整 noise 网格上的自适应 PoD 交点 ----
NOISE_POD = {n: NOISE_GRID[n].copy() for n in N_SHOTS_LIST}
POD_WORKERS = 1  # ★ 内层串行，避免 4×5 嵌套线程 + GIL 吃不满
N_POD_COARSE = 15        # ★ v30：11→15，粗网格间距从 0.4 缩到 ~0.29 decade
N_MC_POD_COARSE = 600    # ★ v30：300→600，粗交点定位误差减半
N_POD_LOCAL_PER_ROOT = 7 # ★ v30：5→7
N_MC_POD_LOCAL = 800     # 每个局部能量点的 MC 次数
N_MC_POD_VERIFY = 5000   # 每个最终临界点的独立验证次数
POD_LEVELS = [0.50, 0.90]
POD_LOG_BOOST_MIN = -6.0
POD_LOG_BOOST_MAX = 2.0
POD_LOCAL_HALF_DECADE = 0.35   # ★ v30：0.22→0.35，粗交点偏一点也还能罩住真根
POD_VERIFY_TOL = 0.02
POD_VERIFY_ROUNDS = 6          # ★ v30：临界点迭代验证的最大轮数（原实现只有 1 步）
SIG_PRE_NS = 3.0         # 信号窗：ToF 之前 [ns]
SIG_POST_NS= 12.0        # 信号窗：ToF 之后 [ns]（覆盖 8 ns 过阈窗 + 余量）
POD_WARM_NS= 60.0        # PoD 子窗的暖机长度 [ns]（≫ 3τ_RC + T_OVER ≈ 34 ns）

# ---- ★ v05 缓存：主文件 + 兼容读取旧版 + 增量检查点 ----
USE_CACHE = True
CACHE_NOISE = "pod_esti_v30_cache_noise.npz"
CACHE_POD   = "pod_esti_v30_cache_pod.npz"
CACHE_SIG   = "pod_esti_v30_cache_signal.npz"  # 模块 9.3 / 15 固定信号
# ★ v30：FAR 列表与 res 结构都变了（新增 10% 与 hist_std），旧缓存一律不复用。
CACHE_NOISE_FALLBACK = []
CACHE_POD_FALLBACK   = []
CACHE_NOISE_CKPT = "pod_esti_v30_cache_noise.partial.npz"
CACHE_POD_CKPT   = "pod_esti_v30_cache_pod.partial.npz"
CACHE_SIG_CKPT   = "pod_esti_v30_cache_signal.partial.npz"

print(f"单光子能量 E_photon = {E_PHOTON:.3e} J")
print(f"目标 D = {D_TARGET} m → ToF = {2*D_TARGET/C_LIGHT*1e9:.2f} ns")
print(f"采集窗 {WIN_LO_NS:.0f}–{WIN_HI_NS:.0f} ns，bin 宽 {PARAMS['hist']['bin_width']*1e9:.0f} ns，"
      f"掐头去尾各 {TRIM_NS:.0f} ns")
print(f"宏像元 {MACRO_BX}×{MACRO_BY} = {MACRO_BX*MACRO_BY} 个 SPAD；N_SHOTS_MAX={N_SHOTS_MAX}；分析 N={N_SHOTS_LIST}")
print(f"  ★ v20 统一 BG_GRID：{BG_GRID[0]:g}→{BG_GRID[-1]:g}，共 {BG_GRID.size} 档，步长 0.25；noise_amb=bg/N")
for _ns, _g in NOISE_GRID.items():
    print(f"  ★ N_shots={_ns} 目标 bg = {_g[0]:g} → {_g[-1]:g}（noise_amb=bg/{_ns}），"
          f"步长 {_g[1]-_g[0]:g}，共 {_g.size} 档（线性等间距）")
print(f"  ★ FAR 目标：{[FAR_LABEL[f] for f in TARGET_FARS]}，每档 {N_MC_NOISE:,} 次 MC")
print(f"  ★ 并行：N_WORKERS={N_WORKERS}，噪声分块={MC_CHUNK:,}；PoD 外层×内层={POD_BIN_WORKERS}×{POD_WORKERS}，POD_MC_CHUNK={POD_MC_CHUNK}；每 {CHECKPOINT_EVERY} 档增量落盘")
print(f"  ★ 缓存主文件：{CACHE_NOISE} / {CACHE_POD}")
print(f"  ★ v30 缓存与 v11/v20 不互通（FAR 列表与 res 字段都变了），本版全量重算")

def _run_cmd_stream(cmd):
    """运行外部命令并实时打印 stdout（耗时多进程扫描用）。"""
    import subprocess, sys
    env = os.environ.copy()
    env["PYTHONUNBUFFERED"] = "1"
    env["PYTHONIOENCODING"] = "utf-8"
    env.setdefault("POD_CORE_QUIET", "1")
    print(f"[run] {' '.join(str(c) for c in cmd)}", flush=True)
    p = subprocess.Popen(
        [str(c) for c in cmd],
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        text=True, encoding="utf-8", errors="replace",
        bufsize=1, env=env,
    )
    assert p.stdout is not None
    for line in p.stdout:
        print(line, end="", flush=True)
    return int(p.wait())

# ---- ★ v30：全局绘图约定与开关 ----
# 三个 N_shots 在全项目所有图里用同一套颜色，方便跨图对照
_COLORS_N = {1: "tab:blue", 2: "tab:green", 4: "tab:red"}
# 模块 3c 的精确引擎 vs 快速引擎比对很费时，且已在 v20 验证为 bit 级一致。
# 需要重新验证时改成 True。
RUN_ENGINE_CHECK = False
# 模块 3d 的 noise → E_lambda 反解校验单次约 2 分钟，闭合误差已验证 <0.2%。
RUN_INVERSE_CHECK = False

print(f"  ★ v30：FAR 共 {len(FAR_SPECS)} 条；PoD 临界能量只解 "
      f"{[FAR_LABEL[f] for f in POD_FARS]}")
print(f"  ★ v30：信号扫描 {N_MC_SIG:,} MC/档；RUN_ENGINE_CHECK={RUN_ENGINE_CHECK}")
# ===== 源自 PoD_esti_v30.ipynb cell 4 =====
# ---- 激光脉冲 ----
def _pulse_norm(p=PARAMS):
    tr, tf_ = p["laser"]["tau_r"], p["laser"]["tau_f"]
    t_peak = np.log(tf_ / tr) / (1.0 / tr - 1.0 / tf_)
    s_peak = np.exp(-t_peak / tf_) - np.exp(-t_peak / tr)
    area_shape = tf_ - tr
    A = (p["laser"]["E_pulse"] * s_peak / area_shape
         if p["laser"]["amp_mode"] == "energy" else p["laser"]["P_peak"])
    return A, s_peak, t_peak, area_shape

def pulse_temporal(t, p=PARAMS):
    """激光时域光功率 P(t) [W]（双指数）。"""
    tr, tf_ = p["laser"]["tau_r"], p["laser"]["tau_f"]
    A, s_peak, _, _ = _pulse_norm(p)
    tpos = np.clip(t, 0.0, None)
    shape = np.where(t >= 0, np.exp(-tpos / tf_) - np.exp(-tpos / tr), 0.0)
    return A * np.clip(shape, 0.0, None) / s_peak

def pulse_energy(p=PARAMS):
    """单脉冲发射能量 [J] = ∫P(t)dt（双指数解析值 = A/s_peak·(τ_f−τ_r)）。"""
    A, s_peak, _, area_shape = _pulse_norm(p)
    return A / s_peak * area_shape

# ---- 发射光学 / 信道 / 目标 ----
def tx_derived(p=PARAMS):
    x_L, y_L = p["laser"]["x_L"], p["laser"]["y_L"]
    f_TX = p["tx"]["f_TX"]; lam = p["laser"]["wavelength"]; w0 = p["tx"]["D_TX"] / 2.0
    tgx, tgy = x_L / (2 * f_TX), y_L / (2 * f_TX)
    tdf = lam / (np.pi * w0)
    tx_, ty_ = np.hypot(tgx, tdf), np.hypot(tgy, tdf)
    return {"theta_x": tx_, "theta_y": ty_, "zR_x": w0 / tx_, "zR_y": w0 / ty_, "w0": w0,
            "eta_TX": p["tx"]["eta_coupling"] * p["tx"]["T_TX"]}

def atm_transmission(D, p=PARAMS):
    return np.exp(-p["channel"]["alpha"] * D)

def beam_spot_size(D, p=PARAMS):
    tx = tx_derived(p); w0 = tx["w0"]
    return (2 * w0 * np.sqrt(1 + (D / tx["zR_x"])**2),
            2 * w0 * np.sqrt(1 + (D / tx["zR_y"])**2))

def time_of_flight(D):
    return 2.0 * D / C_LIGHT

def lambertian_brdf(rho):
    return rho / np.pi

def echo_range_broadening_sigma(D, tilt_deg, p=PARAMS):
    """目标倾斜带来的几何测距展宽 σ [s]（本文件 tilt=0，返回 0）。"""
    x_D, y_D = beam_spot_size(D, p)
    w = 0.5 * np.hypot(x_D, y_D) / np.sqrt(2)
    dz = w * np.tan(np.deg2rad(tilt_deg))
    return (2.0 * dz / C_LIGHT) / 2.0

# ---- 接收光学 ----
def rx_area(p=PARAMS):
    return np.pi * (p["rx"]["D_RX"] / 2.0)**2

def airy_diameter(p=PARAMS):
    fnum = p["rx"]["f_RX"] / p["rx"]["D_RX"]
    return 2.44 * p["laser"]["wavelength"] * fnum

def rx_image_spot_size(D, p=PARAMS):
    """RX 像面光斑 1/e² 全宽 [m] = 几何成像 ⊕ 衍射艾里斑（平方和）。"""
    x_D, y_D = beam_spot_size(D, p)
    s_airy = airy_diameter(p)
    return (np.hypot((x_D / D) * p["rx"]["f_RX"], s_airy),
            np.hypot((y_D / D) * p["rx"]["f_RX"], s_airy))

def link_factor(echo, p=PARAMS):
    """从发射光功率到单位收集比例下接收光功率的总链路因子（无量纲）。"""
    tx = tx_derived(p); D = echo["D"]
    Omega = rx_area(p) / D**2
    eta_rx_total = p["rx"]["eta_RX"] * p["rx"]["T_RX"] * p["rx"]["T_filter"]
    return (tx["eta_TX"] * echo["frac"] * lambertian_brdf(echo["rho"]) * Omega
            * atm_transmission(D, p)**2 * eta_rx_total)

# ---- 像元收集比例 & 环境光 ----
def gaussian_kernel(sigma, dt, n_sigma=5):
    if sigma <= 0:
        return np.array([1.0 / dt])
    half = max(1, int(np.ceil(n_sigma * sigma / dt)))
    tk = np.arange(-half, half + 1) * dt
    k = np.exp(-0.5 * (tk / sigma)**2)
    return k / (k.sum() * dt)

def pixel_grid(p=PARAMS):
    a = p["spad_array"]; pitch = a["pitch"]
    xi = (np.arange(a["Nx"]) - (a["Nx"] - 1) / 2.0) * pitch
    yj = (np.arange(a["Ny"]) - (a["Ny"] - 1) / 2.0) * pitch
    return xi, yj

def pixel_collection_matrix(D, p=PARAMS):
    """每个像元在椭圆高斯像斑上的空间收集比例 f_pix[i,j]（∑ ≤ 1，其余漏到阵列外）。"""
    sx, sy = rx_image_spot_size(D, p)
    sig_x, sig_y = sx / 4.0, sy / 4.0
    xi, yj = pixel_grid(p); pitch = p["spad_array"]["pitch"]
    def _frac(centers_, sig):
        lo = (centers_ - pitch / 2.0) / (np.sqrt(2) * sig)
        hi = (centers_ + pitch / 2.0) / (np.sqrt(2) * sig)
        return 0.5 * (erf(hi) - erf(lo))
    fx = _frac(xi, sig_x); fy = _frac(yj, sig_y)
    return np.outer(fx, fy), fx, fy

def ambient_photon_rate_per_pixel(p=PARAMS, e_lambda=None):
    """单像元环境光【光子到达率】[ph/s]，不含 PDE。
    e_lambda 不给则用 PARAMS 里的基准值 0.68 W/m²/nm（≈100 klux）。"""
    if not p["ambient"]["enable"]:
        return 0.0
    E_l = p["ambient"]["E_lambda"] if e_lambda is None else e_lambda
    E = E_l * (p["rx"]["filter_bw"] * 1e9)                        # 带内辐照 [W/m²]
    L = p["ambient"]["surface_rho"] * E / np.pi                   # 辐亮度 [W/m²/sr]
    iFOV = p["spad_array"]["pitch"] / p["rx"]["f_RX"]
    P_amb = L * iFOV**2 * rx_area(p)
    return P_amb / E_PHOTON * p["rx"]["T_RX"] * p["rx"]["T_filter"]

def signal_photon_rate_fine(echo, f_pix_ij, tf_grid, p=PARAMS):
    """单 SPAD 信号【光子到达率】[ph/s]（不含 PDE），在精细网格 tf_grid 上。"""
    t0 = time_of_flight(echo["D"])
    r = (pulse_temporal(tf_grid - t0, p) * link_factor(echo, p) / E_PHOTON * f_pix_ij)
    sig_b = echo_range_broadening_sigma(echo["D"], echo.get("tilt_deg", 0.0), p)
    if sig_b > 0:
        dtf = tf_grid[1] - tf_grid[0]
        r = np.convolve(r, gaussian_kernel(sig_b, dtf), mode="same") * dtf
    return r

# ---- 自检 ----
ECHO0 = {"D": D_TARGET, "rho": RHO_TARGET, "frac": 1.00, "tilt_deg": 0.0}
E_PULSE_BASE = pulse_energy()
_sx, _sy = rx_image_spot_size(D_TARGET)
print("="*78)
print(f"单脉冲发射能量 E_pulse = {E_PULSE_BASE*1e9:.1f} nJ（P_peak={PARAMS['laser']['P_peak']:.0f} W，"
      f"双指数 τ_r={PARAMS['laser']['tau_r']*1e9:.1f} / τ_f={PARAMS['laser']['tau_f']*1e9:.1f} ns）")
print(f"D={D_TARGET} m: 像面光斑 1/e² 全宽 x={_sx*1e6:.2f} µm, y={_sy*1e6:.2f} µm")
print(f"链路因子 link_factor = {link_factor(ECHO0):.3e}（ρ={RHO_TARGET}）")
print(f"环境光基准 E_lambda = {PARAMS['ambient']['E_lambda']} W/m²/nm（≈100 klux）："
      f"单像元 r_amb = {ambient_photon_rate_per_pixel():.3e} ph/s  "
      f"→ 探测率 r_det = r_amb·PDE = {ambient_photon_rate_per_pixel()*PARAMS['spad']['PDE']:.3e} cps")
# ===== 源自 PoD_esti_v30.ipynb cell 6 =====
# ---- SPAD 器件参数（从 PARAMS["spad"] 解出全局量，命名沿用 v45） ----
_sp = PARAMS["spad"]
PDE        = _sp["PDE"]
JIT        = _sp["jitter_sigma"]
TAU_RC     = _sp["tau_rc"]
VTH_FRAC   = _sp["Vth_frac"]
RESP_SHAPE = _sp["resp_shape"]
RESP_K     = _sp["resp_k"]
jit = JIT

T_OVER = -TAU_RC * np.log(1.0 - VTH_FRAC)          # 过阈窗宽 ≈ 8.00 ns（导出量）

# ---- 时间窗与 bin ----
BIN_W  = PARAMS["hist"]["bin_width"]
WIN_LO, WIN_HI = WIN_LO_NS * 1e-9, WIN_HI_NS * 1e-9
NBINS  = int(round((WIN_HI - WIN_LO) / BIN_W))
CENTERS = WIN_LO + (np.arange(NBINS) + 0.5) * BIN_W        # bin 中心 [s]
TC_NS   = CENTERS * 1e9

# 统计窗（掐头去尾）：bin 中心落在 [TRIM, WIN_HI-TRIM] 内
_keep = (TC_NS >= TRIM_NS) & (TC_NS <= WIN_HI_NS - TRIM_NS)
IDX_STAT = np.where(_keep)[0]
I_STAT0, I_STAT1 = int(IDX_STAT[0]), int(IDX_STAT[-1] + 1)   # 切片 [I_STAT0, I_STAT1)
N_STAT = I_STAT1 - I_STAT0

# 光子生成网格（左扩护带，采样点不动）
TF_GEN = np.arange(WIN_LO - WARM_NS*1e-9, WIN_HI, DT_FINE)

# ---- 宏像元 ----
N_PIX_MACRO = MACRO_BX * MACRO_BY
FPIX, FX, FY = pixel_collection_matrix(D_TARGET)
_n_macro = PARAMS["spad_array"]["Ny"] // MACRO_BY
_macro_fsum = np.array([FPIX[:, m*MACRO_BY:(m+1)*MACRO_BY].sum() for m in range(_n_macro)])
M_PEAK = int(_macro_fsum.argmax())
F_VALS = FPIX[:, M_PEAK*MACRO_BY:(M_PEAK+1)*MACRO_BY].ravel()   # 峰值宏像元 27 个 SPAD 的收集比例

# ---- 环境光基准速率 ----
R_AMB_BASE = ambient_photon_rate_per_pixel()       # E_lambda = 0.68 时的单像元光子率 [ph/s]

# ---- 信号率形状（单位收集比例、boost=1）----
T0_SIG   = time_of_flight(D_TARGET)
T0_SIG_NS= T0_SIG * 1e9
R_SIG_UNIT_GEN = signal_photon_rate_fine(ECHO0, 1.0, TF_GEN)    # [ph/s]，f_pix=1

print("="*78)
print(f"过阈窗宽 T_OVER = -τ_RC·ln(1-Vth_frac) = {T_OVER*1e9:.3f} ns")
print(f"采集窗 {WIN_LO_NS:.0f}–{WIN_HI_NS:.0f} ns → {NBINS} 个 1 ns bin（中心 "
      f"{TC_NS[0]:.1f} … {TC_NS[-1]:.1f} ns）")
print(f"统计窗（掐头去尾 {TRIM_NS:.0f} ns）→ bin 下标 [{I_STAT0}, {I_STAT1}) 共 {N_STAT} 个，"
      f"中心 {TC_NS[I_STAT0]:.1f} … {TC_NS[I_STAT1-1]:.1f} ns")
print(f"生成网格 TF_GEN：{TF_GEN[0]*1e9:.1f} … {TF_GEN[-1]*1e9:.1f} ns，"
      f"步长 {DT_FINE*1e12:.0f} ps，共 {TF_GEN.size} 点（左扩护带 {WARM_NS:.0f} ns）")
print(f"宏像元 {MACRO_BX}×{MACRO_BY} = {N_PIX_MACRO} SPAD；峰值宏像元 m={M_PEAK}"
      f"（Σf_pix={_macro_fsum[M_PEAK]:.4f}）")
for ns_ in N_SHOTS_LIST:
    print(f"  N_shots={ns_} → 轨迹数 n_tr = 27×{ns_} = {N_PIX_MACRO*ns_}，"
          f"二值硬上限 macro_cap = {N_PIX_MACRO*ns_}")
print(f"信号峰位 ToF = {T0_SIG_NS:.2f} ns（bin 下标 {int(T0_SIG_NS)}）")
print(f"环境光基准 r_amb = {R_AMB_BASE:.3e} ph/s，r_det = {R_AMB_BASE*PDE:.3e} cps")
# ===== 源自 PoD_esti_v30.ipynb cell 8 =====
def spad_response_g(vov_frac, shape="linear", k=3.0):
    """SPAD 对光子的响应函数 g(vov_frac) ∈ [0,1]；触发概率 = PDE_max·g。
    约束 g(0)=0（Vov=0 完全不响应）、g(1)=1（满 Vov 触发概率 = PDE_max）。"""
    x = np.clip(vov_frac, 0.0, 1.0)
    if shape == "exp":
        return (1.0 - np.exp(-k * x)) / (1.0 - np.exp(-k))
    return x


def spad_binary_trace(r_sig_fine, r_amb_ph, tf, centers, PDE_max, tau_rc, Vth_frac,
                      jitter_sigma, rng, t_over, t_laser=0.0,
                      resp_shape="linear", resp_k=3.0):
    """【精确引擎，v45 cell 32 原样移植（去掉归因分支）】
    单 SPAD、单次 shot 的二值采样，返回长度 len(centers) 的 int8 数组（每 bin 0/1）。"""
    dt = tf[1] - tf[0]
    nbn = len(centers)
    out = np.zeros(nbn, dtype=np.int8)
    inv_tau = 1.0 / tau_rc
    mu = (r_sig_fine + r_amb_ph) * dt
    n_ph = rng.poisson(mu)
    if n_ph.sum() == 0:
        return out
    t_arr = np.repeat(tf, n_ph)                 # 逐光子到达时刻（升序）
    u = rng.random(t_arr.size)
    last = -1e30
    av = []
    for k in range(t_arr.size):
        t = t_arr[k]
        d = (t - last) * inv_tau
        vov_frac = 1.0 - np.exp(-d) if d < 700 else 1.0
        if u[k] < PDE_max * spad_response_g(vov_frac, resp_shape, resp_k):
            av.append(t); last = t             # 二值：雪崩即复位
    if not av:
        return out
    av = np.asarray(av)
    if jitter_sigma > 0:
        av = av + rng.normal(0.0, jitter_sigma, av.size)   # IRF 抖动
    for tt in av:
        lo = tt + t_laser
        out[(centers >= lo) & (centers < lo + t_over)] = 1
    return out


def macro_hist_exact(n_real, f_arr, r_sig_unit, r_amb, tf, centers, rng, boost=1.0):
    """用精确引擎堆出 (n_real, nbins) 宏像元二值累加直方图（f_arr 已含 N_shots 平铺）。"""
    nb = len(centers)
    h = np.zeros((n_real, nb))
    for i in range(n_real):
        acc = np.zeros(nb, dtype=np.int32)
        for fij in f_arr:
            acc += spad_binary_trace(r_sig_unit * fij * boost, r_amb, tf, centers,
                                     PDE, TAU_RC, VTH_FRAC, JIT, rng, T_OVER, 0.0,
                                     RESP_SHAPE, RESP_K)
        h[i] = acc
    return h
# ===== 源自 PoD_esti_v30.ipynb cell 9 =====
# ============================================================================
# 更新过程工具：累积强度 H、生存函数 S、平衡态点亮概率 p_bin、以及 noise -> r_det 反解
# ============================================================================
U_C    = 20.0 * TAU_RC     # 超过 20τ 后 g≈1，hazard 恒为 r_det ⇒ 尾部可解析积分
E_MAX  = 30.0              # -ln(U) 的截断：P(E>30)=e^-30≈9e-14，远低于任何 MC 规模
N_ETAB = 262144            # 反函数直查表的点数（均匀 E 网格）


def build_renewal_table(r_det, tau_rc=TAU_RC, resp_shape=RESP_SHAPE, resp_k=RESP_K,
                        d_max=None, n_grid=400001):
    """累积强度 H(Δ) = ∫₀^Δ r_det·g(1−e^{−s/τ}) ds，在 [0, d_max] 上数值积分。"""
    if d_max is None:
        d_max = E_MAX / r_det + 40 * tau_rc
    d = np.linspace(0.0, d_max, n_grid)
    h = r_det * spad_response_g(1.0 - np.exp(-d / tau_rc), resp_shape, resp_k)
    H = np.concatenate([[0.0], np.cumsum(0.5 * (h[1:] + h[:-1]) * np.diff(d))])
    return d, H


def p_bin_equilibrium(r_det, t_over=T_OVER, tau_rc=TAU_RC,
                      resp_shape=RESP_SHAPE, resp_k=RESP_K, n=200001):
    """更新过程【平衡态】下，某个 1 ns 采样点被点亮的概率 p_bin，以及平均雪崩间隔 mu。

        点亮 ⟺ 回溯时间 B < T_OVER；平衡态密度 f_B(u) = S(u)/mu
        ⇒ p_bin = (1/mu)·∫₀^{T_OVER} S(u) du,   S(u)=exp(−H(u)),  mu=∫₀^∞ S

    积分网格只覆盖 [0, U_C=20τ]（步长 ≈0.9 ps，足以分辨 8 ns 的过阈窗）；
    尾部 u > U_C 时 hazard 恒为 r_det ⇒ ∫_{U_C}^∞ S = S(U_C)/r_det，解析补上。
    ⚠️ v01 曾把网格设成 d_max=40/r_det 且点数固定，低噪声档步长会粗到上百 ns，
       连 T_OVER 都分辨不了 —— 这里固定 [0, U_C] 就是为了避开这个坑。
    """
    d = np.linspace(0.0, U_C, n)
    h = r_det * spad_response_g(1.0 - np.exp(-d / tau_rc), resp_shape, resp_k)
    H = np.concatenate([[0.0], np.cumsum(0.5 * (h[1:] + h[:-1]) * np.diff(d))])
    S = np.exp(-H)
    mu = np.trapezoid(S, d) + S[-1] / r_det
    m = d <= t_over
    return float(np.trapezoid(S[m], d[m]) / mu), float(mu)


def r_det_for_noise(noise_target, n_tr, lo=1e2, hi=1e11, iters=45):
    """二分反解：使平衡态底噪 = noise_target（宏像元每 bin 平均计数）的 r_det [cps]。"""
    p_t = noise_target / n_tr
    if not (0.0 < p_t < 1.0):
        return np.nan
    for _ in range(iters):
        mid = np.sqrt(lo * hi)
        if p_bin_equilibrium(mid, n=20001)[0] < p_t:
            lo = mid
        else:
            hi = mid
    return float(np.sqrt(lo * hi))


def e_lambda_for_r_det(r_det, p=PARAMS):
    """由所需的单 SPAD 环境探测率反推 E_lambda [W/m²/nm]（链路对 E_lambda 线性）。"""
    return p["ambient"]["E_lambda"] * (r_det / PDE) / R_AMB_BASE


# ============================================================================
# 快速引擎 A（★ v02 优化版）：更新过程反函数【O(1) 直查表】+ float32
# ============================================================================
def build_inv_table(r_det, tau_rc=TAU_RC, resp_shape=RESP_SHAPE, resp_k=RESP_K,
                    e_max=E_MAX, n_etab=N_ETAB):
    """把 Δ = H⁻¹(E) 重采样到【均匀 E 网格】，使采样变成 O(1) 直查 + 线性插值。
    只在建表时做一次二分（np.interp），之后每次采样都是常数时间。"""
    d, H = build_renewal_table(r_det, tau_rc, resp_shape, resp_k)
    e_grid = np.linspace(0.0, e_max, n_etab)
    inv = np.interp(e_grid, H, d).astype(np.float32)
    return inv, np.float32((n_etab - 1) / e_max)


def noise_macro_hist_fast(n_real, n_tr, r_det, rng, win_lo=None, win_hi=None,
                          nbins=None, bin_w=None, t_over=T_OVER, jitter=JIT,
                          warm=None, max_round=20000, inv_tab=None):
    """纯环境光下的 (n_real, nbins) 宏像元二值累加直方图（int32）。

    n_tr = 宏像元 SPAD 数 × N_shots。纯噪声时各 SPAD、各 shot 独立同分布，
    可以直接折进"轨迹数"这一个维度。

    覆盖技巧：bin 被点亮 ⟺ 它之前最近一次雪崩距它 < T_OVER。
    于是把每次雪崩的窗口【在下一次雪崩处截断】，得到互不重叠的区间，
    就能用「差分数组 + bincount + cumsum」直接累进 (实现数, bin 数) 矩阵，无需逐轨迹展开。
    """
    win_lo = WIN_LO if win_lo is None else win_lo
    win_hi = WIN_HI if win_hi is None else win_hi
    nbins  = NBINS  if nbins  is None else nbins
    bin_w  = BIN_W  if bin_w  is None else bin_w
    warm   = WARM_NS*1e-9 if warm is None else warm
    if inv_tab is None:
        inv_tab = build_inv_table(r_det)
    inv, scale = inv_tab
    n_e = inv.size

    t_start = np.float32(win_lo - warm)
    wl = np.float32(win_lo); wh = np.float32(win_hi)
    tov = np.float32(t_over); bw = np.float32(bin_w)
    D_FLOOR = np.float32(1e-13)        # 等待时间下限（0.1 ps），兜底防止时间不前进而空转

    N = n_real * n_tr
    ri = np.repeat(np.arange(n_real, dtype=np.int32), n_tr)
    # 首个间隔：起始时刻完全恢复 (vov=1) ⇒ Exp(r_det)，与精确引擎 last=-1e30 的初值对应
    t_j = (t_start + rng.exponential(1.0 / r_det, N)).astype(np.float32)
    j_j = (rng.normal(0.0, jitter, N).astype(np.float32) if jitter > 0
           else np.zeros(N, np.float32))

    lo_all, hi_all, ri_all = [], [], []
    for _ in range(max_round):
        keep = t_j < wh
        if not keep.any():
            break
        t_j = t_j[keep]; j_j = j_j[keep]; ri = ri[keep]
        m_ = t_j.size
        # Δ = H⁻¹(E)，E ~ Exp(1)。用 standard_exponential 而不是 -log(random())：
        # 后者在 float32 下有 2^-24 概率取到精确 0 → -log(0)=inf → NaN 污染 + 死循环。
        E = rng.standard_exponential(m_, dtype=np.float32)
        x = np.minimum(E, np.float32(E_MAX)) * scale
        i0 = x.astype(np.int32)
        np.clip(i0, 0, n_e - 2, out=i0)
        fr = x - i0
        delta = np.maximum(inv[i0] * (1.0 - fr) + inv[i0 + 1] * fr, D_FLOOR)
        t_n = t_j + delta
        j_n = (rng.normal(0.0, jitter, m_).astype(np.float32) if jitter > 0
               else np.zeros(m_, np.float32))

        lo_all.append(t_j + j_j)
        hi_all.append(np.minimum(t_j + j_j + tov, t_n + j_n))   # 被下一次雪崩截断
        ri_all.append(ri.copy())
        t_j, j_j = t_n, j_n

    diff = np.zeros(n_real * (nbins + 1), dtype=np.int32)
    if lo_all:
        lo_t = np.concatenate(lo_all); hi_t = np.concatenate(hi_all)
        rr = np.concatenate(ri_all).astype(np.int64)
        b_lo = np.clip(np.ceil((lo_t - wl) / bw - 0.5), 0, nbins).astype(np.int64)
        b_hi = np.clip(np.ceil((hi_t - wl) / bw - 0.5), 0, nbins).astype(np.int64)
        m = b_hi > b_lo
        if m.any():
            base = rr[m] * (nbins + 1)
            diff = (np.bincount(base + b_lo[m], minlength=diff.size)
                    - np.bincount(base + b_hi[m], minlength=diff.size)).astype(np.int32)
    return np.cumsum(diff.reshape(n_real, nbins + 1), axis=1)[:, :nbins]


# ============================================================================
# 快速引擎 B：同步时间步进（信号 + 环境），在「实现 × SPAD」维度向量化
# ============================================================================
def _cover_commit(c, tcov, cend, t_over):
    """判断 bin 中心 c 是否落在过阈覆盖区内，并就地推进覆盖终点 `cend`。

    ★ v30 修（2026-08-10）：原实现只保留「最近一次雪崩时刻」`tcov`，判据写成
      `0 <= c - tcov < t_over`。这在**饱和**时是错的 ——
      SPAD 在同一个 bin 内可能连续雪崩多次，`tcov` 被改写成带 IRF 抖动的**更晚**时刻，
      于是 `c - tcov < 0`，本该被前一次雪崩的 8 ns 窗口点亮的 bin 被判成灭的。
      后果：深饱和时平均波形压不出削顶平台（实测 27 个 SPAD 只到 19，精确引擎给 27）。

    正确判据（与 `spad_binary_trace` 取全部窗口并集完全等价）：

        c 被点亮 ⟺ max{ ta + t_over : 雪崩时刻 ta <= c } > c

    因为 bin 是按时间顺序出的，只要把「确实早于 c 的雪崩」并进单调不减的 `cend`，
    再比一次 `cend > c` 即可。`tcov > c`（抖动把它甩到 c 之后）的雪崩本轮不提交，
    留到下一个 bin 再算。
    """
    upd = (tcov <= c) & (tcov + t_over > cend)
    np.copyto(cend, tcov + t_over, where=upd)
    return cend > c


def binary_macro_stepping(n_real, f_arr, r_sig_unit, tgrid, r_amb, centers, rng,
                          boost=1.0, tau_rc=TAU_RC, t_over=T_OVER, pde=PDE,
                          jitter=JIT, resp_shape=RESP_SHAPE, resp_k=RESP_K):
    """返回 (n_real, len(centers)) 宏像元二值累加直方图（int32）。f_arr 已含 N_shots 平铺。

    与 spad_binary_trace 逐步等价：一个细网格步内到达 n~Poisson(μ) 个光子、各自以 φ 触发，
    故该步「至少触发一次」的概率 = 1 − e^{−μ·φ}，且步内至多一次雪崩
    （首个触发后 Vov=0 ⇒ g(0)=0 ⇒ 同步内后续光子不可能再触发）。
    因 t − t_last 恒为步长整数倍，φ 可预先做成查表 phi[age]。

    覆盖判据见 `_cover_commit` 的说明（v30 修：原来只用「最近一次雪崩」会在饱和时丢覆盖）。
    """
    dt = tgrid[1] - tgrid[0]
    n_tr = f_arr.size
    nb = len(centers)
    k_max = int(np.ceil(20.0 * tau_rc / dt))
    phi = pde * spad_response_g(1.0 - np.exp(-np.arange(k_max + 1) * dt / tau_rc),
                                resp_shape, resp_k)
    age  = np.full((n_real, n_tr), k_max, dtype=np.int32)
    tcov = np.full((n_real, n_tr), -1e30)
    cend = np.full((n_real, n_tr), -1e30)
    hist = np.zeros((n_real, nb), dtype=np.int32)
    mu_all = (r_sig_unit[:, None] * f_arr[None, :] * boost + r_amb) * dt
    ib = 0
    for i in range(tgrid.size):
        t = tgrid[i]
        while ib < nb and centers[ib] < t:      # 先出 bin，再处理本步雪崩
            hist[:, ib] = _cover_commit(centers[ib], tcov, cend, t_over).sum(axis=1)
            ib += 1
        p = -np.expm1(-mu_all[i][None, :] * phi[age])
        fire = rng.random((n_real, n_tr)) < p
        age = np.minimum(age + 1, k_max)
        if fire.any():
            age[fire] = 1
            nf = int(fire.sum())
            tcov[fire] = t + (rng.normal(0.0, jitter, nf) if jitter > 0 else 0.0)
    while ib < nb:
        hist[:, ib] = _cover_commit(centers[ib], tcov, cend, t_over).sum(axis=1)
        ib += 1
    return hist




# ============================================================================
# ★ v10：per-shot hist_i 与前缀和
# ============================================================================
def noise_hists_per_shot(n_real, n_shots, r_det, rng, inv_tab=None):
    """纯噪声 hist_i：(n_real, n_shots, NBINS)，每 shot 计数 ∈[0,27]。"""
    if inv_tab is None:
        inv_tab = build_inv_table(r_det)
    out = np.zeros((n_real, n_shots, NBINS), dtype=np.int32)
    for s in range(n_shots):
        out[:, s, :] = noise_macro_hist_fast(
            n_real, N_PIX_MACRO, r_det, rng, inv_tab=inv_tab)
    return out


def binary_macro_stepping_per_shot(n_real, f_pix, n_shots, r_sig_unit, tgrid, r_amb,
                                   centers, rng, boost=1.0, tau_rc=TAU_RC, t_over=T_OVER,
                                   pde=PDE, jitter=JIT, resp_shape=RESP_SHAPE, resp_k=RESP_K):
    """信号+环境 hist_i：(n_real, n_shots, len(centers))。"""
    f_arr = np.tile(np.asarray(f_pix, float), int(n_shots))
    dt = tgrid[1] - tgrid[0]
    n_tr = f_arr.size
    n_pix = int(np.asarray(f_pix).size)
    nb = len(centers)
    k_max = int(np.ceil(20.0 * tau_rc / dt))
    phi = pde * spad_response_g(1.0 - np.exp(-np.arange(k_max + 1) * dt / tau_rc),
                                resp_shape, resp_k)
    age = np.full((n_real, n_tr), k_max, dtype=np.int32)
    tcov = np.full((n_real, n_tr), -1e30)
    cend = np.full((n_real, n_tr), -1e30)
    hist_i = np.zeros((n_real, n_shots, nb), dtype=np.int32)
    mu_all = (r_sig_unit[:, None] * f_arr[None, :] * boost + r_amb) * dt
    ib = 0
    for i in range(tgrid.size):
        t = tgrid[i]
        while ib < nb and centers[ib] < t:
            lit = _cover_commit(centers[ib], tcov, cend, t_over)
            hist_i[:, :, ib] = lit.reshape(n_real, n_shots, n_pix).sum(axis=2)
            ib += 1
        p = -np.expm1(-mu_all[i][None, :] * phi[age])
        fire = rng.random((n_real, n_tr)) < p
        age = np.minimum(age + 1, k_max)
        if fire.any():
            age[fire] = 1
            nf = int(fire.sum())
            tcov[fire] = t + (rng.normal(0.0, jitter, nf) if jitter > 0 else 0.0)
    while ib < nb:
        lit = _cover_commit(centers[ib], tcov, cend, t_over)
        hist_i[:, :, ib] = lit.reshape(n_real, n_shots, n_pix).sum(axis=2)
        ib += 1
    return hist_i


def hist_add_from_prefix(hist_i, n_shots):
    return hist_i[:, :n_shots, :].sum(axis=1)


def stats_from_hist_i(hist_i, n_shots_list=None, i0=None, i1=None):
    """由 N_SHOTS_MAX 发 hist_i 派生各 N 的 noise/bg/peak 充分统计。"""
    if n_shots_list is None:
        n_shots_list = N_SHOTS_LIST
    if i0 is None:
        i0 = I_STAT0
    if i1 is None:
        i1 = I_STAT1
    n_real, n_max, _ = hist_i.shape
    shot_nz = hist_i[:, :, i0:i1].mean(axis=2)
    out = {}
    for n in n_shots_list:
        hadd = hist_add_from_prefix(hist_i, n)
        a = hadd[:, i0:i1]
        bg = a.mean(axis=1)
        pk = a.max(axis=1)
        nz = shot_nz[:, :n].mean(axis=1)
        n_tr = N_PIX_MACRO * n
        out[n] = dict(
            n=n_real,
            noise_sum=float(nz.sum()), noise_sumsq=float((nz*nz).sum()),
            bg_sum=float(bg.sum()), bg_sumsq=float((bg*bg).sum()),
            peak_cnt=np.bincount(pk, minlength=n_tr + 2).astype(np.int64),
            # ★ v30：单条 hist_add 在统计窗内 152 个 bin 上的样本 std，
            #        累加后除以条数就是模块 10 的「hist 内 std 均值」
            hist_std_sum=float(a.std(axis=1).sum()),
        )
    return out


print("引擎与更新过程工具就绪：")
print("  · spad_binary_trace     —— 精确（v45 原样），基准")
print("  · noise_macro_hist_fast —— 快速 A（★v02 O(1) 直查表 + float32），仅纯环境光")
print("  · binary_macro_stepping —— 快速 B，同步时间步进，含信号")
print("  · noise_hists_per_shot / binary_macro_stepping_per_shot —— ★v10 hist_i")
print("  · hist_add_from_prefix / stats_from_hist_i —— ★v10 前缀和")
print("  · p_bin_equilibrium / r_det_for_noise / e_lambda_for_r_det —— noise 与环境光的精确互换")
# ===== 源自 PoD_esti_v30.ipynb cell 17 =====
def _noise_chunk_stats(m, n_tr, r_det, inv_tab, seed):
    """单个纯噪声分块；返回可直接归并的充分统计量。"""
    h = noise_macro_hist_fast(
        m, n_tr, r_det, np.random.default_rng(seed), inv_tab=inv_tab,
    )
    a = h[:, I_STAT0:I_STAT1]
    nz = a.mean(axis=1)
    return (float(nz.sum()), float((nz*nz).sum()),
            np.bincount(a.max(axis=1), minlength=n_tr + 2),
            float(a.std(axis=1).sum()))   # ★ v30：hist 内 std


def peak_stats_from_cnt(cnt):
    """由 peak 的 bincount 精确算出均值/标准差/各分位数。"""
    v = np.arange(cnt.size, dtype=float)
    n = cnt.sum()
    mean = (v * cnt).sum() / n
    var = (v*v * cnt).sum() / n - mean**2
    cum = np.cumsum(cnt) / n
    q = lambda p: float(np.searchsorted(cum, p))
    return dict(n=int(n), mean=mean, std=np.sqrt(max(var, 0.0)),
                p01=q(0.01), p50=q(0.50), p99=q(0.99),
                p999=q(0.999), p9999=q(0.9999),
                pmax=float(np.nonzero(cnt)[0].max()) if np.any(cnt) else 0.0)


def _atomic_savez(path, **kwargs):
    """先写临时文件再替换，避免中断写出半截缓存。
    注意：np.savez_compressed 若路径不以 .npz 结尾会自动追加 .npz，
    因此临时文件必须自带 .npz 后缀，否则 os.replace 会找不到文件（WinError 2）。
    """
    path = str(path)
    tmp = path + ".tmp.npz"
    np.savez_compressed(tmp, **kwargs)
    os.replace(tmp, path)

def _try_load_noise_cache(path, grid_key):
    if not os.path.exists(path):
        return None
    z = np.load(path, allow_pickle=True)
    if (int(z["n_mc"]) == N_MC_NOISE and list(z["n_shots_list"]) == list(N_SHOTS_LIST)
            and z["grid_key"].shape == grid_key.shape
            and np.allclose(z["grid_key"], grid_key)):
        _r = z["res"].item()
        # ★ v30：没有 hist_std 字段的一律当作旧缓存丢弃
        if any("hist_std" not in _r.get(_n, {}) for _n in N_SHOTS_LIST):
            return None
        return _r
    return None


def _save_noise_cache(path, res, grid_key):
    _atomic_savez(path, res=np.array(res, dtype=object),
                  n_mc=N_MC_NOISE, n_shots_list=np.array(N_SHOTS_LIST),
                  grid_key=grid_key)


def run_noise_scan(n_shots, noise_grid, n_mc, chunk, seed0=2000, verbose_every=5,
                   res=None, start_k=0, on_progress=None):
    """对一组【目标 noise】跑纯噪声 MC；支持从 start_k 断点续跑。"""
    n_tr = N_PIX_MACRO * n_shots
    ng = len(noise_grid)
    if res is None:
        res = {"n_shots": n_shots, "n_tr": n_tr,
               "noise_target": np.asarray(noise_grid, float),
               "r_det": np.zeros(ng), "e_lambda": np.zeros(ng),
               "p_eq": np.zeros(ng), "noise_mc": np.zeros(ng), "noise_std": np.zeros(ng),
               "hist_std": np.zeros(ng),   # ★ v30
               "peak_cnt": np.zeros((ng, n_tr + 2), dtype=np.int64),
               "done": np.zeros(ng, dtype=bool)}
    elif "done" not in res:
        # 兼容旧缓存：peak_cnt 有样本即视为已完成
        res["done"] = np.array([int(c.sum()) > 0 for c in res["peak_cnt"]], dtype=bool)

    t_start = time.time()
    for k, nt in enumerate(noise_grid):
        if k < start_k or bool(res["done"][k]):
            continue
        r_det = r_det_for_noise(float(nt), n_tr)
        res["r_det"][k] = r_det
        res["e_lambda"][k] = e_lambda_for_r_det(r_det)
        res["p_eq"][k] = p_bin_equilibrium(r_det)[0]
        inv_tab = build_inv_table(r_det)
        s1 = s2 = s3 = 0.0
        res["peak_cnt"][k][:] = 0
        specs = [
            (min(chunk, n_mc - s), n_tr, r_det, inv_tab, seed0 + 1000*k + s)
            for s in range(0, n_mc, chunk)
        ]
        if NOISE_WORKERS <= 1:
            parts = [_noise_chunk_stats(*spec) for spec in specs]
        else:
            with ThreadPoolExecutor(max_workers=NOISE_WORKERS) as pool:
                parts = list(pool.map(lambda x: _noise_chunk_stats(*x), specs))
        for p1, p2, pcnt, p3 in parts:
            s1 += p1; s2 += p2; s3 += p3
            res["peak_cnt"][k] += pcnt
        res["noise_mc"][k] = s1 / n_mc
        res["noise_std"][k] = np.sqrt(max(s2/n_mc - (s1/n_mc)**2, 0.0))
        res["hist_std"][k] = s3 / n_mc
        res["done"][k] = True
        if on_progress is not None:
            on_progress(res, k)
        if k == 0 or k == ng-1 or (k+1) % verbose_every == 0:
            el = time.time() - t_start
            remain = max(int((~res["done"]).sum()), 0)
            done_n = int(res["done"].sum())
            eta = el / max(done_n - start_k, 1) * remain
            pk = peak_stats_from_cnt(res["peak_cnt"][k])
            print(f"  [N_shots={n_shots} {done_n:>3d}/{ng}] 目标 noise={nt:>6.2f} → "
                  f"实测 {res['noise_mc'][k]:>6.3f}（E_λ={res['e_lambda'][k]:.4f}，"
                  f"≈{res['e_lambda'][k]/0.68*100:>5.0f} klux）  "
                  f"peak 中位={pk['p50']:>5.1f} 99.99%={pk['p9999']:>5.1f}  "
                  f"[已用 {el:.0f}s, 剩约 {eta:.0f}s]")
    return res





def run_noise_scan_v20_bg(bg_grid, n_mc, chunk, seed0=2000, verbose_every=1,
                          res_all=None, on_progress=None):
    """★ v20：按统一 bg 网格扫；对每个 N 单独设 noise_amb=bg/N。

    返回 {N: res_dict}；res_dict 字段兼容旧作图（noise_mc=实测 bg，noise_target=目标 bg）。
    """
    grid = np.asarray(bg_grid, float)
    ng = len(grid)
    if res_all is None:
        res_all = {}
    for n in N_SHOTS_LIST:
        n_tr = N_PIX_MACRO * n
        if n not in res_all:
            res_all[n] = {
                "n_shots": n, "n_tr": n_tr,
                "noise_target": grid.copy(),                 # 目标 bg
                "noise_amb_target": np.round(grid / n, 6),  # 对应单次 noise
                "r_det": np.zeros(ng), "e_lambda": np.zeros(ng), "p_eq": np.zeros(ng),
                "noise_mc": np.zeros(ng), "noise_std": np.zeros(ng),  # = bg
                "noise_amb_mc": np.zeros(ng), "noise_amb_std": np.zeros(ng),
                "peak_cnt": np.zeros((ng, n_tr + 2), dtype=np.int64),
                "done": np.zeros(ng, dtype=bool),
            }
        elif "done" not in res_all[n]:
            res_all[n]["done"] = np.array(
                [int(c.sum()) > 0 for c in res_all[n]["peak_cnt"]], dtype=bool)

    t_start = time.time()
    jobs = [(n, k) for n in N_SHOTS_LIST for k in range(ng)
            if not bool(res_all[n]["done"][k])]
    n_jobs = len(jobs)
    print(f"v20 噪声扫描：{ng} bg × N={list(N_SHOTS_LIST)} = {ng*len(N_SHOTS_LIST)} 档，"
          f"待算 {n_jobs}，每档 {n_mc:,} MC", flush=True)
    for ji, (n, k) in enumerate(jobs):
        bg_t = float(grid[k])
        nt_amb = bg_t / n
        r_det = float(r_det_for_noise(float(nt_amb), N_PIX_MACRO))
        e_lam = float(e_lambda_for_r_det(r_det))
        p_eq = float(p_bin_equilibrium(r_det)[0])
        inv_tab = build_inv_table(r_det)

        acc = dict(noise_sum=0.0, noise_sumsq=0.0, bg_sum=0.0, bg_sumsq=0.0,
                   peak_cnt=np.zeros(N_PIX_MACRO * n + 2, dtype=np.int64), nn=0)
        done_m, part = 0, 0
        while done_m < n_mc:
            m = min(chunk, n_mc - done_m)
            seeds = [seed0 + 10007 * (n * 1000 + k) + 104729 * part + 17 * t
                     for t in range(NOISE_WORKERS)]
            ms = [m // NOISE_WORKERS + (1 if t < m % NOISE_WORKERS else 0)
                  for t in range(NOISE_WORKERS)]

            def _one(args, _n=n, _rd=r_det, _it=inv_tab):
                mm, sd = args
                if mm <= 0:
                    return None
                rng = np.random.default_rng(sd)
                hi = noise_hists_per_shot(mm, _n, _rd, rng, inv_tab=_it)
                return stats_from_hist_i(hi, n_shots_list=[_n])

            with ThreadPoolExecutor(max_workers=NOISE_WORKERS) as pool:
                parts = list(pool.map(_one, zip(ms, seeds)))
            for st in parts:
                if st is None:
                    continue
                b = st[n]
                acc["noise_sum"] += b["noise_sum"]; acc["noise_sumsq"] += b["noise_sumsq"]
                acc["bg_sum"] += b["bg_sum"]; acc["bg_sumsq"] += b["bg_sumsq"]
                acc["peak_cnt"] += b["peak_cnt"]; acc["nn"] += b["n"]
            done_m += m; part += 1

        R = res_all[n]; nn = max(acc["nn"], 1)
        R["r_det"][k] = r_det; R["e_lambda"][k] = e_lam; R["p_eq"][k] = p_eq
        R["noise_amb_mc"][k] = acc["noise_sum"] / nn
        R["noise_amb_std"][k] = float(np.sqrt(max(
            acc["noise_sumsq"]/nn - (acc["noise_sum"]/nn)**2, 0.0)))
        R["noise_mc"][k] = acc["bg_sum"] / nn
        R["noise_std"][k] = float(np.sqrt(max(
            acc["bg_sumsq"]/nn - (acc["bg_sum"]/nn)**2, 0.0)))
        R["peak_cnt"][k] = acc["peak_cnt"]
        R["done"][k] = True
        if on_progress is not None:
            on_progress(res_all, n, k)
        if (ji % verbose_every) == 0 or ji == n_jobs - 1:
            el = time.time() - t_start
            eta = el / (ji + 1) * (n_jobs - ji - 1)
            pk = peak_stats_from_cnt(R["peak_cnt"][k])
            print(f"  [{ji+1}/{n_jobs}] N={n} bg={bg_t:.2f}（amb={nt_amb:.3f}）→ "
                  f"bg_mc={R['noise_mc'][k]:.3f} peakμ={pk['mean']:.2f}  "
                  f"已用 {el/60:.1f} min，预计剩余 {eta/60:.1f} min", flush=True)
    return res_all


# 兼容旧名
run_noise_scan_v20_amb = run_noise_scan_v20_bg


# ---- 估算总耗时并开跑（主缓存 + fallback + 增量检查点）----

# ===== 源自 PoD_esti_v30.ipynb cell 18 =====
def far_threshold_from_cnt(cnt, target_far):
    """由 peak 的 bincount 求满足 P(peak ≥ T) < target_far 的最小整数 T。

    全程用【整数计数】比较（n_ge < target_far·n），避免浮点边界误判。
    返回 (T, 该 T 处实测 FAR, 该 T 处越阈次数, 生存函数数组)。
    """
    n = int(cnt.sum())
    n_ge = np.concatenate([[n], n - np.cumsum(cnt)[:-1]])
    lim = target_far * n
    ok = np.where(n_ge < lim)[0]
    sf = n_ge / n
    if ok.size == 0:
        return int(cnt.size), 0.0, 0, sf
    T = int(ok[0])
    return T, float(sf[T]), int(n_ge[T]), sf


def far_threshold_binom_indep(n_tr, p_bin, n_bins, target_far):
    """独立 Binomial 近似阈值（保守对照）。"""
    a_bin = 1.0 - (1.0 - target_far) ** (1.0 / n_bins)
    T = 0
    while T <= n_tr and _binom.sf(T - 1, n_tr, p_bin) > a_bin:
        T += 1
    return T

# ===== 源自 PoD_esti_v30.ipynb cell 22 =====
# ---- PoD 专用子窗（只计算信号附近，前方保留暖机）----
POD_T_LO = T0_SIG - POD_WARM_NS * 1e-9
POD_T_HI = T0_SIG + SIG_POST_NS * 1e-9
TF_POD = np.arange(POD_T_LO, POD_T_HI, DT_FINE)
_sigmask = (TC_NS >= T0_SIG_NS - SIG_PRE_NS) & (TC_NS <= T0_SIG_NS + SIG_POST_NS)
IDX_SIG = np.where(_sigmask)[0]
CENTERS_SIG = CENTERS[IDX_SIG]
R_SIG_UNIT_POD = signal_photon_rate_fine(ECHO0, 1.0, TF_POD)
_NPH_BASE = np.trapezoid(R_SIG_UNIT_POD, TF_POD) * F_VALS.sum()

print(f"PoD 子窗：{POD_T_LO*1e9:.1f}–{POD_T_HI*1e9:.1f} ns，{TF_POD.size} 个细网格步")
print(f"每种 N_shots 对统一 BG_GRID 求解："
      f"{[len(NOISE_GRID[n]) for n in N_SHOTS_LIST]} 档，bg 步长 0.25；noise_amb=bg/N")
print(f"并行：外层 POD_BIN_WORKERS={POD_BIN_WORKERS} × 内层 POD_WORKERS={POD_WORKERS}；"
      f"MC 分块 POD_MC_CHUNK={POD_MC_CHUNK}；临界验证 {N_MC_POD_VERIFY:,} 次/点")


def _peaks_chunk(boost, n_shots, r_amb, n_real, seed):
    """★ v20：按当前 N 仿 n_shots 发（r_amb 已对应 noise=bg/N）。"""
    rng = np.random.default_rng(seed)
    hist_i = binary_macro_stepping_per_shot(
        n_real, F_VALS, n_shots, R_SIG_UNIT_POD, TF_POD, r_amb, CENTERS_SIG,
        rng, boost=boost,
    )
    return hist_i.sum(axis=1).max(axis=1)


def sig_peaks(boost, n_shots, r_amb, n_real, seed):
    """兼容入口：大 n_real 时按 POD_MC_CHUNK 切开，用 POD_WORKERS 并行。"""
    n_real = int(n_real)
    if n_real <= POD_MC_CHUNK or POD_WORKERS <= 1:
        return _peaks_chunk(boost, n_shots, r_amb, n_real, seed)
    specs = []
    for s in range(0, n_real, POD_MC_CHUNK):
        m = min(POD_MC_CHUNK, n_real - s)
        specs.append((boost, n_shots, r_amb, m, seed + 104729 * s))
    with ThreadPoolExecutor(max_workers=POD_WORKERS) as pool:
        parts = list(pool.map(lambda sp: _peaks_chunk(*sp), specs))
    return np.concatenate(parts)


def _eval_mc_jobs(job_specs, n_shots, r_amb):
    """统一并行入口。
    job_specs: [(boost, n_real, seed), ...] → 与输入等长的 peak 数组列表。
    每个 job 再按 POD_MC_CHUNK 切开，全部丢进同一个线程池，避免嵌套池。
    """
    flat, owners = [], []
    for j, (boost, n_real, seed) in enumerate(job_specs):
        n_real = int(n_real)
        if n_real <= 0:
            continue
        for s in range(0, n_real, POD_MC_CHUNK):
            m = min(POD_MC_CHUNK, n_real - s)
            flat.append((float(boost), n_shots, r_amb, m, int(seed) + 104729 * s))
            owners.append(j)
    out = [None] * len(job_specs)
    if not flat:
        return [np.zeros(0, dtype=int) for _ in job_specs]
    if POD_WORKERS <= 1:
        parts = [_peaks_chunk(*sp) for sp in flat]
    else:
        with ThreadPoolExecutor(max_workers=POD_WORKERS) as pool:
            parts = list(pool.map(lambda sp: _peaks_chunk(*sp), flat))
    buckets = [[] for _ in job_specs]
    for own, pk in zip(owners, parts):
        buckets[own].append(pk)
    for j, segs in enumerate(buckets):
        out[j] = np.concatenate(segs) if segs else np.zeros(0, dtype=int)
    return out


def _eval_boost_grid(boosts, n_shots, r_amb, n_real, seed0):
    """并行评估若干独立能量点；返回每点的 peak 样本。"""
    boosts = np.asarray(boosts, float)
    jobs = [(float(b), n_real, seed0 + 1009 * i) for i, b in enumerate(boosts)]
    return _eval_mc_jobs(jobs, n_shots, r_amb)


def _isotonic(p):
    """简单保序：消除有限 MC 导致的局部 PoD 下降。"""
    return np.maximum.accumulate(np.asarray(p, float))


def _crossing_logboost(boosts, pod, level):
    """在 log10(boost) 上找首次跨越；越界返回 NaN。"""
    order = np.argsort(boosts)
    x = np.log10(np.asarray(boosts)[order])
    p = _isotonic(np.asarray(pod)[order])
    if p[0] >= level:
        return float(x[0])
    if p[-1] < level:
        return np.nan
    i = int(np.searchsorted(p, level))
    dp = p[i] - p[i - 1]
    w = 0.5 if dp <= 0 else (level - p[i - 1]) / dp
    return float(x[i - 1] + w * (x[i] - x[i - 1]))


def _probit_fit(boosts, pod, n_real):
    """拟合 Phi^-1(PoD) = slope*log10(boost) + intercept。"""
    boosts = np.asarray(boosts, float)
    success = np.rint(np.asarray(pod, float) * n_real)
    p = (success + 0.5) / (n_real + 1.0)
    transition = (p > 0.01) & (p < 0.99)
    if transition.sum() < 3:
        transition = np.argsort(np.abs(p - 0.5))[:min(5, len(p))]
    x = np.log10(boosts[transition])
    z = _norm.ppf(p[transition])
    slope, intercept = np.polyfit(x, z, 1)
    return float(slope), float(intercept)


def _pk_to_record(boost, pk, T, n_shots):
    return {
        "boost": float(boost),
        "pod": float((pk >= T).mean()),
        "peak_mean": float(pk.mean()),
        "peak_std": float(pk.std()),
        "peak_cnt": np.bincount(pk, minlength=N_PIX_MACRO * n_shots + 1),
        "n_verify": int(pk.size),
    }


def _probit_fit_local(boosts, pod, n_real, level, half_decade=0.6):
    """只用经验交点附近的点做 probit 拟合。

    ★ v30 修复：全域拟合会被 4 个数量级上的饱和点（PoD≈0 与 PoD≈1）拽偏，
    5% FAR 档的初值经常偏半个数量级以上。
    """
    x0 = _crossing_logboost(boosts, pod, level)
    if not np.isfinite(x0):
        return _probit_fit(boosts, pod, n_real)
    x = np.log10(np.asarray(boosts, float))
    sel = np.abs(x - x0) <= half_decade
    if sel.sum() >= 3:
        return _probit_fit(np.asarray(boosts, float)[sel],
                           np.asarray(pod, float)[sel], n_real)
    return _probit_fit(boosts, pod, n_real)


def _next_root_guess(hist, level, slope, max_step=0.5):
    """由已验证的 (log10 boost, PoD) 历史给出下一个试探点。

    有括号（一点低于目标、一点高于目标）就在 probit 空间做割线，
    割线跑出括号则退回二分；没有括号就用拟合斜率做 Newton 步并主动向外扩。
    """
    n = float(N_MC_POD_VERIFY)
    _clip = lambda p: min(max(p, 0.5 / n), 1.0 - 0.5 / n)
    below = [h for h in hist if h[1] < level]
    above = [h for h in hist if h[1] > level]
    if below and above:
        lo = max(below, key=lambda h: h[0])
        hi = min(above, key=lambda h: h[0])
        if hi[0] > lo[0]:
            zl, zh, zt = _norm.ppf(_clip(lo[1])), _norm.ppf(_clip(hi[1])), _norm.ppf(level)
            x = (lo[0] + (zt - zl) / (zh - zl) * (hi[0] - lo[0])
                 if zh > zl else 0.5 * (lo[0] + hi[0]))
            if not (lo[0] < x < hi[0]):
                x = 0.5 * (lo[0] + hi[0])
            return float(x)
    x0, p0 = hist[-1][0], hist[-1][1]
    s = slope if (slope and slope > 0) else 2.0
    dx = float(np.clip((_norm.ppf(level) - _norm.ppf(_clip(p0))) / s, -max_step, max_step))
    if dx == 0.0:
        dx = max_step if p0 < level else -max_step
    return float(x0 + dx)


def _verify_critical_batch(cands, n_shots, r_amb, seed0):
    """多轮批量迭代求根，把每个 (FAR, PoD 等级) 临界点解到验证 PoD 落进容差。

    ★ v30 修复：v20 只做一次 Newton 步、步长夹在 ±0.25 decade，初值偏 0.5 decade
    以上时根本追不回来，却仍然无条件接受结果。表现是模块 7 的临界能量曲线出现
    3–5 倍的毛刺，验证 PoD 实测 0.68 或 1.000 而不是 0.90。
    现在每轮把所有活跃候选一起并行评估（保持吞吐），最多 POD_VERIFY_ROUNDS 轮，
    最终取历史上最接近目标的那个点，并把 pod_err 一并存进记录备查。
    """
    if not cands:
        return {}
    state = [{"c": c, "i": i, "x": float(np.log10(c["boost"])), "hist": [], "done": False}
             for i, c in enumerate(cands)]

    for rnd in range(POD_VERIFY_ROUNDS):
        act = [s for s in state if not s["done"]]
        if not act:
            break
        pks = _eval_mc_jobs(
            [(float(10.0 ** s["x"]), N_MC_POD_VERIFY,
              seed0 + 7919 * s["i"] + 1_000_003 * rnd) for s in act],
            n_shots, r_amb,
        )
        for s, pk in zip(act, pks):
            c = s["c"]
            rec = _pk_to_record(10.0 ** s["x"], pk, c["T"], n_shots)
            s["hist"].append((s["x"], rec["pod"], rec))
            if abs(rec["pod"] - c["level"]) <= POD_VERIFY_TOL:
                s["done"] = True
                continue
            nx = _next_root_guess(s["hist"], c["level"], c["slope"])
            if not np.isfinite(nx):
                s["done"] = True
            else:
                s["x"] = nx

    finals = {}
    for s in state:
        c = s["c"]
        best = min(s["hist"], key=lambda h: abs(h[1] - c["level"]))
        rec = dict(best[2])
        rec["verify_rounds"] = len(s["hist"])
        rec["pod_err"] = float(best[1] - c["level"])
        finals[(c["tag"], f"{c['level']:.2f}")] = rec
    return finals


def solve_pod_noise(n_shots, k, seed0):
    """求一个 noise 档、全部 FAR 阈值下的 PoD50/90 临界点。"""
    R, Tr = NOISE_RES[n_shots], THRESH[n_shots]
    nt = float(R["noise_target"][k])
    n_tr = int(R["n_tr"])
    r_amb = float(R["r_det"][k] / PDE)
    # ★ v30：只对 POD_FARS 求 PoD 交点（阈值本身七条都在 THRESH 里）
    T_map = {FAR_TAG[far]: int(Tr["T" + FAR_TAG[far]][k]) for far in POD_FARS}
    if max(T_map.values()) > n_tr:
        return (n_shots, nt), {
            "noise": float(R["noise_mc"][k]), "e_lambda": float(R["e_lambda"][k]),
            "n_tr": n_tr, "T_map": T_map, "critical": {}, "invalid": "阈值超过二值硬上限",
        }

    coarse_boost = np.logspace(POD_LOG_BOOST_MIN, POD_LOG_BOOST_MAX, N_POD_COARSE)
    coarse_pk = _eval_boost_grid(
        coarse_boost, n_shots, r_amb, N_MC_POD_COARSE, seed0,
    )
    coarse_pod = {
        tag: np.array([(pk >= T).mean() for pk in coarse_pk])
        for tag, T in T_map.items()
    }

    roots0 = []
    for tag in T_map:
        for level in POD_LEVELS:
            x0 = _crossing_logboost(coarse_boost, coarse_pod[tag], level)
            if np.isfinite(x0):
                roots0.append(x0)
    if roots0:
        local_x = np.unique(np.concatenate([
            np.linspace(x0 - POD_LOCAL_HALF_DECADE, x0 + POD_LOCAL_HALF_DECADE,
                        N_POD_LOCAL_PER_ROOT)
            for x0 in roots0
        ]))
        local_boost = 10.0 ** local_x
        local_pk = _eval_boost_grid(
            local_boost, n_shots, r_amb, N_MC_POD_LOCAL, seed0 + 500_000,
        )
    else:
        local_boost = np.array([], float)
        local_pk = []

    critical = {tag: {} for tag in T_map}
    curve = {}
    cands = []
    for tag, T in T_map.items():
        boosts_fit = np.concatenate([coarse_boost, local_boost])
        pod_fit = np.concatenate([
            coarse_pod[tag],
            np.array([(pk >= T).mean() for pk in local_pk]) if len(local_pk) else np.array([], float),
        ])
        order = np.argsort(boosts_fit)
        boosts_fit, pod_fit = boosts_fit[order], pod_fit[order]
        curve[tag] = {"boost": boosts_fit, "pod": pod_fit}
        for level in POD_LEVELS:
            # ★ v30：逐 level 做局部 probit 拟合，别让远处的饱和点拽偏初值
            slope, intercept = _probit_fit_local(
                boosts_fit, pod_fit, N_MC_POD_LOCAL, level)
            x_root = (_norm.ppf(level) - intercept) / slope if slope > 0 else np.nan
            x_emp = _crossing_logboost(boosts_fit, pod_fit, level)
            if np.isfinite(x_emp) and (not np.isfinite(x_root)
                                       or abs(x_root - x_emp) > 0.5):
                x_root = x_emp   # 拟合外推得离谱时，经验交点更可信
            if not np.isfinite(x_root):
                critical[tag][f"{level:.2f}"] = None
                continue
            cands.append({
                "tag": tag, "level": level, "T": T, "slope": slope,
                "boost": float(10.0 ** x_root),
            })

    verified = _verify_critical_batch(cands, n_shots, r_amb, seed0 + 700_000)
    for (tag, lk), rec in verified.items():
        critical[tag][lk] = rec
    for tag in T_map:
        for level in POD_LEVELS:
            critical[tag].setdefault(f"{level:.2f}", None)

    return (n_shots, nt), {
        "noise": float(R["noise_mc"][k]),
        "noise_target": nt,
        "e_lambda": float(R["e_lambda"][k]),
        "n_tr": n_tr,
        "T_map": T_map,
        "curve": curve,
        "critical": critical,
    }


# ---- 对完整 0.25-bg 网格求解；主缓存 + fallback + 增量检查点 ----
_pod_grid_key = np.concatenate([NOISE_GRID[n] for n in N_SHOTS_LIST])


def _try_load_pod_cache(path, grid_key):
    """仅当噪声网格、MC 精度与 FAR 标签集合全部一致时才接受缓存。"""
    if not os.path.exists(path):
        return None
    z = np.load(path, allow_pickle=True)
    if "far_tags" not in z.files:
        return None  # 旧缓存无 FAR 标签 → 一律作废
    if (
        np.array_equal(z["grid_key"], grid_key)
        and list(z["far_tags"]) == list(FAR_TAGS)
        and int(z["n_coarse"]) == N_MC_POD_COARSE
        and int(z["n_local"]) == N_MC_POD_LOCAL
        and int(z["n_verify"]) == N_MC_POD_VERIFY
    ):
        return z["res"].item()
    return None


def _save_pod_cache(path, res, grid_key):
    _atomic_savez(
        path, res=np.array(res, dtype=object),
        grid_key=grid_key,
        far_tags=np.array(FAR_TAGS),
        n_coarse=N_MC_POD_COARSE, n_local=N_MC_POD_LOCAL,
        n_verify=N_MC_POD_VERIFY,
    )


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
