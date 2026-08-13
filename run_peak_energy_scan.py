# -*- coding: utf-8 -*-
"""信号能量扫描（工作名 `peak_vs_energy`）——多进程 + 落盘缓存 + 断点续跑。

回答的问题：**在没有环境光（bg=0）的条件下，二值 SPAD 宏像元累加直方图的
peak 与半高全宽（FWHM）如何随信号能量变化**，从单光子稀疏区一路到深饱和。

物理内核整个复用 `pod_esti_v30_core.py`，本脚本不定义任何新物理。

口径（与 PoD_esti 全项目一致）：
* 能量用倍率 `boost` 表示，`boost=1` = 默认 ρ=0.1 场景的回波强度。
  用户只关心比例，所以横轴就用 boost（另附「等效 ρ」与「峰值光子率」换算）。
* `hist_add` = N 发累加直方图，取值 0…n_tr，n_tr = 27 × N_shots。
* `peak` = `hist_add` 在统计窗 152 个 bin 内的最大计数。
* `FWHM` = **单次实现的 `hist_add` 波形**的半高全宽（线性插值到亚 bin，
  1 bin = 1 ns），再对多次蒙卡求平均。取「包含峰位的那一段连续过半高区间」，
  这样即使将来把 bg 开到非零、远处出现噪声尖峰也不会串进来。

一次 4 发仿真通过前缀和白拿 N=1/2/4 三档，故 N_shots 维度几乎不额外花钱。

落盘的是**充分统计量**，不是原始样本：
* `peak_cnt_{N}`   (nB, n_tr+1)  peak 的完整分布
* `fwhm_cnt_{N}`   (nB, NFW)     FWHM 的分布（0.25 ns 一格）
* `fwhm_sum/sumsq/nval_{N}`      FWHM 的均值与标准差
* `wave_sum_{N}`   (nB, NBINS)   平均波形（画总览图用）

用法：
    python run_peak_energy_scan.py --workers 20 --n-mc 20000
"""
from __future__ import annotations

import argparse
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

os.environ.setdefault("POD_CORE_QUIET", "1")
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

import numpy as np

import pod_esti_v30_core as core

# ---------------------------------------------------------------- 扫描配置
# 像斑在宏像元 27 个 SPAD 上的分配方式。两种模式各存各的缓存，互不覆盖。
#   "uniform" —— 27 个 SPAD 平分总收集比例 Σf_pix。各轨迹点亮概率严格相同，
#                单 bin 计数是标准二项分布，方差极大值落在占比 0.5。
#   "real"    —— 真实像斑轮廓，f_pix 最大/最小相差约 241 倍。各轨迹 p_t 悬殊，
#                单 bin 计数是泊松二项分布，方差被压低、极大值左移到占比 0.2–0.4。
# 两种模式的 Σf_pix 相同，所以低能段（远未饱和时）的 peak 均值几乎一致，
# 差别只在饱和过程的快慢与涨落的大小。
# 临时切换不必改这里：set PVE_FPIX_MODE=real 后再运行即可。
F_PIX_MODE = os.environ.get("PVE_FPIX_MODE", "uniform")

CACHE = f"peak_vs_energy_cache_{F_PIX_MODE}.npz"
CKPT = f"peak_vs_energy_cache_{F_PIX_MODE}.partial.npz"

BG_SIG = 0.0            # 环境光底噪（用户要求 bg=0；改成非 0 会让缓存键失配、自动重算）
N_SHOTS_MAX = 4
N_SHOTS_LIST = (1, 2, 4)

MC_CHUNK = 2500         # 每块实现数，控制峰值内存（20 进程时约 3 MB/块）
CHECKPOINT_EVERY = 6    # 每完成几个能量点落一次盘

FWHM_MIN_PEAK = 4       # peak 低于这个计数就不给 FWHM（半高只有 1~2 个计数，无意义）
FWHM_MAX_NS = 64.0      # FWHM 分布直方图的上界
FWHM_BIN_NS = 0.25
NFW = int(round(FWHM_MAX_NS / FWHM_BIN_NS))


def build_boost_grid(mode: str = F_PIX_MODE) -> np.ndarray:
    """非均匀采样网格：坐标轴仍然线性，只是把 MC 点放在曲线走得快的地方。

    两种模式的饱和位置差了近两个数量级，所以网格必须分别设计：
    均匀分配时所有 SPAD 同时逼近饱和，boost≈0.3 就封顶；
    真实分配里最暗的 SPAD 要等到 boost≈20 才点亮，压缩过程拖长约 3 个数量级。
    两种模式的末段都要一直铺到 1e4 —— peak 早就封顶了，但 FWHM 每十倍能量还在涨约 4.4 ns。
    """
    if mode == "uniform":
        g = np.concatenate([
            np.linspace(0.0, 0.005, 21),        # peak 0→11，查线性
            np.linspace(0.005, 0.05, 19)[1:],   # 占比 0.12→0.69，σ 拱顶（占比 0.5）在这一段
            np.linspace(0.05, 0.15, 11)[1:],    # 占比 0.69→0.96
            np.linspace(0.15, 0.5, 8)[1:],      # 占比 0.96→0.997，封顶
            np.geomspace(0.5, 1.0e4, 20)[1:],   # 只剩 FWHM 在长
        ])
    elif mode == "real":
        g = np.concatenate([
            np.linspace(0.0, 0.005, 21),
            np.linspace(0.005, 0.05, 19)[1:],
            np.linspace(0.05, 0.25, 11)[1:],
            np.linspace(0.25, 1.0, 7)[1:],
            np.linspace(1.0, 5.0, 6)[1:],
            np.geomspace(5.0, 1.0e4, 17)[1:],
        ])
    else:
        raise ValueError(f"未知的 F_PIX_MODE：{mode!r}（只支持 'uniform' / 'real'）")
    return np.unique(np.round(g, 10))


BOOSTS = build_boost_grid()
NB = BOOSTS.size

# 峰值宏像元 27 个 SPAD 的真实空间收集比例（像斑轮廓，最大/最小相差约 241 倍）
F_PIX_REAL = core.FPIX[:, core.M_PEAK * core.MACRO_BY:(core.M_PEAK + 1) * core.MACRO_BY].ravel()
# 均匀模式：总量不变，27 个 SPAD 平分
F_PIX_UNIFORM = np.full(core.N_PIX_MACRO, F_PIX_REAL.sum() / core.N_PIX_MACRO)
F_PIX = F_PIX_UNIFORM if F_PIX_MODE == "uniform" else F_PIX_REAL

# 引擎版本号。2026-08-10 修了 binary_macro_stepping* 在饱和时丢过阈覆盖的 bug
# （详见 pod_esti_v30_core._cover_commit）。改这个标签就是为了让修复前跑的缓存自动失配，
# 不会被当成有效数据继续用。物理引擎再动，这里必须跟着加版本。
ENGINE_TAG = "cov2"

GRID_KEY = (f"eng={ENGINE_TAG}|fpix={F_PIX_MODE}|bg={BG_SIG:g}|nsmax={N_SHOTS_MAX}|nb={NB}|"
            f"b0={BOOSTS[0]:.6g}|b1={BOOSTS[-1]:.6g}|"
            f"fmin={FWHM_MIN_PEAK}|fw={FWHM_BIN_NS:g}/{FWHM_MAX_NS:g}|"
            f"stat=[{core.I_STAT0},{core.I_STAT1})")


# ---------------------------------------------------------------- FWHM
def fwhm_ns_batch(h: np.ndarray, min_peak: int) -> np.ndarray:
    """一批 hist_add 波形的半高全宽 [ns]，全向量化。

    取**包含峰位的那一段连续过半高区间**，两端用线性插值补到亚 bin。
    peak < min_peak 的返回 NaN。
    """
    h = np.asarray(h, dtype=np.float64)
    nr, nb = h.shape
    pk = h.max(axis=1)
    half = 0.5 * pk
    i_pk = h.argmax(axis=1)
    rows = np.arange(nr)

    below = h < half[:, None]                       # 半高以下
    idx = np.arange(nb)[None, :]

    # 峰位左侧最近的一个「半高以下」bin；没有则 -1
    L = np.maximum.accumulate(np.where(below, idx, -1), axis=1)[rows, i_pk]
    # 峰位右侧最近的一个「半高以下」bin；没有则 nb
    R = np.minimum.accumulate(np.where(below, idx, nb)[:, ::-1], axis=1)[:, ::-1][rows, i_pk]

    okL = L >= 0
    hL = h[rows, np.clip(L, 0, nb - 1)]
    hL1 = h[rows, np.clip(L + 1, 0, nb - 1)]
    xl = np.where(okL, L + (half - hL) / np.maximum(hL1 - hL, 1e-12), 0.0)

    okR = R <= nb - 1
    hR = h[rows, np.clip(R, 0, nb - 1)]
    hR1 = h[rows, np.clip(R - 1, 0, nb - 1)]
    xr = np.where(okR, R - (half - hR) / np.maximum(hR1 - hR, 1e-12), float(nb - 1))

    w = (xr - xl) * float(core.BIN_W * 1e9)         # bin → ns
    return np.where(pk >= min_peak, w, np.nan)


# ---------------------------------------------------------------- 单个能量点
def _one_boost(args):
    ib, boost, n_mc, seed0 = args
    rng = np.random.default_rng(seed0)
    n_pix = core.N_PIX_MACRO
    nbins = core.NBINS

    acc = {}
    for n in N_SHOTS_LIST:
        n_tr = n_pix * n
        acc[n] = dict(
            peak_cnt=np.zeros(n_tr + 1, dtype=np.int64),
            fwhm_cnt=np.zeros(NFW, dtype=np.int64),
            fwhm_sum=0.0, fwhm_sumsq=0.0, fwhm_nval=0,
            wave_sum=np.zeros(nbins, dtype=np.float64),
        )

    done = 0
    while done < n_mc:
        m = int(min(MC_CHUNK, n_mc - done))
        hist_i = core.binary_macro_stepping_per_shot(
            m, F_PIX, N_SHOTS_MAX, core.R_SIG_UNIT_GEN, core.TF_GEN,
            BG_SIG, core.CENTERS, rng, boost=float(boost))
        for n in N_SHOTS_LIST:
            hadd = core.hist_add_from_prefix(hist_i, n)          # (m, nbins)
            a = hadd[:, core.I_STAT0:core.I_STAT1]
            n_tr = n_pix * n
            acc[n]["peak_cnt"] += np.bincount(a.max(axis=1), minlength=n_tr + 1)
            acc[n]["wave_sum"] += hadd.sum(axis=0)

            w = fwhm_ns_batch(a, FWHM_MIN_PEAK)
            v = w[np.isfinite(w)]
            if v.size:
                acc[n]["fwhm_sum"] += float(v.sum())
                acc[n]["fwhm_sumsq"] += float((v * v).sum())
                acc[n]["fwhm_nval"] += int(v.size)
                k = np.clip((v / FWHM_BIN_NS).astype(np.int64), 0, NFW - 1)
                acc[n]["fwhm_cnt"] += np.bincount(k, minlength=NFW)
        done += m

    return ib, {n: {k: (val.copy() if isinstance(val, np.ndarray) else val)
                    for k, val in acc[n].items()} for n in N_SHOTS_LIST}


# ---------------------------------------------------------------- 缓存
def _empty(n_mc):
    R = {"boosts": BOOSTS, "done": np.zeros(NB, dtype=bool),
         "n_mc": np.int64(n_mc), "bg": np.float64(BG_SIG),
         "f_pix_mode": F_PIX_MODE, "f_pix": F_PIX,
         "grid_key": GRID_KEY, "fwhm_min_peak": np.int64(FWHM_MIN_PEAK),
         "fwhm_bin_ns": np.float64(FWHM_BIN_NS), "fwhm_max_ns": np.float64(FWHM_MAX_NS)}
    for n in N_SHOTS_LIST:
        n_tr = core.N_PIX_MACRO * n
        R[f"peak_cnt_{n}"] = np.zeros((NB, n_tr + 1), dtype=np.int64)
        R[f"fwhm_cnt_{n}"] = np.zeros((NB, NFW), dtype=np.int64)
        R[f"fwhm_sum_{n}"] = np.zeros(NB)
        R[f"fwhm_sumsq_{n}"] = np.zeros(NB)
        R[f"fwhm_nval_{n}"] = np.zeros(NB, dtype=np.int64)
        R[f"wave_sum_{n}"] = np.zeros((NB, core.NBINS))
    return R


def _load(n_mc):
    for path in (CACHE, CKPT):
        if not os.path.exists(path):
            continue
        try:
            z = np.load(path, allow_pickle=False)
            if str(z["grid_key"]) != GRID_KEY or int(z["n_mc"]) != int(n_mc):
                print(f"  {path} 的缓存键与当前配置不符，忽略")
                continue
            R = {k: z[k] for k in z.files}
            R["grid_key"] = str(z["grid_key"])
            R["f_pix_mode"] = str(z["f_pix_mode"])
            print(f"  从 {path} 载入 {int(R['done'].sum())}/{NB} 档，断点续跑")
            return R
        except Exception as e:                                  # noqa: BLE001
            print(f"  {path} 读取失败（{e}），忽略")
    return _empty(n_mc)


def _save(R, path):
    tmp = path + ".tmp.npz"
    np.savez_compressed(tmp, **R)
    try:
        os.replace(tmp, path)
    except PermissionError:
        # Windows 下目标文件被别的进程占用（多半是还开着这份缓存、没重启的 Jupyter 内核）时，
        # os.replace 会被直接拒绝。绝不能让十几分钟的扫描结果就这样丢掉：
        # 退一步写到 <path>.new.npz 并大声提示，关掉占用方后手动改名即可。
        alt = path + ".new.npz"
        os.replace(tmp, alt)
        print(f"  [警告] 无法覆盖 {path}：文件被其他进程占用。")
        print(f"         结果已完整写入 {alt}，关掉占用它的进程（通常是 Jupyter 内核）后改名即可。")


def _apply(R, ib, res):
    for n in N_SHOTS_LIST:
        d = res[n]
        R[f"peak_cnt_{n}"][ib] = d["peak_cnt"]
        R[f"fwhm_cnt_{n}"][ib] = d["fwhm_cnt"]
        R[f"fwhm_sum_{n}"][ib] = d["fwhm_sum"]
        R[f"fwhm_sumsq_{n}"][ib] = d["fwhm_sumsq"]
        R[f"fwhm_nval_{n}"][ib] = d["fwhm_nval"]
        R[f"wave_sum_{n}"][ib] = d["wave_sum"]
    R["done"][ib] = True


# ---------------------------------------------------------------- 主流程
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=20)
    ap.add_argument("--n-mc", type=int, default=20000)
    a = ap.parse_args()

    n_tr_max = core.N_PIX_MACRO * N_SHOTS_MAX
    print("=" * 92)
    print("信号能量扫描（peak_vs_energy）：bg=0，看 peak 与 FWHM 随能量从单光子到深饱和")
    print(f"  能量网格 {NB} 档：boost {BOOSTS[0]:.4g} → {BOOSTS[-1]:.4g}（非均匀采样，坐标轴仍线性）")
    print(f"  每档 {a.n_mc:,} 次 MC；一次 {N_SHOTS_MAX} 发仿真前缀和白拿 N={list(N_SHOTS_LIST)}")
    print(f"  宏像元 {core.N_PIX_MACRO} SPAD（Σf_pix={F_PIX.sum():.4f}），"
          f"n_tr 最大 {n_tr_max}；统计窗 bin [{core.I_STAT0}, {core.I_STAT1})")
    print(f"  FWHM：逐次实现测 hist_add 波形，peak ≥ {FWHM_MIN_PEAK} 才计入")
    print("=" * 92)

    R = _load(a.n_mc)
    todo = [i for i in range(NB) if not R["done"][i]]
    if not todo:
        print("全部能量档均已完成")
        _save(R, CACHE)
        return

    print(f"待算 {len(todo)} 档；workers={a.workers}\n")
    t0 = time.time()
    jobs = [(i, float(BOOSTS[i]), int(a.n_mc), 90210 + 7919 * i) for i in todo]
    fin = 0
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        futs = {ex.submit(_one_boost, j): j[0] for j in jobs}
        for fu in as_completed(futs):
            ib, res = fu.result()
            _apply(R, ib, res)
            fin += 1
            n4 = N_SHOTS_LIST[-1]
            cnt = res[n4]["peak_cnt"]
            tot = max(int(cnt.sum()), 1)
            mu = float((np.arange(cnt.size) * cnt).sum()) / tot
            nv = res[n4]["fwhm_nval"]
            fw = res[n4]["fwhm_sum"] / nv if nv else float("nan")
            el = (time.time() - t0) / 60.0
            eta = el / fin * (len(todo) - fin)
            print(f"  [{fin}/{len(todo)} {100*fin/len(todo):5.1f}%] "
                  f"boost={BOOSTS[ib]:>10.4g}　N={n4}: peakμ={mu:6.2f}/{core.N_PIX_MACRO*n4} "
                  f"FWHM={fw:5.2f} ns（有效 {100*nv/max(tot,1):5.1f}%）"
                  f"　已用 {el:.1f} min，预计剩 {eta:.1f} min")
            if fin % CHECKPOINT_EVERY == 0:
                _save(R, CKPT)
                print(f"    …检查点已写入 {CKPT}")

    _save(R, CACHE)
    if os.path.exists(CKPT):
        os.remove(CKPT)
    print(f"\n[能量扫描完成] → {CACHE}，{(time.time()-t0)/60:.1f} min")


if __name__ == "__main__":
    main()
