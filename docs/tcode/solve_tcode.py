# -*- coding: utf-8 -*-
"""
solve_tcode.py —— 求【最省时间】的 tcode 码表
================================================================
v22 约束（比 v21 更强）：
  A. 散开：同一条鬼影源在 N_ACC 次 shot 的落点两两间隔 ≥ SEP
  B. 避真：同 kick 串扰的每次落点都离真峰 ≥ SEP
           （真峰相对码差 = 0，即要求 |c[a][k] - c[b][k]| ≥ SEP）

用法：
    python docs/tcode/solve_tcode.py [SEP] [BUDGET]
    若 BUDGET 下无解，会自动往上找最小可行预算。

缩写：XM（XtalkMark，串扰标记）、TOF（Time of Flight，飞行时间）、
      IRF（Instrument Response Function，仪器响应函数）。
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gen_tcode_figures as G

LASER_IDS = G.LASER_IDS
KICKS_OF  = G.KICKS_OF
N_ACC     = G.N_ACC
VARS      = [(l, k) for l in LASER_IDS for k in KICKS_OF[l]]
VIDX      = {v: i for i, v in enumerate(VARS)}

CROSSTALK_MAX_GAP = 2
DK_SET            = (0, 1)


# ---------------------------------------------------------------------------
# 约束 A：散开（鬼影类内两两分开）
# ---------------------------------------------------------------------------
def build_classes(gap=CROSSTALK_MAX_GAP, dks=DK_SET):
    out = []
    for b in LASER_IDS:
        for a in LASER_IDS:
            if abs(a - b) > gap:
                continue
            for dk in dks:
                if a == b and dk == 0:
                    continue
                pairs = [(VIDX[(a, kb - dk)], VIDX[(b, kb)])
                         for kb in KICKS_OF[b] if (kb - dk) in KICKS_OF[a]]
                if len(pairs) >= 2:
                    kind = ("自身混叠" if a == b else
                            ("同kick串扰" if dk == 0 else "跨kick串扰"))
                    out.append((f"L{a}->L{b} Δk={dk} {kind}", pairs))
    return out


CLASSES = build_classes()
VAR2CLS = [[] for _ in VARS]
for ci, (_, prs) in enumerate(CLASSES):
    for ia, ib in prs:
        VAR2CLS[ia].append(ci)
        VAR2CLS[ib].append(ci)
VAR2CLS = [sorted(set(c)) for c in VAR2CLS]


# ---------------------------------------------------------------------------
# 约束 B：避真（同 kick 串扰不得落在真峰附近）
#   真峰码差 = 0；同 kick 鬼影码差 = c[a][k]-c[b][k]
#   要求每个共同 kick 上都有 |c[a]-c[b]| ≥ SEP
#   （跨 kick / 自身混叠相对真峰还要再减 T_kick≈2200ns，码预算内不可能撞真）
# ---------------------------------------------------------------------------
def build_avoid_true(gap=CROSSTALK_MAX_GAP):
    """返回 [(ia, ib), ...]，每对必须 |u[ia]-u[ib]| ≥ SEP。"""
    pairs = []
    seen = set()
    for b in LASER_IDS:
        for a in LASER_IDS:
            if a == b or abs(a - b) > gap:
                continue
            for k in KICKS_OF[b]:
                if k not in KICKS_OF[a]:
                    continue
                ia, ib = VIDX[(a, k)], VIDX[(b, k)]
                key = (min(ia, ib), max(ia, ib))
                if key not in seen:
                    seen.add(key)
                    pairs.append((ia, ib))
    return pairs


AVOID_PAIRS = build_avoid_true()
VAR2AVOID = [[] for _ in VARS]
for pi, (ia, ib) in enumerate(AVOID_PAIRS):
    VAR2AVOID[ia].append(pi)
    VAR2AVOID[ib].append(pi)


SEP = 12


def set_sep(w):
    global SEP
    SEP = w


def class_cost(u, ci):
    prs = CLASSES[ci][1]
    d = sorted(int(u[ia]) - int(u[ib]) for ia, ib in prs)
    return sum(1 for i in range(len(d) - 1) if d[i + 1] - d[i] < SEP)


def avoid_cost_one(u, pi):
    ia, ib = AVOID_PAIRS[pi]
    return int(abs(int(u[ia]) - int(u[ib])) < SEP)


def total_cost(u):
    return (sum(class_cost(u, ci) for ci in range(len(CLASSES)))
            + sum(avoid_cost_one(u, pi) for pi in range(len(AVOID_PAIRS))))


def _class_cost_with(u, ci, vi, val):
    """把 u[vi] 临时换成 val 时，第 ci 组的散开代价（不改 u）。"""
    prs = CLASSES[ci][1]
    d = []
    for ia, ib in prs:
        a = val if ia == vi else int(u[ia])
        b = val if ib == vi else int(u[ib])
        d.append(a - b)
    d.sort()
    return sum(1 for i in range(len(d) - 1) if d[i + 1] - d[i] < SEP)


def solve(budget, rng, max_steps=80000, plateau=1500):
    """min-conflicts：散开 + 避真。
       plateau：连续这么多步代价不降就放弃（短重启友好）。"""
    u = np.asarray(rng.integers(0, budget + 1, len(VARS)), dtype=np.int32)
    cost = total_cost(u)
    best_cost, stale = cost, 0
    vals = np.arange(budget + 1, dtype=np.int32)

    for _ in range(max_steps):
        if cost == 0:
            return u
        # 只扫一遍建坏约束列表（比每步对每个候选重算全表便宜）
        bad_c = [ci for ci in range(len(CLASSES)) if class_cost(u, ci) > 0]
        bad_a = [pi for pi in range(len(AVOID_PAIRS)) if avoid_cost_one(u, pi)]
        if not bad_c and not bad_a:
            return u

        if bad_a and (not bad_c or rng.random() < 0.45):
            pi = int(bad_a[rng.integers(0, len(bad_a))])
            ia, ib = AVOID_PAIRS[pi]
            vi = int(ia if rng.random() < 0.5 else ib)
        else:
            prs = CLASSES[int(rng.choice(bad_c))][1]
            vi = int(prs[int(rng.integers(0, len(prs)))][int(rng.integers(0, 2))])

        touched_c = VAR2CLS[vi]
        touched_a = VAR2AVOID[vi]
        base = (cost
                - sum(class_cost(u, ci) for ci in touched_c)
                - sum(avoid_cost_one(u, pi) for pi in touched_a))

        # 避真项可向量化：对每个候选 val，涉及 vi 的 avoid 对
        avoid_extra = np.zeros(budget + 1, dtype=np.int32)
        for pi in touched_a:
            ia, ib = AVOID_PAIRS[pi]
            if ia == vi:
                other = int(u[ib])
                avoid_extra += (np.abs(vals - other) < SEP).astype(np.int32)
            else:
                other = int(u[ia])
                avoid_extra += (np.abs(other - vals) < SEP).astype(np.int32)

        best, bestc = [], 1 << 30
        old = int(u[vi])
        for val in range(budget + 1):
            c = (base
                 + sum(_class_cost_with(u, ci, vi, val) for ci in touched_c)
                 + int(avoid_extra[val]))
            if c < bestc:
                best, bestc = [val], c
            elif c == bestc:
                best.append(val)

        pick = old if (bestc == cost and rng.random() < 0.2) else int(best[rng.integers(0, len(best))])
        u[vi] = pick
        cost = (base
                + sum(class_cost(u, ci) for ci in touched_c)
                + sum(avoid_cost_one(u, pi) for pi in touched_a))

        if cost < best_cost:
            best_cost, stale = cost, 0
        else:
            stale += 1
            if stale >= plateau:
                return None
    return None


def make_code_fn(u, w=1):
    tbl = {v: int(u[VIDX[v]]) * w for v in VARS}
    return lambda l, k: tbl.get((l, k), 0)


def evaluate(code_fn, ratio, step=10.0):
    fr = G.build_firings(code_fn)
    t = dict(gb=0, ga=0, tb=0, kill=0, res_same=0, res_cross=0, res_mix=0)
    for D in np.arange(5.0, 601.0, step):
        r = G.simulate(D, fr, ratio=ratio)
        for k in t:
            t[k] += r[k]
    return t


def check_avoid_true(code_fn, sep=None):
    """返回 (是否全过, 违规列表)。"""
    sep = SEP if sep is None else sep
    fails = []
    for a in LASER_IDS:
        for b in LASER_IDS:
            if a >= b or abs(a - b) > CROSSTALK_MAX_GAP:
                continue
            for k in KICKS_OF[a]:
                if k not in KICKS_OF[b]:
                    continue
                d = abs(code_fn(a, k) - code_fn(b, k))
                if d < sep:
                    fails.append((a, b, k, d))
    return len(fails) == 0, fails


def dump_table(u, path, sep, budget, note=""):
    tbl = {f"{l},{k}": int(u[VIDX[(l, k)]]) for (l, k) in VARS}
    with open(path, "w", encoding="utf-8") as f:
        f.write('# -*- coding: utf-8 -*\n')
        f.write('"""本文件由 docs/tcode/solve_tcode.py 自动生成，请勿手改。\n\n')
        f.write(f'   约束：散开 + 避真峰（|码差|≥SEP）\n')
        f.write(f'   最小落点间隔 / 避真间隔 = {sep} ns\n')
        f.write(f'   码预算                   = {budget} ns\n')
        f.write(f'   实际用到的最大码值       = {max(tbl.values())} ns\n')
        f.write(f'   {note}\n"""\n')
        f.write(f'TCODE_SEP_NS    = {sep}\n')
        f.write(f'TCODE_BUDGET_NS = {budget}\n')
        f.write('# (laser_id, kick) -> tx_trig_dly [ns]\n')
        f.write('TCODE_TABLE = {\n')
        for l in LASER_IDS:
            items = ", ".join(f"({l},{k}): {tbl[f'{l},{k}']:>2d}" for k in KICKS_OF[l])
            f.write(f"    {items},\n")
        f.write('}\n')
    print(f"  码表已写入 {path}")


if __name__ == "__main__":
    SEP_ARG = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    BUDGET0 = int(sys.argv[2]) if len(sys.argv) > 2 else 64
    OUT_NAME = sys.argv[3] if len(sys.argv) > 3 else "tcode_table_v22.py"
    rng = np.random.default_rng(1)
    set_sep(SEP_ARG)

    print(f"峰宽/避真间隔 {SEP_ARG} ns | 起始预算 {BUDGET0} ns")
    print(f"  散开约束组 {len(CLASSES)} 个 | 避真约束对 {len(AVOID_PAIRS)} 个")

    # 只从用户给的预算往上找（不先试更小的——更小预算失败极慢）
    search_order = [BUDGET0] + [b for b in
                                (40, 48, 56, 64, 72, 80, 96, 112, 128, 144, 160, 176, 192)
                                if b > BUDGET0]

    feasible = None
    for B in search_order:
        n_try = 12 if B <= 96 else 8
        sols = []
        for _ in range(n_try):
            s = solve(B, rng)
            if s is not None:
                sols.append(s)
        print(f"  预算 {B:>3d} ns -> {'可行（%d/%d）' % (len(sols), n_try) if sols else '搜不到'}")
        if sols and feasible is None:
            feasible = (B, sols)
            # 找到最小可行后，若后面还有更大的就不必继续往上
            # 但若当前来自 below，继续找更小的已经按升序，可 break 于首次成功
            break
    if feasible is None:
        raise SystemExit("给定范围内无解")

    B, sols = feasible
    best = min(sols, key=lambda s: (
        evaluate(make_code_fn(s), 1.6, 20.0)["ga"],
        max(s),
    ))
    fn = make_code_fn(best)
    ok_av, fails = check_avoid_true(fn, SEP_ARG)
    print(f"\n  最小可行预算 = {B} ns | 实际 max(tx) = {max(best)} ns | 避真检查 = {'通过' if ok_av else '失败'}")
    if fails[:5]:
        print("  违规样例:", fails[:5])

    print("  实测（1~600m，步长 2m）：")
    for ratio in (1.6, 2.5):
        t = evaluate(fn, ratio, step=2.0)
        print(f"    ratio={ratio}: 鬼影 {t['gb']} -> 残留 {t['ga']}"
              f" ({t['ga']/max(t['gb'],1):.3%})  误杀 {t['kill']}/{t['tb']}")

    out_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), OUT_NAME)
    dump_table(best, out_path, SEP_ARG, B,
               note="v22：散开 + 避真峰；单光子脉宽 8ns，SEP=峰宽裕度")
    print("\n  tx_trig_dly [ns]")
    for l in LASER_IDS:
        print(f"    L{l:<3d} " + "  ".join(f"K{k:<2d}={int(best[VIDX[(l,k)]]):<3d}"
                                           for k in KICKS_OF[l]))
