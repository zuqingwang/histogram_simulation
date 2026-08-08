# -*- coding: utf-8 -*-
import sys, numpy as np, time
sys.path.insert(0, '.')
import docs.tcode.solve_tcode as S
import docs.tcode.gen_tcode_figures as G

LASER_IDS, KICKS_OF = G.LASER_IDS, G.KICKS_OF

print("理论下界分析（4 个差值两两相差 >= SEP，差值范围 [-B, +B]，需 2B >= 3*SEP）")
for SEP in (8, 10, 12):
    print(f"  SEP={SEP}ns  最小预算 >= {1.5*SEP:.0f}ns")
print()


def quick_search(SEP, BUDGET, restarts=60, steps=25000):
    S.set_sep(SEP)
    rng = np.random.default_rng(42)
    t0 = time.time()
    for i in range(restarts):
        u = S.solve(BUDGET, rng, max_steps=steps)
        if u is not None:
            elapsed = time.time() - t0
            fn = S.make_code_fn(u)
            C = max(fn(l, k) for l in LASER_IDS for k in KICKS_OF[l])
            ev = S.evaluate(fn, 2.5, step=5.0)
            ev16 = S.evaluate(fn, 1.6, step=5.0)
            return True, u, C, elapsed, ev, ev16
        if time.time() - t0 > 45:
            break
    return False, None, None, time.time() - t0, None, None


print(f"{'SEP':>5} {'budget':>7} {'result':>8} {'actual_C':>10} "
      f"{'r2.5%':>8} {'r1.6%':>8} {'time':>7}")
print("-" * 60)
for SEP in (8, 10, 12):
    for B in (24, 32, 40, 48):
        ok, u, C, elapsed, ev, ev16 = quick_search(SEP, B)
        if ok:
            r25 = ev['ga'] / max(ev['gb'], 1) * 100
            r16 = ev16['ga'] / max(ev16['gb'], 1) * 100
            print(f"{SEP:>5} {B:>7} {'OK':>8} {C:>10} {r25:>7.2f}% {r16:>7.2f}% {elapsed:>6.1f}s")
        else:
            print(f"{SEP:>5} {B:>7} {'FAIL':>8} {'-':>10} {'-':>8} {'-':>8} {elapsed:>6.1f}s")
        sys.stdout.flush()
