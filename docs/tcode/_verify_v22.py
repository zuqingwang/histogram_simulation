# -*- coding: utf-8 -*-
import sys
sys.path.insert(0, ".")
import solve_tcode as S
from tcode_table_v22 import TCODE_TABLE as T, TCODE_SEP_NS as SEP, TCODE_BUDGET_NS as B

S.set_sep(SEP)
fn = lambda l, k: T[(l, k)]
ok, fails = S.check_avoid_true(fn, SEP)
spread_bad = 0
for desc, prs in S.CLASSES:
    d = sorted(fn(*S.VARS[ia]) - fn(*S.VARS[ib]) for ia, ib in prs)
    gaps = [d[i+1]-d[i] for i in range(len(d)-1)]
    if gaps and min(gaps) < SEP:
        spread_bad += 1
print(f"budget={B} max={max(T.values())} avoid={ok} spread_bad={spread_bad}")
for ratio in (1.6, 2.5):
    t = S.evaluate(fn, ratio, step=2.0)
    print(f"ratio={ratio}: {t['gb']}->{t['ga']} ({t['ga']/max(t['gb'],1):.3%}) "
          f"kill {t['kill']}/{t['tb']}")
