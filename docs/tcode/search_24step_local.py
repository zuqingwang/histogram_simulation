# -*- coding: utf-8 -*-
"""
search_24step_local.py —— 本地长时间随机搜索零残留码表

字母表：24×(0..4) = {0, 24, 48, 72, 96}
约束  ：ratio=1.5，SEP=12 ns，gap=2，预算 96 ns
验收  ：1~600 m 扫描，XM 后鬼影残留=0、真峰误杀=0

用法（在项目根目录）：

  # 跑 8 小时，8 进程，找到零残留自动写码表并退出
  python docs/tcode/search_24step_local.py --hours 8 --jobs 8

  # 一直循环多轮，直到找到零残留（推荐挂机）
  python docs/tcode/search_24step_local.py --loop --jobs 8

  # 每轮 30 分钟，共跑 10 轮（即使没找到也停）
  python docs/tcode/search_24step_local.py --minutes 30 --rounds 10 --jobs 8

  # 断点续跑（读取 checkpoint 里的 round / seed_base）
  python docs/tcode/search_24step_local.py --loop --jobs 8 --resume

找到零残留后写入：
  docs/tcode/tcode_table_v40_r1.5_L5_24step_96ns.py

进度与最优 near-miss 写入 checkpoint（默认 docs/tcode/search_24step_ckpt.json）。
"""
from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import os
import sys
import time
from datetime import datetime

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import gen_tcode_figures as G
import solve_tcode as S

SEP, GAP, RATIO, BUDGET = 12, 2, 1.5, 96
ALPH = np.array([0, 24, 48, 72, 96], dtype=np.int32)
OUT_TABLE = os.path.join(HERE, "tcode_table_v40_r1.5_L5_24step_96ns.py")
DEFAULT_CKPT = os.path.join(HERE, "search_24step_ckpt.json")
SEED_SPAN = 1_000_000_000  # 每个 worker 独占 1e9 的 seed 段


def rebuild_module(mod):
    mod.CROSSTALK_MAX_GAP = GAP
    mod.set_sep(SEP)
    mod.CLASSES = mod.build_classes(gap=GAP)
    mod.AVOID_PAIRS = mod.build_avoid_true(gap=GAP)
    mod.VAR2CLS = [[] for _ in mod.VARS]
    for ci, (_, prs) in enumerate(mod.CLASSES):
        for ia, ib in prs:
            mod.VAR2CLS[ia].append(ci)
            mod.VAR2CLS[ib].append(ci)
    mod.VAR2CLS = [sorted(set(c)) for c in mod.VAR2CLS]
    mod.VAR2AVOID = [[] for _ in mod.VARS]
    for pi, (ia, ib) in enumerate(mod.AVOID_PAIRS):
        mod.VAR2AVOID[ia].append(pi)
        mod.VAR2AVOID[ib].append(pi)


def solve_alphabet(sloc, alph, rng, max_steps=40000, plateau=2000):
    u = alph[rng.integers(0, len(alph), len(sloc.VARS))].astype(np.int32)
    cost = sloc.total_cost(u)
    best_cost, stale = cost, 0
    for _ in range(max_steps):
        if cost == 0:
            return u
        bad_c = [ci for ci in range(len(sloc.CLASSES)) if sloc.class_cost(u, ci) > 0]
        bad_a = [pi for pi in range(len(sloc.AVOID_PAIRS)) if sloc.avoid_cost_one(u, pi)]
        if not bad_c and not bad_a:
            return u
        if bad_a and (not bad_c or rng.random() < 0.45):
            pi = int(bad_a[rng.integers(0, len(bad_a))])
            ia, ib = sloc.AVOID_PAIRS[pi]
            vi = int(ia if rng.random() < 0.5 else ib)
        else:
            prs = sloc.CLASSES[int(rng.choice(bad_c))][1]
            vi = int(prs[int(rng.integers(0, len(prs)))][int(rng.integers(0, 2))])
        touched_c = sloc.VAR2CLS[vi]
        touched_a = sloc.VAR2AVOID[vi]
        base = (cost
                - sum(sloc.class_cost(u, ci) for ci in touched_c)
                - sum(sloc.avoid_cost_one(u, pi) for pi in touched_a))
        best, bestc = [], 1 << 30
        old = int(u[vi])
        for val in alph:
            av = 0
            for pi in touched_a:
                ia, ib = sloc.AVOID_PAIRS[pi]
                other = int(u[ib if ia == vi else ia])
                av += int(abs(int(val) - other) < SEP)
            c = (base
                 + sum(sloc._class_cost_with(u, ci, vi, int(val)) for ci in touched_c)
                 + av)
            if c < bestc:
                best, bestc = [int(val)], c
            elif c == bestc:
                best.append(int(val))
        pick = old if (bestc == cost and rng.random() < 0.2) else int(best[rng.integers(0, len(best))])
        u[vi] = pick
        cost = (base
                + sum(sloc.class_cost(u, ci) for ci in touched_c)
                + sum(sloc.avoid_cost_one(u, pi) for pi in touched_a))
        if cost < best_cost:
            best_cost, stale = cost, 0
        else:
            stale += 1
            if stale >= plateau:
                return None
    return None


def residual(gloc, sloc, u, step):
    fn = sloc.make_code_fn(u)
    fr = gloc.build_firings(fn)
    ga = gb = kill = tb = 0
    for d in np.arange(5.0, 601.0, step):
        r = gloc.simulate(d, fr, max_gap=GAP, ratio=RATIO)
        ga += r["ga"]
        gb += r["gb"]
        kill += r["kill"]
        tb += r["tb"]
        if ga or kill:
            return ga, gb, kill, tb, False
    return ga, gb, kill, tb, True


def worker(payload):
    """子进程：从 seed0 起顺序递增 seed，直到超时。"""
    worker_id, seed0, seconds, log_every = payload
    sys.path.insert(0, HERE)
    import gen_tcode_figures as gloc
    import solve_tcode as sloc

    rebuild_module(sloc)
    alph = np.array([0, 24, 48, 72, 96], dtype=np.int32)

    t0 = time.time()
    n_feas = n_try = 0
    best_ga = None
    best_u = None
    best_seed = None
    seed = seed0
    seen_feas = set()

    while time.time() - t0 < seconds:
        n_try += 1
        cur_seed = seed
        seed += 1
        u = solve_alphabet(sloc, alph, np.random.default_rng(cur_seed))
        if u is None:
            continue
        key = tuple(int(x) for x in u)
        if key in seen_feas:
            continue
        seen_feas.add(key)
        n_feas += 1

        ga, gb, kill, tb, clean = residual(gloc, sloc, u, 10.0)
        if best_ga is None or ga < best_ga or (ga == best_ga and kill < 1):
            best_ga = ga
            best_u = [int(x) for x in u]
            best_seed = cur_seed

        if not clean:
            if n_feas % log_every == 0:
                print(f"  [W{worker_id:02d}] seed={cur_seed} 可行={n_feas}/{n_try} "
                      f"唯一={len(seen_feas)} 最好残留={best_ga} "
                      f"elapsed={time.time()-t0:.0f}s", flush=True)
            continue

        ga, gb, kill, tb, clean = residual(gloc, sloc, u, 2.0)
        if clean:
            return {
                "ok": True,
                "u": [int(x) for x in u],
                "seed": cur_seed,
                "worker_id": worker_id,
                "feas": n_feas,
                "tries": n_try,
                "unique": len(seen_feas),
                "secs": time.time() - t0,
                "gb": gb,
                "tb": tb,
            }

    return {
        "ok": False,
        "worker_id": worker_id,
        "seed_end": seed - 1,
        "feas": n_feas,
        "tries": n_try,
        "unique": len(seen_feas),
        "secs": time.time() - t0,
        "best_ga": best_ga,
        "best_u": best_u,
        "best_seed": best_seed,
    }


def load_ckpt(path):
    if not os.path.isfile(path):
        return None
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def save_ckpt(path, data):
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(data, f, ensure_ascii=False, indent=2)
    os.replace(tmp, path)


def write_table(u, seed, note_extra=""):
    note = f"v40 alphabet 24*(0..4); ratio={RATIO} zero; seed={seed}"
    if note_extra:
        note += f"; {note_extra}"
    S.dump_table(u, OUT_TABLE, SEP, BUDGET, note=note)
    with open(OUT_TABLE, "a", encoding="utf-8") as f:
        f.write("\nTCODE_ALPHABET=[0, 24, 48, 72, 96]\nTCODE_N_LEVELS=5\n")


def print_table_preview(u):
    for l in S.LASER_IDS:
        print(f"  L{l:<2d} " + " ".join(
            f"K{k}={int(u[S.VIDX[(l, k)]]):>2d}" for k in S.KICKS_OF[l]))


def run_round(round_idx, seed_base, jobs, seconds, log_every):
    seeds = [seed_base + k * SEED_SPAN for k in range(jobs)]
    print(f"\n--- 第 {round_idx} 轮  seed_base={seed_base}  "
          f"worker起点={seeds}  每进程{seconds:.0f}s ---", flush=True)

    payloads = [(k, seeds[k], seconds, log_every) for k in range(jobs)]
    t0 = time.time()
    found = None
    round_feas = 0
    round_best_ga = None
    round_best_u = None
    round_best_seed = None

    with mp.Pool(jobs) as pool:
        for res in pool.imap_unordered(worker, payloads):
            round_feas += res["feas"]
            if res["ok"]:
                found = res
                pool.terminate()
                break
            ga = res.get("best_ga")
            if ga is not None and (round_best_ga is None or ga < round_best_ga):
                round_best_ga = ga
                round_best_u = res.get("best_u")
                round_best_seed = res.get("best_seed")
            print(f"  [round {round_idx}] worker={res['worker_id']:02d} "
                  f"seed→{res.get('seed_end', '?')}  "
                  f"可行={res['feas']}/{res['tries']}  "
                  f"唯一={res['unique']}  "
                  f"最好残留={res.get('best_ga')}  "
                  f"({res['secs']:.0f}s)", flush=True)

    elapsed = time.time() - t0
    return {
        "found": found,
        "round_feas": round_feas,
        "elapsed": elapsed,
        "round_best_ga": round_best_ga,
        "round_best_u": round_best_u,
        "round_best_seed": round_best_seed,
        "next_seed_base": seed_base + jobs * SEED_SPAN,
    }


def main():
    ap = argparse.ArgumentParser(description="本地随机搜索 {0,24,48,72,96} 零残留码表")
    ap.add_argument("--minutes", type=float, default=0.0, help="每轮分钟数（与 --hours 二选一）")
    ap.add_argument("--hours", type=float, default=0.0, help="每轮小时数")
    ap.add_argument("--rounds", type=int, default=1, help="轮数；--loop 时忽略")
    ap.add_argument("--loop", action="store_true", help="无限循环直到找到零残留")
    ap.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    ap.add_argument("--seed-base", type=int, default=0,
                    help="第 0 轮 worker0 起点；0=自动（时间戳）")
    ap.add_argument("--checkpoint", default=DEFAULT_CKPT, help="断点文件路径")
    ap.add_argument("--resume", action="store_true", help="从 checkpoint 续跑")
    ap.add_argument("--log-every", type=int, default=100, help="每 N 个可行解打印一次进度")
    args = ap.parse_args()

    if args.hours > 0:
        seconds_per_round = args.hours * 3600.0
    elif args.minutes > 0:
        seconds_per_round = args.minutes * 60.0
    else:
        seconds_per_round = 30.0 * 60.0  # 默认每轮 30 分钟

    rebuild_module(S)
    print("=" * 64)
    print("search_24step_local — 零残留随机搜索")
    print(f"  字母表 = {list(ALPH)}  ratio={RATIO}  SEP={SEP}  gap={GAP}  预算={BUDGET}ns")
    print(f"  散开类={len(S.CLASSES)}  避真={len(S.AVOID_PAIRS)}")
    print(f"  进程={args.jobs}  每轮={seconds_per_round/60:.1f}min  "
          f"模式={'loop' if args.loop else f'{args.rounds}轮'}")
    print(f"  checkpoint = {args.checkpoint}")
    print(f"  输出码表   = {OUT_TABLE}")
    print("=" * 64, flush=True)

    ckpt = load_ckpt(args.checkpoint) if args.resume else None
    if ckpt:
        round_start = ckpt.get("round", 0) + 1
        seed_base = ckpt.get("next_seed_base", ckpt.get("seed_base", 0))
        total_feas = ckpt.get("total_feas", 0)
        global_best_ga = ckpt.get("best_ga")
        global_best_u = ckpt.get("best_u")
        global_best_seed = ckpt.get("best_seed")
        print(f"  续跑：从第 {round_start} 轮开始，累计可行≈{total_feas}，"
              f"历史最好残留={global_best_ga}", flush=True)
    else:
        round_start = 0
        seed_base = args.seed_base if args.seed_base else int(time.time()) % 1_000_000_000
        total_feas = 0
        global_best_ga = None
        global_best_u = None
        global_best_seed = None

    t_global = time.time()
    round_idx = round_start

    while True:
        if not args.loop and round_idx >= round_start + args.rounds:
            break

        result = run_round(round_idx, seed_base, args.jobs, seconds_per_round, args.log_every)
        total_feas += result["round_feas"]

        if result["round_best_ga"] is not None:
            if global_best_ga is None or result["round_best_ga"] < global_best_ga:
                global_best_ga = result["round_best_ga"]
                global_best_u = result["round_best_u"]
                global_best_seed = result["round_best_seed"]

        ckpt_data = {
            "updated": datetime.now().isoformat(timespec="seconds"),
            "round": round_idx,
            "seed_base": seed_base,
            "next_seed_base": result["next_seed_base"],
            "total_feas": total_feas,
            "best_ga": global_best_ga,
            "best_seed": global_best_seed,
            "best_u": global_best_u,
            "elapsed_total_s": round(time.time() - t_global, 1),
        }
        save_ckpt(args.checkpoint, ckpt_data)

        print(f"\n  [round {round_idx}] 本轮可行={result['round_feas']}  "
              f"累计≈{total_feas}  全局最好残留={global_best_ga}  "
              f"用时={result['elapsed']:.0f}s", flush=True)

        if result["found"]:
            found = result["found"]
            u = np.array(found["u"], dtype=np.int32)
            write_table(u, found["seed"],
                        note_extra=f"round={round_idx} worker={found['worker_id']}")
            print("\n" + "=" * 64)
            print(f"★ 找到零残留！round={round_idx}  worker={found['worker_id']}  "
                  f"seed={found['seed']}")
            print(f"  可行={found['feas']}/{found['tries']}  唯一={found['unique']}  "
                  f"gb={found['gb']}  tb={found['tb']}")
            print(f"  使用值={sorted(set(map(int, u)))}")
            print_table_preview(u)
            print(f"  累计可行≈{total_feas}  总用时={time.time()-t_global:.0f}s")
            print("=" * 64, flush=True)
            sys.exit(0)

        seed_base = result["next_seed_base"]
        round_idx += 1

    print("\n" + "=" * 64)
    print(f"搜索结束：累计可行≈{total_feas}  全局最好残留={global_best_ga}")
    if global_best_u is not None:
        print(f"  历史最优 seed={global_best_seed}  "
              f"（已存 checkpoint，残留={global_best_ga}）")
    print(f"  总用时={time.time()-t_global:.0f}s")
    print("  未找到零残留。可用 --loop --resume 继续跑。", flush=True)
    print("=" * 64)
    sys.exit(2)


if __name__ == "__main__":
    mp.freeze_support()
    main()
