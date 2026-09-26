#!/usr/bin/env python3
"""
Pre-registered analysis for the R-Aevol sweep (v2). Pre-registration: PAPER.md
section 11. Written and committed BEFORE any v2 data existed; do not change the
measures, windows, tests or decision rules after the data arrive. Anything not
specified here is exploratory and must be reported as such.

UNIT OF ANALYSIS: one run. Each run is reduced to means over generations
25,001-30,000 of stats/stats_simd_best.csv. Runs that did not reach generation
30,000 are excluded and listed.

MEASURES (per generation, then averaged over the window)
  S      total regulatory strength = nb_link * avg_influence, and 0 when
         nb_link = 0. (avg_influence is sum/count per link -- Stats_7.cpp --
         and is NaN when there are no links. The exploratory v1 script DROPPED
         those generations, which biases strength upward in sparse arms; a
         network with no links has zero strength, not missing strength.)
  genes  nb_functional_genes
  links  nb_link
  perlink avg_influence, over generations with nb_link > 0 only
  err    metabolic_error

PRIMARY PREDICTIONS, tested at each switch rate r in {0.05, 0.10, 0.20}, pop 1024
  P1  S(fluc_info, r)  >  S(froz_info, r)    trackable fluctuation strengthens regulation
  P2  S(fluc_info, r)  >  S(fluc_yoked, r)   information, with stimuli matched
  P3  S(fluc_none, r)  <  S(froz_none)       untrackable fluctuation weakens regulation
  P4  genes(froz_info, r) > genes(fluc_info, r)   fluctuation shrinks the genome

  froz_none has identical targets and no signals, so it is rate-independent and
  is run once (at r = 0.10); it is the P3 baseline at every rate.

TEST: two-sided permutation test on the difference of means, 100,000
permutations, seed fixed below. Holm correction across the 12 primary tests
(4 predictions x 3 rates).

MANIPULATION CHECK (M), per rate: err(fluc_info, r) < err(fluc_yoked, r) with
unadjusted p < 0.05 -- i.e. the organisms demonstrably use the informative cue.
P1 and P2 are only EVALUABLE at rates where M holds; where the cue is not used,
they are not tests of the hypothesis. P3 and P4 are evaluable at every rate.

VERDICT per prediction, over its evaluable rates:
  SUPPORTED     predicted sign at EVERY evaluable rate, and Holm-adjusted
                p < 0.05 at >= 2 evaluable rates
  REFUTED       opposite sign with unadjusted p < 0.05 at ANY evaluable rate,
                or Holm-adjusted p >= 0.05 at EVERY evaluable rate
  INCONCLUSIVE  anything else
  (P1/P2 with fewer than 2 evaluable rates: INCONCLUSIVE -- the design failed
   to produce cue use, which is not evidence against the hypothesis.)

SECONDARY (reported, no verdict): per-link strength and raw link count for the
P1 contrast (does count track strength?); froz_info vs froz_yoked at r = 0.10
(expected null); the P1-P4 contrasts repeated at population 512, r = 0.10;
rank-based permutation p-values alongside every primary test.

Usage:  analyse_v2.py <root>     (root contains primary/ and secondary/)
"""
import csv, math, os, random, sys
from collections import defaultdict

PERM_SEED = 20260926
NPERM = 100_000
WIN_LO, WIN_HI = 25_000, 30_000
RATES = [50, 100, 200]          # rate tags (per-mille) as they appear in run names


def summarise(d):
    p = os.path.join(d, "stats", "stats_simd_best.csv")
    if not os.path.exists(p):
        return None, "no stats"
    rows, last = [], -1
    with open(p) as f:
        for r in csv.DictReader(f):
            try:
                g = int(float(r["Generation"]))
            except (ValueError, KeyError):
                continue
            last = max(last, g)
            if not (WIN_LO < g <= WIN_HI):
                continue
            nl = float(r["nb_link"])
            ai = float(r["avg_influence"]) if nl > 0 else float("nan")
            rows.append({
                "S": nl * ai if nl > 0 else 0.0,
                "genes": float(r["nb_functional_genes"]),
                "links": nl,
                "perlink": ai,
                "err": float(r["metabolic_error"]),
            })
    if last < WIN_HI:
        return None, "reached only generation %d" % last
    if not rows:
        return None, "empty window"
    out = {}
    for k in ("S", "genes", "links", "err"):
        out[k] = sum(x[k] for x in rows) / len(rows)
    pl = [x["perlink"] for x in rows if not math.isnan(x["perlink"])]
    out["perlink"] = sum(pl) / len(pl) if pl else float("nan")
    out["frac_zero_link"] = sum(1 for x in rows if x["links"] == 0) / len(rows)
    return out, None


def parse(name):
    """fluc_info_r100_w1024_s3 -> (arm, rate, pop)"""
    base = name.rsplit("_s", 1)[0]
    arm, rate, pop = base.rsplit("_", 2)
    return arm, int(rate[1:]), int(pop[1:])


def mean(v):
    return sum(v) / len(v)


def sd(v):
    m = mean(v)
    return math.sqrt(sum((x - m) ** 2 for x in v) / (len(v) - 1)) if len(v) > 1 else 0.0


def perm_p(a, b, rng, n=NPERM):
    obs = abs(mean(a) - mean(b))
    pool = list(a) + list(b); na = len(a); hits = 0
    for _ in range(n):
        rng.shuffle(pool)
        if abs(mean(pool[:na]) - mean(pool[na:])) >= obs - 1e-12:
            hits += 1
    return (hits + 1) / (n + 1)


def rank_p(a, b, rng, n=NPERM):
    pool = sorted([(x, 0) for x in a] + [(x, 1) for x in b])
    ranks = list(range(len(pool)))
    lab = [g for _, g in pool]
    return perm_p([r for r, g in zip(ranks, lab) if g == 0],
                  [r for r, g in zip(ranks, lab) if g == 1], rng, n)


def holm(ps):
    order = sorted(range(len(ps)), key=lambda i: ps[i])
    adj = [0.0] * len(ps); run = 0.0
    for k, i in enumerate(order):
        run = max(run, min(1.0, (len(ps) - k) * ps[i]))
        adj[i] = run
    return adj


def main(root):
    rng = random.Random(PERM_SEED)
    data = defaultdict(list)          # (arm, rate, pop) -> [summary]
    excluded = []
    for tier in ("primary", "secondary"):
        td = os.path.join(root, tier)
        if not os.path.isdir(td):
            continue
        for name in sorted(os.listdir(td)):
            full = os.path.join(td, name)
            if not os.path.isdir(full):
                continue
            s, why = summarise(full)
            if s is None:
                excluded.append("%s (%s)" % (name, why)); continue
            data[parse(name)].append(s)

    print("CELLS")
    for k in sorted(data):
        print("  %-12s r=%3d pop=%4d  n=%2d" % (k[0], k[1], k[2], len(data[k])))
    if excluded:
        print("EXCLUDED: " + "; ".join(excluded))

    def v(arm, rate, pop, m):
        return [x[m] for x in data.get((arm, rate, pop), [])]

    # manipulation check
    print("\nMANIPULATION CHECK  err(fluc_info) < err(fluc_yoked)")
    M = {}
    for r in RATES:
        a, b = v("fluc_info", r, 1024, "err"), v("fluc_yoked", r, 1024, "err")
        if len(a) < 2 or len(b) < 2:
            M[r] = False; print("  r=%3d  insufficient data" % r); continue
        p = perm_p(a, b, rng)
        M[r] = mean(a) < mean(b) and p < 0.05
        print("  r=%3d  %.4f vs %.4f   p=%.5f   -> %s" %
              (r, mean(a), mean(b), p, "HOLDS" if M[r] else "FAILS"))

    # primary tests
    spec = [
        ("P1", "S",     lambda r: ("fluc_info", r),  lambda r: ("froz_info", r),  +1),
        ("P2", "S",     lambda r: ("fluc_info", r),  lambda r: ("fluc_yoked", r), +1),
        ("P3", "S",     lambda r: ("fluc_none", r),  lambda r: ("froz_none", 100), -1),
        ("P4", "genes", lambda r: ("froz_info", r),  lambda r: ("fluc_info", r),  +1),
    ]
    tests = []
    for (pid, m, fa, fb, sign) in spec:
        for r in RATES:
            a = v(*fa(r), 1024, m); b = v(*fb(r), 1024, m)
            if len(a) < 2 or len(b) < 2:
                tests.append((pid, r, None)); continue
            d = mean(a) - mean(b)
            tests.append((pid, r, dict(a=a, b=b, diff=d, sign=sign,
                                       p=perm_p(a, b, rng), prank=rank_p(a, b, rng),
                                       lab=(fa(r), fb(r)))))
    ps = [t[2]["p"] for t in tests if t[2]]
    adj = iter(holm(ps))
    for t in tests:
        if t[2]:
            t[2]["padj"] = next(adj)

    print("\nPRIMARY TESTS (Holm across %d)" % len(ps))
    print("  %-3s %4s  %-26s %9s %9s %9s %9s %9s %9s" %
          ("id", "r", "contrast", "mean A", "mean B", "diff", "p", "p_holm", "p_rank"))
    for pid, r, t in tests:
        if not t:
            print("  %-3s %4d  insufficient data" % (pid, r)); continue
        (aa, ra), (ab, rb) = t["lab"]
        print("  %-3s %4d  %-26s %9.4f %9.4f %9.4f %9.5f %9.5f %9.5f" %
              (pid, r, "%s vs %s" % (aa, ab), mean(t["a"]), mean(t["b"]),
               t["diff"], t["p"], t["padj"], t["prank"]))

    print("\nVERDICTS")
    for (pid, *_rest) in spec:
        rows = [(r, t) for (q, r, t) in tests if q == pid and t]
        if pid in ("P1", "P2"):
            rows = [(r, t) for (r, t) in rows if M.get(r)]
        if pid in ("P1", "P2") and len(rows) < 2:
            print("  %s  INCONCLUSIVE  (evaluable rates: %d -- needs >= 2; a rate is "
                  "evaluable when it has data and the manipulation check holds)" % (pid, len(rows)))
            continue
        right = [t["diff"] * t["sign"] > 0 for _, t in rows]
        sig = [t["padj"] < 0.05 for _, t in rows]
        opp = [t["diff"] * t["sign"] < 0 and t["p"] < 0.05 for _, t in rows]
        if any(opp) or not any(sig):
            verdict = "REFUTED"
        elif all(right) and sum(sig) >= 2:
            verdict = "SUPPORTED"
        else:
            verdict = "INCONCLUSIVE"
        print("  %s  %-12s  evaluable rates %s, predicted sign %d/%d, Holm-significant %d"
              % (pid, verdict, [r for r, _ in rows], sum(right), len(rows), sum(sig)))

    # secondary
    print("\nSECONDARY (exploratory, no verdicts)")
    for r in RATES:
        for m in ("links", "perlink", "S"):
            a, b = v("fluc_info", r, 1024, m), v("froz_info", r, 1024, m)
            if len(a) > 1 and len(b) > 1:
                print("  r=%3d  %-7s fluc_info %.4f  froz_info %.4f  ratio %.2f  p=%.5f"
                      % (r, m, mean(a), mean(b), mean(a) / mean(b) if mean(b) else float("nan"),
                         perm_p(a, b, rng)))
    a, b = v("froz_info", 100, 1024, "S"), v("froz_yoked", 100, 1024, "S")
    if len(a) > 1 and len(b) > 1:
        print("  frozen yoked control (expected null): S %.4f vs %.4f  p=%.5f"
              % (mean(a), mean(b), perm_p(a, b, rng)))
    for pid, m, fa, fb, sign in spec:
        A, B = fa(100), fb(100)
        a, b = v(A[0], 100, 512, m), v(B[0], 100, 512 if B[0] != "froz_none" else 1024, m)
        if len(a) > 1 and len(b) > 1:
            print("  pop512 %s  %s %s vs %s: %.4f vs %.4f  p=%.5f"
                  % (pid, m, A[0], B[0], mean(a), mean(b), perm_p(a, b, rng)))
    zl = [(k, mean([x["frac_zero_link"] for x in xs])) for k, xs in sorted(data.items())]
    print("  fraction of window generations with zero links, by cell:")
    for k, z in zl:
        print("    %-12s r=%3d pop=%4d  %.3f" % (k[0], k[1], k[2], z))


if __name__ == "__main__":
    main(sys.argv[1])
