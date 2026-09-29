#!/usr/bin/env python3
"""
Pre-registered analysis for the R-Aevol CONFIRMATION of the signal-link
decomposition. Pre-registration: PAPER.md section 11.19. Written and committed
BEFORE any confirmation run existed; do not change the measures, window, tests,
decision rules or interpretations after the data arrive. Anything not specified
here is exploratory.

RUNS: gen_params_confirm.py (population 1024, 12 new seeds per cell), run from
generation 0 with the instrumented binary (instrument_signal_links.py), which
appends nb_signal_link and sum_signal_influence at the END of every row of
stats/stats_simd_best.csv (read by position; the header is unchanged).

UNIT OF ANALYSIS: one run, reduced to means over generations 25,001-30,000.
Runs that did not reach generation 30,000 are excluded and listed.

MEASURES (per generation, then averaged over the window)
  S       total regulatory strength = nb_link * avg_influence, 0 when nb_link = 0
  S_sig   sum_signal_influence: the part exerted by signal proteins
  S_int   S - S_sig: the genome-encoded ("internal") part
  err     metabolic_error

PRIMARY PREDICTIONS on S_int, at each rate q in {0.10, 0.20}
  C1  S_int(fluc_none, q)  < S_int(froz_none)       no cue: fluctuation erodes the network
  C2  S_int(fluc_yoked, q) < S_int(froz_yoked, q)   uninformative cue: fluctuation erodes it
  C3  S_int(fluc_info, q)  > S_int(fluc_yoked, q)   information preserves it
  C4  S_int(fluc_info, q) / S_int(froz_info, q) is not below 0.75 (non-inferiority):
      the lower bound of the 95% bootstrap CI of the ratio of means > 0.75

TESTS: C1-C3 two-sided permutation test on the difference of means, 100,000
permutations, seed fixed below; Holm correction across the 6 tests (3 x 2 rates).
C4: bootstrap CI of the ratio of means, 10,000 resamples of runs within each arm,
seed fixed below; not part of the Holm family.

MANIPULATION CHECK (M), per rate: err(fluc_info, q) < err(fluc_yoked, q) with
unadjusted p < 0.05. C3 and C4 are EVALUABLE only at rates where M holds; C1 and
C2 at both rates.

VERDICT per prediction, over its evaluable rates
  C1-C3  SUPPORTED     predicted sign and Holm-adjusted p < 0.05 at every evaluable
                       rate, with at least 2 evaluable rates
         REFUTED       opposite sign with unadjusted p < 0.05 at any evaluable rate,
                       or Holm-adjusted p >= 0.05 at every evaluable rate
         INCONCLUSIVE  anything else, including fewer than 2 evaluable rates
  C4     SUPPORTED     lower CI bound > 0.75 at every evaluable rate (at least 2)
         REFUTED       upper CI bound < 0.75 at any evaluable rate
         INCONCLUSIVE  anything else

WHAT EACH OUTCOME MEANS FOR THE PAPER (fixed in advance)
  C1-C4 supported  -> the headline is confirmed: without information fluctuation
                      erodes the genome-encoded network, with information the
                      network keeps its frozen-environment strength.
  C2 not supported -> an uninformative signal may partly protect the network; the
                      claim is restricted to fluctuation without any cue.
  C3 refuted       -> information does not preserve the internal network; the P2
                      effect is signal wiring only, and the headline changes.
  C4 refuted       -> information does not maintain the network; fluctuation
                      erodes it regardless, and the title must change.
  C4 inconclusive  -> maintenance is reported as not established.

SECONDARY (reported, no verdicts): cell means of S, S_sig, S_int, links, genes,
err; share of S from signal proteins; the P1 contrast (fluc_info vs froz_info) on
S and on S_int at each rate; rank-based permutation p beside every primary test.

Usage:  analyse_confirm.py <root>     (root contains primary/)
"""
import csv, math, os, random, sys
from collections import defaultdict

PERM_SEED = 20260929
BOOT_SEED = 20260930
NPERM = 100_000
NBOOT = 10_000
WIN_LO, WIN_HI = 25_000, 30_000
MARGIN = 0.75
RATES = [100, 200]          # rate tags (per-mille) in run names


def mean(x):
    return sum(x) / len(x)


def parse(name):
    """fluc_info_r100_w1024_s3 -> (arm, rate_tag, pop)"""
    base = name.rsplit("_s", 1)[0]
    arm, r, w = base.rsplit("_", 2)
    return arm, int(r[1:]), int(w[1:])


def summarise(d):
    p = os.path.join(d, "stats", "stats_simd_best.csv")
    if not os.path.exists(p):
        return None, "no stats"
    with open(p) as f:
        rows = list(csv.reader(f))
    head = rows[0]; H = {k: i for i, k in enumerate(head)}; n0 = len(head)
    acc, last = [], -1
    for r in rows[1:]:
        try:
            g = int(float(r[0]))
        except (ValueError, IndexError):
            continue
        last = max(last, g)
        if not (WIN_LO < g <= WIN_HI):
            continue
        extra = [x for x in r[n0:] if x != ""]
        if len(extra) < 2:
            return None, "row without instrumented fields at generation %d" % g
        nl = float(r[H["nb_link"]])
        S = nl * float(r[H["avg_influence"]]) if nl > 0 else 0.0
        ssig = float(extra[-1])
        acc.append({"S": S, "S_sig": ssig, "S_int": S - ssig, "links": nl,
                    "links_sig": float(extra[-2]), "genes": float(r[H["nb_functional_genes"]]),
                    "err": float(r[H["metabolic_error"]])})
    if last < WIN_HI:
        return None, "reached only generation %d" % last
    if not acc:
        return None, "empty window"
    return {k: mean([a[k] for a in acc]) for k in acc[0]}, None


def perm_p(a, b, rng, n=NPERM):
    obs = abs(mean(a) - mean(b)); pool = list(a) + list(b); na = len(a); hits = 0
    for _ in range(n):
        rng.shuffle(pool)
        if abs(mean(pool[:na]) - mean(pool[na:])) >= obs - 1e-12:
            hits += 1
    return (hits + 1) / (n + 1)


def rank_p(a, b, rng, n=NPERM):
    pool = list(a) + list(b)
    order = sorted(range(len(pool)), key=lambda i: pool[i])
    ranks = [0.0] * len(pool); i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and pool[order[j + 1]] == pool[order[i]]:
            j += 1
        for k in range(i, j + 1):
            ranks[order[k]] = (i + j) / 2 + 1
        i = j + 1
    return perm_p(ranks[:len(a)], ranks[len(a):], rng, n)


def holm(ps):
    idx = sorted(range(len(ps)), key=lambda i: ps[i]); out = [0.0] * len(ps); run = 0.0
    for k, i in enumerate(idx):
        run = max(run, min(1.0, (len(ps) - k) * ps[i])); out[i] = run
    return out


def boot_ratio(a, b, rng, n=NBOOT):
    rs = []
    for _ in range(n):
        ma = mean([a[rng.randrange(len(a))] for _ in a]); mb = mean([b[rng.randrange(len(b))] for _ in b])
        rs.append(ma / mb if mb > 0 else float("inf"))
    rs.sort()
    return mean(a) / mean(b), rs[int(0.025 * n)], rs[int(0.975 * n) - 1]


def main(root):
    cells = defaultdict(list); excluded = []
    pdir = os.path.join(root, "primary")
    for name in sorted(os.listdir(pdir)):
        d = os.path.join(pdir, name)
        if not os.path.isdir(d):
            continue
        s, why = summarise(d)
        if s is None:
            excluded.append((name, why)); continue
        arm, r, w = parse(name)
        cells[(arm, r)].append(s)
    print("CELLS")
    for k in sorted(cells):
        c = cells[k]
        print("  %-11s r=%3d n=%2d  S %6.2f  S_sig %6.2f  S_int %6.2f  links %6.1f  genes %6.2f  err %.4f" % (
            k[0], k[1], len(c), mean([x["S"] for x in c]), mean([x["S_sig"] for x in c]), mean([x["S_int"] for x in c]),
            mean([x["links"] for x in c]), mean([x["genes"] for x in c]), mean([x["err"] for x in c])))
    if excluded:
        print("EXCLUDED:", excluded)
    v = lambda arm, r, k: [x[k] for x in cells[(arm, r)]]
    rng = random.Random(PERM_SEED)

    print("\nMANIPULATION CHECK  err(fluc_info) < err(fluc_yoked)")
    M = {}
    for r in RATES:
        a, b = v("fluc_info", r, "err"), v("fluc_yoked", r, "err")
        p = perm_p(a, b, rng); M[r] = (mean(a) < mean(b)) and p < 0.05
        print("  r=%3d  %.4f vs %.4f  p=%.5f  -> %s" % (r, mean(a), mean(b), p, "HOLDS" if M[r] else "FAILS"))

    spec = {"C1": (("fluc_none", None), ("froz_none", 100), -1),
            "C2": (("fluc_yoked", None), ("froz_yoked", None), -1),
            "C3": (("fluc_info", None), ("fluc_yoked", None), +1)}
    tests = []
    for cid, ((aa, ar), (ba, br), sign) in spec.items():
        for r in RATES:
            a, b = v(aa, ar or r, "S_int"), v(ba, br or r, "S_int")
            tests.append({"id": cid, "r": r, "A": "%s r%d" % (aa, ar or r), "B": "%s r%d" % (ba, br or r),
                          "ma": mean(a), "mb": mean(b), "sign": sign,
                          "p": perm_p(a, b, rng), "p_rank": rank_p(a, b, rng)})
    for t, ph in zip(tests, holm([t["p"] for t in tests])):
        t["p_holm"] = ph
    print("\nPRIMARY TESTS on S_int (Holm across %d)" % len(tests))
    print("  %-3s %4s  %-18s %-18s %8s %8s %9s %9s %9s" % ("id", "r", "A", "B", "mean A", "mean B", "p", "p_holm", "p_rank"))
    for t in tests:
        print("  %-3s %4d  %-18s %-18s %8.3f %8.3f %9.5f %9.5f %9.5f" % (t["id"], t["r"], t["A"], t["B"], t["ma"], t["mb"], t["p"], t["p_holm"], t["p_rank"]))

    print("\nC4  non-inferiority  S_int(fluc_info)/S_int(froz_info), margin %.2f" % MARGIN)
    brng = random.Random(BOOT_SEED); c4 = {}
    for r in RATES:
        ratio, lo, hi = boot_ratio(v("fluc_info", r, "S_int"), v("froz_info", r, "S_int"), brng)
        c4[r] = (ratio, lo, hi)
        print("  r=%3d  ratio %.3f  95%% CI [%.3f, %.3f]" % (r, ratio, lo, hi))

    print("\nVERDICTS")
    for cid in ("C1", "C2", "C3"):
        rows = [t for t in tests if t["id"] == cid and (cid in ("C1", "C2") or M[t["r"]])]
        ok_sign = all((t["ma"] - t["mb"]) * t["sign"] > 0 for t in rows)
        sig = [t for t in rows if t["p_holm"] < 0.05]
        opposite = any((t["ma"] - t["mb"]) * t["sign"] < 0 and t["p"] < 0.05 for t in rows)
        if len(rows) < 2:
            verdict = "INCONCLUSIVE"
        elif opposite or not sig:
            verdict = "REFUTED"
        elif ok_sign and len(sig) == len(rows):
            verdict = "SUPPORTED"
        else:
            verdict = "INCONCLUSIVE"
        print("  %s  %-12s evaluable rates %s, predicted sign %d/%d, Holm-significant %d" % (
            cid, verdict, [t["r"] for t in rows], sum((t["ma"] - t["mb"]) * t["sign"] > 0 for t in rows), len(rows), len(sig)))
    ev = [r for r in RATES if M[r]]
    if len(ev) < 2:
        v4 = "INCONCLUSIVE"
    elif any(c4[r][2] < MARGIN for r in ev):
        v4 = "REFUTED"
    elif all(c4[r][1] > MARGIN for r in ev):
        v4 = "SUPPORTED"
    else:
        v4 = "INCONCLUSIVE"
    print("  C4  %-12s evaluable rates %s" % (v4, ev))

    print("\nSECONDARY (no verdicts)")
    for r in RATES:
        S, Ss = mean(v("fluc_info", r, "S")), mean(v("fluc_info", r, "S_sig"))
        print("  r=%3d  signal share of S, fluc_info: %.2f" % (r, Ss / S if S else float("nan")))
        for k in ("S", "S_int"):
            a, b = v("fluc_info", r, k), v("froz_info", r, k)
            print("  r=%3d  P1 on %-5s fluc_info %.3f vs froz_info %.3f  ratio %.2f  p=%.5f" % (r, k, mean(a), mean(b), mean(a) / mean(b), perm_p(a, b, rng)))


if __name__ == "__main__":
    main(sys.argv[1])
