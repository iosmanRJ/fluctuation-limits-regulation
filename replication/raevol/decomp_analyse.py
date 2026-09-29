#!/usr/bin/env python3
"""EXPLORATORY (not pre-registered): how much of R-Aevol's regulatory strength is
exerted by the signal proteins, and do the pre-registered contrasts hold for the
genome-encoded ("internal") part alone?

Input: paper/data/decomp_runs.csv (decomp_extract.py; continuations of every run
from its saved generation-30,000 population to 32,000, instrumented statistics)
and paper/data/raevol_runs.csv (the original analysis window, 25,001-30,000).

1. Validation: the continuation should reproduce the original cell means of S.
2. Decomposition of S into signal-derived and internal parts, per cell.
3. The P1-P3 contrasts and cross-cue comparisons on internal strength, with the
   same two-sided permutation test (100,000 permutations) as the pre-registered
   analysis. No verdicts: these are exploratory.

Usage: decomp_analyse.py <paper/data dir>
"""
import csv, math, os, random, sys
from statistics import mean, stdev

D = sys.argv[1]
C = list(csv.DictReader(open(os.path.join(D, "decomp_runs.csv"))))
O = list(csv.DictReader(open(os.path.join(D, "raevol_runs.csv"))))
for r in C + O:
    r["rate"] = round(float(r["rate"]), 4); r["pop"] = int(r["pop"])


def cell(rows, exp, arm, rate, pop, k):
    return [float(r[k]) for r in rows if r["experiment"] == exp and r["arm"] == arm
            and r["rate"] == rate and r["pop"] == pop]


def perm(a, b, n=100_000, seed=20260929):
    rng = random.Random(seed); obs = abs(mean(a) - mean(b)); z = list(a) + list(b); h = 0
    for _ in range(n):
        rng.shuffle(z)
        if abs(mean(z[:len(a)]) - mean(z[len(a):])) >= obs - 1e-12:
            h += 1
    return (h + 1) / (n + 1)


def d(a, b):
    sp = math.sqrt(((len(a) - 1) * stdev(a) ** 2 + (len(b) - 1) * stdev(b) ** 2) / (len(a) + len(b) - 2))
    return (mean(a) - mean(b)) / sp if sp > 0 else float("nan")


CELLS = [("v2", a, q, 1024) for a in ("fluc_info", "froz_info", "fluc_yoked", "fluc_none") for q in (0.05, 0.1, 0.2)] \
      + [("v2", "froz_yoked", 0.1, 1024), ("v2", "froz_none", 0.1, 1024)] \
      + [("v2", a, 0.1, 512) for a in ("fluc_info", "fluc_yoked", "froz_info")] \
      + [("v1", a, 0.1, 1024) for a in ("fluc_info", "froz_info", "fluc_yoked", "froz_yoked")]

print("1. VALIDATION: continuation (30,001-32,000) against the original window (25,001-30,000)")
print("%-4s %-11s %5s %5s  %8s %8s %6s   %8s %8s" % ("exp", "arm", "q", "pop", "S orig", "S cont", "ratio", "genes o", "genes c"))
pairs = []
for exp, arm, q, pop in CELLS:
    so, sc = cell(O, exp, arm, q, pop, "S"), cell(C, exp, arm, q, pop, "S")
    go, gc = cell(O, exp, arm, q, pop, "genes"), cell(C, exp, arm, q, pop, "genes")
    if not so or not sc:
        continue
    pairs.append((mean(so), mean(sc)))
    print("%-4s %-11s %5.2f %5d  %8.2f %8.2f %6.2f   %8.2f %8.2f" % (exp, arm, q, pop, mean(so), mean(sc), mean(sc) / mean(so) if mean(so) else float("nan"), mean(go), mean(gc)))
x, y = [p[0] for p in pairs], [p[1] for p in pairs]
mx, my = mean(x), mean(y)
r = sum((a - mx) * (b - my) for a, b in zip(x, y)) / math.sqrt(sum((a - mx) ** 2 for a in x) * sum((b - my) ** 2 for b in y))
print("correlation of cell means, original vs continuation: r = %.3f over %d cells" % (r, len(pairs)))

print("\n2. DECOMPOSITION of S in the continuation (cell means)")
print("%-4s %-11s %5s %5s  %7s %7s %7s %6s   %7s %7s %7s" % ("exp", "arm", "q", "pop", "S", "S_sig", "S_int", "sig%", "links", "l_sig", "l_int"))
for exp, arm, q, pop in CELLS:
    S = cell(C, exp, arm, q, pop, "S")
    if not S:
        continue
    g = lambda k: mean(cell(C, exp, arm, q, pop, k))
    print("%-4s %-11s %5.2f %5d  %7.2f %7.2f %7.2f %5.0f%%   %7.1f %7.1f %7.1f" % (
        exp, arm, q, pop, g("S"), g("S_sig"), g("S_int"), 100 * g("S_sig") / g("S") if g("S") else 0,
        g("links"), g("links_sig"), g("links_int")))
neg = [r for r in C if float(r["min_S_int"]) < -1e-3]
print("runs with a negative internal strength in any generation (rounding check): %d" % len(neg))

print("\n3. CONTRASTS ON INTERNAL STRENGTH (exploratory; two-sided permutation, 100,000)")
print("%-44s %5s  %8s %8s %6s %9s %6s   | same contrast on total S in the continuation" % ("contrast", "q", "mean A", "mean B", "ratio", "p", "d"))


def contrast(label, A, B, k="S_int"):
    a, b = cell(C, *A, k), cell(C, *B, k)
    ta, tb = cell(C, *A, "S"), cell(C, *B, "S")
    print("%-44s %5.2f  %8.2f %8.2f %6.2f %9.5f %6.2f   | %6.2f vs %6.2f, ratio %.2f, p=%.5f" % (
        label, A[2], mean(a), mean(b), mean(a) / mean(b) if mean(b) else float("nan"), perm(a, b), d(a, b),
        mean(ta), mean(tb), mean(ta) / mean(tb) if mean(tb) else float("nan"), perm(ta, tb)))


for q in (0.05, 0.1, 0.2):
    contrast("P1-int  fluc_info vs froz_info", ("v2", "fluc_info", q, 1024), ("v2", "froz_info", q, 1024))
for q in (0.05, 0.1, 0.2):
    contrast("P2-int  fluc_info vs fluc_yoked", ("v2", "fluc_info", q, 1024), ("v2", "fluc_yoked", q, 1024))
for q in (0.05, 0.1, 0.2):
    contrast("P3      fluc_none vs froz_none (no signals)", ("v2", "fluc_none", q, 1024), ("v2", "froz_none", 0.1, 1024))
for q in (0.05, 0.1, 0.2):
    contrast("cross   fluc_yoked vs fluc_none", ("v2", "fluc_yoked", q, 1024), ("v2", "fluc_none", q, 1024))
contrast("cross   froz_info vs froz_none", ("v2", "froz_info", 0.1, 1024), ("v2", "froz_none", 0.1, 1024))
contrast("cross   froz_yoked vs froz_none", ("v2", "froz_yoked", 0.1, 1024), ("v2", "froz_none", 0.1, 1024))
contrast("frozen  froz_info vs froz_yoked", ("v2", "froz_info", 0.1, 1024), ("v2", "froz_yoked", 0.1, 1024))
contrast("pop512  fluc_info vs froz_info", ("v2", "fluc_info", 0.1, 512), ("v2", "froz_info", 0.1, 512))
contrast("pop512  fluc_info vs fluc_yoked", ("v2", "fluc_info", 0.1, 512), ("v2", "fluc_yoked", 0.1, 512))
contrast("v1      fluc_info vs froz_info", ("v1", "fluc_info", 0.1, 1024), ("v1", "froz_info", 0.1, 1024))
contrast("v1      fluc_info vs fluc_yoked", ("v1", "fluc_info", 0.1, 1024), ("v1", "fluc_yoked", 0.1, 1024))

print("\n   internal LINK counts for the P1 contrast (counts vs strength, internal part only)")
for q in (0.05, 0.1, 0.2):
    contrast("links_int fluc_info vs froz_info", ("v2", "fluc_info", q, 1024), ("v2", "froz_info", q, 1024), k="links_int")
