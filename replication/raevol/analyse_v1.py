#!/usr/bin/env python3
"""
Analyse the R-Aevol replication.

Each run writes stats/stats_simd_best.csv, one row per generation, for the
best individual. We summarise each run by the MEAN over a window of
generations ending at a COMMON cutoff shared by every run, which gives one
independent number per run -- the unit of analysis is a RUN, not a
generation. Pooling generations as if independent would inflate n by ~10^4
and is exactly the mistake that produced a spurious result earlier in this
project.

Metrics
  genes    nb_functional_genes            genome SIZE
  links    nb_link                        number of regulatory edges
  density  nb_link / nb_functional_genes  edges per gene  <- the DENSITY claim
  dna      amount_of_dna
  err      metabolic_error (lower = better adapted)

Primary contrasts (pre-specified):
  C1  fluc_none  vs froz_none    does fluctuation raise size and density?
  C2  fluc_info  vs fluc_yoked   does the cue's INFORMATION CONTENT matter?

C2 is the yoked control: both arms receive the same two signal species at the
same marginal frequencies and the same per-age target- and signal-change
rates; only the contingency between signal and target differs.
"""
import csv, math, os, random, sys
from collections import defaultdict

random.seed(20260925)

KEEP = ("Generation", "amount_of_dna", "nb_functional_genes", "nb_link",
        "metabolic_error", "fitness")


def read_rows(d):
    p = os.path.join(d, "stats", "stats_simd_best.csv")
    if not os.path.exists(p):
        return []
    rows = []
    with open(p) as f:
        for r in csv.DictReader(f):
            try:
                rows.append({k: float(r[k]) for k in KEEP})
            except (ValueError, TypeError, KeyError):
                continue
    rows.sort(key=lambda x: x["Generation"])
    return rows


def summarise(rows, cutoff, window):
    """Mean of each metric over (cutoff-window, cutoff]."""
    tail = [x for x in rows if cutoff - window < x["Generation"] <= cutoff]
    if len(tail) < 10:
        return None

    def mean(k):
        return sum(x[k] for x in tail) / len(tail)

    live = [x for x in tail if x["nb_functional_genes"] > 0]
    return {
        "genes": mean("nb_functional_genes"),
        "links": mean("nb_link"),
        "density": (sum(x["nb_link"] / x["nb_functional_genes"] for x in live)
                    / len(live)) if live else float("nan"),
        "dna": mean("amount_of_dna"),
        "err": mean("metabolic_error"),
    }


def mean_sd(v):
    m = sum(v) / len(v)
    if len(v) < 2:
        return m, 0.0
    return m, math.sqrt(sum((x - m) ** 2 for x in v) / (len(v) - 1))


def welch(a, b):
    ma, sa = mean_sd(a); mb, sb = mean_sd(b)
    na, nb = len(a), len(b)
    se2 = sa * sa / na + sb * sb / nb
    if se2 <= 0:
        return ma - mb, float("nan")
    return ma - mb, (ma - mb) / math.sqrt(se2)


def perm_p(a, b, n=50000):
    obs = abs(sum(a) / len(a) - sum(b) / len(b))
    pool = list(a) + list(b)
    na = len(a)
    hits = 0
    for _ in range(n):
        random.shuffle(pool)
        if abs(sum(pool[:na]) / na - sum(pool[na:]) / (len(pool) - na)) >= obs - 1e-12:
            hits += 1
    return (hits + 1) / (n + 1)


def cohend(a, b):
    ma, sa = mean_sd(a); mb, sb = mean_sd(b)
    na, nb = len(a), len(b)
    sp = math.sqrt(((na - 1) * sa ** 2 + (nb - 1) * sb ** 2) / (na + nb - 2))
    return (ma - mb) / sp if sp > 0 else float("nan")


METRICS = ["genes", "links", "density", "dna", "err"]


def contrast(name, arms, A, B):
    if A not in arms or B not in arms:
        return
    print("\n" + "=" * 80)
    print("%s\n    %s (n=%d)  vs  %s (n=%d)"
          % (name, A, len(arms[A]), B, len(arms[B])))
    print("=" * 80)
    print("%-8s %16s %16s %11s %7s %9s %7s"
          % ("metric", A, B, "diff", "t", "perm p", "d"))
    for m in METRICS:
        a = [r[m] for r in arms[A]]
        b = [r[m] for r in arms[B]]
        ma, sa = mean_sd(a); mb, sb = mean_sd(b)
        diff, t = welch(a, b)
        print("%-8s %8.3f+-%-7.3f %8.3f+-%-7.3f %11.4f %7.2f %9.5f %7.2f"
              % (m, ma, sa, mb, sb, diff, t, perm_p(a, b), cohend(a, b)))


if __name__ == "__main__":
    root = sys.argv[1]
    window = int(sys.argv[2]) if len(sys.argv) > 2 else 2000

    raw = {}
    for d in sorted(os.listdir(root)):
        full = os.path.join(root, d)
        if os.path.isdir(full):
            rows = read_rows(full)
            if rows:
                raw[d] = rows

    if not raw:
        sys.exit("no data under %s" % root)

    cutoff = min(r[-1]["Generation"] for r in raw.values())
    print("runs with data: %d" % len(raw))
    print("common cutoff generation: %d   (window %d)" % (cutoff, window))
    lag = sorted((r[-1]["Generation"], d) for d, r in raw.items())
    print("slowest: %s@%d   fastest: %s@%d"
          % (lag[0][1], lag[0][0], lag[-1][1], lag[-1][0]))

    arms = defaultdict(list)
    for d, rows in raw.items():
        s = summarise(rows, cutoff, window)
        if s:
            arms[d.rsplit("_s", 1)[0]].append(s)

    print("\nARMS")
    for a in sorted(arms):
        print("  %-12s n=%d" % (a, len(arms[a])))

    contrast("C1  FLUCTUATION EFFECT", arms, "fluc_none", "froz_none")
    contrast("C2  YOKED CONTROL -- does information content matter?",
             arms, "fluc_info", "fluc_yoked")
    for nm, A, B in [("S1  fluctuation, informative cue", "fluc_info", "froz_info"),
                     ("S2  fluctuation, yoked cue", "fluc_yoked", "froz_yoked"),
                     ("S3  cue present vs absent (fluctuating)", "fluc_info", "fluc_none"),
                     ("S4  yoked control under frozen env", "froz_info", "froz_yoked")]:
        contrast(nm, arms, A, B)
