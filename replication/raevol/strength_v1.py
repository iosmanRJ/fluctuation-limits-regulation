#!/usr/bin/env python3
"""
Is the frozen arm's higher link count FUNCTIONAL, or just incidental binding
that scales with having more genes?

nb_link counts every protein->RNA influence, however weak. The frozen arms
carry ~42% more genes, so they have more proteins and more RNAs and therefore
more opportunities for negligible binding. Three sharper measures:

  wlinks   nb_link * avg_influence        strength-weighted connectivity
  perpair  nb_link / (genes * nb_TF)      fraction of POSSIBLE regulator->
                                          target pairs that are realised
  tf       nb_TF                          how many proteins actually regulate

If the fluctuation effect on connectivity survives all three, it is real.
If it disappears under perpair, the raw nb_link difference was a by-product
of gene count and the density claim was measuring size twice.
"""
import csv, math, os, random, sys
from collections import defaultdict

random.seed(11)

KEEP = ("Generation", "nb_functional_genes", "nb_link", "avg_influence",
        "nb_TF", "nb_pure_TF", "metabolic_error")


def rows_of(d):
    p = os.path.join(d, "stats", "stats_simd_best.csv")
    out = []
    if not os.path.exists(p):
        return out
    with open(p) as f:
        for r in csv.DictReader(f):
            try:
                v = {k: float(r[k]) for k in KEEP}
            except (ValueError, TypeError, KeyError):
                continue
            if any(math.isnan(x) for x in v.values()):
                continue
            out.append(v)
    out.sort(key=lambda x: x["Generation"])
    return out


def summarise(rows, cutoff=30000, window=5000):
    tail = [x for x in rows if cutoff - window < x["Generation"] <= cutoff
            and x["nb_functional_genes"] > 0]
    if len(tail) < 10:
        return None
    n = len(tail)
    wl = sum(x["nb_link"] * x["avg_influence"] for x in tail) / n
    pp = sum(x["nb_link"] / (x["nb_functional_genes"] * x["nb_TF"])
             for x in tail if x["nb_TF"] > 0)
    npp = sum(1 for x in tail if x["nb_TF"] > 0)
    return {
        "genes": sum(x["nb_functional_genes"] for x in tail) / n,
        "links": sum(x["nb_link"] for x in tail) / n,
        "density": sum(x["nb_link"] / x["nb_functional_genes"] for x in tail) / n,
        "wlinks": wl,
        "perpair": pp / npp if npp else float("nan"),
        "tf": sum(x["nb_TF"] for x in tail) / n,
        "infl": sum(x["avg_influence"] for x in tail) / n,
    }


def msd(v):
    m = sum(v) / len(v)
    return m, (math.sqrt(sum((x - m) ** 2 for x in v) / (len(v) - 1))
               if len(v) > 1 else 0.0)


def perm(a, b, n=50000):
    obs = abs(sum(a) / len(a) - sum(b) / len(b))
    pool = list(a) + list(b); na = len(a); h = 0
    for _ in range(n):
        random.shuffle(pool)
        if abs(sum(pool[:na]) / na - sum(pool[na:]) / (len(pool) - na)) >= obs - 1e-12:
            h += 1
    return (h + 1) / (n + 1)


root = sys.argv[1]
arms = defaultdict(list)
for d in sorted(os.listdir(root)):
    p = os.path.join(root, d)
    if os.path.isdir(p):
        r = rows_of(p)
        if r:
            s = summarise(r)
            if s:
                arms[d.rsplit("_s", 1)[0]].append(s)

MET = ["genes", "tf", "links", "density", "perpair", "wlinks", "infl"]

print("%-11s %s" % ("arm", " ".join("%10s" % m for m in MET)))
for a in sorted(arms):
    print("%-11s %s" % (a, " ".join("%10.4f" % msd([r[m] for r in arms[a]])[0]
                                    for m in MET)))

for nm, A, B in [("C1  fluc_none vs froz_none", "fluc_none", "froz_none"),
                 ("S1  fluc_info vs froz_info", "fluc_info", "froz_info"),
                 ("C2  fluc_info vs fluc_yoked", "fluc_info", "fluc_yoked"),
                 ("S4  froz_info vs froz_yoked", "froz_info", "froz_yoked")]:
    print("\n=== %s ===" % nm)
    for m in MET:
        a = [r[m] for r in arms[A]]; b = [r[m] for r in arms[B]]
        ma, sa = msd(a); mb, sb = msd(b)
        se = math.sqrt(sa ** 2 / len(a) + sb ** 2 / len(b))
        t = (ma - mb) / se if se > 0 else float("nan")
        print("  %-8s %9.4f+-%-8.4f %9.4f+-%-8.4f  diff %9.4f  t %7.2f  p %.5f"
              % (m, ma, sa, mb, sb, ma - mb, t, perm(a, b)))
