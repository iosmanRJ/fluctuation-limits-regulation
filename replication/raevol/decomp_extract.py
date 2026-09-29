#!/usr/bin/env python3
"""One row per continued R-Aevol run (decomp_launch.sh, generations 30,001-32,000,
instrumented statistics), written to decomp_runs.csv. EXPLORATORY: this analysis
was not pre-registered.

Per generation, from stats/stats_simd_best.csv:
  S      = nb_link * avg_influence (0 without links)     total regulatory strength
  S_sig  = sum_signal_influence  (instrumented column)   part exerted by signal proteins
  S_int  = S - S_sig                                     genome-encoded ("internal") part
  links, links_sig, links_int likewise from nb_link and nb_signal_link.
The two instrumented fields follow the original columns (header unchanged on
resume), so they are read by position. Per-run values are means over the window.

Usage: decomp_extract.py <decomp_root> <out_csv>
"""
import csv, os, sys

root, out = sys.argv[1], sys.argv[2]
LO, HI = 30000, 32000
w = csv.writer(open(out, "w", newline=""))
w.writerow(["experiment", "arm", "rate", "pop", "seed", "gens", "S", "S_sig", "S_int",
            "links", "links_sig", "links_int", "genes", "err", "tfs", "min_S_int"])
for tier in sorted(os.listdir(root)):
    for name in sorted(os.listdir(os.path.join(root, tier))):
        d = os.path.join(root, tier, name)
        p = os.path.join(d, "stats", "stats_simd_best.csv")
        if not os.path.exists(p):
            print("no stats", d); continue
        with open(p) as f:
            rows = list(csv.reader(f))
        head = rows[0]; H = {k: i for i, k in enumerate(head)}
        n0 = len(head)
        acc = []
        for r in rows[1:]:
            try:
                g = int(float(r[0]))
            except (ValueError, IndexError):
                continue
            if not (LO < g <= HI):
                continue
            extra = [x for x in r[n0:] if x != ""]
            if len(extra) < 2:
                sys.exit("row without instrumented fields in %s gen %d" % (p, g))
            nl = float(r[H["nb_link"]])
            S = nl * float(r[H["avg_influence"]]) if nl > 0 else 0.0
            ls, ss = float(extra[-2]), float(extra[-1])
            acc.append((S, ss, S - ss, nl, ls, nl - ls, float(r[H["nb_functional_genes"]]),
                        float(r[H["metabolic_error"]]), float(r[H["nb_TF"]])))
        if len(acc) < HI - LO:
            print("incomplete", d, len(acc)); continue
        m = [sum(x[i] for x in acc) / len(acc) for i in range(9)]
        # arm / rate / pop / seed from the directory name
        if tier == "main":
            arm, seed = name.rsplit("_s", 1); rate, pop, exp = 0.10, 1024, "v1"
        else:
            base, seed = name.rsplit("_s", 1)
            arm, rtag, wtag = base.rsplit("_", 2)
            rate, pop, exp = int(rtag[1:]) / 1000, int(wtag[1:]), "v2"
        w.writerow([exp, arm, rate, pop, seed, len(acc)] + ["%.6g" % v for v in m]
                   + ["%.6g" % min(x[2] for x in acc)])
print("wrote", out)
