#!/usr/bin/env python3
"""Per-run trajectories of the environment-switch experiment to
paper/data/switch_trajectories.csv: means over 500-generation blocks (block b
covers generations 500(b-1)+1 .. 500b) of S = nb_link * avg_influence (0 without
links), nb_TF, nb_link, nb_functional_genes and metabolic_error, plus the source
run's values over generations 25,001-30,000 (block 0). No signals in any cell, so
S is the genome-encoded network. Runs on the Mac mini."""
import csv, os
E = "/Users/ivanosman/raevol_exp"
ROOT = E + "/switch"
OUT = E + "/switch_trajectories.csv"


def blocks(p, lo, hi, width):
    rows = list(csv.reader(open(p)))
    H = {k: i for i, k in enumerate(rows[0])}
    acc = {}
    for r in rows[1:]:
        try:
            g = int(float(r[0]))
        except (ValueError, IndexError):
            continue
        if not (lo < g <= hi):
            continue
        b = (g - lo - 1) // width + 1
        nl = float(r[H["nb_link"]])
        v = (nl * float(r[H["avg_influence"]]) if nl > 0 else 0.0, float(r[H["nb_TF"]]), nl,
             float(r[H["nb_functional_genes"]]), float(r[H["metabolic_error"]]))
        a = acc.setdefault(b, [0.0] * 5 + [0])
        for i in range(5):
            a[i] += v[i]
        a[5] += 1
    return {b: [a[i] / a[5] for i in range(5)] for b, a in acc.items()}


man = list(csv.DictReader(open(os.path.join(ROOT, "manifest.csv"))))
with open(OUT, "w", newline="") as fo:
    w = csv.writer(fo, lineterminator="\n")
    w.writerow(["run", "cell", "rep", "block", "gen_hi", "S", "tfs", "links", "genes", "err"])
    for m in man:
        src = blocks(os.path.join(m["source_dir"], "stats", "stats_simd_best.csv"), 25000, 30000, 5000)
        w.writerow([m["run"], m["cell"], m["rep"], 0, 0] + ["%.6g" % x for x in src[1]])
        p = os.path.join(ROOT, "primary", m["run"], "stats", "stats_simd_best.csv")
        for b, v in sorted(blocks(p, 0, 10000, 500).items()):
            w.writerow([m["run"], m["cell"], m["rep"], b, b * 500] + ["%.6g" % x for x in v])
print("wrote", OUT)
