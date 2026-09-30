#!/usr/bin/env python3
"""Exploratory: trajectories of regulatory strength in the NO-CUE arms of every
experiment (no signals, so S is the genome-encoded network). Cell means of S,
genes and nb_TF over generation windows."""
import csv, os, re
from collections import defaultdict
E = "/Users/ivanosman/raevol_exp"
SETS = [("v1", E + "/main"), ("v2", E + "/v2/primary"), ("grain", E + "/grain/primary")]
WINS = [(0, 500), (500, 1000), (1000, 2000), (2000, 5000), (5000, 10000), (10000, 20000), (20000, 30000)]
cells = defaultdict(list)
for tag, root in SETS:
    for name in sorted(os.listdir(root)):
        if "none" not in name:
            continue
        p = os.path.join(root, name, "stats", "stats_simd_best.csv")
        if not os.path.exists(p):
            continue
        cell = tag + ":" + re.sub(r"_s\d+$", "", name)
        rows = list(csv.reader(open(p)))
        H = {k: i for i, k in enumerate(rows[0])}
        acc = {w: defaultdict(list) for w in WINS}
        for r in rows[1:]:
            try:
                g = int(float(r[0]))
            except (ValueError, IndexError):
                continue
            nl = float(r[H["nb_link"]])
            S = nl * float(r[H["avg_influence"]]) if nl > 0 else 0.0
            for lo, hi in WINS:
                if lo < g <= hi:
                    acc[(lo, hi)]["S"].append(S)
                    acc[(lo, hi)]["genes"].append(float(r[H["nb_functional_genes"]]))
                    acc[(lo, hi)]["tfs"].append(float(r[H["nb_TF"]]))
        cells[cell].append({w: {k: sum(v) / len(v) for k, v in acc[w].items()} for w in WINS})
print("%-34s" % "cell (n)" + "".join("%16s" % ("%d-%dk" % (lo // 1000, hi // 1000) if lo >= 1000 else "%d-%d" % (lo, hi)) for lo, hi in WINS))
for c in sorted(cells):
    runs = cells[c]
    line = "%-34s" % ("%s (%d)" % (c, len(runs)))
    for w in WINS:
        s = sum(r[w]["S"] for r in runs) / len(runs); t = sum(r[w]["tfs"] for r in runs) / len(runs)
        line += "%16s" % ("%.2f / %.1f" % (s, t))
    print(line)
print("(values: S / nb_TF)")
