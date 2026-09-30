#!/usr/bin/env python3
"""Exploratory: trajectories of the genome-encoded network in the confirmation
runs (instrumented from generation 0). Per run, means over generation windows of
S_int, S_sig, genes and nb_TF; per cell, the mean over runs. Also each run's peak
S_int over 1,000-generation blocks, to see whether networks were built and lost."""
import csv, os, sys
from collections import defaultdict
ROOT = "/Users/ivanosman/raevol_exp/confirm/primary"
WINS = [(0, 500), (500, 1000), (1000, 2000), (2000, 5000), (5000, 10000), (10000, 20000), (20000, 30000)]

def parse(name):
    base = name.rsplit("_s", 1)[0]
    arm, r, w = base.rsplit("_", 2)
    return arm, int(r[1:])

cells = defaultdict(list)
for name in sorted(os.listdir(ROOT)):
    p = os.path.join(ROOT, name, "stats", "stats_simd_best.csv")
    rows = list(csv.reader(open(p)))
    H = {k: i for i, k in enumerate(rows[0])}; n0 = len(rows[0])
    acc = {w: defaultdict(list) for w in WINS}
    blocks = defaultdict(list)
    for r in rows[1:]:
        try:
            g = int(float(r[0]))
        except (ValueError, IndexError):
            continue
        extra = [x for x in r[n0:] if x != ""]
        nl = float(r[H["nb_link"]])
        S = nl * float(r[H["avg_influence"]]) if nl > 0 else 0.0
        ssig = float(extra[-1]) if len(extra) >= 2 else 0.0
        vals = {"S_int": S - ssig, "S_sig": ssig, "genes": float(r[H["nb_functional_genes"]]),
                "tfs": float(r[H["nb_TF"]]) if "nb_TF" in H else float("nan")}
        for lo, hi in WINS:
            if lo < g <= hi:
                for k, v in vals.items():
                    acc[(lo, hi)][k].append(v)
        blocks[(g - 1) // 1000].append(S - ssig)
    run = {w: {k: sum(v) / len(v) for k, v in acc[w].items()} for w in WINS}
    peak = max(sum(v) / len(v) for v in blocks.values())
    peak_block = max(blocks, key=lambda b: sum(blocks[b]) / len(blocks[b]))
    final = sum(blocks[29]) / len(blocks[29])
    cells[parse(name)].append((run, peak, peak_block, final))

print("cell mean per window: S_int / genes / TFs")
print("%-18s" % "cell" + "".join("%18s" % ("%d-%dk" % (lo // 1000, hi // 1000) if lo >= 1000 else "%d-%d" % (lo, hi)) for lo, hi in WINS))
for c in sorted(cells):
    runs = cells[c]
    line = "%-18s" % ("%s q=%.2f" % (c[0], c[1] / 1000))
    for w in WINS:
        s = sum(r[0][w]["S_int"] for r in runs) / len(runs)
        g = sum(r[0][w]["genes"] for r in runs) / len(runs)
        t = sum(r[0][w]["tfs"] for r in runs) / len(runs)
        line += "%18s" % ("%.2f/%.1f/%.1f" % (s, g, t))
    print(line)
print("\npeak 1,000-generation block of S_int vs final block (29-30k), per cell")
for c in sorted(cells):
    runs = cells[c]
    pk = sum(r[1] for r in runs) / len(runs); fin = sum(r[3] for r in runs) / len(runs)
    early = sum(1 for r in runs if r[2] < 10)
    print("%-18s mean peak %.2f  mean final %.2f  ratio final/peak %.2f  runs peaking before 10k: %d/%d" % ("%s q=%.2f" % (c[0], c[1] / 1000), pk, fin, fin / pk if pk else float("nan"), early, len(runs)))
