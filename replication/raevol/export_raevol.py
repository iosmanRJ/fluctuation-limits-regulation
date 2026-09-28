#!/usr/bin/env python3
"""Export one row per R-Aevol run (v1 exploratory + v2 pre-registered) to
paper/data/raevol_runs.csv, using the SAME summary function as the
pre-registered analysis (analyse_v2.summarise: mean over generations
25,001-30,000; S = nb_link*avg_influence, 0 if no links).

tfs (mean nb_TF over the same window) is not part of the pre-registered
summary; it is computed here with the same window and reported only for the
exploratory experiment. Runs on the Mac mini, where the run directories live."""
import csv, os, sys
sys.path.insert(0, "/Users/ivanosman/raevol_exp")
from analyse_v2 import summarise, parse, WIN_LO, WIN_HI

R = "/Users/ivanosman/raevol_exp"


def mean_tfs(d):
    v = []
    with open(os.path.join(d, "stats", "stats_simd_best.csv")) as f:
        for r in csv.DictReader(f):
            try:
                g = int(float(r["Generation"]))
            except (ValueError, KeyError):
                continue
            if WIN_LO < g <= WIN_HI:
                v.append(float(r["nb_TF"]))
    return sum(v) / len(v)


out = open(os.path.join(R, "raevol_runs.csv"), "w", newline="")
w = csv.writer(out)
w.writerow(["experiment", "arm", "rate", "pop", "seed", "S", "genes", "links", "perlink", "err", "tfs"])

# v1: main/<arm>_s<seed>, rate 0.10, pop 1024
for name in sorted(os.listdir(os.path.join(R, "main"))):
    d = os.path.join(R, "main", name)
    if not os.path.isdir(d):
        continue
    s, why = summarise(d)
    if s is None:
        print("v1 skip", name, why); continue
    arm, seed = name.rsplit("_s", 1)
    w.writerow(["v1", arm, 0.10, 1024, seed, s["S"], s["genes"], s["links"], s["perlink"], s["err"], mean_tfs(d)])

# v2: v2/{primary,secondary}/<arm>_r<rate>_w<pop>_s<seed>
for tier in ("primary", "secondary"):
    td = os.path.join(R, "v2", tier)
    for name in sorted(os.listdir(td)):
        d = os.path.join(td, name)
        if not os.path.isdir(d):
            continue
        s, why = summarise(d)
        if s is None:
            print("v2 skip", name, why); continue
        arm, rate, pop = parse(name)
        seed = name.rsplit("_s", 1)[1]
        w.writerow(["v2", arm, rate / 1000, pop, seed, s["S"], s["genes"], s["links"], s["perlink"], s["err"], mean_tfs(d)])
out.close()
print("wrote", os.path.join(R, "raevol_runs.csv"))
