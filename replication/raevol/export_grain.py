#!/usr/bin/env python3
"""One row per grain-experiment run to paper/data/grain_runs.csv, summarised
with the pre-registered function (analyse_v2.summarise: generations
25,001-30,000; S = nb_link * avg_influence, 0 when there are no links).
Runs on the Mac mini, where the run directories live."""
import csv, os, sys
sys.path.insert(0, "/Users/ivanosman/raevol_exp")
from analyse_v2 import summarise

R = "/Users/ivanosman/raevol_exp/grain/primary"
with open("/Users/ivanosman/raevol_exp/grain_runs.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["arm", "q", "dwell_ages", "seed", "S", "genes", "links", "perlink", "err"])
    for name in sorted(os.listdir(R)):
        d = os.path.join(R, name)
        if not os.path.isdir(d):
            continue
        s, why = summarise(d)
        if s is None:
            print("skip", name, why); continue
        base, seed = name.rsplit("_s", 1)
        arm, tag = base.rsplit("_g", 1)
        q = int(tag) / 1e5
        w.writerow([arm, q, 1 / q, seed, s["S"], s["genes"], s["links"], s["perlink"], s["err"]])
print("wrote grain_runs.csv")
