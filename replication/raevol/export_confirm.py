#!/usr/bin/env python3
"""One row per confirmation run (paper/data/confirm_runs.csv), summarised with the
pre-registered function analyse_confirm.summarise (generations 25,001-30,000;
S = nb_link * avg_influence, 0 without links; S_sig from the instrumented columns;
S_int = S - S_sig). Runs on the Mac mini, where the run directories live."""
import csv, os, sys
sys.path.insert(0, "/Users/ivanosman/raevol_exp")
from analyse_confirm import summarise, parse

R = "/Users/ivanosman/raevol_exp/confirm/primary"
with open("/Users/ivanosman/raevol_exp/confirm_runs.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["arm", "rate", "pop", "seed", "S", "S_sig", "S_int", "links", "links_sig", "genes", "err"])
    for name in sorted(os.listdir(R)):
        s, why = summarise(os.path.join(R, name))
        if s is None:
            print("skip", name, why); continue
        arm, r, pop = parse(name)
        seed = name.rsplit("_s", 1)[1]
        w.writerow([arm, r / 1000, pop, seed] + ["%.6g" % s[k] for k in ("S", "S_sig", "S_int", "links", "links_sig", "genes", "err")])
print("wrote confirm_runs.csv")
