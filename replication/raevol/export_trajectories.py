#!/usr/bin/env python3
"""Per-run trajectories to paper/data/trajectories.csv: means over 500-generation
blocks (block b covers generations 500(b-1)+1 .. 500b) of the best individual's
statistics, for every run of the exploratory experiment (v1), the pre-registered
sweep at population 1024 (v2), the grain experiment and the confirmation.
S = nb_link * avg_influence (0 without links), as in the pre-registered analyses.
S_sig and links_sig exist only for the confirmation (instrumented from generation
0); for cue-free arms S is the genome-encoded network in every experiment.
Runs on the Mac mini, where the run directories live."""
import csv, os, re
E = "/Users/ivanosman/raevol_exp"
SETS = [("v1", E + "/main"), ("v2", E + "/v2/primary"), ("grain", E + "/grain/primary"),
        ("confirm", E + "/confirm/primary")]
OUT = E + "/trajectories.csv"
BLOCK = 500

def cell_of(tag, name):
    base = re.sub(r"_s\d+$", "", name)
    seed = name.rsplit("_s", 1)[1]
    if tag == "v1":
        return base, 0.1, seed
    if tag == "grain":
        arm, g = base.rsplit("_g", 1)
        return arm, int(g) / 1e5, seed
    arm, r, w = base.rsplit("_", 2)
    return arm, int(r[1:]) / 1000, seed

with open(OUT, "w", newline="") as fo:
    w = csv.writer(fo, lineterminator="\n")
    w.writerow(["experiment", "arm", "rate", "seed", "block", "gen_hi", "S", "S_sig", "links", "links_sig", "tfs", "genes", "err"])
    for tag, root in SETS:
        for name in sorted(os.listdir(root)):
            p = os.path.join(root, name, "stats", "stats_simd_best.csv")
            if not os.path.exists(p):
                continue
            arm, rate, seed = cell_of(tag, name)
            rows = list(csv.reader(open(p)))
            H = {k: i for i, k in enumerate(rows[0])}; n0 = len(rows[0])
            acc = {}
            for r in rows[1:]:
                try:
                    g = int(float(r[0]))
                except (ValueError, IndexError):
                    continue
                if g < 1 or g > 30000:
                    continue
                b = (g - 1) // BLOCK + 1
                nl = float(r[H["nb_link"]])
                S = nl * float(r[H["avg_influence"]]) if nl > 0 else 0.0
                extra = [x for x in r[n0:] if x != ""]
                ssig = float(extra[-1]) if (tag == "confirm" and len(extra) >= 2) else float("nan")
                lsig = float(extra[-2]) if (tag == "confirm" and len(extra) >= 2) else float("nan")
                v = (S, ssig, nl, lsig, float(r[H["nb_TF"]]), float(r[H["nb_functional_genes"]]), float(r[H["metabolic_error"]]))
                a = acc.setdefault(b, [0.0] * 7 + [0])
                for i in range(7):
                    a[i] += v[i]
                a[7] += 1
            for b in sorted(acc):
                a = acc[b]; n = a[7]
                w.writerow([tag, arm, rate, seed, b, b * BLOCK] + ["%.6g" % (a[i] / n) for i in range(7)])
print("wrote", OUT)
