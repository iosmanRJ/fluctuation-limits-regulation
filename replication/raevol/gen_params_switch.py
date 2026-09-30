#!/usr/bin/env python3
"""
R-Aevol ENVIRONMENT-SWITCH experiment (pre-registration: PAPER.md section 11.24).

Source populations: the cue-free cells of the confirmation (population 1024,
q = 0.10), best individual of each run at generation 30,000. Each new run is a
clonal population of that genome (aevol_create -C), created and run with the
instrumented binary for 10,000 generations. Four cells, 12 runs each, run s
drawn from source replicate s:

  FF  source froz_none_r100_w1024_s<s>  ->  frozen (same target as the source)   control
  FU  source froz_none_r100_w1024_s<s>  ->  cue-free fluctuation, q = 0.10        switch
  UU  source fluc_none_r100_w1024_s<s>  ->  cue-free fluctuation, q = 0.10        control
  UF  source fluc_none_r100_w1024_s<s>  ->  frozen                               switch

Parameter files come from gen_params_v2.make_param, so every setting is that of
the sweep and the confirmation. Frozen targets follow the seed's parity (A odd,
B even); SEED = 500000 + 1000*c + s has the parity of s, as the source froz_none
seeds (410000 + s) do, so FF keeps the source's target.

Usage:  gen_params_switch.py <root>     writes <root>/primary/<run>/ and manifest.csv
"""
import csv, os, sys, shutil
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gen_params_v2 import make_param, BINDING_MATRIX, NSEEDS

SEED_BASE = 500000          # v1 1e5+, v2 2e5+, grain 3e5+, confirmation 4e5+
SOURCE = "/Users/ivanosman/raevol_exp/confirm/primary"
CELLS = [  # (cell, source arm, target arm)
    ("FF", "froz_none", "froz_none"),
    ("FU", "froz_none", "fluc_none"),
    ("UU", "fluc_none", "fluc_none"),
    ("UF", "fluc_none", "froz_none"),
]
RATE = 0.10

if __name__ == "__main__":
    root = sys.argv[1]
    os.makedirs(os.path.join(root, "primary"), exist_ok=True)
    rows = []
    for ci, (cell, src_arm, tgt_arm) in enumerate(CELLS):
        for s in range(1, NSEEDS + 1):
            seed = SEED_BASE + 1000 * ci + s
            name = "%s_s%d" % (cell, s)
            d = os.path.join(root, "primary", name)
            os.makedirs(d, exist_ok=True)
            with open(os.path.join(d, "param.in"), "w") as f:
                f.write(make_param(tgt_arm, seed, RATE, 32, 32, name))
            shutil.copy(BINDING_MATRIX, os.path.join(d, "binding_matrix.rae"))
            src = os.path.join(SOURCE, "%s_r100_w1024_s%d" % (src_arm, s))
            assert os.path.exists(os.path.join(src, "world", "world_000030000.ae")), src
            rows.append([name, cell, s, seed, src_arm, tgt_arm, src])
    with open(os.path.join(root, "manifest.csv"), "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["run", "cell", "rep", "seed", "source_arm", "target_arm", "source_dir"])
        w.writerows(rows)
    print("wrote %d runs under %s" % (len(rows), root))
