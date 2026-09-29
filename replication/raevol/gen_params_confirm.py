#!/usr/bin/env python3
"""
R-Aevol pre-registered CONFIRMATION of the signal-link decomposition
(pre-registration: PAPER.md section 11.19).

Parameter files are built by gen_params_v2.make_param, so every setting is
identical to the pre-registered sweep (v2) except the cells and the seeds.
The runs use the instrumented binary (instrument_signal_links.py) from
generation 0, so the genome-encoded part of regulation is measured directly in
the analysis window, without continuations.

Cells, population 1024 (32x32), 12 new seeds each (SEED_BASE 400000):
  q = 0.10 and q = 0.20:  fluc_info, fluc_yoked, fluc_none, froz_info, froz_yoked
  q = 0.10 only:          froz_none (identical targets, no signals: rate-independent)
froz_yoked is rate-dependent (its signals switch at 1.5q), so it is run at both
rates and is the like-for-like frozen baseline for the yoked cue.

Usage:  gen_params_confirm.py <root>
"""
import os, sys, shutil
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gen_params_v2 import make_param, run_name, BINDING_MATRIX, NSEEDS

SEED_BASE = 400000          # v1 100000+, v2 200000+, grain 300000+; no overlap
RATES = [0.10, 0.20]


def cells():
    """(arm, rate) for every cell, in a fixed order."""
    out = []
    for r in RATES:
        for arm in ("fluc_info", "fluc_yoked", "fluc_none", "froz_info", "froz_yoked"):
            out.append((arm, r))
    out.append(("froz_none", 0.10))
    return out


if __name__ == "__main__":
    root = sys.argv[1]
    n = 0
    for ci, (arm, rate) in enumerate(cells()):
        for s in range(1, NSEEDS + 1):
            seed = SEED_BASE + 1000 * ci + s
            name = run_name(arm, rate, 32, 32, s)
            d = os.path.join(root, "primary", name)
            os.makedirs(d, exist_ok=True)
            with open(os.path.join(d, "param.in"), "w") as f:
                f.write(make_param(arm, seed, rate, 32, 32, name))
            shutil.copy(BINDING_MATRIX, os.path.join(d, "binding_matrix.rae"))
            n += 1
    print("wrote %d runs under %s" % (n, root))
