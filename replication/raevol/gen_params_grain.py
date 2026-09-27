#!/usr/bin/env python3
"""
R-Aevol grain experiment (G). Pre-registration: PAPER.md section 11.8.

Question: does the SIGN of the effect of cue-free fluctuation on genome size and
regulatory strength depend on environmental GRAIN? Tsuda & Kawata 2010 report
that cue-free fluctuation of the optimum ACROSS generations promotes GRN
complexity; P3/P4 (section 11.7) found that cue-free switching WITHIN a lifetime
weakens regulation and shrinks the genome.

Grain is one parameter here, with no change to the code. In the SIMD engine
(SIMD_PhenotypicTargetHandler_R::ApplyVariation) the environment is a single
Markov chain over ages that continues from one generation to the next: at every
one of the 20 ages it switches with probability p. Mean dwell = 1/p ages.

    p        mean dwell        grain
    0.1      10 ages           fine: switches within almost every lifetime
    0.01     100 ages = 5 gen  intermediate
    0.001    1000 ages = 50 gen      coarse
    0.0002   5000 ages = 250 gen     very coarse: Tsuda & Kawata's regime

The environment is population-wide (one PRNG), as in Tsuda & Kawata. Unlike
theirs, it alternates between two fixed targets rather than random-walking to
new optima -- their CYCLIC case, which they found promotes complexity less than
the random walk. That limit travels with the result.

Everything else is identical to gen_params_v2.py (imported from it): targets,
mutation, regulation, binding matrix, population 1024, 30,000 generations.
Cells: fluc_none at each p; froz_none once (identical targets, rate-independent).
12 seeds per cell, new seeds (base 300000; v1 used 100000, v2 200000).

Usage:  gen_params_grain.py <root>
"""
import os, sys, shutil
from gen_params_v2 import make_param, BINDING_MATRIX

NSEEDS = 12
SEED_BASE = 300000
RATES = [0.1, 0.01, 0.001, 0.0002]


def cells():
    out = [("fluc_none", r) for r in RATES]
    out.append(("froz_none", 0.1))
    return out


def run_name(arm, rate, s):
    # rate tag = p x 1e5, so 0.1 -> 10000, 0.0002 -> 20
    return "%s_g%d_s%d" % (arm, round(rate * 1e5), s)


if __name__ == "__main__":
    root = sys.argv[1]
    n = 0
    for ci, (arm, rate) in enumerate(cells()):
        for s in range(1, NSEEDS + 1):
            seed = SEED_BASE + 1000 * ci + s
            name = run_name(arm, rate, s)
            d = os.path.join(root, name)
            os.makedirs(d, exist_ok=True)
            with open(os.path.join(d, "param.in"), "w") as f:
                f.write(make_param(arm, seed, rate, 32, 32, name))
            shutil.copy(BINDING_MATRIX, os.path.join(d, "binding_matrix.rae"))
            n += 1
    print("wrote %d runs under %s" % (n, root))
