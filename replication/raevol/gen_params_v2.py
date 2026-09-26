#!/usr/bin/env python3
"""
R-Aevol pre-registered sweep (v2). Pre-registration: PAPER.md section 11.

Identical to gen_params_v1.py (the exploratory 72-run experiment) in every
respect EXCEPT that the switch rate and the population size are swept and the
seeds are new. Targets, signals, mutation rates, regulation parameters, the
binding matrix and the yoked-control construction are unchanged, so any
difference from v1 is attributable to rate, population, or seed.

Cells
  PRIMARY    switch rate p in {0.05, 0.10, 0.20}, population 1024 (32x32)
               fluc_info, fluc_yoked, fluc_none, froz_info      at every rate
               froz_none, froz_yoked                            at p = 0.10 only
  SECONDARY  population 512 (32x16), p = 0.10
               fluc_info, fluc_yoked, fluc_none, froz_info

Why froz_none runs at one rate only: its two targets are identical and it has
no signals, so the switch rate changes nothing the organism can experience --
it is the same frozen baseline at every rate. froz_info is NOT rate-independent
(its signal still switches), so it is run at every rate.

Switch-rate matching (see v1): a switch rejection-samples among the OTHER
environments, so 4-environment arms (yoked) use 1.5x the 2-environment rate to
equalise per-age target-change and signal-change rates.

Usage:  gen_params_v2.py <root>
"""
import os, sys, shutil

SIG1 = "h0 w0 h1 m1 w0 h1 m0 h0 h1 w0 h0 m1 h1 h1 m1 m0 h0 w0 h1 m1 w0 h1 m0 h0 h1 w0 h0 m1 h1 h1 m1 m0 h1 m0 m1"
SIG2 = "m0 h0 m1 h1 m1 w0 m0 m1 m0 h0 m1 h1 w0 h0 h0 h1 m1 m0 h1 w0 h1 h0 m1 h1 m0 w0 w0 m0 m1 w0 w0 h1 h0 w0 h1 h0 h0 m0 h0 w0 h0 m1 m0 w0 h1 w0 w0 h1 m0"

TARGET_A = [(1.2, 0.52, 0.12), (-1.4, 0.50, 0.07), (0.3, 0.80, 0.03)]
TARGET_B = [(1.2, 0.48, 0.12), (-1.4, 0.50, 0.07), (0.3, 0.20, 0.03)]

NSEEDS = 12
SEED_BASE = 200000          # v1 used 100000+; no overlap

RATES = [0.05, 0.10, 0.20]  # 2-environment switch probability per age

HEADER = """###################################
#  R-Aevol pre-registered sweep v2
#  arm={arm} rate={rate} world={wx}x{wy} seed={seed}
###################################
STRAIN_NAME             {name}
SEED                    {seed}
INIT_POP_SIZE           {npop}
WORLD_SIZE              {wx} {wy}
INIT_METHOD             ONE_GOOD_GENE CLONE
CHROMOSOME_INITIAL_LENGTH   5000
FUZZY_FLAVOR            3

SELECTION_SCHEME        fitness_proportionate   750

POINT_MUTATION_RATE     5e-5
SMALL_INSERTION_RATE    5e-5
SMALL_DELETION_RATE     5e-5
MAX_INDEL_SIZE          6

WITH_ALIGNMENTS         false
DUPLICATION_RATE        5e-5
DELETION_RATE           5e-5
TRANSLOCATION_RATE      5e-5
INVERSION_RATE          5e-5

BACKUP_STEP             10000
TREE_STEP               10000
RECORD_TREE             false
MORE_STATS              false

####### Regulation #######
HILL_SHAPE_N                   4
HILL_SHAPE_THETA               0.5
DEGRADATION_RATE               1
NB_DEGRADATION_STEP            10
NB_INDIV_AGE                   20
INDIVIDUAL_EVALUATION_AGES     1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
BINDING_ZEROS_PERCENTAGE       0.75
RANDOM_BINDING_MATRIX          false
WITH_HEREDITY                  false
PROTEIN_PRESENCE_LIMIT         1e-2

####### Environment #######
ENV_SAMPLING            300
MAX_TRIANGLE_WIDTH      0.033333333
ENV_AXIS_FEATURES       METABOLISM
ENV_VARIATION           switch_in_a_list {switch_p}

NB_ENVIRONMENTS {nenv}
"""


def env_layout(arm, seed, rate):
    fluctuating = arm.startswith("fluc")
    cue = arm.split("_")[1]
    if fluctuating:
        a, b = TARGET_A, TARGET_B
    else:
        t = TARGET_A if (seed % 2 == 1) else TARGET_B
        a, b = t, t
    if cue == "none":
        return [a, b], [[], []], rate
    if cue == "info":
        return [a, b], [[1], [2]], rate
    if cue == "yoked":
        return [a, a, b, b], [[1], [2], [1], [2]], round(1.5 * rate, 6)
    raise ValueError(cue)


def make_param(arm, seed, rate, wx, wy, name):
    targets, signals, p = env_layout(arm, seed, rate)
    out = [HEADER.format(arm=arm, rate=rate, wx=wx, wy=wy, seed=seed,
                         name=name, npop=wx * wy, switch_p=p, nenv=len(targets))]
    for i, tgt in enumerate(targets, start=1):
        for (h, m, s) in tgt:
            out.append("ENV_ADD_GAUSSIAN  %d  %s  %s  %s\n" % (i, h, m, s))
    out.append("\n")
    if any(signals):
        out.append("CREATE_SIGNAL %s\n" % SIG1)
        out.append("CREATE_SIGNAL %s\n" % SIG2)
        for i, sigs in enumerate(signals, start=1):
            for sid in sigs:
                out.append("ENV_ADD_SIGNAL %d %d\n" % (i, sid))
    return "".join(out)


def cells():
    """(tier, arm, rate, wx, wy) for every cell, in a fixed order."""
    out = []
    for r in RATES:
        for arm in ("fluc_info", "fluc_yoked", "fluc_none", "froz_info"):
            out.append(("primary", arm, r, 32, 32))
    for arm in ("froz_none", "froz_yoked"):
        out.append(("primary", arm, 0.10, 32, 32))
    for arm in ("fluc_info", "fluc_yoked", "fluc_none", "froz_info"):
        out.append(("secondary", arm, 0.10, 32, 16))
    return out


def run_name(arm, rate, wx, wy, s):
    return "%s_r%03d_w%d_s%d" % (arm, round(rate * 1000), wx * wy, s)


BINDING_MATRIX = "/Users/ivanosman/raevol7/examples/dynamic_env/binding_matrix.rae"

if __name__ == "__main__":
    root = sys.argv[1]
    n = 0
    for ci, (tier, arm, rate, wx, wy) in enumerate(cells()):
        for s in range(1, NSEEDS + 1):
            seed = SEED_BASE + 1000 * ci + s
            name = run_name(arm, rate, wx, wy, s)
            d = os.path.join(root, tier, name)
            os.makedirs(d, exist_ok=True)
            with open(os.path.join(d, "param.in"), "w") as f:
                f.write(make_param(arm, seed, rate, wx, wy, name))
            shutil.copy(BINDING_MATRIX, os.path.join(d, "binding_matrix.rae"))
            n += 1
    print("wrote %d runs under %s (primary + secondary)" % (n, root))
