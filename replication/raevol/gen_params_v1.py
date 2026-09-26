#!/usr/bin/env python3
"""
Generate R-Aevol parameter files for an independent replication of the EvoSim
fluctuation experiment, on a completely different model of the genome
(sequence -> RNA -> protein -> Hill-function regulation).

Baseline is the published example examples/dynamic_env/param.in from the
aevol-7 distribution: NB_ENVIRONMENTS 2, ENV_VARIATION switch_in_a_list 0.1,
signals attached to environments. The binding matrix is the one shipped with
that example, used unchanged and identically in every arm, so the "chemistry"
is a fixed shared constant rather than something that varies between arms.

Design, 2 x 3:

  fluctuation
    FLUC  two DIFFERENT phenotypic targets (A and B)
    FROZ  two IDENTICAL targets, so the switching machinery, the PRNG draws
          and the number of target models are all preserved but the organism
          experiences no environmental change

  cue
    NONE   no signal proteins exist
    INFO   2 environments, (A,sig1) and (B,sig2): exactly one signal present
           at any age, and its identity determines the target. I(sig;target)=1 bit
    YOKED  4 environments, (A,sig1) (A,sig2) (B,sig1) (B,sig2): exactly one
           signal present at any age, same two species, same marginal
           frequencies -- but independent of the target. I(sig;target)=0 bits

YOKED is the operant-conditioning yoked control: stimuli matched, contingency
broken.

MATCHING THE SWITCH RATES.  A switch rejection-samples uniformly among the
OTHER environments (PhenotypicTargetHandler_R::ApplyVariation). With 2
environments a switch always changes the target; with 4 it changes the target
with probability 2/3 and the signal with probability 2/3. So to give INFO and
YOKED the same per-age target-change rate AND the same per-age signal-change
rate, YOKED uses switch probability 0.15:

    INFO  (2 env, p=0.10): target changes 0.100/age, signal changes 0.100/age
    YOKED (4 env, p=0.15): target changes 0.15*2/3 = 0.100/age,
                           signal changes 0.15*2/3 = 0.100/age

Without this correction the yoked arm would fluctuate one third less often
than the informative arm, and any difference between them would be
confounded with the fluctuation rate -- the very thing being measured.

Targets A and B are mirror images about x=0.5, so their geometric areas are
equal by construction. The genotype->phenotype map is NOT mirror-symmetric,
so the frozen arms alternate between (A,A) and (B,B) across seeds to balance
any residual asymmetry.
"""
import os, sys, shutil

# The two signal proteins are taken verbatim from the distribution's own
# examples (complex_env signals 1 and 2) rather than invented here.
SIG1 = "h0 w0 h1 m1 w0 h1 m0 h0 h1 w0 h0 m1 h1 h1 m1 m0 h0 w0 h1 m1 w0 h1 m0 h0 h1 w0 h0 m1 h1 h1 m1 m0 h1 m0 m1"
SIG2 = "m0 h0 m1 h1 m1 w0 m0 m1 m0 h0 m1 h1 w0 h0 h0 h1 m1 m0 h1 w0 h1 h0 m1 h1 m0 w0 w0 m0 m1 w0 w0 h1 h0 w0 h1 h0 h0 m0 h0 w0 h0 m1 m0 w0 h1 w0 w0 h1 m0"

# gaussians: (height, mean, sigma)
TARGET_A = [(1.2, 0.52, 0.12), (-1.4, 0.50, 0.07), (0.3, 0.80, 0.03)]
TARGET_B = [(1.2, 0.48, 0.12), (-1.4, 0.50, 0.07), (0.3, 0.20, 0.03)]

P_TWO_ENV  = 0.10   # published value, examples/dynamic_env
P_FOUR_ENV = 0.15   # matches per-age target- and signal-change rates

HEADER = """###################################
#  R-Aevol replication of EvoSim
#  arm={arm}  seed={seed}
###################################
STRAIN_NAME             {arm}_s{seed}
SEED                    {seed}
INIT_POP_SIZE           1024
WORLD_SIZE              32 32
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

BACKUP_STEP             {backup}
TREE_STEP               {backup}
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


def env_layout(arm, seed):
    """Return (list_of_targets, list_of_signal_ids_per_env, switch_p).

    signal id 0 means 'no signal in this environment'.
    """
    fluctuating = arm.startswith("fluc")
    cue = arm.split("_")[1]

    if fluctuating:
        a, b = TARGET_A, TARGET_B
    else:
        # frozen: every environment carries the SAME target
        t = TARGET_A if (seed % 2 == 1) else TARGET_B
        a, b = t, t

    if cue == "none":
        return [a, b], [[], []], P_TWO_ENV
    if cue == "info":
        # (A,sig1) (B,sig2)  -- signal identity determines target
        return [a, b], [[1], [2]], P_TWO_ENV
    if cue == "yoked":
        # (A,sig1) (A,sig2) (B,sig1) (B,sig2) -- signal independent of target
        return [a, a, b, b], [[1], [2], [1], [2]], P_FOUR_ENV
    raise ValueError(cue)


def make_param(arm, seed, backup):
    targets, signals, switch_p = env_layout(arm, seed)
    out = [HEADER.format(arm=arm, seed=seed, switch_p=switch_p,
                         nenv=len(targets), backup=backup)]
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


ARMS = ["fluc_none", "fluc_info", "fluc_yoked",
        "froz_none", "froz_info", "froz_yoked"]

BINDING_MATRIX = "/Users/ivanosman/raevol7/examples/dynamic_env/binding_matrix.rae"

if __name__ == "__main__":
    root, nseeds, backup = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
    for arm in ARMS:
        for s in range(1, nseeds + 1):
            seed = 100000 + 1000 * ARMS.index(arm) + s
            d = os.path.join(root, "%s_s%d" % (arm, s))
            os.makedirs(d, exist_ok=True)
            with open(os.path.join(d, "param.in"), "w") as f:
                f.write(make_param(arm, seed, backup))
            shutil.copy(BINDING_MATRIX, os.path.join(d, "binding_matrix.rae"))
    print("wrote %d runs under %s" % (len(ARMS) * nseeds, root))
