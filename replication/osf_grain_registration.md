# OSF registration text — grain experiment (G1–G4)

For pasting into an OSF registration ("OSF Preregistration" or the Simulation
Studies template). The value of registering is an independent timestamp, so it
should be submitted **before** `grain_results.txt` exists on the Mac mini
(expected around 01:30, 2026-09-28).

---

**Title.** Does environmental grain decide the sign of cue-free fluctuation's
effect on genome size and regulatory strength? A pre-registered R-Aevol
experiment.

**Status of data at registration.** The 60 simulation runs began on
2026-09-27 at 18:24 and are in progress. No outcome data have been inspected.
The analysis script runs automatically, once, when the last run finishes. The
design, analysis code and decision rules were committed to version control
before any run (git commit fa00777; SHA-256 prefixes: `analyse_grain.py`
89a8613289731999, `gen_params_grain.py` 2048f1443cc9367e, `analyse_v2.py`
d0609140a970de48).

**Background.** Tsuda & Kawata (2010, PLoS Comput Biol 6:e1000873) report that
cue-free fluctuation of the optimum across generations promotes gene-regulatory
complexity. In a pre-registered R-Aevol sweep we found that cue-free switching
within a lifetime weakens regulation and shrinks the genome. We test whether
environmental grain reconciles the two.

**Hypothesis (H_grain).** The sign of the effect of cue-free fluctuation depends
on grain: negative when fine-grained, non-negative when coarse-grained.

**Design.** R-Aevol (Aevol branch `aevol-7`, commit 33277374, built with
`-Dwith-raevol=ON`), code unchanged. The environment is a Markov chain over the
20 ages of a lifetime, continuing across generations, switching with
probability p per age between two mirror-image target functions of equal area.
Arms: cue-free fluctuating at p = 0.1, 0.01, 0.001 and 0.0002 (mean dwell 10,
100, 1,000 and 5,000 ages, i.e. half a lifetime to ~250 generations), and a
frozen arm (identical targets). 12 seeds per arm, population 1024, 30,000
generations. All other parameters as in the earlier sweep.

**Measures.** Per run, mean over generations 25,001–30,000 of the best
individual: genes = number of functional genes; S = total regulatory strength =
nb_link × avg_influence (0 when there are no links).

**Predictions** (coarse = p 0.0002, fine = p 0.1):
- G1: genes(coarse) > genes(fine)
- G2: genes(coarse) > genes(frozen)
- G3: S(coarse) > S(fine)
- G4: S(coarse) > S(frozen)

**Analysis.** Two-sided permutation test on the difference of means, 100,000
permutations; Holm correction across the four. Verdict per prediction:
SUPPORTED if the predicted sign holds with Holm-adjusted p < 0.05; REFUTED if the
opposite sign holds with unadjusted p < 0.05; INCONCLUSIVE otherwise. Runs not
reaching generation 30,000 are excluded and listed.

**Interpretation fixed in advance.** If both predictions in a pair (genes G1/G2,
regulation G3/G4) are supported, grain reconciles the conflict. If only the
first is supported, coarse grain removes the fine-grain penalty without a gain;
the gain Tsuda & Kawata report probably requires novel optima, which this
two-target design cannot produce. If the first is not supported, grain is not
the reconciling variable here.

**Secondary, exploratory.** Cell means; Spearman correlation of genes and S
with log10(dwell); fine-grain sanity check (p = 0.1 against frozen);
metabolic error by cell.
