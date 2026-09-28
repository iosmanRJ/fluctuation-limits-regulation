# Replication

Code behind PAPER.md §11. Everything here runs on the Mac mini; paths inside
the scripts point there.

## `raevol/` — independent replication in R-Aevol

R-Aevol is the regulation-enabled build of Aevol: branch `aevol-7` of
`https://gitlab.inria.fr/aevol/aevol.git`, compiled with `-Dwith-raevol=ON`.
Build notes for Apple Silicon are in the project memory (`raevol-replication`).

| file | role |
|---|---|
| `gen_params_v1.py` | exploratory design: 6 arms × 12 seeds, switch rate 0.10 |
| `strength_v1.py` | exploratory analysis that found the count/strength divergence |
| `results_v1_analysis.txt` | v1 output, first-pass analysis |
| `gen_params_v2.py` | **pre-registered** sweep: 3 switch rates, 2 population sizes, new seeds |
| `analyse_v2.py` | **pre-registered** measures, tests and decision rules — do not edit after launch |
| `launch_v2.sh` | launcher: primary cells first, seed-major, `DONE_<exit>` markers |
| `results_v2_preregistered.txt` | output of the single pre-registered run of `analyse_v2.py`: P1–P4 all SUPPORTED (PAPER.md §11.7) |
| `gen_params_grain.py`, `analyse_grain.py` | **pre-registered** grain experiment (G1–G4), committed as `fa00777` before any run |
| `results_grain_preregistered.txt` | its single analysis run: G1, G3 inconclusive; G2, G4 refuted (PAPER.md §11.10) |

The yoked control uses four environments against the informative arm's two,
at 1.5× the switch probability, so per-age target- and signal-change rates are
equal. See the docstring in `gen_params_v2.py`.

## `evosim_strength/` — EvoSim, edge count versus wiring strength

Exploratory. `mean_abs_w` is Σ|w| / ng² (all slots, zeros included), i.e. a
strength *density*; see the header of `analyse.jl` before deriving anything
from it. The pre-registered confirmation is P5 in PAPER.md §11.5.
