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
| `instrument_signal_links.py` | patches a COPY of the aevol-7 source so the best-individual stats also report links and summed influence of SIGNAL proteins (two columns appended at the end of each row) |
| `decomp_launch.sh`, `decomp_extract.py`, `decomp_analyse.py` | EXPLORATORY decomposition: continue every run from its generation-30,000 population to 32,000 with the instrumented binary, summarise per run (`paper/data/decomp_runs.csv`), and test contrasts on the genome-encoded part → `results_decomp_exploratory.txt`. Resuming is not bit-reproducible, hence continuation rather than re-run; validated against the original window (r = 0.995) |
| `gen_params_confirm.py`, `analyse_confirm.py` | **pre-registered** confirmation of the decomposition (C1–C4 on the genome-encoded strength), committed as `59ff365` before any run; runs use the instrumented binary from generation 0 |
| `launch_confirm.sh`, `chain_confirm.sh` | launch the 132 confirmation runs under `caffeinate`, then verify the checksum of `analyse_confirm.py` and run it once |
| `results_confirm_preregistered.txt` | its single analysis run: C1–C4 all SUPPORTED (PAPER.md §11.20) |
| `export_confirm.py` | writes `paper/data/confirm_runs.csv` with the pre-registered `summarise` |
| `export_raevol.py`, `export_grain.py` | write `paper/data/raevol_runs.csv` and `grain_runs.csv`, one row per run, with the pre-registered summary function (plus `tfs`, mean `nb_TF`, not pre-registered) |

The yoked control uses four environments against the informative arm's two,
at 1.5× the switch probability, so per-age target- and signal-change rates are
equal. See the docstring in `gen_params_v2.py`.

## `evosim_strength/` — EvoSim, edge count versus wiring strength

Exploratory. `mean_abs_w` is Σ|w| / ng² (all slots, zeros included), i.e. a
strength *density*; see the header of `analyse.jl` before deriving anything
from it. The pre-registered confirmation is P5 in PAPER.md §11.5.
`export_evosim.jl` writes `paper/data/evosim_runs.csv`; `gen_time.jl` gives
steps per generation (fluctuating 23.6–24.3, frozen 30.7).
