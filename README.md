# Environmental fluctuation limits gene regulation unless it can be tracked

Code, data and manuscript for the paper of the same title (Ivan Osman,
manuscript in preparation for *PLOS Computational Biology*).

In R-Aevol, a model of genome evolution in which regulation arises from
transcription factors encoded in a DNA-like sequence, we separated
environmental fluctuation from the information available for tracking it. Under
fluctuation, populations received an informative cue, no cue, or a *yoked* cue
that supplies the same signal proteins at the same rates but independently of
the environment. Without a usable cue, fluctuation kept the genome-encoded
regulatory network from growing and removed networks that already existed. With
an informative cue the network was at least as strong as in a constant
environment, and the cue itself was wired into it. The four R-Aevol
experiments that test predictions were pre-registered by committing the
predictions and analysis code before the runs.

## Contents

| Path | What it holds |
|---|---|
| `paper/manuscript.tex`, `manuscript.pdf` | the manuscript (PLOS LaTeX template) |
| `paper/supplement/` | S1–S6 (ODD description of EvoSim, R-Aevol build and configuration, result tables, earlier EvoSim results, prior work, the EvoSim experiment) |
| `paper/data/` | one row per simulation run (or per run and 500-generation block) for every experiment; all figures and tables are computed from these files |
| `paper/make_figures.py` | draws every figure from `paper/data/` |
| `replication/raevol/` | R-Aevol parameter generators, launch scripts, the pre-registered analysis scripts and their single outputs (`results_*_preregistered.txt`), exports and exploratory analyses; see `replication/README.md` |
| `replication/evosim_strength/` | the EvoSim experiment of S6 Text |
| `src/`, `scripts/`, `test/` | EvoSim, the Julia simulator in which the study began (S1, S4 and S6 Text); its documentation is `docs/EvoSim.md` |
| `PAPER.md` | the full research log, in order, including retracted analyses |

## Reproducing the figures, tables and PDFs

From the per-run data in `paper/data/`:

```sh
python3 paper/make_figures.py                     # needs numpy and matplotlib
(cd paper/supplement && python3 make_s3.py)       # S3 Table, standard library only
(cd paper && tectonic manuscript.tex)             # or any LaTeX with the PLOS template files included
(cd paper/supplement && for f in S*_Text.tex S*_Table.tex; do tectonic $f; done)
```

## Re-running the simulations

**R-Aevol.** Branch `aevol-7` of <https://gitlab.inria.fr/aevol/aevol>, commit
`33277374`, built with regulation enabled and per-generation statistics; the
exact build options, the parameter template, targets, signals and seeds are in
S2 Text. `replication/raevol/instrument_signal_links.py` patches a copy of the
source so the statistics also report the signal proteins' links (used for the
confirmation and the environment switch). The parameter generators and
launchers contain the absolute paths of the machine the runs were made on;
change them before use. The analysis scripts need only the Python standard
library.

**EvoSim.** Julia 1.10 or later:

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. -t auto scripts/run.jl --help
```

## Pre-registration record

Each experiment's design, predictions, decision rules and analysis code were
committed before its runs; each analysis was run once, after its SHA-256
checksum was verified (checksums in S2 Text).

| Experiment | Commit |
|---|---|
| Sweep, predictions P1–P4 | `e6173aa` |
| EvoSim confirmation P5 (rule; outcome: not run) | `cbca829`, `15d0418` |
| Environmental grain, G1–G4 | `aafb136` |
| Confirmation of the signal-link decomposition, C1–C4 | `e15104a` |
| Environment switch, W1–W2 (design; analysis script) | `5bd4fee`; `baeddca` |

The commit timestamps are local and cannot be verified independently.

## Data and archive

This repository is archived on Zenodo, doi:
[10.5281/zenodo.24054019](https://doi.org/10.5281/zenodo.24054019) (all
versions; version 1.0 was used for the paper). `paper/data/` holds the per-run
summaries behind every result. The raw per-generation outputs of all runs
(several GB) will be deposited on Zenodo separately and linked here.

## Licence

MIT; see `LICENSE`.
