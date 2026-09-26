# Paper preparation

Everything needed to write up this work: the claims with their support, what is
new against what reproduces published results, the caveats that must travel with
each number, and what is still open.

**Status of the evidence base:** ~1,835 simulation runs across 25 experiments,
plus a completed 8-seed × 25M-tick long run (~630,000 generations per seed).
Two independent code reviews. One defect found and fixed.

> **SUPERSEDED IN PART — 2026-09-26.** The independent R-Aevol replication did
> not reproduce §1.1 or §1.2 as stated, and an EvoSim re-analysis by wiring
> *strength* rather than edge *count* reverses the sign of §1.1's density
> effect. §§1–10 are kept as the record of what was claimed; the current state
> of the evidence, the prior-art position and the pre-registered test that
> decides the revised claim are in **§11**. Do not draft from §1.1, §1.2 or §10
> without reading §11 first.

---

## 1. What to claim, ranked by strength

### 1.1 Fluctuation raises connection density, not just network size

The central result. **7 of 7 independent comparisons positive on node count,
every one significant.**

| measure | comparisons | mean effect | range |
|---|---|---|---|
| node count `ng` | 7/7 positive | **+6.32** (~+32%) | +2.64 to +10.15 |
| edge count | 6/6 positive | **+104%** | +63% to +122% |
| **connection density** | **5/6 significant** | **+0.09 to +0.14** | see below |

**State the claim as density, not as raw edge count.** Edges scale as `ng²`, so a
32% rise in nodes mechanically yields ~75% more edges *at constant density*. The
raw "+104% edges against +32% nodes" comparison therefore overstates the case,
and a reviewer will say so immediately. The measured excess over that arithmetic
baseline:

| experiment | `ng` % | predicted from `ng²` | observed | excess | density frozen → fluctuating | p |
|---|---|---|---|---|---|---|
| SP | +32% | +75% | +112% | +37% | 0.827 → 0.916 | <0.0001 |
| SWA | +42% | +102% | +98% | **−4%** | 0.987 → 0.988 | 0.64 |
| SWB | +35% | +82% | +122% | +41% | 0.823 → 0.922 | <0.0001 |
| SWC | +16% | +36% | +63% | +27% | 0.618 → 0.722 | 0.009 |
| N2-4 | +34% | +78% | +116% | +38% | 0.606 → 0.750 | <0.0001 |
| N2-8 | +32% | +75% | +112% | +37% | 0.592 → 0.730 | <0.0001 |

Mean excess **+29 percentage points**, and density rises significantly in 5 of 6.

**The exception is the evidence.** SWA is the saturated regime (density 0.987),
where density *cannot* rise — and there the excess is −4%, i.e. the edge increase
is exactly what node growth predicts and nothing more. Where density can move it
does; where it cannot, the effect reduces to arithmetic. That internal
consistency is worth reporting explicitly.

**It is not a population artefact.** Frozen populations are larger (3,238 vs
2,583), and this work's own results have population driving genome size, so the
confound must be addressed. Regressing on both:

| response | fluctuating | log₂(population) |
|---|---|---|
| `ng` | **+4.76 (t = +4.76)** | −2.80 (t = −2.15) |
| density (SP) | **+0.084 (t = +9.88)** | −0.013 (t = −1.15, n.s.) |
| density (N2-8) | **+0.127 (t = +10.89)** | −0.028 (t = −1.85, n.s.) |

The fluctuation term survives, and the population term carries the **expected
negative sign**, consistent with Lynch and with §1.4. The density effect is
essentially population-independent.

**It is robust to development time.** `dev_iters` controls the phenotype
completely — changing it 8 → 9 or 8 → 10 flips the top uptake channel in 100% of
organisms (§3) — so it is the sharpest available test of whether the result is an
artefact of one setting. Re-run at 20 seeds per arm with sparsity on:

| `dev_iters` | Δ`ng` | p | Δ density | p | excess over `ng²` |
|---|---|---|---|---|---|
| 4 | +6.28 | <0.0001 | **+0.144** | <0.0001 | +38% |
| 8 | +5.81 | <0.0001 | **+0.138** | <0.0001 | +37% |
| 12 | +5.71 | <0.0001 | **+0.118** | <0.0001 | +31% |
| 20 | +5.36 | <0.0001 | **+0.120** | <0.0001 | +30% |

`ng` spans +5.4 to +6.3, density +0.118 to +0.144, excess +30% to +38% — all at
p < 0.0001. **Across a parameter that reshuffles every phenotype, the effect
barely moves.** This is also the answer to "development never converges, so why
should the phenotype mean anything": the result does not depend on where in the
oscillation the state is sampled.

*(An earlier version of this test appeared to show instability, +23% to +112%.
That run set no sparsity parameters, so networks saturated and "edges" merely
tracked `ng²`. The corrected version is above; the flawed one is recorded in §3
because the failure mode — measuring a proxy instead of the quantity — is worth
stating.)*

*Caveat to carry:* the largest effects come from the mildest sparsity regimes. At
the sparsest setting the magnitude is measured on survivors (see §1.3), so treat
existence as established and magnitude at high extinction as approximate.

### 1.2 The information content of the fluctuation is irrelevant

**15 comparisons pooled: 9 positive, 6 negative, mean −0.171, SE 0.562, t = −0.30.**
Sign test p = 0.61. Two of fifteen reach p < 0.05 against a chance expectation of
0.8 — **and they point in opposite directions** (−3.19 and +2.98).

Established with a **yoked control**: identical world, identical rotation,
identical ecology, with only the cue's mutual information with the environmental
state removed (2.585 bits → ~0.067). A third arm with no cue at all anchors the
baseline, and the scrambled arm is statistically indistinguishable from it.

**It also holds across development time.** The one comparison that had suggested
otherwise (`dev_iters` = 4, +2.98, p = 0.009) does not replicate at 8, 12 or 20
(p = 0.47, 0.73, 0.48). On density the pattern is the same and the isolated hit
is *negative* (`dev_iters` = 8, p = 0.021). Two isolated hits in eight
comparisons, pointing opposite ways, is noise.

| `dev_iters` | Δ`ng` | p |
|---|---|---|
| 4 | +2.98 | **0.009** |
| 8 | +0.90 | 0.47 |
| 12 | +0.54 | 0.73 |
| 20 | +0.70 | 0.48 |

*This is the most robust finding in the work.* Complexity responds to change, not
to *solvable* change.

### 1.3 Complexity is what surviving a changing environment costs

**Frozen environment: 0 of 60 populations extinct. Fluctuating: 37 and 33 of 60.
Fisher exact p = 3.6 × 10⁻¹³.**

Arms differ only in the seasonal period; same grid, sparsity, inflow fraction and
routing. Information-independent (37 vs 33).

This gives 1.1 a mechanism: organisms meet fluctuation by adding wiring, so
taxing wiring removes the response and they die.

### 1.4 Task difficulty does not set the complexity plateau

Four distinct manipulations, all null, with detection floors of 2.21, 3.41 and 4.19 nodes:
zoned environments, an abiotic cue, independent (non-shared) task structure, and
a five-point heterogeneity gradient (H = 2–8). In the one regime where organisms
demonstrably *did* evolve to read the cue (switchers 0.645 vs 0.229, p < 10⁻⁵),
genome size still did not move (p = 0.836, 95% CI [−2.13, +2.64]).

### 1.5 Reversible plasticity is what made the experiment possible

**The one prediction in this work that came out right.** With development fixed
at birth, the regime theory points at — 1–4 generations per environmental state —
is uninhabitable. Reversible re-expression removes exactly that failure:

| period | generations/state | fixed at birth | redevelop 5 | redevelop 25 |
|---|---|---|---|---|
| 600 | 3.8 | **0/3** | **3/3** | **3/3** |
| 1200 | 7.7 | 1/3 | 2/3 | **3/3** |

Predicted in advance from Levins' grain argument — irreversible developmental
plasticity suits coarse-grained environments, fine-grained ones call for
reversible responses — and confirmed. Worth reporting because the four preceding
task manipulations failed for reasons that were only diagnosable afterwards; this
one was diagnosed first and then fixed.

### 1.6 Long-run outcome: the floor, then nothing

8 seeds × 25M ticks ≈ 630,000 generations each. Final `ng`: 13, 14, 12, 12, 14,
12, 12, 12 — six of eight exactly at `NG_MIN` = 12, with every available weight
non-zero. Reached the floor by generation ~14,000; last-quarter slope
**−0.0002 per million ticks** over 6.2M ticks.

### 1.7 Long-run dynamics beyond genome size

Three findings from the completed run that qualify or extend other results:

- **Cross-feeding collapses at this timescale.** `top_share` for channel 1 is
  **0.999** throughout, Shannon diversity 0.006–0.014, waste channels drained
  (res5 20,903 → 1,082). The 0.26–0.30 cross-feeding figure that this model
  produces on 150k-tick runs **does not survive 630,000 generations.** Whatever
  maintains ecotype structure on short timescales does not hold on long ones.
- **A ~70× population swing that is imposed, not emergent.** Period exactly
  15,000 ticks, sd 0 across seeds — because `season_period` is a parameter. Worth
  reporting only for the **amplification**: ±60% inflow → ~70-fold population
  swing.
- **Motility tracks the season genetically, not plastically.** It swings
  0.002 ↔ 0.994 in antiphase with population (r = −0.938), which looks like
  evolved plasticity. Developing the *same* genomes across resource levels
  0.05–4.0 gives a flat 0.94–0.998. Selection cycles the population; nothing
  plastic happens. Consistent with 1.4.

### 1.8 The no-cue baseline that anchors the yoked control

A third arm where nothing varies (`cond_states` = 1), probed identically:

| arm | switchers | spread | matched |
|---|---|---|---|
| no cue | 0.133 | 0.060 | **1.00/6** |
| scrambled | 0.250 | 0.093 | 1.11/6 |
| informative | **0.704** | **0.415** | 1.75/6 |

Scrambled is **statistically indistinguishable from no-cue** (p = 0.26, 0.34) —
correct, since neither carries usable information — while informative departs
from it at p = 0.0004 / 0.0003 / 0.0014. And `matched` for the no-cue arm lands
on **exactly 1.00/6**, the analytically derived non-responsive baseline (an
unresponsive population answers one fixed channel, and each channel is
best-paying in exactly one of six states). The measure is calibrated, not merely
relative.

---

## 2. Novelty: what is new against what reproduces

**Five literature checks were run; four cut the claim back.** Be explicit about
this in the paper — it is more defensible than discovering it in review.

### Reproduces published work
- **Fluctuation raises genome size.** Bentkowski, van Oosterhout & Mock 2015
  (*GBE*): variability → larger genomes, 2 → 15–16 genes, mechanism *relaxed
  selection on metabolic efficiency*. Also Tsuda & Kawata 2010 (*PLOS CB*),
  *Size Matters* 2011.
- **Population size drives genome size down.** Luiselli et al. 2024 (*GBE*,
  Aevol): coding density ~80% at 4× population, ~88% at 16×; compound parameter
  **N×μ** governs.
- **Genome size bounded by mutational robustness, not function.** Knibbe et al.
  2007 (*MBE*): indirect selection for robustness to duplications/deletions
  bounds genome size; Drake's rule emerges spontaneously.
- **"Rewiring not adding."** Established empirically: rewired bacterial networks
  (*Nature* 2008), Ndt80 regulon, echinoid GRNs.
- **Complexity as byproduct.** van Nimwegen's quadratic TF scaling; *GBE* 2026.
- **Environmental grain and plasticity.** Levins; Moran 1992; Walasek et al. 2025.

### Appears to be new
1. **The yoked control applied to this question** — separating a fluctuation's
   *information content* from the fluctuation itself. No precedent found across
   searches of cue reliability, sensing, network reliability, evolvability,
   digital evolution and bet-hedging.
2. **The edge/node decomposition under fluctuation.** Tsuda & Kawata's model
   *could* have measured this — both dimensions evolve there — but they report
   node duplication as primary and never quantify edges separately.
3. **The extinction cost** (1.3), which follows from 2 and has no obvious
   precedent.

### Position against, do not ignore
- **Aevol reaches the same conclusion by a mechanism this model lacks.** Knibbe
  is explicit that it *"would collapse without rearrangements."* This model has
  none — and their **N×μ law was tested here and fails** (matched-product cells
  differ, p = 0.022; a 4× product cell does not, p = 0.77). Frame as convergent
  evidence from different mutational architectures.
- **Tsuda & Kawata emphasise nodes; this finds edges dominate.** Their mutational
  operators differ. Address directly.
- **LaBar & Adami 2016 find large populations evolve *larger* genomes** (U-shape);
  this finds the largest world at the floor. Unresolved; state it.

---

## 3. Caveats that must travel with the numbers

| caveat | detail |
|---|---|
| **Development never converges** | Activations saturate (mean \|e\| 0.846) and ~92% of nodes flip sign every iteration. The phenotype is a snapshot of an oscillation, not a fixed point, and `dev_iters` changes the top uptake channel in 100% of organisms. **But the result is robust to it**: density effect +0.118 to +0.144 across `dev_iters` ∈ {4,8,12,20}, all p < 0.0001 (§1.1). An earlier version of this test appeared unstable (+23% to +112%) because it ran *without* sparsity, so "edges" tracked `ng²` rather than connectivity — a proxy measured in place of the quantity. |
| **Scramble leaked 0.067 bits** | The zero-information control was verified with a spatial gradient but used without one, where an integer multiplier locks to the cycle. 2.585 vs 0.067 is still 39×, and residual signal makes the null *conservative*. Fixed (irrational multiplier); pre-fix runs not bit-reproducible. |
| **Survivorship in sparse regimes** | Frozen 0/60 extinct vs fluctuating 37/60. Effect size *falls* as extinction rises (+122% → +55%), so survivorship does not manufacture it — but magnitudes at high extinction are approximate. |
| **Modularity is not properly testable** | Q is measurable only at low density (0.146 at 0.60 vs 0.06 at saturation), and that regime kills 55–62% of fluctuating populations. Q null has a detection floor of 0.022 (~15%). |
| **`nzw` was a bad measure before sparsity existed** | Point mutations are additive, so density ratchets with lineage depth (0.25 → 0.86) regardless of selection. One earlier inference retracted. |
| **Single model, single architecture** | No independent replication. Aevol is the obvious platform. |
| **Edge count is not independent of node count** | Edges scale as `ng²`. Always report density alongside, and show the excess over the `ng²` baseline (mean +29 points, 5/6 significant). |
| **Arms differ in population size** | Frozen 3,238 vs fluctuating 2,583. Addressed by regression: the fluctuation term survives and population carries the expected negative sign. Report the regression, not just the means. |

---

## 4. Methods summary

- **Model:** individual-based, spatially explicit consumer–resource dish.
  Gene-regulatory-network genome (weight matrix + bias), developed once at birth
  (or periodically with `redevelop_every`), phenotype cached.
- **Evolvable dimensions:** node count (duplication/deletion), connectivity
  (`p_edge_silence` + `maint_edge`), and all phenotype traits.
- **Determinism:** bit-identical results at 1, 4 and 8 threads, verified with the
  full feature set active. Per-child RNG seeds drawn serially.
- **Controls:** yoked (scrambled) cue, no-cue baseline, frozen environment,
  no-op graft, sham-free ablation, and 9 "bit-identical when off" feature checks.
- **Statistics:** permutation tests throughout, validated against null data
  (slightly conservative: 0.033 of nulls below p = 0.05), a known 2-SD effect,
  identical samples, and unequal group sizes.
- **Tests:** 35,297 + 1,304 covering the feature set.

---

## 5. Figures to produce

1. **Fluctuation raises density** — observed edge increase against the `ng²`
   prediction at each sparsity level, with the saturated regime (SWA) shown as
   the case where the excess vanishes. This is the figure that answers the
   "edges scale as `ng²`" objection pre-emptively.
2. **The yoked control** — three arms (no cue / scrambled / informative) on
   cue usage *and* on `ng`, showing behaviour responds and complexity does not.
3. **Extinction cost** — survival by arm at each sparsity level, with the
   0/60 vs 37/60 endpoint.
4. **Long run** — `ng` against generations to 630,000, flat after ~14,000.
5. **Meta-analysis forest plot** — all 15 information comparisons, centred on
   zero.
6. *(supplementary)* Development non-convergence — step size by iteration,
   showing no decay.

---

## 6. Still open

- The LaBar & Adami discrepancy.
- Independent replication, ideally in Aevol.
- Modularity, which needs a regime that is both sparse and survivable.

---

## 7. Pre-registration record

Every substantive prediction was committed to git *before* the run that tested
it, so the history is a usable pre-registration record rather than a
reconstruction. **Five of six were wrong.**

| prediction, made in advance | outcome |
|---|---|
| `ng` rises with task complexity → plateau is task-limited | **wrong** — flat |
| rotating yields alone makes the task conditional | **wrong** — inflow monopoly |
| an abiotic cue will be read because it cannot be depleted | **wrong** — legible but unread |
| the conditional genome is adaptive but mutationally unreachable | **wrong** — selected against |
| reversible plasticity rescues the uninhabitable regime | **correct** — 0/3 → 3/3 |
| conditional response, once evolved, raises genome size | **wrong** — response evolves, complexity does not |

## 8. Data and code availability

- **Code:** the simulator is a stdlib-only Julia package; analysis and plotting
  live in a separate environment so the core has no plotting dependencies.
- **Reproducibility:** every run is determined by its `Params`, including the
  seed. Results are **bit-identical at any thread count** (verified at 1, 4 and
  8). Nine feature flags are bit-identical when off, so historical results
  reproduce under current code — with one documented exception: the scrambled cue
  multiplier was corrected, so pre-fix scrambled runs are not bit-reproducible.
- **Tests:** 35,297 + 1,304, covering determinism, backward compatibility, the
  information content of the yoked control in both gradient modes, and the
  measurement functions.
- **Scale:** ~1,835 runs across 25 experiments, plus 8 × 25M ticks.

## 9. Anticipated reviewer objections

Written out so they can be answered in the paper rather than in a rebuttal.

| objection | answer | where |
|---|---|---|
| "Edges scale as `ng²` — your edge effect is just node growth." | Density is measured directly and rises +0.09–0.14, significant in 5/6, a mean +29 points above the `ng²` baseline. The saturated regime shows −4% excess, exactly as it should. | §1.1 |
| "Your arms differ in population size, and you say population drives genome size." | Regression: fluctuation +4.76 (t = +4.76) on `ng` and +0.084/+0.127 (t ≈ 10) on density; population carries the expected negative sign and is n.s. for density. | §1.1 |
| "Absence of evidence is not evidence of absence — what is your power?" | Detection floors stated for every null: 2.21, 3.41, 4.19 nodes; Q floor 0.022 (~15%); and the information null pooled to t = −0.30 over 15 comparisons. | §1.2, §1.4 |
| "Survivorship — 37/60 of your fluctuating runs died." | Effect size *falls* as extinction rises (+122% → +55%), the opposite of what survivorship would produce. Magnitudes at high extinction reported as approximate. | §1.3, §3 |
| "Your organisms only reach `matched` 1.55/6 — maybe they never solved the task well enough for information to matter." | Granted, and quantified: they exploit ~10.9% of available information, so the realised advantage is ≤0.28 bits against a 5.9-point detection floor. The fitness null is *consistent with* Kussell–Leibler, not evidence against it. | §3 |
| "Development never reaches a fixed point — why does the phenotype mean anything?" | It is deterministic, heritable, and robust to small weight perturbations; evolution acts on a well-defined trait. Direction of every result is `niter`-robust; magnitude is not, and that is stated. | §3 |
| "Aevol already showed this." | Agreed for the population/mutation result, and cited. Their mechanism requires chromosomal rearrangements, which this model lacks, and their N×μ law is tested here and fails. Convergent evidence from different architectures. | §2 |
| "Single model." | Conceded. Independent replication in Aevol named as the priority next step. | §2, §6 |
| "Why these parameter values?" | Every regime was chosen by a documented viability pilot, not tuned to a result; extinction boundaries are reported alongside effects. | §4 |

## 10. Suggested framing

**Venue:** ALIFE conference or a short contribution to *Artificial Life*. Not a
general evolution journal — the science is two modest refinements of established
results.

**Title direction:** something on *what* environmental fluctuation selects for —
connectivity rather than gene count — with the information-irrelevance as the
second claim.

**The honest strength of this paper is its controls**, not its findings. Five of
six pre-registered predictions were wrong; fourteen measurement traps are
documented with the control that caught each; four novelty claims were cut back
by literature checks; three small-sample signals evaporated on replication. A
methods-forward framing is both more defensible and more useful to the field
than presenting the findings as discoveries.

---

## 11. Independent replication, prior art, and the pre-registered test (2026-09-26)

### 11.1 Where the evidence now stands

**R-Aevol, exploratory (v1).** 72 runs: 6 arms × 12 seeds × 30,000 generations,
all completed. Independent platform (sequence → RNA → protein → Hill-function
regulation; Beslon group, `aevol-7` built with `-Dwith-raevol=ON`). Design,
code and output: `replication/raevol/*_v1*`. Per-run means over generations
25,001–30,000; n = 12 per arm.

| arm | genes | TFs | raw links | total strength | strength / link | metabolic error |
|---|---|---|---|---|---|---|
| fluctuating, informative cue | 19.0 | 11.2 | 69.6 | **12.68** | 0.190 | 0.029 |
| fluctuating, yoked cue | 17.9 | 6.2 | 39.6 | 5.66 | 0.141 | 0.058 |
| fluctuating, no cue | 18.9 | 7.4 | 15.0 | **1.42** | 0.099 | 0.063 |
| frozen, informative cue | 26.1 | 13.9 | 80.0 | 4.69 | 0.089 | 0.008 |
| frozen, no cue | 26.7 | 16.6 | 68.8 | 8.21 | 0.129 | 0.008 |
| frozen, yoked cue | 24.6 | 12.1 | 76.8 | 6.17 | 0.101 | 0.009 |

- Genome size: **frozen > fluctuating in every cue condition** (no cue: −7.9
  genes, d = −2.35, p = 4e-5). Opposite sign to §1.
- Raw link counts misled: similar counts carry very different total strength.
  By strength, fluctuation *with a usable cue* beats its frozen counterpart
  2.7× (t = 6.88, p = 2e-5); fluctuation *without* one gives the weakest
  networks of all six arms.
- Information matters: informative vs yoked, total strength 2.2× (p = 2e-5),
  metabolic error halved — the organisms demonstrably use the cue. Contradicts
  §1.2, which was measured in a model whose organisms never read the cue.
- Frozen informative vs frozen yoked: complete null (p = 0.30–0.90) — a real
  control here, because the signal still switches under a fixed target.

These numbers are **exploratory**: the trackable/untrackable interpretation was
formed after seeing them, so they cannot also confirm it. §11.4 is the test.

**EvoSim, strength vs count (exploratory).** 180 runs (6 arms × 30 seeds,
150k ticks, grid 128, `p_edge_silence` 0.2, `maint_edge` 5e-5). grid 48 was
abandoned after 6/6 fluctuating runs and 0/6 frozen runs went extinct. The
three frozen arms are bit-identical (with one state the cue is zero whatever
its settings), so EvoSim has one frozen arm and no frozen yoked control.
`mean_abs_w` is Σ|w| / ng² — a strength *density*, the direct counterpart of
count density nzw / ng². Interim at n ≈ 13 per arm (to be replaced by n = 30):

| measure, fluctuating (no cue) vs frozen | change | p |
|---|---|---|
| count density nzw / ng² | **+28%** | 2e-5 |
| strength density Σ\|w\| / ng² | **−13%** | 5e-5 |
| total strength Σ\|w\| | −10% | 0.17 |
| node count | +1% (at the floor, ~14 vs minimum 13) | 0.65 |

§1.1's count-density effect replicates (+0.129 here vs +0.118–0.144 before) but
wiring strength moves the other way: more edges, each weaker. Population is 35%
smaller under fluctuation and the ranges do not overlap between arms; within
arms, population tracks count density (r = −0.55) but not strength density in
the direction needed to explain it.

### 11.2 Prior-art check, sixth (2026-09-26)

| claim | status | source |
|---|---|---|
| Cue reliability decides whether regulatory plasticity evolves in a fluctuating environment | **Established** | Kussell & Leibler 2005 (*Science*); Botero et al. 2015 (*PNAS*); Jorritsma & van den Berg 2026 (*Evolution* 80:230) — a GRN model: plasticity evolves under fast, erratic change "especially if cues are reliable" |
| Fluctuating environments can select *for or against* gene regulation | **Established, experimentally** | Pollack, Nozoe & Kussell 2025 (bioRxiv 2025.06.23.661119) — determinant there is metabolic flux, not cue information |
| Fluctuation without cues *promotes* GRN complexity | **Published, conflicts with v1** | Tsuda & Kawata 2010 (*PLOS CB*): random-walk optimum promotes complex GRNs via duplication. Their fluctuation is across generations and tracked genetically (coarse-grained); ours switches within a lifetime (fine-grained). Grain is the likely reconciling variable — must be addressed, not ignored |
| Sensing and non-sensing populations diverge genomically under fluctuation | **Published, different measures** | Lalejini, Ferguson, Grant & Ofria 2021 (*Front. Ecol. Evol.*, Avida): trait loss and deleterious accumulation in non-plastic populations; no network architecture |
| Weighted network size (Σ edge strength) as a GRN variable | **Standard method** | Iwasaki, Tsuda & Kawata 2013 (*BMC Evol Biol*) |
| Neutral accumulation of interactions inflates counts | **Established** | Soyer & Bonhoeffer 2006 (*PNAS*) |

**Not found:** (a) edge count and edge strength moving in *opposite* directions
under fluctuation; (b) a stimulus-matched (yoked) control separating a cue's
information from its presence in GRN-architecture evolution; (c) total
regulatory strength and genome size contrasted across informative, yoked and
absent cues.

**Unread, and must be read before claiming (c):** Jorritsma & van den Berg 2026
full text (paywalled). Their abstract says they "zoom in on network
architecture"; if they report connection strengths by cue reliability, (c) is
theirs.

### 11.3 The revised claim

> In a fine-grained fluctuating environment, whether selection strengthens or
> weakens regulatory wiring depends on whether a usable cue exists — and edge
> counts do not track wiring strength, so count-based connectivity claims can
> have the wrong sign.

The phenotype-level half (cue reliability → plasticity) is established and is
cited, not claimed. What is potentially new is the architectural consequence,
the yoked control that isolates information from stimulation, and the
count/strength divergence.

### 11.4 Pre-registration: R-Aevol sweep (v2)

Committed before any v2 run. Code: `replication/raevol/gen_params_v2.py`
(design), `replication/raevol/analyse_v2.py` (measures, tests and decision
rules, fixed). Anything not specified there is exploratory and will be
labelled as such.

**Design.** Identical to v1 except switch rate, population and seeds.
- Primary: switch rate p ∈ {0.05, 0.10, 0.20} per age, population 1024;
  arms fluc_info, fluc_yoked, fluc_none, froz_info at every rate;
  froz_none and froz_yoked at p = 0.10 (froz_none has identical targets and no
  signals, so it is rate-independent). 168 runs.
- Secondary: population 512, p = 0.10, the four fluctuating/frozen-info arms.
  48 runs.
- 12 seeds per cell, all new (base 200000; v1 used 100000). 30,000 generations.
- Four-environment (yoked) arms use 1.5 × p, which equalises per-age target-
  and signal-change rates with the two-environment arms.
- Launched seed-major, so a partial batch is still balanced.

**Measures.** Per run, mean over generations 25,001–30,000 of the best
individual. Total regulatory strength S = nb_link × avg_influence, **counted
as 0 when nb_link = 0** (v1's script dropped those generations, which biases
sparse arms upward).

**Predictions** (each at each of the three rates, population 1024):

| id | prediction | meaning |
|---|---|---|
| P1 | S(fluc_info) > S(froz_info) | trackable fluctuation strengthens regulation |
| P2 | S(fluc_info) > S(fluc_yoked) | the cue's information, stimuli matched, is what does it |
| P3 | S(fluc_none) < S(froz_none) | untrackable fluctuation weakens regulation |
| P4 | genes(froz_info) > genes(fluc_info) | fine-grained fluctuation shrinks the genome |

**Test.** Two-sided permutation test on the difference of means, 100,000
permutations; Holm correction across all 12 tests.

**Manipulation check.** At each rate, metabolic error lower in fluc_info than
fluc_yoked (p < 0.05). P1 and P2 are evaluable only at rates where it holds —
where the cue is not used they do not test the hypothesis.

**Verdicts.** SUPPORTED: predicted sign at every evaluable rate and
Holm-adjusted p < 0.05 at ≥ 2. REFUTED: opposite sign with p < 0.05 at any
evaluable rate, or no Holm-significant rate. INCONCLUSIVE: otherwise, or
P1/P2 with fewer than two evaluable rates.

**What each outcome means for the paper.**
- P1–P3 supported → the revised claim stands on a pre-registered test across
  a threefold range of switch rates.
- P2 refuted → information is not what matters; the v1 difference came from
  something else in the 2-vs-4-environment construction, despite the frozen
  control.
- P1 or P3 rate-dependent → report the boundary; the claim becomes
  conditional on switch rate.
- P4 refuted → the genome-size sign flip in v1 was not robust, and the
  precision-is-gene-limited mechanism has nothing to explain.

**Secondary, no verdicts.** Raw link count and per-link strength for the P1
contrast; frozen informative vs frozen yoked (expected null); P1–P4 at
population 512; rank-based p-values beside every primary test.

**Expected duration.** ~22 h on the Mac mini at 10 concurrent runs (v1: 72
runs in 8.2 h), primary cells first.

### 11.5 Pre-registration: EvoSim confirmation (P5)

The EvoSim batch in §11.1 is exploratory — the strength measure was chosen, and
an interim look taken, while it ran. Confirmation in a second regime, where node
count sits above the floor:

- **Regime.** Same settings as §11.1 except grid. Grid chosen by a viability
  pilot before any confirmatory run: the smallest of {80, 96, 112} at which
  6/6 fluctuating (no cue) and 6/6 frozen pilot runs survive 150k ticks. Pilot
  seeds are not reused.
- **Arms.** fluc_none vs froz_none, 30 seeds each, new seeds.
- **Measures.** Per run, mean over the final 20% of ticks: count density
  nzw / ng², strength density `mean_abs_w`.
- **Predictions.** P5a: count density higher under fluctuation.
  P5b: strength density lower under fluctuation.
- **Test.** Two-sided permutation, 100,000 permutations, Holm across the two.
  Both significant in the predicted direction → SUPPORTED; either significant
  in the opposite direction → REFUTED; otherwise INCONCLUSIVE.
- **Exclusions.** Extinct runs excluded and counted; if survival differs
  between arms by more than 10 percentage points, the result is reported as
  survivorship-confounded rather than given a verdict.
- **Scheduling.** After the R-Aevol sweep (the Mac mini is fully used by it).

### 11.6 Deviations from the pre-registration

None yet. Any change to design, measures, tests or decision rules after
launch is logged here with its reason.
