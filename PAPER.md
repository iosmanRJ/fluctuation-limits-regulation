# Paper preparation

Everything needed to write up this work: the claims with their support, what is
new against what reproduces published results, the caveats that must travel with
each number, and what is still open.

**Status of the evidence base:** ~1,835 simulation runs across 25 experiments,
plus a completed 8-seed × 25M-tick long run (~630,000 generations per seed).
Two independent code reviews. One defect found and fixed.

---

## 1. What to claim, ranked by strength

### 1.1 Fluctuation raises complexity, and acts far more on connectivity than on gene count

The central result. **7 of 7 independent comparisons positive, every one significant.**

| measure | comparisons | mean effect | range |
|---|---|---|---|
| node count `ng` | 7/7 positive | **+6.32** | +2.64 to +10.15 |
| edge count | 6/6 positive | **+104%** | +63% to +122% |

Nodes rise ~31%; edges **more than double**. Holds across edge densities from
0.99 to 0.60, across two development-time settings, and across separate
experiments.

*Caveat to carry:* the largest effects come from the mildest sparsity regimes.
At the sparsest setting the magnitude is measured on survivors (see 1.3), so
treat existence as established and magnitude at high extinction as approximate.

### 1.2 The information content of the fluctuation is irrelevant

**15 comparisons pooled: 9 positive, 6 negative, mean −0.171, SE 0.562, t = −0.30.**
Sign test p = 0.61. Two of fifteen reach p < 0.05 against a chance expectation of
0.8 — **and they point in opposite directions** (−3.19 and +2.98).

Established with a **yoked control**: identical world, identical rotation,
identical ecology, with only the cue's mutual information with the environmental
state removed (2.585 bits → ~0.067). A third arm with no cue at all anchors the
baseline, and the scrambled arm is statistically indistinguishable from it.

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

Four distinct manipulations, all null, with detection floors of 2–4 nodes:
zoned environments, an abiotic cue, independent (non-shared) task structure, and
a five-point heterogeneity gradient (H = 2–8). In the one regime where organisms
demonstrably *did* evolve to read the cue (switchers 0.645 vs 0.229, p < 10⁻⁵),
genome size still did not move (p = 0.836, 95% CI [−2.13, +2.64]).

### 1.5 Long-run outcome: the floor, then nothing

8 seeds × 25M ticks ≈ 630,000 generations each. Final `ng`: 13, 14, 12, 12, 14,
12, 12, 12 — six of eight exactly at `NG_MIN` = 12, with every available weight
non-zero. Reached the floor by generation ~14,000; last-quarter slope
**−0.0002 per million ticks** over 6.2M ticks.

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
| **Development never converges** | Activations saturate (mean \|e\| 0.846) and ~92% of nodes flip sign every iteration. The phenotype is a snapshot of an oscillation, not a fixed point. `dev_iters` changes the top uptake channel in 100% of organisms. Direction of results is `niter`-robust; magnitude is not (edges +23% to +112% across `niter` ∈ {4,8,12,20}). |
| **Scramble leaked 0.067 bits** | The zero-information control was verified with a spatial gradient but used without one, where an integer multiplier locks to the cycle. 2.585 vs 0.067 is still 39×, and residual signal makes the null *conservative*. Fixed (irrational multiplier); pre-fix runs not bit-reproducible. |
| **Survivorship in sparse regimes** | Frozen 0/60 extinct vs fluctuating 37/60. Effect size *falls* as extinction rises (+122% → +55%), so survivorship does not manufacture it — but magnitudes at high extinction are approximate. |
| **Modularity is not properly testable** | Q is measurable only at low density (0.146 at 0.60 vs 0.06 at saturation), and that regime kills 55–62% of fluctuating populations. Q null has a detection floor of 0.022 (~15%). |
| **`nzw` was a bad measure before sparsity existed** | Point mutations are additive, so density ratchets with lineage depth (0.25 → 0.86) regardless of selection. One earlier inference retracted. |
| **Single model, single architecture** | No independent replication. Aevol is the obvious platform. |

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

1. **Fluctuation effect, nodes vs edges** — paired bars across the sparsity
   sweep, showing edges responding ~3× more.
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

- `niter` 12 and 20 in the sparse regime (running at time of writing) — decides
  whether the edge decomposition is development-time-robust.
- The LaBar & Adami discrepancy.
- Independent replication, ideally in Aevol.
- Modularity, which needs a regime that is both sparse and survivable.

---

## 7. Suggested framing

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
