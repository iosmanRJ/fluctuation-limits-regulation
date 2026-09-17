# Environmental fluctuation selects for connection density, and does not care whether it is predictable

**Draft v1.** Target: ALIFE conference / *Artificial Life* short paper.

---

## Abstract

Gene regulatory network complexity is known to increase under fluctuating
selection, but "complexity" in such studies is usually a single quantity — the
number of genes — and the *information* the fluctuation carries is rarely
separated from the fluctuation itself. Using an individual-based
consumer–resource model in which both network size and connectivity evolve
freely, we show that environmental fluctuation raises **connection density**
beyond what the increase in gene number alone would produce (+0.12 to +0.14 in
density, *p* < 0.0001, robust across six sparsity regimes and four
development-time settings), and that this occurs **regardless of whether the
fluctuation carries usable information about the environment**. Using a yoked
control that holds the world fixed while destroying only the cue's mutual
information with the environmental state, we find no effect of information on
network size across 15 independent comparisons (*t* = −0.30, sign test
*p* = 0.61), even in populations that demonstrably evolved to read the cue
(switchers 0.645 vs 0.229, *p* < 10⁻⁵). When connectivity is made costly, the
same fluctuation becomes lethal — 0 of 60 populations die in a frozen
environment against 37 of 60 in a fluctuating one (*p* = 3.6 × 10⁻¹³) — placing
the added connectivity as the price of persistence rather than an incidental
byproduct. The fluctuation effect reproduces established results; the
contribution is the decomposition into size and density, and the demonstration
that predictability is irrelevant to both.

---

## 1. Introduction

That fluctuating environments select for larger genomes and more complex
regulatory networks is well established, in simulation and in comparative
genomics. Bentkowski et al. (2015) find environmental variability driving genome
expansion from 2 to 15–16 genes; Tsuda and Kawata (2010) find complex GRNs
promoted by duplication under unpredictable fluctuation; Luiselli et al. (2024)
find population size driving genome size down through indirect selection for
replicative robustness. Empirically, terrestrial bacteria — in environments
varying two to three orders of magnitude more than marine ones — carry larger
genomes disproportionately enriched in regulatory genes.

Two things are usually left unresolved.

**First, "complexity" is measured as one number.** Most models either fix network
size and evolve connectivity, or evolve size and report it alone. Where both can
move, the decomposition is rarely reported: Tsuda and Kawata's model permits both
and emphasises gene duplication, but does not quantify edge density separately.
Yet the two are not independent — edges scale as the square of nodes — so an
apparent connectivity effect can be nothing more than the arithmetic consequence
of network growth.

**Second, fluctuation and information are confounded.** A fluctuating environment
that an organism can *sense* offers two distinct things: variation, and
predictive information about that variation. Theory says information has
quantifiable fitness value (Kussell and Leibler, 2005; Rivoire and Leibler,
2011), and plasticity theory says cue reliability governs whether plastic
responses evolve (Moran, 1992). But in a simulation where the cue is derived from
the environment, the two cannot be separated by observation alone.

We address both. The model evolves node count and connectivity independently, and
a **yoked control** — matched stimuli, different contingency — holds the
environment, the ecology and the population fixed while removing only the cue's
information content.

## 2. Model

An individual-based, spatially explicit consumer–resource dish. Organisms carry a
gene-regulatory-network genome (a weight matrix plus per-node bias) that develops
into a phenotype at birth, optionally re-expressed during life. The phenotype
allocates a fixed uptake budget across six resource channels and sets size,
motility, division threshold, repair fidelity, stress tolerance and leakage.
Resources arrive by abiotic inflow, diffuse, decay, and are converted along a
fixed excretion ladder that makes one organism's waste another's substrate.

**Evolvable dimensions.** Node count changes by duplication and deletion.
Connectivity changes by point mutation and by an edge-silencing operator; a
per-edge maintenance cost makes sparsity payable. Without both a silencing
operator *and* a cost, connectivity cannot evolve: point mutations are additive,
so a weight that becomes non-zero never returns to zero, and density ratchets
from 0.25 at founding to 0.86 by lineage depth 1,400 irrespective of selection.

**The environment.** An abiotic condition cycles through six states, each
rotating which resource channel pays best and receives inflow. The condition is
not consumable, so unlike resource concentration it remains legible however hard
the population exploits the dish. It reaches the genome through a dedicated
sensor node with no phenotypic readout.

**Arms.** *Frozen*: the rotation machinery runs but the state never changes.
*Informative*: the cue predicts the state (2.585 bits, the full entropy of a
six-state environment). *Scrambled*: identical world, identical rotation,
identical ecology, with the cue driven from an incommensurate phase so its mutual
information with the state is ~0.067 bits. A fourth arm with no cue at all
anchors the measurement scale.

**Reproducibility.** Every run is determined by its parameters including seed, and
results are bit-identical at any thread count.

## 3. Results

### 3.1 Fluctuation raises density, not merely size

At the reference regime (20 seeds per arm):

| arm | *n* | nodes | edges | density |
|---|---|---|---|---|
| frozen | 20 | 18.01 | 201.3 | 0.592 |
| informative | 18 | 24.26 | 436.3 | 0.717 |
| scrambled | 17 | 23.36 | 415.2 | 0.744 |

Fluctuation raises node count by +5.81 (*p* < 0.0001) and density by +0.138
(*p* < 0.0001). The node increase alone would predict ~+75% more edges; the
observed increase is +112%, an excess of **+37 percentage points**. Across six
sparsity regimes the density rise is significant in five, with a mean excess of
+29 points.

The sixth is the informative exception: in the saturated regime (density 0.987),
where density *cannot* rise, the excess is −4% — the edge increase is exactly
what node growth predicts and nothing more. Where density can move it does;
where it cannot, the effect reduces to arithmetic.

The effect is not a population artefact. Frozen populations are larger (3,238 vs
2,583), so with population entered alongside: fluctuation +4.76 (*t* = +4.76) on
node count and +0.084 to +0.127 (*t* ≈ 10) on density, with population carrying
the expected negative sign and going non-significant for density.

It is also robust to development time, which in this model controls the phenotype
completely — changing the iteration count reassigns the top uptake channel in
100% of organisms. Across `dev_iters` ∈ {4, 8, 12, 20} the density effect spans
only +0.118 to +0.144, every cell at *p* < 0.0001.

### 3.2 The information content is irrelevant

Organisms in the informative arm demonstrably evolve to read the cue: switchers
0.645 against 0.229 scrambled (*p* < 10⁻⁵), profile spread 0.374 against 0.107,
and matched responses 1.55/6 against 0.98/6. A no-cue arm gives switchers 0.133
and matched **exactly 1.00/6**, the analytically derived non-responsive baseline,
and the scrambled arm is statistically indistinguishable from it.

Network size does not follow. Pooling all 15 informative-versus-scrambled
comparisons: 9 positive, 6 negative, mean −0.171 (SE 0.562, *t* = −0.30), sign
test *p* = 0.61. Two reach *p* < 0.05 against a chance expectation of 0.8, and
they point in opposite directions. The null holds across three parameter regimes,
a five-point heterogeneity gradient, a structure manipulation, six sparsity
levels and four development times, with detection floors of 2.21 to 4.19 nodes.

**Behavioural complexity and genomic complexity dissociate.** The environment
selects for regulatory capacity in response to *change*; it does not select for
more of it when the change is *predictable*.

### 3.3 The added connectivity is the price of persistence

With connectivity taxed, the same fluctuation becomes lethal. Across 60 seeds per
arm: **0 of 60 frozen populations go extinct, against 37 of 60 informative and 33
of 60 scrambled** (Fisher exact *p* = 3.6 × 10⁻¹³). The arms differ only in
whether the state changes.

This gives §3.1 a mechanism. Organisms meet fluctuation by adding wiring; price
the wiring and they cannot mount the response, so they die. The extinction is
information-independent (37 vs 33), consistent with §3.2.

### 3.4 Task difficulty does not substitute for fluctuation

Four manipulations designed to make the environment *demand* regulatory
complexity — spatial zoning, an abiotic conditional cue, task states sharing no
structure, and a five-point gradient in the number of distinct states — all
returned null effects on network size, at detection floors of 2.21 to 4.19 nodes.
In the one regime where organisms demonstrably solved the conditional task, size
still did not move (*p* = 0.836, 95% CI [−2.13, +2.64]).

## 4. Discussion

**What reproduces.** Fluctuation raising genome size reproduces Bentkowski et al.
(2015) and Tsuda and Kawata (2010). Population size driving it down reproduces
Luiselli et al. (2024). Genome size bounded by mutational rather than functional
pressure reproduces Knibbe et al. (2007). Adaptation by rewiring rather than gene
addition is established empirically. We claim none of these.

**What is new.** The decomposition of the fluctuation response into size and
density, with the `ng²` baseline subtracted; and the demonstration, via a yoked
control, that neither responds to the cue's information content.

**Relation to Aevol.** Knibbe et al. reach our conclusion about functional demand
by a mechanism — indirect selection for robustness to chromosomal rearrangements
— that this model lacks entirely, and which they state "would collapse without
rearrangements." Consistent with that, their compound parameter N×μ does not
govern here: cells matched on N×μ differ (*p* = 0.022) while a cell at four times
the product does not (*p* = 0.77). Two models with different mutational
architectures reaching the same conclusion is convergent evidence, not
duplication.

**Why information does not help.** Organisms exploit only ~10.9% of the available
information (matched 1.55 of 6 against a baseline of 1), so the realised fitness
advantage is bounded at ≤0.28 bits per generation — far below what this design
could detect (5.9 percentage points in extinction rate). Our fitness null is
therefore *consistent with* Kussell and Leibler rather than evidence against them.

## 5. Limitations

Single model, single architecture; independent replication, ideally in Aevol, is
the priority. Development does not reach a fixed point — activations saturate and
~92% of nodes flip sign each iteration — so the phenotype is a snapshot of an
oscillation; results are robust to where it is sampled but the mechanism is not
the "settling to a stable expression state" the design language implies. In the
sparse regimes the fluctuating arms lose 55–62% of runs, so magnitudes there are
measured on survivors, though effect size *falls* as extinction rises, the
opposite of what survivorship would produce. Modularity could not be tested: the
regime sparse enough for it to be measurable is the regime where fluctuating
populations die.

## 6. Conclusion

Environmental fluctuation selects for regulatory connectivity, not merely for
more genes, and does so whether or not the fluctuation can be predicted. When
connectivity is made expensive the same fluctuation kills most populations,
identifying the added wiring as the cost of persisting through change rather than
an incidental accumulation. Complexity here responds to change, not to *solvable*
change.

## Data and code availability

~2,000 runs across 25 experiments, plus 8 × 25 × 10⁶ ticks (≈630,000 generations
per seed). Every run reproducible from its parameters and seed; bit-identical at
any thread count; 36,601 automated tests covering determinism, backward
compatibility and the measurement functions. All predictions were committed to
version control before the runs that tested them — five of six were wrong.
