# Environmental fluctuation strengthens gene regulation only when it can be sensed

**Short title:** Sensed fluctuation strengthens gene regulation

**Target:** PLOS Computational Biology, Research Article.
**Status:** second draft, 2026-09-27. Replaces the first draft (git `9e00da9`),
whose two central claims did not survive independent replication. Sections
marked ⟨PENDING⟩ wait on the grain experiment (pre-registered, running; results
expected 2026-09-28). References marked [verify] need checking against the
source before submission. Numbers are taken from `PAPER.md` §11 and the result
files under `replication/`.

---

## Abstract

Environmental variability is widely expected to favour more elaborate gene
regulation, yet models disagree even on the direction of the effect, and most
measure regulatory complexity by counting network connections. We asked what
fluctuation does to regulatory architecture, separating the fluctuation itself
from the information an organism has with which to track it. In two independent
individual-based models of genome evolution — EvoSim, whose genomes are
regulatory networks of evolvable size, and R-Aevol, in which regulation emerges
from sequence-encoded transcription factors — we compared frozen and
fluctuating environments under an informative cue, no cue, and a yoked control
that delivers the same signal molecules at the same frequencies but decoupled
from the environment. In EvoSim, fluctuation raised the density of regulatory
connections by 26% while lowering wiring strength per possible connection by
11%, leaving total wiring strength unchanged. In R-Aevol, a pre-registered test
across three switching rates (two for the cue-dependent predictions) supported
all four predictions: with an
informative cue, fluctuation increased total regulatory strength 1.8–2.7-fold
over a frozen environment; the effect required the cue's information rather
than its mere presence (informative versus yoked, 2.3–2.9-fold); without a cue,
fluctuation produced the weakest regulation of any condition; and fluctuation
reduced genome size by 8–10 genes. The pattern replicated at a second
population size. In both models, raw connection counts moved in the opposite
direction to regulatory strength. ⟨PENDING: one sentence on grain.⟩
Environmental fluctuation therefore does not select for regulatory complexity as
such: it strengthens regulation when a cue makes the environment trackable and
weakens it when none does, and counting connections can report the wrong sign.

*(≈265 words; limit 300.)*

## Author summary

Organisms that live in changing environments often carry elaborate gene
regulation, and it is tempting to conclude that change itself builds regulatory
complexity. We tested this in two very different computer models of evolving
genomes. In both, we let populations evolve in environments that were either
constant or switching, and we controlled what the organisms could know about
the switching: some received a signal that announced which environment they
were in, some received the same signal scrambled so that it announced nothing,
and some received no signal at all. We found that change alone did not build
stronger regulation. When the environment switched and a meaningful signal was
available, regulatory connections became much stronger; when it switched with
no usable signal, they became the weakest we saw; and a scrambled signal did
not help, even though it was physically identical. We also found that simply
counting regulatory connections, a common way to measure network complexity,
gave the opposite answer to measuring how strong those connections were. We
tested our main predictions after writing them down in advance.

*(≈177 words; limit 150–200.)*

---

## Introduction

Organisms in variable environments are often thought to need, and therefore to
evolve, more elaborate gene regulation. The intuition is reasonable: regulation
is how a single genome produces different phenotypes in different conditions,
and a constant environment offers nothing for it to do. It has also received
support from comparative genomics, where the number of transcription factors
scales with genome size and with the range of environments a lineage
encounters [refs].

Theory is precise about one level of this question and silent about another.
At the level of the phenotype, whether variability favours a plastic response
depends on grain and on information. A response to the environment within a
lifetime is favoured when the environment changes on the timescale of an
individual's life and a reliable cue predicts the state it will face
[1–4]; without such a cue, fluctuation favours fixed generalist phenotypes,
bet-hedging, or genetic tracking instead [3,4]. A recent gene-regulatory-network
model confirms this for the evolution of plasticity itself: plasticity evolves
mostly under fast, erratic change, especially when cues are reliable, while
evolvability evolves under intermediate variability and less reliable cues [5].

What this means for regulatory *architecture* — how many genes a genome carries,
and how strongly they regulate one another — is much less clear, and models
disagree. Tsuda and Kawata found that fluctuating selection without any
environmental cue [verify: their networks respond to external stimuli; confirm
these carry no information about the optimum] promotes complex regulatory networks, largely through the
fixation of beneficial gene duplications [6]; related models report larger
genomes in fluctuating environments [7]. Work on the Aevol platform instead
finds genome size bounded by selection for mutational robustness and reduced by
large population size and high mutation rate [8,9], and models of modularity
find that its evolution depends on sparseness and on the structure of the
variation rather than on variability as such [10,11].

Two methodological problems stand between these results and a general answer.
First, fluctuation and information are usually confounded. An environment that
changes more often also offers a different amount of predictable structure, and
adding a cue changes both what an organism can know and what it experiences. The
standard remedy from the study of learning is the yoked control [12]: a second
subject receives exactly the same stimuli as the first, but with the contingency
between stimulus and outcome broken. We are not aware of it having been used in
the evolution of regulatory networks. Second, regulatory complexity is usually
measured by counting connections [6,10,11]. Counts can grow by near-neutral
accumulation [13] and say nothing about how strongly a connection regulates.
Weighted measures exist [14] but are rarely the primary readout.

Here we address both problems in two independent models of genome evolution
that share nothing but the question. EvoSim represents a genome directly as a
regulatory network whose size and connectivity evolve; R-Aevol derives
regulation from transcription factors encoded in a DNA-like sequence and binding
to promoter regions [15,16]. In each we compared a frozen environment with a
fluctuating one, and within fluctuation we compared an informative cue, no cue,
and a yoked cue matched to the informative one in signal identity, frequency and
switching rate. We measured regulatory connectivity both as a count and as a
strength. After an exploratory replication in R-Aevol overturned two of our
earlier conclusions from EvoSim, we stated four predictions in advance and
tested them across three switching rates and two population sizes, and we
tested a fifth, about grain, to address the conflict with Tsuda and Kawata.

We find that fluctuation strengthens regulation only when a cue makes the
environment trackable, weakens it when none does, and that the cue's
information — not its presence — is what matters. Fine-grained fluctuation also
shrinks the genome. Throughout, counting connections gives the opposite answer
to measuring their strength.

---

## Results

### Two models and one design

**EvoSim** is an individual-based model of single-celled organisms on a
two-dimensional lattice of diffusing, cross-fed resources (Methods; full ODD
description in S1 Text). Each genome is a gene regulatory network — a weight
matrix and node biases — of 13 to 32 nodes. Thirteen nodes form a fixed
interface: twelve are read out as the phenotype (resource-uptake allocation,
body size, motility, division threshold, DNA repair, stress tolerance and
leakage) and one senses an abiotic cue. The remaining hidden nodes are gained by
duplication and lost by deletion, so genome size evolves. The network is
developed once at birth and the phenotype is fixed for life. Connectivity evolves
because mutation can both create and remove connections and each non-zero
connection carries a small maintenance cost.

**R-Aevol** is the regulation-enabled version of the Aevol platform [15–17].
Genomes are circular sequences in which promoters, genes and binding sites are
found by motif search; genes are translated into proteins that both contribute
to a phenotype (a function matched against a target function) and, by binding
upstream of promoters, enhance or inhibit transcription of other genes.
Protein concentrations are integrated over 20 ages of an individual's life, and
fitness is evaluated against the target at every age. Environmental signals are
signal proteins injected at a fixed concentration during the ages in which their
environment is active, so an organism can evolve to read them.

**The design** (Fig 1) crosses environment (frozen or fluctuating) with cue
(informative, none, or yoked). In R-Aevol the fluctuating environment switches
between two target functions, A and B, mirror images of one another with
identical geometric area, so neither is intrinsically harder. The frozen
environment keeps the switching machinery but makes both targets identical.
The informative arm pairs each target with its own signal protein; the yoked
arm uses four environments — each target with each signal — so that the same
two signals occur at the same frequencies but carry no information about the
target. Because a switch in R-Aevol moves to a different environment chosen at
random, a four-environment arm switches its target less often than a
two-environment arm at the same switching probability; we raised the yoked
arm's probability by a factor of 1.5, which equalises both the rate at which the
target changes and the rate at which the signal changes (Methods). EvoSim uses
the same logic with an abiotic cue whose yoked version is decoupled from the
environment by an incommensurate period.

> **Fig 1. Design.** (A) The two models: EvoSim's network genome and R-Aevol's
> sequence-encoded regulation. (B) The yoked control in R-Aevol: informative arm
> (two environments, each target with its own signal, switching probability p)
> and yoked arm (four environments, each target with each signal, 1.5p), with
> the resulting target- and signal-change rates shown equal. (C) Environmental
> grain as mean dwell time of an environment, from half a lifetime to 250
> generations.

### In EvoSim, fluctuation adds connections but not wiring strength

We first measured connectivity in EvoSim, where earlier work had shown that
fluctuation raises connection density (S4 Text). In 180 populations (six arms ×
30 seeds, 150,000 time steps), fluctuation without a cue raised count density —
the fraction of possible connections that are non-zero — from 0.467 to 0.587
(+25.5%, p < 10⁻⁴, d = 2.2; survivors only, n = 28 versus 30; Fig 2A). Over the
same networks, strength density — summed absolute connection weight divided by
the number of possible connections — fell from 0.207 to 0.184 (−10.9%,
p < 10⁻⁴, d = −1.4; Fig 2B). Total wiring strength did not change (38.0 versus
37.5, p = 0.86). Fluctuation added connections without adding strength: the mean
strength of a non-zero connection fell by about 30%. The same pattern held with an
informative cue and with a yoked one. Node count was near the minimum of 13 in
this regime and rose only slightly (+0.6 nodes, p = 0.047).

Two features of this regime need stating. First, fluctuating populations were
40% smaller than frozen ones, and population size can affect genome
architecture on its own. The population ranges of the two arms did not overlap,
so the two cannot be separated statistically; but within every arm, larger
populations had *lower* strength density (r = −0.39), so population size, taken
at face value, predicts higher strength density in the smaller fluctuating
populations — the opposite of what we observed. Second, in this regime the cue
went unused: informative and yoked arms had identical strength density
(p = 0.99). EvoSim organisms read the cue only when development is made
reversible, and even then exploit little of its information (S4 Text), so
EvoSim speaks mainly to fluctuation that cannot be tracked.

> **Fig 2. EvoSim: more connections, less strength.** Count density (A) and
> strength density (B) for fluctuating (no cue) and frozen populations, one
> point per population, with means and 95% CIs. (C) Total wiring strength.
> (D) Within-arm relation between population size and strength density.

### In R-Aevol, whether fluctuation strengthens regulation depends on a usable cue

We next asked whether an independent model reproduced EvoSim's earlier
conclusions — that fluctuation raises genome size and connectivity, and that
the information content of a cue is irrelevant. It did not. In an exploratory
experiment (six arms × 12 seeds, 30,000 generations), fluctuation *reduced*
genome size in every cue condition, and whether it strengthened regulation
depended on the cue (Fig 3). With an informative cue, fluctuating populations
evolved total regulatory strength of 12.7 against 4.7 in the frozen arm with the
same signals; with no cue, fluctuating populations had the weakest regulation of
all six arms (1.4).

Because this interpretation was formed after seeing the data, we tested it with
new seeds. Before running anything, we committed four predictions, their
measures, statistical tests and decision rules (Table 2; Methods). The test
varied the switching probability over a fourfold range (0.05, 0.10 and 0.20 per
age) at a population of 1024, with 12 new seeds per cell, and repeated the
central contrasts at a population of 512. For the two predictions that depend on
cue use, we also specified a manipulation check: that the informative cue
lowers metabolic error relative to the yoked one, i.e. that organisms actually
use it. The check held at switching probabilities 0.10 and 0.20 (p < 10⁻⁴) but
not at 0.05 (p = 0.088), so by the pre-registered rule those predictions were
judged at the two faster rates only.

All four predictions were supported (Table 2, Fig 4). With an informative cue,
fluctuation increased total regulatory strength over the frozen environment
1.8-fold at p = 0.10 and 2.7-fold at p = 0.20 (P1; Holm-adjusted p = 0.005 and
≈ 10⁻⁴). Without a cue, fluctuation reduced it roughly fivefold at every rate
(1.1–1.2 against 6.3; P3; Holm-adjusted p ≈ 10⁻⁴ at each rate). The
slowest-switching informative cells, excluded from the verdict, pointed the same
way (11.0 against 5.9). One of these populations lost every regulatory
connection for the whole analysis window, the only such case in the experiment.

> **Fig 3. R-Aevol, six arms (exploratory).** Total regulatory strength (A), raw
> link count (B) and gene number (C) for each environment × cue combination,
> one point per population.

### The cue's information, not its presence, is what matters

The yoked control isolates information from stimulation: the informative and
yoked arms receive the same two signal proteins, at the same concentration and
frequency, switching at the same rate; only the contingency between signal and
target differs. Informative-cue populations evolved 2.3-fold (p = 0.10) and
2.9-fold (p = 0.20) greater regulatory strength than yoked ones (P2;
Holm-adjusted p ≈ 10⁻⁴ at both) and matched the target roughly twice as well
(metabolic error 0.029 versus 0.058 at p = 0.10); in the exploratory experiment
they also carried more transcription factors (11.2 versus 6.2). In a frozen environment, where the signal still
switches but announces nothing because the target never changes, informative and
yoked cues were indistinguishable (7.5 versus 5.4, p = 0.18). The difference
between cue arms therefore arises only when the cue carries information about
something that changes — the property the yoked design is built to isolate.

This reverses EvoSim's earlier conclusion that information content is
irrelevant. The two results are compatible once cue use is taken into account:
in the EvoSim regime where that conclusion was drawn, organisms exploited only
about a tenth of the information available to them (S4 Text), whereas R-Aevol
organisms halved their error when given an informative cue.

### Fine-grained fluctuation shrinks the genome

Fluctuation reduced the number of functional genes at every switching rate, by
10, 10 and 8 genes (26.3, 28.3 and 25.6 in frozen against 16.3, 18.1 and 18.0 in
fluctuating populations; P4; Holm-adjusted p ≤ 3 × 10⁻⁴; Fig 4), and at a
population of 512 (25.3 against 17.5, p < 10⁻⁴). This is the opposite of the
genome enlargement reported for fluctuating selection in other models [6,7] and
of EvoSim's own earlier finding (S4 Text). One plausible mechanism is that in
Aevol precision is bought with genes: the phenotype is a sum of contributions
from individual genes, so a constant target that can be matched closely rewards
additional genes, whereas a target that switches within a lifetime cannot be
matched closely and leaves a coarser, cheaper phenotype. Frozen populations
reached metabolic errors near 0.008, fluctuating ones without a usable cue near
0.06. We did not test this mechanism.

> **Fig 4. Pre-registered test.** Effect of fluctuation on total regulatory
> strength (P1, P3), of information (P2), and on gene number (P4), at each
> switching probability, as differences with 95% CIs. Open symbols: cells
> excluded from the verdict by the manipulation check. Right panel: the same
> contrasts at population 512.

### Counting connections gives the wrong sign

In both models, connection counts and connection strength moved in opposite
directions. In EvoSim, fluctuation raised the count and lowered strength per
possible connection. In R-Aevol, informative-cue populations in a fluctuating
environment had *fewer* raw regulatory links than their frozen counterparts at
every switching rate (0.60, 0.78 and 0.83 times as many), while their total
regulatory strength was 1.8–2.8 times higher and their strength per link roughly
twice as high (Fig 5). A connectivity measure based on counting links would
therefore have concluded that fluctuation with a usable cue *simplifies*
regulatory networks. This comparison was specified in advance as a secondary
analysis and carries no verdict.

> **Fig 5. Counting versus weighing.** Ratio of fluctuating to frozen raw link
> count against ratio of total regulatory strength, per switching rate
> (R-Aevol) and for EvoSim; the dashed diagonal marks agreement between the two
> measures.

### ⟨PENDING⟩ Grain

To address the conflict with Tsuda and Kawata [6], who found that cue-free
fluctuation across generations promotes regulatory complexity, we varied grain
in R-Aevol without changing the model. The environment is a single Markov chain
over ages that continues from one generation to the next, so the switching
probability sets the mean time an environment persists: from half a lifetime at
0.1 to about 250 generations at 0.0002. We predicted that coarse-grained
cue-free fluctuation would produce larger genomes and stronger regulation than
fine-grained fluctuation (G1, G3) and than a frozen environment (G2, G4) —
Tsuda and Kawata's direction — and fixed in advance how each outcome would be
interpreted (Methods). ⟨Results: 60 runs, expected 2026-09-28.⟩

> **Fig 6. ⟨PENDING⟩ Grain.** Gene number and total regulatory strength against
> mean environmental dwell time for cue-free fluctuating populations, with the
> frozen mean for reference.

---

## Discussion

Environmental fluctuation did not select for regulatory complexity as such. In
R-Aevol, whether it strengthened regulation depended on whether the organism
could tell which environment it was in: with an informative cue, fluctuation
produced the strongest regulation we observed; with none, the weakest. A yoked
cue that delivered the same molecules at the same rate did not substitute for
an informative one, so it was the information, not the stimulation, that
mattered. These conclusions held across switching rates and population sizes,
in a test whose predictions and analysis were fixed before the data existed.

This extends a well-established result at the level of the phenotype to the
level of architecture. That reliable cues favour plastic responses and their
absence favours generalism or bet-hedging is classical [1–4], and has recently
been shown for the evolution of plasticity in a gene regulatory network [5].
Our contribution is to show the corresponding consequence for the network that
implements the response — its regulatory strength and its size — and to
separate a cue's information from its presence, which manipulations of cue
reliability alone cannot do.

⟨PENDING: paragraph on grain and Tsuda & Kawata. EvoSim's fluctuation is
coarse relative to a lifetime — development occurs once at birth and each
environmental state lasts on the order of twenty generations — and EvoSim
earlier found genome enlargement under fluctuation, as Tsuda and Kawata do;
R-Aevol's fluctuation is fine-grained and shrinks the genome. The grain
experiment tests whether this is the reconciling variable within one model.⟩

Counting connections gave the opposite answer to weighing them in both models,
for different reasons. In EvoSim, duplication and node addition copy a node's
connections with small added noise, so each copy creates many weak connections
where there were none; these persist when purifying selection is weak, as it is
in the smaller populations of a fluctuating environment, and connection number
rises while strength does not. In R-Aevol, frozen populations carried more but
weaker regulatory links, while informative-cue populations carried fewer links
of roughly twice the strength; why frozen populations accumulate weak links we
did not establish. In both cases
the count tracks the history of mutation more closely than the function of the
network. Since connection counts are the usual measure of regulatory complexity
in evolutionary models [6,10,11], we suggest that claims about connectivity
report a weighted measure alongside or instead.

**Limitations.**
- *One model for the trackable case.* The positive result — that a usable cue
  lets fluctuation strengthen regulation — comes from R-Aevol. EvoSim organisms
  read a cue only when development is reversible and exploit little of it; we
  have not measured regulatory strength in that regime.
- *Replication depth.* Each R-Aevol cell has 12 populations, fewer than the 60–100
  common in comparable work [6,18]. The effects are large (most d > 2) and survive
  correction for multiple testing, but modest effects would be missed.
- *Measures.* R-Aevol statistics describe the best individual in each population
  at each generation, as the platform records them.
- *Environments.* R-Aevol's fluctuation alternates between two fixed targets. It
  does not include the random walk to novel optima that Tsuda and Kawata found
  most favourable to complexity [6].
- *Population size.* In EvoSim, fluctuation and population size are confounded
  across arms; the within-arm evidence argues against population as the
  explanation of the strength decline, but cannot exclude it.
- *Genome size.* P4 confirms the direction of the genome-size effect, not its
  mechanism.
- *A pre-registered confirmation that could not run.* We pre-registered a
  confirmation of the EvoSim count–strength divergence in a regime with more
  hidden nodes, conditional on fluctuating populations surviving there. They did
  not (2 to 5 of 6 survived, against 6 of 6 frozen), so it was not run and the
  EvoSim result remains exploratory.
- *Timestamps.* Pre-registrations were committed to version control before the
  corresponding runs; for the first sweep these timestamps are local and not
  independently verifiable.

**Conclusion.** A changing world does not in itself build stronger regulation.
It does so when the organism has the information to track the change, and when
it has not, fluctuation leaves regulation weaker than a constant world would.
Whether a model says otherwise may depend on whether it counts connections or
weighs them.

---

## Methods

### EvoSim

EvoSim is described in full, following the ODD protocol [19], in S1 Text; we
summarise the parts relevant here. Organisms live on a G × G lattice (G = 128)
with six diffusing resource channels. Only the first receives external inflow;
consuming channel r releases a byproduct into channel r + 1, with declining
yields along the chain. Each organism carries an energy store, divides when it
exceeds a heritable threshold, and dies when energy runs out.

*Genome and development.* A genome is a weight matrix W and bias vector b over
ng nodes, 13 ≤ ng ≤ 32. Nodes 1–12 are read out as the phenotype and node 13
receives the abiotic cue; hidden nodes 14…ng carry no fixed meaning.
Development runs once, at birth: expression is initialised at 0.1 and updated
eight times as e_k ← tanh(b_k + Σ_m W_mk e_m + input_k), where input_k is the
local concentration of resource k for k ≤ 6 and the cue for k = 13. Expression is
then mapped to the twelve traits.

*Mutation.* At each birth, each weight receives an additive Gaussian perturbation
(σ = 0.25) with probability 10⁻³, reduced by the organism's heritable repair
trait by up to 85%. Structural mutations: gene silencing (probability 0.002;
zeroes a node's outgoing weights), gene duplication (0.002; copies a node's
regulatory role onto another with N(0, 0.05²) noise on every copied weight), and
addition or deletion of a hidden node (0.003 each; an added node duplicates an
existing one in the same way). Connection silencing (probability 0.2 per birth)
sets one randomly chosen weight to exactly zero; without it, additive mutation
makes density ratchet towards one. Each hidden node costs 0.0015 and each
non-zero weight 5 × 10⁻⁵ energy units per time step.

*Environment and cue.* An abiotic condition cycles through `cond_states` discrete
states with a period of 2,000 time steps; each state rotates which resource
channel pays best, and 40% of external inflow follows the favoured channel. The
condition is fed to node 13. Fluctuating arms use six states, frozen arms one.
The yoked cue replaces the condition's value with one generated on an
incommensurate period (a multiplier of 5 + φ), so that it cycles through the
same values without predicting the state; the no-cue arm sets the cue's gain to
zero. With a single state the cue is identically zero, so the three frozen arms
are identical runs and were analysed as one arm.

*Measures.* Count density is the number of weights with |w| > 10⁻⁴ divided by
ng²; strength density is Σ|w| / ng² over all ng² weights, zeros included; total
strength is Σ|w|. Each is computed per organism at birth, averaged over the
population at each logged time step, and averaged over the final 20% of the run.
Populations that went extinct were excluded and counted.

*Runs.* 180 runs: six arms × 30 seeds, 150,000 time steps each. Runs are
deterministic given the seed and bit-identical across thread counts.

### R-Aevol

We used the `aevol-7` branch of Aevol (commit 33277374,
https://gitlab.inria.fr/aevol/aevol) compiled with regulation enabled
(`-Dwith-raevol=ON`) and per-generation statistics enabled
(`-Dwithout-stats=OFF`); building on Apple Silicon required correcting a
hard-coded OpenMP library path (S2 Text). Parameters are in Table 1 and follow
the platform's own two-environment example (`examples/dynamic_env`), including
its binding matrix, used unchanged in every run.

*Environments.* Each environment is a target function built from three
Gaussians. Target B is target A reflected about the midpoint of the phenotypic
axis; both have geometric area 0.152468. Fluctuating arms switch between A and B;
frozen arms use two copies of the same target (A for odd seeds, B for even) so
that the switching machinery runs identically but the organism experiences no
change. The environment is a Markov chain over ages that continues across
generations: at each of the 20 ages it moves, with probability p, to a different
environment chosen uniformly among the others.

*Signals and the yoked control.* Two signal proteins, taken unchanged from the
platform's `complex_env` example, are injected at concentration 0.9 during ages
when their environment is active and are otherwise absent. Informative arm: two
environments, (A, s1) and (B, s2). Yoked arm: four environments, (A, s1),
(A, s2), (B, s1), (B, s2). No-cue arm: no signals. In the two-environment arm
every switch changes both target and signal, so both change at rate p per age.
In the four-environment arm a switch changes the target with probability 2/3 and
the signal with probability 2/3, so at switching probability 1.5p both again
change at rate p per age.

*Measures.* R-Aevol writes, each generation, statistics of the best individual.
Genome size is the number of functional genes. Total regulatory strength
S = nb_link × avg_influence, where nb_link is the number of protein–RNA
influences and avg_influence their mean strength (the sum of enhancing and
operating factors divided by their number); S = 0 when there are no links. Each
run is summarised by its mean over generations 25,001–30,000.

*Runs.* Exploratory: six arms × 12 seeds at p = 0.10, population 1024
(72 runs). Pre-registered sweep: fluctuating informative, yoked and no-cue arms
and frozen informative arms at p ∈ {0.05, 0.10, 0.20}, frozen no-cue and frozen
yoked arms at p = 0.10 (frozen no-cue is rate-independent), 12 new seeds each,
population 1024 (168 runs), plus the fluctuating and frozen informative arms at
population 512 (48 runs). Grain experiment: cue-free fluctuating arms at
p ∈ {0.1, 0.01, 0.001, 0.0002} and a frozen arm, 12 new seeds each (60 runs).
All runs reached 30,000 generations.

> **Table 1. R-Aevol parameters.** Population 1024 (32 × 32 grid) or 512;
> fitness-proportionate selection, pressure 750; point mutation, small
> insertion, small deletion, duplication, deletion, translocation and inversion
> rates 5 × 10⁻⁵ per base per generation; maximum indel 6; initial genome 5,000
> bases with one good gene; Hill function n = 4, θ = 0.5; degradation rate 1,
> 10 degradation steps per age; 20 ages, fitness evaluated at all 20; binding
> matrix fixed (the platform example's file, identical in every run); target
> sampling 300 points;
> signal concentration 0.9.

### Statistics

The unit of analysis is a population. Each run contributes one number per
measure, so replicate generations are never treated as independent. Differences
between arms were tested with two-sided permutation tests on the difference of
means (100,000 permutations), with rank-based permutation tests reported
alongside. Multiple predictions were corrected with the Holm procedure. Effect
sizes are Cohen's d.

### Pre-registration and deviations

Predictions P1–P4, their measures, tests and decision rules, and the analysis
code were committed to version control (commit f985ccd) before the
pre-registered sweep was launched; the analysis was run once, on complete data,
with the committed code (verified by checksum). The decision rules are given in
Table 2. The only deviation was operational: the host machine's idle sleep
suspended runs overnight, which lengthened wall time but cannot affect results,
since runs are deterministic given their seed. The EvoSim confirmation (commit
4bc16c3) required fluctuating and frozen populations to survive in a
viability pilot; none of the three candidate regimes qualified and it was not
run (commit 56fe3a0). The grain experiment (commit fa00777) ⟨PENDING: OSF
registration if made; outcome⟩.

> **Table 2. Pre-registered predictions and verdicts.**
>
> | id | prediction | verdict |
> |---|---|---|
> | P1 | fluctuation with an informative cue raises total regulatory strength over frozen | supported (p = 0.10, 0.20) |
> | P2 | informative cue raises it over a yoked cue | supported (p = 0.10, 0.20) |
> | P3 | fluctuation without a cue lowers it below frozen | supported (all rates) |
> | P4 | fluctuation lowers gene number | supported (all rates) |
> | P5 | EvoSim: count density up, strength density down, in a regime above the node floor | not run (no viable regime) |
> | G1–G4 | coarse grain reverses the effect of cue-free fluctuation | ⟨PENDING⟩ |
>
> Rules: supported if the predicted sign holds at every evaluable rate and the
> Holm-adjusted p < 0.05 at two or more; refuted if the opposite sign is
> significant at any rate or no rate is Holm-significant; inconclusive
> otherwise. P1 and P2 are evaluable only at rates where the informative cue
> measurably lowers metabolic error relative to the yoked cue.

### Software, hardware and availability

EvoSim is written in Julia (1.12–1.13) with no dependencies beyond the standard
library for simulation. R-Aevol runs were made on a single 12-core machine.
All code — both models' configurations, launchers, pre-registered analysis
scripts and the full history of the analyses — and the per-run statistics are
available at ⟨repository URL and Zenodo DOI⟩.

---

## Supporting information

- **S1 Text.** Full ODD description of EvoSim.
- **S2 Text.** R-Aevol build, configuration files and the yoked-control derivation.
- **S3 Table.** Per-arm and per-run results for every experiment.
- **S4 Text.** Earlier EvoSim results and the audit trail: claims that were
  superseded, the measurement traps that produced them, and the control that
  caught each.
- **S5 Table.** Prior-art comparison.

---

## References

1. Levins R. Evolution in changing environments. Princeton: Princeton University Press; 1968.
2. Moran NA. The evolutionary maintenance of alternative phenotypes. Am Nat. 1992;139(5):971–989.
3. Kussell E, Leibler S. Phenotypic diversity, population growth, and information in fluctuating environments. Science. 2005;309(5743):2075–2078.
4. Botero CA, Weissing FJ, Wright J, Rubenstein DR. Evolutionary tipping points in the capacity to adapt to environmental change. Proc Natl Acad Sci USA. 2015;112(1):184–189.
5. Jorritsma A, van den Berg P. The evolution of plasticity and evolvability in a simple gene regulatory network. Evolution. 2026;80(1):230–239. doi:10.1093/evolut/qpaf219
6. Tsuda ME, Kawata M. Evolution of gene regulatory networks by fluctuating selection and intrinsic constraints. PLoS Comput Biol. 2010;6(8):e1000873.
7. Bentkowski P, van Oosterhout C, Mock T. A model of genome size evolution for prokaryotes in stable and fluctuating environments. Genome Biol Evol. 2015;7(8):2344–2351.
8. Knibbe C, Coulon A, Mazet O, Fayard JM, Beslon G. A long-term evolutionary pressure on the amount of noncoding DNA. Mol Biol Evol. 2007;24(10):2344–2353. [verify]
9. Luiselli J, Rouzaud-Cornabas J, Lartillot N, Beslon G. Genome streamlining: effect of mutation rate and population size on genome size reduction. Genome Biol Evol. 2024;16(12):evae250.
10. Espinosa-Soto C. On the role of sparseness in the evolution of modularity in gene regulatory networks. PLoS Comput Biol. 2018;14(5):e1006172.
11. Kashtan N, Alon U. Spontaneous evolution of modularity and network motifs. Proc Natl Acad Sci USA. 2005;102(39):13773–13778.
12. Church RM. Systematic effect of random error in the yoked control design. Psychol Bull. 1964;62:122–131.
13. Soyer OS, Bonhoeffer S. Evolution of complexity in signaling pathways. Proc Natl Acad Sci USA. 2006;103(44):16337–16342.
14. Iwasaki WM, Tsuda ME, Kawata M. Genetic and environmental factors affecting cryptic variations in gene regulatory networks. BMC Evol Biol. 2013;13:91.
15. Sanchez-Dehesa Y, Parsons D, Peña JM, Beslon G. Modelling evolution of regulatory networks in artificial bacteria. Math Model Nat Phenom. 2008;3(2):27–66.
16. Beslon G, Sanchez-Dehesa Y, Parsons DP, Peña JM, Knibbe C. Scaling laws in bacterial genomes: a side-effect of selection of mutational robustness? BioSystems. 2010;102(1):32–40. [verify]
17. Aevol-9: a simulation platform to decipher the evolution of genome architecture. bioRxiv. 2025. doi:10.1101/2025.04.10.648095 [verify authors]
18. Lalejini A, Ferguson AJ, Grant NA, Ofria C. Adaptive phenotypic plasticity stabilizes evolution in fluctuating environments. Front Ecol Evol. 2021;9:715381.
19. Grimm V, Railsback SF, Vincenot CE, Berger U, Gallagher C, DeAngelis DL, et al. The ODD protocol for describing agent-based and other simulation models: a second update to improve clarity, replication, and structural realism. J Artif Soc Soc Simul. 2020;23(2):7.
20. Pollack D, Nozoe T, Kussell E. Gene regulation, phenotypic memory, and selection in fluctuating environments. bioRxiv. 2025. doi:10.1101/2025.06.23.661119
21. LaBar T, Adami C. Different evolutionary paths to complexity for small and large populations of digital organisms. PLoS Comput Biol. 2016;12(12):e1005066. [verify]
22. Crombach A, Hogeweg P. Evolution of evolvability in gene regulatory networks. PLoS Comput Biol. 2008;4(7):e1000112.

*Refs 20–22 are for the Discussion paragraphs still to be written (regulation
selected for or against by environment structure; population size and
complexity; evolution of evolvability). The comparative-genomics sentence in
the first Introduction paragraph needs its own citations.*
