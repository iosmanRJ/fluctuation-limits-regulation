# EvoSim

An agent-based evolution simulator. Single-celled organisms with
gene-regulatory-network genomes live in a 2D petri dish of diffusing,
cross-fed resources, under seasonal cycles and mass-extinction events, for
hundreds of thousands of generations.

```sh
# core simulator — stdlib only, no plotting stack, runs on a headless box
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. -t auto scripts/run.jl --gens 200000 --out runs/first

# charting — separate environment, pulls in Makie
julia --project=scripts -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=scripts scripts/dashboard.jl runs/first
julia --project=scripts -t auto scripts/live.jl     # watch one evolve
```

**Two environments on purpose.** The core package depends on nothing but the
standard library, so a headless machine can run simulations without OpenGL or
any plotting stack. `scripts/Project.toml` adds Makie for charting and depends
on the core by path.

---

## Ticks are not generations

**This is the single most important thing to know about reading results.**

A *tick* is one simulation step. A *generation* is one round of descent. One
generation takes **10–30 ticks**, and the ratio is **not constant** — it grows
with population, because a crowded dish means less energy per organism and
slower division. Measured: ~15 ticks/generation at grid 128, ~29 at grid 256.

So a million-tick run is tens of thousands of generations, not a million, and
you cannot convert between them with a fixed factor.

Every run logs `mean_depth` — the mean lineage depth of living organisms —
which is the real generation count. Never infer it from the tick number. The
runner prints the measured ticks-per-generation ratio when it finishes.

## The model

### Genome

A gene-regulatory network: a weight matrix plus per-node bias. At birth the
network is run to an approximate fixed point, and the resulting expression
vector is read off as a **12-trait phenotype** which is then cached for life.

**Genome size is variable.** Nodes `1:NPHEN` are a fixed interface — the first
`NRES` receive environmental input and the first `NPHEN` are read off as the
phenotype — so they cannot be added or removed without breaking the mapping.
Nodes `NPHEN+1 : ng` are **hidden units**, gained by duplication and lost by
deletion, so an organism's `ng` evolves over `NG_MIN`(12) to `NG_MAX`(32).
`maint_gene` charges for each hidden node, since more genes cost more to
replicate and express.

This exists to make the open-endedness question answerable. With a fixed node
count it was not: connectivity saturated at ~392 of 400 possible weights in
every run, and "complexity plateaued" could not be told apart from "complexity
hit the ceiling".

Development happens **once, at birth** — not every tick. Re-expressing every
tick measured 15× slower (44.3 vs 2.9 ms/tick, on an earlier version of the
model) for no modelling benefit, and "develop, then live as an adult" is the
more faithful model anyway.

Local resource concentrations feed the first six network nodes, so the same
genome develops differently in different conditions (phenotypic plasticity).
The `plasticity` parameter scales that input; `0.0` makes allocation purely
genetic, which is useful for asking whether plastic responses are masking
heritable variation (here, they are not).

Mutations are point perturbations of weights, plus gene **silencing** and gene
**duplication** — the latter is how regulatory networks gain new functions.

### Phenotype

| Slots | Trait | Notes |
|---|---|---|
| 1–6 | uptake affinity per channel | an **allocation budget**, sums to 1 |
| 7 | body size | scales uptake *and* upkeep (charged on size²) |
| 8 | motility | movement probability, costed |
| 9 | division threshold | energy needed to divide |
| 10 | repair fidelity | suppresses mutation rate, costed |
| 11 | stress tolerance | shock survival **and** waste-toxicity protection (charged on stress²) |
| 12 | leakage | fraction of intake released as byproduct, above a floor |

There is deliberately **no** "which byproduct do I emit" trait — see
[Cross-feeding](#cross-feeding-how-it-is-modelled-and-why).

Dietary breadth can additionally be charged a concave proteome cost via
`maint_channel` / `channel_cost_exp` (off by default — it costs more than it
buys; see [Still open](#still-open-coexisting-ecotypes)).

**Every trait that buys something also costs something.** This is not
decoration — without a cost, evolution pins a trait to 1.0 and the model stops
being interesting. Size is charged quadratically because its benefit is
multiplicative and frequency-dependent (being bigger than your neighbour pays),
which under a linear cost is an arms race with only a boundary solution.

Uptake affinity is an *allocation* rather than six independent dials. With
independent dials, breadth is merely expensive and the optimum is always
"specialise on the highest-yield channel", which kills cross-feeding. With a
budget, breadth costs opportunity, so an under-exploited channel is always
worth something to somebody.

### Ecology

Six resource channels diffuse across the dish with no-flux edges. **Only
channel 1 has external inflow.** Channels 2–6 exist purely as metabolic waste,
so every other niche must be *created* by some organism before anything can
live in it.

Byproducts run strictly downstream: eating channel `r` leaks part of the intake
into channel `r+1`, and channel 6's byproduct leaves the system. Yields decline
along the ladder, so each rung is a living but a poorer one than the rung above.
That ladder is fixed chemistry, not an evolvable trait — see
[Cross-feeding](#cross-feeding-how-it-is-modelled-and-why) for why that matters
more than it sounds.

Three pressures keep selection from settling:

- **Cross-feeding** — the fitness landscape is built by the population itself.
- **Seasons and shocks** — sinusoidal inflow, plus `cull` (indiscriminate
  bottleneck), `heat` (selects on stress tolerance) and `famine` (inflow
  collapses for a while) events.
- **Mutation load** — most mutations are deleterious, some beneficial, so load
  can fall as well as rise; Muller's ratchet is a dynamic, not a countdown.

### Dormancy

Below `dorm_threshold` energy an organism goes dormant: upkeep drops to
`dorm_upkeep` of normal, **ageing pauses**, it cannot divide, and it keeps
taking up resources at `dorm_uptake` of normal so it can climb back out. This is
sporulation.

It is not decoration — it is what makes long runs possible at all. **Every
extinction measured was caused by famine**: 8 of 8 across a 24-run batch, each
within 24–211 ticks of one, with no cull or heat event killing anything. At
~0.53 mean energy against ~0.3 upkeep per tick an organism holds under two
ticks of reserve against a 400-tick famine, so surviving one was impossible by
construction rather than by ecology.

Measured under a deliberately brutal regime — famines every ~8k ticks, about 3×
normal, over 60k ticks, 8 seeds:

| `dorm_upkeep` | extinctions |
|---|---|
| none | **7/8** |
| 0.200 | **1/8** |
| 0.050 | **1/8** |
| 0.020 | 2/8 |
| 0.005 | 3/8 |

**Shallow dormancy beats deep**, which is the opposite of the intuition it was
built on. Going from 20% upkeep to 0.5% makes things *worse*: dormant cells take
up at only 35% rate, so digging too deep means surviving the famine but never
climbing back out, and the population lingers dormant instead of recovering.
**Recovery is the binding constraint, not survival.** The default sits at 0.05,
in the middle of the flat region.

### Clades

A full per-organism phylogeny is not representable over a million ticks. A
clade is instead a **phenotypic neighbourhood**: a child that has drifted
further than `clade_radius2` from its clade's founder phenotype founds a new
one. Clades with no living members or descendants are garbage-collected
periodically — without that the table grew ~40 entries per tick, several GB
over a long run.

## Architecture

| File | Contents |
|---|---|
| [src/genome.jl](src/genome.jl) | GRN development, mutation, inheritance |
| [src/params.jl](src/params.jl) | every knob, mutable so the live view can retune |
| [src/population.jl](src/population.jl) | struct-of-arrays population, slot allocator, clade table |
| [src/environment.jl](src/environment.jl) | resource fields, diffusion, replenishment |
| [src/sim.jl](src/sim.jl) | the per-tick loop |
| [src/stats.jl](src/stats.jl) | stats log and snapshot IO |
| [src/run.jl](src/run.jl) | headless driver |
| [scripts/run.jl](scripts/run.jl) | CLI |
| [scripts/dashboard.jl](scripts/dashboard.jl) | post-hoc charts |
| [scripts/live.jl](scripts/live.jl) | live interactive window |

Performance decisions that are load-bearing, and should not be undone casually:

- **Struct-of-arrays.** The hot loop reads a few fields of every organism; an
  array-of-structs layout would drag whole 1.6 KB genomes through cache to read
  one float. Genomes live in their own array and are touched only at birth.
- **O(1) slot allocation** via a free-list and a watermark. Scanning for a free
  slot is O(capacity) per birth and dominates everything else.
- **Parameters hoisted into locals** at the top of the loop. `Params` is
  mutable, so the compiler must reload every field access unless you pull them
  out first.
- **The hot loop allocates nothing.** There is a test asserting this; if it
  starts failing, something in the loop began allocating and throughput will
  have fallen by a large factor.

## Output format

A run directory holds:

- `stats.bin` — flat `Float32`, `NSTAT` values per record, no framing. Read it
  with `read_stats(dir)`, which returns an `NSTAT × nticks` matrix. Column
  names are in `STAT_NAMES`; index by name with `stat_index("pop")`.
- `snap_*.bin` — full resource grid plus a sample of organisms, on a
  **log-spaced** schedule. Evolution is roughly log-uniform in interest: the
  first thousand ticks deserve as many samples as the last million, and a
  linear schedule spends all its resolution on the boring end.
- `clades.tsv`, `shocks.tsv`, `meta.txt` — phylogeny, event log, and the full
  parameter set the run used.

## Performance

Re-measured on the current model. Population is resource-limited, so it is
reported alongside throughput — ticks/s depends on how many organisms are alive
far more than on grid size alone.

| grid | capacity | population | organism arrays | ticks/s |
|---|---|---|---|---|
| 128 | 65,536 | ~14,500 | 282 MB | **~900** |
| 256 | 262,144 | ~62,000 | 1.13 GB | **~240** |
| 512 | 1,048,576 | ~272,000 | 4.5 GB | **~55** |

Memory scales as `NG_MAX²` because every organism's weight block is dimensioned
for the largest genome it could have, not the one it has. Raising `NG_MAX` from
20 to 32 for variable genomes multiplied the arrays by 2.6×. Throughput figures
predate that change and are **due a re-measurement on an idle machine** — a
benchmark taken while other runs were competing for cores showed every phase
slowing by ~4×, including `diffuse!`, which never touches a genome and so
cannot have been affected. Do not benchmark this under load.

Grids 128 and 256 measured on an Apple M5 (10 threads); grid 512 is the median
of 57 samples from a live run on an Apple M4 Pro (12 threads), where the range
is 44–132 ticks/s — throughput swings by ~3× across the seasonal population
cycle, so a single instantaneous reading means little.

Reaching a million *generations* needs well over 15M ticks at grid 512 — use
the measured `mean_depth`, not an assumed ratio. Run long jobs headless with
`scripts/run.jl`: the core has no plotting dependency, so it runs fine on a
machine with no display, and `stats.bin` flushes on a 10s timer so a long run
can be charted or stopped at any point without losing what it has done.

For scale, the same loop in pure Python measured **440 ns per organism-tick
versus 8.3 ns** in Julia — and that comparison excluded the diffusion step,
which Julia also had to do.

### Threading

Run with `-t auto` (or `-t 12`). Two phases are parallel: resource diffusion
across channels, and the expensive part of the birth phase.

Current split, grid 256 on 10 threads:

| phase | ms/tick | share | |
|---|---|---|---|
| `organism_step!` | 2.80 | **67%** | serial — the bottleneck |
| `settle_births!` | 0.95 | 23% | parallel |
| `diffuse!` | 0.26 | 6% | parallel across channels |
| `replenish!` | 0.14 | 3% | serial, negligible |

Parallelising the birth phase was worth **2.1×** when it was measured against
the serial version (20.76 → 9.87 ms/tick at grid 512, 48 → 101 ticks/s). It
used to be 72% of runtime on its own; it is now 23%, which is why
`organism_step!` dominates the table above.

Births dominate because each one develops a fresh regulatory network — tens of
thousands of developments per tick at this scale. The parallel phase has no
shared writes: each child touches only its own genome and phenotype columns,
and the environment is read-only while it runs.

**Parallelism does not cost reproducibility.** Each child's RNG seed is drawn
serially in a fixed order, so its random stream depends on its own seed and
never on how the work was divided. Verified identical population, phenotype sum
and clade count at 1, 4 and 10 threads. Per-thread RNGs would have been simpler
and would have made results thread-count dependent.

**The remaining opportunity is `organism_step!`**, now ~48% of runtime and
still serial. It cannot be threaded naively: organisms in the same cell deplete
a shared resource pool by sequential read-modify-write, so concurrent access
would double-spend. Doing it properly means partitioning organisms into grid
row-bands so no two threads touch the same cell, which needs a counting sort by
band each tick (~0.5 ms) and deferred movement, since an organism crossing a
band boundary would otherwise race on the occupancy grid. Estimated gain is a
further ~1.5×; it has not been done because a correct 2.1× beats a racy 3×.

## Development

```sh
git config core.hooksPath .githooks     # run tests before every commit
julia --project=. -t auto test/runtests.jl
```

- Work on a branch; `main` is the integration branch.
- CI runs the suite on Julia 1.10 and 1.12, Linux and macOS.
- The suite is fast (a few seconds) deliberately — it is meant to be run
  constantly, not nightly.

### What the tests actually protect

They are not smoke tests. Each one guards something that has already broken or
would break silently:

- **Diffusion conserves mass** under no-flux boundaries. A boundary bug leaks
  resources and quietly changes every ecological result.
- **The occupancy grid matches reality** across births, deaths, moves and mass
  extinctions. It is maintained incrementally; if it drifts, crowding costs
  decouple from the actual population and nothing downstream means anything.
- **Clade pruning preserves the exact partition** of living organisms. Pruning
  renumbers every id in the table and the population at once.
- **The uptake budget sums to 1** for every organism, including under
  plasticity.
- **Raising a trait's cost lowers its evolved value** — the property that keeps
  traits off the boundary.
- **Same seed reproduces exactly**, and identically at any thread count.
  Without this no result is ever checkable.
- **Allocation does not scale with population.** The test deliberately asserts
  the scaling property, not a fixed byte budget: a budget passes while
  per-organism garbage creeps back in. Two per-birth allocations were found
  this way — `Random.seed!(::Xoshiro, ::UInt64)` costs 368 bytes a call, and a
  `view(phen, :, child)` escaping into a call allocates.

### Reproducibility, precisely

Same machine and same Julia version: bit-identical, at any thread count.

*Across* Julia versions the discrete trajectory holds for a long time but is
not guaranteed forever — LLVM reassociates the SIMD dot product in `develop!`
differently between versions, so some phenotypes differ by ~1 ULP. Measured on
1.12.7 vs 1.13.0 after 500 ticks: population and lineage depth identical, float
sums differing in the 8th significant figure. Eventually a 1-ULP difference
flips a comparison and the runs diverge. Pin the Julia version when a result
needs to be exactly reproducible.

### A warning from experience

Three separate rounds of this model looked fine and were degenerate. The
failure mode is always the same shape: some quantity has a benefit with no
matching cost, or a constraint that was meant to be ecological is actually an
implementation limit. Symptoms to watch for on the dashboard:

- Traits pinned at exactly 0.0 or 1.0 → that trait has no interior optimum.
- Population pinned at exactly `capacity` → the array bound has become the
  selective force, births are being dropped, and selection has gone slack.
- Trait spread collapsing to zero → the population is monomorphic and nothing
  further will happen no matter how long you run it.
- An ecological *mechanism* quietly fading out (cross-feeding did this) → look
  for a degree of freedom the organisms can exploit to short-circuit it. Making
  the excretion destination evolvable let them excrete waste back into the
  channel they had just eaten, which silently disabled the whole waste economy.
  If something is chemistry rather than strategy, do not make it a trait.

### Two ways to measure the wrong thing

Both of these produced confident, completely wrong conclusions here. They are
general, not specific to cross-feeding.

**A mean cannot tell you whether a population has split.** Mean allocation of
ch1 0.70 / ch2 0.19 is *identical* whether every organism is a 0.70/0.19
generalist or 70% are pure-ch1 and 25% are pure-ch2. Standard deviation is not
much better. For anything that might be bimodal, plot a histogram or count how
many organisms fall in each category — the trait-spread panel on the dashboard
is a hint, not an answer.

**Judge a resource by its renewal flux, never by its standing pool.** Computing
the value of each channel from how much of it was lying around said channels 3
and 4 were worth ~1.8× channel 1, and an invasion calculation built on that said
every specialist should sweep. Both were wrong: a waste channel's pool is large
precisely *because nothing is eating it*, while its renewal flux is several
leakage steps down the ladder. A specialist on it eats the accumulated pool once
and then starves.

Charts are the debugger here. Run a few tens of thousands of ticks and look at
them before committing to a long run. The full worked example of both mistakes
is in [Still open: coexisting ecotypes](#still-open-coexisting-ecotypes).

## Cross-feeding: how it is modelled, and why

The waste channels exist so organisms can build niches for each other. Getting
that to persist rather than decay took three specific things, each of which the
microbial consumer-resource literature is explicit about.

### 1. A byproduct never returns to its own channel

Channel `r`'s byproduct is channel `r+1`, always; the last channel's byproduct
leaves the system. **Where a byproduct goes is chemistry, not a trait.**

This was the bug that broke everything. When the excretion destination *was* an
evolvable trait, organisms found a futile cycle — they evolved to excrete
**36% of their waste straight back into the channel they had just eaten**
(against 17% for a neutral profile). Eat channel 1, excrete a third of it back
as channel 1, re-eat it. The waste channels were never fed, so cross-feeding
had nothing to run on and allocation to it decayed to zero.

The consumer-resource models exclude exactly this: the metabolic matrix `D` is
column-stochastic **with the diagonal excluded**, so a resource cannot be
reconverted into itself, on thermodynamic grounds. A strictly downstream ladder
enforces it and additionally rules out longer cycles.

### 2. Leakage is a constant, not a trait

`leak_min` is 0.45, and it is a floor rather than a target. Left free to evolve,
leakage went to ~0.003 — uptake saturates at whatever is locally available, so
the rate/yield bonus cannot pay for the yield it costs, and the waste channels
starved. Leakage is an environmental constant in this literature (Marsland et
al. run `l = 0.8`) because biochemistry sets it. The trait can still push
leakage *above* the floor, which keeps the rate/yield trade-off alive.

### 3. Waste is mildly toxic, and tolerance is quadratically costly

Accumulated byproducts inhibit growth locally, so eating a waste channel is a
*service* to whoever produced it, not mere scavenging. That is what converts
one-way cross-feeding into a mutualism — the Black Queen mechanism, where
partitioning a detoxifying function stabilises coexistence even for organisms
competing for a single resource.

Tolerance is charged on `stress²`. Under a linear cost this trait has no
interior optimum: it pins at 1.0. **`maint_toxic` is not monotonic** — around
0.015 the dish can die during the ramp-up, before tolerance has evolved. The
0.008 default was checked across four seeds.

### Measured outcome

| | waste-channel uptake |
|---|---|
| before | 0.180 → **0.076** over 80k ticks, still falling |
| after | **0.26–0.30**, flat across 60k ticks and four seeds |

Things that were tried and did **not** work, so nobody repeats them:

- *Starving the dish.* Worse, not better (0.090 → 0.019 as inflow went
  0.25 → 0.08): a smaller population makes less waste, thinning the niches.
- *Raising waste yields.* Slows the decay without stopping it.
- *A convex allocation trade-off* (`uptake_exponent` > 1). Adaptive-dynamics
  theory says evolutionary branching needs a trade-off with positive curvature,
  and the exponent does drive specialisation up (0.987 → 0.999) — but every
  organism specialises on the *same* channel, so there is no branching, and
  small dishes starve. Curvature was not the binding constraint; the futile
  cycle was.

### Still open: coexisting ecotypes

Cross-feeding is stable, but the population remains **one type**. A histogram of
channel-1 allocation shows a single sharp peak with **99.9% of organisms topping
out on channel 1** — not two coexisting specialists. Distinct ecotypes do not
reliably form.

This was investigated properly and the diagnosis is solid even though the fix
is not. Recorded here because three plausible explanations turned out to be
wrong, and re-testing them would waste days.

**Measure branching with a histogram, never with a mean.** A mean allocation of
ch1 0.70 / ch2 0.19 is *identical* whether every organism is a 0.70/0.19
generalist or 70% are pure-ch1 and 25% pure-ch2. `std` is also a poor guide.
Count how many organisms have each channel as their top channel.

**Falsified:**

- *Canalisation.* The regulatory network is not stuck — re-develop the same
  evolved genomes in an environment with channel 4 abundant and channel 1 absent
  and 96–100% of them produce channel-4-dominant phenotypes. The machinery works.
- *Plasticity shielding genotypes from selection.* Setting `plasticity = 0`
  (allocation purely genetic) changes nothing.
- *Mutation supply.* A 10× mutation rate changes nothing.

**A measurement trap worth knowing about.** Marginal value computed from
*standing stock* says channels 3 and 4 are worth ~1.8× channel 1, and an
invasion calculation on that basis says every pure specialist beats the
resident by up to +0.21 energy/tick. **Both are wrong.** A waste channel's stock
is high precisely *because nothing eats it*, while its renewal flux is several
leakage steps down the ladder. A channel-4 specialist would eat the accumulated
pool once and then starve. Judge a strategy by the resource's **renewal flux**,
never by its standing pool.

Working the flux cascade through, the energy carried by each channel goes
roughly `1 : 0.27 : 0.09 : 0.03`, and the observed allocation (0.78 / 0.17 /
0.05 / ~0) *is* that optimum. The population is right; there is simply no
advantage to specialising, because uptake is **linear** in allocation on a
budget that sums to 1. At the equilibrium where resources equalise, every
allocation earns the same. Invasion fitness is flat, its second derivative is
zero, and the singular strategy is degenerate — neutral rather than disruptive.
Adaptive dynamics requires a branching point to be convergence stable *and*
evolutionarily unstable, i.e. a fitness **minimum**.

**What partly works — measured properly with a controlled A/B.** Three
mechanisms together shift the outcome: a concave proteome cost for dietary
breadth (`maint_channel`), spatial structure (low `diff_rate` plus
`diff_falloff` so byproducts diffuse slower than the substrate), and enough
`inflow` to afford the first two.

The diagnostic that explains why is the **spatial correlation between channel 1
and total waste**. Niche separation needs them *anti*-correlated. By default
they are positively correlated, because waste is produced where organisms are
and organisms are where channel 1 is — so every organism sees the same mixture
and there is nothing to partition:

| config | non-ch1 specialists | ρ(ch1, waste) |
|---|---|---|
| baseline | 0.1% | **+0.04** |
| proteome cost only | 5.4% | −0.05 |
| all three | 19.8% | **−0.13** |

Paired A/B, same 12 seeds per arm, grid 128, 150k ticks:

| | control | combo |
|---|---|---|
| surviving runs | 8/12 | 8/12 |
| **extinctions** | **4/12** | **4/12** |
| non-ch1 specialists, median (range) | **0.0%** (0.0–0.1) | **1.7%** (0.3–6.4) |
| runs ≥5% non-ch1 | **0/8** | **3/8** |
| runs ≥10% non-ch1 | 0/8 | 0/8 |
| **clade diversity (Shannon)** | **0.03–1.9** | **4.9–8.2** |
| generations reached | ~6,500 | ~11,600–15,700 |
| genome connectivity (of 400) | 380 | 392 |

Read that carefully, because the headline is not the one being aimed at.

**Dietary branching is real but small.** The control is not merely low, it is
flat zero — eight runs, maximum 0.1%. The combo reaches a median of 1.7% and
gets three of eight runs past 5%. So the mechanism works and the effect is
unambiguous, but it never produced a stable co-dominant ecotype: nothing
cleared 10%.

**It costs nothing in extinction risk** — 4/12 in both arms, identical.

**The large effect is on clade diversity, not diet.** Control runs collapse to
one lineage (Shannon 0.03–1.9); combo runs sustain many (4.9–8.2). The
mechanism preserves lineage diversity far more strongly than it partitions
diets. That was not what it was built for, and it is the more interesting
result: what these three changes buy is a population that does not collapse to
a single clone, even while all those clades keep eating much the same thing.

None of this is the default. The combo needs inflow raised 2.4x, which changes
the energy economy of everything else, and the payoff is a minority ecotype
rather than a split population.

### Genome complexity does not grow — measured

Connectivity saturating at ~392 of 400 weights was originally read here as
evidence that complexity had stopped growing. That reading was not available:
`NG` was a hard cap, so saturation was forced. With genome size now variable the
question is answerable, and the answer is negative.

36 replicates, five arms, 200k ticks (~10,000 generations) at grid 128:

| arm | gene cost | started at | peak | **final** | final connectivity |
|---|---|---|---|---|---|
| cost0 | **0** | 16 | 20.8 | **15.1** | 223 |
| cost15 | 0.0015 | 16 | 17.0 | **13.0** | 161 |
| cost60 | 0.006 | 16 | 14.2 | **12.0** | 141 |
| hi28 | 0.0015 | **28** | 28.2 | **14.1** | 189 |
| lo12 | 0.0015 | **12** | 13.5 | **12.1** | 142 |

**Genome size decays in every arm, monotonically, and connectivity decays with
it** (416 → 218 weights in a typical run). Two results carry the argument:

- **It decays with the gene cost set to zero.** `cost0` still falls from a peak
  of 20.8 to 15.1, so this is not the explicit `maint_gene` charge. Adding a
  node duplicates an existing one and perturbs a network that already works, so
  additions are deleterious on average and selection removes them. Complexity
  is costly here even when nothing charges for it.
- **The floor is absorbing.** Starting at `NG_MIN`, genome size never rises —
  12.1 → 12.0 across 200k ticks. If hidden units conferred any benefit at all, a
  population sitting at the minimum would acquire them. It does not.

Starting high converges downward (28 → 14.1) while starting low stays put, so
this is not a common interior attractor. It is a ratchet toward the minimum.

The honest reading: **this environment does not reward regulatory complexity.**
Mapping six environmental inputs to twelve traits is a task simple enough that
the fixed interface layer suffices, and every hidden unit is overhead. Open-
ended complexity growth would need a task that actually demands it — conditional
responses with memory, organism-organism interaction, or an environment whose
structure keeps changing faster than a minimal network can track. That is the
next thing worth building, and it is a much larger change than adding nodes.

**Where to go next.** The same work finds cross-feeding lineages emerge only
under strong spatial structure, self-organising into interleaved patterns, and
that in well-mixed conditions autonomous mutants take over. The fields here are
close to uniform (`avail@organisms ≈ avail@global`), which is the real blocker:
every organism sees the same resource mixture, so there is one optimal strategy
and nothing to partition. Lowering `diff_rate` raises field heterogeneity
(CV 0.51 → 0.89) and lifts branching to ~2% but no further on its own, and
making motility expensive enough to matter starves the dish. Per-channel
diffusion (`diff_falloff`) is now implemented and alone does nothing —
heterogeneity rises (waste CV 0.125 → 0.302) with branching still at 0%. It
only pays combined with the proteome cost, as the A/B above shows.

What has *not* been tried: making organisms' offspring disperse less (birth
placement is currently a neighbouring cell, which is already local), or giving
the dish a spatial gradient in inflow so that different regions genuinely
favour different channels rather than relying on the population to generate
the heterogeneity itself.

### References

- Marsland et al., *Available energy fluxes drive a transition in the diversity,
  stability, and functional structure of microbial communities* (MiCRM) —
  [bioRxiv](https://www.biorxiv.org/content/10.1101/335893v2.full)
- *Stability and feasibility of the Microbial Consumer-Resource Model* —
  [arXiv:2607.19982](https://arxiv.org/abs/2607.19982)
- Doebeli, *A model for the evolutionary dynamics of cross-feeding polymorphisms
  in microorganisms*, Population Ecology (2002) —
  [Springer](https://link.springer.com/content/pdf/10.1007/s101440200008.pdf)
- van Dijk et al., *Contingent evolution of alternative metabolic network
  topologies determines whether cross-feeding evolves*, Communications Biology —
  [PMC7391776](https://pmc.ncbi.nlm.nih.gov/articles/PMC7391776/)
- Pakzad & DeAngelis, *Traveling wave formation enables strain coexistence in a
  spatial model of bacterial cross-feeding*, Mathematical Biosciences (2026) —
  [PubMed](https://pubmed.ncbi.nlm.nih.gov/41667035/)
- Morris et al., *Coexistence of evolving bacteria stabilized by a shared Black
  Queen function*, Evolution —
  [Wiley](https://onlinelibrary.wiley.com/doi/full/10.1111/evo.12485)
- *Stochastic trade-offs and the emergence of diversification in E. coli
  evolution experiments* — [arXiv:2307.11033](https://arxiv.org/abs/2307.11033)
- Geritz et al., *Evolutionarily singular strategies and the adaptive growth and
  branching of the evolutionary tree* —
  [Springer](https://link.springer.com/article/10.1023/A:1006554906681)
- *Continuously stable strategies as evolutionary branching points* —
  [arXiv:1005.3862](https://arxiv.org/pdf/1005.3862)
- *Shaping of microbial phenotypes by trade-offs* (proteome allocation and the
  specialist/generalist trade-off) —
  [Nature Communications](https://www.nature.com/articles/s41467-024-48591-9)
