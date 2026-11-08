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

#### Dormancy does not fix long-run extinction

The table above is a 60k-tick stress test. Over **200k ticks at the normal
famine rate**, dormancy makes no difference to extinction at all:

| arm | extinct | population | generations | clade Shannon | dormant |
|---|---|---|---|---|---|
| no dormancy | **4/8** | 16,594 | 8,629 | 1.18 | 7.6% |
| dormancy | **4/8** | 30,877 | 6,261 | **4.49** | 35.2% |

Both tests saw a similar *number* of famines (~2.5 vs ~2.7), so the difference
is not famine exposure — it is **run length**. Dormancy buys survival over 60k
ticks and has stopped mattering by 200k. Whatever kills long runs is therefore
not simply an inability to sit out a famine, and that cause is not yet
identified.

What dormancy *does* do is substantial and was not the goal:

- **population ×1.9**, because organisms that would have starved persist
- **clade diversity ×3.8** (Shannon 1.18 → 4.49) — dormant cells act as a
  **seed bank**, so lineages that would have been lost sit out the bad period
  and return. This is a well-established ecological principle and it fell out
  of the model unprompted.
- **generations −28%**, because dormant cells do not divide

Beware of generalising a short stress test. The 7/8 → 1/8 result is real for the
regime it was measured in and does not transfer to longer runs.

#### What did reduce long-run extinction: variable genomes

Not dormancy — **variable genome size**, and by accident rather than design.
Re-running the same 200k-tick configuration after that change:

| code | extinctions at 200k ticks |
|---|---|
| fixed `NG = 20` | 4/8 = **50%** |
| variable genome, seeds 1–8 | 0/8 |
| variable genome, seeds 9–20 | 3/12 |
| **variable genome, combined** | **3/20 = 15%** |

The 0/8 was a lucky draw and is shown here only because reporting it alone
would have overstated the effect — 15% across 20 seeds is the number to use.
Extinction is reduced, not eliminated, and 15% per 200k ticks is still far from
what a million-generation run would need.

**Two observations about why, one falsified and one open:**

- *Falsified:* genome size does **not** shrink during famine. The obvious story
  — shrink the genome to cut upkeep exactly when food is scarce — is wrong.
  Measured across 10 runs, mean `ng` during famine was lower in 5 and higher in
  5. No systematic response.
- *Open:* **larger genomes correlate with dying.** Every run that went extinct
  had mean `ng` above 14.7 (14.9, 15.2, 15.3); every survivor sat below 14.6,
  with most between 12.4 and 13.2. Clean separation, no overlap across 12 runs.
  A permanent `maint_gene` plus development tax leaving less margin when inflow
  collapses would explain it. But those sizes arose by drift, so a run could
  drift large *and* be fragile for a shared upstream reason — correlation only
  until an arm with `ng` pinned at 12 versus 20 settles it.

#### Tested: genome size does not cause fragility

`ng` pinned at 12 versus 20, node add/delete switched off so size never moves,
same seeds, only genome size differing.

**One decisive result:** a 20-node genome costs **40% of carrying capacity**
(median surviving population 32,645 vs 17,846, and the ng=12 populations sit
within ±2% of each other across seeds). Extra genes are a real metabolic tax.

**But no fragility effect.** Extinction was 1/12 versus 2/12 — Fisher exact
p = 1.0, and separating 8% from 17% would need ~250 seeds per arm. Re-run with
a continuous measure instead — minimum population reached as a fraction of that
run's own median, so it is not merely restating the carrying-capacity gap:

| arm | min-pop ratio, median | range | extinct |
|---|---|---|---|
| ng = 12 | 0.0100 | 0.0009 – 0.0520 | 0/6 |
| ng = 20 | 0.0311 | 0.0004 – 0.0377 | 1/6 |

The distributions overlap almost entirely and the ng=20 median is nominally
*higher*, i.e. less fragile — the opposite direction. **There is no effect.**
So "every extinct run had `ng` > 14.7" was coincidence in a small sample, and
the population-level-force-against-complexity idea it suggested has no support.

#### What the bottleneck numbers actually show

| arm | minimum population reached |
|---|---|
| ng = 12 | 21, 91, 159, 303, 780, 1188 |
| ng = 20 | 0, 5, 179, 466, 478, 497 |

**Every run crashes from ~20,000 to double or single digits.** One ng=12 run
fell to **21 organisms and recovered**; an ng=20 run fell to **5** and
recovered; another reached 0 and did not. Survival at these bottlenecks looks
close to a coin flip, and is independent of genome size.

That reframes extinction here. The model is not fragile because of any
particular trait — it is fragile because the seasonal cycle and shocks
routinely drive the population to within tens of organisms of zero, at which
point demographic stochasticity decides. **Bottleneck depth, not genome size,
is what a long run would have to address.** Pinning `ng` small buys carrying
capacity, not safety.

#### What drives the bottlenecks: season sets *when*, famine sets *how deep*

Every deep trough across 6 seeds × 200k ticks, recorded with the season phase,
whether a famine was active, and time since the last shock.

**Season is the clock.** Inflow bottoms at phase 0.75, and **16 of 19 troughs
fall in phase 0.6–1.1** — the half-cycle around and after that minimum, against
~10 expected if phase were irrelevant. The median trough sits at phase 0.91,
about 16% of a cycle after the low, which is the lag of a population tracking a
declining resource.

**Famine sets the depth.** Of the six deepest troughs (<100 organisms), their
ticks-since-last-shock were 428, 436, 429 and 413 — against a `famine_len` of
400. The bottoms land as a famine expires on top of the seasonal minimum.

#### But famine DURATION is not the lever — total severity is

That timing invited an obvious mechanism: a 400-tick famine spans a whole
seasonal descent, so shortening it should stop the two stacking. Shortening
`famine_len` 400 → 150 did give a **23× deeper margin and 0/6 extinctions**.

**That result was confounded and the hypothesis is false.** Integrated severity
is `famine_len × (1 − famine_factor)`, and the arms were not matched: 400 × 0.88
= **352** tick-equivalents against 150 × 0.95 = **142**. The short arm simply
carried 40% of the total food loss — and that cannot be fixed by deepening,
since at 150 ticks the maximum possible deficit is 150, still under 352.

Matching the deficit by making the *long* famine shallower instead:

| arm | deficit | median min-pop ratio | extinct |
|---|---|---|---|
| long + shallow, 400 × 0.65 | 140 | **0.0600** | 0/8 |
| short + deep, 150 × 0.05 | 142 | **0.0469** | 0/8 |

**Indistinguishable.** Paired by seed the short arm was better in only 3 of 8,
and long+shallow is nominally the *safer* of the two — the opposite of the
prediction. Duration does not matter; only the product does.

The deepest troughs landing 413–436 ticks after a shock was a *consequence* of
famine duration, not evidence that duration caused the depth. Suggestive
timing, wrong mechanism.

#### Softening famine removes extinction *and* raises diversity

`famine_factor` therefore defaults to **0.65**, not 0.12, cutting the integrated
deficit from 352 to 140. Validated over 6 seeds × 200k ticks against the old
default:

| metric | harsh 0.12 | soft 0.65 | |
|---|---|---|---|
| extinctions | 1/6 | **0/6** | |
| minimum population | 31 | **1,222** | **39× safer** |
| cross-feeding | 0.277 | 0.270 | unchanged |
| **clade diversity (Shannon)** | 1.15 | **3.48** | **3× higher** |
| trait spread | 0.0235 | 0.0261 | slightly higher |
| generations reached | 5,965 | 5,632 | unchanged |

The expected cost did not appear. Less famine *should* mean weaker selection and
a duller model — instead **diversity tripled**. Harsh famines bottleneck the
population to 6–31 organisms, which purges nearly every lineage; soft ones
bottom out at 720–4,885, so far more clades survive. **The harsh setting was not
generating diversity through selection, it was destroying it through
bottlenecks.**

This is the first configuration under which a million-generation run is
plausible rather than doomed.

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
| [src/genome.jl](../src/genome.jl) | GRN development, mutation, inheritance |
| [src/params.jl](../src/params.jl) | every knob, mutable so the live view can retune |
| [src/population.jl](../src/population.jl) | struct-of-arrays population, slot allocator, clade table |
| [src/environment.jl](../src/environment.jl) | resource fields, diffusion, replenishment |
| [src/sim.jl](../src/sim.jl) | the per-tick loop |
| [src/stats.jl](../src/stats.jl) | stats log and snapshot IO |
| [src/run.jl](../src/run.jl) | headless driver |
| [scripts/run.jl](../scripts/run.jl) | CLI |
| [scripts/dashboard.jl](../scripts/dashboard.jl) | post-hoc charts |
| [scripts/live.jl](../scripts/live.jl) | live interactive window |

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
| 128 | 65,536 | ~15,800 | 282 MB | **886** |
| 256 | 262,144 | ~66,000 | 1.13 GB | **234** |
| 512 | 1,048,576 | ~272,000 | 4.5 GB | ~55 |

Grids 128 and 256 measured on an idle Apple M5 (10 threads) after the move to
variable genomes; 512 is from a live run on an M4 Pro and is the median of 57
samples spanning 44–132 ticks/s, since throughput swings ~3× over the seasonal
population cycle.

Memory scales as `NG_MAX²`, because every organism's weight block is dimensioned
for the largest genome it could have rather than the one it has — raising
`NG_MAX` from 20 to 32 multiplied the arrays by 2.6×. **Speed was not affected**
(886 vs ~900 tick/s at grid 128 before the change).

**Do not benchmark this under load.** A profile taken while other runs competed
for cores showed every phase ~4× slower, including `diffuse!`, which never
touches a genome and therefore cannot have been affected by the genome change.
That artifact was very nearly written up as a regression.

Current phase split, grid 256, idle, 10 threads: `organism_step!` 73%,
`settle_births!` 18%, `diffuse!` 6%, `replenish!` 3%.

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

### Controls, or: three hypotheses that died the same way

Three mechanisms were proposed, each fitted the data, and each was false. In
every case the giveaway only appeared because a control existed.

| hypothesis | why it looked right | what killed it |
|---|---|---|
| large genomes cause fragility | every extinct run had `ng` > 14.7, every survivor below 14.6, no overlap in 12 runs | pinning `ng` at 12 vs 20: extinction 1/12 vs 2/12, and on a continuous measure the ng=20 arm was nominally *less* fragile |
| genomes shrink during famine | an adaptive story that would neatly explain improved survival | mean `ng` during famine was lower in 5 of 10 runs and higher in 5 |
| famine *duration* drives bottleneck depth | deepest troughs landed 413–436 ticks after a shock, against `famine_len` 400 | matching total deficit: 400 × 0.65 and 150 × 0.05 are indistinguishable |

The pattern is the same each time: a **correlation in a small sample plus a
plausible mechanism**, which together feel like an explanation. The control is
what separates them, and in two of the three cases the control was cheap and
took twenty minutes.

Two specific traps worth naming:

- **Check that a "control" actually controls.** The famine-duration arms were
  run as a matched comparison and were not matched at all — integrated severity
  differed by 2.5×, which nobody would have noticed without writing out
  `len × (1 − factor)`. The confounded version gave a clean-looking 23× effect.
- **Binary extinction is the wrong statistic at these sample sizes.** 1/12 vs
  2/12 is Fisher p = 1.0, and separating 8% from 17% needs ~250 seeds per arm.
  Minimum population reached, as a fraction of the run's own median, measures
  the same fragility continuously and settles it at n=6.

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
- Lukeš et al., *How a neutral evolutionary ratchet can build cellular
  complexity* (Constructive Neutral Evolution) —
  [IUBMB Life](https://iubmb.onlinelibrary.wiley.com/doi/10.1002/iub.489)
- *Constructive Neutral Evolution 20 Years Later* —
  [PMC7982386](https://pmc.ncbi.nlm.nih.gov/articles/PMC7982386/)
- Lynch, *Drift and genome complexity* (mutational hazard hypothesis) —
  [PLoS Genetics](https://journals.plos.org/plosgenetics/article?id=10.1371%2Fjournal.pgen.1002092)
- Lynch, *Non-adaptive origins of interactome complexity* —
  [Nature](https://www.nature.com/articles/nature09992)
- *Evolutionary simulations reveal a role for genomic recombination in the
  evolution of gene regulatory network complexity and robustness* (2026) —
  [PMC13455014](https://pmc.ncbi.nlm.nih.gov/articles/PMC13455014/)
- *Open-Ended Evolution* (ALife encyclopedia) —
  [alife.org](https://alife.org/encyclopedia/introduction/open-ended-evolution/)
- Geritz et al., *Evolutionarily singular strategies and the adaptive growth and
  branching of the evolutionary tree* —
  [Springer](https://link.springer.com/article/10.1023/A:1006554906681)
- *Continuously stable strategies as evolutionary branching points* —
  [arXiv:1005.3862](https://arxiv.org/pdf/1005.3862)
- *Shaping of microbial phenotypes by trade-offs* (proteome allocation and the
  specialist/generalist trade-off) —
  [Nature Communications](https://www.nature.com/articles/s41467-024-48591-9)

## Genome complexity: does it grow?

It does not, and the experiment that shows it is the cleanest in the project.

Connectivity saturating at ~392 of 400 weights was originally read as
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

**Genome size decays in every arm here, monotonically, and connectivity decays
with it** (416 → 218 weights in a typical run). But these five arms do not
support the conclusion they look like they support, because **four of the five
charge for genes** and the fifth (`cost0`) settles at ~15.1 — about three nodes
*above* the floor, not at it.

A sixth and seventh arm settle it. Zero gene cost, started at the floor and
started high:

| arm | gene cost | started at | at 10k ticks | **final** |
|---|---|---|---|---|
| free12 | **0** | **12** | **14.5–18.1** | 12.2–15.1 |
| free20 | **0** | 20 | 20.1–26.1 | 16.0–21.0 |

**Started at the floor with free genes, genome size rises** — 12 to 14.5–18.1
within 10,000 ticks. So the floor is *not* absorbing; it only looked that way
because `lo12` carried a gene cost, which makes sitting at the minimum close to
tautological.

Two caveats that keep this from being a clean positive result either:

- **At `ng = NG_MIN` deletion cannot fire** — only additions are possible at the
  boundary. Part of the rise from 12 is therefore a mutational boundary effect
  rather than evidence that hidden units pay.
- **The two zero-cost arms do not converge** (final medians 13.3 from below,
  16.5 from above), so there is no demonstrated interior attractor. With
  symmetric add/delete rates and no cost, drift with a reflecting boundary would
  also produce roughly this, and 200k ticks may simply be too short.

### What can actually be claimed

- **Complexity does not grow open-endedly.** No arm, at any gene cost, from any
  starting size, showed sustained growth. Everything either decays or wanders in
  the mid-teens.
- **With any real gene cost, size decays to the minimum** and stays there.
- **With genes free, size does not collapse** — it sits in the mid-teens,
  reached from both directions but not convergently within 200k ticks.
- **Whether hidden units carry positive value is UNRESOLVED.** Distinguishing a
  weak selective benefit from drift-with-a-reflecting-boundary needs an arm that
  breaks the boundary asymmetry — for instance allowing deletion to fail
  harmlessly at the floor so add and delete stay symmetric there, or measuring
  the fitness of matched genotypes differing only in hidden-node count.

The earlier version of this section claimed complexity "decays toward an
absorbing floor even without a gene cost". That was wrong on both halves: it
conflated the zero-cost arm (which stops well above the floor) with the
start-at-the-floor arm (which had a cost). It is recorded here rather than
quietly deleted because the mistake is instructive — the arms were not a
factorial design, and the two factors got read off each other.

### Entrenchment, measured directly

Constructive Neutral Evolution rests on **entrenchment**: neutrally-added
components become indispensable over time, because the rest of the system comes
to depend on them. In real organisms that is inferred from phylogeny, because
you cannot delete a component and read off the cost as a function of how long it
has been there. Here you can — delete each hidden node, re-develop the network,
and measure how far the phenotype moves.

**Result: entrenchment is real.** 519,839 measurements across 24 seeds, 4,229
nodes tracked individually by unique id:

| | |
|---|---|
| nodes qualifying (≥5 observations spanning ≥20,000 ticks) | 248 |
| **positive slopes** | **173/248 (69.8%)** |
| median slope | +1.02×10⁻⁶ disruption per tick |
| effect over a 100k-tick lifetime | **+0.102 disruption units** |
| sign test | z = +6.22, **p = 4.9×10⁻¹⁰** |

Nodes present longer are measurably harder to remove.

#### Two earlier versions of this measurement were wrong

Both are recorded because the failure modes are the instructive part.

**Attempt 1 said entrenchment was absent** (r = −0.052). It sampled nodes only
at the *end* of a run, so every node observed was a survivor and **96% were
already older than 60,000 ticks**. With almost no young nodes to compare
against, the age axis was nearly constant and the correlation was noise around a
flat line. A further 20% of the sample were founder nodes, structurally
different from duplication-added ones, which dominated the oldest bin.

**Attempt 2 said entrenchment was strong** (r = +0.211, n = 67,068). It measured
during the run, fixing the age coverage, but grouped nodes into 5,000-tick birth
cohorts — mixing different nodes together — and pooled every measurement as if
independent. Most observations are repeat measurements of the *same* node, so
the correlation was inflated roughly by the number of re-measurements. Its own
paired cohort test was far weaker: 12/20 rising, p = 0.25.

**The unit of analysis is a node, not a measurement.** Attempt 3 gives each node
a unique id, follows it as it ages, and reduces it to a single regression slope.
That yields 248 independent points instead of half a million correlated ones,
and a clean answer. The inclusion rule — ≥5 observations spanning ≥20,000 ticks
— matters too: slopes fitted over a narrow age window are dominated by noise,
which is what flattened attempt 2's cohort test.

Caveat: only 248 of 4,229 nodes qualified, because most do not survive long
enough to yield a slope. That is survivorship again, though here it is
conservative — short-lived nodes are the ones least likely to have entrenched.

#### What entrenchment does and does not explain

It does **not** rescue open-endedness. Complexity still plateaus at ~22 while
entrenchment is demonstrably operating, so the ratchet exists but is **too weak
to overcome drift–selection balance** at this population size.

#### None of this is novel — checked

This was briefly written up here as a new finding. A literature check says
otherwise, and the check is recorded because the temptation to believe it was
new was strong.

**The method is standard.** Protein evolution has measured entrenchment this way
for years: *"entrenchment is measured by reverting a mutation to its original
amino acid at times posterior to the mutation's occurrence, and observing
whether this reversion is typically deleterious at subsequent times."* That is
exactly the delete-and-measure-against-age design above. See [Emergent time
scales of epistasis in protein
evolution](https://www.pnas.org/doi/10.1073/pnas.2406807121) and [Hsp90 over a
billion years](https://pmc.ncbi.nlm.nih.gov/articles/PMC5924896/).

**The result is established.** Entrenchment rising with evolutionary time is the
expected finding, and even the "weaker than expected" nuance is published —
[*Epistasis plays a limited role in driving entrenchment during neutral protein
evolution*](https://link.springer.com/article/10.1186/s13059-026-04157-0)
(Genome Biology, 2026).

**The conclusion is published too.** Walsh (2026), [*Neutrally Evolving
Interlocking Complexity in the Quandary Den*](https://arxiv.org/html/2604.18361),
simulates CNE directly and concludes the ratchet is self-limiting and complexity
plateaus. Their limiting mechanism differs — shrinking solution space rather
than drift–selection balance — and they run a population of **one**, so drift
plays no part. But the headline is theirs.

What is left unclaimed is narrow: Walsh does not examine population size or
drift, and the protein work measures sites within a fixed-length sequence rather
than components whose *number* evolves. Age-resolved entrenchment under evolving
component count with real population dynamics does not appear to have been done.
That is a gap-filling increment, not a finding, and it does not support a paper.

### What makes complexity accumulate: a 2×2×2 factorial

The absence of growth is not mysterious. Constructive Neutral Evolution names
four conditions under which complexity fails to accumulate, and this model
violated two outright: **selection too strong** (populations of 20,000–30,000
purge neutral additions) and **mutation unbiased** (`p_node_add` equalled
`p_node_del` exactly, giving no directional ratchet). Lynch's mutational hazard
hypothesis says the same from the other side — genome expansion needs small
effective population size so drift can overwhelm purifying selection.

Both are one-parameter changes. Full factorial, 3 seeds per cell, 150k ticks,
reporting final mean genome size:

| small pop | mutation bias | gene cost | final `ng` |
|---|---|---|---|
| — | — | yes | **13.0** |
| — | ✓ | yes | 13.4 |
| ✓ | — | yes | 14.8 |
| ✓ | ✓ | yes | 15.0 |
| — | — | no | 15.1 |
| — | ✓ | no | 17.0 |
| ✓ | — | no | 20.0 |
| **✓** | **✓** | **no** | **21.5** |

**Complexity is sustained at 21.5 against a baseline of 13.0** — 1.7×, held flat
across 150k ticks rather than decaying. This is the first configuration in which
genome size does not collapse.

**The gene cost gates everything.** With a real per-gene cost, selection purges
additions regardless of anything else: small populations buy **+1.8** and
mutation bias buys **+0.4**. Remove the cost so additions are neutral and the
same two factors become potent: small populations buy **+4.9**, bias **+1.9**.

That gating is the whole result, and it is mechanistically obvious in hindsight.
**Mutational bias changes the supply of additions, not their fate.** If every
extra gene carries a fitness cost, selection removes them however fast they
arrive — which is why the CNE condition about insertion bias does nothing here
until the additions are actually neutral. Drift can only preserve what selection
is indifferent to.

So the earlier framing — that the environment is too simple to reward complexity
— was at best half the story. The population was also far too large and every
gene was taxed. Under **large populations with efficient purifying selection,
streamlining is the predicted outcome**, and the literature says so explicitly:
prokaryotic ribosomes are streamlined precisely because prokaryotic populations
are huge. This model's organisms are prokaryote-like with enormous populations,
so shedding genome was never a failure — it was the textbook result.

The likeliest explanation for the absence of growth remains that **this
environment does not demand regulatory complexity**: mapping six inputs to
twelve traits is simple enough that the interface layer largely suffices. But
that is now a hypothesis the data is consistent with, not a result the data
establishes.

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

## Making the environment demand complexity — and why it could not be done

Complexity plateaus at ~22 nodes. Walsh (2026) attributes that to a shrinking
solution space, Lynch to population size. A third possibility appears untested:
maybe the plateau is set by the *task*, and a job that actually requires a
conditional response would move it. `yield_zones` was built to test that — the
dish splits into vertical bands and each band rotates which channel pays best,
so the *set* of yields is identical everywhere and only the context changes.

The answer is that **the hypothesis remains untested, because the environment
could not be made to demand complexity.** Genome size did move, substantially
and reproducibly, but nothing ties that to the task.

### The first version failed for an obvious reason

Rotating the yields alone does nothing. Channel 1 is the sole recipient of
external inflow, so it is far more abundant everywhere, and a yield rotation
spanning 1.0 → 0.22 cannot outweigh that. Organisms went on eating channel 1 in
every zone and simply paid for it: genome size unchanged, population 958 → 505.
The dish became uniformly worse rather than context-dependent. `rotate_inflow`
fixes the alignment so the abundant channel and the best-paying channel coincide.

### With the alignment fixed, genome size rises — monotonically

Six seeds per condition, 150k ticks, complexity-permissive regime (grid 32,
`maint_gene = 0`, `p_node_add = 0.006`, `p_node_del = 0.002`):

| zones | final `ng` | SE | population |
|---|---|---|---|
| 1 | 21.81 | ±0.51 | 958 |
| 2 | 23.96 | ±0.58 | 1013 |
| 4 | 24.82 | ±0.96 | 1263 |
| 6 | **26.69** | ±1.17 | 1297 |

`ng ~ log2(zones)`: **+1.74 per doubling** (t = 4.17, R² = 0.44). Spearman
ρ = +0.66, permutation p = 0.00066. Controlling for population the slope is
+1.86 (t = 2.90) while `log2(pop)` contributes nothing (t = −0.24). Weights per
node rise too, 21.0 → 25.5.

> **Retracted.** That last inference does not hold. `nzw` counts weights that are
> non-zero, and point mutations are *additive* (`W[m,k] += d`), so a weight that
> once becomes non-zero never returns to exactly zero. Density therefore ratchets
> upward with lineage depth regardless of selection — measured: 0.25 at founding,
> 0.777 by depth 304, 0.862 by depth 1396. "Weights per node" tracks how long a
> lineage has been mutating, not how densely it is wired, so it cannot support
> "denser rather than padded". See the connectivity section below.

It is a **higher plateau, not open-ended growth** — every condition flattens by
25–50% through the run and holds. And it is not a ceiling artifact: re-running
zones = 6 with `NG_MAX` raised 32 → 64 gives 26.01 vs 26.69 (t = −0.43), so the
bound never bound.

### Two findings that take it all back

**Organisms do not respond conditionally.** The decisive test takes the *same*
organism and develops it under each zone's environment in turn. One encoding a
fixed optimum gives the same uptake profile whatever it is shown; one encoding a
conditional map shifts. Averaged over 3 seeds:

| evolved in | final `ng` | switchers | profile spread | max shift |
|---|---|---|---|---|
| 1-zone world | 21.70 | 40.7% | 0.185 | 0.491 |
| 6-zone world | 27.12 | 11.9% | 0.086 | 0.354 |

Organisms from the zoned world are **less** plastic, not more, on every measure
and in 2 of 3 seeds each. In one seed the 6-zone population's top channel was
`ch6 ch6 ch6 ch6 ch6 ch6` — identical across all six probes. The bigger genomes
are not doing conditional work.

**And the control arm produced the largest genomes of all.** Holding
`rotate_inflow` on and switching off only `rotate_yield` was meant to isolate the
payoff rotation. It does not: without the yield rotation the abundant channel in
each zone pays 0.22–0.35 instead of 1.0, so most of the dish is poor.

| cell | population | standing resource | final `ng` |
|---|---|---|---|
| zones 2, yield rot **no** | 516 | 5342 | **28.72** |
| zones 2, yield rot yes | 1013 | 2347 | 23.96 |
| zones 4, yield rot **no** | 371 | 7599 | **27.10** |
| zones 4, yield rot yes | 1263 | 1968 | 24.82 |
| zones 6, yield rot **no** | 452 | 7174 | 25.23 |
| zones 6, yield rot yes | 1297 | 2749 | 26.69 |

Resource piles up unconsumed, population falls 2.5–3.5×, and the arm with the
*less* conditional task carries the **largest genomes in the dataset** — with a
trend running *downward* in zones (28.7 → 27.1 → 25.2). Smaller population, more
drift, more complexity. That is Lynch again, arriving as a confound. Pooled
across both arms nothing survives: `log2(zones)` t = −0.34, `yield_rot` t = −0.95,
R² = 0.09.

### Why the task could not be built

The zones never present distinct readable optima, and the reason is the
conversion ladder. Channel *r*'s byproduct is *r+1*, always. An organism eats the
inflow channel down to nothing and excretes one step along, so **the realised
concentration peak sits systematically downstream of the channel actually being
supplied**, and waste accumulates in ch5/ch6 everywhere. Expected top channel by
zone is `ch1 ch6 ch5 ch4 ch3 ch2`; what the field actually shows is
`ch3 ch6 ch5 ch6 ch6 ch6`.

This is not a diffusion problem. Lowering `diff_rate` 6× (0.12 → 0.02),
tripling the grid (32 → 96), and adding per-channel falloff all leave it at 1–2
zones out of 6. Raising `decay` 13× (0.006 → 0.08) reaches 3 of 6 — at the cost
of cross-feeding (downstream uptake share 0.846 → 0.704) and population
(1608 → 1115); at 0.20 the dish goes extinct.

**So in a model with a fixed downstream excretion ladder, a spatially varying
optimum cannot be made legible from local concentration.** Consumption plus
excretion displace the peak away from the supply. Fixing that means dismantling
the cross-feeding ladder — which exists because an evolvable excretion
destination produced a futile cycle and killed cross-feeding outright.

The `ng` rise is therefore real, reproducible, and **unattributable**: most
plausibly a side effect of the zoned dish changing selection intensity, not
evidence that complexity was demanded. `rotate_inflow` and `rotate_yield` are
kept so the experiment can be redone if the ladder is ever replaced by a
metabolic matrix that does not pile waste into a terminal sink.

### A cue that cannot be eaten

The zoned task failed because organisms had to read a field they were actively
flattening. So the cue was moved out of the resource system entirely: node
`P_SENSE` (index 13) receives an **abiotic condition** that rotates which
channel pays best, is not consumable, and is fed straight to the genome. With
`cond_states = 1` the cue is constant and development is bit-identical to the
model without the feature, so the control arm is the old model rather than a
differently-tuned one.

Development happens once at birth, so an individual is committed to the state it
was born into. The *genome* therefore has to hold a cue → phenotype map that
works in every state; one that ignores the cue is right 1/`cond_states` of the
time.

**First attempt: the cue was legible but the reward was not.** `ng` came out
21.05 / 22.44 / 21.07 / **27.64** for 1/2/3/6 states — flat, then a jump at 6.
But the diagnostic showed `matched = 1/6` in *every* arm: the 6-state
populations answered `ch1 ch1 ch1 ch1 ch1 ch1` against a target of
`ch1 ch6 ch5 ch4 ch3 ch2`. Channel 1 holds the inflow monopoly, so it stays the
most abundant channel whatever the cue says and "always eat ch1" remains good
enough. Only 6 states moved anything because it is the sole setting where every
channel receives every yield exactly once, so no fixed diet keeps an edge.

**Second attempt: routing inflow to the cue channel killed the dish.** 48 of 48
runs extinct, median at the first or second switch (~2.5–4.6k ticks). The cause
is architectural: with development cached at birth, a discrete global switch
mis-adapts every adult simultaneously and intake collapses before new births can
adapt. Fractional inflow did not rescue it — at `cond_inflow_frac = 0.2` only
2 of 4 survived, and the survivors sat at `ng` = 32.00, exactly `NG_MAX`, which
is a tiny population drifting into the cap rather than a result.

**What worked: a travelling wave.** With `cond_gradient` the cue varies in space
as well as time, so the dish is never globally starved — some region always
suits somebody — while the drift means no lineage can park in a favourable band.
At `cond_inflow_frac = 0.4` this survives 4/4.

### The control that finally controls

Every comparison in this project had been spoiled the same way: the
manipulation changed population size as well as the task, and population size
drives genome size on its own. The fix is a **scrambled cue** — same rotation,
same inflow, same energy, same ecology, same population, with only the
*information* removed. Measured mutual information between cue and state is
**2.585 bits informative against 0.0000 bits scrambled**, at identical range and
smoothness. (A golden-ratio scramble was tried first and leaked 0.72 bits,
because cue and state then share the spatial term; 7 temporal turns against 5
spatial ones gives a clean zero.)

16 seeds per arm, 150k ticks:

| arm | final `ng` | SE | population |
|---|---|---|---|
| flat — no rotation | 21.36 | ±1.40 | 1020 |
| **informative cue** — 2.585 bits | 25.44 | ±1.03 | 497 |
| **scrambled cue** — 0.000 bits | 26.09 | ±1.02 | 469 |

Informative against scrambled: **−0.65, t = −0.45, permutation p = 0.65**, with
populations matched at 497 and 469.

Both wave arms sit well above flat — but **the arm whose cue predicts nothing
rises just as much as the arm whose cue predicts everything.** The increase is
the ecology, not the task. Information contributes zero.

Caveat: 3/16, 3/16 and 1/16 seeds finished at or above `ng` 31 against a cap of
32, so the top of the distribution is mildly compressed. That is symmetric
across arms and does not bias the comparison.

**The control is the only one in this project that holds ecology constant**, and
on its own terms it is clean: the rise is ecology, not task.

### But the null does not mean what it looks like — checked, and it failed

The conclusion above was written first as "a task that genuinely requires a
conditional response does not move the plateau." That claim needs organisms in
the informative arm to have actually *read* the cue. They did not.

Probing each organism with all six cue values, resources held uniform:

| arm | switchers | profile spread | matched | top channel by state |
|---|---|---|---|---|
| informative | 0.0 / 0.2 / 0.0 / 0.4 % | 0.157 / 0.047 / 0.019 / 0.003 | 1/6 all seeds | `ch1 ch1 ch1 ch1 ch1 ch1` |
| scrambled | 6.1 / 0.0 / 0.0 / 8.7 % | 0.018 / 0.000 / 0.037 / 0.099 | 1/6 all seeds | `ch1 ch1 ch1 ch1 ch1 ch1` |
| flat (no cue at all) | 96.8 / 15.7 / 0.1 / 0.0 % | 0.360 / 0.113 / 0.076 / 0.031 | 1/6 all seeds | mixed |

**Nobody reads the cue.** The informative arm shifts its uptake profile by
essentially nothing across the full range of the signal, and answers the same
channel in all six states. It is in fact *less* cue-sensitive than the flat arm,
whose organisms never experienced a varying cue at all — what little response
the flat arm shows is an unused node being perturbed, not adaptation.

So the comparison was **an ignored cue against a useless one**, which are the
same situation. The null is real but uninformative about the hypothesis: it
shows that an unread signal changes nothing, which was never in doubt.

What the experiment does establish is narrower and still worth recording:
**conditional response did not evolve at all** — not with a cue that is abiotic,
undepleteable, and carries 2.585 bits about the rotation, and not when tracking
it should pay. A rough accounting says the conditional strategy is worth ~38%
more than the fixed one at `cond_inflow_frac = 0.4` (40% of inflow at yield 1.0,
against 60% of inflow at a mean yield of ~0.48), so the advantage is there and
goes untaken.

Two candidate reasons, neither tested. The GRN may have no incremental path to a
conditional map — partial responsiveness could be worse than none, making it a
fitness valley rather than a slope. Or the useful fraction of inflow may be
unreachable: raising `cond_inflow_frac` past 0.4 is what makes the advantage
decisive, and 0.4 is already near the edge where the dish dies (at 1.0, 48 of 48
runs went extinct). If that window is genuinely empty, the model cannot pose
this question at all.

### Why it goes unread: conditional specialisation is deleterious

Rather than guess, the conditional genome was **built by hand and competed**.

**The architecture can represent it perfectly.** Five hidden units on top of the
sensor — threshold units at the band boundaries, each uptake channel driven by
the difference of two adjacent steps — give a **6/6 match**, with the allocation
essentially one-hot on the correct channel at every cue value. Representability
is not the obstacle.

**But it loses, and not because mutation erodes it.** Seeded at 20% into an
established population, the conditional type is gone within ~500 ticks. With
**mutation switched off entirely** it still disappears, and its energy runs at
**0.52–0.64 of its neighbours'**. That rules out mutational erosion and
unreachability alike: the behaviour is simply worse.

A tilt sweep says the same thing monotonically — 50/50 competition, mutation off,
three seeds, final carrier frequency per seed:

| tilt `A` | per-seed outcome | reading |
|---|---|---|
| 0.0 (no-op control) | 1.0 / 0.0 / 0.0 | drifts to fixation or loss — neutral |
| 0.5 | 0.0 / 1.0 / 0.0 | drifts — no detectable benefit |
| 1.5 | 0.0 / 0.0 / 0.07 | losing |
| 4.0 | 0.0 / 0.0 / 0.0 | eliminated |

The stronger the conditional response, the worse it does. No tilt is beneficial.

**The mechanism.** Uptake is an *allocation* across channels and capture is
proportional to local availability: `got = avail × uptake_rate × aff^exp`.
Concentrating the budget on the best-paying channel forfeits access to the other
five — including the waste channels the ladder keeps stocked — and the yield
premium (1.0 against a mean of ~0.48) does not cover the abundance lost. Breadth
beats rate. The earlier "conditional is worth ~38% more" estimate counted the
yield gain and ignored the opportunity cost, and was simply wrong.

#### Three grafts failed their own control before one passed

Worth recording, because each failure was invisible until the control was
checked. The first graft zeroed the uptake nodes' incoming weights — which also
changed what those nodes *output* to nodes 7:12, wrecking every other trait
(`repair` 0.995 → 0.127, `leak` 0.002 → 0.723, `divide` 0.002 → 0.695). It was
measuring a broken organism. The second cut that coupling and re-tuned the
downstream biases, preserving all other traits to within 0.001 — that one is
sound, and it is the basis of the 0.52–0.64 result. The third, for the tilt
sweep, cannibalised five evolved hidden nodes, and its A = 0 control went to
**zero** when a no-op must hold at 0.5. Only a graft that *extends* the genome
with five fresh nodes is bit-identical at A = 0 (max phenotype change
0.000000), and only then is the sweep interpretable.

An index-based lineage tag also had to be thrown away: dead slots are recycled,
so the tag leaked to unrelated organisms and the A = 0 control reached fixation
in every arm. A genome marker fixed it. **Every one of these was caught by a
control that had a known correct answer** — which is the argument for always
including one.

### Rebuilding uptake, and the trade-off that closes the question

The obvious fix was to stop specialisation being dominated by breadth.
`uptake_half` adds Monod saturation — off by default, bit-identical to
proportional capture at 0, with a test pinning that.

**It is the wrong lever, in an instructive direction.** Saturation compresses
channel values *toward each other*: a channel ten times richer stops being worth
ten times more, so spreading across many thin channels gains relative value. It
rewards generalism. The seeded conditional genome did *worse* with saturation on
than off — energy ratio 0.388 against 0.833 — and measured specialisation moved
the same way (top-share 0.76 → 0.70).

What decides this is the shape of the **allocation trade-off**, not the shape of
the concentration response. Levins: a convex trade-off selects specialists, a
concave one generalists. `aff^uptake_exponent` with the exponent above 1 makes
concentrating the budget pay superlinearly, and that knob already existed.
At 1.5 it produces genuinely specialist populations (top-share 0.955) that still
survive the travelling wave.

**And in that regime, information finally has a measurable fitness effect.**
Carrier energy against non-carriers, averaged over ticks 250–2000, 12 seeds:

| cue | energy ratio |
|---|---|
| informative (2.585 bits) | **0.896** |
| scrambled (0.000 bits) | **0.492** |

**+0.404, t = 3.88, permutation p = 0.0006**, ecology matched. Reading the cue
nearly doubles a carrier's energy. (The scrambled *n* is 9 rather than 12 because
three seeds lost their carriers before the first sample, which biases the
scrambled mean *upward* — so the true gap is wider.)

It is still below 1.0, so the carrier remains net deleterious and is eliminated.

### The window is empty

Mapping payoff against viability across 18 cells — `uptake_exponent` ×
`cond_inflow_frac` × grid, 2 seeds each:

| `uptake_exponent` | burn-in survival | mean energy ratio |
|---|---|---|
| 1.0 (original) | 12/12 | 0.707 |
| 1.5 (convex) | 12/12 | 0.929 |
| 1.7 | **7/12** | 0.898 |

Best single cell 0.976. **Not one of the eighteen reaches parity.** The payoff
climbs steeply with convexity and then viability collapses before it can cross:
`exp` 1.9 survives 0/8 burn-ins. The conditions that make specialisation pay and
the conditions that impose a moving target do not overlap — the more
specialisation pays, the more lethal a moving target becomes.

A tilt-strength sweep under the convex trade-off says the same, 12 seeds,
carriers seeded at 0.5, frequency after 12,000 ticks:

| tilt `A` | final frequency | vs no-op |
|---|---|---|
| 0.0 (verified no-op) | 0.286 | — |
| 0.5 | 0.446 | +0.160, p = 0.27 |
| 1.0 | 0.460 | +0.173, p = 0.26 |
| 1.5 | 0.225 | −0.061, p = 0.57 |

Nothing reaches significance, and the effect is non-monotonic — 1.5 goes the
other way. Drift dominates: most seeds fix or lose outright, so detecting a
difference this size would need on the order of 60 seeds per arm.

Two caveats recorded against these numbers. The no-op control reads an energy
ratio of **1.070** rather than 1.000, because carriers are seeded by slot index
and slot order correlates with age and lineage — a non-random sample. The
comparison is still fair (every arm carries the same bias) but the absolute
ratios are not. And a 3-seed pilot showed A = 1.0 climbing to 0.688, which did
not survive 12 seeds — the second time in this work that a small-*n* signal
evaporated, after a 4-seed invasion result did the same.

### The axis that was never varied: environmental grain

All of the above searched the *payoff* axis — uptake trade-offs, yield premiums,
inflow fractions. A literature check found the binding constraint was somewhere
else entirely.

Levins' framework turns on the **grain** of environmental variation relative to
generation time. The [Chlamydomonas selection
experiments](https://www.nature.com/articles/6883290) state the rule directly:
*versatile generalists capable of reversible responses evolve when the
environment is fine-grained, whereas plastic generalists that respond
irreversibly to the conditions of growth evolve when it is coarse-grained.*
Development here is cached at birth and never updated — irreversible by
construction — so it needs a coarse-grained world. Measuring whether it had one:

| quantity | measured |
|---|---|
| ticks per environmental state | 333 |
| median lifespan | 346 ticks |
| **fraction of life spent in the birth state** | **0.334 mean, 0.250 median** |
| organisms never leaving their birth state | **6.2%** |
| organisms spending under half their life in it | **80.8%** |
| evolved motility / band width | 0.960 / 10.7 cells |

Organisms streamed across bands continuously. [Moran's
condition](https://pmc.ncbi.nlm.nih.gov/articles/PMC6369965/) says the same in
other words — plasticity evolves only when the environment of *development*
predicts the environment of *selection*, and when they decouple *the evolved
plastic response is shallower than the perfect response*. The tilt sweep found
exactly that shape: shallow tilts best at ~1.017, steep ones falling to 0.719.
That was read as noise at the time; it is the signature of an unreliable cue.
High motility with narrow bands is also textbook [gene flow swamping local
adaptation](https://onlinelibrary.wiley.com/doi/full/10.1111/mec.13606).

### Both failure modes, bracketed

Removing the spatial gradient and slowing the wave produced a genuinely
coarse-grained world — **96% of life in the birth state**. The conditional
genome was still eliminated, 12/12 and 9/9. The reason is the *other* failure
mode: at ~21 generations per state the environment holds still long enough for
ordinary genetic adaptation to track it, so plasticity is unnecessary. Kussell
and colleagues measure this directly — [selection for regulation weakens as
generations per environment rises](https://www.biorxiv.org/content/10.1101/2025.06.23.661119v2.full).

So plasticity is bracketed: **too fine and the cue cannot predict; too coarse
and genetics suffices.** Mapping the axis between them, 3 seeds per cell:

| period | viable | grain | generations/state | payoff |
|---|---|---|---|---|
| 3000 | 5/6 | 0.960 | 15.1 | 0.954 |
| 6000 | 4/6 | 0.982 | 13.2 | 0.705 |
| 12000 | 6/6 | 0.987 | 26.9 | 0.967 |
| 24000 | 6/6 | 1.000 | 92.7 | 0.845 |
| 48000 | 6/6 | 1.000 | 47.4 | **0.991** |

No systematic dependence on grain — the payoff hovers just under parity
throughout, and the fast cells start dying.

### The strongest result available — and an audit that cut it down

The best cell (period 48000) re-run at 12 seeds per arm. **Four errors in the
first write-up were found by auditing it; all four are recorded here because
three of them changed the number and one changed the conclusion.**

| | informative | scrambled | **no-op (marker only)** |
|---|---|---|---|
| energy ratio, mean | 1.103 | 0.429 | **0.902** |
| energy ratio, median | 1.099 | 0.257 | **0.999** |
| **carriers retained** | **7/9** | **3/12** | **12/12** |

**Error 1 — the Fisher table counted extinctions as retentions.** Reported as
9/12 vs 3/12; the informative arm actually had 2 burn-in extinctions and 1
post-seeding population collapse, so "12 − 3 carrier-losses = 9" was wrong.
Restricted to runs whose population survived, the correct table is 7/9 vs 3/12:
**Fisher exact p = 0.0300**, slightly *stronger* than the 0.0391 first reported.
The retention result survives.

**Error 2 — the energy mean rests on one seed.** Values are 0.983, 1.105, 1.001,
1.149, 0.937, **1.444**, 1.100. Dropping the outlier gives 1.046. Median 1.099.

**Error 3 — generations-per-state was wrong by a factor of ~3.** It was computed
from a forward-tracked median lifespan of 95 ticks, which oversamples long-lived
survivors. Measuring generation time properly, as the slope of mean lineage depth
against tick, gives **26.1 ticks**. So the real figures are 14 generations per
state at period 2000, **77** at 12000, and **307** at 48000 — not 21 and 47.
Every configuration tested was temporally *coarse*.

**Error 4 — the no-op control was never run for this configuration.** It matters
because carriers are seeded by slot index, which correlates with age and lineage.
Running it (additive graft at A = 0, phenotype change 3.8×10⁻⁶): mean 0.902,
median 0.999, and **12/12 carrier retention**.

That last one changes the reading. Retention orders as **no-op 12/12 >
informative 7/9 > scrambled 3/12**. The conditional genome is retained *less
well than a marker carrying no function at all*.

So the honest statement is narrower than "information helps":

> Carrying the conditional machinery costs something. An informative cue recovers
> most of that cost — decisively more than a scrambled cue (energy permutation
> p = 0.0084; retention p = 0.0300) — but it does not make the machinery better
> than not having it.

### The predicted window is lethal — now tested rather than assumed

With generation time known to be 26 ticks, the regime theory actually points at —
one to four generations per environmental state, fast enough that genetic
tracking cannot keep up — turns out never to have been tested. It is:

| period | generations/state | viable | payoff |
|---|---|---|---|
| 300 | 1.9 | **0/6** | — |
| 600 | 3.8 | **0/6** | — |
| 1200 | 7.7 | 3/6 | 0.53–0.75 |
| 2000 | 12.8 | 3/6 | 0.92–0.98 |
| 12000 | 76.6 | 6/6 | 0.967 |
| 48000 | 306.9 | 6/6 | 0.991 |

**The window where plasticity should win is uninhabitable.** Everything
survivable is 13+ generations per state, where genetic tracking suffices. And the
best payoff sits at the *slowest* end — 307 generations per state — which is
where plasticity should matter least, so the payoff gradient is not about
plasticity at all.

One further correction to the grain diagnosis: holding position fixed, the
*temporal* cue agreement is **0.907**. The 0.334 figure is produced almost
entirely by **movement** — motility evolved to 0.99 across bands 10.7 cells wide.
The environment is fine-grained spatially and coarse-grained temporally, which is
the worst combination: an unreliable cue *and* slow enough change that genetics
can track it.

What is *not* closed, and would be the honest next step: every result here tests
irreversible development. Theory says a fine-grained world of this kind calls for
**reversible** plasticity — re-development during life — which this architecture
cannot express, and which the 2025 *Evolution Letters*
[treatment](https://academic.oup.com/evlett/article/9/4/491/8193125) explicitly
flags as unaddressed for within-lifetime change. That is a change to the
development model, not a parameter sweep.

## The hypothesis, finally tested properly

Five attempts failed before the experiment could even run: the cue was
illegible, then unread, then the behaviour it enabled was selected against, then
the regime where it should pay was uninhabitable. Reversible plasticity
(`redevelop_every`) removed the last of those, and the test finally ran.

**Regime A: 100 seeds per arm**, informative against scrambled cue, ecologically
identical, 150k ticks, `NG_MAX` = 64 so nothing is censored.

| measure | scrambled | informative | p |
|---|---|---|---|
| **`ng`** | **30.55** | **30.81** | **0.836** |
| switchers | 0.229 | **0.645** | **< 10⁻⁵** |
| profile spread | 0.107 | **0.374** | **< 10⁻⁵** |
| matched | 0.976 | **1.545** | **< 10⁻⁵** |

n = 88 and 83 after extinctions. **Organisms decisively evolve to read the cue,
and genome size does not follow.**

The cue-usage result is anchored by a third arm with no cue at all
(`cond_states = 1`, nothing varies): switchers 0.133, spread 0.060, matched
**exactly 1.00/6** — which is the non-responsive baseline derived analytically,
since an unresponsive population answers one fixed channel and each channel is
best-paying in exactly one of six states. The scrambled arm is statistically
indistinguishable from that no-cue baseline (p = 0.26, 0.34), which is correct:
neither carries usable information. Only the informative arm departs from it.

### The null is now a measurement, not an absence

The first version of this test used 12 seeds, and a power analysis showed it
could only have detected an **8.7-node (27%) rise** — a 95% CI of [−6.4, +5.7].
That is not evidence of absence. At 100 seeds per arm:

| | n = 12 | **n = 100** |
|---|---|---|
| minimum detectable rise (80% power) | 8.67 nodes (27%) | **3.41 nodes (11.1%)** |
| 95% CI on the difference | [−6.44, +5.70] | **[−2.13, +2.64]** |

**Any rise larger than 2.64 nodes is excluded by the data.** Sample size also
turned the cue-usage effects from marginal (p = 0.006–0.09 at n = 12) into
decisive — they were real all along and merely underpowered, which is the
argument for doing this properly the first time.

### The sharpest result: solving the task does not require a bigger genome

Between-arm comparisons can always be argued about. The within-arm test cannot:
among the 88 populations that all received the *same* informative cue, does
responding more strongly go with carrying more regulatory nodes?

| arm | corr(switchers, `ng`) | n | p |
|---|---|---|---|
| informative | **−0.189** | 88 | 0.077 |
| scrambled | +0.049 | 83 | 0.665 |

If anything it is negative. So the finding is not merely that the task failed to
raise complexity — **the populations that solved it best did not need extra
machinery to do so.** Conditional response here is achieved by rewiring, not by
adding parts.

At n = 12 this test was uninformative, giving +0.481 in one batch and −0.040 in
another. It took the larger sample to resolve.

**Conclusion: task difficulty does not set the complexity plateau.** Population
size and gene cost do, which is where every independent line of evidence in this
README arrives.

## Fluctuation raises complexity — a reproduction, plus one new control

After five failed attempts to show that a *task* demanding conditional response
raises genome complexity, the controls that were built to test that hypothesis
answered a different question, and answered it cleanly.

**Read the novelty audit below before quoting any of this as new.** The main
effect here is published — three times over — and what this work adds is one
control the prior studies did not run.

### The result

| condition | mean `ng` | n |
|---|---|---|
| **frozen** — rotation machinery on, state never changes | **23.77** | 20 |
| no cue at all (`cond_states` = 1) | 23.43 | 12 |
| **fluctuating, unpredictable** (scrambled cue, 0.000 bits) | **30.55** | 83 |
| **fluctuating, predictable** (informative cue, 2.585 bits) | **30.81** | 88 |

**Environmental fluctuation raises genome complexity by ~7 nodes (~30%).
Whether that fluctuation is predictable makes no difference at all.**

- frozen → scrambled: **+6.78, p = 0.0013**
- frozen → informative: **+7.04, p = 0.0002**
- scrambled → informative: +0.25, **p = 0.836**

### Why the frozen control is what makes this a result

Every earlier comparison in this project was spoiled the same way: the
manipulation moved the ecology as well as the thing being tested. The
no-fluctuation arm (`cond_states` = 1) has that flaw — it also changes how
inflow is routed.

The frozen arm removes it. `cond_period` = 10⁹ leaves the rotation machinery
running and routing inflow exactly as the fluctuating arms do, while the state
never changes within a run (verified: one state visited, cue constant at −1.0).
It lands at 23.77 — **statistically identical to the no-cue arm (p = 0.88)** and
~7 nodes below both fluctuating arms. So the effect is fluctuation itself, not
the inflow structure.

### It holds across a five-point heterogeneity gradient

Reproducing the manipulation from [*Adaptive Dynamics of Regulatory Networks:
Size Matters*](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3171436/), whose
variable H = 2–8 is the direct analogue of `cond_states`:

| states | `ng` informative | `ng` scrambled | p |
|---|---|---|---|
| 2 | 26.75 | 27.79 | 0.61 |
| 3 | 26.59 | 28.03 | 0.57 |
| 4 | 29.51 | 26.76 | 0.30 |
| 6 | 31.15 | 30.45 | 0.74 |
| 8 | 29.81 | 29.06 | 0.81 |

`ng` rises with the number of states — **+2.13 nodes per doubling informative,
+1.04 scrambled** — which is the direction that paper reports. But the
informative-scrambled difference is null at **every single H**. The rise tracks
how much the environment varies, not how solvable the variation is.

That reconciles this work with the prior result rather than contradicting it.
The 2011 study varied heterogeneity and saw optimal network size rise; with no
zero-information control, reading that as the task demanding complexity was
natural. Add the control and the same rise appears when the cue predicts nothing.

Caveat: the 8-state cell has switchers at 0.014 — organisms stopped tracking the
cue there entirely, and `ng` fell from 31.15 to 29.81. The gradient is not
monotone at the top, and the slope rests on the 2→6 range.

### The task result, restated correctly

Across three independent regimes, an informative cue reliably produces
conditional *behaviour* and no extra genome:

| regime | `ng` scram → info | p | switchers | p |
|---|---|---|---|---|
| A: exp 1.5, per 600, 6 states | 30.55 → 30.81 | 0.836 | 0.229 → 0.645 | <10⁻⁵ |
| B: exp 1.0, per 1200, rd 5 | 22.46 → 21.73 | 0.357 | 0.103 → 0.433 | <10⁻⁵ |
| C: exp 1.5, per 600, 3 states | 29.62 → 26.43 | 0.036 | 0.109 → 0.226 | 0.0004 |

Detection floors 3.41, 2.21 and 4.19 nodes, so these are measured nulls rather
than underpowered ones. Regime C is nominally significant and *negative*, but
does not survive Bonferroni across three regimes (needs p < 0.017). Within-arm
correlations between how well a population tracks the cue and how large its
genome is: −0.189, +0.066, −0.062 — scattered around zero.

**So: complexity responds to environmental change, not to solvable environmental
change.** Behavioural complexity and genomic complexity dissociate.

### Two weaker results, stated as such

**The sensor knockout does not isolate a mechanism.** Ablating `P_SENSE` costs
~15% energy and drives carriers from 50% to 5.8% (informative) and 0.6%
(scrambled), so it is clearly deleterious — but in *both* arms, with no
significant difference between them (p = 0.64). Cutting a node's outgoing
weights disrupts the network whether or not the cue matters. A sham ablation on
a random non-sensor node is the control this needed and it was not built.

**A bigger world barely rescues the fast regime.** At 1.9 generations per state:
0/3 viable at grid 64, 0/3 at 128, 2/3 at 192. Tripling the world helps only
marginally. And `ng` *falls* at grid 192 (27.8 → 23.1), which is the Lynch effect
reasserting itself — larger populations, more efficient selection, smaller
genomes.

## Did the states have to share structure? No.

The most plausible explanation for the task null was that the six environmental
states were **cyclic rotations of one yield vector** — six views of a single
mapping. Knowing one gives you all six, and "find the abundant channel and eat
it" covers every case, so no extra genome would be needed however many states
there are. Kashtan & Alon and the continual-learning literature both point at
shared structure as the relevant variable.

`independent_states` removes the sharing: each state gets its own random
permutation, preserving the yield set (and therefore the energy on offer)
exactly. Under rotation state 2 is state 1 shifted; under independent
permutations they are unrelated.

**2×2, structure × information, 30 seeds per cell:**

| cell | n | extinct | `ng` | switchers | distinct optima |
|---|---|---|---|---|---|
| rotation / informative | 28 | 2 | 31.12 | 0.642 | 6.0 |
| rotation / scrambled | 24 | 6 | 29.73 | 0.313 | 6.0 |
| independent / informative | 18 | **12** | 32.59 | 0.275 | 4.0 |
| independent / scrambled | 19 | **11** | 34.00 | 0.246 | 3.9 |

Structure: **+1.47 (p = 0.52)** informative, **+4.27 (p = 0.10)** scrambled.
Information: null under both structures (p = 0.45, 0.66).

**The hypothesis was wrong.** Making the states genuinely independent does not
raise complexity.

### Three caveats, and why the null survives them

**Extinction is heavy in the independent arms** — 12/30 and 11/30 against 2/30
and 6/30 — so those cells are a selected subset and the +4.27 in particular
should not be leaned on.

**The manipulation was diluted.** Random permutations collide, so the
independent arms presented **4.0 distinct optima on average, not 6**, ranging as
low as 2.

**But dilution does not explain the null.** Within the independent arms,
`corr(distinct optima, ng) = −0.059` across 37 runs. Runs that genuinely
presented six distinct tasks show no complexity gain either, which rules out
"the manipulation was too weak."

### What the organisms did instead

Switchers **fell** in the independent arms, 0.642 → 0.275. When states stopped
sharing structure, organisms tracked the cue *less*, not more, while extinction
roughly doubled. They did not learn six mappings; they stopped trying and paid
the cost.

### Standing of the task null

It now survives rotation vs independent structure, informative vs scrambled
information, three parameter regimes, a five-point heterogeneity gradient, up to
100 seeds per arm, and detection floors of 2–4 nodes. **Genome size here responds
to how much the environment changes, and to population size and gene cost — and
to nothing about the task.** Four distinct task manipulations have now failed,
the last two with adequate power.

## Two follow-ups from adjacent literatures

### Modularity cannot be measured here — the networks are saturated

[Kashtan & Alon (*PNAS* 2005)](https://www.pnas.org/doi/10.1073/pnas.0503610102)
report that modularity does not evolve under a fixed goal, does not evolve under
*randomly* varying goals, and does emerge under goals built from shared subgoals.
Every measurement in this project had been of network *size*; their result says
the quantity that responds to the structure of environmental variation is
**modularity**. The arms here map onto their contrast directly.

Newman's Q on the weight graph, 20 seeds per arm:

| arm | n | Q | edge density | `ng` |
|---|---|---|---|---|
| frozen (fixed goal) | 20 | 0.0587 | 0.9958 | 23.77 |
| scrambled (randomly varying) | 17 | 0.0607 | 0.9980 | 30.45 |
| informative (structured) | 20 | 0.0635 | 0.9947 | 31.15 |

All comparisons null (p = 0.57, 0.80, 0.75).

**This is a failed measurement, not a negative result.** Edge density is **0.9961
across every arm**, and Newman's Q is defined by the *absence* of edges between
communities. At that density there are almost no absent edges, so Q is pinned
near zero by construction and cannot register a difference whatever the
environment does — Q spans only 0.024–0.144 across all 60 runs. Reporting "no
modularity effect" would have been wrong.

This is consistent with something already visible in the long runs, where
`mean_nonzero_w` ≈ `ng²`: these networks saturate. Genome size changes by *node
count* while connectivity stays near complete. Kashtan & Alon's circuits were
sparse by construction; these are not, so testing their claim would need sparsity
to be achievable — a per-edge cost or a hard connectivity cap — which is a model
change rather than another run.

### The fitness null is consistent with theory, not a contradiction of it

The extinction comparison (6.9% informative against 8.8% scrambled, p = 0.46) was
written up earlier as "no detectable fitness advantage to information," which
reads as a finding. Putting numbers on the theory shows it is not one.

[Kussell & Leibler (2005)](https://www.science.org/doi/10.1126/science.1114383)
and [Rivoire & Leibler
(2011)](https://link.springer.com/article/10.1007/s10955-011-0166-2) bound the
growth-rate gain from environmental information by the mutual information between
cue and state:

| quantity | value |
|---|---|
| H(environment) = log₂(6) | 2.585 bits |
| mutual information, informative cue | 2.585 bits (perfect) |
| mutual information, scrambled cue | 0.000 bits |
| **theoretical ceiling on the advantage** | **2.585 bits/generation** |

But that ceiling assumes an organism that reads the cue perfectly *and* converts
it into a perfectly matched phenotype. Organisms here match the correct channel
in **1.55 of 6 states** against a non-responsive baseline of 1/6, so they exploit
**(1.55 − 1)/(6 − 1) = 10.9%** of the available information:

> **realised advantage ≤ 0.282 bits per generation**

Against a detection floor of **5.9 percentage points** in extinction rate, with an
observed difference of 1.9 points. An edge that small, spread over runs of ~5,700
generations with populations in the thousands, is absorbed by demographic
stochasticity long before it appears as differential extinction.

**So the experiment was never powered to see the effect theory predicts.** The
null is consistent with Kussell–Leibler rather than evidence against it, and the
earlier phrasing is withdrawn.

## Connectivity could not evolve — and that invalidates a measure

The modularity attempt above found every arm at ~99.6% edge density. Chasing
that turned up a structural fault rather than a tuning problem.

### Why the networks saturate

Point mutations are **additive**: `W[m, k, i] += d`. Nothing ever sets an
individual weight back to exactly zero, and `p_silence` knocks out an entire
*row* — a node — not an edge. Founders start at 25% density by construction, and
from there density only ever climbs.

Measured against lineage depth, with no per-edge selection at all:

| mean lineage depth | edge density |
|---|---|
| 0 (founders) | 0.25 |
| 304 | 0.777 |
| 569 | 0.828 |
| 1,122 | 0.846 |
| 1,396 | 0.862 |

**Saturation is mutationally inevitable, not a selective outcome.**

### The measure this invalidates

`nzw` counts non-zero weights, so it counts *weights touched since the founder* —
which tracks how long a lineage has been mutating, not how densely it is wired.
Any inference from `nzw` or from "weights per node" is therefore unsound, and one
such inference earlier in this README (that zoned networks were "denser rather
than padded") is retracted where it appears.

It also means **complexity in this model has only ever had one evolvable
dimension: node count.** Connectivity was fixed by construction, so every result
above — the fluctuation effect included — concerns `ng` alone.

### Making sparsity achievable

Two parameters, both off by default and bit-identical when off, because a cost
with no mechanism to produce zeros would only tax everyone equally with no escape:

- `p_edge_silence` — per birth, set one random weight to exactly zero. The only
  route back from non-zero, given additive point mutations.
- `maint_edge` — upkeep per non-zero weight, so sparsity can pay.

Measured effect, grid 48, 30k ticks:

| `p_edge_silence` | `maint_edge` | density | `ng` | pop |
|---|---|---|---|---|
| 0 | 0 | 0.824 | 19.7 | 1872 |
| 0.05 | 0 | 0.629 | 18.0 | 1848 |
| 0.05 | 0.0002 | 0.614 | 17.0 | 1557 |
| 0.2 | 0.0005 | **0.426** | 17.0 | 1370 |

Density now responds, and populations stay healthy. Connectivity becomes a second
evolvable dimension alongside node count, which is what the modularity question
and any topology question needs.

## The completed long run: 8 seeds x 25M ticks x ~630,000 generations

The run finished. Final genome sizes: **13.0, 14.0, 12.0, 12.0, 14.0, 12.0, 12.0,
12.0** — six of eight exactly at `NG_MIN` = 12, with `nzw` ≈ 144 ≈ 12². The
trajectory ran 12.78 at generation 13,875, 12.75 through the middle, 12.63 at
generation 634,225, last-quarter slope **−0.0002 per million ticks**. It reached
the floor by generation ~14,000 and did nothing for the following 620,000.

Digging past the genome size turned up four things.

### Cross-feeding collapses at this timescale

`top_share` for channel 1 is **0.999 throughout**, Shannon diversity 0.006–0.014.
Essentially the whole population specialises on the primary resource, and the
waste channels drain with it — res5 falls 20,903 → 1,082 and res6 3,495 → 10.

**This qualifies a headline result above.** The cross-feeding figure reported
earlier (stable at 0.26–0.30) comes from runs of 150k ticks. Over 630,000
generations the ecotype structure does not merely weaken, it disappears. Whatever
maintains cross-feeding on short timescales does not hold on long ones.

### The population swings ~70×, but the cycle is imposed, not emergent

Population ranges from ~500 to ~37,000 with CV 0.48–0.54, and the pattern is
near-identical in all eight seeds. Autocorrelation gives a period of **exactly
15,000 ticks, standard deviation 0 across seeds, r ≈ 0.995**.

That exactness is the giveaway: `season_period = 15_000` is a parameter, and
inflow is driven by `1 + 0.6·sin(2π·gen/15000)`. **This is the imposed seasonal
forcing, not a discovered limit cycle** — it was nearly written up as one.

What is worth noting is the **amplification**: a ±60% swing in inflow produces a
~70-fold swing in population. The demographic response is strongly nonlinear.

### Motility tracks the season — genetically, not plastically

Motility swings between 0.002 and 0.994, in near-perfect antiphase with
population (**r = −0.938 at zero lag**). Low population, high motility, larger
size, resource-rich; high population, motility off, smaller size, resources drawn
down.

That looks like evolved plasticity, so it was tested directly: take the standing
population and develop the *same genomes* across resource levels spanning
0.05–4.0 per cell. Motility comes out **0.980, 0.994, 0.997, 0.998, 0.987,
0.941** — essentially flat.

**The swing is genetic turnover between seasons, not plasticity.** At ~27 ticks
per generation there are ~550 generations per season, ample for selection to
cycle the population. This is the same conclusion the grain analysis reached from
the other direction: when the environment changes slowly relative to generation
time, genetic tracking does the work and plasticity is never needed.

### Genome size is decoupled from all of it

`corr(pop, ng)` is **−0.015, −0.018, −0.030** in the first three seeds. The
demographic cycle, the motility turnover and the resource swings leave genome
size untouched. Repair also sits at ~0.001 — effectively switched off — while
`mean_load` stays flat at 0.25–0.34, so Muller's ratchet is not running despite
repair being abandoned.

## With connectivity evolvable: fluctuation acts on edges more than nodes

Every result above concerns `ng` alone, because connectivity was saturated by
construction. With `p_edge_silence` and `maint_edge` making sparsity achievable,
the three arms were re-run — 20 seeds each, 150k ticks:

| arm | n | `ng` | edges | density | Q | pop |
|---|---|---|---|---|---|---|
| frozen (no fluctuation) | 20 | 18.01 | **201.3** | 0.827 | 0.101 | 3238 |
| fluctuating, 0 bits | 17 | 23.36 | **415.2** | 0.923 | 0.085 | 2284 |
| fluctuating, 2.585 bits | 18 | 24.26 | **436.3** | 0.909 | 0.090 | 2865 |

Effect of fluctuation (frozen → fluctuating):

| measure | frozen → informative | p | frozen → scrambled | p |
|---|---|---|---|---|
| `ng` | **+6.25** | <0.0001 | +5.35 | <0.0001 |
| **edges** | **+235.0** | **<0.0001** | +213.9 | <0.0001 |
| density | +0.083 | <0.0001 | +0.096 | <0.0001 |
| Q | −0.011 | 0.148 | **−0.017** | **0.016** |

**Fluctuation acts on both, and far more strongly on connectivity.** Node count
rises ~31%; edge count **more than doubles**, +108%. The effect reported earlier
as "fluctuation raises complexity by ~7 nodes" was measuring the only dimension
that could move; the larger response was in wiring all along.

Information remains irrelevant on every measure — scrambled against informative
gives p = 0.47 (`ng`), 0.68 (edges), 0.079 (density), 0.53 (Q). That conclusion
survives the richer measurement.

### Modularity falls under fluctuation

Q drops from 0.101 frozen to 0.085 scrambled (**p = 0.016**), with the
informative arm between. That is the *opposite* of what Kashtan & Alon report for
structured variation, and it follows directly from the density result: fluctuation
drives networks denser, and denser networks are less modular.

Two caveats. Density is still high in absolute terms (0.83–0.92), because the
sparsity settings had to be mild — the wave arm goes extinct at
`maint_edge` = 0.0005, so 0.00005 was the largest survivable cost. And these Q
values, while now responsive, remain small. A regime with genuinely sparse
networks would test Kashtan & Alon properly; this one only shows the direction.

## Development does not converge, and `niter` was a hidden parameter

`develop!` iterates the network and the code described the result as *"an
approximate fixed point"*. Prompted by the criticality literature — Kauffman's
ordered/critical/chaotic classification turns on exactly this, and these networks
run at 0.83–0.99 edge density — that assumption was tested. It is false.

| measurement | value |
|---|---|
| mean \|activation\| after settling | **0.846** (1.0 = fully saturated) |
| nodes flipping sign **every** iteration | **92.1%** |
| mean \|weight\| | 0.183 |
| step size, iteration 24 ÷ iteration 8 | **1.03 — no decay at all** |

Activations saturate and then the network oscillates almost globally: nearly
every node flips sign on every iteration, indefinitely. Step size *grows* from
0.365 at iteration 1 to ~0.9 and stays there. `develop!` stops at iteration 8 and
reads off whatever phase of that oscillation it lands on.

### The consequence

Changing the iteration count from 8 to 9 — **or to 10** — flips the top uptake
channel in **100% of organisms**. Mean phenotype distance from the niter-8 result:
1.78 at 9, 2.30 at 10, 0.97 at 11, 0.23 at 12, 0.39 at 20, 1.82 at 21. Not even a
clean parity effect; the trajectory is complex enough that any change reshuffles
every phenotype.

So the iteration count is a parameter with **total effect on the phenotype**, and
it had been hard-coded as a keyword default and never varied. It is now exposed
as `dev_iters` and swept like anything else.

### What this does and does not invalidate

It does **not** invalidate the comparative results. `niter` = 8 throughout, the
genotype→phenotype map is deterministic and heritable, and it is robust to small
weight perturbations — sensitivity ≈ 0, because saturation absorbs them.
Evolution optimises "the state at iteration 8", which is a well-defined trait.

It does mean three things. The description of development in this codebase was
**wrong** and is corrected. The biological metaphor the design rests on —
development settling to a stable expression state — is **not what happens**. And
a parameter with complete control over the phenotype went unexamined through
every experiment above, which is why the headline result is being re-run across
`dev_iters` ∈ {4, 8, 12, 20} before anything is written up.

This is the same regime as the saturation finding: 99.6% edge density and 92%
per-iteration sign-flipping are one phenomenon. In Kauffman's terms these
networks sit firmly on the chaotic side rather than at criticality, which is
where dense networks are expected to be.

## Code review: one real defect, found and fixed

A systematic review of the paths the claims rest on. Most held; one did not.

### The scrambled cue was not at zero bits in the runs that used it

The scramble was tuned and verified **with a spatial gradient**, where it
measured 0.0000 bits. But `cond_gradient = false` drops the spatial term, and the
remaining integer temporal multiplier (7 turns per cycle) is commensurate with
the cycle, so the phases lock:

| multipliers | no gradient | with gradient |
|---|---|---|
| 7.0 / 5.0 (as used) | **0.0666 bits** | 0.0000 bits |
| 6.1803399 / 5.0 (fixed) | **0.0004 bits** | 0.0014 bits |

**Every experiment after the grain work used `cond_gradient = false`**, so their
scrambled arms carried ~0.067 bits rather than the 0.000 claimed here and in the
commit messages.

**Does it change any conclusion? No.** The comparison remains 2.585 bits against
0.067 — a **39× difference** — so "informative versus effectively uninformative"
still holds, and any residual signal in the scrambled arm makes the information
null *conservative* rather than inflated. But the figure quoted was wrong, and
statements of "0.000 bits" for non-gradient runs are corrected to 0.067.

The multiplier is now **irrational** (5 + φ), which cannot lock in either mode,
with a test covering gradient and non-gradient. Note this changes the scrambled
cue stream, so scrambled runs from before this fix are not bit-reproducible under
the current code.

### What the review checked and cleared

- **Determinism under threading.** Identical population, mean `ng` to six
  decimals, total `nzw` and phenotype sum at **1, 4 and 8 threads**, on the full
  feature set — sparsity, reversible plasticity, abiotic cue and all.
- **Sparsity is real disconnection, not an accounting artefact.** `maint_edge`
  charges for `|w| > 1e-4`, so weights could in principle hide below the
  threshold. Measured: **0.0%** of weights fall in that band, 28.7% are exactly
  zero, and density by exact-zero and by threshold agree to four decimals
  (0.7133 both ways).
- **The frozen arm is genuinely frozen** — one state visited across 150,000
  ticks, cue constant to within 3×10⁻⁴.
- **Inflow follows the best-paying channel 6/6** under both rotation and
  independent permutations.
- **`redevelop!` keeps `chan_cost` consistent** with the phenotype it derives
  from — zero stale entries across a run.
- **Maintenance accounting** applies `maint_edge` to the cached `nzw`, which is
  refreshed at birth and stays valid because the genome does not change during
  life.

## Meta-analysis: pooling every comparison in this project

Before any write-up, every comparison run here was pooled rather than reported
cell by cell. Reporting sixteen experiments individually invites picking the
flattering ones; pooling shows whether the pattern holds.

### The fluctuation effect — 7 comparisons, all positive

| experiment | `ng` frozen → fluctuating | p | edges | p |
|---|---|---|---|---|
| edges experiment | 18.01 → 23.82 (+5.81) | <0.0001 | +112% | <0.0001 |
| sparse level A | 24.04 → 34.20 (+10.15) | 0.010 | +98% | 0.022 |
| sparse level B | 18.60 → 25.06 (+6.47) | 0.0001 | +122% | 0.0001 |
| sparse level C | 16.04 → 18.68 (+2.64) | 0.005 | +63% | 0.002 |
| modularity experiment | 23.77 → 30.83 (+7.06) | 0.0003 | — | — |
| `niter` 4 (sparse) | 18.71 → 24.99 (+6.28) | <0.0001 | +116% | <0.0001 |
| `niter` 8 (sparse) | 18.01 → 23.82 (+5.81) | <0.0001 | +112% | <0.0001 |

**Every comparison positive, every one significant.** `ng` mean +6.32 (range
+2.64 to +10.15); edges mean **+104%** (range +63% to +122%). Edges respond
roughly **3.3× more than nodes**, across sparsity levels, `niter` values and
separate experiments.

### The information null — 15 comparisons, centred on zero

Pooling every informative-vs-scrambled comparison on `ng`:

- **9 positive, 6 negative**
- mean difference **−0.171**, SE 0.562 → **t = −0.30** against zero
- sign test: 9/15 positive, **p = 0.61**
- significant at p < 0.05: **2/15**, against a chance expectation of 0.8

And the two significant results point in **opposite directions** — regime C at
−3.19 (p = 0.037) and `niter` 4 at +2.98 (p = 0.009). That is what noise looks
like, not a suppressed effect. It also defuses the `niter` = 4 result flagged
above: one of two hits in fifteen, with its counterpart going the other way.

### Code-level claims, verified together

All nine "off by default, bit-identical when off" claims checked against a common
baseline on population, phenotypes **and** weight matrices: `uptake_half`,
`cond_states`, `redevelop_every`, `independent_states`, `p_edge_silence`,
`maint_edge`, `dev_iters`, `cond_scramble`, `rotate_yield` — **9/9 pass**, with
every default matching what this README states, and a negative control
confirming the check can detect divergence.

The first version of that check silently passed regardless of failures: the
accumulator hit Julia's soft-scope rule, became a new local, and the script
printed the untouched global. It is the trap already documented in this
README's development notes, written again anyway — which is the argument for the
negative control that now accompanies it.

## The cost of complexity: fluctuation is lethal when wiring is expensive

The sparse regime was re-run at 60 seeds per arm on a second machine, because the
first pass gave only 3 and 2 survivors — the sample size that has produced three
false signals in this work.

| arm | n | extinct | `ng` | edges | density | Q | pop |
|---|---|---|---|---|---|---|---|
| frozen | 60 | **0/60** | 15.99 | 98.9 | 0.604 | 0.146 | 3091 |
| informative | 23 | **37/60** | 17.42 | **153.4** | 0.735 | 0.138 | 2296 |
| scrambled | 27 | **33/60** | 16.97 | **147.6** | 0.729 | 0.140 | 2225 |

### The strongest single result in this project

**Frozen dies in 0 of 60 runs. Fluctuating dies in 37 and 33 of 60. Fisher exact
p = 3.6×10⁻¹³.**

The arms differ *only* in `cond_period` (10⁹ against 600) — same grid, same
sparsity, same inflow fraction, and the frozen arm still routes 20% of inflow by
the rotation, merely to a fixed state. There is no confound; the comparison is
exactly what it claims.

And it is information-independent: 37/60 against 33/60 is indistinguishable
whether or not the cue predicts anything. **What kills these populations is
change itself, not unpredictable change.**

### It also explains the edge result mechanistically

Organisms meet fluctuation by **adding wiring** — edges +55% and +49%
(p < 0.0001) against nodes +1.4 and +1.0 (p = 0.001, 0.008). Put a price on
wiring and they cannot mount that response, so they die. The complexity increase
is not decorative: **it is what surviving a changing environment costs**, and
under a per-edge tax most populations cannot pay it.

Information remains null on every measure — `ng` p = 0.37, edges p = 0.62,
density p = 0.76, Q p = 0.92.

### Audit of these numbers

**Survivorship was the obvious threat and it does not explain the edge effect.**
Frozen is unselected at 0/60 while the fluctuating survivors are a filtered
sample, so the comparison could in principle be an artefact of selection. But the
pattern runs the wrong way:

| level | density | info extinct | edge effect |
|---|---|---|---|
| 0.05 / 0 | 0.987 | 1/10 | +62% |
| 0.2 / 5e-5 | 0.823 | 1/10 | **+122%** |
| 0.3 / 2e-4 | 0.618 | 7/10 | +66% |
| 0.3 / 2e-4, n = 60 | 0.604 | **37/60** | **+55%** |

Extinction rises from 10% to 62% while the effect *falls* from +122% to +55%. If
survivorship were manufacturing it, the effect would be largest where selection
is hardest; it is smallest there.

The qualification stands, though: **the +55% is measured on survivors** and is
not directly comparable to the +122% at lower extinction. Treat the existence of
the effect as established and its magnitude at high extinction as approximate.

**The Q null needs its detection floor attached.** Q is now measurable — 0.146 at
density 0.60, against 0.06 at saturation — and shows no fluctuation effect
(p = 0.31, 0.48). But with n = 23 in the informative arm the minimum detectable
difference is **0.022, roughly 15% of the frozen mean**. So the correct statement
is "no effect larger than 15%", not "no effect". Kashtan & Alon's reported effects
are larger than that, so it remains informative.

### What this closes

The regime sparse enough for modularity to be meaningful is the regime where
fluctuating environments kill 55–62% of populations. Testing Kashtan & Alon
properly here is therefore not a matter of tuning — it is blocked by a measured
trade-off, backed by 180 runs.

## Aevol: the same conclusion by a different mechanism

A late sweep of genome-architecture literature — searched with different terms
than the GRN work above — found the Aevol platform (Knibbe, Beslon and
colleagues), which has been asking this question with explicit sequence-level
genomes since 2007.

**Their conclusion is this project's conclusion.** [Luiselli et al.
2024](https://academic.oup.com/gbe/article/16/12/evae250/7905804) find larger
populations drive genome size down (coding density ~80% at 4× population, ~88% at
16×) through a trade-off between direct selection for adaptation and **indirect
selection for replicative robustness**. [Knibbe et al.
2007](https://academic.oup.com/mbe/article/24/10/2344/1074262) show that
robustness to large duplications and deletions *bounds* genome size, with Drake's
rule emerging spontaneously and log-log fits at r² = 0.99.

So "genome size is set by mutational and population-genetic pressure rather than
functional demand" is established, quantified, and older than this work. The long
run here reproduces their streamlining result.

### But their mechanism cannot operate in this model

Knibbe et al. are explicit: the mechanism *"would collapse without rearrangements
— point mutations alone wouldn't couple genome size to mutational variability."*
Their organisms tune noncoding DNA to hold roughly **one neutral offspring per
generation**, and that works only because larger genomes suffer proportionally
more deleterious rearrangements.

This model has **no chromosomal rearrangements**. Node duplication copies a node
with its weights; it creates no size-proportional mutational exposure. So the
Aevol explanation is unavailable here, and the two models reach the same
conclusion by different routes — convergent evidence rather than duplication.

It also explains a difference: the floor here is **hard** (`ng` = 12 exactly, six
of eight seeds) while theirs is a stabilised distribution, because without the
equilibrium-producing pressure selection simply runs to the minimum.

### Their compound parameter does not govern here — tested

Luiselli et al. report that **N×μ** governs coding density: constant N×μ gives
similar architecture. Tested directly, scaling μ with inverse grid area so
population × mutation rate is held constant:

| cell | grid | μ | pop | `ng` | N×μ |
|---|---|---|---|---|---|
| N/4, μ×4 | 32 | 0.004 | 1,257 | 20.72 | 5.03 |
| baseline | 64 | 0.001 | 6,778 | 19.98 | 6.78 |
| N×4, μ/4 | 128 | 0.00025 | 30,616 | 22.87 | 7.65 |
| **N, μ×4** | 64 | 0.004 | 4,651 | 19.68 | **18.60** |

Constant-N×μ cells should match and the 4× cell should differ. **Both fail**:
`large` differs from baseline despite matched N×μ (+2.89, p = 0.022), while
`refhigh` at nearly 4× the product is indistinguishable (−0.30, p = 0.77).

**N×μ does not govern architecture in this model**, exactly as Knibbe's own
caveat predicts for a model without rearrangements. Their framework is
rearrangement-specific rather than general.

## `dev_iters` sensitivity, corrected: the result barely moves

The first version of this test set neither sparsity parameter, so networks
saturated and "edges" merely tracked `ng²` — it probed the saturated regime
rather than the sparse one where connectivity is a real measurement. That run is
retained below as a methods note; this is the corrected one, sparsity on, 20
seeds per arm.

`dev_iters` controls the phenotype completely — changing it 8 → 9, or 8 → 10,
flips the top uptake channel in **100%** of organisms — so it is the sharpest
available test of whether these results are an artefact of one setting.

| `dev_iters` | Δ`ng` | p | Δ density | p | excess over `ng²` |
|---|---|---|---|---|---|
| 4 | +6.28 | <0.0001 | **+0.144** | <0.0001 | +38% |
| 8 | +5.81 | <0.0001 | **+0.138** | <0.0001 | +37% |
| 12 | +5.71 | <0.0001 | **+0.118** | <0.0001 | +31% |
| 20 | +5.36 | <0.0001 | **+0.120** | <0.0001 | +30% |

**The effect barely moves.** `ng` spans +5.4 to +6.3, density +0.118 to +0.144,
excess over the `ng²` baseline +30% to +38%, every cell at p < 0.0001. Across a
parameter that reshuffles every phenotype in the population, the result is
stable.

That is also the answer to the obvious objection raised by the non-convergence
finding above — if development is a snapshot of an oscillation, why should the
phenotype mean anything? Because the result does not depend on where in the
oscillation it is sampled.

### And the information null holds here too

One cell had suggested otherwise. It does not replicate:

| `dev_iters` | Δ`ng` informative − scrambled | p |
|---|---|---|
| 4 | +2.98 | **0.009** |
| 8 | +0.90 | 0.47 |
| 12 | +0.54 | 0.73 |
| 20 | +0.70 | 0.48 |

On density the pattern repeats and the single hit is *negative* (`dev_iters` = 8,
p = 0.021). Two isolated hits in eight comparisons pointing opposite ways is
noise, and it matches what the pooled meta-analysis already indicated.

### The methods note worth keeping

The flawed first version reported the fluctuation effect on "edges" swinging
+23% to +112% across `dev_iters`, and that instability was entirely an artefact
of measuring `ng²` under another name. **A proxy substituted for the quantity of
interest produced an apparent result about robustness that was not about
robustness at all.**

## Novelty audit: what is actually new here## Novelty audit: what is actually new here

Run *after* writing the result up, which was the wrong order and is the fourth
time in this project a novelty claim has had to be cut back.

### The main effect is published, three times over

[Bentkowski, van Oosterhout & Mock 2015,
*GBE*](https://pmc.ncbi.nlm.nih.gov/articles/PMC4558865/) is essentially this
result: greater environmental variability produces larger genomes, expanding from
2 genes in a constant environment to 15–16 under turbulence, with both
transcription-factor count and metabolic-network size rising. Their proposed
mechanism is **relaxed selection on metabolic efficiency** rather than selection
for plasticity — which is what the data here show too, since the scrambled arm
rises just as much as the informative one.

[Fluctuating selection and intrinsic constraints
(2010)](https://journals.plos.org/ploscompbiol/article?id=10.1371%2Fjournal.pcbi.1000873)
reports the same direction, as does [*Size
Matters* (2011)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3171436/), whose
heterogeneity manipulation the H sweep above reproduces.

**So "fluctuation raises complexity" is not a finding of this work.** It is a
reproduction, in a different model, of an established result.

### The control has a name, and a known critique

The scrambled-cue design is a **[yoked
control](https://en.wikipedia.org/wiki/Yoked_control_design)** — matched
stimuli, different contingency — standard in operant conditioning since the
1950s. It is not a new design, only a new application of an old one.

It carries **Church's 1964 critique**: in the classic form the yoked subject
receives a *copy* of its partner's stimulus schedule, so the two are not
statistically independent, and that correlated-error structure can produce
apparent group differences from nothing but between-subject noise.

**That critique does not apply to the design used here**, and the reason should
be stated up front rather than left for a reviewer to raise. The scrambled arm
does not replay a paired informative run: its cue is generated independently
(7 temporal turns against 5 spatial, chosen to drive mutual information to
0.0000 bits), every run is an independent replicate, and arms are compared
between independent samples rather than as matched pairs. There is no
master–yoke coupling to induce correlated errors.

### What does appear to be new

None of those studies **separates the variation from its information content**.
Bentkowski et al. drive the environment with a bounded random walk and never
manipulate predictability; the others likewise have no zero-information control.

Two controls here do separate them, and both are what any claim rests on:

- the **scrambled cue** — identical world, identical rotation, identical
  ecology, mutual information with the environmental state driven from 2.585
  bits to 0.0000;
- the **frozen arm** — `cond_period` = 10⁹, so the rotation machinery still
  routes inflow exactly as in the fluctuating arms while the state never changes.

Together they show the complexity rise **does not require the fluctuation to
carry usable information**. That is a refinement of a published effect, not a
discovery.

### A dissociation that does NOT hold, and was nearly claimed

[Kussell & Leibler 2005](https://www.science.org/doi/10.1126/science.1114383)
established that information about a fluctuating environment has quantifiable
fitness value — the KL divergence is the penalty paid for poor information. That
made an attractive story available: information buys fitness and behaviour but
not genome complexity.

Testing it on this data kills it. Across **640 runs**, extinction was **6.9% with
an informative cue and 8.8% with a scrambled one — Fisher exact p = 0.46**. There
is no detectable fitness advantage to information here.

So the defensible statement is narrower: an informative cue reliably changes
**behaviour** (switchers, p < 10⁻⁵ in three independent regimes) and changes
neither genome size nor survival.

### A population-size result this may conflict with

[LaBar & Adami 2016, *PLOS Computational
Biology*](https://journals.plos.org/ploscompbiol/article?id=10.1371%2Fjournal.pcbi.1005066)
report a **U-shape**: both small *and* large populations evolve larger genomes,
by different routes — small populations fix slightly deleterious insertions
(drift), large populations fix rare beneficial ones — concluding that *purifying
selection is not powerful enough to prevent the evolution of complexity in large
populations*.

The data here are **monotonic, not U-shaped**: the largest world tested sits at
the floor (`ng` = 12.0, grid 128, populations 6,800–34,400, flat for 285,000
generations) and raising grid 64 → 192 *reduced* `ng` from 27.8 to 23.1. That is
the opposite of their large-population arm.

This is recorded as an open discrepancy rather than resolved. Their populations,
genome architecture and mutational regime all differ substantially, and nothing
here tests their mechanism directly. But it is the one place where this work
appears to disagree with a published result, and it should not be glossed.

### The empirical literature agrees with the effect too

Beyond simulation: terrestrial environments show two to three orders of magnitude
greater variation than marine ones, and terrestrial cyanobacteria carry genomes
at the larger end (6.0–8.0 Mb) *enriched in regulatory, transport and motility
genes*. Across bacteria generally, [regulatory genes are disproportionately
enriched in larger
genomes](https://www.pnas.org/doi/10.1073/pnas.0308653100), and [genome size
tracks ecosystem type and trophic
strategy](https://www.frontiersin.org/journals/microbiology/articles/10.3389/fmicb.2021.761869/full),
with streamlining in nutrient-limited settings.

So the fluctuation effect holds in real organisms as well as in models — which
strengthens the reproduction and further weakens any novelty claim attached to it.

### Also already established, and not claimed here

- **The architecture.** [Sakamoto & Innan
  2024](https://academic.oup.com/g3journal/article/14/9/jkae144/7710408) evolve
  discrete plasticity via a dedicated sensor gene feeding a GRN — the `P_SENSE`
  design — with gene number fixed. [Jorritsma & van den Berg, *Evolution*
  2026](https://academic.oup.com/evolut/article-abstract/80/1/230/8297131)
  occupy the same parameter space.
- **"Rewiring rather than adding"**: [rewired bacterial
  networks](https://www.nature.com/articles/nature06847), [Ndt80
  regulon](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5391208/), [echinoid
  GRNs](https://www.pnas.org/doi/10.1073/pnas.1509845112).
- **"Complexity is not directly selected"**: van Nimwegen's [quadratic TF
  scaling](https://academic.oup.com/nar/article/45/13/7615/3865434), explained by
  addition/deletion rates rather than environmental demand; [GBE
  2026](https://academic.oup.com/gbe/article/18/8/evag196/8753180).

### Standing

A published effect reproduced, with one new control showing it does not depend on
the information content of the variation, and a fitness dissociation explicitly
checked and rejected. That is a short note refining Bentkowski et al., not a
headline result.

## Predictions, recorded before the runs that tested them

Every substantive prediction in this work was written into the repository
*before* the experiment that tested it, so the git history is a usable
pre-registration record rather than a reconstruction. The claims and the commits
that predate their tests:

| prediction, made in advance | outcome |
|---|---|
| `ng` rises with task complexity → the plateau is task-limited; `ng` flat → drift–selection sets it regardless of task | **flat** |
| rotating yields alone will make the task conditional | **wrong** — channel 1's inflow monopoly kept a fixed diet good enough |
| an abiotic cue will be read because it cannot be depleted | **wrong** — legible but unread, for the same inflow reason |
| the conditional genome is adaptive but unreachable by mutation | **wrong** — it is adaptive-looking but *selected against* |
| reversible plasticity will rescue the regime that irreversible development cannot survive | **correct** — 0/3 → 3/3 viable |
| conditional response, once it evolves, will raise genome size | **wrong** — response evolves, complexity does not follow |

Five of six predictions were wrong, which is the honest summary of the work and
the reason the controls matter more than the hypotheses.

## The long run: what actually happens after 120,000 generations

8 seeds, 25M ticks each, grid 128, default parameters — which means the gene
cost is **on** (`maint_gene = 0.0015`) and populations are large (6,800–34,400).
At ~4.6–5.0M ticks and ~120,000 generations:

| seed | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|---|---|
| mean `ng` | 13.0 | 15.0 | 12.0 | 12.0 | 14.0 | 12.0 | 12.0 | 12.0 |

`NG_MIN` is 12. **Five of eight seeds sit exactly at the minimum legal genome
size**, with `mean_nonzero_w` = 144.1 ≈ 12² — every available weight nonzero.
The networks minimised their node count and then saturated connectivity inside
the floor.

| tick | generations | mean `ng` | population |
|---|---|---|---|
| 92,875 | 2,958 | 13.94 | 29,679 |
| 464,500 | 12,937 | 12.77 | 18,178 |
| 1,858,000 | 49,816 | 12.76 | 8,841 |
| 3,251,625 | 85,077 | 12.76 | 6,290 |
| 4,645,125 | 120,355 | 12.75 | 7,044 |

It fell to the floor by generation ~13,000 and has not moved in the 107,000
generations since. Last-quarter slope: **−0.0005 `ng` per million ticks** over
1.2M ticks — a dead-flat line.

So the answer to "what happens if you run it for millions of generations" is
**streamlining to the minimum, then nothing**. Note this is not the ~21 plateau
reported above: that figure belongs to the complexity-permissive regime (small
population, no gene cost). Under default settings, with a gene cost and a large
population, complexity does not plateau at a moderate value — it is driven to
the absolute floor and pinned there. Large populations with efficient purifying
selection streamline, exactly as the mutational hazard hypothesis predicts, now
confirmed at roughly ten times the scale of anything else tested here.

### Putting it together

Every result in this project resolves to the same two variables:

| condition | outcome |
|---|---|
| large population + gene cost | floor: `ng` 12, fully connected |
| small population, no gene cost | plateau ~21 |
| harsher ecology → smaller population | ~25–26 |
| task demanding conditional response | **untestable — it never pays in any survivable regime** |

Genome size is set by population size and by whether genes are taxed — that much
is measured, and Lynch's mutational hazard hypothesis accounts for all of it.
Whether task difficulty *would* enter is unknown, because no arm of any
experiment here ever evolved the conditional behaviour that would test it. The
scrambled-cue control remains the useful piece: holding the world fixed and
destroying only the information is the right shape for this question, and it is
what exposed the fact that the cue was being ignored.
