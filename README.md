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
node rise too, 21.0 → 25.5, so the networks are denser rather than padded.

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

**This is the cleanest negative in the project, because it is the only one whose
control holds ecology constant.** A task that genuinely requires a conditional
response does not move the complexity plateau. What moves it is how small and
how hard-pressed the population is.

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
| task genuinely demanding conditional response | **no additional effect** |

Genome size is set by population size and by whether genes are taxed. Task
difficulty does not enter. That is Lynch's mutational hazard hypothesis, and the
scrambled-cue control is what makes it a measurement here rather than an
assumption.
