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

A *tick* is one simulation step. A *generation* is one round of descent. At
default settings one generation takes roughly **15 ticks**, so a million-tick
run is about 65,000 generations, not a million.

Every run logs `mean_depth` — the mean lineage depth of living organisms —
which is the real generation count. Never infer it from the tick number. The
runner prints the measured ticks-per-generation ratio when it finishes.

## The model

### Genome

A gene-regulatory network: a 20×20 weight matrix plus per-node bias. At birth
the network is run to an approximate fixed point, and the resulting expression
vector is read off as a **17-trait phenotype** which is then cached for life.

Development happens **once, at birth** — not every tick. Re-expressing every
tick measured 15× slower (44.3 vs 2.9 ms/generation) for no modelling benefit,
and "develop, then live as an adult" is the more faithful model anyway.

Local resource concentrations feed the first six network nodes, so the same
genome develops differently in different conditions (phenotypic plasticity).

Mutations are point perturbations of weights, plus gene **silencing** and gene
**duplication** — the latter is how regulatory networks gain new functions.

### Phenotype

| Slots | Trait | Notes |
|---|---|---|
| 1–6 | uptake affinity per channel | an **allocation budget**, sums to 1 |
| 7–12 | excretion profile | which wastes this organism emits, sums to 1 |
| 13 | body size | scales uptake *and* upkeep (charged on size²) |
| 14 | motility | movement probability, costed |
| 15 | division threshold | energy needed to divide |
| 16 | repair fidelity | suppresses mutation rate, costed |
| 17 | stress tolerance | survival odds in heat shocks, costed |

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
live in it. Waste channels yield progressively less energy, so eating them is a
living but not a good one.

Three pressures keep selection from settling:

- **Cross-feeding** — the fitness landscape is built by the population itself.
- **Seasons and shocks** — sinusoidal inflow, plus `cull` (indiscriminate
  bottleneck), `heat` (selects on stress tolerance) and `famine` (inflow
  collapses for a while) events.
- **Mutation load** — most mutations are deleterious, some beneficial, so load
  can fall as well as rise; Muller's ratchet is a dynamic, not a countdown.

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

Measured on an Apple M5 (10 threads), grid 256, ~30k organisms:

| | ticks/s | 1M ticks | ≈ generations |
|---|---|---|---|
| default settings | ~200 | ~1.4 h | ~65,000 |

Reaching a million *generations* needs roughly 15M ticks. Run that headless
with `scripts/run.jl` — the core has no plotting dependency, so it runs fine on
a machine with no display. `stats.bin` is flushed as it goes, so a long run can
be charted (or stopped) at any point without losing what it has done.

Memory scales with dish area, since capacity is derived from it:

| grid | capacity | organism arrays | approx ticks/s |
|---|---|---|---|
| 128 | 65,536 | 29 MB | ~250 |
| 256 | 262,144 | 465 MB | ~50 |
| 512 | 1,048,576 | 1.8 GB | ~13 |

For scale, the same loop in pure Python measured **440 ns per organism-tick
versus 8.3 ns** in Julia — and that comparison excluded the diffusion step,
which Julia also had to do.

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
- **Same seed reproduces exactly.** Without this no result is ever checkable.
- **The hot loop does not allocate.**

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

Charts are the debugger here. Run a few tens of thousands of ticks and look at
them before committing to a long run.
