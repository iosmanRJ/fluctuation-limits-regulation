# =============================================================================
# Every knob in one place. Mutable so the live window can tune during a run.
# =============================================================================

Base.@kwdef mutable struct Params
    # --- world ---------------------------------------------------------------
    grid::Int            = 256
    # Capacity is an ARRAY bound, not an ecological one. It must stay well
    # above the resource-determined carrying capacity: if the dish fills, births
    # get dropped, reproductive success stops mattering, and selection collapses
    # into "minimise every cost" — every trait then walks to zero.
    #
    # 0 means "scale with the dish" (grid^2 * capacity_per_cell). A fixed default
    # silently became the binding constraint at grid 256, so it is now derived.
    capacity::Int          = 0
    capacity_per_cell::Int = 4          # headroom; max_per_cell is the real limit
    max_per_cell::Int32  = Int32(8)    # hard local packing limit
    # Founders, 0 = scale with dish area. A fixed count is a density trap: 200
    # founders is 0.012/cell at grid 128 but 0.0008/cell at grid 512, 16x
    # sparser, which makes runs at different grid sizes incomparable.
    #
    # It was ALSO hypothesised to explain the early extinctions seen at grid
    # 256. It does not: scaling founders moved extinctions only 3/8 -> 2/8 over
    # 8 seeds, which is nothing at that sample size, and a scaled run still died
    # at tick 10,986. The actual cause turned out to be famine (8 of 8 deaths
    # followed one within 24-211 ticks), addressed by dormancy. Keep the
    # scaling because comparability across grid sizes is worth having on its
    # own; do not credit it with fixing extinction.
    n_founders::Int      = 0
    seed::Int            = 42

    # --- resources -----------------------------------------------------------
    # Channel 1 is the only abiotic input. Channels 2..NRES exist ONLY as
    # metabolic waste, so every one of them is a niche that has to be created
    # by some other organism before anything can live in it.
    inflow::Float32      = 0.250f0
    yields::Vector{Float32} = Float32[1.0, 0.60, 0.45, 0.35, 0.28, 0.22]  # thermodynamic decline

    # --- how hard is the task? ------------------------------------------------
    # The dish is split into `yield_zones` vertical bands, and each band ROTATES
    # the yield vector, so which channel pays best depends on where you are.
    # With 1 zone the task is the original one: a single fixed optimum, which a
    # minimal network can encode. With several zones an organism must sense its
    # context and respond differently — a conditional map, which is what hidden
    # units are for.
    #
    # This is the untested question. Complexity is known to plateau; Walsh 2026
    # attributes that to a shrinking solution space and Lynch to population
    # size, but nobody appears to have asked whether a task that actually
    # REQUIRES complexity moves the plateau.
    # NOTE on the first attempt: rotating only the YIELDS does not make the task
    # harder in the intended way. Channel 1 is the sole channel with external
    # inflow, so it is far more abundant everywhere; a yield rotation spanning
    # 1.0 to 0.22 cannot overcome that. Organisms kept eating channel 1 in every
    # zone and simply paid for it (population 958 -> 505 with no change in
    # genome size). Rotating the INFLOW channel too is what actually moves the
    # optimum between zones.
    #
    # The two halves of the manipulation are separable so the effect can be
    # decomposed: `rotate_inflow` moves WHICH CHANNEL IS ABUNDANT between zones,
    # `rotate_yield` moves WHICH CHANNEL PAYS BEST. Turning either off while
    # holding `yield_zones` fixed keeps the spatial structure — and therefore
    # the population size and the standing resource stock — and removes only
    # one source of context-dependence.
    yield_zones::Int        = 1
    rotate_inflow::Bool     = true     # zone r also receives its inflow on the rotated channel
    rotate_yield::Bool      = true     # zone r also rotates which channel pays best

    # --- the abiotic condition -----------------------------------------------
    # The zoned task failed for a structural reason: an organism had to infer
    # its context from resource concentrations, but a consumer FLATTENS exactly
    # the field it would need to read. Consumption plus the excretion ladder put
    # the realised peak downstream of whatever was actually supplied, so the
    # zones were never legible and organisms solved the dish with a fixed
    # strategy. No amount of diffusion, grid or decay tuning fixed that.
    #
    # An abiotic condition cannot be eaten, so it stays readable however hard
    # the population exploits the dish. It cycles through `cond_states` discrete
    # states, each rotating which channel pays best exactly as the zones did,
    # and it is fed to the genome through the sensor node P_SENSE. With
    # `cond_states = 1` nothing rotates and the signal is constant: the control.
    #
    # It varies in TIME, not just space, so there is no region to sit in and no
    # lineage that can opt out. Development happens once at birth, so an
    # organism is committed to the state it was born into; the GENOME therefore
    # has to encode a map from cue to phenotype that works in every state, which
    # is the conditional task hidden units are supposed to be for. A genome that
    # ignores the cue is right 1/cond_states of the time.
    cond_states::Int       = 1
    cond_period::Int       = 2000     # ticks for a full cycle through all states
    # Adds a spatial ramp on top of the temporal cycle, making the condition a
    # travelling wave rather than uniform. Off by default: uniform-in-space is
    # the escape-proof version, since a gradient lets a lineage track its
    # preferred band instead of learning to respond.
    cond_gradient::Bool    = false
    cond_strength::Float32 = 1.0f0    # gain of the cue into the sensor node
    # Making the cue legible was not enough. Channel 1 is the sole recipient of
    # external inflow, so it is the most ABUNDANT channel whatever the cue says,
    # and "always eat channel 1" stayed good enough: at 6 states the evolved
    # populations still answered ch1 in every state (1/6 matched). The reward has
    # to move with the cue, not just the exchange rate -- the same lesson the
    # zones taught, arrived at from the other side.
    #
    # Routing ALL of the inflow to the cue channel kills the dish: development
    # is fixed at birth, so at the instant the cue moves every adult is
    # mis-adapted at once and intake collapses before new births can adapt.
    # Median extinction came at the first or second switch (~2.5-4.6k ticks) in
    # 48 of 48 runs. So only a FRACTION follows the cue and the rest keeps
    # feeding channel 1: no channel is ever starved outright, the transition is
    # survivable, and tracking the cue still pays.
    cond_inflow_frac::Float32 = 0.0f0   # 0 = all inflow to channel 1, 1 = all to the cue
    # THE CONTROL THAT MATTERS. Every comparison in this project has been spoiled
    # by the manipulation changing population size as well as the task, and
    # population size drives genome size on its own. Scrambling decouples the cue
    # from the rotation it is supposed to predict, using an incommensurate period
    # so the two never lock. The dish is then ecologically IDENTICAL -- same
    # rotation, same inflow, same energy, same population -- and the only thing
    # removed is the information. If genome size rises with an informative cue
    # and not with a useless one, the task is what did it.
    cond_scramble::Bool    = false

    # --- reversible plasticity -----------------------------------------------
    # Development is cached at birth and never revisited, which makes plasticity
    # here IRREVERSIBLE. Levins' rule is that irreversible developmental
    # plasticity needs a COARSE-grained environment while fine-grained ones call
    # for REVERSIBLE responses instead, and every experiment in this project used
    # the irreversible form.
    #
    # It is also why the regime theory points at is lethal. At 1-4 generations
    # per environmental state nothing survived (0/6): an organism committed at
    # birth simply starves when the world switches under it. Re-developing during
    # life removes exactly that failure, so this is the one change with a
    # specific prediction behind it rather than a guess.
    #
    # 0 keeps the cached-at-birth behaviour and reproduces every earlier result
    # bit-for-bit. A positive value re-expresses every living organism from its
    # CURRENT surroundings every that-many ticks. Generation time is ~26 ticks,
    # so 25 is roughly one re-expression per lifetime and 5 is several.
    #
    # Deliberately FREE of any cost. Plasticity is normally costly (DeWitt et
    # al.), but giving the hypothesis its most favourable case makes a negative
    # result decisive rather than a matter of how the cost was tuned.
    redevelop_every::Int   = 0
    diff_rate::Float32   = 0.12f0
    # Per-channel diffusion falloff: channel r diffuses at
    # diff_rate * diff_falloff^(r-1). 1.0 means every channel diffuses alike.
    #
    # Biologically, metabolic byproducts are larger molecules than the primary
    # substrate and diffuse more slowly. Mechanically this is an attempt at the
    # blocker described under 'Still open': the resource fields are near-uniform,
    # so every organism sees the same mixture and there is one optimal strategy
    # and nothing to partition. Slow-diffusing waste should form persistent
    # patches that a specialist could actually live in.
    diff_falloff::Float32 = 1.0f0
    res_cap::Float32     = 6.0f0
    decay::Float32       = 0.0060f0    # abiotic loss, stops waste accumulating forever

    # --- seasonality ---------------------------------------------------------
    season_period::Int   = 15_000
    season_amp::Float32  = 0.60f0

    # --- shocks --------------------------------------------------------------
    shock_interval::Int  = 25_000      # mean generations between mass-extinction events
    shock_kill::Float32  = 0.92f0      # baseline fraction killed
    famine_len::Int      = 400         # generations of suppressed inflow after a famine
    # 0.65, not 0.12. What matters is the integrated deficit
    # famine_len * (1 - famine_factor); duration alone does nothing (measured:
    # 400x0.65 and 150x0.05 are indistinguishable). Cutting the deficit
    # 352 -> 140 buys a 28x deeper safety margin and removes extinction.
    #
    # The obvious worry — less famine means less selection and a duller model —
    # is wrong. Measured over 6 seeds x 200k ticks, softening famine left
    # cross-feeding and generation rate unchanged and RAISED clade diversity 3x
    # (Shannon 1.15 -> 3.48). Harsh famines bottleneck to 6-31 organisms, which
    # purges nearly every lineage; soft ones bottom out at 720-4885, so far more
    # clades survive. The harsh setting was not generating diversity by
    # selection, it was destroying it by bottleneck.
    famine_factor::Float32 = 0.65f0

    # --- metabolism ----------------------------------------------------------
    # The whole energy economy is scaled so a tick is roughly one cell cycle.
    # Slow division makes ticks and generations diverge badly: at max_age 3000
    # a 50k-tick run contained only ~66 generations of actual descent, far too
    # few for anything to evolve.
    # Uptake is PROPORTIONAL to what is locally available rather than capped by
    # it. A hard cap means affinity above the ambient concentration buys
    # nothing, so under scarcity affinity becomes pure cost and specialisation
    # can never pay for itself.
    # Every trait that buys something must also cost something, or evolution
    # simply pins it to the boundary. Size is charged QUADRATICALLY: its benefit
    # is multiplicative and frequency-dependent (being bigger than your
    # neighbour pays), so a linear cost just yields an arms race to 1.0.
    uptake_rate::Float32    = 0.60f0   # fraction of local resource taken at full allocation
    # Uptake scales as allocation^uptake_exponent.
    #
    # 1.0 is the tested default. Raising it was an attempt to make the
    # allocation trade-off convex, on the theory that a specialist would then
    # out-eat a generalist on its own channel and selection would become
    # disruptive. MEASURED: it does the opposite — at exponent 2 mean
    # specialisation fell to 0.32 from 0.70, i.e. the population became MORE
    # generalist, and at small grid sizes the dish starved to extinction.
    # Plausibly because aff^2 <= aff shrinks total uptake for everyone, so the
    # system is energy-poorer and a backup channel is worth more than
    # committing. Left exposed as a knob, but the hypothesis it was added for
    # is disproved; do not raise it expecting niche partitioning.
    # --- saturating uptake ----------------------------------------------------
    # Capture was strictly PROPORTIONAL to local concentration, and that is what
    # made conditional specialisation unprofitable: concentrating the allocation
    # budget on the best-paying channel forfeited access to the other five, and
    # the yield premium (1.0 against a mean of ~0.48) never covered the abundance
    # lost. Breadth beat rate, so a hand-built conditional genome was eliminated
    # even with mutation switched off (energy 0.52-0.64 of its neighbours).
    #
    # Monod saturation compresses that abundance advantage: capture goes as
    # avail / (uptake_half + avail), so a channel ten times richer is no longer
    # worth ten times as much and the yield premium can decide. This is also the
    # standard consumer-resource form; strict proportionality was the unusual
    # choice.
    #
    # 0 disables it and restores exact proportional capture, so every earlier
    # result in this README reproduces bit-for-bit.
    uptake_half::Float32    = 0.0f0
    uptake_exponent::Float32 = 1.0f0
    maint_base::Float32     = 0.0800f0
    maint_size::Float32     = 0.2200f0 # charged on size^2
    maint_motil::Float32    = 0.0350f0 # movement is not free
    maint_stress::Float32   = 0.0400f0 # maintaining shock tolerance is not free
    maint_repair::Float32   = 0.0200f0 # repair fidelity is not free
    # --- proteome cost of breadth -------------------------------------------
    # Charged as maint_channel * sum(allocation_r ^ channel_cost_exp).
    #
    # This is the term that makes ecotypes possible at all. Uptake is LINEAR in
    # allocation on a budget that sums to 1, so at the equilibrium where the
    # resources equalise, every allocation earns exactly the same — invasion
    # fitness is flat, its second derivative is zero, and the singular strategy
    # is neutral rather than disruptive. A population there just drifts as one
    # cloud, which is precisely what was observed: a single peak with 99.9% of
    # organisms topping out on channel 1.
    #
    # A concave cost (exponent < 1) breaks the tie. A specialist pays
    # maint_channel * 1; an even six-way generalist pays maint_channel * sqrt(6)
    # = 2.45x that. Uptake is untouched, so specialists lose nothing — unlike
    # raising uptake_exponent, which starved everyone. The interior becomes a
    # fitness MINIMUM, which is the textbook condition for evolutionary
    # branching: convergence stable and evolutionarily unstable.
    #
    # Biologically this is proteome allocation: every transporter you express
    # occupies membrane and ribosome capacity whether or not you are using it,
    # so breadth has a fixed overhead that specialisation avoids.
    # DEFAULT OFF (0.0), because measuring it honestly showed it costs more than
    # it buys: at 0.04 the population fell from ~12,700 to ~4,600, cross-feeding
    # fell 0.26-0.30 -> 0.22, and 1 seed in 5 went extinct — while branching went
    # from 0% to at most 2%. Tuned hard (inflow 0.60, cost 0.10) it does produce
    # a real second ecotype, but only 0-17% of the time across seeds with 1 in 6
    # dying. See the README for the full result.
    # Setting channel_cost_exp to 1.0 also disables it (the sum is then always 1).
    maint_channel::Float32       = 0.0000f0
    channel_cost_exp::Float32    = 0.50f0

    maint_crowd::Float32    = 0.0500f0 # per co-located neighbour -> local competition
    maint_load::Float32     = 0.0150f0 # per accumulated deleterious mutation
    # --- rate / yield trade-off (Pfeiffer & Bonhoeffer 2004) -----------------
    # Leaking more of what you take gives a faster uptake rate but less energy
    # per unit consumed — the reason real cells ferment to acetate instead of
    # respiring fully. P_LEAK is the evolvable dial; this is its rate bonus.
    leak_rate_bonus::Float32 = 1.20f0
    # Leakage FLOOR, and it is deliberately high. In the microbial
    # consumer-resource literature leakage is an environmental constant, not a
    # trait (Marsland et al. run l = 0.8), because it is set by biochemistry.
    # When it was free to evolve down it went to ~0.003: uptake saturates at
    # whatever is locally available, so the rate bonus cannot pay for the yield
    # loss, and the waste channels were starved of supply. The trait can still
    # raise leakage above this floor.
    leak_min::Float32        = 0.45f0

    # --- waste toxicity (Black Queen / detoxification) -----------------------
    # Accumulated byproducts inhibit growth locally, so an organism that eats a
    # waste channel is doing everyone upwind a service. This is what turns
    # one-way cross-feeding into a mutualism with negative frequency dependence;
    # without it a cross-feeder confers nothing and simply gets outcompeted.
    # Stress tolerance protects against it, which also gives that trait a
    # continuous role instead of mattering only during rare shocks.
    # 0.008 is the tested default: it leaves stress tolerance with an INTERIOR
    # optimum (~0.35 across seeds) instead of pinning it at 1.0, which a linear
    # tolerance cost did. The response is not monotonic — around 0.015 the dish
    # can die during the ramp-up, before tolerance has had time to evolve. Sweep
    # this one carefully.
    maint_toxic::Float32     = 0.0080f0
    div_energy::Float32     = 2.0f0
    e_init::Float32         = 0.8f0
    e_max::Float32          = 4.0f0    # cells cannot hoard energy indefinitely

    # --- dormancy (sporulation) ----------------------------------------------
    # Every single extinction measured was caused by a famine — 8 of 8, each
    # within 24-211 ticks of one. The reason is that organisms carry no buffer:
    # at ~0.53 mean energy against ~0.3 upkeep per tick, they hold under two
    # ticks of reserve against a 400-tick famine, so it is unsurvivable by
    # construction rather than by ecology.
    #
    # Below `dorm_threshold` energy an organism goes dormant: upkeep drops to
    # `dorm_upkeep` of normal, ageing pauses, and it cannot divide. It still
    # takes up resources (at `dorm_uptake` of normal) so it can climb back out
    # when food returns. This is sporulation, which is precisely the adaptation
    # to famine, and it makes surviving one a strategy rather than an
    # impossibility.
    dorm_threshold::Float32 = 0.15f0
    dorm_upkeep::Float32    = 0.05f0   # fraction of normal upkeep while dormant
    dorm_uptake::Float32    = 0.35f0   # fraction of normal uptake while dormant
    max_age::Int32          = Int32(60)

    # How strongly local resource concentrations feed the regulatory network.
    # 1.0 = full plasticity, 0.0 = allocation is purely genetic. Plasticity can
    # SUPPRESS diversification: if every genotype plastically converges on the
    # same behaviour, selection sees no heritable variation to split.
    plasticity::Float32     = 1.0f0

    # --- mutation ------------------------------------------------------------
    mut_rate::Float32       = 1.0f-3   # per weight site per birth, before repair
    mut_sigma::Float32      = 0.25f0
    repair_strength::Float32 = 0.85f0  # max suppression a perfect repairer achieves
    p_silence::Float32      = 0.0020f0
    p_duplicate::Float32    = 0.0020f0
    # Most mutations hurt, a few help. Load can therefore fall as well as rise,
    # which is what makes Muller's ratchet a dynamic (small populations click
    # upward, large ones purge) rather than a one-way count to extinction.
    # --- genome size ---------------------------------------------------------
    # Hidden nodes are gained by duplication and lost by deletion, so genome
    # SIZE evolves rather than only its wiring. maint_gene charges for each
    # hidden node, since more genes cost more to replicate and express: set it
    # to 0 to ask whether complexity grows by drift alone.
    ng_init::Int            = 16       # founders start here, between NG_MIN and NG_MAX
    p_node_add::Float32     = 0.0030f0
    p_node_del::Float32     = 0.0030f0
    maint_gene::Float32     = 0.0015f0 # per hidden node, per tick

    p_deleterious::Float32  = 0.70f0
    p_beneficial::Float32   = 0.05f0
    clade_radius2::Float32  = 0.16f0   # squared phenotype distance that names a new clade
    # Collect often: at grid 512 the table reached 1.4M entries between
    # collections, ~95 MB of founder phenotypes churned through the garbage
    # collector every cycle, almost all for clades with no living members.
    clade_gc_interval::Int  = 1_000    # ticks between clade-table garbage collections

    # --- logging -------------------------------------------------------------
    snapshot_decades::Int   = 12       # log-spaced snapshots per decade of generations
end

season(p::Params, gen::Int) =
    1.0f0 + p.season_amp * sin(2.0f0 * Float32(pi) * Float32(gen) / Float32(p.season_period))
