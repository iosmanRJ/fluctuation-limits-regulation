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
    # Founders, 0 = scale with dish area. A fixed count is a trap: 200 founders
    # is 0.012/cell at grid 128 but 0.0008/cell at grid 512, 16x sparser, and
    # sparse founders simply fail to establish. Measured: 3 of 8 grid-256 runs
    # died inside 90k ticks, two of them inside 27k — establishment failures,
    # not ecology.
    n_founders::Int      = 0
    seed::Int            = 42

    # --- resources -----------------------------------------------------------
    # Channel 1 is the only abiotic input. Channels 2..NRES exist ONLY as
    # metabolic waste, so every one of them is a niche that has to be created
    # by some other organism before anything can live in it.
    inflow::Float32      = 0.250f0
    yields::Vector{Float32} = Float32[1.0, 0.60, 0.45, 0.35, 0.28, 0.22]  # thermodynamic decline
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
    famine_factor::Float32 = 0.12f0

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
