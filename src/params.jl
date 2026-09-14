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
    seed::Int            = 42

    # --- resources -----------------------------------------------------------
    # Channel 1 is the only abiotic input. Channels 2..NRES exist ONLY as
    # metabolic waste, so every one of them is a niche that has to be created
    # by some other organism before anything can live in it.
    inflow::Float32      = 0.250f0
    yields::Vector{Float32} = Float32[1.0, 0.60, 0.45, 0.35, 0.28, 0.22]  # thermodynamic decline
    diff_rate::Float32   = 0.12f0
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
    maint_base::Float32     = 0.0800f0
    maint_size::Float32     = 0.2200f0 # charged on size^2
    maint_motil::Float32    = 0.0350f0 # movement is not free
    maint_stress::Float32   = 0.0400f0 # maintaining shock tolerance is not free
    maint_repair::Float32   = 0.0200f0 # repair fidelity is not free
    maint_crowd::Float32    = 0.0500f0 # per co-located neighbour -> local competition
    maint_load::Float32     = 0.0150f0 # per accumulated deleterious mutation
    excrete_frac::Float32   = 0.35f0   # share of intake returned as waste
    div_energy::Float32     = 2.0f0
    e_init::Float32         = 0.8f0
    e_max::Float32          = 4.0f0    # cells cannot hoard energy indefinitely
    max_age::Int32          = Int32(60)

    # --- mutation ------------------------------------------------------------
    mut_rate::Float32       = 1.0f-3   # per weight site per birth, before repair
    mut_sigma::Float32      = 0.25f0
    repair_strength::Float32 = 0.85f0  # max suppression a perfect repairer achieves
    p_silence::Float32      = 0.0020f0
    p_duplicate::Float32    = 0.0020f0
    # Most mutations hurt, a few help. Load can therefore fall as well as rise,
    # which is what makes Muller's ratchet a dynamic (small populations click
    # upward, large ones purge) rather than a one-way count to extinction.
    p_deleterious::Float32  = 0.70f0
    p_beneficial::Float32   = 0.05f0
    clade_radius2::Float32  = 0.16f0   # squared phenotype distance that names a new clade
    clade_gc_interval::Int  = 5_000    # ticks between clade-table garbage collections

    # --- logging -------------------------------------------------------------
    snapshot_decades::Int   = 12       # log-spaced snapshots per decade of generations
end

season(p::Params, gen::Int) =
    1.0f0 + p.season_amp * sin(2.0f0 * Float32(pi) * Float32(gen) / Float32(p.season_period))
