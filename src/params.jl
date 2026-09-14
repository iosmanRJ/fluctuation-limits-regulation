# =============================================================================
# Every knob in one place. Mutable so the live window can tune during a run.
# =============================================================================

Base.@kwdef mutable struct Params
    # --- world ---------------------------------------------------------------
    grid::Int            = 256
    capacity::Int        = 60_000      # hard ceiling on simultaneous organisms
    seed::Int            = 42

    # --- resources -----------------------------------------------------------
    # Channel 1 is the only abiotic input. Channels 2..NRES exist ONLY as
    # metabolic waste, so every one of them is a niche that has to be created
    # by some other organism before anything can live in it.
    inflow::Float32      = 0.020f0
    yields::Vector{Float32} = Float32[1.0, 0.60, 0.45, 0.35, 0.28, 0.22]  # thermodynamic decline
    diff_rate::Float32   = 0.12f0
    res_cap::Float32     = 4.0f0
    decay::Float32       = 0.0015f0    # abiotic loss, stops waste accumulating forever

    # --- seasonality ---------------------------------------------------------
    season_period::Int   = 20_000
    season_amp::Float32  = 0.60f0

    # --- shocks --------------------------------------------------------------
    shock_interval::Int  = 50_000      # mean generations between mass-extinction events
    shock_kill::Float32  = 0.92f0      # baseline fraction killed
    famine_len::Int      = 400         # generations of suppressed inflow after a famine
    famine_factor::Float32 = 0.12f0

    # --- metabolism ----------------------------------------------------------
    uptake_scale::Float32   = 0.35f0
    maint_base::Float32     = 0.012f0
    maint_size::Float32     = 0.050f0
    maint_uptake::Float32   = 0.018f0  # cost per unit total affinity -> forces specialisation
    maint_repair::Float32   = 0.030f0  # repair fidelity is not free
    excrete_frac::Float32   = 0.35f0   # share of intake returned as waste
    div_energy::Float32     = 2.0f0
    e_init::Float32         = 0.8f0
    max_age::Int32          = Int32(3000)

    # --- mutation ------------------------------------------------------------
    mut_rate::Float32       = 2.0f-4   # per weight site per birth, before repair
    mut_sigma::Float32      = 0.25f0
    repair_strength::Float32 = 0.85f0  # max suppression a perfect repairer achieves
    p_silence::Float32      = 0.0020f0
    p_duplicate::Float32    = 0.0020f0
    clade_threshold::Int    = 8        # big-effect mutations before a new clade is named

    # --- logging -------------------------------------------------------------
    snapshot_decades::Int   = 12       # log-spaced snapshots per decade of generations
end

season(p::Params, gen::Int) =
    1.0f0 + p.season_amp * sin(2.0f0 * Float32(pi) * Float32(gen) / Float32(p.season_period))
