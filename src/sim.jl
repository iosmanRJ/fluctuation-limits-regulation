# =============================================================================
# The simulation loop.
#
# Order within one generation:
#   replenish -> maybe shock -> organism step (collect births) -> settle births
#   -> diffuse -> log
#
# Births are collected during the organism loop and instantiated afterwards, so
# a newborn never acts in the generation it was born in. Instantiating inline
# would also let a newborn land in a slot the loop has not yet reached, which
# silently gives some offspring a free extra tick.
# =============================================================================

mutable struct Sim
    p::Params
    env::Environment
    pop::Population
    clades::CladeTable
    rng::Xoshiro
    gen::Int

    # scratch, reused every generation so the hot loop allocates nothing
    birth_parents::Vector{Int32}
    e_scratch::Vector{Float32}
    buf_scratch::Vector{Float32}
    phen_scratch::Vector{Float32}
    env_scratch::Vector{Float32}

    births::Int
    deaths::Int
    dropped::Int          # births lost to a full organism array — must stay ~0
    last_shock::Int
    shock_log::Vector{Tuple{Int,Symbol}}
end

function Sim(p::Params = Params())
    rng = Xoshiro(p.seed)
    env = Environment(p)
    pop = Population(p)
    sim = Sim(p, env, pop, CladeTable(), rng, 0,
              Int32[], zeros(Float32, NG), zeros(Float32, NG),
              zeros(Float32, NPHEN), zeros(Float32, NRES),
              0, 0, 0, 0, Tuple{Int,Symbol}[])
    seed_population!(sim, 200)
    return sim
end

"""Drop `n` founder organisms with random regulatory networks into the dish."""
function seed_population!(sim::Sim, n::Int)
    p, pop, env, rng = sim.p, sim.pop, sim.env, sim.rng
    G = env.G
    for _ in 1:n
        i = claim!(pop)
        i == 0 && break
        seed_genome!(pop.W, pop.bias, i, rng)
        pop.x[i] = rand(rng, Int32(1):Int32(G))
        pop.y[i] = rand(rng, Int32(1):Int32(G))
        pop.occ[pop.x[i], pop.y[i]] += Int32(1)
        pop.energy[i] = p.e_init
        pop.age[i] = 0
        pop.load[i] = 0
        pop.depth[i] = 0
        pop.clade[i] = UInt32(0)
        @inbounds for r in 1:NRES
            sim.env_scratch[r] = env.res[r, pop.x[i], pop.y[i]]
        end
        develop!(sim.phen_scratch, sim.e_scratch, sim.buf_scratch,
                 pop.W, pop.bias, i, sim.env_scratch)
        @inbounds for k in 1:NPHEN
            pop.phen[k, i] = sim.phen_scratch[k]
        end
    end
    return nothing
end

# -----------------------------------------------------------------------------

function step!(sim::Sim)
    p, env, pop, rng = sim.p, sim.env, sim.pop, sim.rng
    sim.gen += 1
    sim.births = 0
    sim.deaths = 0

    replenish!(env, p, sim.gen)
    maybe_shock!(sim)
    organism_step!(sim)
    settle_births!(sim)
    diffuse!(env, p)

    if p.clade_gc_interval > 0 && sim.gen % p.clade_gc_interval == 0
        prune_clades!(sim.clades, pop)
    end
    return sim
end

function organism_step!(sim::Sim)
    p, env, pop, rng = sim.p, sim.env, sim.pop, sim.rng
    res = env.res
    G = env.G
    hi = high_water(pop)
    empty!(sim.birth_parents)

    # hoist every parameter: these are fields of a mutable struct, so the
    # compiler must reload them on each access unless we pull them out first.
    yields        = p.yields
    uptake_rate   = p.uptake_rate
    maint_base    = p.maint_base
    maint_size    = p.maint_size
    maint_motil   = p.maint_motil
    maint_stress  = p.maint_stress
    maint_repair  = p.maint_repair
    maint_crowd   = p.maint_crowd
    maint_load    = p.maint_load
    excrete_frac  = p.excrete_frac
    div_energy    = p.div_energy
    res_cap       = p.res_cap
    e_max         = p.e_max
    max_age       = p.max_age
    occ           = pop.occ
    deaths        = 0

    @inbounds for i in 1:hi
        pop.alive[i] || continue

        xi = pop.x[i]
        yi = pop.y[i]

        # --- uptake across channels, weighted by evolved affinity ------------
        gained     = 0.0f0
        raw_intake = 0.0f0
        affinity_total = 0.0f0
        # Size scales uptake as well as upkeep. Without a benefit, size has
        # only a cost and evolution trivially drives it to zero; with both it
        # has an interior optimum that shifts as conditions change.
        size_gain = 0.5f0 + pop.phen[P_SIZE, i]
        for r in 1:NRES
            aff = pop.phen[r, i]
            affinity_total += aff
            avail = res[r, xi, yi]
            # proportional, not capped: more affinity always captures more, and
            # co-located organisms therefore compete for shares of the same pool
            got = avail * uptake_rate * aff * size_gain
            if got > 0.0f0
                res[r, xi, yi] = avail - got
                raw_intake += got
                gained += got * yields[r]      # lower-quality wastes yield less
            end
        end

        # --- excretion: this organism's waste is someone else's niche --------
        if raw_intake > 0.0f0
            waste = raw_intake * excrete_frac
            for r in 1:NRES
                frac = pop.phen[NRES + r, i]
                if frac > 0.0f0
                    v = res[r, xi, yi] + waste * frac
                    res[r, xi, yi] = v > res_cap ? res_cap : v
                end
            end
        end

        # --- maintenance: broad diets and good repair both cost energy -------
        crowd = Float32(occ[xi, yi] - Int32(1))     # neighbours sharing this cell
        sz = pop.phen[P_SIZE, i]
        upkeep = maint_base +
                 maint_size   * sz * sz +             # quadratic: bounds the size arms race
                 maint_motil  * pop.phen[P_MOTIL, i] +
                 maint_stress * pop.phen[P_STRESS, i] +
                 maint_repair * pop.phen[P_REPAIR, i] +
                 maint_crowd  * crowd +
                 maint_load   * Float32(pop.load[i])   # genetic deterioration bites here

        e = pop.energy[i] + gained * (1.0f0 - excrete_frac) - upkeep
        e = e > e_max ? e_max : e                      # no indefinite hoarding
        pop.energy[i] = e
        pop.age[i] += Int32(1)

        # --- movement ---------------------------------------------------------
        if rand(rng, Float32) < pop.phen[P_MOTIL, i]
            nx = xi + rand(rng, Int32(-1):Int32(1))
            ny = yi + rand(rng, Int32(-1):Int32(1))
            nx = nx < Int32(1) ? Int32(1) : (nx > Int32(G) ? Int32(G) : nx)
            ny = ny < Int32(1) ? Int32(1) : (ny > Int32(G) ? Int32(G) : ny)
            place!(pop, i, nx, ny)
        end

        # --- death and division ------------------------------------------------
        if e <= 0.0f0 || pop.age[i] > max_age
            kill!(pop, i)
            deaths += 1
        else
            # division threshold is itself an evolved trait
            thresh = div_energy * (0.5f0 + pop.phen[P_DIVIDE, i])
            if e >= thresh
                push!(sim.birth_parents, Int32(i))
            end
        end
    end

    sim.deaths = deaths
    return nothing
end

function settle_births!(sim::Sim)
    p, env, pop, rng = sim.p, sim.env, sim.pop, sim.rng
    G = env.G
    births = 0

    @inbounds for pi_ in sim.birth_parents
        parent = Int(pi_)
        pop.alive[parent] || continue          # parent died to a shock meanwhile
        child = claim!(pop)
        if child == 0
            # The organism array is full. This is an IMPLEMENTATION limit, not
            # an ecological one, and if it happens at scale the run is invalid:
            # dropped births mean reproductive success stops being selected on.
            sim.dropped += length(sim.birth_parents)
            break
        end

        # child appears next to the parent, but only if there is room there
        nx = pop.x[parent] + rand(rng, Int32(-1):Int32(1))
        ny = pop.y[parent] + rand(rng, Int32(-1):Int32(1))
        nx = nx < Int32(1) ? Int32(1) : (nx > Int32(G) ? Int32(G) : nx)
        ny = ny < Int32(1) ? Int32(1) : (ny > Int32(G) ? Int32(G) : ny)
        if pop.occ[nx, ny] >= p.max_per_cell
            kill!(pop, child)                  # no room; the division fails
            continue
        end

        repair = pop.phen[P_REPAIR, parent]
        nmut, _big, dload = inherit!(pop.W, pop.bias, child, parent, rng, p, repair)

        pop.x[child] = nx
        pop.y[child] = ny
        pop.occ[nx, ny] += Int32(1)

        # parent splits its energy with the child
        half = pop.energy[parent] * 0.5f0
        pop.energy[parent] = half
        pop.energy[child]  = half
        pop.age[child]     = Int32(0)
        pop.depth[child]   = pop.depth[parent] + Int32(1)
        # load can fall as well as rise: beneficial mutations undo damage
        pop.load[child]    = max(pop.load[parent] + Int32(dload), Int32(0))

        # --- develop the child's network in ITS OWN local conditions ---------
        for r in 1:NRES
            sim.env_scratch[r] = env.res[r, pop.x[child], pop.y[child]]
        end
        develop!(sim.phen_scratch, sim.e_scratch, sim.buf_scratch,
                 pop.W, pop.bias, child, sim.env_scratch)
        for k in 1:NPHEN
            pop.phen[k, child] = sim.phen_scratch[k]
        end

        # --- clade bookkeeping -------------------------------------------------
        # A child that has drifted far enough from its clade's founder
        # phenotype founds a new clade. Radius, not mutation count: the same
        # number of mutations can mean nothing or everything depending on where
        # they land in the network.
        parent_clade = pop.clade[parent]
        if phen_dist2(sim.clades, parent_clade, view(pop.phen, :, child)) > p.clade_radius2
            pop.clade[child] = new_clade!(sim.clades, parent_clade, sim.gen,
                                          view(pop.phen, :, child))
        else
            pop.clade[child] = parent_clade
        end
        births += 1
    end

    sim.births = births
    return nothing
end

# -----------------------------------------------------------------------------
# Mass-extinction events. Three flavours so the survivors are not always the
# same kind of organism.
# -----------------------------------------------------------------------------

function maybe_shock!(sim::Sim)
    p, pop, rng = sim.p, sim.pop, sim.rng
    p.shock_interval <= 0 && return nothing
    rand(rng, Float32) < (1.0f0 / p.shock_interval) || return nothing

    kind = rand(rng, (:cull, :heat, :famine))
    hi = high_water(pop)

    if kind === :cull
        # indiscriminate: a genetic bottleneck, nothing to adapt to
        @inbounds for i in 1:hi
            pop.alive[i] || continue
            if rand(rng, Float32) < p.shock_kill
                kill!(pop, i)
            end
        end
    elseif kind === :heat
        # selective: stress tolerance decides who lives
        @inbounds for i in 1:hi
            pop.alive[i] || continue
            if rand(rng, Float32) < p.shock_kill * (1.0f0 - pop.phen[P_STRESS, i])
                kill!(pop, i)
            end
        end
    else
        # famine: nobody dies now, but inflow collapses for a while
        sim.env.famine = p.famine_len
    end

    sim.last_shock = sim.gen
    push!(sim.shock_log, (sim.gen, kind))
    return nothing
end
