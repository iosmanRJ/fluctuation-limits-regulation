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

    # Pending births: slots and seeds are assigned serially, then the expensive
    # per-birth work (genome copy, mutation, network development) runs in
    # parallel. Each child writes only its OWN genome slot, so the parallel
    # phase has no shared writes at all.
    pend_parent::Vector{Int32}
    pend_child::Vector{Int32}
    pend_seed::Vector{UInt64}
    pend_dload::Vector{Int32}
    # one scratch set per worker chunk
    tl_e::Vector{Vector{Float32}}
    tl_buf::Vector{Vector{Float32}}
    tl_phen::Vector{Vector{Float32}}
    tl_env::Vector{Vector{Float32}}
    tl_rng::Vector{Xoshiro}        # re-seeded in place; constructing one per birth allocates

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
    nt = Threads.nthreads()
    sim = Sim(p, env, pop, CladeTable(), rng, 0,
              Int32[], zeros(Float32, NG), zeros(Float32, NG),
              zeros(Float32, NPHEN), zeros(Float32, NRES),
              Int32[], Int32[], UInt64[], Int32[],
              [zeros(Float32, NG)    for _ in 1:nt],
              [zeros(Float32, NG)    for _ in 1:nt],
              [zeros(Float32, NPHEN) for _ in 1:nt],
              [zeros(Float32, NRES)  for _ in 1:nt],
              [Xoshiro(0) for _ in 1:nt],
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
    uptake_exp    = p.uptake_exponent
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

        # --- uptake across channels, by evolved allocation --------------------
        gained     = 0.0f0
        raw_intake = 0.0f0
        # Size scales uptake as well as upkeep. Without a benefit, size has
        # only a cost and evolution trivially drives it to zero; with both it
        # has an interior optimum that shifts as conditions change.
        size_gain = 0.5f0 + pop.phen[P_SIZE, i]
        for r in 1:NRES
            aff = pop.phen[r, i]          # allocation share; these sum to 1
            avail = res[r, xi, yi]
            # Proportional, not capped: more allocation always captures more,
            # and co-located organisms compete for shares of the same pool.
            # Raised to uptake_exp so the trade-off is CONVEX — a specialist
            # out-eats a generalist on its own channel. Under a linear
            # trade-off spreading costs nothing in aggregate, there is no
            # disruptive selection, and the population collapses to one
            # monomorphic generalist whose waste allocation then erodes away.
            got = avail * uptake_rate * (aff^uptake_exp) * size_gain
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

    # ---- phase 1 (serial, cheap): claim slots, place children, draw seeds ---
    # Everything here touches shared state — the free list, the occupancy grid,
    # the RNG — so it stays serial. It is only a few operations per birth.
    empty!(sim.pend_parent); empty!(sim.pend_child); empty!(sim.pend_seed)
    @inbounds for pi_ in sim.birth_parents
        parent = Int(pi_)
        pop.alive[parent] || continue          # parent died to a shock meanwhile

        nx = pop.x[parent] + rand(rng, Int32(-1):Int32(1))
        ny = pop.y[parent] + rand(rng, Int32(-1):Int32(1))
        nx = nx < Int32(1) ? Int32(1) : (nx > Int32(G) ? Int32(G) : nx)
        ny = ny < Int32(1) ? Int32(1) : (ny > Int32(G) ? Int32(G) : ny)
        pop.occ[nx, ny] >= p.max_per_cell && continue   # no room; division fails

        child = claim!(pop)
        if child == 0
            # The organism array is full. This is an IMPLEMENTATION limit, not
            # an ecological one, and if it happens at scale the run is invalid:
            # dropped births mean reproductive success stops being selected on.
            sim.dropped += length(sim.birth_parents)
            break
        end

        pop.x[child] = nx
        pop.y[child] = ny
        pop.occ[nx, ny] += Int32(1)

        # Each child gets its own RNG seed, drawn serially in a fixed order.
        # That is what keeps results EXACTLY reproducible regardless of how
        # many threads run phase 2 — per-thread RNGs would make the stream
        # depend on the thread count.
        push!(sim.pend_parent, Int32(parent))
        push!(sim.pend_child, Int32(child))
        push!(sim.pend_seed, rand(rng, UInt64))
    end

    np = length(sim.pend_child)
    resize!(sim.pend_dload, np)

    # ---- phase 2 (parallel, expensive): mutate and develop -----------------
    # ~72% of total runtime lived here and ran on one core. Each task writes
    # only to its own children's genome and phenotype columns, so there are no
    # shared writes; env is read-only for the duration.
    if np > 0
        nchunk = min(Threads.nthreads(), np)
        Threads.@threads for c in 1:nchunk
            e    = sim.tl_e[c]
            buf  = sim.tl_buf[c]
            phen = sim.tl_phen[c]
            envs = sim.tl_env[c]
            crng = sim.tl_rng[c]
            lo = ((c - 1) * np) ÷ nchunk + 1
            hi = (c * np) ÷ nchunk
            @inbounds for k in lo:hi
                parent = Int(sim.pend_parent[k])
                child  = Int(sim.pend_child[k])
                # allocation-free re-seed; seeding from the child's OWN seed
                # keeps the random stream independent of the thread count
                seed_fast!(crng, sim.pend_seed[k])

                repair = pop.phen[P_REPAIR, parent]
                _nmut, _big, dload = inherit!(pop.W, pop.bias, child, parent,
                                              crng, p, repair)
                sim.pend_dload[k] = Int32(dload)

                for r in 1:NRES
                    envs[r] = env.res[r, pop.x[child], pop.y[child]]
                end
                develop!(phen, e, buf, pop.W, pop.bias, child, envs)
                for kk in 1:NPHEN
                    pop.phen[kk, child] = phen[kk]
                end
            end
        end
    end

    # ---- phase 3 (serial, cheap): energy split, lineage and clade bookkeeping
    # The clade table is shared and grows, so assignment stays serial.
    @inbounds for k in 1:np
        parent = Int(sim.pend_parent[k])
        child  = Int(sim.pend_child[k])

        half = pop.energy[parent] * 0.5f0
        pop.energy[parent] = half
        pop.energy[child]  = half
        pop.age[child]     = Int32(0)
        pop.depth[child]   = pop.depth[parent] + Int32(1)
        # load can fall as well as rise: beneficial mutations undo damage
        pop.load[child]    = max(pop.load[parent] + sim.pend_dload[k], Int32(0))

        # A child that has drifted far enough from its clade's founder
        # phenotype founds a new clade. Radius, not mutation count: the same
        # number of mutations can mean nothing or everything depending on where
        # they land in the network.
        parent_clade = pop.clade[parent]
        if phen_dist2(sim.clades, parent_clade, pop.phen, child) > p.clade_radius2
            pop.clade[child] = new_clade!(sim.clades, parent_clade, sim.gen,
                                          pop.phen, child)
        else
            pop.clade[child] = parent_clade
        end
    end

    sim.births = np
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
