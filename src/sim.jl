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
    pend_nid::Vector{Int32}       # pre-reserved node ids, one per pending birth
    # one scratch set per worker chunk
    tl_e::Vector{Vector{Float32}}
    tl_buf::Vector{Vector{Float32}}
    tl_phen::Vector{Vector{Float32}}
    tl_env::Vector{Vector{Float32}}
    tl_rng::Vector{Xoshiro}        # re-seeded in place; constructing one per birth allocates

    births::Int
    deaths::Int
    dropped::Int          # births lost to a full organism array — must stay ~0
    next_node_id::Int     # hands out unique ids so individual nodes can be tracked
    last_shock::Int
    shock_log::Vector{Tuple{Int,Symbol}}
end

function Sim(p::Params = Params())
    rng = Xoshiro(p.seed)
    env = Environment(p)
    pop = Population(p)
    nt = Threads.nthreads()
    # Independent per-state yield permutations, drawn from a stream of their own
    # so they do not perturb the main RNG and a run stays reproducible.
    if p.independent_states && size(p.yield_perm, 2) != max(1, p.cond_states)
        prng = Xoshiro(p.seed ⊻ 0x5EEDBEEF)
        ns = max(1, p.cond_states)
        perm = zeros(Int32, NRES, ns)
        for s_ in 1:ns
            perm[:, s_] = Int32.(randperm(prng, NRES))
        end
        p.yield_perm = perm
    end
    sim = Sim(p, env, pop, CladeTable(), rng, 0,
              Int32[], zeros(Float32, NG_MAX), zeros(Float32, NG_MAX),
              zeros(Float32, NPHEN), zeros(Float32, NRES),
              Int32[], Int32[], UInt64[], Int32[], Int32[],
              [zeros(Float32, NG_MAX) for _ in 1:nt],
              [zeros(Float32, NG_MAX) for _ in 1:nt],
              [zeros(Float32, NPHEN) for _ in 1:nt],
              [zeros(Float32, NRES)  for _ in 1:nt],
              [Xoshiro(0) for _ in 1:nt],
              0, 0, 0, 0, 1, Tuple{Int,Symbol}[])
    seed_population!(sim, founder_count(p))
    return sim
end

"""Founders to seed: explicit if set, else scaled to dish area."""
founder_count(p::Params) =
    p.n_founders > 0 ? p.n_founders : max(200, (p.grid * p.grid) ÷ 80)

"""Drop `n` founder organisms with random regulatory networks into the dish."""
function seed_population!(sim::Sim, n::Int)
    p, pop, env, rng = sim.p, sim.pop, sim.env, sim.rng
    G = env.G
    for _ in 1:n
        i = claim!(pop)
        i == 0 && break
        pop.ng[i] = Int32(clamp(p.ng_init, NG_MIN, NG_MAX))
        for q in 1:NG_MAX
            pop.node_born[q, i] = Int32(0)
            pop.node_id[q, i] = Int32(q)      # founders share ids 1..NG_MAX
        end
        seed_genome!(pop.W, pop.bias, i, rng, Int(pop.ng[i]))
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
                 pop.W, pop.bias, i, sim.env_scratch, Int(pop.ng[i]);
                 plasticity = p.plasticity,
                 cond = condition_at(p, pop.x[i], env.G, sim.gen)[1])
        @inbounds for k in 1:NPHEN
            pop.phen[k, i] = sim.phen_scratch[k]
        end
        pop.chan_cost[i] = channel_cost(sim.phen_scratch, p.channel_cost_exp)
        nz, act, aw = genome_structure(pop.W, pop.bias, i, Int(pop.ng[i]))
        pop.nzw[i] = nz; pop.active[i] = act; pop.absw[i] = aw
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
    redevelop!(sim)
    organism_step!(sim)
    settle_births!(sim)
    diffuse!(env, p)

    if p.clade_gc_interval > 0 && sim.gen % p.clade_gc_interval == 0
        prune_clades!(sim.clades, pop)
    end
    return sim
end

"""
    redevelop!(sim)

Re-express every living organism from its CURRENT surroundings, making plasticity
reversible instead of fixed at birth.

Development is otherwise run once and cached, so an organism is committed for
life to the conditions it was born into. That is the irreversible form, and
Levins' rule is that it suits coarse-grained environments while fine-grained ones
call for reversible responses. `redevelop_every = 0` keeps the cached behaviour
and reproduces earlier results exactly.

`develop!` draws no randomness, so this changes no RNG stream and determinism is
preserved.
"""
function redevelop!(sim::Sim)
    p = sim.p
    every = p.redevelop_every
    (every <= 0 || sim.gen % every != 0) && return nothing
    env, pop = sim.env, sim.pop
    G = env.G
    plast = p.plasticity
    chan_exp = p.channel_cost_exp
    ph = sim.phen_scratch; e = sim.e_scratch; buf = sim.buf_scratch
    es = sim.env_scratch
    @inbounds for i in 1:high_water(pop)
        pop.alive[i] || continue
        xi = pop.x[i]; yi = pop.y[i]
        for r in 1:NRES
            es[r] = env.res[r, xi, yi]
        end
        cue = condition_at(p, xi, G, sim.gen)[1]
        develop!(ph, e, buf, pop.W, pop.bias, i, es, Int(pop.ng[i]);
                 plasticity = plast, cond = cue)
        for k in 1:NPHEN
            pop.phen[k, i] = ph[k]
        end
        pop.chan_cost[i] = channel_cost(ph, chan_exp)
    end
    return nothing
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
    nzones        = p.yield_zones
    # zones for the YIELD rotation specifically; 1 disables it while leaving the
    # inflow rotation (and hence the spatial structure) untouched
    yzones        = p.rotate_yield ? nzones : 1
    indep         = p.independent_states && !isempty(p.yield_perm)
    yperm         = p.yield_perm
    cond_states   = p.cond_states
    uptake_rate   = p.uptake_rate
    uptake_exp    = p.uptake_exponent
    up_half       = p.uptake_half
    maint_base    = p.maint_base
    maint_size    = p.maint_size
    maint_motil   = p.maint_motil
    maint_stress  = p.maint_stress
    maint_repair  = p.maint_repair
    maint_gene    = p.maint_gene
    maint_chan    = p.maint_channel
    maint_crowd   = p.maint_crowd
    maint_load    = p.maint_load
    maint_toxic   = p.maint_toxic
    leak_bonus    = p.leak_rate_bonus
    leak_min      = p.leak_min
    div_energy    = p.div_energy
    res_cap       = p.res_cap
    e_max         = p.e_max
    max_age       = p.max_age
    dorm_thresh   = p.dorm_threshold
    dorm_upkeep   = p.dorm_upkeep
    dorm_uptake   = p.dorm_uptake
    occ           = pop.occ
    deaths        = 0

    @inbounds for i in 1:hi
        pop.alive[i] || continue

        xi = pop.x[i]
        yi = pop.y[i]

        # --- uptake, conversion and leakage ----------------------------------
        gained     = 0.0f0
        toxic      = 0.0f0
        # Dormancy is decided from the energy carried INTO the tick, so it
        # cannot flicker within one.
        dormant = pop.energy[i] < dorm_thresh
        # Size scales uptake as well as upkeep. Without a benefit, size has
        # only a cost and evolution trivially drives it to zero; with both it
        # has an interior optimum that shifts as conditions change.
        zone = zone_of(nzones, xi, G)
        # The zone and the abiotic condition both rotate the yield vector; their
        # offsets add. The condition cannot be depleted, so unlike the zone it
        # stays legible however hard the dish is exploited.
        cstate = condition_at(p, xi, G, sim.gen)[2]
        yoff = (yzones > 1 ? zone - 1 : 0) + (cond_states > 1 ? cstate - 1 : 0)
        size_gain = 0.5f0 + pop.phen[P_SIZE, i]
        # Rate/yield trade-off: leaking more lets a cell pull resource in faster
        # but extract less from each unit of it.
        leak = leak_min + (1.0f0 - leak_min) * pop.phen[P_LEAK, i]
        rate_mult = size_gain * (1.0f0 + leak_bonus * leak)
        dormant && (rate_mult *= dorm_uptake)

        for r in 1:NRES
            aff = pop.phen[r, i]          # allocation share; these sum to 1
            avail = res[r, xi, yi]
            # waste channels inhibit growth wherever they accumulate
            r > 1 && (toxic += avail)
            aff <= 0.0f0 && continue
            # More allocation always captures more, and co-located organisms
            # compete for shares of the same pool. With uptake_half > 0 the
            # concentration term saturates, so being on the richest channel
            # stops being decisive and yield can matter instead.
            # Normalised so the two forms AGREE at avail = uptake_half: below it
            # capture is near-proportional, above it flattens. Without this the
            # cap sits at uptake_rate and the dish starves (population 434 -> 5).
            conc = up_half > 0.0f0 ?
                   2.0f0 * up_half * avail / (up_half + avail) : avail
            got = conc * uptake_rate * (aff^uptake_exp) * rate_mult
            got <= 0.0f0 && continue
            got > avail && (got = avail)
            res[r, xi, yi] = avail - got
            # Rotations share structure across states; independent permutations
            # do not, so each state must be handled on its own terms.
            yv = indep ? @inbounds(yields[yperm[r, cstate]]) :
                         rot_yield(yields, yoff, r)
            gained += got * yv * (1.0f0 - leak)

            # --- the conversion ladder ---------------------------------------
            # Channel r's byproduct is channel r+1, always. A resource is never
            # converted back into itself or into anything upstream: that is the
            # excluded diagonal of the metabolic matrix, and the reason is
            # thermodynamic — you cannot climb back up the energy ladder. The
            # last channel's byproduct leaves the system entirely.
            if r < NRES
                v = res[r + 1, xi, yi] + got * leak
                res[r + 1, xi, yi] = v > res_cap ? res_cap : v
            end
        end

        # --- maintenance: broad diets and good repair both cost energy -------
        crowd = Float32(occ[xi, yi] - Int32(1))     # neighbours sharing this cell
        sz = pop.phen[P_SIZE, i]
        st = pop.phen[P_STRESS, i]
        upkeep = maint_base +
                 maint_size   * sz * sz +             # quadratic: bounds the size arms race
                 maint_motil  * pop.phen[P_MOTIL, i] +
                 maint_stress * st * st +            # quadratic, like size: a linear
                                                       # cost gives tolerance no interior optimum

                 maint_repair * pop.phen[P_REPAIR, i] +
                 maint_gene   * Float32(pop.ng[i] - NG_MIN) +  # cost of extra genes
                 maint_chan   * pop.chan_cost[i] +   # proteome cost of breadth
                 maint_crowd  * crowd +
                 maint_load   * Float32(pop.load[i]) +  # genetic deterioration bites here
                 # Waste toxicity, which stress tolerance buys protection from.
                 # This is what makes eating a waste channel a SERVICE to
                 # whoever produced it, rather than pure scavenging.
                 maint_toxic  * toxic * (1.0f0 - pop.phen[P_STRESS, i])

        # a spore costs almost nothing to maintain
        dormant && (upkeep *= dorm_upkeep)

        e = pop.energy[i] + gained - upkeep
        e = e > e_max ? e_max : e                      # no indefinite hoarding
        pop.energy[i] = e
        # ageing pauses while dormant: a spore waits out the famine
        dormant || (pop.age[i] += Int32(1))

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
            # division threshold is itself an evolved trait; a dormant cell
            # does not divide however much it happens to be holding
            thresh = div_energy * (0.5f0 + pop.phen[P_DIVIDE, i])
            if !dormant && e >= thresh
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
    chan_exp = p.channel_cost_exp
    plast    = p.plasticity

    # ---- phase 1 (serial, cheap): claim slots, place children, draw seeds ---
    # Everything here touches shared state — the free list, the occupancy grid,
    # the RNG — so it stays serial. It is only a few operations per birth.
    empty!(sim.pend_parent); empty!(sim.pend_child); empty!(sim.pend_seed)
    empty!(sim.pend_nid)
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
        # reserve an id per pending birth in the SERIAL phase; the parallel
        # phase must not share a mutable counter
        push!(sim.pend_nid, Int32(sim.next_node_id))
        sim.next_node_id += 1
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
                nid = sim.pend_nid[k]
                ng_p = Int(pop.ng[parent])
                _nmut, _big, dload, ng_c, added, moved =
                    inherit!(pop.W, pop.bias, child, parent, crng, p, repair, ng_p)
                pop.ng[child] = Int32(ng_c)
                # carry each node's birth tick across, then record what changed
                for q in 1:max(ng_p, ng_c)
                    pop.node_born[q, child] = pop.node_born[q, parent]
                    pop.node_id[q, child]   = pop.node_id[q, parent]
                end
                if moved != 0                       # a deletion swapped the last
                    pop.node_born[moved, child] = pop.node_born[ng_c + 1, parent]
                    pop.node_id[moved, child]   = pop.node_id[ng_c + 1, parent]
                end
                if added != 0                       # a brand new node
                    pop.node_born[added, child] = Int32(sim.gen)
                    pop.node_id[added, child]   = Int32(nid)
                end
                sim.pend_dload[k] = Int32(dload)

                for r in 1:NRES
                    envs[r] = env.res[r, pop.x[child], pop.y[child]]
                end
                # pure in (p, x, G, tick), so safe to call from the parallel phase
                ccue = condition_at(p, pop.x[child], G, sim.gen)[1]
                # An organism develops under the condition it is BORN into and
                # is committed to it for life, so the genome -- not the
                # individual -- has to hold a map that works in every state.
                develop!(phen, e, buf, pop.W, pop.bias, child, envs, ng_c;
                         plasticity = plast, cond = ccue)
                for kk in 1:NPHEN
                    pop.phen[kk, child] = phen[kk]
                end
                pop.chan_cost[child] = channel_cost(phen, chan_exp)
                nz, act, aw = genome_structure(pop.W, pop.bias, child, ng_c)
                pop.nzw[child] = nz; pop.active[child] = act; pop.absw[child] = aw
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
