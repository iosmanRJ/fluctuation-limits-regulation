# =============================================================================
# Logging.
#
# Per-generation stats go to a flat Float32 binary file: NSTAT values per
# record, no framing. 1e6 generations is ~190 MB, which is fine, and it maps
# straight into a Matrix on read with no parsing.
#
# Full state is far too big to keep every generation, so snapshots are taken on
# a LOG-SPACED schedule. Evolution is roughly log-uniform in interest: the
# first thousand generations deserve the same number of samples as the last
# million, and a linear schedule wastes all its resolution on the boring end.
# =============================================================================

const NSTAT = 10 + NRES + 2 * NPHEN + NRES + 4   # 50

const STAT_NAMES = vcat(
    ["gen", "pop", "births", "deaths", "mean_energy", "mean_age",
     "mean_load", "clade_shannon", "clades_alive", "mean_depth"],
    ["res$(r)_total" for r in 1:NRES],
    ["mean_phen$(k)" for k in 1:NPHEN],
    ["std_phen$(k)" for k in 1:NPHEN],
    # Fraction of the population whose LARGEST uptake allocation is channel r.
    # This is the branching measure: a mean cannot distinguish one generalist
    # from two coexisting specialists, but these shares can.
    ["top_share$(r)" for r in 1:NRES],
    # Network structure: is the genome gaining or losing complexity?
    ["mean_nonzero_w", "mean_active_genes", "mean_abs_w", "mean_ng"],
)

mutable struct StatLog
    io::IO
    buf::Matrix{Float32}
    k::Int
    counts::Dict{UInt32,Int}
    dir::String
    last_flush::Float64
    flush_after::Float64
end

"""
    StatLog(dir; chunk, flush_after)

`chunk` records are buffered before writing, but the buffer is also flushed
every `flush_after` seconds. The time bound matters: a long run with a large
`stat_every` can take hundreds of thousands of ticks to fill a chunk, during
which the file stays empty and the run is indistinguishable from a stalled one.
"""
function StatLog(dir::String; chunk::Int = 512, flush_after::Real = 10.0)
    mkpath(dir)
    io = open(joinpath(dir, "stats.bin"), "w")
    StatLog(io, zeros(Float32, NSTAT, chunk), 0, Dict{UInt32,Int}(), dir,
            time(), Float64(flush_after))
end

function record!(sl::StatLog, sim::Sim)
    pop = sim.pop
    hi  = high_water(pop)

    n = 0
    sum_e = 0.0; sum_age = 0.0; sum_load = 0.0; sum_depth = 0.0
    sum_nzw = 0.0; sum_active = 0.0; sum_absw = 0.0; sum_ng = 0.0
    top = ntuple(_ -> 0, NRES)
    mean_phen = zeros(Float64, NPHEN)
    sq_phen   = zeros(Float64, NPHEN)
    counts = sl.counts
    empty!(counts)

    @inbounds for i in 1:hi
        pop.alive[i] || continue
        n += 1
        sum_e    += pop.energy[i]
        sum_age  += pop.age[i]
        sum_load += pop.load[i]
        sum_depth += pop.depth[i]
        sum_nzw += pop.nzw[i]; sum_active += pop.active[i]; sum_absw += pop.absw[i]
        sum_ng += pop.ng[i]
        # which channel does this organism allocate the most to?
        best = 1; bv = pop.phen[1, i]
        for r in 2:NRES
            v = pop.phen[r, i]
            v > bv && (bv = v; best = r)
        end
        top = Base.setindex(top, top[best] + 1, best)
        counts[pop.clade[i]] = get(counts, pop.clade[i], 0) + 1
        for k in 1:NPHEN
            v = Float64(pop.phen[k, i])
            mean_phen[k] += v
            sq_phen[k]   += v * v
        end
    end

    invn = n > 0 ? 1.0 / n : 0.0

    # Shannon entropy over clade abundances: the standard diversity index, and
    # unlike a raw clade count it is not inflated by singletons.
    H = 0.0
    if n > 0
        for (_, c) in counts
            pr = c * invn
            pr > 0 && (H -= pr * log(pr))
        end
    end

    sl.k += 1
    col = sl.k
    B = sl.buf
    @inbounds begin
        B[1, col] = sim.gen
        B[2, col] = n
        B[3, col] = sim.births
        B[4, col] = sim.deaths
        B[5, col] = sum_e * invn
        B[6, col] = sum_age * invn
        B[7, col] = sum_load * invn
        B[8, col] = H
        B[9, col]  = length(counts)
        B[10, col] = sum_depth * invn        # generations of descent, not ticks
        for r in 1:NRES
            B[10 + r, col] = total_resource(sim.env, r)
        end
        for k in 1:NPHEN
            m = mean_phen[k] * invn
            v = max(sq_phen[k] * invn - m * m, 0.0)
            B[10 + NRES + k, col]         = m
            B[10 + NRES + NPHEN + k, col] = sqrt(v)
        end
        off = 10 + NRES + 2 * NPHEN
        for r in 1:NRES
            B[off + r, col] = top[r] * invn
        end
        B[off + NRES + 1, col] = sum_nzw * invn
        B[off + NRES + 2, col] = sum_active * invn
        B[off + NRES + 3, col] = sum_absw * invn
        B[off + NRES + 4, col] = sum_ng * invn
    end

    if sl.k == size(B, 2) || time() - sl.last_flush >= sl.flush_after
        flush_stats!(sl)
    end
    return nothing
end

function flush_stats!(sl::StatLog)
    sl.last_flush = time()
    sl.k == 0 && return nothing
    write(sl.io, @view sl.buf[:, 1:sl.k])
    flush(sl.io)
    sl.k = 0
    return nothing
end

function Base.close(sl::StatLog)
    flush_stats!(sl)
    close(sl.io)
end

"""Read a stats.bin back as an NSTAT x ngen matrix."""
function read_stats(dir::String)
    path = isdir(dir) ? joinpath(dir, "stats.bin") : dir
    raw = read(path)
    ncol = length(raw) ÷ (4 * NSTAT)
    reshape(reinterpret(Float32, raw)[1:(NSTAT * ncol)], NSTAT, ncol)
end

stat_index(name::AbstractString) = findfirst(==(name), STAT_NAMES)

# -----------------------------------------------------------------------------
# Snapshots
# -----------------------------------------------------------------------------

"""
    snapshot_schedule(ngen, per_decade) -> Vector{Int}

Log-spaced generations at which to dump full state, always including gen 1 and
the final generation.
"""
function snapshot_schedule(ngen::Int, per_decade::Int)
    gens = Int[1]
    g = 1.0
    step = 10.0^(1 / per_decade)
    while g < ngen
        g *= step
        gi = round(Int, g)
        gi > last(gens) && gi <= ngen && push!(gens, gi)
    end
    last(gens) == ngen || push!(gens, ngen)
    gens
end

"""
Dump the resource grid plus a random sample of live organisms. Sampling rather
than dumping every genome keeps a snapshot near 2 MB instead of 100 MB.
"""
function write_snapshot(sim::Sim, dir::String; nsample::Int = 1500)
    mkpath(dir)
    pop = sim.pop
    hi = high_water(pop)

    live = Int32[]
    @inbounds for i in 1:hi
        pop.alive[i] && push!(live, Int32(i))
    end
    sample = length(live) <= nsample ? live :
             live[randperm(sim.rng, length(live))[1:nsample]]
    ns = length(sample)

    open(joinpath(dir, "snap_$(lpad(sim.gen, 9, '0')).bin"), "w") do io
        write(io, Int32(sim.gen), Int32(sim.env.G), Int32(NRES),
                  Int32(NPHEN), Int32(NG_MAX), Int32(ns))
        write(io, sim.env.res)
        # The TRUE occupancy grid, not one rebuilt from the organism sample.
        # A density map derived from a 1500-of-113893 sample reads as empty
        # almost everywhere and is actively misleading. This costs G^2 Int32
        # (256 KB at grid 512) — nothing next to the resource fields.
        write(io, pop.occ)
        for i in sample
            write(io, pop.x[i], pop.y[i], pop.energy[i], pop.age[i],
                      pop.load[i], pop.depth[i], pop.clade[i])
            write(io, @view pop.phen[:, i])
        end
        # genomes of the first 200 sampled organisms, for post-hoc network analysis
        ng = min(ns, 200)
        write(io, Int32(ng))
        for j in 1:ng
            write(io, @view pop.W[:, :, sample[j]])
            write(io, @view pop.bias[:, sample[j]])
        end
    end
    return nothing
end

function read_snapshot(path::String)
    open(path, "r") do io
        gen  = read(io, Int32); G = read(io, Int32); nres = read(io, Int32)
        nphen = read(io, Int32); ng = read(io, Int32); ns = read(io, Int32)
        res = Array{Float32,3}(undef, nres, G, G)
        read!(io, res)
        occ = Array{Int32,2}(undef, G, G)
        read!(io, occ)
        x = zeros(Int32, ns); y = zeros(Int32, ns); e = zeros(Float32, ns)
        age = zeros(Int32, ns); load = zeros(Int32, ns)
        depth = zeros(Int32, ns); clade = zeros(UInt32, ns)
        phen = zeros(Float32, nphen, ns)
        for i in 1:ns
            x[i] = read(io, Int32); y[i] = read(io, Int32)
            e[i] = read(io, Float32); age[i] = read(io, Int32)
            load[i] = read(io, Int32); depth[i] = read(io, Int32)
            clade[i] = read(io, UInt32)
            read!(io, @view phen[:, i])
        end
        nw = read(io, Int32)
        Ws = [Array{Float32,2}(undef, ng, ng) for _ in 1:nw]
        biases = [Vector{Float32}(undef, ng) for _ in 1:nw]
        for j in 1:nw
            read!(io, Ws[j]); read!(io, biases[j])
        end
        (; gen = Int(gen), G = Int(G), res, occ, x, y, energy = e, age, load,
           depth, clade, phen, W = Ws, bias = biases)
    end
end
