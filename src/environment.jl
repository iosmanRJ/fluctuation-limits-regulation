# =============================================================================
# The dish: NRES stacked resource fields on a G x G grid with no-flux edges.
# =============================================================================

mutable struct Environment
    res::Array{Float32,3}     # NRES x G x G
    buf::Array{Float32,3}
    G::Int
    famine::Int               # generations of famine still to run
end

function Environment(p::Params)
    G = p.grid
    res = zeros(Float32, NRES, G, G)
    @inbounds for j in 1:G, i in 1:G
        res[1, i, j] = 1.0f0          # seed the dish with primary resource only
    end
    Environment(res, similar(res), G, 0)
end

"""
    replenish!(env, p, gen)

Add abiotic inflow to channel 1 (seasonally modulated, suppressed during
famine) and apply decay to every channel.
"""
function replenish!(env::Environment, p::Params, gen::Int)
    G = env.G
    nz = p.yield_zones
    # The abiotic cue's rotation, if it is also steering the inflow. This is
    # uniform in space, so it shifts the whole dish at once.
    cfrac = (p.cond_states > 1) ? clamp(p.cond_inflow_frac, 0.0f0, 1.0f0) : 0.0f0
    coff  = cfrac > 0.0f0 ? condition_at(p, Int32(1), G, gen)[2] - 1 : 0
    indep = p.independent_states && !isempty(p.yield_perm)
    infl = p.inflow * season(p, gen)
    if env.famine > 0
        infl *= p.famine_factor
        env.famine -= 1
    end
    infl = max(infl, 0.0f0)
    keep = 1.0f0 - p.decay
    res = env.res
    @inbounds for j in 1:G
        # Which channel receives the abiotic inflow HERE. With rotate_inflow the
        # zone's best-paying channel is also the one being fed, so the resource
        # landscape genuinely differs between zones and an organism has to
        # respond to where it is rather than eat channel 1 everywhere.
        for i in 1:G
            for r in 1:NRES
                res[r, i, j] *= keep
            end
        end
        for i in 1:G
            # Feed whichever channel the current rotation makes best-paying, so
            # yield and abundance move together and a fixed diet loses ground.
            # The offsets compose exactly as rot_yield composes them.
            zoff = (nz > 1 && p.rotate_inflow) ? zone_of(nz, Int32(i), G) - 1 : 0
            base = 1 + mod(-zoff, NRES)
            if cfrac > 0.0f0
                # Which channel pays best in the current state. Under rotation
                # that is the inverse shift; under independent permutations it is
                # whichever channel the permutation maps to yields[1].
                cue_src = if indep
                    best = 1
                    @inbounds for r in 1:NRES
                        p.yield_perm[r, coff + 1] == Int32(1) && (best = r; break)
                    end
                    best
                else
                    1 + mod(-(zoff + coff), NRES)
                end
                v = res[base, i, j] + infl * (1.0f0 - cfrac)
                res[base, i, j] = v > p.res_cap ? p.res_cap : v
                v = res[cue_src, i, j] + infl * cfrac
                res[cue_src, i, j] = v > p.res_cap ? p.res_cap : v
                continue
            end
            src = base
            v = res[src, i, j] + infl
            res[src, i, j] = v > p.res_cap ? p.res_cap : v
        end
    end
    return nothing
end

"""
    diffuse!(env, p)

Five-point Laplacian on every channel with no-flux (reflective) boundaries.
Threaded across channels; the interior is the SIMD fast path and the four
borders are fixed up afterwards.
"""
function diffuse!(env::Environment, p::Params)
    G = env.G
    src = env.res
    dst = env.buf
    cap = p.res_cap

    Threads.@threads for r in 1:NRES
        # larger byproducts diffuse more slowly than the primary substrate
        rate = p.diff_rate * p.diff_falloff^(r - 1)
        @inbounds for j in 2:(G - 1)
            @simd for i in 2:(G - 1)
                lap = src[r, i-1, j] + src[r, i+1, j] +
                      src[r, i, j-1] + src[r, i, j+1] - 4.0f0 * src[r, i, j]
                dst[r, i, j] = min(src[r, i, j] + rate * lap, cap)
            end
        end
        # no-flux borders: a missing neighbour is replaced by the centre cell,
        # so nothing leaks out of the dish.
        @inbounds for j in 1:G, i in 1:G
            (1 < i < G && 1 < j < G) && continue
            c = src[r, i, j]
            lap = (i > 1 ? src[r, i-1, j] : c) + (i < G ? src[r, i+1, j] : c) +
                  (j > 1 ? src[r, i, j-1] : c) + (j < G ? src[r, i, j+1] : c) - 4.0f0 * c
            dst[r, i, j] = min(c + rate * lap, cap)
        end
    end
    env.res, env.buf = env.buf, env.res
    return nothing
end

"""
    zone_of(p, x, G) -> 1:yield_zones

Which yield regime applies at column `x`. Vertical bands so that an organism's
zone changes as it moves, making the context genuinely something it has to
track rather than a fixed property of its lineage.
"""
@inline function zone_of(nzones::Int, x::Int32, G::Int)
    nzones <= 1 && return 1
    1 + (Int(x) - 1) * nzones ÷ G
end

"""
    zone_yield(yields, nzones, zone, r) -> Float32

Yield of channel `r` in `zone`. Each zone rotates the yield vector, so the
best-paying channel differs by zone while the SET of yields is identical
everywhere — the task changes, the total energy on offer does not.
"""
@inline function zone_yield(yields::Vector{Float32}, nzones::Int, zone::Int, r::Int)
    nzones <= 1 && return yields[r]
    @inbounds yields[1 + mod(r - 1 + (zone - 1), NRES)]
end

"""
    rot_yield(yields, off, r) -> Float32

Yield of channel `r` under a rotation of `off` steps. The zone and the abiotic
condition are both rotations of the same vector, so they compose by ADDING their
offsets rather than by multiplying yields — which would change the total energy
on offer and stop being the same manipulation.
"""
@inline function rot_yield(yields::Vector{Float32}, off::Int, r::Int)
    off == 0 && return @inbounds yields[r]
    @inbounds yields[1 + mod(r - 1 + off, NRES)]
end

"""
    condition_at(p, x, G, t) -> (cue, state)

The abiotic condition at column `x` and tick `t`. `cue` is the signal the genome
senses through `P_SENSE`, in [-1, 1]; `state` in 1:cond_states selects which
rotation of the yield vector applies, exactly as a zone did.

Unlike a resource field this is not consumed, so it stays legible no matter how
hard the population exploits the dish — which is the whole point, since the
zoned task failed precisely because consumers flatten the field they must read.
"""
@inline function condition_at(p::Params, x::Int32, G::Int, t::Int)
    ns = p.cond_states
    ns <= 1 && return (0.0f0, 1)
    phase = Float32(mod(t, p.cond_period)) / Float32(p.cond_period)
    if p.cond_gradient
        phase += Float32(Int(x) - 1) / Float32(G)
    end
    phase -= floor(phase)
    state = 1 + min(ns - 1, unsafe_trunc(Int, phase * ns))
    if p.cond_scramble
        # Drive the cue from an incommensurate period so it stays bounded and
        # smoothly varying -- the sensor sees a signal of the same shape and
        # magnitude -- but carries no usable information about `state`.
        # The temporal multiplier must be IRRATIONAL. An integer one is
        # commensurate with the cycle: 7.0 measured 0.000 bits WITH a spatial
        # gradient but leaked 0.067 bits without one, because with `cond_gradient
        # = false` the spatial term drops out and the two phases lock. 10/phi
        # (6.1803399 below) never locks, and measures 0.0004 bits without a gradient and 0.0014
        # with one. (A pure golden-ratio pair was tried earlier and leaked 0.72
        # bits, because cue and state then SHARE the spatial term.)
        q = Float32(mod(t, p.cond_period) * 6.1803399f0) / Float32(p.cond_period)
        if p.cond_gradient
            q += Float32(Int(x) - 1) * 5.0f0 / Float32(G)
        end
        q -= floor(q)
        return (p.cond_strength * (2.0f0 * q - 1.0f0), state)
    end
    # a signed cue, so a genome that ignores it is not merely mis-scaled
    (p.cond_strength * (2.0f0 * phase - 1.0f0), state)
end

total_resource(env::Environment, r::Int) = begin
    s = 0.0
    @inbounds for j in 1:env.G, i in 1:env.G
        s += env.res[r, i, j]
    end
    s
end
