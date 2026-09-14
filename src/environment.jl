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
    infl = p.inflow * season(p, gen)
    if env.famine > 0
        infl *= p.famine_factor
        env.famine -= 1
    end
    infl = max(infl, 0.0f0)
    keep = 1.0f0 - p.decay
    res = env.res
    @inbounds for j in 1:G
        @simd for i in 1:G
            res[1, i, j] = min(res[1, i, j] * keep + infl, p.res_cap)
        end
        for r in 2:NRES
            @simd for i in 1:G
                res[r, i, j] *= keep
            end
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

total_resource(env::Environment, r::Int) = begin
    s = 0.0
    @inbounds for j in 1:env.G, i in 1:env.G
        s += env.res[r, i, j]
    end
    s
end
