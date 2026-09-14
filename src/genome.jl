# =============================================================================
# Gene-regulatory-network genome.
#
# A genome is a weight matrix W (node -> node) plus per-node bias. Development
# runs the network to an approximate fixed point ONCE at birth; the resulting
# expression vector is read off as the phenotype and CACHED for life. That is
# ~15x cheaper than re-expressing every tick and is the biologically sensible
# model anyway: develop, then live as an adult.
#
# W is stored (from, to, organism) so the inner product loop is contiguous.
# =============================================================================

const NG    = 20                  # regulatory nodes
const NRES  = 6                   # resource channels
const NPHEN = 2 * NRES + 5        # 17 phenotype outputs

# --- phenotype layout --------------------------------------------------------
const P_UPTAKE  = 1:NRES              # 1:6   affinity for each resource channel
const P_EXCRETE = (NRES + 1):(2NRES)  # 7:12  which wastes this organism emits
const P_SIZE    = 2NRES + 1           # 13    body size -> maintenance cost
const P_MOTIL   = 2NRES + 2           # 14    movement probability
const P_DIVIDE  = 2NRES + 3           # 15    division threshold
const P_REPAIR  = 2NRES + 4           # 16    DNA repair fidelity (heritable, costly)
const P_STRESS  = 2NRES + 5           # 17    tolerance to environmental shocks

"""
    develop!(phen, e, buf, W, bias, i, env; niter)

Run organism `i`'s regulatory network to an approximate steady state and write
its phenotype into `phen`. `env` holds local resource concentrations, so the
same genome can express differently in different conditions (plasticity).

`e` and `buf` are caller-owned scratch vectors of length `NG`; this function
allocates nothing.
"""
function develop!(phen::Vector{Float32}, e::Vector{Float32}, buf::Vector{Float32},
                  W::Array{Float32,3}, bias::Matrix{Float32}, i::Int,
                  env::Vector{Float32}; niter::Int = 8)
    @inbounds for k in 1:NG
        e[k] = 0.1f0
    end
    @inbounds for _ in 1:niter
        for k in 1:NG
            s = bias[k, i]
            @simd for m in 1:NG
                s += W[m, k, i] * e[m]      # contiguous in m
            end
            # resource concentrations feed the first NRES nodes -> plasticity
            if k <= NRES
                s += env[k]
            end
            buf[k] = tanh(s)
        end
        for k in 1:NG
            e[k] = buf[k]
        end
    end

    # map expression (-1,1) -> phenotype (0,1)
    @inbounds for k in 1:NPHEN
        phen[k] = 0.5f0 * (e[k] + 1.0f0)
    end

    # excretion profile is a proportion, so normalise it
    @inbounds begin
        tot = 0.0f0
        for k in P_EXCRETE
            tot += phen[k]
        end
        if tot > 1.0f-6
            for k in P_EXCRETE
                phen[k] /= tot
            end
        else
            for k in P_EXCRETE
                phen[k] = 1.0f0 / NRES
            end
        end
    end
    return nothing
end

"""
    inherit!(W, bias, child, parent, rng, p) -> (nmut, big_effect)

Copy parent's genome to the child slot and mutate it. Returns how many
mutations landed and whether any was large enough to count as a clade-founding
event.

Mutation classes:
  * point      - perturb one regulatory weight (the overwhelming majority)
  * silencing  - zero every outgoing weight of a node (gene knockout)
  * duplication- copy one node's regulatory role onto another (raw material
                 for novelty; this is how GRNs gain functions)

Most point mutations are deleterious simply because a population sitting near
a local optimum has more ways to get worse than better. That emergent load is
the honest version; `nmut` tracks accumulation so Muller's ratchet is visible.
"""
function inherit!(W::Array{Float32,3}, bias::Matrix{Float32},
                  child::Int, parent::Int, rng, p, repair::Float32)
    @inbounds for k in 1:NG, m in 1:NG
        W[m, k, child] = W[m, k, parent]
    end
    @inbounds for k in 1:NG
        bias[k, child] = bias[k, parent]
    end
    return mutate!(W, bias, child, rng, p, repair)
end

function mutate!(W::Array{Float32,3}, bias::Matrix{Float32}, i::Int, rng, p, repair::Float32)
    # repair fidelity is itself heritable: high fidelity suppresses the rate
    rate = p.mut_rate * (1.0f0 - p.repair_strength * repair)
    rate = max(rate, 1.0f-6)

    nsites = NG * NG
    nmut = 0
    big  = false

    # expected point mutations this birth; sampled as a small Poisson draw
    lam = rate * nsites
    npoint = pois_small(rng, lam)
    @inbounds for _ in 1:npoint
        m = rand(rng, 1:NG)
        k = rand(rng, 1:NG)
        d = p.mut_sigma * randn(rng, Float32)
        W[m, k, i] += d
        nmut += 1
        if abs(d) > 2.0f0 * p.mut_sigma
            big = true
        end
    end

    # bias mutations
    if rand(rng, Float32) < rate * NG
        k = rand(rng, 1:NG)
        bias[k, i] += p.mut_sigma * randn(rng, Float32)
        nmut += 1
    end

    # structural: gene silencing
    if rand(rng, Float32) < p.p_silence
        k = rand(rng, 1:NG)
        @inbounds for m in 1:NG
            W[k, m, i] = 0.0f0
        end
        nmut += 1
        big = true
    end

    # structural: gene duplication (donor's regulatory role copied onto target)
    if rand(rng, Float32) < p.p_duplicate
        src = rand(rng, 1:NG)
        dst = rand(rng, 1:NG)
        if src != dst
            @inbounds for m in 1:NG
                W[m, dst, i] = W[m, src, i] + 0.05f0 * randn(rng, Float32)
                W[dst, m, i] = W[src, m, i] + 0.05f0 * randn(rng, Float32)
            end
            bias[dst, i] = bias[src, i]
            nmut += 1
            big = true
        end
    end

    return nmut, big
end

"""Poisson draw, fine for the small lambda (<~10) we always use here."""
@inline function pois_small(rng, lam::Float32)
    lam <= 0.0f0 && return 0
    L = exp(-lam)
    k = 0
    pr = 1.0f0
    while true
        pr *= rand(rng, Float32)
        pr <= L && return k
        k += 1
        k > 64 && return k          # guard against pathological lambda
    end
end

"""Random founder genome: sparse, weakly connected, mostly quiet."""
function seed_genome!(W::Array{Float32,3}, bias::Matrix{Float32}, i::Int, rng, density = 0.25f0)
    @inbounds for k in 1:NG, m in 1:NG
        W[m, k, i] = rand(rng, Float32) < density ? 0.6f0 * randn(rng, Float32) : 0.0f0
    end
    @inbounds for k in 1:NG
        bias[k, i] = 0.2f0 * randn(rng, Float32)
    end
    return nothing
end
