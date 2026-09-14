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
const NPHEN = NRES + 6            # 12 phenotype outputs

# --- phenotype layout --------------------------------------------------------
const P_UPTAKE  = 1:NRES          # 1:6  allocation across resource channels
const P_SIZE    = NRES + 1        # 7    body size -> maintenance cost
const P_MOTIL   = NRES + 2        # 8    movement probability
const P_DIVIDE  = NRES + 3        # 9    division threshold
const P_REPAIR  = NRES + 4        # 10   DNA repair fidelity (heritable, costly)
const P_STRESS  = NRES + 5        # 11   tolerance to shocks AND to waste toxicity
const P_LEAK    = NRES + 6        # 12   fraction of intake released as byproduct

# NOTE: there is deliberately no evolvable "which byproduct do I emit" trait.
# Where a byproduct goes is chemistry (see the conversion ladder in sim.jl), not
# a decision. When it WAS evolvable, organisms discovered a futile cycle — they
# evolved to excrete 36% of their waste straight back into the channel they had
# just eaten (vs 17% for a neutral profile), which meant the waste channels were
# never fed and cross-feeding had nothing to run on. Excluding self-conversion
# is also what the microbial consumer-resource literature does: the metabolic
# matrix D has its diagonal excluded on thermodynamic grounds.

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
                  env::Vector{Float32}; niter::Int = 8, plasticity::Float32 = 1.0f0)
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
                s += plasticity * env[k]
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

    # Uptake affinity is an ALLOCATION, not a set of independent dials: the
    # organism divides one fixed budget across channels. Independent affinities
    # make breadth merely expensive, which collapses to "specialise on the best
    # channel"; a budget makes breadth genuinely costly in opportunity terms, so
    # an under-exploited channel is always worth something to somebody.
    @inbounds normalise!(phen, P_UPTAKE)
    return nothing
end

"""
    inherit!(W, bias, child, parent, rng, p, repair) -> (nmut, big_effect, dload)

Copy parent's genome to the child slot and mutate it. Returns how many
mutations landed, whether any was large enough to count as a clade-founding
event, and the net change in deleterious load.

Mutation classes:
  * point      - perturb one regulatory weight (the overwhelming majority)
  * silencing  - zero every outgoing weight of a node (gene knockout)
  * duplication- copy one node's regulatory role onto another (raw material
                 for novelty; this is how GRNs gain functions)

Point mutations perturb the network, and a population near a local optimum has
more ways to get worse than better — that part is emergent. On top of it each
mutation is classified deleterious/neutral/beneficial so an explicit load
counter can carry a direct metabolic penalty; because load can fall as well as
rise, Muller's ratchet behaves like a ratchet rather than a countdown.
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
    dload = 0

    # expected point mutations this birth; sampled as a small Poisson draw
    lam = rate * nsites
    npoint = pois_small(rng, lam)
    @inbounds for _ in 1:npoint
        m = rand(rng, 1:NG)
        k = rand(rng, 1:NG)
        d = p.mut_sigma * randn(rng, Float32)
        W[m, k, i] += d
        nmut += 1
        dload += classify_mutation(rng, p)
        if abs(d) > 2.0f0 * p.mut_sigma
            big = true
        end
    end

    # bias mutations
    if rand(rng, Float32) < rate * NG
        k = rand(rng, 1:NG)
        bias[k, i] += p.mut_sigma * randn(rng, Float32)
        nmut += 1
        dload += classify_mutation(rng, p)
    end

    # structural: gene silencing
    if rand(rng, Float32) < p.p_silence
        k = rand(rng, 1:NG)
        @inbounds for m in 1:NG
            W[k, m, i] = 0.0f0
        end
        nmut += 1
        dload += 1          # knocking out a gene is reliably damaging
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

    return nmut, big, dload
end

"""Deleterious (+1), beneficial (-1) or neutral (0)."""
@inline function classify_mutation(rng, p)
    u = rand(rng, Float32)
    u < p.p_deleterious && return 1
    u < p.p_deleterious + p.p_beneficial && return -1
    return 0
end

@inline function normalise!(phen::Vector{Float32}, rng_idx)
    tot = 0.0f0
    @inbounds for k in rng_idx
        tot += phen[k]
    end
    @inbounds if tot > 1.0f-6
        for k in rng_idx
            phen[k] /= tot
        end
    else
        for k in rng_idx
            phen[k] = 1.0f0 / length(rng_idx)
        end
    end
    return nothing
end

"""
    seed_fast!(rng, seed)

Re-seed an `Xoshiro` in place, allocation-free.

`Random.seed!(::Xoshiro, ::UInt64)` allocates 368 bytes per call, which at tens
of thousands of births per tick is pure garbage. This does the same SplitMix64
state expansion by assigning the state words directly.
"""
@inline function splitmix64(z::UInt64)
    x = z
    x = (x ⊻ (x >> 30)) * 0xbf58476d1ce4e5b9
    x = (x ⊻ (x >> 27)) * 0x94d049bb133111eb
    return x ⊻ (x >> 31)
end

@inline function seed_fast!(rng::Xoshiro, seed::UInt64)
    # Written out rather than with an inner helper closure: a closure that
    # mutates a captured variable gets boxed, which allocates more than the
    # Random.seed! call this replaces.
    z = seed
    z += 0x9e3779b97f4a7c15; a = splitmix64(z)
    z += 0x9e3779b97f4a7c15; b = splitmix64(z)
    z += 0x9e3779b97f4a7c15; c = splitmix64(z)
    z += 0x9e3779b97f4a7c15; d = splitmix64(z)
    z += 0x9e3779b97f4a7c15; e = splitmix64(z)
    # a zero state is absorbing for xoshiro; vanishingly unlikely here, but one
    # occurrence would be silent, so guard it
    if (a | b | c | d | e) == 0x0000000000000000
        a = 0x9e3779b97f4a7c15
    end
    setfield!(rng, :s0, a); setfield!(rng, :s1, b); setfield!(rng, :s2, c)
    setfield!(rng, :s3, d); setfield!(rng, :s4, e)
    return rng
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

"""
    channel_cost(phen, exp) -> Float32

Proteome overhead of the uptake allocation: `sum(a_r ^ exp)`. Computed once at
birth and cached, because the allocation never changes afterwards and a `pow`
per channel per organism per tick would be a real cost at scale.

With `exp < 1` this is concave, so spreading the same total allocation across
more channels costs more — which is what makes the interior of the allocation
simplex a fitness minimum rather than a flat ridge.
"""
@inline function channel_cost(phen::AbstractVector{Float32}, ex::Float32)
    c = 0.0f0
    @inbounds if ex == 0.5f0          # by far the common case; avoids a pow
        for r in P_UPTAKE
            c += sqrt(phen[r])
        end
    else
        for r in P_UPTAKE
            c += phen[r]^ex
        end
    end
    c
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
