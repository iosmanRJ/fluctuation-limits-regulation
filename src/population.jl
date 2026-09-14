# =============================================================================
# Population in struct-of-arrays layout with a free-slot stack.
#
# SoA rather than Vector{Organism} because the hot loop touches a few fields of
# every organism: AoS would drag whole 1.6 KB genomes through cache for the
# sake of reading one float. Genomes live in their own array and are only
# touched at birth.
# =============================================================================

mutable struct Population
    # --- state touched every tick (keep these compact) ------------------------
    x::Vector{Int32}
    y::Vector{Int32}
    energy::Vector{Float32}
    age::Vector{Int32}
    alive::Vector{Bool}
    phen::Matrix{Float32}        # NPHEN x cap   cached at birth

    # --- genome, touched only at birth ---------------------------------------
    W::Array{Float32,3}          # NG x NG x cap
    bias::Matrix{Float32}        # NG x cap

    # --- bookkeeping ----------------------------------------------------------
    depth::Vector{Int32}         # lineage depth = generations of descent from a founder
    load::Vector{Int32}          # accumulated mutations in this lineage (Muller's ratchet)
    bigmut::Vector{Int32}        # big-effect mutations since this clade was founded
    clade::Vector{UInt32}
    occ::Matrix{Int32}           # per-cell occupancy, for local competition
    free::Vector{Int32}          # stack of reusable slots
    n::Int                       # number currently alive
    watermark::Int               # highest slot ever handed out
    capacity::Int
end

function Population(p::Params)
    cap = p.capacity
    Population(
        zeros(Int32, cap), zeros(Int32, cap),
        zeros(Float32, cap), zeros(Int32, cap), fill(false, cap),
        zeros(Float32, NPHEN, cap),
        zeros(Float32, NG, NG, cap), zeros(Float32, NG, cap),
        zeros(Int32, cap), zeros(Int32, cap), zeros(Int32, cap), zeros(UInt32, cap),
        zeros(Int32, p.grid, p.grid),
        Int32[], 0, 0, cap,
    )
end

"""
    claim!(pop) -> slot or 0

Reuse a freed slot if one exists, otherwise extend the watermark. Returns 0
when the dish is at capacity. O(1) — never scan for a free slot, at 1e6
generations that cost would dominate everything else.
"""
@inline function claim!(pop::Population)
    i = 0
    if !isempty(pop.free)
        i = Int(pop!(pop.free))
    elseif pop.watermark < pop.capacity
        pop.watermark += 1
        i = pop.watermark
    else
        return 0                  # dish is full
    end
    @inbounds pop.alive[i] = true
    pop.n += 1
    return i
end

@inline function kill!(pop::Population, i::Int)
    @inbounds begin
        pop.alive[i] = false
        pop.occ[pop.x[i], pop.y[i]] -= Int32(1)
    end
    push!(pop.free, Int32(i))
    pop.n -= 1
    return nothing
end

"""Move organism `i`, keeping the occupancy grid in step."""
@inline function place!(pop::Population, i::Int, nx::Int32, ny::Int32)
    @inbounds begin
        pop.occ[pop.x[i], pop.y[i]] -= Int32(1)
        pop.x[i] = nx
        pop.y[i] = ny
        pop.occ[nx, ny] += Int32(1)
    end
    return nothing
end

"""Upper bound for iteration: slots above this were never occupied."""
@inline high_water(pop::Population) = pop.watermark

# =============================================================================
# Phylogeny: a full parent-of-every-organism tree is impossible over 1e6
# generations. Instead organisms carry a clade id that only changes when enough
# large-effect mutations have accumulated in a lineage. That yields a coarse
# but genuine tree of manageable size.
# =============================================================================

struct CladeTable
    parent::Vector{UInt32}
    born::Vector{Int32}
    founder_phen::Vector{Vector{Float32}}
end

CladeTable() = CladeTable(UInt32[0], Int32[0], [zeros(Float32, NPHEN)])

nclades(ct::CladeTable) = length(ct.parent)

function new_clade!(ct::CladeTable, parent::UInt32, gen::Int, phen::AbstractVector{Float32})
    push!(ct.parent, parent)
    push!(ct.born, Int32(gen))
    push!(ct.founder_phen, copy(Vector{Float32}(phen)))
    return UInt32(length(ct.parent) - 1)
end
