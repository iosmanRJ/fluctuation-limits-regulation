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
        zeros(Int32, cap), zeros(Int32, cap), zeros(UInt32, cap),
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
    # Founder phenotypes kept FLAT (NPHEN values per clade) rather than as a
    # vector of vectors: one allocation per clade founding is ruinous when
    # clades are founded hundreds of thousands of times in a run.
    founder::Vector{Float32}
end

CladeTable() = CladeTable(UInt32[0], Int32[0], zeros(Float32, NPHEN))

nclades(ct::CladeTable) = length(ct.parent)

@inline founder_offset(c::UInt32) = Int(c) * NPHEN

function new_clade!(ct::CladeTable, parent::UInt32, gen::Int, phen::AbstractVector{Float32})
    push!(ct.parent, parent)
    push!(ct.born, Int32(gen))
    @inbounds for k in 1:NPHEN
        push!(ct.founder, phen[k])
    end
    return UInt32(length(ct.parent) - 1)
end

"""
    phen_dist2(ct, clade, phen) -> squared Euclidean distance

How far a phenotype has drifted from its clade's founder. A clade is defined as
a phenotypic neighbourhood, which is both more meaningful than counting
"large-effect" mutations and self-limiting: drift has to actually go somewhere
to found a new clade, so the table cannot explode.
"""
@inline function phen_dist2(ct::CladeTable, clade::UInt32, phen::AbstractVector{Float32})
    off = founder_offset(clade)
    d = 0.0f0
    @inbounds for k in 1:NPHEN
        dk = phen[k] - ct.founder[off + k]
        d += dk * dk
    end
    d
end

"""
    prune_clades!(ct, pop) -> nkept

Garbage-collect the clade table. A clade is kept if it has living members or is
an ancestor of one; everything else is unreachable history. Without this the
table grows unboundedly (measured: ~40 new clades per tick, which is several GB
over a million-tick run) purely from lineages that flicker across the radius
boundary and die out.

Returns the number of clades kept. Organism clade ids are remapped in place.
"""
function prune_clades!(ct::CladeTable, pop::Population)
    n = nclades(ct)
    keep = falses(n)
    keep[1] = true                       # the root founder clade always survives

    @inbounds for i in 1:pop.watermark
        pop.alive[i] || continue
        keep[Int(pop.clade[i]) + 1] = true
    end

    # walk each live clade up to the root; parents always precede children in
    # the table, so a single reverse pass marks every ancestor
    @inbounds for c in n:-1:2
        keep[c] && (keep[Int(ct.parent[c]) + 1] = true)
    end

    remap = zeros(UInt32, n)
    nk = 0
    @inbounds for c in 1:n
        keep[c] || continue
        nk += 1
        remap[c] = UInt32(nk - 1)
    end
    nk == n && return nk                 # nothing to do

    new_parent = Vector{UInt32}(undef, nk)
    new_born   = Vector{Int32}(undef, nk)
    new_found  = Vector{Float32}(undef, nk * NPHEN)
    j = 0
    @inbounds for c in 1:n
        keep[c] || continue
        j += 1
        new_parent[j] = remap[Int(ct.parent[c]) + 1]
        new_born[j]   = ct.born[c]
        src = (c - 1) * NPHEN
        dst = (j - 1) * NPHEN
        for k in 1:NPHEN
            new_found[dst + k] = ct.founder[src + k]
        end
    end

    resize!(ct.parent, nk); copyto!(ct.parent, new_parent)
    resize!(ct.born, nk);   copyto!(ct.born, new_born)
    resize!(ct.founder, nk * NPHEN); copyto!(ct.founder, new_found)

    @inbounds for i in 1:pop.watermark
        pop.alive[i] || continue
        pop.clade[i] = remap[Int(pop.clade[i]) + 1]
    end
    return nk
end
