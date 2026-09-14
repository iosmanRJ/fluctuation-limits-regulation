#!/usr/bin/env julia
# =============================================================================
# Tabulate a batch of runs, one row each, and summarise per group.
#
#   julia --project=. scripts/compare.jl runs/ab/A-s1 runs/ab/A-s2 ...
#   julia --project=. scripts/compare.jl runs/ab/*        # shell glob
#
# Groups are taken from the directory name up to the first '-', so
# runs/ab/A-s3 belongs to group "A". Reports extinction, the top-channel shares
# that measure branching, clade diversity and genome structure.
# =============================================================================

using EvoSim, Printf, Statistics

dirs = filter(isdir, ARGS)
isempty(dirs) && error("usage: compare.jl <rundir> [<rundir> ...]")

meta(d, k, default) = begin
    f = joinpath(d, "meta.txt")
    isfile(f) || return default
    for l in eachline(f)
        parts = split(l, '=', limit = 2)
        length(parts) == 2 && parts[1] == k && return parts[2]
    end
    default
end

group(d) = begin
    b = basename(rstrip(d, '/'))
    i = findfirst('-', b)
    i === nothing ? b : b[1:i-1]
end

rows = NamedTuple[]
for d in dirs
    isfile(joinpath(d, "stats.bin")) || continue
    M = read_stats(d)
    size(M, 2) == 0 && continue
    e = size(M, 1)   # guard against older runs with fewer stats
    last = size(M, 2)
    tick = M[stat_index("gen"), last]
    target = parse(Int, meta(d, "gens", "0"))
    pop = M[stat_index("pop"), last]
    # a run that stopped short of its target with pop 0 went extinct
    extinct = pop <= 0
    shares = e >= stat_index("top_share$(NRES)") ?
             [M[stat_index("top_share$(r)"), last] for r in 1:NRES] : fill(NaN, NRES)
    push!(rows, (; dir = d, grp = group(d), tick = Int(tick), target,
                 extinct, pop,
                 depth = M[stat_index("mean_depth"), last],
                 nonch1 = 100 * (1 - shares[1]),
                 shares,
                 shannon = M[stat_index("clade_shannon"), last],
                 nzw = e >= stat_index("mean_nonzero_w") ? M[stat_index("mean_nonzero_w"), last] : NaN,
                 waste = sum(M[stat_index("mean_phen$(r)"), last] for r in 2:NRES)))
end

isempty(rows) && error("no runs with stats found")

@printf("%-18s %9s %9s %8s %8s %8s %8s %8s\n",
        "run", "tick", "gens", "pop", "nonch1%", "waste", "shannon", "nz_w")
for r in sort(rows, by = x -> (x.grp, x.dir))
    if r.extinct
        @printf("%-18s %9d %9s %8s %8s %8s %8s %8s   EXTINCT\n",
                basename(r.dir), r.tick, "-", "-", "-", "-", "-", "-")
    else
        @printf("%-18s %9d %9.0f %8.0f %8.1f %8.3f %8.3f %8.1f\n",
                basename(r.dir), r.tick, r.depth, r.pop, r.nonch1, r.waste,
                r.shannon, r.nzw)
    end
end

println("\n", "="^92)
@printf("%-8s %6s %10s %22s %18s %10s\n",
        "group", "n", "extinct", "non-ch1 % (median, range)", "waste alloc", "nz weights")
for g in sort(unique(r.grp for r in rows))
    gr = filter(r -> r.grp == g, rows)
    surv = filter(r -> !r.extinct, gr)
    ne = count(r -> r.extinct, gr)
    if isempty(surv)
        @printf("%-8s %6d %10s %22s %18s %10s\n", g, length(gr), "$ne/$(length(gr))", "-", "-", "-")
        continue
    end
    nc = [r.nonch1 for r in surv]
    wa = [r.waste for r in surv]
    nz = [r.nzw for r in surv]
    @printf("%-8s %6d %10s %8.1f  (%.1f-%.1f) %12.3f %10.1f\n",
            g, length(gr), "$ne/$(length(gr))", median(nc), minimum(nc), maximum(nc),
            median(wa), median(nz))
end

# how many runs in each group produced a substantial second ecotype?
println()
for g in sort(unique(r.grp for r in rows))
    surv = filter(r -> r.grp == g && !r.extinct, rows)
    isempty(surv) && continue
    for thr in (5.0, 10.0)
        n = count(r -> r.nonch1 >= thr, surv)
        @printf("group %-4s: %2d/%2d surviving runs have >=%.0f%% non-ch1 specialists\n",
                g, n, length(surv), thr)
    end
end
