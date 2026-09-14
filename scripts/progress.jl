#!/usr/bin/env julia
# Report how far a run has got, without disturbing it.
#   julia --project=. scripts/progress.jl runs/flagship
using EvoSim, Printf, Dates

dir = length(ARGS) >= 1 ? ARGS[1] : error("usage: progress.jl <rundir>")
meta = Dict(split(l, '=', limit = 2)[1] => split(l, '=', limit = 2)[2]
            for l in eachline(joinpath(dir, "meta.txt")) if occursin('=', l))
target = parse(Int, meta["gens"])
every  = parse(Int, meta["stat_every"])

M = read_stats(dir)
n = size(M, 2)
n == 0 && (println("no records flushed yet"); exit(0))

tick  = Int(M[stat_index("gen"), end])
depth = M[stat_index("mean_depth"), end]
@printf("%s\n", dir)
@printf("  tick        %d / %d  (%.2f%%)\n", tick, target, 100tick/target)
@printf("  generations %.0f  (%.1f ticks/gen)\n", depth, tick / max(depth, 1))
@printf("  population  %.0f\n", M[stat_index("pop"), end])
@printf("  clades      %.0f alive, Shannon %.3f\n",
        M[stat_index("clades_alive"), end], M[stat_index("clade_shannon"), end])
@printf("  load        %.3f\n", M[stat_index("mean_load"), end])
@printf("  records     %d (every %d ticks), %.1f MB\n", n, every,
        filesize(joinpath(dir, "stats.bin")) / 1e6)
snaps = filter(startswith("snap_"), readdir(dir))
@printf("  snapshots   %d\n", length(snaps))

mt = mtime(joinpath(dir, "stats.bin"))
age = time() - mt
@printf("  last write  %.0fs ago%s\n", age, age > 600 ? "  ** stalled? **" : "")
