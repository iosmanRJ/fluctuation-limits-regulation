# Ticks per generation in the 180-run strength batch, from the slope of mean
# lineage depth over the second half of each run (unit = run).
using EvoSim, Statistics, Printf
root = "/Users/ivanosman/evosim_exp/main"
out = Dict{String,Vector{Float64}}()
for d in sort(readdir(root))
    full = joinpath(root, d); isdir(full) || continue
    isfile(full * ".DONE") || continue
    occursin("extinct", read(full * ".err", String)) && continue
    S = read_stats(full); n = size(S, 2)
    ig, id = stat_index("gen"), stat_index("mean_depth")
    a = n ÷ 2
    ticks = S[ig, n] - S[ig, a]; gens = S[id, n] - S[id, a]
    gens > 0 || continue
    push!(get!(out, join(split(d, "_s")[1:end-1], "_s"), Float64[]), ticks / gens)
end
for (k, v) in sort(collect(out))
    @printf("%-11s n=%2d  ticks/generation mean %.1f  (range %.1f-%.1f)\n", k, length(v), mean(v), minimum(v), maximum(v))
end
