# One row per EvoSim run of the strength-vs-count batch (paper/data/evosim_runs.csv),
# summarised as replication/evosim_strength/analyse.jl does: mean over the final
# 20% of logged ticks, per-tick values averaged. Extinct = .err reports extinction.
#
# Each organism's number of non-zero weights (|w| > 1e-4), its Σ|w|/ng² and its ng
# are computed at birth; the log holds their population means at each tick. So
#   edges            = mean number of non-zero weights
#   count_density    = edges / (mean ng)²
#   strength_density = mean Σ|w|/ng²
#   total_strength   = strength_density * (mean ng)²
#   per_edge         = total_strength / edges
# Runs on the Mac mini, where the run directories live.
using EvoSim, Statistics, Printf

root = "/Users/ivanosman/evosim_exp/main"
io = open("/Users/ivanosman/evosim_exp/evosim_runs.csv", "w")
println(io, "arm,seed,extinct,count_density,strength_density,total_strength,per_edge,ng,pop,edges")
for d in sort(readdir(root))
    full = joinpath(root, d)
    isdir(full) || continue
    isfile(full * ".DONE") || continue
    parts = split(d, "_s")
    arm = join(parts[1:end-1], "_s"); seed = parts[end]
    ext = isfile(full * ".err") && occursin("extinct", read(full * ".err", String))
    S = read_stats(full); n = size(S, 2)
    ing, inz, iaw, ip = stat_index("mean_ng"), stat_index("mean_nonzero_w"),
                        stat_index("mean_abs_w"), stat_index("pop")
    w = max(10, round(Int, n * 0.2))
    keep = [c for c in (n - w + 1):n if S[ip, c] > 0]
    if isempty(keep)
        @printf(io, "%s,%s,%s,NaN,NaN,NaN,NaN,NaN,NaN,NaN\n", arm, seed, ext); continue
    end
    g(c) = Float64(S[ing, c]); z(c) = Float64(S[inz, c]); a(c) = Float64(S[iaw, c])
    @printf(io, "%s,%s,%s,%.6f,%.6f,%.6f,%.6f,%.4f,%.1f,%.4f\n", arm, seed, ext,
            mean(z(c) / g(c)^2 for c in keep), mean(a(c) for c in keep),
            mean(a(c) * g(c)^2 for c in keep),
            mean(a(c) * g(c)^2 / z(c) for c in keep if z(c) > 0),
            mean(g(c) for c in keep), mean(Float64(S[ip, c]) for c in keep),
            mean(z(c) for c in keep))
end
close(io)
println("wrote evosim_runs.csv")
