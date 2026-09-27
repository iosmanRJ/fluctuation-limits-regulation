#!/usr/bin/env julia
# =============================================================================
# P5 -- pre-registered EvoSim confirmation (PAPER.md section 11.5).
# Committed before any P5 run. Do not change measures, tests or rules after the
# confirmatory runs exist; anything not specified here is exploratory.
#
#   p5_analyse.jl <root>        root holds <arm>_s<seed>/ run dirs, with
#                               <arm>_s<seed>.err and <arm>_s<seed>.DONE beside them
#
# Arms: fluc_none (cond_states 6, cue gain 0) vs froz_none (cond_states 1, gain 0).
#
# Measures, per run, mean over the final 20% of logged ticks:
#   count density     nzw / ng^2        (per tick, then averaged)
#   strength density  mean_abs_w        (= sum|w| / ng^2, already a density)
#
# Predictions:
#   P5a  count density    higher under fluctuation
#   P5b  strength density lower  under fluctuation
#
# Test: two-sided permutation on the difference of means, 100,000 permutations,
# Holm across the two.
# Verdict: SUPPORTED if both Holm-significant in the predicted direction;
#          REFUTED if either is significant (unadjusted p < 0.05) in the
#          opposite direction; INCONCLUSIVE otherwise.
# Exclusions: a run is EXTINCT if its .err reports extinction or it has no
# living population in the final window; extinct runs are excluded and counted.
# If survival differs between arms by more than 10 percentage points, the
# result is reported as SURVIVORSHIP-CONFOUNDED and given no verdict.
#
# Secondary (no verdict): node count, total strength sum|w|, population.
# =============================================================================

using EvoSim, Printf, Statistics, Random

const NPERM = 100_000
const WINDOW_FRAC = 0.2

function extinct(dir)
    e = dir * ".err"
    isfile(e) && occursin("extinct", read(e, String))
end

function summarise(dir)
    S = read_stats(dir)
    n = size(S, 2)
    ing, inz, iaw, ip = stat_index("mean_ng"), stat_index("mean_nonzero_w"),
                        stat_index("mean_abs_w"), stat_index("pop")
    w = max(10, round(Int, n * WINDOW_FRAC))
    keep = [c for c in (n - w + 1):n if S[ip, c] > 0]
    isempty(keep) && return nothing
    g(c) = Float64(S[ing, c]); z(c) = Float64(S[inz, c]); a(c) = Float64(S[iaw, c])
    (; count_density    = mean(z(c) / g(c)^2 for c in keep),
       strength_density = mean(a(c) for c in keep),
       ng               = mean(g(c) for c in keep),
       total_strength   = mean(a(c) * g(c)^2 for c in keep),
       pop              = mean(Float64(S[ip, c]) for c in keep))
end

function permp(a, b, rng)
    obs = abs(mean(a) - mean(b)); pool = vcat(a, b); na = length(a); h = 0
    for _ in 1:NPERM
        shuffle!(rng, pool)
        abs(mean(@view pool[1:na]) - mean(@view pool[na+1:end])) >= obs - 1e-12 && (h += 1)
    end
    (h + 1) / (NPERM + 1)
end

function main(root)
    rng = MersenneTwister(20260927)
    runs = Dict("fluc_none" => Any[], "froz_none" => Any[])
    total = Dict("fluc_none" => 0, "froz_none" => 0)
    lost = Dict("fluc_none" => String[], "froz_none" => String[])
    for d in sort(readdir(root))
        full = joinpath(root, d)
        isdir(full) || continue
        isfile(full * ".DONE") || continue
        arm = join(split(d, "_s")[1:end-1], "_s")
        haskey(runs, arm) || continue
        total[arm] += 1
        s = extinct(full) ? nothing : summarise(full)
        if s === nothing
            push!(lost[arm], d); continue
        end
        push!(runs[arm], s)
    end

    println("SURVIVAL")
    surv = Dict{String,Float64}()
    for a in ("fluc_none", "froz_none")
        surv[a] = total[a] == 0 ? NaN : length(runs[a]) / total[a]
        @printf("  %-10s %d/%d survived (%.1f%%)%s\n", a, length(runs[a]), total[a],
                100surv[a], isempty(lost[a]) ? "" : "   extinct: " * join(lost[a], ", "))
    end
    gap = abs(surv["fluc_none"] - surv["froz_none"])
    @printf("  survival gap %.1f percentage points (limit 10)\n", 100gap)

    F, Z = runs["fluc_none"], runs["froz_none"]
    tests = [("P5a", :count_density, +1), ("P5b", :strength_density, -1)]
    res = []
    for (id, m, sign) in tests
        a = Float64[getproperty(r, m) for r in F]
        b = Float64[getproperty(r, m) for r in Z]
        push!(res, (id, m, sign, mean(a), mean(b), mean(a) - mean(b), permp(a, b, rng)))
    end
    ps = [r[7] for r in res]
    o = sortperm(ps); adj = similar(ps); run = 0.0
    for (k, i) in enumerate(o)
        run = max(run, min(1.0, (length(ps) - k + 1) * ps[i])); adj[i] = run
    end

    println("\nPRIMARY (fluctuating vs frozen)")
    for (k, (id, m, sign, ma, mb, d, p)) in enumerate(res)
        @printf("  %-4s %-17s %.4f vs %.4f  diff %+.4f (%+.1f%%)  p=%.5f  p_holm=%.5f\n",
                id, m, ma, mb, d, 100d / mb, p, adj[k])
    end

    print("\nVERDICT  ")
    if gap > 0.10
        println("SURVIVORSHIP-CONFOUNDED (no verdict)")
    else
        right = [r[6] * r[3] > 0 for r in res]
        sig = adj .< 0.05
        opp = [r[6] * r[3] < 0 && r[7] < 0.05 for r in res]
        println(any(opp) ? "REFUTED" : all(right .& sig) ? "SUPPORTED" : "INCONCLUSIVE")
    end

    println("\nSECONDARY (exploratory)")
    for m in (:ng, :total_strength, :pop)
        a = Float64[getproperty(r, m) for r in F]; b = Float64[getproperty(r, m) for r in Z]
        @printf("  %-15s %.4f vs %.4f  p=%.5f\n", m, mean(a), mean(b), permp(a, b, rng))
    end
end

main(ARGS[1])
