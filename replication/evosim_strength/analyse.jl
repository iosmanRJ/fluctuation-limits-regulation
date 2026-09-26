#!/usr/bin/env julia
# =============================================================================
# The strength-vs-count test.
#
# EvoSim's headline connectivity result was measured as an EDGE COUNT
# (mean_nonzero_w) and a count-based density nzw/ng^2. The R-Aevol replication
# showed that counting links and weighing them can give OPPOSITE answers about
# the same networks.
#
# DEFINITIONS -- checked against the source, because getting them wrong is how
# the first version of this script produced meaningless columns:
#
#   genome_structure() (src/genome.jl) returns
#       nz   = number of weights with |w| > 1e-4
#       absw = sum(|w|) / ng^2        <- averaged over ALL ng^2 slots, zeros
#                                        included. It is already a DENSITY.
#
#   so, per run:
#       dens  = nzw / ng^2            count density
#       absw  = sum|w| / ng^2         strength density   (count density's twin)
#       totw  = absw * ng^2           total wiring strength = sum|w|
#       perw  = absw * ng^2 / nzw     mean strength per non-zero edge
#
#   R-Aevol for comparison (Stats_7.cpp): avg_influence = sum(factors)/nb_link,
#   i.e. per-LINK mean, so nb_link * avg_influence = total strength. The two
#   platforms' "total strength" and "per-edge strength" are therefore the
#   comparable quantities; their raw averages are not.
#
# One number per RUN (mean over the final 20% of ticks). Only runs with a DONE
# marker are read, so an in-progress run can never leak mid-run values into the
# tail window.
#
# The three FROZEN arms are bit-identical: with cond_states = 1 the cue is
# identically zero whatever its gain or scramble setting, so the same seed gives
# the same run. They are ONE arm, and only froz_none is used as the baseline.
# =============================================================================

using EvoSim, Printf, Statistics, Random, LinearAlgebra

Random.seed!(20260926)

const WINDOW_FRAC = 0.2

function summarise(dir::String)
    isfile(joinpath(dir, "stats.bin")) || return nothing
    S = read_stats(dir)
    size(S, 2) < 50 && return nothing
    ing, inz, iaw = stat_index("mean_ng"), stat_index("mean_nonzero_w"), stat_index("mean_abs_w")
    ip, ie = stat_index("pop"), stat_index("mean_energy")
    n = size(S, 2)
    w = max(10, round(Int, n * WINDOW_FRAC))
    keep = [c for c in (n - w + 1):n if S[ip, c] > 0]
    isempty(keep) && return nothing
    ng(c) = Float64(S[ing, c]); nz(c) = Float64(S[inz, c]); aw(c) = Float64(S[iaw, c])
    (; ng    = mean(ng(c) for c in keep),
       nzw   = mean(nz(c) for c in keep),
       dens  = mean(nz(c) / ng(c)^2 for c in keep),
       absw  = mean(aw(c) for c in keep),
       totw  = mean(aw(c) * ng(c)^2 for c in keep),
       perw  = mean(aw(c) * ng(c)^2 / nz(c) for c in keep if nz(c) > 0),
       pop   = mean(Float64(S[ip, c]) for c in keep),
       energy = mean(Float64(S[ie, c]) for c in keep))
end

msd(v) = (mean(v), length(v) < 2 ? 0.0 : std(v))

function permp(a, b; n = 100_000)
    obs = abs(mean(a) - mean(b)); pool = vcat(a, b); na = length(a); hits = 0
    for _ in 1:n
        shuffle!(pool)
        abs(mean(@view pool[1:na]) - mean(@view pool[na+1:end])) >= obs - 1e-12 && (hits += 1)
    end
    (hits + 1) / (n + 1)
end

function cohend(a, b)
    ma, sa = msd(a); mb, sb = msd(b); na, nb = length(a), length(b)
    sp = sqrt(((na - 1) * sa^2 + (nb - 1) * sb^2) / (na + nb - 2))
    sp > 0 ? (ma - mb) / sp : NaN
end

const MEASURES = [:ng, :nzw, :dens, :absw, :totw, :perw, :pop, :energy]

function contrast(name, arms, A, B)
    (haskey(arms, A) && haskey(arms, B)) || return
    println("\n", "="^92)
    @printf("%s\n    %s (n=%d)  vs  %s (n=%d)\n", name, A, length(arms[A]), B, length(arms[B]))
    println("="^92)
    @printf("%-7s %20s %20s %11s %7s %6s %9s %7s\n",
            "measure", A, B, "diff", "%", "t", "perm p", "d")
    for m in MEASURES
        a = Float64[getproperty(r, m) for r in arms[A]]
        b = Float64[getproperty(r, m) for r in arms[B]]
        ma, sa = msd(a); mb, sb = msd(b)
        se = sqrt(sa^2 / length(a) + sb^2 / length(b))
        @printf("%-7s %10.4f+-%-9.4f %10.4f+-%-9.4f %11.4f %6.1f%% %6.2f %9.5f %7.2f\n",
                m, ma, sa, mb, sb, ma - mb, 100 * (ma - mb) / mb,
                se > 0 ? (ma - mb) / se : NaN, permp(a, b), cohend(a, b))
    end
end

# ---------------------------------------------------------------------------
root = ARGS[1]
arms = Dict{String,Vector{Any}}()
for d in sort(readdir(root))
    full = joinpath(root, d)
    isdir(full) || continue
    isfile(full * ".DONE") || continue
    s = summarise(full)
    s === nothing && continue
    push!(get!(arms, join(split(d, "_s")[1:end-1], "_s"), Any[]), s)
end

println("ARMS (completed runs only)")
for a in sort(collect(keys(arms)))
    @printf("  %-12s n=%2d\n", a, length(arms[a]))
end

# the frozen arms must be identical -- verify rather than assume
if all(haskey(arms, a) for a in ("froz_none", "froz_info", "froz_yoked"))
    k = minimum(length(arms[a]) for a in ("froz_none", "froz_info", "froz_yoked"))
    same = all(arms["froz_none"][i].dens == arms["froz_info"][i].dens ==
               arms["froz_yoked"][i].dens for i in 1:k)
    println("\n  frozen arms bit-identical over first $k seeds: ", same)
end

println("""

  dens  count density   nzw/ng^2        absw  strength density  sum|w|/ng^2
  totw  total strength  sum|w|          perw  strength per non-zero edge
""")

contrast("C1   FLUCTUATION  (no cue)",                 arms, "fluc_none",  "froz_none")
contrast("C1b  FLUCTUATION  (informative cue)",        arms, "fluc_info",  "froz_none")
contrast("C1c  FLUCTUATION  (yoked cue)",              arms, "fluc_yoked", "froz_none")
contrast("C2   YOKED CONTROL: does information matter?", arms, "fluc_info", "fluc_yoked")
contrast("S3   cue present vs absent (fluctuating)",   arms, "fluc_info",  "fluc_none")

# ---------------------------------------------------------------------------
# Population confound. Fluctuating dishes carry fewer organisms, and population
# size drives genome architecture on its own (weaker purifying selection lets
# small near-neutral weights persist). Fluctuation and population are almost
# perfectly separated between arms, so an ANCOVA cannot pull them apart -- but
# the WITHIN-arm correlation says whether population can move these measures
# at all over the range it actually spans.
println("\n", "="^92)
println("POPULATION CONFOUND -- within-arm correlation with population size")
println("="^92)
@printf("%-12s %4s %14s %9s %9s %9s %9s\n", "arm", "n", "pop range", "r(dens)", "r(absw)", "r(perw)", "r(ng)")
for a in ["froz_none", "fluc_none", "fluc_info", "fluc_yoked"]
    haskey(arms, a) || continue
    R = arms[a]
    pop = Float64[r.pop for r in R]
    cr(m) = length(R) > 2 ? cor(pop, Float64[getproperty(r, m) for r in R]) : NaN
    @printf("%-12s %4d %6.0f-%-7.0f %9.3f %9.3f %9.3f %9.3f\n",
            a, length(R), minimum(pop), maximum(pop), cr(:dens), cr(:absw), cr(:perw), cr(:ng))
end
# pooled fluctuating arms: they share the environment, differ only in the cue,
# so pooling gives the widest population range with the treatment held fixed
fl = vcat([get(arms, a, Any[]) for a in ("fluc_none", "fluc_info", "fluc_yoked")]...)
if length(fl) > 3
    pop = Float64[r.pop for r in fl]
    @printf("%-12s %4d %6.0f-%-7.0f %9.3f %9.3f %9.3f %9.3f   (pooled fluctuating)\n",
            "fluc_all", length(fl), minimum(pop), maximum(pop),
            cor(pop, Float64[r.dens for r in fl]), cor(pop, Float64[r.absw for r in fl]),
            cor(pop, Float64[r.perw for r in fl]), cor(pop, Float64[r.ng for r in fl]))
end
