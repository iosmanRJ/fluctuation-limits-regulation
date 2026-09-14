#!/usr/bin/env julia
# =============================================================================
# Headless runner.
#
#   julia --project=. -t auto scripts/run.jl --gens 1000000 --out runs/big
#
# Any Params field can be overridden with --<field> <value>:
#   --grid 512 --capacity 200000 --mut-rate 2e-3 --inflow 0.3
# (dashes in the flag map to underscores in the field name)
# =============================================================================

using EvoSim, Dates, Printf

function parse_args(args)
    opts = Dict{Symbol,Any}()
    gens = 100_000
    out  = joinpath("runs", Dates.format(now(), "yyyymmdd-HHMMSS"))
    snapshots = true
    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--gens"
            gens = parse(Int, args[i+1]); i += 2
        elseif a == "--out"
            out = args[i+1]; i += 2
        elseif a == "--no-snapshots"
            snapshots = false; i += 1
        elseif a == "--help" || a == "-h"
            print_help(); exit(0)
        elseif startswith(a, "--")
            key = Symbol(replace(a[3:end], '-' => '_'))
            hasfield(Params, key) || error("unknown option $a (no Params field `$key`)")
            T = fieldtype(Params, key)
            opts[key] = T <: AbstractVector ?
                        parse.(eltype(T), split(args[i+1], ',')) :
                        parse(T, args[i+1])
            i += 2
        else
            error("unexpected argument: $a")
        end
    end
    (gens, out, snapshots, opts)
end

function print_help()
    println("""
    usage: run.jl [--gens N] [--out DIR] [--no-snapshots] [--<param> VALUE ...]

      --gens N          ticks to simulate (default 100000)
      --out DIR         output directory (default runs/<timestamp>)
      --no-snapshots    skip spatial snapshots, stats only

    Any Params field may be set with --<field> <value>, dashes for underscores:
    """)
    for f in fieldnames(Params)
        @printf("      --%-22s %s\n", replace(String(f), '_' => '-'),
                fieldtype(Params, f))
    end
    println("""
    Note: ticks are not generations. One generation of descent takes roughly 15
    ticks at default settings; the run reports mean lineage depth so you never
    have to guess. Chart a finished run with:

      julia --project=. scripts/dashboard.jl <outdir>
    """)
end

gens, out, snapshots, opts = parse_args(ARGS)
p = Params(; opts...)

@printf("EvoSim — %d ticks, grid %d, capacity %d, %d threads\n",
        gens, p.grid, p.capacity, Threads.nthreads())
isempty(opts) || println("overrides: ", join(("$k=$v" for (k, v) in opts), ", "))
println("output: ", out)

r = run_sim(; gens, p, outdir = out, snapshots, verbose = true)

live = [i for i in 1:high_water(r.sim.pop) if r.sim.pop.alive[i]]
if isempty(live)
    @printf("\nEXTINCT at tick %d\n", r.extinct_at)
else
    depth = sum(r.sim.pop.depth[i] for i in live) / length(live)
    @printf("\nfinal: %d organisms · %.0f generations of descent (%.1f ticks/gen) · %d clades\n",
            length(live), depth, r.gens_run / max(depth, 1), nclades(r.sim.clades))
end
println("\nchart it:  julia --project=. scripts/dashboard.jl $out")
