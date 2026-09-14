# =============================================================================
# Headless driver.
# =============================================================================

"""
    run_sim(; gens, p, outdir, stat_every, verbose) -> RunResult

Run `gens` generations, writing per-generation stats and log-spaced snapshots
into `outdir`. Stops early if the population goes extinct — with shocks and
mutation load enabled that is a real outcome, not a bug, so it is reported
rather than hidden.
"""
function run_sim(; gens::Int,
                   p::Params = Params(),
                   outdir::String = joinpath("runs", Dates.format(now(), "yyyymmdd-HHMMSS")),
                   stat_every::Int = max(1, gens ÷ 200_000),
                   snapshots::Bool = true,
                   verbose::Bool = true)

    mkpath(outdir)
    sim = Sim(p)
    slog = StatLog(outdir)
    sched = snapshots ? Set(snapshot_schedule(gens, p.snapshot_decades)) : Set{Int}()

    write_meta(outdir, p, gens, stat_every)

    t0 = time()
    extinct_at = 0
    last_report = t0
    warned_full = false

    for g in 1:gens
        step!(sim)

        if sim.pop.n == 0
            extinct_at = g
            record!(slog, sim)
            verbose && @warn "population went extinct" generation=g
            break
        end

        if !warned_full && sim.dropped > 0
            warned_full = true
            @warn """organism array is full — births are being dropped.
                   Capacity is an implementation limit, not an ecological one; while
                   it binds, reproductive success stops being selected on and every
                   trait drifts toward whatever is cheapest. Raise `capacity` (or
                   `capacity_per_cell`) and rerun.""" capacity=sim.pop.capacity population=sim.pop.n generation=g
        end

        (g % stat_every == 0) && record!(slog, sim)
        (g in sched) && write_snapshot(sim, outdir)

        if verbose && time() - last_report > 5.0
            el = time() - t0
            eta = el * (gens - g) / g
            @printf("[%s] tick %9d/%d  pop %8d  clades %6d  %.0f tick/s  ETA %s\n",
                    Dates.format(now(), "HH:MM:SS"), g, gens, sim.pop.n,
                    nclades(sim.clades), g / el, fmt_dur(eta))
            # stdout is block-buffered when redirected to a file, so a long
            # headless run would show no progress at all without this.
            flush(stdout)
            last_report = time()
        end
    end

    flush_stats!(slog)
    close(slog)
    write_clades(outdir, sim.clades)
    write_shocks(outdir, sim.shock_log)

    elapsed = time() - t0
    if verbose
        @printf("done: %d ticks in %s (%.0f tick/s)%s\n",
                sim.gen, fmt_dur(elapsed), sim.gen / elapsed,
                extinct_at > 0 ? "  [EXTINCT at $extinct_at]" : "")
        flush(stdout)
    end

    if sim.dropped > 0
        @warn "run finished with dropped births; treat its evolutionary results as unreliable" dropped=sim.dropped
    end
    return (; sim, outdir, elapsed, extinct_at, gens_run = sim.gen)
end

function fmt_dur(s::Real)
    s < 60  && return @sprintf("%.0fs", s)
    s < 3600 && return @sprintf("%.0fm%02.0fs", s ÷ 60, s % 60)
    @sprintf("%.0fh%02.0fm", s ÷ 3600, (s % 3600) ÷ 60)
end

function write_meta(dir, p::Params, gens, stat_every)
    open(joinpath(dir, "meta.txt"), "w") do io
        println(io, "gens=", gens)
        println(io, "stat_every=", stat_every)
        println(io, "nstat=", NSTAT)
        println(io, "nres=", NRES)
        println(io, "nphen=", NPHEN)
        println(io, "ng_max=", NG_MAX)
        println(io, "started=", now())
        for f in fieldnames(Params)
            println(io, f, "=", getfield(p, f))
        end
    end
end

function write_clades(dir, ct::CladeTable)
    open(joinpath(dir, "clades.tsv"), "w") do io
        println(io, "id\tparent\tborn")
        for i in 1:nclades(ct)
            println(io, i - 1, "\t", ct.parent[i], "\t", ct.born[i])
        end
    end
end

function write_shocks(dir, shocks)
    open(joinpath(dir, "shocks.tsv"), "w") do io
        println(io, "gen\tkind")
        for (g, k) in shocks
            println(io, g, "\t", k)
        end
    end
end
