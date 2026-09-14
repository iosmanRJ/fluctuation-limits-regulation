#!/usr/bin/env julia
# =============================================================================
# Post-hoc dashboard. Reads a run directory and writes PNG charts.
#
#   julia --project=. scripts/dashboard.jl runs/my-run [--dark]
#
# Produces:
#   timeseries.png  — the run's history
#   fields.png      — resource fields and organism distribution at the last snapshot
# =============================================================================

using CairoMakie, EvoSim, Printf, Statistics
include(joinpath(@__DIR__, "theme.jl"))

const RUNDIR = length(ARGS) >= 1 ? ARGS[1] : error("usage: dashboard.jl <rundir> [--dark]")
const DARK   = "--dark" in ARGS

P = palette(DARK)
set_theme!(evo_theme(P))

M = read_stats(RUNDIR)
ix(n) = stat_index(n)
gens = M[ix("gen"), :]
ngen = length(gens)
@printf("loaded %d records from %s\n", ngen, RUNDIR)

# shock times, drawn as recessive vertical rules so events line up across panels
shocks = Int[]
shockfile = joinpath(RUNDIR, "shocks.tsv")
if isfile(shockfile)
    for (i, ln) in enumerate(eachline(shockfile))
        i == 1 && continue
        push!(shocks, parse(Int, split(ln, '\t')[1]))
    end
end

"""Mark mass-extinction events on a panel. Recessive: chrome, not data."""
function mark_shocks!(ax)
    isempty(shocks) && return
    vlines!(ax, shocks; color = (P.muted, 0.45), linewidth = 1, linestyle = :dash)
end

xs(v) = decimate(v, 3000)

# =============================================================================
fig = Figure(size = (1500, 1450))

Label(fig[0, 1:3], "Evolution in a petri dish — $(basename(RUNDIR))";
      fontsize = 20, color = P.primary, font = :bold, halign = :left)
Label(fig[1, 1:3],
      "$(ngen) ticks logged · dashed rules mark mass-extinction events";
      fontsize = 11, color = P.muted, halign = :left)

# --- row 1: the three headline single-series panels ---------------------------
# One series each, so no legend: the title names the measure.
ax1 = Axis(fig[2, 1], title = "Population", xlabel = "tick", ylabel = "organisms")
let (x, lo, hi) = envelope(M[ix("pop"), :])
    band!(ax1, x, lo, hi; color = (P.series[1], 0.25))
    lines!(ax1, x, (lo .+ hi) ./ 2; color = P.series[1], linewidth = 2)
end
mark_shocks!(ax1)

ax2 = Axis(fig[2, 2], title = "Generations of descent",
           xlabel = "tick", ylabel = "mean lineage depth")
let (x, y) = xs(M[ix("mean_depth"), :])
    lines!(ax2, x, y; color = P.series[3], linewidth = 2)
end
mark_shocks!(ax2)

ax3 = Axis(fig[2, 3], title = "Mean energy per organism",
           xlabel = "tick", ylabel = "energy")
let (x, y) = xs(M[ix("mean_energy"), :])
    lines!(ax3, x, y; color = P.series[2], linewidth = 2)
end
mark_shocks!(ax3)

# --- row 2: the cross-feeding story ------------------------------------------
# Channel r uses series slot r in BOTH panels, so identity carries across.
chan_labels = ["ch$(r)" * (r == 1 ? " (abiotic)" : " (waste)") for r in 1:NRES]

ax4 = Axis(fig[3, 1:2], title = "Resource pools by channel — only channel 1 has external inflow",
           xlabel = "tick", ylabel = "total in dish")
for r in 1:NRES
    (x, y) = xs(M[ix("res$(r)_total"), :])
    lines!(ax4, x, y; color = P.series[r], linewidth = 2, label = chan_labels[r])
end
mark_shocks!(ax4)
axislegend(ax4; position = :lt, nbanks = 2)

ax5 = Axis(fig[3, 3], title = "Mean uptake affinity",
           xlabel = "tick", ylabel = "affinity")
for r in 1:NRES
    (x, y) = xs(M[ix("mean_phen$(r)"), :])
    lines!(ax5, x, y; color = P.series[r], linewidth = 2, label = chan_labels[r])
end
mark_shocks!(ax5)
axislegend(ax5; position = :lt, nbanks = 2)

# --- row 3: diversity and genetic load ---------------------------------------
ax6 = Axis(fig[4, 1], title = "Clade diversity (Shannon)",
           xlabel = "tick", ylabel = "H")
let (x, y) = xs(M[ix("clade_shannon"), :])
    lines!(ax6, x, y; color = P.series[7], linewidth = 2)
end
mark_shocks!(ax6)

# separate panel rather than a second y-axis on the Shannon plot
ax7 = Axis(fig[4, 2], title = "Coexisting clades", xlabel = "tick", ylabel = "count")
let (x, y) = xs(M[ix("clades_alive"), :])
    lines!(ax7, x, y; color = P.series[6], linewidth = 2)
end
mark_shocks!(ax7)

ax8 = Axis(fig[4, 3], title = "Genetic load — Muller's ratchet",
           xlabel = "tick", ylabel = "mean deleterious mutations")
let (x, y) = xs(M[ix("mean_load"), :])
    lines!(ax8, x, y; color = P.series[8], linewidth = 2)
end
mark_shocks!(ax8)

# --- row 4: life-history traits ----------------------------------------------
life = [("body size", P_SIZE), ("motility", P_MOTIL), ("division threshold", P_DIVIDE),
        ("repair fidelity", P_REPAIR), ("stress tolerance", P_STRESS)]
ax9 = Axis(fig[5, 1:2], title = "Life-history traits (population mean)",
           xlabel = "tick", ylabel = "trait value (0–1)")
for (j, (nm, k)) in enumerate(life)
    (x, y) = xs(M[ix("mean_phen$(k)"), :])
    lines!(ax9, x, y; color = P.series[j], linewidth = 2, label = nm)
end
mark_shocks!(ax9)
axislegend(ax9; position = :lt, nbanks = 3)

# Within-population spread: rising spread is the signature of the population
# splitting into distinct ecotypes rather than tracking one moving optimum.
ax10 = Axis(fig[5, 3], title = "Trait spread (std dev of uptake affinity)",
            xlabel = "tick", ylabel = "std dev")
for r in 1:NRES
    (x, y) = xs(M[ix("std_phen$(r)"), :])
    lines!(ax10, x, y; color = P.series[r], linewidth = 2, label = chan_labels[r])
end
mark_shocks!(ax10)
axislegend(ax10; position = :lt, nbanks = 2)

for (r, h) in [(0, 26), (1, 18), (2, 220), (3, 220), (4, 220), (5, 220)]
    rowsize!(fig.layout, r, Fixed(h))
end
rowgap!(fig.layout, 14)
colgap!(fig.layout, 18)

out1 = joinpath(RUNDIR, DARK ? "timeseries-dark.png" : "timeseries.png")
save(out1, fig)
println("wrote ", out1)

# =============================================================================
# Spatial fields from the final snapshot.
# Sequential magnitude -> ONE hue, light to dark. Never a rainbow.
# =============================================================================
snaps = sort(filter(startswith("snap_"), readdir(RUNDIR)))
if !isempty(snaps)
    snap = read_snapshot(joinpath(RUNDIR, last(snaps)))
    seqmap = cgrad(P.seq)

    fig2 = Figure(size = (1700, 1000))
    Label(fig2[0, 1:4], "Dish at tick $(snap.gen)"; fontsize = 20,
          color = P.primary, font = :bold, halign = :left)
    Label(fig2[1, 1:4],
          "Each channel on its own scale — the point is spatial pattern, not cross-channel magnitude";
          fontsize = 11, color = P.muted, halign = :left)

    """
    One panel: heatmap plus its own colorbar, in a nested grid.

    A colorbar placed in the parent cell's Right() protrusion overlaps whatever
    sits in the next column, and gets its tick labels clipped at the figure
    edge. Nesting gives the colorbar a real column of its own.
    """
    function field_panel!(fig, row, col, title, data, cmap; scatter_pts = nothing)
        gl = GridLayout(fig[row, col])
        ax = Axis(gl[1, 1]; title = title, aspect = DataAspect())
        hidedecorations!(ax)
        pl = if scatter_pts === nothing
            heatmap!(ax, data; colormap = cmap)
        else
            sc = scatter!(ax, scatter_pts[1], scatter_pts[2];
                          color = data, colormap = cmap, markersize = 4)
            xlims!(ax, 1, snap.G); ylims!(ax, 1, snap.G)
            sc
        end
        Colorbar(gl[1, 2], pl; width = 9, ticklabelsize = 9, ticklabelpad = 3)
        colgap!(gl, 6)
        ax
    end

    seqmap = cgrad(P.seq)
    # Organism attributes get the SECOND sequential hue so they never read as
    # another resource channel.
    orgmap = cgrad(P.seq2)

    for r in 1:NRES
        row = 2 + (r - 1) ÷ 3
        col = 1 + (r - 1) % 3
        field_panel!(fig2, row, col, chan_labels[r], snap.res[r, :, :], seqmap)
    end

    field_panel!(fig2, 2, 4, "Organisms, shaded by genetic load", snap.load, orgmap;
                 scatter_pts = (snap.x, snap.y))

    # the true occupancy grid, recorded in the snapshot
    field_panel!(fig2, 3, 4, "Organism density (all organisms)",
                 Float32.(snap.occ), orgmap)

    rowsize!(fig2.layout, 0, Fixed(26)); rowsize!(fig2.layout, 1, Fixed(18))
    rowgap!(fig2.layout, 14); colgap!(fig2.layout, 26)

    out2 = joinpath(RUNDIR, DARK ? "fields-dark.png" : "fields.png")
    save(out2, fig2)
    println("wrote ", out2)
end

# =============================================================================
# Plain-text summary, so the run is readable without opening an image.
# =============================================================================
@printf("\n%-26s %12s %12s\n", "metric", "start", "end")
for nm in ["pop", "mean_depth", "clades_alive", "clade_shannon", "mean_load", "mean_energy"]
    i = ix(nm)
    @printf("%-26s %12.3f %12.3f\n", nm, M[i, 1], M[i, end])
end
