#!/usr/bin/env julia
# =============================================================================
# Live interactive dish.
#
#   julia --project=. -t auto scripts/live.jl            # interactive window
#   julia --project=. -t auto scripts/live.jl --frames 60 --dark
#
# Watch the dish evolve in real time and tune parameters while it runs. For
# million-tick runs use scripts/run.jl headless and scripts/dashboard.jl — this
# view is for building intuition and for seeing what a parameter change does.
#
# --frames N runs N frames and exits (used to smoke-test the viewer).
# =============================================================================

using GLMakie, EvoSim, Printf, Statistics
include(joinpath(@__DIR__, "theme.jl"))

const DARK = "--dark" in ARGS
const MAXFRAMES = let i = findfirst(==("--frames"), ARGS)
    i === nothing ? typemax(Int) : parse(Int, ARGS[i + 1])
end
const HEADLESS = MAXFRAMES != typemax(Int)

P = palette(DARK)
set_theme!(evo_theme(P))

p   = Params(grid = 192, capacity = 60_000, seed = 2024)
sim = Sim(p)

const HIST = 1200                     # rolling chart window, in frames
const MAXPTS = 12_000                 # scatter cap; drawing 60k points stalls the UI

# --- observables -------------------------------------------------------------
field    = Observable(Matrix{Float32}(sim.env.res[1, :, :]))
orgpts   = Observable(Point2f[])
orgcol   = Observable(Float32[])
hist_t   = Observable(Float32[])
hist_pop = Observable(Float32[])
hist_cl  = Observable(Float32[])
hist_ld  = Observable(Float32[])
hist_up1 = Observable(Float32[])
hist_upW = Observable(Float32[])
status   = Observable("initialising")

running  = Observable(true)
channel  = Observable(1)
spf      = Observable(4)              # sim ticks per rendered frame

# =============================================================================
# Layout: the dish is the point of this view, so it gets the left half at full
# height; charts stack down the right; controls span the bottom. Nesting each
# region in its own GridLayout keeps the DataAspect dish from being squeezed by
# the chart rows.
# =============================================================================
fig = Figure(size = (1720, 1020))

head = GridLayout(fig[1, 1:2])
Label(head[1, 1], "Petri dish — live"; fontsize = 20, color = P.primary,
      font = :bold, halign = :left)
Label(head[1, 2], status; fontsize = 12, color = P.muted, halign = :right)
colsize!(head, 1, Relative(0.42))

# --- the dish ----------------------------------------------------------------
dish = GridLayout(fig[2, 1])
axd = Axis(dish[1, 1], aspect = DataAspect(), title = "Resource field + organisms")
hidedecorations!(axd)
hm = heatmap!(axd, field; colormap = cgrad(P.seq))
# Organisms sit ON TOP of the resource field, so they need their own hue —
# shading them with the field's blue ramp makes them disappear into it.
sc = scatter!(axd, orgpts; color = orgcol, colormap = cgrad(P.seq2),
              markersize = 4, strokewidth = 0)
# two sequential scales, two colorbars, in their own columns so they cannot
# overlap each other's labels
Colorbar(dish[1, 2], hm; width = 10, ticklabelpad = 4, label = "resource")
Colorbar(dish[1, 3], sc; width = 10, ticklabelpad = 4, label = "genetic load")
xlims!(axd, 1, p.grid); ylims!(axd, 1, p.grid)
colsize!(dish, 1, Relative(0.82))
colgap!(dish, 26)

# --- live charts (each its own panel; never two scales on one axis) ----------
charts = GridLayout(fig[2, 2])
axp = Axis(charts[1, 1], title = "Population", ylabel = "organisms")
lines!(axp, hist_t, hist_pop; color = P.series[1], linewidth = 2)

axc = Axis(charts[2, 1], title = "Coexisting clades", ylabel = "count")
lines!(axc, hist_t, hist_cl; color = P.series[6], linewidth = 2)

axl = Axis(charts[3, 1], title = "Genetic load", ylabel = "mean")
lines!(axl, hist_t, hist_ld; color = P.series[8], linewidth = 2)

axu = Axis(charts[4, 1], title = "Uptake allocation", xlabel = "tick", ylabel = "share")
lines!(axu, hist_t, hist_up1; color = P.series[1], linewidth = 2, label = "channel 1")
lines!(axu, hist_t, hist_upW; color = P.series[3], linewidth = 2, label = "waste channels")
axislegend(axu; position = :rt, orientation = :horizontal, patchsize = (12, 2))
# only the bottom panel carries tick labels; the three above share its axis
for ax in (axp, axc, axl)
    hidexdecorations!(ax; grid = false)
end
rowgap!(charts, 6)

# --- controls ----------------------------------------------------------------
ctrl = GridLayout(fig[3, 1:2])

btn_run   = Button(ctrl[1, 1]; label = "pause", width = 90)
btn_reset = Button(ctrl[1, 2]; label = "reset", width = 90)
Label(ctrl[1, 3], "channel"; color = P.secondary, fontsize = 11, halign = :right)
menu_ch = Menu(ctrl[1, 4]; options = ["ch$(r)" for r in 1:NRES], default = "ch1", width = 90)
Label(ctrl[1, 5], "ticks/frame"; color = P.secondary, fontsize = 11, halign = :right)
sl_spf = Slider(ctrl[1, 6]; range = 1:32, startvalue = 4, width = 200)
Label(ctrl[1, 7], ""; width = 10)

sliders = SliderGrid(ctrl[2, 1:7],
    (label = "mutation rate (x1e-3)", range = 0.0:0.05:5.0,  startvalue = 1.0),
    (label = "resource inflow",       range = 0.02:0.01:1.0, startvalue = p.inflow),
    (label = "diffusion rate",        range = 0.0:0.01:0.24, startvalue = p.diff_rate),
    (label = "crowding cost",         range = 0.0:0.005:0.2, startvalue = p.maint_crowd),
    (label = "load penalty",          range = 0.0:0.005:0.1, startvalue = p.maint_load),
    (label = "shock interval (x1e3)", range = 1:1:100,       startvalue = 25),
)
for s in sliders.labels;      s.color = P.secondary; s.fontsize = 11; end
for s in sliders.valuelabels; s.color = P.primary;   s.fontsize = 11; end
rowgap!(ctrl, 10)

# --- wiring ------------------------------------------------------------------
on(btn_run.clicks) do _
    running[] = !running[]
    btn_run.label = running[] ? "pause" : "resume"
end
on(btn_reset.clicks) do _
    global sim = Sim(p)
    for h in (hist_t, hist_pop, hist_cl, hist_ld, hist_up1, hist_upW)
        empty!(h[]); notify(h)
    end
end
on(menu_ch.selection) do s
    channel[] = parse(Int, s[3:end])
end
on(sl_spf.value) do v
    spf[] = v
end
# Params is mutable precisely so the live view can retune a run in flight.
on(sliders.sliders[1].value) do v; sim.p.mut_rate      = Float32(v * 1e-3); end
on(sliders.sliders[2].value) do v; sim.p.inflow        = Float32(v); end
on(sliders.sliders[3].value) do v; sim.p.diff_rate     = Float32(v); end
on(sliders.sliders[4].value) do v; sim.p.maint_crowd   = Float32(v); end
on(sliders.sliders[5].value) do v; sim.p.maint_load    = Float32(v); end
on(sliders.sliders[6].value) do v; sim.p.shock_interval = Int(v) * 1000; end

# --- per-frame refresh -------------------------------------------------------
function refresh!()
    pop = sim.pop
    r = channel[]
    field[] = Matrix{Float32}(sim.env.res[r, :, :])

    hw = high_water(pop)
    live = Int[]
    @inbounds for i in 1:hw
        pop.alive[i] && push!(live, i)
    end

    if isempty(live)
        orgpts[] = Point2f[]; orgcol[] = Float32[]
        status[] = @sprintf("tick %d — EXTINCT", sim.gen)
        running[] = false
        return
    end

    sel = length(live) <= MAXPTS ? live : live[1:cld(length(live), MAXPTS):end]
    orgpts[] = [Point2f(pop.x[i], pop.y[i]) for i in sel]
    orgcol[] = Float32[pop.load[i] for i in sel]

    n = length(live)
    meanload = sum(pop.load[i] for i in live) / n
    up1 = sum(pop.phen[1, i] for i in live) / n
    upW = sum(sum(pop.phen[r2, i] for r2 in 2:NRES) for i in live) / n
    ncl = length(unique(pop.clade[i] for i in live))
    depth = sum(pop.depth[i] for i in live) / n

    for (h, v) in ((hist_t, Float32(sim.gen)), (hist_pop, Float32(n)),
                   (hist_cl, Float32(ncl)), (hist_ld, Float32(meanload)),
                   (hist_up1, Float32(up1)), (hist_upW, Float32(upW)))
        push!(h[], v)
        length(h[]) > HIST && popfirst!(h[])
        notify(h)
    end
    for ax in (axp, axc, axl, axu)
        autolimits!(ax)
    end

    status[] = @sprintf("tick %d · %d organisms · %d clades · %.0f generations · load %.2f",
                        sim.gen, n, ncl, depth, meanload)
    return
end

rowsize!(fig.layout, 1, Fixed(30))
rowsize!(fig.layout, 3, Fixed(200))
colsize!(fig.layout, 1, Relative(0.54))
rowgap!(fig.layout, 12); colgap!(fig.layout, 20)

refresh!()
display(fig)

# The loop lives in a function: at top level Julia scopes loop-body
# assignments as local, and a bare `while` here would not see its own counter.
function mainloop!()
    frames = 0
    while frames < MAXFRAMES && (HEADLESS || isopen(fig.scene))
        if running[]
            for _ in 1:spf[]
                step!(sim)
            end
            refresh!()
        end
        frames += 1
        sleep(HEADLESS ? 0.0 : 0.001)
    end
    frames
end

nframes = mainloop!()

if HEADLESS
    out = joinpath(@__DIR__, "..", "runs", DARK ? "live-check-dark.png" : "live-check.png")
    save(out, fig)
    println("ran $nframes frames to tick $(sim.gen); wrote ", normpath(out))
end
