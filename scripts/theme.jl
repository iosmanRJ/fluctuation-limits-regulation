# =============================================================================
# Shared palette and Makie theme.
#
# Colours come from a validated categorical palette and are used in FIXED slot
# order — resource channel r is always slot r, in every panel. Channel 3 is the
# same green on the totals chart and the affinity chart, so a reader who learns
# one mapping keeps it.
#
# No panel here plots two different measures on two y-scales. Where two
# quantities matter they get two panels.
# =============================================================================

const LIGHT = (
    surface   = "#fcfcfb",
    panel     = "#ffffff",
    primary   = "#0b0b0b",
    secondary = "#52514e",
    muted     = "#8a8880",
    grid      = "#e8e7e3",
    series    = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4",
                 "#008300", "#4a3aa7", "#e34948"],
    # single-hue sequential ramp (blue 100 -> 700) for magnitude fields
    seq       = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf",
                 "#184f95", "#0d366b"],
)

const DARK = (
    surface   = "#1a1a19",
    panel     = "#222221",
    primary   = "#ffffff",
    secondary = "#c3c2b7",
    muted     = "#8a8880",
    grid      = "#333331",
    series    = ["#3987e5", "#d95926", "#199e70", "#c98500", "#d55181",
                 "#008300", "#9085e9", "#e66767"],
    seq       = ["#0d366b", "#184f95", "#256abf", "#3987e5", "#6da7ec",
                 "#9ec5f4", "#cde2fb"],
)

palette(dark::Bool) = dark ? DARK : LIGHT

function evo_theme(P)
    Theme(
        backgroundcolor = P.surface,
        fontsize = 13,
        Axis = (
            backgroundcolor = P.panel,
            xgridcolor = P.grid, ygridcolor = P.grid,
            xgridwidth = 1, ygridwidth = 1,
            leftspinevisible = false, rightspinevisible = false,
            topspinevisible = false, bottomspinevisible = false,
            xtickcolor = P.grid, ytickcolor = P.grid,
            xticklabelcolor = P.secondary, yticklabelcolor = P.secondary,
            xlabelcolor = P.secondary, ylabelcolor = P.secondary,
            titlecolor = P.primary, titlesize = 14, titlealign = :left,
            titlefont = :bold, xlabelsize = 11, ylabelsize = 11,
            xticklabelsize = 10, yticklabelsize = 10,
        ),
        Legend = (
            framevisible = false, labelcolor = P.secondary, labelsize = 10,
            patchsize = (14, 2), rowgap = 1, titlecolor = P.secondary,
            titlesize = 10, titlefont = :bold,
        ),
        Colorbar = (
            ticklabelcolor = P.secondary, labelcolor = P.secondary,
            ticklabelsize = 9, labelsize = 10, spinewidth = 0, tickcolor = P.grid,
        ),
    )
end

"""
    decimate(v, target)

Reduce a long series to ~`target` points for plotting. A million-point line is
slower to draw than to simulate and resolves to nothing on screen.
"""
function decimate(v::AbstractVector, target::Int = 3000)
    n = length(v)
    n <= target && return (collect(1:n), collect(v))
    step = cld(n, target)
    idx = collect(1:step:n)
    (idx, v[idx])
end

"""Envelope decimation: keeps the min and max of each bucket so spikes survive."""
function envelope(v::AbstractVector, target::Int = 1500)
    n = length(v)
    n <= 2target && return (collect(1:n), collect(v), collect(v))
    w = cld(n, target)
    lo = Float64[]; hi = Float64[]; xs = Int[]
    for s in 1:w:n
        e = min(s + w - 1, n)
        seg = @view v[s:e]
        push!(xs, (s + e) ÷ 2); push!(lo, minimum(seg)); push!(hi, maximum(seg))
    end
    (xs, lo, hi)
end
