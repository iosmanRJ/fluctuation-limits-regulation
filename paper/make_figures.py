#!/usr/bin/env python3
"""
Figures 2-6 for the manuscript, from the per-run exports in paper/data/.
One point = one population (the unit of analysis). Intervals are bootstrap 95%
CIs (10,000 resamples, fixed seed). Static print figures, PLOS sizing
(max width 7.5 in), Arial, vector PDF.

Colour encodes the CUE, the one categorical identity shared by every figure:
informative = slot 1 blue, yoked = slot 2 orange, none = slot 3 aqua (first
three slots of the reference palette; validated all-pairs, light surface:
worst CVD dE 9.2, normal-vision dE 24.0; aqua is below 3:1 contrast, so every
cue colour is also named in a legend or axis label). ENVIRONMENT is encoded by
marker fill: filled = fluctuating, hollow = frozen.
"""
import csv, os
import numpy as np
import matplotlib
matplotlib.use("pdf")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
OUT = os.path.join(HERE, "figures")
rng = np.random.default_rng(20260927)

CUE = {"info": "#2a78d6", "yoked": "#eb6834", "none": "#1baf7a"}
CUE_LABEL = {"info": "Informative cue", "yoked": "Yoked cue", "none": "No cue"}
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#898781"
GRID, AXIS = "#e1e0d9", "#c3c2b7"

plt.rcParams.update({
    "font.family": "Arial", "font.size": 8, "axes.titlesize": 8.5,
    "axes.labelsize": 8, "xtick.labelsize": 7.5, "ytick.labelsize": 7.5,
    "legend.fontsize": 7.5, "axes.edgecolor": AXIS, "axes.linewidth": 0.6,
    "xtick.color": MUTED, "ytick.color": MUTED, "xtick.labelcolor": INK2,
    "ytick.labelcolor": INK2, "axes.labelcolor": INK, "text.color": INK,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "axes.grid.axis": "y", "grid.color": GRID,
    "grid.linewidth": 0.5, "axes.axisbelow": True, "pdf.fonttype": 42,
    "figure.dpi": 150, "savefig.bbox": "tight", "savefig.pad_inches": 0.03,
})


def load(name):
    with open(os.path.join(DATA, name)) as f:
        return list(csv.DictReader(f))


def boot_mean(x, n=10000):
    x = np.asarray(x, float)
    b = rng.choice(x, size=(n, len(x)), replace=True).mean(1)
    return x.mean(), np.percentile(b, 2.5), np.percentile(b, 97.5)


def boot_diff(a, b, n=10000):
    a, b = np.asarray(a, float), np.asarray(b, float)
    ba = rng.choice(a, size=(n, len(a)), replace=True).mean(1)
    bb = rng.choice(b, size=(n, len(b)), replace=True).mean(1)
    d = ba - bb
    return a.mean() - b.mean(), np.percentile(d, 2.5), np.percentile(d, 97.5)


def strip(ax, x, vals, color, filled, width=0.22):
    """Points (jittered, with a surface ring) plus mean and 95% CI to the right."""
    vals = np.asarray(vals, float)
    jit = (rng.random(len(vals)) - 0.5) * width
    ax.scatter(x - 0.08 + jit, vals, s=16, zorder=3,
               facecolor=color if filled else "white", edgecolor=color if not filled else "white",
               linewidth=0.9 if not filled else 0.5, alpha=0.9 if filled else 1.0)
    m, lo, hi = boot_mean(vals)
    ax.plot([x + 0.2, x + 0.2], [lo, hi], color=INK, lw=1.0, zorder=4, solid_capstyle="round")
    ax.scatter([x + 0.2], [m], s=26, zorder=5, marker="D",
               facecolor=color if filled else "white", edgecolor=INK, linewidth=0.7)


def env_legend(fig, color=MUTED, y=-0.08):
    from matplotlib.lines import Line2D
    h = [Line2D([], [], ls="", marker="o", ms=5, mfc=color, mec="white", label="Fluctuating"),
         Line2D([], [], ls="", marker="o", ms=5, mfc="white", mec=color, mew=0.9, label="Frozen"),
         Line2D([], [], ls="", marker="D", ms=4.5, mfc="white", mec=INK, label="Mean, 95% CI")]
    fig.legend(handles=h, loc="lower center", ncol=3, frameon=False, handletextpad=0.3,
               bbox_to_anchor=(0.5, y))


def panel_label(ax, s):
    ax.text(-0.02, 1.04, s, transform=ax.transAxes, fontsize=10, fontweight="bold",
            ha="right", va="bottom", color=INK)


# ----------------------------------------------------------------------------
# Fig 2 -- EvoSim: more connections, less strength (survivors only)
# ----------------------------------------------------------------------------
def fig2():
    E = [r for r in load("evosim_runs.csv") if r["extinct"] == "false"]
    fl = [r for r in E if r["arm"] == "fluc_none"]
    fz = [r for r in E if r["arm"] == "froz_none"]
    col = CUE["none"]
    fig, axes = plt.subplots(1, 4, figsize=(7.5, 2.3), gridspec_kw={"wspace": 0.55})
    specs = [("count_density", "Count density\n(non-zero / possible)", "A"),
             ("strength_density", "Strength density\n(Σ|w| / possible)", "B"),
             ("total_strength", "Total wiring strength\n(Σ|w|)", "C")]
    for ax, (k, lab, pl) in zip(axes[:3], specs):
        a = [float(r[k]) for r in fl]; b = [float(r[k]) for r in fz]
        strip(ax, 0, a, col, True); strip(ax, 1, b, col, False)
        ax.set_xticks([0, 1]); ax.set_xticklabels(["Fluctuating", "Frozen"])
        ax.set_xlim(-0.55, 1.55); ax.set_ylabel(lab)
        ch = 100 * (np.mean(a) - np.mean(b)) / np.mean(b)
        ax.set_title("%+.1f%%" % ch, color=INK2, fontsize=8)
        panel_label(ax, pl)
    ax = axes[3]
    for rows, filled in ((fl, True), (fz, False)):
        ax.scatter([float(r["pop"]) / 1000 for r in rows], [float(r["strength_density"]) for r in rows],
                   s=16, facecolor=col if filled else "white", edgecolor="white" if filled else col,
                   linewidth=0.5 if filled else 0.9, zorder=3)
    ax.set_xlabel("Population (thousands)"); ax.set_ylabel("Strength density")
    ax.grid(axis="x", color=GRID, lw=0.5)
    panel_label(ax, "D")
    env_legend(fig, color=col, y=-0.1)
    fig.savefig(os.path.join(OUT, "fig2_evosim.pdf"))
    plt.close(fig)


# ----------------------------------------------------------------------------
# Fig 3 -- R-Aevol exploratory six arms
# ----------------------------------------------------------------------------
def fig3():
    R = [r for r in load("raevol_runs.csv") if r["experiment"] == "v1"]
    fig, axes = plt.subplots(1, 3, figsize=(7.5, 2.5), gridspec_kw={"wspace": 0.38})
    specs = [("S", "Total regulatory strength", "A"),
             ("links", "Regulatory links (count)", "B"),
             ("genes", "Functional genes", "C")]
    cues = ["info", "yoked", "none"]
    for ax, (k, lab, pl) in zip(axes, specs):
        for i, c in enumerate(cues):
            a = [float(r[k]) for r in R if r["arm"] == "fluc_" + c]
            b = [float(r[k]) for r in R if r["arm"] == "froz_" + c]
            strip(ax, 2 * i, a, CUE[c], True)
            strip(ax, 2 * i + 1, b, CUE[c], False)
        ax.set_xticks([0.5, 2.5, 4.5]); ax.set_xticklabels([CUE_LABEL[c].replace(" cue", "\ncue") for c in cues])
        ax.set_xlim(-0.6, 5.6); ax.set_ylabel(lab)
        for xx in (1.5, 3.5):
            ax.axvline(xx, color=GRID, lw=0.5, zorder=0)
        panel_label(ax, pl)
    env_legend(fig, y=-0.12)
    fig.savefig(os.path.join(OUT, "fig3_raevol_v1.pdf"))
    plt.close(fig)


# ----------------------------------------------------------------------------
# Fig 4 -- pre-registered contrasts with bootstrap 95% CIs
# ----------------------------------------------------------------------------
def fig4():
    R = [r for r in load("raevol_runs.csv") if r["experiment"] == "v2"]

    def v(arm, rate, pop, k):
        return [float(r[k]) for r in R if r["arm"] == arm and abs(float(r["rate"]) - rate) < 1e-9
                and int(r["pop"]) == pop]

    conds = [(0.05, 1024, "0.05"), (0.10, 1024, "0.10"), (0.20, 1024, "0.20"), (0.10, 512, "0.10\nN = 512")]
    preds = [
        ("P1", "S", "fluc_info", "froz_info", "info", "Informative cue:\nfluctuating − frozen", True),
        ("P2", "S", "fluc_info", "fluc_yoked", "info", "Fluctuating:\ninformative − yoked", True),
        ("P3", "S", "fluc_none", "froz_none", "none", "No cue:\nfluctuating − frozen", False),
        ("P4", "genes", "fluc_info", "froz_info", "info", "Informative cue:\nfluctuating − frozen", False),
    ]
    fig, axes = plt.subplots(1, 4, figsize=(7.5, 2.55), gridspec_kw={"wspace": 0.5})
    for ax, (pid, k, A, B, cue, title, needs_check) in zip(axes, preds):
        for i, (rate, pop, lab) in enumerate(conds):
            a = v(A, rate, pop, k)
            bpop, brate = (1024, 0.10) if B == "froz_none" else (pop, rate)
            b = v(B, brate, bpop, k)
            d, lo, hi = boot_diff(a, b)
            evaluable = not (needs_check and abs(rate - 0.05) < 1e-9)
            secondary = pop == 512
            c = CUE[cue]
            ax.plot([i, i], [lo, hi], color=INK, lw=1.0, zorder=3)
            ax.scatter([i], [d], s=30, zorder=4,
                       facecolor=c if (evaluable and not secondary) else "white",
                       edgecolor=c if not (evaluable and not secondary) else INK,
                       linewidth=1.0 if not (evaluable and not secondary) else 0.6,
                       marker="o" if not secondary else "s")
        ax.axhline(0, color=AXIS, lw=0.8, zorder=1)
        ax.axvline(2.5, color=GRID, lw=0.5, zorder=0)
        ax.set_xticks(range(4)); ax.set_xticklabels([c[2] for c in conds])
        ax.set_xlim(-0.5, 3.5)
        ax.set_xlabel("Switching probability")
        ax.set_ylabel("Δ genes" if k == "genes" else "Δ total regulatory strength")
        ax.set_title(title, fontsize=7.5, color=INK2)
        panel_label(ax, pid)
    from matplotlib.lines import Line2D
    h = [Line2D([], [], ls="", marker="o", ms=5, mfc=MUTED, mec=INK, mew=0.6, label="Pre-registered, evaluable"),
         Line2D([], [], ls="", marker="o", ms=5, mfc="white", mec=MUTED, mew=1.0, label="Excluded by manipulation check"),
         Line2D([], [], ls="", marker="s", ms=5, mfc="white", mec=MUTED, mew=1.0, label="Secondary (N = 512)")]
    fig.legend(handles=h, loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.5, -0.2))
    fig.savefig(os.path.join(OUT, "fig4_preregistered.pdf"))
    plt.close(fig)


# ----------------------------------------------------------------------------
# Fig 5 -- counting versus weighing
# ----------------------------------------------------------------------------
def fig5():
    R = [r for r in load("raevol_runs.csv") if r["experiment"] == "v2"]
    E = [r for r in load("evosim_runs.csv") if r["extinct"] == "false"]

    def mv(arm, rate, k):
        x = [float(r[k]) for r in R if r["arm"] == arm and abs(float(r["rate"]) - rate) < 1e-9
             and int(r["pop"]) == 1024]
        return np.mean(x)

    fig, ax = plt.subplots(figsize=(3.6, 3.3))
    lim = (0.12, 4.0)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(*lim); ax.set_ylim(*lim)
    ax.plot(lim, lim, ls=(0, (3, 3)), color=AXIS, lw=0.8, zorder=1)
    ax.axhline(1, color=GRID, lw=0.6, zorder=0); ax.axvline(1, color=GRID, lw=0.6, zorder=0)
    ax.grid(False)
    pts = []
    for rate in (0.05, 0.10, 0.20):
        pts.append(("info", mv("fluc_info", rate, "links") / mv("froz_info", rate, "links"),
                    mv("fluc_info", rate, "S") / mv("froz_info", rate, "S"), "%.2f" % rate))
        pts.append(("none", mv("fluc_none", rate, "links") / mv("froz_none", 0.10, "links"),
                    mv("fluc_none", rate, "S") / mv("froz_none", 0.10, "S"), "%.2f" % rate))
    for cue, x, y, lab in pts:
        ax.scatter([x], [y], s=34, facecolor=CUE[cue], edgecolor="white", linewidth=0.7, zorder=4)
        if cue == "info":
            ax.annotate(lab, (x, y), xytext=(4, 2), textcoords="offset points", fontsize=6.5, color=INK2)
    nx = max(p[1] for p in pts if p[0] == "none"); ny = max(p[2] for p in pts if p[0] == "none")
    ax.annotate("all three rates", (nx, ny), xytext=(6, 4), textcoords="offset points",
                fontsize=6.5, color=INK2)
    fl = [r for r in E if r["arm"] == "fluc_none"]; fz = [r for r in E if r["arm"] == "froz_none"]
        # raw counts and total strength, as for R-Aevol; the density measures share the
    # n_g^2 denominator, which grows with the extra nodes of fluctuating networks
    ex = np.mean([float(r["edges"]) for r in fl]) / np.mean([float(r["edges"]) for r in fz])
    ey = np.mean([float(r["total_strength"]) for r in fl]) / np.mean([float(r["total_strength"]) for r in fz])
    ax.scatter([ex], [ey], s=40, marker="^", facecolor=CUE["none"], edgecolor=INK, linewidth=0.6, zorder=5)
    ax.annotate("EvoSim", (ex, ey), xytext=(5, -9), textcoords="offset points", fontsize=6.5, color=INK2)
    ax.text(0.14, 3.2, "count down,\nstrength up", fontsize=6.5, color=MUTED, va="top")
    ax.text(3.6, 0.14, "count up,\nstrength down", fontsize=6.5, color=MUTED, ha="right", va="bottom")
    ax.text(2.2, 2.9, "measures\nagree", fontsize=6.5, color=MUTED, ha="center", va="top")
    ax.set_xlabel("Fluctuating / frozen: connection count")
    ax.set_ylabel("Fluctuating / frozen: connection strength")
    from matplotlib.lines import Line2D
    h = [Line2D([], [], ls="", marker="o", ms=5.5, mfc=CUE["info"], mec="white", label="R-Aevol, informative cue"),
         Line2D([], [], ls="", marker="o", ms=5.5, mfc=CUE["none"], mec="white", label="R-Aevol, no cue"),
         Line2D([], [], ls="", marker="^", ms=6, mfc=CUE["none"], mec=INK, mew=0.6, label="EvoSim, no cue")]
    ax.legend(handles=h, loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=1,
              frameon=False, handletextpad=0.3)
    ax.set_xticks([0.25, 0.5, 1, 2, 4]); ax.set_xticklabels(["0.25", "0.5", "1", "2", "4"])
    ax.set_yticks([0.25, 0.5, 1, 2, 4]); ax.set_yticklabels(["0.25", "0.5", "1", "2", "4"])
    ax.minorticks_off()
    fig.savefig(os.path.join(OUT, "fig5_count_vs_strength.pdf"))
    plt.close(fig)
    return pts, (ex, ey)


# ----------------------------------------------------------------------------
# Fig 4 -- signal-protein links versus the genome-encoded network (EXPLORATORY)
# Continuations of every run from its saved generation-30,000 population to
# 32,000 with instrumented statistics (replication/raevol/decomp_*.py).
# ----------------------------------------------------------------------------
def fig_decomp():
    C = load("decomp_runs.csv")

    def vals(arm, rate, k, pop=1024):
        return [float(r[k]) for r in C if r["experiment"] == "v2" and r["arm"] == arm
                and abs(float(r["rate"]) - rate) < 1e-9 and int(r["pop"]) == pop]

    groups = [("q = 0.05", [("froz_info", 0.05), ("fluc_info", 0.05), ("fluc_yoked", 0.05), ("fluc_none", 0.05)]),
              ("q = 0.10", [("froz_info", 0.10), ("fluc_info", 0.10), ("fluc_yoked", 0.10), ("fluc_none", 0.10)]),
              ("q = 0.20", [("froz_info", 0.20), ("fluc_info", 0.20), ("fluc_yoked", 0.20), ("fluc_none", 0.20)]),
              ("frozen,\nother cues", [("froz_yoked", 0.10), ("froz_none", 0.10)])]
    fig, axes = plt.subplots(1, 2, figsize=(7.3, 2.9))
    for ax, (ktot, ksig, kint, ylab, lab) in zip(axes, [
            ("S", "S_sig", "S_int", "Total regulatory strength", "A"),
            ("links", "links_sig", "links_int", "Regulatory links (count)", "B")]):
        x, ticks, tlabels = 0.0, [], []
        for gname, cells in groups:
            start = x
            for arm, rate in cells:
                cue = arm.split("_")[1]; fluc = arm.startswith("fluc")
                col = CUE[cue]
                ti, si = np.mean(vals(arm, rate, kint)), np.mean(vals(arm, rate, ksig))
                face = col if fluc else "white"
                ax.bar(x, ti, width=0.72, color=face, edgecolor=col, linewidth=0.9, zorder=3)
                ax.bar(x, si, width=0.72, bottom=ti, color="white" if not fluc else col, alpha=1.0 if not fluc else 0.35,
                       edgecolor=col, linewidth=0.9, hatch="////", zorder=3)
                m, lo, hi = boot_mean(vals(arm, rate, kint))
                ax.plot([x - 0.12, x - 0.12], [lo, hi], color=INK, lw=0.9, zorder=5)
                m, lo, hi = boot_mean(vals(arm, rate, ktot))
                ax.plot([x + 0.12, x + 0.12], [lo, hi], color=MUTED, lw=0.9, zorder=5)
                x += 0.85
            ticks.append((start + x - 0.85) / 2); tlabels.append(gname)
            x += 0.9
        ax.set_xticks(ticks); ax.set_xticklabels(tlabels)
        ax.tick_params(axis="x", length=0)
        ax.set_ylabel(ylab)
        panel_label(ax, lab)
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    h = [Patch(facecolor=CUE["info"], edgecolor=CUE["info"], label="Informative cue"),
         Patch(facecolor=CUE["yoked"], edgecolor=CUE["yoked"], label="Yoked cue"),
         Patch(facecolor=CUE["none"], edgecolor=CUE["none"], label="No cue"),
         Patch(facecolor="white", edgecolor=MUTED, label="Frozen (hollow)"),
         Patch(facecolor=MUTED, edgecolor=MUTED, label="Genome-encoded"),
         Patch(facecolor="white", edgecolor=MUTED, hatch="////", label="Signal proteins"),
         Line2D([], [], color=INK, lw=0.9, label="95% CI, genome-encoded"),
         Line2D([], [], color=MUTED, lw=0.9, label="95% CI, total")]
    fig.legend(handles=h, loc="lower center", ncol=4, frameon=False, bbox_to_anchor=(0.5, -0.2),
               handletextpad=0.4, columnspacing=1.2)
    fig.savefig(os.path.join(OUT, "fig4_decomposition.pdf"))
    plt.close(fig)


# ----------------------------------------------------------------------------
# Fig 6 -- grain: cue-free fluctuation against mean environmental dwell time
# ----------------------------------------------------------------------------
def fig6():
    G = load("grain_runs.csv")
    col = CUE["none"]
    fl = [r for r in G if r["arm"] == "fluc_none"]
    fz = [r for r in G if r["arm"] == "froz_none"]
    dwells = sorted({float(r["dwell_ages"]) for r in fl})
    fig, axes = plt.subplots(1, 3, figsize=(7.5, 2.45), gridspec_kw={"wspace": 0.42})
    specs = [("genes", "Functional genes", "A"),
             ("S", "Total regulatory strength", "B"),
             ("err", "Metabolic error", "C")]
    for ax, (k, lab, pl) in zip(axes, specs):
        z = [float(r[k]) for r in fz]
        m, lo, hi = boot_mean(z)
        ax.axhspan(lo, hi, color=GRID, zorder=0, lw=0)
        ax.axhline(m, color=MUTED, lw=0.9, ls=(0, (4, 2)), zorder=1)
        ax.text(0.03, m, " frozen", transform=ax.get_yaxis_transform(), fontsize=6.5,
                color=INK2, va="bottom", ha="left")
        for dw in dwells:
            v = [float(r[k]) for r in fl if abs(float(r["dwell_ages"]) - dw) < 1e-6]
            jit = 10 ** ((rng.random(len(v)) - 0.5) * 0.12)
            ax.scatter(np.array([dw] * len(v)) * jit / 1.25, v, s=14, facecolor=col,
                       edgecolor="white", linewidth=0.5, zorder=3)
            mm, l2, h2 = boot_mean(v)
            ax.plot([dw * 1.25, dw * 1.25], [l2, h2], color=INK, lw=1.0, zorder=4)
            ax.scatter([dw * 1.25], [mm], s=24, marker="D", facecolor=col, edgecolor=INK,
                       linewidth=0.6, zorder=5)
        ax.set_xscale("log")
        ax.set_xticks([10, 100, 1000, 5000]); ax.set_xticklabels(["10", "100", "1,000", "5,000"])
        ax.minorticks_off()
        ax.set_xlim(5.5, 11000)
        ax.set_xlabel("Mean dwell of an environment (ages)")
        ax.set_ylabel(lab)
        if k == "err":
            ax.set_ylim(bottom=0)
        else:
            ax.set_ylim(bottom=0 if k == "S" else None)
        panel_label(ax, pl)
    from matplotlib.lines import Line2D
    h = [Line2D([], [], ls="", marker="o", ms=5, mfc=col, mec="white", label="Fluctuating, no cue (one point per population)"),
         Line2D([], [], ls="", marker="D", ms=4.5, mfc=col, mec=INK, label="Mean, 95% CI"),
         Line2D([], [], color=MUTED, ls=(0, (4, 2)), lw=0.9, label="Frozen mean (band: 95% CI)")]
    fig.legend(handles=h, loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.5, -0.13))
    fig.savefig(os.path.join(OUT, "fig6_grain.pdf"))
    plt.close(fig)


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    fig2()        # EvoSim; now shown in S6 Text
    fig3(); fig4(); fig_decomp(); fig6()
    print("figures written to", OUT)
