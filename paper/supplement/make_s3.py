#!/usr/bin/env python3
"""Generate the tables of S3 Table from the per-run exports in paper/data/.
Each cell: mean +/- SD over populations (the unit of analysis), n shown.
Run from paper/supplement/:  python3 make_s3.py  -> s3_tables.tex"""
import csv, math, os

D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")


def load(n):
    return list(csv.DictReader(open(os.path.join(D, n))))


def ms(v, fmt="%.2f"):
    v = [x for x in v if not math.isnan(x)]
    if not v:
        return "--"
    m = sum(v) / len(v)
    sd = math.sqrt(sum((x - m) ** 2 for x in v) / (len(v) - 1)) if len(v) > 1 else 0.0
    return (fmt % m) + r" $\pm$ " + (fmt % sd)


ARM = {"fluc_info": "Fluctuating, informative", "fluc_yoked": "Fluctuating, yoked",
       "fluc_none": "Fluctuating, no cue", "froz_info": "Frozen, informative",
       "froz_none": "Frozen, no cue", "froz_yoked": "Frozen, yoked"}
ORDER = ["fluc_info", "fluc_yoked", "fluc_none", "froz_info", "froz_none", "froz_yoked"]
out = []

# ---- Table A: EvoSim -----------------------------------------------------
E = load("evosim_runs.csv")
out.append(r"""\begin{table}[!ht]\small
\caption{\textbf{EvoSim, 180 runs} (grid 128, 150{,}000 steps, 30 seeds per arm).
Populations that went extinct are counted and excluded from the means. The three
frozen arms are identical runs (the cue is zero with a single state) and are
listed once.}
\begin{tabular}{@{}lcccccc@{}}\toprule
Arm & Survived & Count density & Strength density & Total strength & Nodes & Population\\\midrule""")
for a in ["fluc_info", "fluc_yoked", "fluc_none", "froz_none"]:
    rows = [r for r in E if r["arm"] == a]
    surv = [r for r in rows if r["extinct"] == "false"]
    f = lambda k: [float(r[k]) for r in surv]
    lab = "Frozen (all three arms)" if a == "froz_none" else ARM[a]
    out.append("%s & %d/%d & %s & %s & %s & %s & %s\\\\" % (
        lab, len(surv), len(rows), ms(f("count_density"), "%.3f"), ms(f("strength_density"), "%.3f"),
        ms(f("total_strength"), "%.1f"), ms(f("ng"), "%.2f"), ms(f("pop"), "%.0f")))
out.append(r"\bottomrule\end{tabular}" + "\n" + r"""\par\footnotesize Time steps per generation, from the slope of mean lineage depth
over the second half of each run: fluctuating arms 23.6--24.3, frozen 30.7.
\end{table}""")

# ---- Table B/C: R-Aevol --------------------------------------------------
R = load("raevol_runs.csv")


def raevol_table(rows, caption, cells):
    t = [r"\begin{table}[!ht]\small", r"\caption{" + caption + "}",
         r"\begin{tabular}{@{}lccccc@{}}\toprule",
         r"Cell & $n$ & Total strength $S$ & Links & Genes & Metabolic error\\\midrule"]
    for lab, sel in cells:
        v = [r for r in rows if sel(r)]
        if not v:
            continue
        f = lambda k: [float(r[k]) for r in v]
        t.append("%s & %d & %s & %s & %s & %s\\\\" % (lab, len(v), ms(f("S"), "%.2f"), ms(f("links"), "%.1f"),
                                                     ms(f("genes"), "%.2f"), ms(f("err"), "%.4f")))
    t.append(r"\bottomrule\end{tabular}\end{table}")
    return "\n".join(t)


v1 = [r for r in R if r["experiment"] == "v1"]
out.append(raevol_table(v1, r"\textbf{R-Aevol, exploratory experiment} (72 runs; $q=0.10$, population 1024, 30{,}000 generations).",
                        [(ARM[a], (lambda a: lambda r: r["arm"] == a)(a)) for a in ORDER]))

v2 = [r for r in R if r["experiment"] == "v2"]
cells = []
for q in (0.05, 0.10, 0.20):
    for a in ["fluc_info", "fluc_yoked", "fluc_none", "froz_info"]:
        cells.append(("%s, $q=%.2f$" % (ARM[a], q),
                      (lambda a, q: lambda r: r["arm"] == a and abs(float(r["rate"]) - q) < 1e-9 and r["pop"] == "1024")(a, q)))
for a in ["froz_none", "froz_yoked"]:
    cells.append(("%s, $q=0.10$" % ARM[a], (lambda a: lambda r: r["arm"] == a and r["pop"] == "1024")(a)))
for a in ["fluc_info", "fluc_yoked", "fluc_none", "froz_info"]:
    cells.append(("%s, $q=0.10$, $N=512$" % ARM[a], (lambda a: lambda r: r["arm"] == a and r["pop"] == "512")(a)))
out.append(raevol_table(v2, r"\textbf{R-Aevol, pre-registered sweep} (216 runs, 30{,}000 generations; population 1024 unless marked).", cells))

# ---- Table D: grain ------------------------------------------------------
G = load("grain_runs.csv")
t = [r"\begin{table}[!ht]\small",
     r"\caption{\textbf{R-Aevol, grain experiment} (60 runs, no cue, population 1024, 30{,}000 generations). Mean dwell is $1/q$ ages; 20 ages make one generation.}",
     r"\begin{tabular}{@{}lccccc@{}}\toprule",
     r"Cell & $n$ & Total strength $S$ & Links & Genes & Metabolic error\\\midrule"]
for q in (0.1, 0.01, 0.001, 0.0002):
    v = [r for r in G if r["arm"] == "fluc_none" and abs(float(r["q"]) - q) < 1e-12]
    f = lambda k: [float(r[k]) for r in v]
    t.append("Fluctuating, $q=%g$ (dwell %s ages) & %d & %s & %s & %s & %s\\\\" % (
        q, "{:,}".format(int(round(1 / q))).replace(",", "{,}"), len(v), ms(f("S"), "%.2f"), ms(f("links"), "%.1f"),
        ms(f("genes"), "%.2f"), ms(f("err"), "%.4f")))
v = [r for r in G if r["arm"] == "froz_none"]
f = lambda k: [float(r[k]) for r in v]
t.append("Frozen & %d & %s & %s & %s & %s\\\\" % (len(v), ms(f("S"), "%.2f"), ms(f("links"), "%.1f"),
                                                ms(f("genes"), "%.2f"), ms(f("err"), "%.4f")))
t.append(r"\bottomrule\end{tabular}\end{table}")
out.append("\n".join(t))

# ---- Table E: signal-link decomposition (exploratory continuations) --------
C = load("decomp_runs.csv")
t = [r"\begin{table}[!ht]\small",
     r"\caption{\textbf{R-Aevol, signal-link decomposition} (exploratory; 252 runs continued from their saved generation-30{,}000 populations to 32{,}000, population 1024 unless stated). $S_\mathrm{sig}$ is the part of the regulatory strength exerted by signal proteins, $S_\mathrm{int}$ the genome-encoded part.}",
     r"\begin{tabular}{@{}llcccccc@{}}\toprule",
     r"Experiment & Cell & $n$ & $S$ & $S_\mathrm{sig}$ & $S_\mathrm{int}$ & Signal links & Genome links\\\midrule"]
cells = [("v2", a, q, 1024) for q in (0.05, 0.1, 0.2) for a in ("fluc_info", "fluc_yoked", "fluc_none", "froz_info")] \
      + [("v2", "froz_yoked", 0.1, 1024), ("v2", "froz_none", 0.1, 1024)] \
      + [("v2", a, 0.1, 512) for a in ("fluc_info", "fluc_yoked", "froz_info")] \
      + [("v1", a, 0.1, 1024) for a in ("fluc_info", "froz_info", "fluc_yoked", "froz_yoked")]
for exp, a, q, pop in cells:
    v = [r for r in C if r["experiment"] == exp and r["arm"] == a and abs(float(r["rate"]) - q) < 1e-9 and int(r["pop"]) == pop]
    f = lambda k: [float(r[k]) for r in v]
    lab = "%s, $q=%.2f$%s" % (ARM[a], q, "" if pop == 1024 else ", population 512")
    t.append("%s & %s & %d & %s & %s & %s & %s & %s\\\\" % ("Sweep" if exp == "v2" else "Exploratory", lab, len(v),
             ms(f("S")), ms(f("S_sig")), ms(f("S_int")), ms(f("links_sig"), "%.1f"), ms(f("links_int"), "%.1f")))
t.append(r"\bottomrule\end{tabular}\end{table}")
out.append("\n".join(t))

# ---- Table F: pre-registered confirmation -------------------------------
F = load("confirm_runs.csv")
t = [r"\begin{table}[!ht]\small",
     r"\caption{\textbf{R-Aevol, pre-registered confirmation} (132 runs with the instrumented platform from generation 0, population 1024, 30{,}000 generations; means over generations 25{,}001--30{,}000).}",
     r"\begin{tabular}{@{}lccccccc@{}}\toprule",
     r"Cell & $n$ & $S$ & $S_\mathrm{sig}$ & $S_\mathrm{int}$ & Links & Genes & Metabolic error\\\midrule"]
for q in (0.1, 0.2):
    for a in ("fluc_info", "fluc_yoked", "fluc_none", "froz_info", "froz_yoked"):
        v = [r for r in F if r["arm"] == a and abs(float(r["rate"]) - q) < 1e-9]
        f = lambda k: [float(r[k]) for r in v]
        t.append("%s, $q=%.2f$ & %d & %s & %s & %s & %s & %s & %s\\\\" % (ARM[a], q, len(v), ms(f("S")), ms(f("S_sig")), ms(f("S_int")),
                 ms(f("links"), "%.1f"), ms(f("genes")), ms(f("err"), "%.4f")))
v = [r for r in F if r["arm"] == "froz_none"]
f = lambda k: [float(r[k]) for r in v]
t.append("%s & %d & %s & %s & %s & %s & %s & %s\\\\" % (ARM["froz_none"], len(v), ms(f("S")), ms(f("S_sig")), ms(f("S_int")),
         ms(f("links"), "%.1f"), ms(f("genes")), ms(f("err"), "%.4f")))
t.append(r"\bottomrule\end{tabular}\end{table}")
out.append("\n".join(t))

# ---- Table G: pre-registered environment switch -------------------------
W = load("switch_trajectories.csv")
t = [r"\begin{table}[!ht]\small",
     r"\caption{\textbf{R-Aevol, pre-registered environment switch} (48 runs, 12 per cell, each founded by the best genome of a cue-free confirmation population at $q=0.10$, generation 30{,}000, and run for 10{,}000 generations; source over generations 25{,}001--30{,}000, start over 1--500, end over 8{,}001--10{,}000; no signals, so $S$ is the genome-encoded network).}",
     r"\begin{tabular}{@{}lccccccc@{}}\toprule",
     r"Cell & $n$ & $S$ source & $S$ start & $S$ end & Regulators end & Genes end & Metabolic error end\\\midrule"]
CELLNAME = {"FF": "Frozen, kept frozen", "FU": "Frozen, moved to fluctuation", "UU": "Fluctuating, kept fluctuating", "UF": "Fluctuating, moved to frozen"}
for c in ("FF", "FU", "UU", "UF"):
    reps = sorted({r["rep"] for r in W if r["cell"] == c}, key=int)
    def per_run(block_set, k):
        return [sum(float(r[k]) for r in W if r["cell"] == c and r["rep"] == s and int(r["block"]) in block_set) / len(block_set) for s in reps]
    end = {17, 18, 19, 20}
    t.append("%s & %d & %s & %s & %s & %s & %s & %s\\\\" % (CELLNAME[c], len(reps), ms(per_run({0}, "S")), ms(per_run({1}, "S")),
             ms(per_run(end, "S")), ms(per_run(end, "tfs"), "%.1f"), ms(per_run(end, "genes"), "%.1f"), ms(per_run(end, "err"), "%.4f")))
t.append(r"\bottomrule\end{tabular}\end{table}")
out.append("\n".join(t))

open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "s3_tables.tex"), "w").write("\n\n".join(out) + "\n")
print("wrote s3_tables.tex")
