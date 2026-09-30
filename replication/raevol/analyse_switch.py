#!/usr/bin/env python3
"""
Pre-registered analysis of the ENVIRONMENT-SWITCH experiment (PAPER.md section
11.24, committed as 313d744 before any run). Written after the launch and before
any output was inspected; it implements the committed rules and nothing else is
a verdict.

Per run, from stats/stats_simd_best.csv (no signals, so S is the genome-encoded
network): S = nb_link * avg_influence (0 without links), nb_TF, nb_link,
nb_functional_genes, metabolic_error. START = mean over generations 1-500,
END = mean over generations 8,001-10,000. Runs that did not reach generation
10,000 are excluded and listed.

W1 loss    S_END(FU) vs S_END(FF)   two-sided permutation, 100,000
W2 growth  S_END(UF) vs S_END(UU)   two-sided permutation, 100,000
Holm across W1 and W2. For W1 also S_END(FU)/S_END(FF) with a 10,000-resample
bootstrap 95% CI.

Interpretation (fixed in 11.24):
  W1 FU < FF, Holm p < 0.05                      -> existing networks are LOST
  W1 n.s. and lower CI bound of FU/FF > 0.75     -> existing networks PERSIST
  W1 otherwise                                   -> INCONCLUSIVE about loss
  W2 UF > UU, Holm p < 0.05                      -> networks GROW once fluctuation stops
  W2 n.s.                                        -> history constrains growth within 10,000 generations
Secondary (no verdicts): trajectories in 500-generation blocks per cell; START
against END within cells; fraction lost 1 - S_END(FU)/S_START(FU).

Usage: analyse_switch.py <root>     (root contains primary/ and manifest.csv)
"""
import csv, os, random, sys
from collections import defaultdict

PERM_SEED, BOOT_SEED = 20261001, 20261002
NPERM, NBOOT = 100_000, 10_000
START, END, NGEN = (0, 500), (8000, 10000), 10000


def mean(x):
    return sum(x) / len(x)


def summarise(d):
    p = os.path.join(d, "stats", "stats_simd_best.csv")
    if not os.path.exists(p):
        return None, None, "no stats"
    rows = list(csv.reader(open(p)))
    H = {k: i for i, k in enumerate(rows[0])}
    win = {"START": defaultdict(list), "END": defaultdict(list)}
    blocks = defaultdict(lambda: defaultdict(list))
    last = -1
    for r in rows[1:]:
        try:
            g = int(float(r[0]))
        except (ValueError, IndexError):
            continue
        last = max(last, g)
        if g < 1:
            continue
        nl = float(r[H["nb_link"]])
        v = {"S": nl * float(r[H["avg_influence"]]) if nl > 0 else 0.0, "tfs": float(r[H["nb_TF"]]),
             "links": nl, "genes": float(r[H["nb_functional_genes"]]), "err": float(r[H["metabolic_error"]])}
        for name, (lo, hi) in (("START", START), ("END", END)):
            if lo < g <= hi:
                for k, x in v.items():
                    win[name][k].append(x)
        for k, x in v.items():
            blocks[(g - 1) // 500 + 1][k].append(x)
    if last < NGEN:
        return None, None, "reached only generation %d" % last
    out = {w: {k: mean(x) for k, x in win[w].items()} for w in win}
    traj = {b: {k: mean(x) for k, x in blocks[b].items()} for b in blocks}
    return out, traj, None


def perm_p(a, b, rng, n=NPERM):
    obs = abs(mean(a) - mean(b)); pool = list(a) + list(b); na = len(a); hits = 0
    for _ in range(n):
        rng.shuffle(pool)
        if abs(mean(pool[:na]) - mean(pool[na:])) >= obs - 1e-12:
            hits += 1
    return (hits + 1) / (n + 1)


def boot_ratio(a, b, rng, n=NBOOT):
    rs = []
    for _ in range(n):
        ma = mean([a[rng.randrange(len(a))] for _ in a]); mb = mean([b[rng.randrange(len(b))] for _ in b])
        rs.append(ma / mb if mb > 0 else float("inf"))
    rs.sort()
    return mean(a) / mean(b), rs[int(0.025 * n)], rs[int(0.975 * n) - 1]


def main(root):
    man = list(csv.DictReader(open(os.path.join(root, "manifest.csv"))))
    cells = defaultdict(list); trajs = defaultdict(list)
    for m in man:
        s, t, why = summarise(os.path.join(root, "primary", m["run"]))
        if s is None:
            print("EXCLUDED %s: %s" % (m["run"], why)); continue
        cells[m["cell"]].append(s); trajs[m["cell"]].append(t)
    print("CELLS (mean over runs; START = generations 1-500, END = 8,001-10,000)")
    for c in ("FF", "FU", "UU", "UF"):
        R = cells.get(c, [])
        if not R:
            print("  %s  no runs" % c); continue
        f = lambda w, k: mean([r[w][k] for r in R])
        print("  %s n=%2d  S %6.2f -> %6.2f   TFs %5.1f -> %5.1f   links %6.1f -> %6.1f   genes %5.1f -> %5.1f   err %.4f -> %.4f"
              % (c, len(R), f("START", "S"), f("END", "S"), f("START", "tfs"), f("END", "tfs"), f("START", "links"),
                 f("END", "links"), f("START", "genes"), f("END", "genes"), f("START", "err"), f("END", "err")))
    rng = random.Random(PERM_SEED); brng = random.Random(BOOT_SEED)
    end = lambda c: [r["END"]["S"] for r in cells.get(c, [])]
    ps = []
    for a, b in (("FU", "FF"), ("UF", "UU")):
        ps.append(perm_p(end(a), end(b), rng) if end(a) and end(b) else float("nan"))
    holm = sorted(range(2), key=lambda i: ps[i]); adj = [0.0, 0.0]; run = 0.0
    for k, i in enumerate(holm):
        run = max(run, min(1.0, (2 - k) * ps[i])); adj[i] = run
    print("\nPRIMARY (S_END, two-sided permutation, Holm across 2)")
    r1, lo1, hi1 = boot_ratio(end("FU"), end("FF"), brng)
    print("  W1 FU vs FF: %.3f vs %.3f  p=%.5f  Holm p=%.5f  ratio FU/FF %.3f [%.3f, %.3f]"
          % (mean(end("FU")), mean(end("FF")), ps[0], adj[0], r1, lo1, hi1))
    print("  W2 UF vs UU: %.3f vs %.3f  p=%.5f  Holm p=%.5f" % (mean(end("UF")), mean(end("UU")), ps[1], adj[1]))
    print("\nVERDICTS")
    if mean(end("FU")) < mean(end("FF")) and adj[0] < 0.05:
        v1 = "LOST: existing networks are lost under uninformative fluctuation"
    elif adj[0] >= 0.05 and lo1 > 0.75:
        v1 = "PERSIST: existing networks persist; the gap in the main experiments is a failure to grow"
    else:
        v1 = "INCONCLUSIVE about loss"
    if mean(end("UF")) > mean(end("UU")) and adj[1] < 0.05:
        v2 = "GROW: networks grow once fluctuation stops"
    else:
        v2 = "NO GROWTH within 10,000 generations: history constrains growth"
    print("  W1 " + v1); print("  W2 " + v2)
    print("\nSECONDARY (no verdicts)")
    fu = cells.get("FU", [])
    if fu:
        print("  fraction of FU network lost, 1 - S_END/S_START: %.3f" % (1 - mean([r["END"]["S"] for r in fu]) / mean([r["START"]["S"] for r in fu])))
    for c in ("FF", "FU", "UU", "UF"):
        R = cells.get(c, [])
        if R:
            d = [r["END"]["S"] - r["START"]["S"] for r in R]
            print("  %s S_END - S_START: %+.3f  (permutation p within cell, sign-flip: %.4f)" % (c, mean(d), signflip(d, rng)))
    print("\nTRAJECTORIES (S per 1,000 generations, cell means)")
    for c in ("FF", "FU", "UU", "UF"):
        T = trajs.get(c, [])
        if T:
            line = "  %s " % c
            for k in range(1, 11):
                bs = [2 * k - 1, 2 * k]
                line += "%6.2f" % mean([mean([t[b]["S"] for b in bs]) for t in T])
            print(line)


def signflip(d, rng, n=NPERM):
    obs = abs(mean(d)); hits = 0
    for _ in range(n):
        if abs(mean([x if rng.random() < 0.5 else -x for x in d])) >= obs - 1e-12:
            hits += 1
    return (hits + 1) / (n + 1)


if __name__ == "__main__":
    main(sys.argv[1])
