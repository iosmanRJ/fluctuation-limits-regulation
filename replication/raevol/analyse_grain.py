#!/usr/bin/env python3
"""
Pre-registered analysis for the R-Aevol grain experiment (PAPER.md 11.8).
Committed before any grain run; do not change after the data arrive.

Measures are imported from analyse_v2.py unchanged: one number per run, mean
over generations 25,001-30,000 of the best individual; S = total regulatory
strength (nb_link * avg_influence, 0 when there are no links); genes =
nb_functional_genes. Runs not reaching generation 30,000 are excluded, listed.

HYPOTHESIS (H_grain): the sign of the effect of cue-free fluctuation depends on
grain -- negative when fine-grained (P3/P4), non-negative when coarse-grained,
as Tsuda & Kawata 2010 report.

PRIMARY PREDICTIONS (coarse = p 0.0002, fine = p 0.1):
  G1  genes(fluc_none, coarse) > genes(fluc_none, fine)
  G2  genes(fluc_none, coarse) > genes(froz_none)          Tsuda & Kawata direction
  G3  S(fluc_none, coarse)     > S(fluc_none, fine)
  G4  S(fluc_none, coarse)     > S(froz_none)              Tsuda & Kawata direction

TEST: two-sided permutation on the difference of means, 100,000 permutations;
Holm across the four.
VERDICT per prediction: SUPPORTED if predicted sign and Holm p < 0.05;
REFUTED if opposite sign and unadjusted p < 0.05; INCONCLUSIVE otherwise.

INTERPRETATION, fixed in advance (genes: G1/G2; regulation: G3/G4):
  both supported       grain reconciles the conflict: cue-free fluctuation
                       costs complexity when fine-grained and pays when coarse
  first only           coarse grain removes the fine-grain penalty but gives no
                       gain with two alternating targets; Tsuda & Kawata's gain
                       probably needs novel optima, which this design lacks
  first not supported  grain is not the reconciler here; the conflict stays
                       unresolved and is reported as a model difference

SECONDARY (no verdicts): every cell's means; Spearman correlation of genes and
S with log10(mean dwell) across fluc_none runs, with permutation p; fine-grain
sanity check (fluc_none p 0.1 vs froz_none should reproduce P3/P4); metabolic
error by cell.

Usage:  analyse_grain.py <root>
"""
import math, os, random, sys
from collections import defaultdict
from analyse_v2 import summarise, perm_p, rank_p, holm, mean

PERM_SEED = 20260928
FINE, COARSE = 10000, 20            # rate tags: p x 1e5


def parse(name):
    """fluc_none_g20_s3 -> (arm, rate_tag)"""
    base = name.rsplit("_s", 1)[0]
    arm, g = base.rsplit("_g", 1)
    return arm, int(g)


def spearman(x, y):
    def ranks(v):
        o = sorted(range(len(v)), key=lambda i: v[i]); r = [0.0] * len(v); i = 0
        while i < len(v):
            j = i
            while j + 1 < len(v) and v[o[j + 1]] == v[o[i]]:
                j += 1
            for k in range(i, j + 1):
                r[o[k]] = (i + j) / 2.0
            i = j + 1
        return r
    rx, ry = ranks(x), ranks(y)
    mx, my = mean(rx), mean(ry)
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    den = math.sqrt(sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry))
    return num / den if den else float("nan")


def main(root):
    rng = random.Random(PERM_SEED)
    data = defaultdict(list); excluded = []
    for name in sorted(os.listdir(root)):
        full = os.path.join(root, name)
        if not os.path.isdir(full):
            continue
        s, why = summarise(full)
        if s is None:
            excluded.append("%s (%s)" % (name, why)); continue
        data[parse(name)].append(s)

    print("CELLS  (dwell = 1/p ages; 20 ages per generation)")
    for k in sorted(data, key=lambda k: (k[0], -k[1])):
        p = k[1] / 1e5
        xs = data[k]
        print("  %-10s p=%-7g dwell %6.0f ages  n=%2d  genes %6.2f  S %7.3f  links %6.1f  err %.4f"
              % (k[0], p, 1 / p, len(xs), mean([x["genes"] for x in xs]),
                 mean([x["S"] for x in xs]), mean([x["links"] for x in xs]),
                 mean([x["err"] for x in xs])))
    if excluded:
        print("EXCLUDED: " + "; ".join(excluded))

    def v(arm, tag, m):
        return [x[m] for x in data.get((arm, tag), [])]
    froz = ("froz_none", FINE)

    spec = [
        ("G1", "genes", ("fluc_none", COARSE), ("fluc_none", FINE)),
        ("G2", "genes", ("fluc_none", COARSE), froz),
        ("G3", "S",     ("fluc_none", COARSE), ("fluc_none", FINE)),
        ("G4", "S",     ("fluc_none", COARSE), froz),
    ]
    rows = []
    for gid, m, A, B in spec:
        a, b = v(*A, m), v(*B, m)
        if len(a) < 2 or len(b) < 2:
            rows.append((gid, m, A, B, None)); continue
        rows.append((gid, m, A, B, dict(a=mean(a), b=mean(b), d=mean(a) - mean(b),
                                         p=perm_p(a, b, rng), pr=rank_p(a, b, rng))))
    ps = [r[4]["p"] for r in rows if r[4]]
    adj = iter(holm(ps))
    print("\nPRIMARY (Holm across %d)" % len(ps))
    for gid, m, A, B, t in rows:
        if not t:
            print("  %s  insufficient data" % gid); continue
        t["padj"] = next(adj)
        verdict = ("SUPPORTED" if t["d"] > 0 and t["padj"] < 0.05 else
                   "REFUTED" if t["d"] < 0 and t["p"] < 0.05 else "INCONCLUSIVE")
        print("  %s  %-5s %s p=%g vs %s p=%g:  %.4f vs %.4f  diff %+.4f  p=%.5f  p_holm=%.5f  p_rank=%.5f  -> %s"
              % (gid, m, A[0], A[1] / 1e5, B[0], B[1] / 1e5, t["a"], t["b"], t["d"],
                 t["p"], t["padj"], t["pr"], verdict))

    print("\nSECONDARY (exploratory, no verdicts)")
    xs_dwell, xs_g, xs_s = [], [], []
    for (arm, tag), xs in data.items():
        if arm != "fluc_none":
            continue
        for x in xs:
            xs_dwell.append(math.log10(1e5 / tag)); xs_g.append(x["genes"]); xs_s.append(x["S"])
    for lab, ys in (("genes", xs_g), ("S", xs_s)):
        if len(ys) > 3:
            r = spearman(xs_dwell, ys); hits = 0; N = 20000
            for _ in range(N):
                sh = ys[:]; rng.shuffle(sh)
                if abs(spearman(xs_dwell, sh)) >= abs(r) - 1e-12:
                    hits += 1
            print("  dose-response  Spearman(%s, log10 dwell) = %+.3f  perm p=%.5f  (n=%d)"
                  % (lab, r, (hits + 1) / (N + 1), len(ys)))
    for m in ("genes", "S"):
        a, b = v("fluc_none", FINE, m), v(*froz, m)
        if len(a) > 1 and len(b) > 1:
            print("  fine-grain sanity  %-5s fluc_none p=0.1 %.4f vs froz_none %.4f  p=%.5f"
                  % (m, mean(a), mean(b), perm_p(a, b, rng)))


if __name__ == "__main__":
    main(sys.argv[1])
