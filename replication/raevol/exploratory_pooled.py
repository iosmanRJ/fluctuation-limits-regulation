#!/usr/bin/env python3
"""EXPLORATORY (2026-09-30, PAPER.md section 11.23). From paper/data only:
(1) composition of the genome-encoded network by cell (genes, regulator proteins,
links, strength per regulator and per link) and how much of the effect of
uninformative fluctuation is explained by fewer regulators; (2) power at n = 12;
(3) pooled estimates across the independent experiments for the claims that
single experiments leave imprecise. Needs numpy.
Usage: python3 exploratory_pooled.py"""
import csv, numpy as np
from collections import defaultdict
import os
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "paper", "data") + "/"
rng = np.random.default_rng(1)
def perm(a, b, n=100000):
    a, b = np.asarray(a, float), np.asarray(b, float)
    obs = abs(a.mean() - b.mean()); x = np.concatenate([a, b]); na = len(a)
    idx = np.argsort(rng.random((n, len(x))), axis=1)
    y = x[idx]; d = np.abs(y[:, :na].mean(1) - y[:, na:].mean(1))
    return (np.sum(d >= obs - 1e-12) + 1) / (n + 1)
def cohend(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    sp = np.sqrt(((len(a)-1)*a.var(ddof=1) + (len(b)-1)*b.var(ddof=1))/(len(a)+len(b)-2))
    return (a.mean()-b.mean())/sp
def load(f):
    return list(csv.DictReader(open(D + f)))

# ---- composition
m = lambda v: float(np.mean(v))
X = load("decomp_runs.csv")
def rows(arm, q, exp="v2", pop=1024):
    return [r for r in X if r["experiment"]==exp and r["arm"]==arm and abs(float(r["rate"])-q)<1e-9 and int(r["pop"])==pop]
f = lambda rs, c: np.array([float(r[c]) for r in rs])
print("cell                 genes   tfs  links_int  S_int  S_int/tf  S_int/link_int  S_int/gene")
cells = [("froz_none",0.1),("froz_info",0.05),("froz_info",0.1),("froz_info",0.2),("froz_yoked",0.1),
         ("fluc_none",0.05),("fluc_none",0.1),("fluc_none",0.2),("fluc_yoked",0.05),("fluc_yoked",0.1),("fluc_yoked",0.2),
         ("fluc_info",0.05),("fluc_info",0.1),("fluc_info",0.2)]
for a, q in cells:
    rs = rows(a, q)
    g, t, li, si = f(rs,"genes"), f(rs,"tfs"), f(rs,"links_int"), f(rs,"S_int")
    print(f"{a:11s} q={q:<5} {g.mean():6.1f} {t.mean():5.1f} {li.mean():8.1f} {si.mean():7.2f} {si.mean()/t.mean():8.3f} {si.mean()/max(li.mean(),1e-9):11.3f} {si.mean()/g.mean():10.3f}")
# pooled regression: log S_int ~ fluctuating + log tfs (all v2 pop-1024 cells), to see how much of the fluc effect is composition
import itertools
Y, F, T, G = [], [], [], []
for a, q in cells:
    for r in rows(a, q):
        s = float(r["S_int"])
        if s <= 0: continue
        Y.append(np.log(s)); F.append(1.0 if a.startswith("fluc") and not a.endswith("info") else 0.0); T.append(np.log(float(r["tfs"]))); G.append(np.log(float(r["genes"])))
Y, F, T, G = map(np.array, (Y, F, T, G))
def ols(Xcols):
    Xm = np.column_stack([np.ones(len(Y))] + Xcols)
    b, *_ = np.linalg.lstsq(Xm, Y, rcond=None); return b
b0 = ols([F]); b1 = ols([F, T]); b2 = ols([F, T, G])
print("\nlog S_int on 'uninformative fluctuation' (no cue or yoked):")
print("  unadjusted effect: x%.2f" % np.exp(b0[1]))
print("  adjusting for log TFs: x%.2f   (TF elasticity %.2f)" % (np.exp(b1[1]), b1[2]))
print("  adjusting for log TFs and log genes: x%.2f" % np.exp(b2[1]))
# within-arm: does a stronger internal network go with a better match to the target?
print("\nwithin-cell correlation of S_int with metabolic error (negative = stronger network, better match):")
for a, q in cells:
    rs = rows(a, q)
    print(f"  {a:11s} q={q:<5} r = {np.corrcoef(f(rs,'S_int'), f(rs,'err'))[0,1]:+.2f}")

# ---- power and pooling
m = lambda v: float(np.mean(v))
rng2 = np.random.default_rng(7)
def power(d, n, sims=20000):
    a = rng2.normal(0, 1, (sims, n)); b = rng2.normal(d, 1, (sims, n))
    sp = np.sqrt((a.var(1, ddof=1) + b.var(1, ddof=1)) / 2)
    t = (b.mean(1) - a.mean(1)) / (sp * np.sqrt(2 / n))
    crit = {12: 2.074, 18: 2.032, 24: 2.013, 30: 2.002, 36: 1.994}[n]
    return float(np.mean(np.abs(t) > crit))
print("power of a two-sided test at alpha 0.05 (per test, before Holm):")
for d in (0.5, 0.8, 1.0, 1.2, 1.5, 2.0):
    print("  d=%.1f:" % d, "  ".join("n=%d %.2f" % (n, power(d, n)) for n in (12, 24, 36)))
X = load("decomp_runs.csv"); C = load("confirm_runs.csv"); R = load("raevol_runs.csv")
def col(rows, exp, arm, q, c, pop=1024):
    return [float(r[c]) for r in rows if r.get("experiment", exp) == exp and r["arm"] == arm and abs(float(r["rate"]) - q) < 1e-9 and int(r["pop"]) == pop]
def strat(pairs, B=20000):
    def stat(ps):
        da = np.mean([np.mean(a) for a, b in ps]); db = np.mean([np.mean(b) for a, b in ps]); return da - db, da / db
    d0, r0 = stat(pairs)
    boots = []
    for _ in range(B):
        ps = [(rng2.choice(a, len(a)), rng2.choice(b, len(b))) for a, b in pairs]
        boots.append(stat(ps))
    boots = np.array(boots)
    hits = 0; P = 20000
    for _ in range(P):
        ps = []
        for a, b in pairs:
            x = np.concatenate([a, b]); rng2.shuffle(x); ps.append((x[:len(a)], x[len(a):]))
        if abs(stat(ps)[0]) >= abs(d0) - 1e-12: hits += 1
    return d0, np.percentile(boots[:, 0], [2.5, 97.5]), r0, np.percentile(boots[:, 1], [2.5, 97.5]), (hits + 1) / (P + 1)
print("\nsecondary claim: genome-encoded strength, informative-cue fluctuating vs frozen informative")
for q, strata in ((0.10, [("decomp v1", col(X,"v1","fluc_info",0.1,"S_int"), col(X,"v1","froz_info",0.1,"S_int")),
                          ("decomp v2", col(X,"v2","fluc_info",0.1,"S_int"), col(X,"v2","froz_info",0.1,"S_int")),
                          ("confirm", col(C,None,"fluc_info",0.1,"S_int"), col(C,None,"froz_info",0.1,"S_int"))]),
                  (0.20, [("decomp v2", col(X,"v2","fluc_info",0.2,"S_int"), col(X,"v2","froz_info",0.2,"S_int")),
                          ("confirm", col(C,None,"fluc_info",0.2,"S_int"), col(C,None,"froz_info",0.2,"S_int"))])):
    for name, a, b in strata:
        print("  q=%.2f %-10s %.2f vs %.2f ratio %.2f d %.2f" % (q, name, m(a), m(b), m(a)/m(b), cohend(a, b)))
    d0, dci, r0, rci, p = strat([(np.array(a), np.array(b)) for _, a, b in strata])
    print("  q=%.2f POOLED (n=%d per arm): ratio %.2f [%.2f, %.2f], p=%.4f" % (q, 12 * len(strata), r0, rci[0], rci[1], p))
print("\n'information does not raise genome size': informative minus yoked, genes (fluctuating)")
strata = [(np.array(col(R,"v1","fluc_info",0.1,"genes")), np.array(col(R,"v1","fluc_yoked",0.1,"genes")))]
for q in (0.05, 0.1, 0.2):
    strata.append((np.array(col(R,"v2","fluc_info",q,"genes")), np.array(col(R,"v2","fluc_yoked",q,"genes"))))
for q in (0.1, 0.2):
    strata.append((np.array(col(C,None,"fluc_info",q,"genes")), np.array(col(C,None,"fluc_yoked",q,"genes"))))
d0, dci, r0, rci, p = strat(strata)
print("  pooled over 6 experiments x rates (n=72 per arm): %+.2f genes [%+.2f, %+.2f], p=%.3f" % (d0, dci[0], dci[1], p))
G = load("grain_runs.csv")
gc = lambda q, c: np.array([float(r[c]) for r in G if r["arm"] == "fluc_none" and abs(float(r["q"]) - q) < 1e-12])
for c in ("genes", "S"):
    d0, dci, r0, rci, p = strat([(gc(0.0002, c), gc(0.1, c))])
    print("  grain coarse minus fine, %s: %+.2f [%+.2f, %+.2f]" % (c, d0, dci[0], dci[1]))
