"""What is the ~0.033 t-independent excess in the empirical temporal F_c? (CLAUDE.md 7.9)

THE FINDING THIS EXPLAINS. fc_chrom_check.py TEST 4 measured that the empirical F_c sits about
0.033 above the simulated one by a margin that does NOT grow with the generation gap, while the
gap CONTRAST -- which cancels any t-independent term -- implies POPMULT ~ 6,000, inside the prior
and compatible with F_st. So the level carries a contaminant roughly three times the size of the
whole POPMULT signal. This script asks what it is made of.

THE HYPOTHESIS, and it is not a fishing expedition -- 7.2 named the mechanism before F_c existed.
Two sites are already known to look like samples of close relatives rather than draws from a deme:
Mortensen9-2015 and H41-2023, each simultaneously the LOWEST-pi and (near-)HIGHEST
self-relatedness site of its year, with mean F_st 3-5.7x the next site. A sample of relatives has
allele frequencies far from its deme's, and that displacement enters F_c as a t-INDEPENDENT
variance -- exactly the shape measured. If that is the cause, the per-pair F_c excess should track
the site's relatedness and run opposite to its pi.

WHAT MAKES THIS ACTIONABLE RATHER THAN INTERESTING. The simulation cannot produce a sample of
relatives -- it draws n_i diploids at random from a deme of hundreds. So any excess attributable
to family structure is a mis-specification the model structurally cannot match, which is precisely
what retired LD (7.5). The difference is that LD's contaminant was spread over every distance bin
while this one is concentrated in identifiable FIELD PAIRS, so it can be measured, and the
leave-one-out section says whether the pooled statistic and its slope survive dropping them.

Reads only empirical files; simulates nothing; runs in seconds.

Run from diagnostics/:  python fc_attribute.py
"""
import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "Python_Code"))
import ABCAnalysisNoRedis as ABC  # noqa: E402  -- the REAL readers and specifier ordering

POOLED = Path("../data/empiricalStats/averaged_temporalFc.csv")
EMP = Path("../data/empiricalStats")
N_PERM = 9999


def site_year(name):
    tail = name.rsplit("-", 1)[-1]
    if not (len(tail) == 4 and tail.isdigit()):
        raise ValueError(f"cannot read a year off the site name {name!r}")
    return tail


def year_tables(year):
    """(pi, self-relatedness, mean off-diagonal F_st) per specifier row, for one year."""
    pi = ABC._read_vector(EMP / f"averaged_pi_{year}.csv")
    rel = ABC._read_matrix(EMP / f"averaged_genRel_{year}.csv")
    fst = ABC._read_matrix(EMP / f"averaged_fst_{year}.csv")
    off = fst.copy()
    np.fill_diagonal(off, np.nan)
    with np.errstate(invalid="ignore"):
        mean_fst = np.nanmean(off, axis=1)
    return pi, np.diag(rel).copy(), mean_fst


def perm_corr(x, y, groups, rng, n_perm=N_PERM):
    """Pearson r with a label-permutation p, permuting WITHIN group.

    The group is the generation gap. Permuting across gaps would let a real t effect masquerade
    as an association with whatever predictor happens to differ between the two field sets, and
    the two sets are largely disjoint (7.9 TEST 4). Same reasoning as 7.1's Mantel test: permute
    the exchangeable unit, which here is the field pair within its own gap.
    """
    x, y, groups = np.asarray(x), np.asarray(y), np.asarray(groups)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y, groups = x[ok], y[ok], groups[ok]
    if len(x) < 4 or np.std(x) == 0 or np.std(y) == 0:
        return float("nan"), float("nan"), len(x)
    obs = float(np.corrcoef(x, y)[0, 1])
    idx_by_g = {g: np.where(groups == g)[0] for g in np.unique(groups)}
    hits = 0
    for _ in range(n_perm):
        yp = y.copy()
        for g, idx in idx_by_g.items():
            yp[idx] = y[rng.permutation(idx)]
        if abs(float(np.corrcoef(x, yp)[0, 1])) >= abs(obs):
            hits += 1
    return obs, (hits + 1) / (n_perm + 1), len(x)


def spearman(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 4:
        return float("nan")

    def rank(v):
        order = np.argsort(v, kind="mergesort")
        r = np.empty(len(v), float)
        r[order] = np.arange(len(v), dtype=float)
        # average ties
        _, inv, cnt = np.unique(v, return_inverse=True, return_counts=True)
        for i, c in enumerate(cnt):
            if c > 1:
                r[inv == i] = r[inv == i].mean()
        return r

    return float(np.corrcoef(rank(x[ok]), rank(y[ok]))[0, 1])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pooled", default=str(POOLED))
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)

    rows = list(csv.DictReader(open(a.pooled, newline="", encoding="utf-8")))
    tables, names = {}, {}
    for r in rows:
        for s in (r["site_a"], r["site_b"]):
            y = site_year(s.strip())
            if y not in tables:
                tables[y] = year_tables(y)
                names[y] = ABC._specifier_site_names(y)

    recs = []
    for r in rows:
        sa, sb = r["site_a"].strip(), r["site_b"].strip()
        vals = {}
        for tag, s in (("a", sa), ("b", sb)):
            y = site_year(s)
            if s not in names[y]:
                raise SystemExit(f"{s} is not a row of specifier_matrix_{y}.csv")
            i = names[y].index(s)
            pi, rel, fst = tables[y]
            vals[tag] = (pi[i], rel[i], fst[i])
        fc, ped = float(r["fc"]), float(r["pedestal"])
        na, nb = int(r["n_a"]), int(r["n_b"])
        recs.append(dict(
            a=sa, b=sb, t=int(r["t"]), n_a=na, n_b=nb, fc=fc, ped=ped, drift=fc - ped,
            pi_mean=np.mean([vals["a"][0], vals["b"][0]]),
            pi_min=min(vals["a"][0], vals["b"][0]),
            rel_mean=np.mean([vals["a"][1], vals["b"][1]]),
            rel_max=max(vals["a"][1], vals["b"][1]),
            fst_mean=np.mean([vals["a"][2], vals["b"][2]]),
            fst_max=max(vals["a"][2], vals["b"][2]),
            harm_n=2.0 / (1.0 / na + 1.0 / nb),
            n_loci=int(r["n_loci"])))

    # ---- the per-pair table, ranked by drift -------------------------------
    print("PER-PAIR DRIFT (F_c - pedestal), ranked. The pedestal decomposition is NOT reliable in")
    print("absolute terms at n=7 (7.9.4 -- the MAF cut conditions on F_c's own denominator), but")
    print("it is a fair CONTRAST across pairs, and most pairs here sit at the identical n=7,7.\n")
    print(f"{'pair':<46} {'t':>3} {'n':>7} {'F_c':>8} {'drift':>8} "
          f"{'selfRel':>9} {'pi':>9} {'meanFst':>8}")
    for r in sorted(recs, key=lambda r: -r["drift"]):
        print(f"{r['a'] + ' / ' + r['b']:<46} {r['t']:>3} {r['n_a']:>3},{r['n_b']:<3} "
              f"{r['fc']:>8.5f} {r['drift']:>8.5f} {r['rel_max']:>9.6f} "
              f"{r['pi_min']:>9.6f} {r['fst_max']:>8.5f}")
    print()

    # ---- correlations ------------------------------------------------------
    print(f"ASSOCIATION with the drift term, {N_PERM} within-gap label permutations")
    print(f"{'predictor':<22} {'Pearson':>9} {'p':>8} {'Spearman':>10} {'n':>4}   prediction")
    preds = [("rel_max", "self-relatedness (max of the two endpoints)", "+ if family structure"),
             ("rel_mean", "self-relatedness (mean)", "+ if family structure"),
             ("pi_min", "pi (min of the two endpoints)", "- if family structure"),
             ("pi_mean", "pi (mean)", "- if family structure"),
             ("fst_max", "site mean F_st (max)", "+ if family structure"),
             ("harm_n", "harmonic sample size", "- if a sampling artifact"),
             ("n_loci", "loci in the pair", "0 -- precision, not level")]
    drift = [r["drift"] for r in recs]
    gaps = [r["t"] for r in recs]
    for key, label, pred in preds:
        x = [r[key] for r in recs]
        rr, p, n = perm_corr(x, drift, gaps, rng)
        print(f"{label:<22.22} {rr:>9.3f} {p:>8.4f} {spearman(x, drift):>10.3f} {n:>4}   {pred}")
    print()

    # ---- who dominates the pooled statistic --------------------------------
    print("LEAVE-ONE-PAIR-OUT on the two quantities the fit actually uses.")
    print("The gap-pooled mean is what _fc_loss differences; the CONTRAST is TEST 4's slope.\n")

    def pooled(sub):
        by = defaultdict(list)
        for r in sub:
            by[r["t"]].append(r["fc"])
        return {t: float(np.mean(v)) for t, v in by.items()}

    base = pooled(recs)
    gs = sorted(base)
    lo, hi = gs[0], gs[-1]
    base_slope = base[hi] - base[lo]
    print(f"{'dropped pair':<46} {'t' + str(lo):>9} {'t' + str(hi):>9} {'contrast':>10} "
          f"{'d contrast':>11}")
    print(f"{'(none)':<46} {base[lo]:>9.6f} {base[hi]:>9.6f} {base_slope:>10.6f} {'':>11}")
    deltas = []
    for r in recs:
        p = pooled([q for q in recs if q is not r])
        deltas.append((abs(p[hi] - p[lo] - base_slope), r, p))
    for d, r, p in sorted(deltas, reverse=True)[:6]:
        print(f"{r['a'] + ' / ' + r['b']:<46} {p[lo]:>9.6f} {p[hi]:>9.6f} "
              f"{p[hi] - p[lo]:>10.6f} {p[hi] - p[lo] - base_slope:>+11.6f}")
    print("   (six most influential of %d; a contrast that flips sign when ONE pair is dropped is"
          % len(recs))
    print("    not a measurement of drift.)")


if __name__ == "__main__":
    main()
