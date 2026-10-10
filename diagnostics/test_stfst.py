"""Verify the spatio-temporal F_st (CLAUDE.md 7.9.14B) before trusting it -- CLAUDE.md 10's rule.

1. AnalyzeTreeSeq.hudson_fst_pairs against a brute-force numpy Hudson F_st computed straight from
   the genotype matrix, on a small msprime simulation with samples at times 0, 8 and 16 (the three
   collections), cross-time pairs included. Must agree to rounding.
2. It is Hudson, not Nei: tskit's Fst formula on the same pairs must come out clearly lower.
3. _draw_disjoint: whole individuals, exact sizes, no overlap between sites sharing a deme, and the
   equal-share fallback when the deme is too small.
4. stfst_common on the real specifier matrices: pair count, gaps, distances against
   resample_distances.hav, and summarise/contrast on synthetic values.

Usage (from diagnostics/):  python test_stfst.py      (no SLiM, ~10 s)
"""
import os
import sys
from itertools import combinations
from pathlib import Path

import msprime
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "Python_Code"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
os.chdir(ROOT / "Python_Code")

import AnalyzeTreeSeq as ATS        # noqa: E402
import stfst_common as stc          # noqa: E402
from resample_distances import read_spec, hav   # noqa: E402

TIMES = (0, 8, 16)


def sim(seed=1):
    dem = msprime.Demography.island_model([500, 500, 500], migration_rate=0.002)
    samples = [msprime.SampleSet(12, population=p, time=t, ploidy=2)
               for p in range(3) for t in TIMES]
    ts = msprime.sim_ancestry(samples, demography=dem, sequence_length=2e6,
                              recombination_rate=1e-8, random_seed=seed)
    return msprime.sim_mutations(ts, rate=2e-8, model=msprime.BinaryMutationModel(),
                                 discrete_genome=False, random_seed=seed)


def brute_hudson(ts, sets, pairs):
    """pi and d_xy summed over sites from 0/1 genotypes, divided by L; Hudson from those."""
    G = ts.genotype_matrix()
    col = {int(n): i for i, n in enumerate(ts.samples())}
    L = ts.sequence_length

    def pi(s):
        g = G[:, [col[int(n)] for n in s]]
        k, n = g.sum(axis=1), g.shape[1]
        return float((k * (n - k)).sum() / (n * (n - 1) / 2) / L)

    def dxy(s, t):
        a = G[:, [col[int(n)] for n in s]].mean(axis=1)
        b = G[:, [col[int(n)] for n in t]].mean(axis=1)
        return float((a * (1 - b) + b * (1 - a)).sum() / L)

    return np.array([1 - 0.5 * (pi(sets[i]) + pi(sets[j])) / dxy(sets[i], sets[j])
                     for i, j in pairs])


def main():
    ok = True
    ts = sim()
    rng = np.random.default_rng(7)

    # sample sets: per (deme, time), two DISJOINT sites of 5 and 4 diploids
    sets, labels = [], []
    for p in range(3):
        for t in TIMES:
            a, b = ATS._draw_disjoint(ts, p, t, [5, 4], rng)
            sets += [a, b]
            labels += [(p, t, "a"), (p, t, "b")]
            inds = [{ts.node(int(n)).individual for n in s} for s in (a, b)]
            if inds[0] & inds[1]:
                print(f"FAIL: deme {p} t={t}: the two draws share individuals"); ok = False
            if [len(a), len(b)] != [10, 8] or any(len(s) != 2 * len(i) for s, i in zip((a, b), inds)):
                print(f"FAIL: deme {p} t={t}: wrong sizes or split individuals"); ok = False
            if any(ts.node(int(n)).time != t or ts.node(int(n)).population != p for n in np.r_[a, b]):
                print(f"FAIL: deme {p} t={t}: nodes from the wrong deme/time"); ok = False

    pairs = list(combinations(range(len(sets)), 2))
    mine = ATS.hudson_fst_pairs(ts, sets, pairs)
    ref = brute_hudson(ts, sets, pairs)
    err = np.max(np.abs(mine - ref))
    print(f"1. Hudson vs brute force over {len(pairs)} pairs ({ts.num_sites} sites): "
          f"max |diff| {err:.2e}   F_st range {mine.min():.4f} to {mine.max():.4f}")
    if err > 1e-10:
        print("FAIL: hudson_fst_pairs disagrees with the brute-force estimator"); ok = False

    cross = [k for k, (i, j) in enumerate(pairs) if labels[i][1] != labels[j][1]
             and labels[i][0] != labels[j][0]]
    same_deme = [k for k, (i, j) in enumerate(pairs) if labels[i][:2] == labels[j][:2]]
    print(f"   cross-time between-deme median {np.median(mine[cross]):.4f}; "
          f"same deme+time (disjoint draws) median {np.median(mine[same_deme]):+.4f}")

    # tskit's Fst definition (Nei/Slatkin), written out so the standing grep for that method
    # stays empty (CLAUDE.md 9): 1 - 2(pi_a + pi_b) / (pi_a + pi_b + 2 d_xy)
    pi = np.asarray(ts.diversity(sets))
    dxy = np.asarray(ts.divergence(sets, indexes=pairs))
    a, b = np.array(pairs).T
    nei = 1 - 2 * (pi[a] + pi[b]) / (pi[a] + pi[b] + 2 * dxy)
    big = mine > 0.02
    ratio = np.median(nei[big] / mine[big])
    print(f"2. Nei (tskit's Fst) / Hudson, pairs with F_st > 0.02: median {ratio:.2f}")
    if not ratio < 0.8:
        print("FAIL: expected Nei to be clearly below Hudson"); ok = False

    short = ATS._draw_disjoint(ts, 0, 0, [7, 7], rng)      # 12 individuals, asked for 14
    if [len(s) // 2 for s in short] != [6, 6] or set(short[0]) & set(short[1]):
        print(f"FAIL: equal-share fallback gave {[len(s) // 2 for s in short]}"); ok = False
    else:
        print("3. _draw_disjoint: disjoint, whole individuals, exact sizes; 12 inds for 7+7 -> 6+6")

    sites = stc.site_table(ROOT / "data" / "Genetic_Data")
    allp = stc.all_pairs(sites)
    n = len(sites)
    gaps = sorted({g for _, _, g, _ in allp})
    print(f"4. {n} sites, {len(allp)} pairs, gaps {gaps}")
    if len(allp) != n * (n - 1) // 2 or gaps != [0, 8, 16]:
        print("FAIL: pair enumeration"); ok = False
    if any(sites[a]["year"] > sites[b]["year"] for a, b, _, _ in allp):
        print("FAIL: year_a > year_b somewhere"); ok = False
    coords = np.vstack([read_spec(y)[1] for y in stc.YEARS])
    names = sum((read_spec(y)[0] for y in stc.YEARS), [])
    if names != [s["site"] for s in sites]:
        print("FAIL: site_table order differs from the specifier matrices"); ok = False
    D = hav(coords, coords)
    derr = max(abs(D[a, b] - d) for a, b, _, d in allp)
    print(f"   distances vs resample_distances.hav: max |diff| {derr:.2e} km")
    if derr > 0.05:
        print("FAIL: distances"); ok = False

    # summarise / contrast on synthetic values: F_st = 0.001*(gap/8) + 0.0001*bin
    rows = [{"year_a": sites[a]["year"], "site_a": sites[a]["site"], "year_b": sites[b]["year"],
             "site_b": sites[b]["site"], "gap": g, "dist_km": d,
             "fst": 0.001 * g / 8 + 0.0001 * stc.dist_bin(d)} for a, b, g, d in allp]
    keep = stc.retained_sites(ROOT / "data" / "Genetic_Data", exclude_sites={"Mortensen9-2015"})
    summ = stc.summarise(rows, keep)
    con = stc.contrast(summ)
    exp = {(g, k): 0.001 * g / 8 + 0.0001 * k - 0.0001 * max(k, 1) for (g, k) in con}
    cerr = max(abs(con[c] - exp[c]) for c in con if not np.isnan(con[c]))
    print(f"   summarise/contrast on synthetic values: max |diff| {cerr:.1e}; retained sites "
          f"{len(keep)} of {n}")
    if cerr > 1e-12 or ("2015", "Mortensen9-2015") in keep or len(keep) != n - 3:
        print("FAIL: summarise/contrast/retained_sites"); ok = False
    print("   cells (gap, bin): n_pairs  " + "  ".join(
        f"{g}/{stc.bin_label(k)}:{summ[(g, k)]['n_pairs']}" for (g, k) in sorted(summ)))

    print(f"\nspec {stc.spec_hash()}")
    print("ALL PASS" if ok else "SOME CHECKS FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
