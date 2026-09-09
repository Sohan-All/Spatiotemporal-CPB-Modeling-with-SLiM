"""Between-chromosome tests on the far-field LD excess (CLAUDE.md 7.6.2, 7.5.1 pt 7).

THE QUESTION. Far past the decay, E[r^2] -> the 1/n_hap sampling floor (Hill 1981). Both sides
sit on that floor, but the OBSERVED side carries ~45% more real LD above it than the simulation
does. That excess is either genuine -- a smaller recent Ne than the model produces -- or it is LD
that Beagle imputation manufactured by copying haplotypes between samples. The empirical files are
post-imputation hard calls with no quality field, so the two cannot be separated directly.

These two tests are what the per-chromosome files (data/ld_per_chr/, 17 chr x 3 years) can say
without the pre-imputation VCFs.

  TEST 1 -- is the gap even significant?  Leave-one-chromosome-out jackknife on the pooled
  excess. If the between-chromosome spread is comparable to the sim-obs gap, there is nothing
  to explain.

  TEST 2 -- is the excess homogeneous across chromosomes?  A small recent Ne is a property of the
  POPULATION, so it should raise long-range LD by about the same amount on every chromosome.
  Imputation accuracy depends on local marker density, which varies (5.4: chr6 6.26% vs 8.72%
  genome-wide). Heterogeneity that tracks density points at imputation; homogeneity points at a
  real demographic cause.

NEITHER IS CONCLUSIVE. Local recombination rate also varies between chromosomes and also acts on
the far field, so heterogeneity has an innocent explanation too. Read this as evidence, not proof.

Run from diagnostics/:  python ld_chrom_check.py
"""
import argparse
import csv
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "Python_Code"))
import ABCAnalysisNoRedis as ABC  # noqa: E402  -- the REAL fitted mask
from AnalyzeTreeSeq import _real_sample_sizes  # noqa: E402  -- the REAL sample sizes

YEARS = (2015, 2019, 2023)
N_CHR = 17
PER_CHR = Path("../data/ld_per_chr")
PROBE = Path("../out/ld_probe")
DEFAULT_FAR = 100_000


def read_curves(path):
    rows = list(csv.DictReader(open(path, newline="", encoding="utf-8")))
    lo = np.array([float(r["bin_lo"]) for r in rows])
    r2c = [c for c in rows[0] if c.startswith("r2_")]
    nc = [c for c in rows[0] if c.startswith("n_")]

    def grab(cols):
        return np.array([[float(r[c]) if r[c] not in ("", "None") else np.nan
                          for c in cols] for r in rows])

    return lo, grab(r2c), grab(nc)


def far_field(path, far, keep):
    """Pair-count-weighted pooled r^2 over far-field bins, and the total pair count.

    Pooling is sum(r2*n)/sum(n) -- the same rule ld_loss uses (invariant 4), never a mean
    of means.
    """
    lo, r2, n = read_curves(path)
    rows = lo >= far
    num = np.nansum(r2[rows][:, keep] * n[rows][:, keep])
    den = np.nansum(n[rows][:, keep])
    per_deme = np.nansum(r2[rows][:, keep] * n[rows][:, keep], axis=0) / \
        np.nansum(n[rows][:, keep], axis=0)
    return num / den, den, per_deme


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--far", type=int, default=DEFAULT_FAR,
                    help="bins with bin_lo >= this are treated as plateau")
    ap.add_argument("--sim-pat", default="ld_{y}_pop2000_thin25_s1.csv",
                    help="simulated per-deme curves in out/ld_probe (POPMULT=2000, r=8e-7)")
    a = ap.parse_args()

    print(f"Far field = bins with bin_lo >= {a.far:,}.  Excess = pooled r^2 - 1/n_hap.")
    print("Pooling is sum(r2*n)/sum(n) over the FITTED demes, per chromosome.\n")

    for y in YEARS:
        keep = np.asarray(ABC.get_keep_mask(y), dtype=bool)
        n_hap = 2 * np.asarray(_real_sample_sizes(y), dtype=float)
        # the floor the pooled curve should sit on: pair-weighted mean of 1/n_hap
        floor_deme = 1.0 / n_hap

        vals, pairs = [], []
        for c in range(1, N_CHR + 1):
            p = PER_CHR / f"chr{c}_{y}_ldDecay.csv"
            if not p.exists():
                print(f"  MISSING {p.name}")
                continue
            plat, npair, per_deme = far_field(p, a.far, keep)
            # weight the floor the same way the r2 was weighted
            w = np.nansum(read_curves(p)[2][read_curves(p)[0] >= a.far][:, keep], axis=0)
            fl = float(np.nansum(floor_deme[keep] * w) / np.nansum(w))
            vals.append(plat - fl)
            pairs.append(npair)
        vals = np.array(vals)
        pairs = np.array(pairs, dtype=float)

        # pooled excess, weighted by pair count (what the genome-wide file reports)
        pooled = float(np.sum(vals * pairs) / np.sum(pairs))
        # leave-one-chromosome-out jackknife
        n = len(vals)
        loo = np.array([np.sum(np.delete(vals * pairs, i)) / np.sum(np.delete(pairs, i))
                        for i in range(n)])
        se = float(np.sqrt((n - 1) / n * np.sum((loo - loo.mean()) ** 2)))

        # simulated side, same far field, same mask
        sim_path = PROBE / a.sim_pat.format(y=y)
        sim_txt = "n/a"
        gap = None
        if sim_path.exists():
            splat, _, _ = far_field(sim_path, a.far, keep)
            lo_s, _, n_s = read_curves(sim_path)
            ws = np.nansum(n_s[lo_s >= a.far][:, keep], axis=0)
            sfl = float(np.nansum(floor_deme[keep] * ws) / np.nansum(ws))
            sim_ex = splat - sfl
            gap = pooled - sim_ex
            sim_txt = f"{sim_ex:+.5f}"

        print(f"=== {y} ===")
        print(f"  per-chromosome excess: min {vals.min():+.5f}  max {vals.max():+.5f}  "
              f"sd {vals.std(ddof=1):.5f}")
        print(f"  pooled excess          {pooled:+.5f}   jackknife SE {se:.5f}")
        print(f"  simulated excess       {sim_txt}")
        if gap is not None:
            print(f"  GAP (obs - sim)        {gap:+.5f}   = {gap/se:.1f} jackknife SE"
                  f"   {'-> SIGNIFICANT' if abs(gap) > 2*se else '-> not distinguishable'}")
        # TEST 2: heterogeneity, and does it track the pair-count density proxy?
        het = vals.std(ddof=1) / se
        rho = float(np.corrcoef(np.log(pairs), vals)[0, 1])
        print(f"  heterogeneity: sd/SE = {het:.1f}   "
              f"({'chromosomes differ far more than sampling allows' if het > 3 else 'consistent with homogeneous'})")
        print(f"  corr(excess, log pair count)  {rho:+.3f}   "
              f"(pair count is a rough SNP-density proxy)")
        order = np.argsort(vals)
        print("  lowest 3 chr: " + ", ".join(f"chr{order[i]+1}={vals[order[i]]:+.5f}" for i in range(3)))
        print("  highest 3 chr: " + ", ".join(f"chr{order[-1-i]+1}={vals[order[-1-i]]:+.5f}" for i in range(3)))
        print()




# ---------------------------------------------------------------- TEST 3
def test_shape(far=DEFAULT_FAR, ref_lo=100.0, frac=0.25):
    """Is a chromosome's extra far-field LD a DISTANCE SHIFT or a LEVEL SHIFT?

    The two candidate causes make different predictions, and this is what separates them:

      LOW RECOMBINATION rescales DISTANCE. Under Sved, r^2 depends on d only through 4*Ne*r*d,
      so halving r is the same as halving d -- the whole curve slides right, keeping its shape.
      A chromosome with more far-field excess must therefore ALSO decay more slowly.

      IMPUTATION adds a LEVEL. Copying haplotypes between samples injects LD without changing
      how fast real associations break down, so the far-field sits higher while the decay
      position stays put.

    So: measure each chromosome's decay distance (where the excess-over-floor falls to `frac` of
    its value at `ref_lo` bp) and correlate it against its far-field excess. Strong positive
    correlation => recombination. No relationship => the level moved on its own => imputation.
    """
    print("=" * 74)
    print("TEST 3: DISTANCE SHIFT (recombination) or LEVEL SHIFT (imputation)?")
    print("=" * 74)
    print(f"  decay distance = where excess-over-floor falls to {frac:.0%} of its {ref_lo:.0f} bp value")
    print("  Strong positive corr with far-field excess => recombination (one cause, both effects).")
    print("  No corr => the far-field level moved independently of decay => imputation.\n")

    for y in YEARS:
        keep = np.asarray(ABC.get_keep_mask(y), dtype=bool)
        floor_deme = 1.0 / (2 * np.asarray(_real_sample_sizes(y), dtype=float))
        dist, exc = [], []
        for c in range(1, N_CHR + 1):
            p = PER_CHR / f"chr{c}_{y}_ldDecay.csv"
            if not p.exists():
                continue
            lo, r2, n = read_curves(p)
            mid = np.array([float(r["bin_mid"]) for r in
                            csv.DictReader(open(p, newline="", encoding="utf-8"))])
            w = n[:, keep]
            pooled = np.nansum(r2[:, keep] * w, axis=1) / np.nansum(w, axis=1)
            fl = np.nansum(floor_deme[keep] * w, axis=1) / np.nansum(w, axis=1)
            e = pooled - fl

            i0 = int(np.flatnonzero(lo == ref_lo)[0])
            ref = e[i0]
            target = frac * ref
            d = np.nan
            for i in range(i0, len(e) - 1):
                if e[i] >= target > e[i + 1]:
                    # log-linear interpolation in distance
                    t = (e[i] - target) / (e[i] - e[i + 1])
                    d = float(np.exp(np.log(mid[i]) + t * (np.log(mid[i + 1]) - np.log(mid[i]))))
                    break
            rows = lo >= far
            ff = float(np.nansum(r2[rows][:, keep] * n[rows][:, keep]) /
                       np.nansum(n[rows][:, keep])
                       - np.nansum(floor_deme[keep] * np.nansum(n[rows][:, keep], axis=0)) /
                       np.nansum(n[rows][:, keep]))
            dist.append(d)
            exc.append(ff)
        dist, exc = np.array(dist), np.array(exc)
        ok = np.isfinite(dist) & np.isfinite(exc)
        rho = float(np.corrcoef(np.log(dist[ok]), exc[ok])[0, 1])
        order = np.argsort(exc)
        print(f"=== {y} ===   (n={ok.sum()} chromosomes)")
        print(f"  decay distance: min {np.nanmin(dist):,.0f} bp   max {np.nanmax(dist):,.0f} bp   "
              f"ratio {np.nanmax(dist)/np.nanmin(dist):.2f}x")
        print(f"  corr(log decay distance, far-field excess) = {rho:+.3f}")
        print("   chr  far-field excess   decay distance")
        for i in list(order[:3]) + list(order[-3:]):
            print(f"   {i+1:>3}  {exc[i]:+.5f}          {dist[i]:>10,.0f} bp")
        print()


if __name__ == "__main__":
    main()
    test_shape()
