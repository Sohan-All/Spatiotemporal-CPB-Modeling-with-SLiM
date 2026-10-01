"""Is the persistent-deme drift misfit a SCALE artifact? (TODO 2026-09-29 step 2, CLAUDE.md 7.9.14A)

Reads batch_signed_features.py output and, for the OFF arm (persistent demes, refound_k = -1),
asks how drift growth and the F_c excess levels move with deme size, and where they extrapolate.
If growth heads to ~0 (or the observed value) at deme sizes above the prior, the misfit is the
POPMULT ceiling, not model structure.

Mean deme size = total N / demes = 3.33 * POPMULT / (33 * numClusters)   (CLAUDE.md 3, 1.1).

Observed values and their field-pair SEs come from fc_growth_jackknife.py:
  growth -0.00151 (SE 0.00279), excess t=8 0.02832 (0.00190), excess t=16 0.02681 (0.00205).

Usage (from diagnostics/):  python offarm_scale.py --features ../out/batch6/signed_features.csv
"""
import argparse

import numpy as np
import pandas as pd

OBS = {"growth": (-0.00151, 0.00279), "ex8": (0.02832, 0.00190), "ex16": (0.02681, 0.00205)}
OBS_PAIR_SD = {"exsd8": 0.00600, "exsd16": 0.00614}


def fit_inv(x, y):
    """y = a + b/x by least squares; bootstrap CI on a (the infinite-deme asymptote)."""
    X = np.column_stack([np.ones(len(x)), 1.0 / x])
    beta = np.linalg.lstsq(X, y, rcond=None)[0]
    rng = np.random.default_rng(1)
    boots = []
    for _ in range(2000):
        i = rng.integers(0, len(x), len(x))
        boots.append(np.linalg.lstsq(X[i], y[i], rcond=None)[0])
    boots = np.array(boots)
    return beta, np.percentile(boots, [2.5, 97.5], axis=0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--features", required=True)
    args = ap.parse_args()
    pd.set_option("display.width", 200)
    d = pd.read_csv(args.features)
    d["deme"] = 3.33 * d["pop"] / (33 * d["numClusters"])
    off = d[d.refound_k < 0].copy()
    print(f"OFF arm: {len(off)} trials; mean deme size {off.deme.min():.0f}-{off.deme.max():.0f} "
          f"(median {off.deme.median():.0f})")

    # 1. binned medians
    off["dbin"] = pd.qcut(off.deme, 6)
    cols = ["growth", "ex8", "ex16", "exsd8", "exsd16"]
    tab = off.groupby("dbin", observed=True)[cols].median()
    tab["n"] = off.groupby("dbin", observed=True).size()
    print("\n1. medians by deme-size sextile (observed: growth -0.00151, ex8 0.02832, ex16 0.02681,"
          " pair sd 0.0060/0.0061)")
    print(tab.round(5).to_string())

    # 2. y = a + b/deme, whole arm and slow kernel
    print("\n2. fit y = a + b/deme  (a = the value at infinite deme size)")
    for label, sub in [("all OFF", off), ("OFF, m < 2e-4", off[off.m < 2e-4])]:
        for c in ["growth", "ex8", "ex16"]:
            (a, b), ci = fit_inv(sub.deme.values, sub[c].values)
            o, se = OBS[c]
            print(f"  {label:14s} {c:6s} a = {a:+.5f} [{ci[0, 0]:+.5f}, {ci[1, 0]:+.5f}]  b = {b:+.3f}"
                  f"   observed {o:+.5f} +/- {se:.5f}")

    # 3. same fit within total_migration terciles (Nm, not N, sets drift in a connected deme)
    print("\n3. growth asymptote by total_migration tercile")
    off["tmq"] = pd.qcut(off.total_migration, 3)
    for q, sub in off.groupby("tmq", observed=True):
        (a, b), ci = fit_inv(sub.deme.values, sub.growth.values)
        print(f"  tm {str(q):22s} n={len(sub):3d}  growth a = {a:+.5f} [{ci[0, 0]:+.5f}, {ci[1, 0]:+.5f}]"
              f"  b = {b:+.3f}   median growth {sub.growth.median():+.5f}")

    # 4. where does the observed growth sit in the simulated distribution, by deme size?
    print("\n4. share of OFF trials with growth <= observed, and within the observed 95% band "
          "[-0.0070, +0.0040]")
    for q, sub in off.groupby("dbin", observed=True):
        print(f"  deme {str(q):22s}  <= obs: {np.mean(sub.growth <= -0.00151):.2f}   "
              f"in band: {np.mean((sub.growth > -0.0070) & (sub.growth < 0.0040)):.2f}")

    # 5. simulated between-pair spread vs observed
    print("\n5. per-pair excess sd within a gap (sim median over OFF trials vs observed)")
    for c in ["exsd8", "exsd16"]:
        print(f"  {c}: sim median {off[c].median():.5f} (IQR {off[c].quantile(.25):.5f}-"
              f"{off[c].quantile(.75):.5f})   observed {OBS_PAIR_SD[c]:.5f}")


if __name__ == "__main__":
    main()
