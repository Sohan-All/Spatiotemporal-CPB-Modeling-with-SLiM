"""
7.4.3 conditioning analysis, as a script (it was ad hoc the first time).

QUESTION. 7.4.2 measured that the pass identifies the dispersal kernel `m` and leaves N at the
prior. 7.4.3 answered, on batch 1, what happens if you PIN the nuisance parameters instead:
`pop` tightens somewhat, but the inferred N becomes a direct function of what you pinned it to.

This script re-asks that on any pooled batch, and -- the reason it exists -- lets you ask it under
TWO weightings at once, so you can see what ADDING a statistic to D actually buys. No simulation;
runs in seconds.

WHAT IT DOES NOT DO. It cannot tell you the right value to pin. That is an external question
(TODO 3). The output is a sensitivity surface, not a posterior.

  python condition_slices.py --results ../out/abc_results_pilot.csv
"""
import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "Python_Code"))
from abc_standardize import robust_sigma  # noqa: E402  -- the EXACT sigma the pass freezes

# Weightings to compare. Both come from collect_batch.py 9 (demographic-R2 rule, 7.4.1);
# BATCH1 is what abc_standardize.py currently holds, PILOT adds ld_loss.
WEIGHTINGS = {
    "pi+fst  (batch-1 weights)": {"pi_loss": 0.125, "fst_loss": 0.875},
    "pi+fst+ld (pilot weights)": {"pi_loss": 0.110, "fst_loss": 0.436, "ld_loss": 0.455},
}

# Subpop size = mean(Average Count) * POPMULT / numSubpops, with numSubpops = 33*numClusters
# and mean(Average Count) ~ 3.33 (CLAUDE.md 3.1). Nm uses the DEME size, not total N.
MEAN_AVG_COUNT = 3.33
CLUSTER_MULT = 33


def deme_n(pop, n_clusters):
    return MEAN_AVG_COUNT * pop / (CLUSTER_MULT * n_clusters)


def load(path):
    import csv
    with open(path, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        raise SystemExit(f"{path}: no rows")
    cols = {}
    for k in rows[0]:
        try:
            cols[k] = np.array([float(r[k]) for r in rows])
        except (ValueError, TypeError):
            continue
    n = len(rows)
    print(f"Loaded {n} rows from {path}")
    missing = [c for c in ("pop", "m", "total_migration", "numClusters") if c not in cols]
    if missing:
        raise SystemExit(f"missing columns: {missing}")
    cols["Nm"] = deme_n(cols["pop"], cols["numClusters"]) * cols["total_migration"]
    return cols, n


def distance(A, weights, sigmas):
    """D = sqrt( sum_j w_j (loss_j/sigma_j)^2 ), weights renormalised to 1."""
    tot = float(sum(weights.values()))
    sq = np.zeros(len(A["pop"]))
    for stat, w in weights.items():
        sq += (w / tot) * np.square(A[stat] / sigmas[stat])
    return np.sqrt(sq)


def iqr(v):
    v = v[np.isfinite(v)]
    if v.size < 4:
        return float("nan")
    return float(np.percentile(v, 75) - np.percentile(v, 25))


def slice_row(A, mask, D, accept, prior_iqr):
    """Rank WITHIN the slice -- this emulates a batch actually run at those pinned values."""
    idx = np.flatnonzero(mask)
    n = idx.size
    if n < 12:
        return n, None
    d = D[idx]
    keep = idx[np.argsort(d, kind="mergesort")[:max(4, int(round(accept * n)))]]
    out = {"n": n, "n_acc": keep.size}
    for c in ("pop", "Nm", "m", "total_migration"):
        out[c + "_med"] = float(np.median(A[c][keep]))
        out[c + "_iqr"] = iqr(A[c][keep]) / prior_iqr[c]
    return n, out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default="../out/batch2/abc_results_pilot.csv")
    ap.add_argument("--accept", type=float, default=0.20, help="acceptance fraction within a slice")
    ap.add_argument("--min-n", type=int, default=12, help="skip slices thinner than this")
    args = ap.parse_args()

    A, n_rows = load(args.results)
    prior_iqr = {c: iqr(A[c]) for c in ("pop", "Nm", "m", "total_migration")}
    print(f"\nPrior IQRs (the batch's own draws, n={n_rows}):")
    for c, v in prior_iqr.items():
        print(f"  {c:<16s} {v:.6g}")
    print("\nAn IQR ratio of 1.0 means the data said NOTHING about that parameter.")
    print(f"Ranking is WITHIN each slice, top {args.accept:.0%} by D.")

    m, tm, nc = A["m"], A["total_migration"], A["numClusters"]
    all_true = np.ones(n_rows, dtype=bool)
    slices = [("(unconditioned)", all_true)]
    for lo, hi in [(3e-5, 5e-5), (5e-5, 8e-5), (8e-5, 1e-4), (1e-4, 3e-4)]:
        slices.append((f"m in [{lo:.0e}, {hi:.0e}]", (m >= lo) & (m < hi)))
    band = (m >= 3e-5) & (m < 8e-5)
    for lo, hi in [(0.02, 0.10), (0.10, 0.20), (0.20, 0.31)]:
        slices.append((f"m[3e-5,8e-5] x tm[{lo:.2f},{hi:.2f}]", band & (tm >= lo) & (tm < hi)))
    for k in (1, 2, 3):
        slices.append((f"m[3e-5,8e-5] x tm[.02,.20] x nC={k}",
                       band & (tm >= 0.02) & (tm < 0.20) & (nc == k)))

    for label, weights in WEIGHTINGS.items():
        missing = [s for s in weights if s not in A]
        if missing:
            print(f"\n--- SKIPPED {label}: missing {missing} ---")
            continue
        sigmas = {s: robust_sigma(A[s]) for s in weights}
        D = distance(A, weights, sigmas)
        print("\n" + "=" * 100)
        print(f"D = {label}")
        print("  sigma: " + "  ".join(f"{s}={sigmas[s]:.6g}" for s in weights))
        print("=" * 100)
        print(f"  {'slice':<36s} {'n':>4s} {'pop med':>9s} {'pop IQR':>8s} "
              f"{'Nm med':>8s} {'Nm IQR':>7s} {'m med':>10s} {'tm med':>7s}")
        for lbl, mask in slices:
            n, r = slice_row(A, mask, D, args.accept, prior_iqr)
            if r is None:
                print(f"  {lbl:<36s} {n:>4d}   -- too thin (<{args.min_n}) --")
                continue
            print(f"  {lbl:<36s} {n:>4d} {r['pop_med']:>9.0f} {r['pop_iqr']:>8.2f} "
                  f"{r['Nm_med']:>8.1f} {r['Nm_iqr']:>7.2f} {r['m_med']:>10.2e} "
                  f"{r['total_migration_med']:>7.3f}")

    print("\nRead the pop column DOWN, not across: if median pop swings with the slice, that is")
    print("the answer changing with the assumption, not an error bar (7.4.3).")


if __name__ == "__main__":
    main()
