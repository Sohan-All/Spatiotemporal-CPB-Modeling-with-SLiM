"""
Offline standardization + ranking for the rejection-ABC pass (CLAUDE.md 7).

Runs AFTER the big pass, on the accumulated results CSV. It does NOT run any simulations.

Offline because sigma can only be measured from the spread of losses across the whole run set,
which is its own pilot batch. This script:

  1. reads the results CSV,
  2. computes sigma_j = 1.4826 * MAD for each FITTED statistic,
  3. combines them into one standardized distance D per run,
  4. writes a ranked CSV and freezes the sigmas to a JSON file.

FITTED (enter D):      pi_loss (log-space), fst_loss, fc_loss   -- from batch 3 (2026-09-12)
DIAGNOSTIC (not in D): ibd_loss, dxy_loss, genrel_loss, ld_loss -- posterior-predictive checks only
                       (ld_loss is retired as a fitted statistic, CLAUDE.md 7.5).
"""
import json
import numpy as np
import pandas as pd
from pathlib import Path

# ----------------------------- CONFIG -----------------------------
# Pooled batch results live under out/<batch>/, NOT at ../out/abc_results.csv -- that path is
# where each CHTC job writes its own output, so a tracked file there is cloned into every job
# and appended to (CLAUDE.md 7.6.0). Point these at the batch you mean to standardize.
RESULTS_CSV   = Path("../out/batch3/abc_results.csv")          # input: the pass results
RANKED_CSV    = Path("../out/batch3/abc_results_ranked.csv")   # output: results + D, sorted
SIGMAS_JSON   = Path("../out/batch3/abc_sigmas.json")          # output: frozen sigmas
FITTED_STATS  = ["pi_loss", "fst_loss", "fc_loss"]      # statistics that enter the distance D

# BATCH 3 WEIGHTS: NOT SET YET, deliberately. They come from batch 3's own variance decomposition
# (python ../diagnostics/collect_batch.py, section 9 prints them), by the 7.4.1 rule: weight each
# fitted statistic by the share of its batch spread that is DEMOGRAPHIC signal -- the unique
# rank-space R2 of pop + total_migration + m + numClusters. Never guessed, and the batch 1/2 values
# below do NOT carry over: both were fitted with mu FREE, and mu is fixed from batch 3
# (scale_constants.py; CLAUDE.md 7.9.3), which removes pi_loss's main nuisance driver.
# main() RAISES while this is None, rather than falling back to equal weights -- equal weights were
# ruled out on both earlier batches, so a silent fallback would be a known-wrong D.
WEIGHTS       = None   # e.g. {"pi_loss": ..., "fst_loss": ..., "fc_loss": ...} from collect_batch

# HISTORY -- batch 1 (2,495 trials, mu free), results in ../out/batch1/:
#   WEIGHTS = {"pi_loss": 0.125, "fst_loss": 0.875}, FITTED_STATS = ["pi_loss", "fst_loss"]
#   statistic   R2_total   demographic   mu-nuisance   unexplained
#   pi_loss        0.240        0.0893        0.1533         0.760
#   fst_loss       0.622        0.6243        0.0004         0.378
# pi_loss's spread was MAJORITY mu-draw (63% of what the parameters explained). The 7.3
# replicate-noise rule gave near-equal 0.488/0.512 and was rejected, because mu-driven spread is
# not replicate noise and so counted as signal under it. To re-rank batch 1, restore these values
# AND the batch1 paths above.
ACCEPT_FRAC   = 0.20            # fraction of runs to flag as 'accepted' (top by smallest D)
# ------------------------------------------------------------------


def robust_sigma(values):
    """sigma = 1.4826 * MAD (median absolute deviation), ignoring NaN. Robust to the
    2-individual-subpop outliers (plan §1). Returns np.nan if <2 finite values or MAD==0."""
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    if v.size < 2:
        return np.nan
    mad = np.median(np.abs(v - np.median(v)))
    return 1.4826 * mad if mad > 0 else np.nan


def main():
    # Config checks first, so a half-edited config is reported as such rather than as missing data.
    if WEIGHTS is None:
        raise ValueError(
            "WEIGHTS is not set. Run diagnostics/collect_batch.py on this batch and copy the "
            "WEIGHTS its section 9 prints (CLAUDE.md 7.4.1). Equal weights are deliberately not "
            "a fallback.")
    if set(WEIGHTS) != set(FITTED_STATS):
        raise ValueError(f"WEIGHTS keys {sorted(WEIGHTS)} != FITTED_STATS {sorted(FITTED_STATS)} "
                         f"-- a statistic without a weight, or a weight without a statistic, "
                         f"means the config was only half updated.")
    if not RESULTS_CSV.exists():
        raise FileNotFoundError(f"Results CSV not found: {RESULTS_CSV.resolve()} "
                                f"(run the pass first).")

    df = pd.read_csv(RESULTS_CSV)

    missing = [s for s in FITTED_STATS if s not in df.columns]
    if missing:
        raise ValueError(f"Results CSV missing fitted-stat columns: {missing}. "
                         f"Columns present: {list(df.columns)}")

    # 1) robust per-statistic scale, frozen
    sigmas = {}
    for stat in FITTED_STATS:
        s = robust_sigma(df[stat].values)
        if not np.isfinite(s):
            raise ValueError(f"Could not compute a finite sigma for '{stat}' "
                             f"(all-NaN, <2 values, or zero MAD). Inspect the pass output.")
        sigmas[stat] = float(s)

    # 2) weights (validated non-None and matching FITTED_STATS above)
    total = float(sum(WEIGHTS[s] for s in FITTED_STATS))
    w = {stat: WEIGHTS[stat] / total for stat in FITTED_STATS}

    # 3) combined standardized distance:  D = sqrt( sum_j w_j * (loss_j / sigma_j)^2 )
    sq = np.zeros(len(df))
    for stat in FITTED_STATS:
        z = df[stat].values / sigmas[stat]
        sq = sq + w[stat] * np.square(z)
    df["D"] = np.sqrt(sq)

    # 4) rank + flag acceptance (smaller D = better). NaN D sinks to the bottom.
    df = df.sort_values("D", kind="mergesort", na_position="last").reset_index(drop=True)
    n_accept = max(1, int(np.floor(ACCEPT_FRAC * np.isfinite(df["D"]).sum())))
    df["accepted"] = False
    df.loc[df.index[:n_accept], "accepted"] = np.isfinite(df["D"].iloc[:n_accept]).values

    # 5) write outputs
    RANKED_CSV.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(RANKED_CSV, index=False)
    SIGMAS_JSON.write_text(json.dumps(
        {"sigmas": sigmas, "weights": w, "fitted_stats": FITTED_STATS,
         "accept_frac": ACCEPT_FRAC, "n_runs": int(len(df)), "n_accepted": int(n_accept)},
        indent=2))

    # 6) report
    print(f"Runs: {len(df)}  |  finite D: {int(np.isfinite(df['D']).sum())}")
    print("Frozen sigmas (1.4826*MAD):")
    for stat in FITTED_STATS:
        print(f"  {stat:10s} sigma={sigmas[stat]:.6g}  weight={w[stat]:.3f}")
    print(f"Accepted (top {ACCEPT_FRAC:.0%} by D): {n_accept} runs")
    print(f"  best D = {df['D'].iloc[0]:.4g}   acceptance threshold eps = {df['D'].iloc[n_accept-1]:.4g}")
    print(f"Wrote: {RANKED_CSV}")
    print(f"Wrote: {SIGMAS_JSON}")
    print("\nNOTE: keep eps above the replicate noise floor (CLAUDE.md 7.9.10C; fc_loss ~0.0018-0.0022, "
          "fst_loss ~0.0003, pi_loss ~0.0017-0.011 by POPMULT). Below it you are selecting on "
          "coalescent noise, and with a misspecified model that biases toward noisy prior regions.")


if __name__ == "__main__":
    main()
