"""Uncertainty on the OBSERVED temporal-F_c levels and drift growth (TODO 2026-09-29 step 1).

The observed growth (-0.00151) has never had an error bar, and CLAUDE.md 7.9.14A's reading -- that
the data show no detectable drift -- turns on it. Two resampling schemes, because they see
different noise:

  1. Delete-one-CHROMOSOME jackknife over data/Fc_per_chr/ (17 chromosomes). Captures genomic
     (coalescent + linkage) noise given the sampled individuals. It does NOT see individual
     sampling that is shared across the genome -- relatedness / family structure -- because the
     same beetles are on every chromosome. Both the plain (Tukey) and the loci-weighted delete-m
     (Busing et al. 1999) versions are printed, since chromosomes differ in size.
  2. Resampling FIELD PAIRS within each gap (10 at t=8, 9 at t=16). Sees everything that varies
     between pairs -- genomic noise, individual sampling and family structure, and real
     between-field differences. This is the wider, more honest bar.

Pooling matches ABCAnalysisNoRedis._fc_loss / fc_common: per pair F_c = sum_fc / n_loci over the
genome; per gap = unweighted mean over pairs. Growth = (F16 - ped16) - (F8 - ped8), pedestals from
the target file (they depend on n only, so they are constants here).

Usage (from diagnostics/):  python fc_growth_jackknife.py
"""
import csv
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
PER_CHR = ROOT / "data" / "Fc_per_chr"
TARGET = ROOT / "data" / "empiricalStats" / "averaged_temporalFc.csv"
N_BOOT = 20000


def read_target():
    rows = list(csv.DictReader(open(TARGET, newline="", encoding="utf-8")))
    return {(r["site_a"], r["site_b"]): (int(r["t"]), float(r["fc"]), float(r["pedestal"]),
                                         int(r["n_loci"])) for r in rows}


def read_per_chr():
    """{chrom: {pair: (sum_fc, n_loci)}}"""
    out = {}
    for p in sorted(PER_CHR.glob("chr*_temporalFc.csv")):
        c = int(p.stem.split("_")[0][3:])
        out[c] = {(r["site_a"], r["site_b"]): (float(r["sum_fc"]), int(r["n_loci"]))
                  for r in csv.DictReader(open(p, newline="", encoding="utf-8"))}
    return out


def stats_from(per_chr, chroms, target):
    """(F8, F16, growth) pooled over the given chromosomes."""
    by_gap = {8: [], 16: []}
    for key, (t, _, ped, _) in target.items():
        s = sum(per_chr[c][key][0] for c in chroms)
        n = sum(per_chr[c][key][1] for c in chroms)
        by_gap[t].append((s / n, ped))
    f8 = np.mean([f for f, _ in by_gap[8]]); p8 = np.mean([p for _, p in by_gap[8]])
    f16 = np.mean([f for f, _ in by_gap[16]]); p16 = np.mean([p for _, p in by_gap[16]])
    return np.array([f8, f16, (f16 - p16) - (f8 - p8), f8 - p8, f16 - p16])


NAMES = ["F_c(t=8)", "F_c(t=16)", "growth", "excess t=8", "excess t=16"]


def main():
    target = read_target()
    per_chr = read_per_chr()
    chroms = sorted(per_chr)
    print(f"{len(chroms)} chromosomes, {len(target)} field pairs")

    # --- guard: the per-chromosome files must re-pool EXACTLY to the fitted target, else they
    # are from a different spec (e.g. before the 7.2.2H sample exclusions) and this means nothing.
    worst = 0.0
    for key, (t, fc, _, nl) in target.items():
        s = sum(per_chr[c][key][0] for c in chroms)
        n = sum(per_chr[c][key][1] for c in chroms)
        if n != nl:
            raise SystemExit(f"{key}: per-chr n_loci {n} != target {nl} -- different spec, stop.")
        worst = max(worst, abs(s / n - fc))
    if worst > 1e-8:
        raise SystemExit(f"per-chr files re-pool to F_c off by {worst:.2e} -- different spec, stop.")
    print(f"re-pool check: per-chr files reproduce the target (max |dF_c| {worst:.1e}, n_loci exact)")

    full = stats_from(per_chr, chroms, target)
    g = len(chroms)
    loo = np.array([stats_from(per_chr, [c for c in chroms if c != d], target) for d in chroms])
    # Tukey delete-one
    se_plain = np.sqrt((g - 1) / g * ((loo - loo.mean(0)) ** 2).sum(0))
    # Busing et al. 1999 weighted delete-m, weights = share of loci on the chromosome
    m = np.array([sum(per_chr[c][k][1] for k in target) for c in chroms], dtype=float)
    h = m.sum() / m
    pseudo = h[:, None] * full - (h[:, None] - 1) * loo
    theta_w = g * full - ((1 - m / m.sum())[:, None] * loo).sum(0)
    se_w = np.sqrt((((pseudo - theta_w) ** 2) / (h[:, None] - 1)).sum(0) / g)

    # --- field-pair bootstrap within gap
    rng = np.random.default_rng(20260929)
    pairs = {8: [], 16: []}
    for key, (t, fc, ped, _) in target.items():
        pairs[t].append(fc - ped)
    e8, e16 = np.array(pairs[8]), np.array(pairs[16])
    b8 = rng.choice(e8, (N_BOOT, len(e8))).mean(1)
    b16 = rng.choice(e16, (N_BOOT, len(e16))).mean(1)
    bg = b16 - b8
    se_pair_analytic = np.sqrt(e8.var(ddof=1) / len(e8) + e16.var(ddof=1) / len(e16))

    print()
    print(f"  {'statistic':12s} {'value':>10s} {'SE chr':>9s} {'SE chr(w)':>10s}")
    for i, nm in enumerate(NAMES):
        print(f"  {nm:12s} {full[i]:+10.5f} {se_plain[i]:9.5f} {se_w[i]:10.5f}")
    print()
    print("  field-pair resampling (sees family structure and real between-field differences):")
    print(f"    excess t=8  mean {e8.mean():+.5f}  sd over pairs {e8.std(ddof=1):.5f}  SE {e8.std(ddof=1)/np.sqrt(len(e8)):.5f}")
    print(f"    excess t=16 mean {e16.mean():+.5f}  sd over pairs {e16.std(ddof=1):.5f}  SE {e16.std(ddof=1)/np.sqrt(len(e16)):.5f}")
    print(f"    growth      {full[2]:+.5f}  SE analytic {se_pair_analytic:.5f}  "
          f"bootstrap 95% [{np.percentile(bg, 2.5):+.5f}, {np.percentile(bg, 97.5):+.5f}]  "
          f"P(growth>0) {np.mean(bg > 0):.3f}")
    print(f"    z vs 0 (pair SE): {full[2] / se_pair_analytic:+.2f}   z vs 0 (chr SE): {full[2] / se_plain[2]:+.2f}")

    # pairs excluding the small/renamed/failed ones, as a robustness line
    flagged = {k for k in target if "Arlington-2015" in k[0] or "H41" in k[0]}
    keep8 = [target[k][1] - target[k][2] for k in target if target[k][0] == 8 and k not in flagged]
    keep16 = [target[k][1] - target[k][2] for k in target if target[k][0] == 16 and k not in flagged]
    print(f"\n  robustness -- dropping Arlington-2015 (n=4) and H41 (failed 2023 samples) pairs: "
          f"growth {np.mean(keep16) - np.mean(keep8):+.5f} ({len(keep8)} + {len(keep16)} pairs)")


if __name__ == "__main__":
    main()
