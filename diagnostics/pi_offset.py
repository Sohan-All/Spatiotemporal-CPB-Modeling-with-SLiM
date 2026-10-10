"""TODO 6b: is the ~1% between-year pi offset (CLAUDE.md 7.9.15) the het deficit acting through n?

pixy's pi averages over all C(2n,2) haplotype pairs, n of which are the two haplotypes of ONE
individual. A per-individual het deficit e (CLAUDE.md 7.2.2 -- miscall OR real inbreeding, the
algebra is the same) removes a fraction e of the within-individual differences and leaves the
between-individual pairs unbiased, so

    E[pi_sample] = pi * (1 - e_bar / (2n - 1)),     e_bar = site mean of the per-sample e.

The simulated pi (whole demes, no deficit) carries no such factor. So the predicted log offset per
site is log(1 - e_bar/(2n-1)): ~-3.5% at n=4, ~-1% at n=11. d_xy has no within-individual pairs
and is the control -- it should show no such offset.

Reference: the accepted trials of the hold-out (out/holdout/holdout_trials.csv, primary). For each
kept site, resid = log pi_obs - median over accepted trials of log pi_sim. Reports year means of
resid before/after removing the predicted factor, and the site-level regression of resid on
x = e_bar/(2n-1) (expected slope ~ -1, since log(1-x) ~ -x).

Usage (from diagnostics/):  python pi_offset.py [--holdout ../out/holdout]
                            (needs holdout_2023.py's output; no sim)
Per-site n and e_bar leave out sample_exclusions.EXCLUDED_SAMPLES, matching the target (7.9.16).
"""
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent / "Python_Code"))

import argparse
import csv
import os
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
YEARS = ["2015", "2019", "2023"]
ISOLATES = {"Mortensen9-2015", "H41-2023"}   # family; failed samples (CLAUDE.md 7.1-7.5, 7.2.2)
SRC = {"batch3": ROOT / "out" / "batch3_raw", "batch6off": ROOT / "out" / "batch6"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--holdout", default=str(ROOT / "out" / "holdout"),
                    help="holdout_2023.py output directory whose accepted trials are the reference")
    args = ap.parse_args()
    holdout = Path(args.holdout).resolve()
    os.chdir(ROOT / "Python_Code")
    import ABCAnalysisNoRedis as ABC
    from sample_exclusions import EXCLUDED_SAMPLES, exclusion_id
    print(f"sample_exclusions: [{exclusion_id()}]   reference: {holdout}")

    obs = ABC.getObservedData()
    acc = [r for r in csv.DictReader(open(holdout / "holdout_trials.csv"))
           if r["accepted_primary"] == "1"]
    dirs = [SRC[r["source"]] / f"detailed_sim_results_{r['job_id']}" / f"run{int(r['iteration']) + 1}"
            for r in acc]

    e_by = {}
    for r in csv.DictReader(open(ROOT / "data" / "empiricalStats" / "miscall_rates.csv")):
        e_by[r["sample_id"]] = float(r["e"])

    allx, allr, ally = [], [], []
    print(f"Reference: {len(acc)} accepted hold-out trials.\n")
    print(f"{'year':5s} {'sites':>5s} {'mean n':>6s} {'pi resid':>9s} {'predicted':>9s} {'corrected':>9s}"
          f" {'dxy resid':>9s}")
    for y in YEARS:
        keep = ABC.get_keep_mask(y)
        names = ABC._specifier_site_names(y)
        members = {}
        for line in open(ROOT / "data" / "Genetic_Data" / f"popFile{y}", encoding="utf-8"):
            p = line.split()
            if len(p) >= 2 and p[0] not in EXCLUDED_SAMPLES:   # same samples as the target (7.9.16)
                members.setdefault(p[1], []).append(p[0])
        miss = [s for v in members.values() for s in v if s not in e_by]
        if miss:
            raise SystemExit(f"{y}: no miscall rate for {miss}")

        simpi = np.array([np.log(ABC._read_vector(d / f"diversities_{y}.csv")) for d in dirs])
        simdx = np.array([ABC._read_matrix(d / f"divergences_{y}.csv") for d in dirs])
        n = np.array([len(members[s]) for s in names])
        eb = np.array([np.mean([e_by[s] for s in members[nm]]) for nm in names])
        x = eb / (2 * n - 1)
        pred = np.log(1 - x)
        resid = np.log(obs[f"{y}_diversity"]) - np.median(simpi, axis=0)

        kk = np.ix_(keep, keep)
        iu = np.triu_indices(int(keep.sum()), 1)
        dres = (np.log(obs[f"{y}_divergence"][kk][iu]) - np.log(np.median(simdx, axis=0)[kk][iu]))

        k = keep
        print(f"{y:5s} {int(k.sum()):5d} {n[k].mean():6.1f} {resid[k].mean():+9.4f} {pred[k].mean():+9.4f}"
              f" {(resid[k] - pred[k]).mean():+9.4f} {dres.mean():+9.4f}")
        allx += list(x[k]); allr += list(resid[k]); ally += [y] * int(k.sum())
        kn = [i for i in np.where(k)[0] if names[i] not in ISOLATES]
        low = sorted(np.where(k)[0], key=lambda i: resid[i])[:3]
        print(f"      without isolates: pi resid {resid[kn].mean():+.4f}   lowest sites: " +
              ", ".join(f"{names[i]} {resid[i]:+.3f} (n={n[i]}, e={eb[i]:.2f})" for i in low))

    allx, allr, ally = np.array(allx), np.array(allr), np.array(ally)
    spread = lambda v: max(v) - min(v)
    yr = {y: allr[ally == y].mean() for y in YEARS}
    yc = {y: (allr - np.log(1 - allx))[ally == y].mean() for y in YEARS}
    print(f"\nBetween-year spread of mean pi resid: {spread(list(yr.values())):.4f} raw -> "
          f"{spread(list(yc.values())):.4f} after the predicted factor")

    # Site-level: resid on x with year fixed effects; slope expected ~ -1
    X = np.column_stack([allx] + [(ally == y).astype(float) for y in YEARS])
    b, res, *_ = np.linalg.lstsq(X, allr, rcond=None)
    dof = len(allr) - X.shape[1]
    s2 = float(np.sum((allr - X @ b) ** 2)) / dof
    se = float(np.sqrt(s2 * np.linalg.inv(X.T @ X)[0, 0]))
    X0 = np.column_stack([allx, np.ones(len(allx))])
    b0, *_ = np.linalg.lstsq(X0, allr, rcond=None)
    print(f"Site regression of pi resid on e_bar/(2n-1), {len(allr)} sites:")
    print(f"  with year effects: slope {b[0]:+.2f} (SE {se:.2f}); expected about -1")
    print(f"  pooled, no year effects: slope {b0[0]:+.2f}")
    r = np.corrcoef(allx, allr)[0, 1]
    print(f"  correlation {r:+.3f}")


if __name__ == "__main__":
    main()
