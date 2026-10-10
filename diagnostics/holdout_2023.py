"""The 2023 hold-out (CLAUDE.md 1.1) on the persistent-deme trials already on disk. No simulation.

Pools batch 3 (2,500 trials, persistent demes) and batch 6's OFF arm (refound_k == -1, 628 trials;
the toggle's OFF arm is verified identical to the pre-toggle model, CLAUDE.md 2). Both ran the
current prior and Q=100 constants. For every trial it recomputes pi_loss and fst_loss PER YEAR from
the trial's detail files, with the production helpers and keep mask, and refuses to go on unless the
mean over years reproduces the CSV's pi_loss / fst_loss.

FIT: 2015 + 2019 only. L_j = mean over those two years of the per-year loss (production's own
year-averaging, minus 2023). sigma_j = 1.4826*MAD over the pool. Weights (primary) = each loss's
demographic unique R2 (rank regression on pop, total_migration, m, numClusters -- collect_batch's
rule, CLAUDE.md 7), measured on the FIT-year losses so 2023 never touches the fit.
D = sqrt(sum_j w_j (L_j/sigma_j)^2). Accept the 5% of trials with smallest D.

PREDICT 2023: per accepted trial, simulated 2023 summaries over the kept subpops --
  fst_med    median off-diagonal F_st            (fst_loss chases the median pair, 7.9.13)
  fst_mean   mean off-diagonal F_st              (pulled up by isolates; secondary)
  fst_med_noiso  median with H41-2023 dropped    (failed samples, 7.2.2)
  logpi      mean log pi
  fst_mean_noiso, logpi_noiso   the same with H41-2023 dropped -- ADDED AFTER the first run
             (2026-10-09, TODO 6b), secondary only; they do not change the pre-registered verdict
The predictive interval is the 2.5-97.5% range across accepted trials.

PRE-REGISTERED CRITERION (written before the first run, 2026-10-09):
  PASS on a statistic iff the OBSERVED 2023 value lies inside the accepted trials' 95% interval.
  The primary statistics are fst_med and logpi; each is reported separately -- which half carried
  the test is part of the result. If either primary falls outside, the report says the reframe
  FAILED on that statistic. Power check: the same interval over ALL trials (prior predictive); if
  it covers the observed value about as well, the hold-out has no power and is reported as such,
  not as a pass. The 2023 per-year losses of accepted vs all trials are the second power check.
Sensitivity (reported, not the verdict): equal weights; 20% acceptance; batch 3 vs batch 6-OFF.

No pop posterior is reported (CLAUDE.md 7): the sensitivity-surface section shows where accepted
trials sit in (pop, total_migration, numClusters), relative to the prior, and nothing more.

Usage (from diagnostics/):  python holdout_2023.py [--out ../out/holdout]
"""
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent / "Python_Code"))

import argparse
import csv
import json
import os
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
SOURCES = [  # (label, dir with abc_results_<job>.csv + detailed_sim_results_<job>/, keep-row filter)
    ("batch3", ROOT / "out" / "batch3_raw", lambda r: True),
    ("batch6off", ROOT / "out" / "batch6", lambda r: int(float(r["refound_k"])) == -1),
]
YEARS = ["2015", "2019", "2023"]
FIT_YEARS = ["2015", "2019"]
FITTED = ["pi_loss", "fst_loss"]
DEMOG = ["pop", "total_migration", "m", "numClusters"]
ACCEPT = 0.05
ISOLATE_2023 = "H41-2023"
SUMMARIES = ["fst_med", "fst_mean", "fst_med_noiso", "logpi", "fst_mean_noiso", "logpi_noiso"]
PRIMARY = ["fst_med", "logpi"]


def offd(M, keep):
    M = M[np.ix_(keep, keep)]
    return M[np.triu_indices(M.shape[0], 1)]


def summaries(fst, pi, keep, iso_idx):
    v = offd(fst, keep)
    keep_ni = keep.copy()
    keep_ni[iso_idx] = False
    return {"fst_med": float(np.median(v)), "fst_mean": float(v.mean()),
            "fst_med_noiso": float(np.median(offd(fst, keep_ni))),
            "logpi": float(np.mean(np.log(pi[keep]))),
            "fst_mean_noiso": float(offd(fst, keep_ni).mean()),
            "logpi_noiso": float(np.mean(np.log(pi[keep_ni])))}


def _rank(x):
    order = np.argsort(x, kind="mergesort")
    r = np.empty(len(x))
    r[order] = np.arange(len(x))
    # average ties (as collect_batch._rank)
    xs = x[order]
    i = 0
    while i < len(xs):
        j = i
        while j + 1 < len(xs) and xs[j + 1] == xs[i]:
            j += 1
        if j > i:
            r[order[i:j + 1]] = 0.5 * (i + j)
        i = j + 1
    return r


def _z(v):
    r = _rank(np.asarray(v, dtype=float))
    return (r - r.mean()) / r.std()


def demographic_r2(y, P):
    """Sum of positive unique rank-R2 over DEMOG -- collect_batch.report_decomposition's rule."""
    n = len(y)
    X = np.column_stack([np.ones(n)] + [_z(P[p]) for p in DEMOG])
    yz = _z(y)
    sst = float(np.sum((yz - yz.mean()) ** 2))
    b, *_ = np.linalg.lstsq(X, yz, rcond=None)
    r2 = 1 - float(np.sum((yz - X @ b) ** 2)) / sst
    tot = 0.0
    for k in range(len(DEMOG)):
        Xr = np.delete(X, k + 1, axis=1)
        b2, *_ = np.linalg.lstsq(Xr, yz, rcond=None)
        tot += max(0.0, r2 - (1 - float(np.sum((yz - Xr @ b2) ** 2)) / sst))
    return tot, r2


def robust_sigma(v):
    v = np.asarray(v, dtype=float)
    return 1.4826 * float(np.median(np.abs(v - np.median(v))))


def load(ABC, obs, masks):
    rows = []
    for label, d, keep_row in SOURCES:
        for f in sorted(d.glob("abc_results_*.csv")):
            tail = f.stem[len("abc_results_"):]
            if not tail.isdigit():
                continue
            jid = int(tail)
            for r in csv.DictReader(open(f, newline="", encoding="utf-8")):
                if r.get("iteration") in (None, "iteration") or not keep_row(r):
                    continue
                it = int(r["iteration"])
                dd = d / f"detailed_sim_results_{jid}" / f"run{it + 1}"
                o = {"source": label, "job_id": jid, "iteration": it}
                o.update({p: float(r[p]) for p in DEMOG})
                o["csv_pi_loss"], o["csv_fst_loss"] = float(r["pi_loss"]), float(r["fst_loss"])
                for y in YEARS:
                    keep = masks[y]
                    kk = np.ix_(keep, keep)
                    fst = ABC._read_matrix(dd / f"fst_{y}.csv")
                    pi = ABC._read_vector(dd / f"diversities_{y}.csv")
                    o[f"pi_{y}"] = ABC._pi_log_mean_abs_diff(obs[f"{y}_diversity"][keep], pi[keep])
                    o[f"fst_{y}"] = ABC._offdiag_mean_abs_diff(obs[f"{y}_fst"][kk], fst[kk])
                    if y == "2023":
                        o.update({f"sim_{k}": v for k, v in
                                  summaries(fst, pi, keep, masks["iso"]).items()})
                rows.append(o)
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(ROOT / "out" / "holdout"))
    args = ap.parse_args()
    out = Path(args.out).resolve()
    out.mkdir(parents=True, exist_ok=True)

    os.chdir(ROOT / "Python_Code")   # ABC reads its inputs by relative path
    import ABCAnalysisNoRedis as ABC

    obs = ABC.getObservedData()
    masks = {y: ABC.get_keep_mask(y) for y in YEARS}
    names23 = ABC._specifier_site_names("2023")
    masks["iso"] = names23.index(ISOLATE_2023)
    obs_s = summaries(obs["2023_fst"], obs["2023_diversity"], masks["2023"], masks["iso"])

    rows = load(ABC, obs, masks)
    n = len(rows)
    T = {k: np.array([r[k] for r in rows]) for k in rows[0] if k != "source"}
    src = np.array([r["source"] for r in rows])
    print(f"Loaded {n} trials: " + ", ".join(f"{s} {int((src == s).sum())}" for s in np.unique(src)))

    # --- self-check: per-year losses must reproduce the CSV (production's mean over years)
    for s in ["pi", "fst"]:
        recon = np.nanmean(np.column_stack([T[f"{s}_{y}"] for y in YEARS]), axis=1)
        bad = ~np.isclose(recon, T[f"csv_{s}_loss"], rtol=1e-9, atol=1e-15)
        if bad.any():
            raise SystemExit(f"{s}_loss: {int(bad.sum())} trials do not reproduce the CSV "
                             f"(max abs diff {np.max(np.abs(recon - T[f'csv_{s}_loss'])):.3g})")
    print("Self-check: per-year pi/fst losses reproduce the CSV losses for every trial.")

    L = {s: np.mean(np.column_stack([T[f"{s.split('_')[0]}_{y}"] for y in FIT_YEARS]), axis=1)
         for s in FITTED}
    sig = {s: robust_sigma(L[s]) for s in FITTED}
    dem = {s: demographic_r2(L[s], T) for s in FITTED}
    tot = sum(d for d, _ in dem.values())
    W_r2 = {s: dem[s][0] / tot for s in FITTED}
    W_eq = {s: 1 / len(FITTED) for s in FITTED}

    print("\nFIT (2015+2019 only)")
    for s in FITTED:
        print(f"  {s:9s} sigma {sig[s]:.6f}  R2_tot {dem[s][1]:.3f}  demographic {dem[s][0]:.3f}"
              f"  -> weight {W_r2[s]:.3f}")
    print(f"  Spearman(pi_fit, fst_fit) = {np.corrcoef(_rank(L['pi_loss']), _rank(L['fst_loss']))[0,1]:+.3f}")

    def D(W):
        return np.sqrt(sum(W[s] * (L[s] / sig[s]) ** 2 for s in FITTED))

    def run(W, frac, sel=None, tag=""):
        sel = np.ones(n, bool) if sel is None else sel
        d = np.where(sel, D(W), np.inf)
        k = max(1, int(round(frac * sel.sum())))
        acc = np.argsort(d)[:k]
        res = {"k": k, "stats": {}}
        for s in SUMMARIES:
            a, p = T[f"sim_{s}"][acc], T[f"sim_{s}"][sel]
            lo, hi = np.percentile(a, [2.5, 97.5])
            plo, phi = np.percentile(p, [2.5, 97.5])
            o = obs_s[s]
            res["stats"][s] = {"obs": o, "acc_lo": lo, "acc_med": float(np.median(a)), "acc_hi": hi,
                               "acc_pct": float((a < o).mean() * 100),
                               "prior_lo": plo, "prior_hi": phi,
                               "prior_pct": float((p < o).mean() * 100),
                               "inside": bool(lo <= o <= hi),
                               "width_ratio": float((hi - lo) / (phi - plo))}
        for s in ["pi_2023", "fst_2023"]:
            res[s] = {"acc_med": float(np.median(T[s][acc])), "prior_med": float(np.median(T[s][sel])),
                      "acc_rank_pct": float(np.mean([(T[s][sel] < v).mean() for v in T[s][acc]]) * 100)}
        res["acc_idx"] = acc
        return res

    def show(res, title):
        print(f"\n{title}  (accepted {res['k']})")
        print(f"  {'summary':14s} {'observed':>10s} {'acc 2.5%':>10s} {'acc med':>10s} {'acc 97.5%':>10s}"
              f" {'obs pct':>8s} {'prior 2.5%':>10s} {'prior 97.5%':>11s} {'pr pct':>7s} {'width':>6s}  verdict")
        for s in SUMMARIES:
            x = res["stats"][s]
            v = "inside" if x["inside"] else "OUTSIDE"
            star = "*" if s in PRIMARY else " "
            print(f" {star}{s:14s} {x['obs']:10.5f} {x['acc_lo']:10.5f} {x['acc_med']:10.5f} "
                  f"{x['acc_hi']:10.5f} {x['acc_pct']:7.1f}% {x['prior_lo']:10.5f} {x['prior_hi']:11.5f}"
                  f" {x['prior_pct']:6.1f}% {x['width_ratio']:6.2f}  {v}")
        for s in ["pi_2023", "fst_2023"]:
            x = res[s]
            print(f"  2023 {s[:-5]}_loss: accepted median {x['acc_med']:.5f} vs all-trial median "
                  f"{x['prior_med']:.5f}; accepted trials sit at the {x['acc_rank_pct']:.0f}th pct "
                  f"of all trials on average (50 = no information)")

    primary = run(W_r2, ACCEPT)
    show(primary, f"PRIMARY: demographic-R2 weights, top {ACCEPT:.0%}")
    sens = {"equal_w_5pct": run(W_eq, ACCEPT), "r2_w_20pct": run(W_r2, 0.20),
            "batch3_only": run(W_r2, ACCEPT, src == "batch3"),
            "batch6off_only": run(W_r2, ACCEPT, src == "batch6off")}
    for k, r in sens.items():
        show(r, f"SENSITIVITY: {k}")

    # --- sensitivity surface: where the accepted trials sit, relative to the prior
    acc = primary["acc_idx"]
    print("\nWHERE THE ACCEPTED TRIALS SIT (primary). Enrichment = accepted share / prior share.")
    q = {"pop": np.quantile(T["pop"], [0, .25, .5, .75, 1]),
         "total_migration": np.quantile(T["total_migration"], [0, .25, .5, .75, 1]),
         "m": np.quantile(T["m"], [0, .25, .5, .75, 1])}
    surf = []
    for p in ["pop", "total_migration", "m"]:
        e = q[p]
        cells = []
        for b in range(4):
            inb = (T[p] >= e[b]) & (T[p] <= e[b + 1])
            enr = np.isin(np.where(inb)[0], acc).sum() / len(acc) / inb.mean()
            cells.append(f"Q{b+1} [{e[b]:.3g},{e[b+1]:.3g}] x{enr:.2f}")
            surf.append({"param": p, "bin": f"Q{b+1}", "lo": e[b], "hi": e[b + 1], "enrichment": enr})
        print(f"  {p:16s} " + "  ".join(cells))
    cells = []
    for c in [1, 2, 3]:
        inb = T["numClusters"] == c
        enr = np.isin(np.where(inb)[0], acc).sum() / len(acc) / inb.mean()
        cells.append(f"{c*33} demes x{enr:.2f}")
        surf.append({"param": "numClusters", "bin": str(c * 33), "lo": c, "hi": c, "enrichment": enr})
    print(f"  {'numClusters':16s} " + "  ".join(cells))
    print("  pop x total_migration (rows pop quartile, cols tm quartile; accepted counts):")
    for i in range(4):
        pi_ = (T["pop"] >= q["pop"][i]) & (T["pop"] <= q["pop"][i + 1])
        cnt = []
        for j in range(4):
            tj = (T["total_migration"] >= q["total_migration"][j]) & (T["total_migration"] <= q["total_migration"][j + 1])
            cnt.append(int(np.isin(np.where(pi_ & tj)[0], acc).sum()))
        print(f"    pop Q{i+1}: " + " ".join(f"{c:4d}" for c in cnt))
    rho = np.corrcoef(_rank(T["pop"][acc]), _rank(T["total_migration"][acc]))[0, 1]
    print(f"  Spearman(pop, total_migration) among accepted: {rho:+.3f} (a ridge shows as negative)")

    # --- write
    with open(out / "holdout_trials.csv", "w", newline="") as fh:
        keys = ["source", "job_id", "iteration"] + DEMOG + [f"{s}_{y}" for s in ["pi", "fst"] for y in YEARS] \
               + [f"sim_{s}" for s in SUMMARIES]
        w = csv.writer(fh)
        w.writerow(keys + ["D_primary", "accepted_primary"])
        dp = D(W_r2)
        accset = set(acc.tolist())
        for i, r in enumerate(rows):
            w.writerow([r[k] for k in keys] + [dp[i], int(i in accset)])
    with open(out / "surface.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(surf[0]))
        w.writeheader()
        w.writerows(surf)
    strip = lambda r: {k: v for k, v in r.items() if k != "acc_idx"}
    with open(out / "holdout_summary.json", "w") as fh:
        json.dump({"n_trials": n, "sources": {s: int((src == s).sum()) for s in np.unique(src)},
                   "fit_years": FIT_YEARS, "sigma": sig, "weights_r2": W_r2,
                   "demographic_r2": {s: dem[s][0] for s in FITTED}, "observed_2023": obs_s,
                   "primary": strip(primary), "sensitivity": {k: strip(v) for k, v in sens.items()}},
                  fh, indent=2)
    print(f"\nWrote {out / 'holdout_trials.csv'}, surface.csv, holdout_summary.json")


if __name__ == "__main__":
    main()
