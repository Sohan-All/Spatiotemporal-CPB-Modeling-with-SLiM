"""Per-trial SIGNED features for a batch, from its raw detail files (CLAUDE.md 7.9.13 recipe).

The losses in abc_results_*.csv are absolute differences, so they cannot say which side of the data
a trial sits on. This reads each trial's detailed_sim_results_<job>/run<iteration+1>/ and writes one
row per trial with:

  fc8, fc16            pooled simulated F_c per gap (unweighted mean over field pairs, as _fc_loss)
  ex8, ex16            the same minus the pooled pedestal
  growth               ex16 - ex8
  exsd8, exsd16        sd of the per-pair excess within the gap -- the simulated between-pair spread,
                       to set against the observed 0.0060 / 0.0061 (fc_growth_jackknife.py)
  fstr_<year>          masked mean off-diagonal sim F_st / observed mean
  fstmed_<year>        masked MEDIAN off-diagonal sim F_st (fst_loss chases the median pair)
  logpi_<year>         mean log(sim pi / obs pi) over kept subpops
plus the parameters and fitted losses. Runs in about a minute for 2,500 trials, no simulation.

Usage (from diagnostics/):
    python batch_signed_features.py --batch-dir ../out/batch6 --out ../out/batch6/signed_features.csv
"""
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent / "Python_Code"))

import argparse
import csv
import os
from pathlib import Path

import numpy as np

YEARS = ["2015", "2019", "2023"]
PARAMS = ["m", "total_migration", "pop", "numClusters", "refound_k"]
LOSSES = ["pi_loss", "fst_loss", "fc_loss"]


def fc_by_gap(path):
    g = {}
    for r in csv.DictReader(open(path, newline="", encoding="utf-8")):
        g.setdefault(int(r["t"]), []).append((float(r["fc"]), float(r["pedestal"])))
    return {t: (np.mean([a for a, _ in v]), np.mean([b for _, b in v]),
                float(np.std([a - b for a, b in v], ddof=1))) for t, v in g.items()}


def offd(M, keep):
    M = M[np.ix_(keep, keep)]
    return M[np.triu_indices(M.shape[0], 1)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--batch-dir", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    B = Path(args.batch_dir).resolve()
    out_path = Path(args.out).resolve()

    # ABCAnalysisNoRedis reads its inputs by paths relative to Python_Code/
    os.chdir(Path(__file__).resolve().parent.parent / "Python_Code")
    import ABCAnalysisNoRedis as ABC

    obs = ABC.getObservedData()
    masks = {y: ABC.get_keep_mask(y) for y in YEARS}
    obs_fst = {y: offd(obs[f"{y}_fst"], masks[y]).mean() for y in YEARS}
    obs_pi = {y: obs[f"{y}_diversity"][masks[y]] for y in YEARS}

    rows = []
    for f in sorted(B.glob("abc_results_*.csv")):
        jid = int(f.stem.split("_")[-1])
        for r in csv.DictReader(open(f, newline="", encoding="utf-8")):
            it = int(r["iteration"])
            d = B / f"detailed_sim_results_{jid}" / f"run{it + 1}"
            fc = fc_by_gap(d / "temporal_fc.csv")
            o = {"job_id": jid, "iteration": it}
            o.update({k: float(r[k]) for k in PARAMS + LOSSES})
            o["fc8"], o["fc16"] = fc[8][0], fc[16][0]
            o["ex8"], o["ex16"] = fc[8][0] - fc[8][1], fc[16][0] - fc[16][1]
            o["growth"] = o["ex16"] - o["ex8"]
            o["exsd8"], o["exsd16"] = fc[8][2], fc[16][2]
            for y in YEARS:
                v = offd(ABC._read_matrix(d / f"fst_{y}.csv"), masks[y])
                o[f"fstr_{y}"] = v.mean() / obs_fst[y]
                o[f"fstmed_{y}"] = float(np.median(v))
                pv = ABC._read_vector(d / f"diversities_{y}.csv")
                o[f"logpi_{y}"] = float(np.mean(np.log(pv[masks[y]] / obs_pi[y])))
            rows.append(o)
    with open(out_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"{len(rows)} trials -> {out_path}")


if __name__ == "__main__":
    main()
