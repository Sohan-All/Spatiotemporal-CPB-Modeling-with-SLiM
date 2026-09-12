"""The fc_loss noise floor, measured through the PRODUCTION path (CLAUDE.md 7.9.9E, TODO 8.5).

7.9.8's floor was on the pooled STATISTIC, from diagnostics/temporal_fc.py, on clean genotypes.
fc_loss is a different object: |pooled F_c,sim - pooled F_c,obs| per generation gap, averaged over
the gaps, with the simulated side subsampled per field-year, carrying imitated het miscalls, and with
fc_common.EXCLUDED_SAMPLES applied (7.2.2H). None of that was in the old floor. So this re-measures
it the only way that cannot drift from production: every replicate is a call to
ABCAnalysisNoRedis.model() -> calculate_losses(), exactly what a CHTC trial runs.

What is re-rolled per replicate: SLiM's forward mating (Main passes no -s), recapitation, the
mutation overlay, the per-field-year subsample and the miscall-rate assignment -- all five dice a
real trial rolls. Clustering and the kernel are deterministic, as in production.

Every loss is recorded, not only fc_loss, so this also gives the Hudson-scale fst_loss floor at the
Q=100 constants that 7.3 deferred. Per replicate it records the pooled simulated F_c per gap too:
fc_loss folds an absolute value, and the pooled level is what moves with POPMULT.

Pooling (--summarize) refuses records that disagree on any spec hash, scale constant or parameter
(7.9.8D). All replicates of one invocation run in ONE process, so the code cannot change under them.

WRITES data/ IN PLACE, like production (CLAUDE.md 10). The files the pipeline overwrites are backed
up to out/fc_loss_floor/_backup/ first and restored on exit, including on error. If the process is
KILLED the restore does not run -- copy them back from there by hand.

Usage (from diagnostics/):
    python fc_loss_floor.py --popmult 2000 --reps 5
    python fc_loss_floor.py --popmult 5000 --reps 3
    python fc_loss_floor.py --summarize
"""
import os
os.environ["COMPUTE_FC"] = "1"   # both modules read it AT IMPORT, so it must precede them

import argparse
import itertools
import json
import shutil
import subprocess
import sys
import threading
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "Python_Code"))
import ABCAnalysisNoRedis as ABC   # noqa: E402
import fc_common as fcc            # noqa: E402
import ld_common as ldc            # noqa: E402
import scale_constants as sc       # noqa: E402

ROOT = HERE.parent
OUTDIR = ROOT / "out" / "fc_loss_floor"
BACKUP = OUTDIR / "_backup"
OVERWRITTEN = [ROOT / "data" / "cluster_data.csv", ROOT / "data" / "cluster_distances.csv",
               ROOT / "data" / "migration_rates.csv", ROOT / "out" / "simTreeSeq.trees"]
OUTPUT_DATA = ROOT / "data" / "Output_Data"

# Records pool only if ALL of these agree.
POOL_KEYS = ("fc_spec", "ld_spec", "fc_empirical_spec", "m", "total_migration", "numClusters",
             "mutation_rate", "recombination_rate", "ancestral_Ne", "Q")


def _backup():
    BACKUP.mkdir(parents=True, exist_ok=True)
    if any(BACKUP.iterdir()):
        raise SystemExit(f"{BACKUP} is not empty -- a previous run was killed before restoring. "
                         f"Restore those files by hand (or confirm they are stale) and delete it.")
    saved = []
    for p in OVERWRITTEN:
        if p.exists():
            shutil.copy2(p, BACKUP / p.name); saved.append(p)
    (BACKUP / "Output_Data").mkdir()
    for p in OUTPUT_DATA.glob("*.csv"):
        shutil.copy2(p, BACKUP / "Output_Data" / p.name); saved.append(p)
    return saved


def _restore(saved):
    had = {p.name for p in (BACKUP / "Output_Data").glob("*.csv")}
    for p in OVERWRITTEN:
        if (BACKUP / p.name).exists():
            shutil.copy2(BACKUP / p.name, p)
    for name in had:
        shutil.copy2(BACKUP / "Output_Data" / name, OUTPUT_DATA / name)
    for p in OUTPUT_DATA.glob("*.csv"):     # files the run CREATED (temporal_fc.csv) go away
        if p.name not in had:
            p.unlink()
    shutil.rmtree(BACKUP)
    print(f"restored {len(saved)} files from backup", flush=True)


class PeakRSS:
    """High-water RSS of this process and of its children (SLiM), sampled every 0.5 s."""
    def __init__(self):
        import psutil
        self.proc, self.me, self.kids, self._stop = psutil.Process(), 0, 0, threading.Event()

    def _run(self):
        import psutil
        while not self._stop.is_set():
            try:
                self.me = max(self.me, self.proc.memory_info().rss)
                k = sum(c.memory_info().rss for c in self.proc.children(recursive=True))
                self.kids = max(self.kids, k)
            except psutil.Error:
                pass
            self._stop.wait(0.5)

    def __enter__(self):
        self._t = threading.Thread(target=self._run, daemon=True); self._t.start(); return self

    def __exit__(self, *exc):
        self._stop.set(); self._t.join()


def _git_head():
    try:
        h = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=ROOT, capture_output=True,
                           text=True, check=True).stdout.strip()
        dirty = subprocess.run(["git", "status", "--porcelain", "--", "Python_Code", "SLiM_Code"],
                               cwd=ROOT, capture_output=True, text=True).stdout.strip()
        return h + ("+dirty" if dirty else "")
    except Exception:
        return None


def _pooled_by_gap(table):
    by = {}
    for t, fc in table.values():
        if np.isfinite(fc):
            by.setdefault(t, []).append(fc)
    return {str(t): float(np.mean(v)) for t, v in sorted(by.items())}


def main(a):
    if fcc.spec_hash() != ABC.FC_EMPIRICAL_SPEC:
        raise SystemExit(f"fc spec {fcc.spec_hash()} != FC_EMPIRICAL_SPEC {ABC.FC_EMPIRICAL_SPEC}")
    obs = ABC.getObservedData()
    obs_gap = _pooled_by_gap(obs["temporal_fc"])
    params = {"m": a.m, "total_migration": a.total_migration, "pop": a.popmult,
              "numClusters": a.num_clusters, "mutation_rate": sc.MUTATION_RATE,
              "recombination_rate": sc.RECOMBINATION_RATE}
    common = {"kind": "replicate", "popmult": a.popmult, "m": a.m,
              "total_migration": a.total_migration, "numClusters": a.num_clusters,
              "mutation_rate": sc.MUTATION_RATE, "recombination_rate": sc.RECOMBINATION_RATE,
              "ancestral_Ne": sc.ANCESTRAL_NE, "Q": sc.Q, "fc_spec": fcc.spec_hash(),
              "ld_spec": ldc.spec_hash(), "fc_empirical_spec": ABC.FC_EMPIRICAL_SPEC,
              "git": _git_head(), "fc_pooled_obs": obs_gap}
    print(f"POPMULT={a.popmult} x{a.reps}  m={a.m} tm={a.total_migration} "
          f"demes={a.num_clusters * 33}  fc spec {common['fc_spec']}  git {common['git']}",
          flush=True)

    saved = _backup()
    try:
        for r in range(a.reps):
            stamp = time.strftime("%Y%m%d-%H%M%S")
            print(f"\n=== rep {r + 1}/{a.reps}  [{stamp}] ===", flush=True)
            t0 = time.perf_counter()
            with PeakRSS() as mem:
                sim = ABC.model(dict(params))
            dt = time.perf_counter() - t0
            losses = ABC.calculate_losses(obs, sim)
            sim_gap = _pooled_by_gap(sim["temporal_fc"])

            repdir = OUTDIR / f"p{int(a.popmult)}_{stamp}"
            repdir.mkdir(parents=True, exist_ok=True)
            for p in OUTPUT_DATA.glob("*.csv"):
                shutil.copy2(p, repdir / p.name)

            rec = {**common, "stamp": stamp, "seconds": round(dt, 1),
                   "peak_python_mb": round(mem.me / 2**20), "peak_children_mb": round(mem.kids / 2**20),
                   "losses": losses, "fc_pooled_sim": sim_gap}
            with open(a.out, "a", encoding="utf-8") as f:
                f.write(json.dumps(rec) + "\n")
            print(f"  {dt / 60:.1f} min, peak python {rec['peak_python_mb']} MB / "
                  f"children {rec['peak_children_mb']} MB", flush=True)
            print("  " + ABC._format_losses(losses), flush=True)
            print("  pooled sim F_c " + "  ".join(f"t{t}={v:.5f}" for t, v in sim_gap.items())
                  + "   obs " + "  ".join(f"t{t}={v:.5f}" for t, v in obs_gap.items()), flush=True)
    finally:
        _restore(saved)


def _stats(v):
    v = np.asarray(v, float)
    pair = [abs(x - y) for x, y in itertools.combinations(v, 2)]
    sd = float(v.std(ddof=1)) if len(v) > 1 else float("nan")
    return float(v.mean()), sd, float(np.mean(pair)) if pair else float("nan")


def summarize(a):
    recs = [json.loads(l) for l in open(a.out, encoding="utf-8") if l.strip()]
    recs = [r for r in recs if r.get("kind") == "replicate"]
    if not recs:
        raise SystemExit(f"no records in {a.out}")
    for k in POOL_KEYS:
        vals = {json.dumps(r[k]) for r in recs}
        if len(vals) > 1:
            raise SystemExit(f"refusing to pool: records disagree on {k} -> {sorted(vals)}")
    print(f"{len(recs)} records; fc spec {recs[0]['fc_spec']}, m={recs[0]['m']}, "
          f"tm={recs[0]['total_migration']}, demes={recs[0]['numClusters'] * 33}, "
          f"git {sorted({r['git'] for r in recs}, key=str)}")
    obs = recs[0]["fc_pooled_obs"]
    gaps = sorted(obs, key=int)
    by_p = {}
    for r in recs:
        by_p.setdefault(r["popmult"], []).append(r)

    table = {}
    for p in sorted(by_p):
        rs = by_p[p]
        print(f"\nPOPMULT={p}  n={len(rs)}  mean {np.mean([r['seconds'] for r in rs]) / 60:.1f} min/rep, "
              f"peak python {max(r['peak_python_mb'] for r in rs)} MB, "
              f"SLiM {max(r['peak_children_mb'] for r in rs)} MB")
        print(f"  {'quantity':<16} {'mean':>10} {'sd':>9} {'mean|diff|':>10} {'CV%':>6}")
        rows = [(k, [r["losses"][k] for r in rs]) for k in ABC.LOSS_NAMES]
        rows += [(f"F_c sim t{t}", [r["fc_pooled_sim"][t] for r in rs]) for t in gaps]
        table[p] = {}
        for name, v in rows:
            m, sd, md = _stats(v)
            table[p][name] = (m, sd)
            print(f"  {name:<16} {m:>10.5f} {sd:>9.5f} {md:>10.5f} {100 * sd / m if m else 0:>6.1f}")
        print("  obs pooled F_c  " + "  ".join(f"t{t}={obs[t]:.5f}" for t in gaps))

    ps = sorted(table)
    for lo, hi in itertools.combinations(ps, 2):
        print(f"\n--- signal POPMULT {lo} -> {hi}, in units of the replicate sd ---")
        for name in ["fc_loss"] + [f"F_c sim t{t}" for t in gaps] + ["fst_loss", "pi_loss"]:
            (m1, s1), (m2, s2) = table[lo][name], table[hi][name]
            sd = float(np.sqrt((s1 ** 2 + s2 ** 2) / 2))
            print(f"  {name:<16} {m1:.5f} -> {m2:.5f}   change {m2 - m1:+.5f}   "
                  f"pooled sd {sd:.5f}   separation {abs(m2 - m1) / sd if sd else float('nan'):.1f} sd"
                  f"   floor/signal {100 * sd / abs(m2 - m1) if m2 != m1 else float('nan'):.0f}%")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--popmult", type=float, default=2000)
    ap.add_argument("--reps", type=int, default=5)
    ap.add_argument("--m", type=float, default=5e-5, help="kernel decay (7.9.6 baseline)")
    ap.add_argument("--total-migration", type=float, default=0.05)
    ap.add_argument("--num-clusters", type=int, default=1, help="RAW draw; demes = 33x this")
    ap.add_argument("--out", default=str(ROOT / "out" / "fc_loss_floor.jsonl"))
    ap.add_argument("--summarize", action="store_true")
    args = ap.parse_args()
    summarize(args) if args.summarize else main(args)
