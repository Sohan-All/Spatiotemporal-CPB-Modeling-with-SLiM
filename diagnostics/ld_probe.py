"""Does the SIMULATED LD decay curve move with POPMULT, in the bins we can afford to fit?

Why this exists (CLAUDE.md 7.5, and the 2026-09-07 empirical run). The empirical side has now
been measured: half-decay of the excess-over-plateau is ~60 bp in all three years, the plateau is
the 1/n_hap floor (corr 0.995/0.984/0.769 against 1/n_hap), and the pooled curve still carries
34% of its dynamic range above 100 bp and 13-14% above 562 bp. Two consequences set this script
up, and they happen to agree on where to cut:

  COST. ~60 bp sits below thin=25's ~98 bp site spacing, so resolving the steep part of the curve
  needs thin=1, measured at ~141 s/deme -- ~3.9 h/trial at 33 demes x 3 years, against a ~1.2 h
  trial. Infeasible.

  > CORRECTED 2026-09-07, mid-run. An earlier draft of this note said the >=562 bp bins were
  > "resolved by the existing STAGES for ~1 min/trial". WRONG: ld_common.STAGES is
  > [(10_000, 1), (1_000_000, 25)], so thin=25 applies only BEYOND 10 kb and every fitted bin
  > from 562 bp to 10 kb is enumerated at full ~3.9 bp density -- which is exactly where the cost
  > lives. The affordable rows of the 7.5c table assume an INPUT-level thin that was not
  > implemented anywhere until `calculate_ld_decay(..., thin=)` was added. Use --thin 25.

  TIME DEPTH. At r = 2.75e-6 the 324-generation forward window controls distances
  d >= 1/(2*324*r) ~ 561 bp; everything shorter reflects coalescence in the FIXED ancestral phase
  (Ne=6700) and therefore carries no POPMULT signal at all (Hayes et al. 2003, CLAUDE.md 11).

So: fit bins with bin_lo >= 562. This script measures whether that region actually moves with
POPMULT. If it does not, LD does not break the N/m confound (7.4.2) and steps 4-6 of 7.5 should
not be built.

WHAT IT DOES NOT SETTLE. LD half-decay measures 4*Ne*r -- structurally the same confound as pi's
4*Ne*mu. It helps only because MIGRATION IS ABSENT FROM IT, not because it is assumption-free:
any N it yields is conditional on r exactly as pi's is on mu, and r itself is disputed (6.8).

One POPMULT per process, appending to out/ld_probe.jsonl, so a failure costs one point.

Usage (run from diagnostics/ -- paths are ../data, ../out):
    # the .trees already on disk is POPMULT=5000; reuse it, no SLiM needed
    python ld_probe.py --popmult 5000 --skip-slim --save-ts ../out/ld_probe_p5000.trees

    # a second point. POPMULT=2000 recapitates in minutes and fits in 16 GB (3.1)
    python ld_probe.py --popmult 2000

    python ld_probe.py --summarize
"""
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent / 'Python_Code'))
import scale_constants as _sc
import recapitate_util as _ru

import argparse
import csv
import gc
import json
import math
import platform
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import msprime
import pyslim
import tskit

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "Python_Code"))
import ld_common as ldc                      # noqa: E402
import ABCAnalysisNoRedis as ABC             # noqa: E402  (real get_keep_mask / readers)
import AnalyzeTreeSeq as ATS                 # noqa: E402  (real calculate_ld_decay)

TIMES = {"2015": 16, "2019": 8, "2023": 0}
ANCESTRAL_NE = _sc.ANCESTRAL_NE
JSONL = Path("../out/ld_probe.jsonl")

# The spec the empirical run was computed under (ldCalcOut.txt, 2026-09-07). A mismatch means the
# two curves are binned differently and any ld_loss is meaningless -- the exact drift ld_common
# exists to prevent, so it is a hard stop, not a warning.
EMPIRICAL_SPEC = "4d1d1d92b25b"

# Bins at or above this bp are (a) resolved by the current STAGES at thin=25 and (b) inside the
# 324-generation forward window at r = 2.75e-6. See the module docstring.
DEFAULT_MIN_BIN = 562


def genome_indices(cluster_data, year):
    """Subpop index -> cluster row. Same construction as AnalyzeTreeSeq.py."""
    a = cluster_data[f"Genome Assignment {year}"]
    idx = [-1] * (int(max(a.dropna())) + 1)
    for i in range(len(a)):
        if not math.isnan(a[i]):
            k = int(a[i])
            if idx[k] == -1:
                idx[k] = i
    return idx


def site_names(year):
    """Site names in specifier-matrix row order (col 0) -- the canonical ordering (4)."""
    names = []
    with open(Path(f"../data/Genetic_Data/specifier_matrix_{year}.csv"), encoding="utf-8") as f:
        for line in f:
            if line.strip():
                names.append(line.split(",")[0].strip())
    return names


def read_empirical_ld(year):
    """averaged_ldDecay_{year}.csv -> (sum_r2, cnt, labels), both (N_BINS, K).

    Returns SUMS, not means: pooling across demes is sum(numerators)/sum(denominators), never a
    mean of means (5.2, invariant 4). write_decay stores the mean and the count, so the numerator
    is recovered as mean*count.
    """
    path = Path(f"../data/empiricalStats/averaged_ldDecay_{year}.csv")
    rows = list(csv.DictReader(open(path, encoding="utf-8")))
    if len(rows) != ldc.N_BINS:
        raise ValueError(f"{path}: {len(rows)} bins vs ld_common's {ldc.N_BINS} -- "
                         f"the empirical run used a different spec")
    labels = [c[3:] for c in rows[0] if c.startswith("r2_")]
    lo = np.array([int(r["bin_lo"]) for r in rows], dtype=np.int64)
    if not np.array_equal(lo, ldc.BIN_EDGES[:-1]):
        raise ValueError(f"{path}: bin edges differ from ld_common's -- specs have drifted")
    cnt = np.array([[float(r[f"n_{p}"]) for p in labels] for r in rows])
    mean = np.array([[float(r[f"r2_{p}"]) if r[f"r2_{p}"] != "" else 0.0 for p in labels]
                     for r in rows])
    return mean * cnt, cnt, labels


def pooled(sum_r2, cnt, keep):
    """Pair-count-weighted mean r^2 per bin over the KEPT demes. NaN where a bin has no pairs."""
    s = np.asarray(sum_r2)[:, keep].sum(axis=1)
    c = np.asarray(cnt)[:, keep].sum(axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(c > 0, s / np.where(c > 0, c, 1), np.nan), c


def run_slim(popmult, recomb, seed, raw_ts=None):
    """Run the forward sim, then move its output to `raw_ts` if one is given.

    SLiM's output path is HARDCODED in CPBSampleSim*.slim (out/simTreeSeq.trees), so
    concurrent SLiM runs would clobber each other. Hence: SLiM stays serial, and each point's
    tree is moved aside immediately so the EXPENSIVE half (recapitation) can then run in
    parallel from per-point files.
    """
    script = Path("../SLiM_Code/CPBSampleSim"
                  + ("Win" if platform.system() == "Windows" else "Linux") + ".slim")
    t0 = time.perf_counter()
    subprocess.run(["slim", "-l", "0", "-s", str(seed),
                    "-d", f"POPMULT={popmult}", "-d", f"RECOMB={recomb!r}",
                    str(script)], check=True)
    dt = time.perf_counter() - t0
    produced = Path("../out/simTreeSeq.trees")
    mb = produced.stat().st_size / 1024 / 1024
    if raw_ts and Path(raw_ts) != produced:
        Path(raw_ts).parent.mkdir(parents=True, exist_ok=True)
        shutil.move(str(produced), str(raw_ts))
    print(f"  [slim] seed={seed} {dt:.1f}s trees={mb:.1f}MB -> {raw_ts or produced}", flush=True)
    return dt


def build_mutated_ts(recomb, mu, seed, save_path=None, raw_ts=None):
    """Recapitate -> simplify -> mutate, mirroring AnalyzeTreeSeq.analyze_tree_sequence."""
    cluster_data = pd.read_csv(Path("../data/cluster_data.csv"))
    gi = {y: genome_indices(cluster_data, y) for y in TIMES}

    t0 = time.perf_counter()
    print(f"[{time.strftime('%H:%M:%S')}] recapitating anc_ne={ANCESTRAL_NE} ...", flush=True)
    ts = tskit.load(Path(raw_ts or "../out/simTreeSeq.trees"))
    ts = _ru.recapitate(ts, recombination_rate=recomb, ancestral_Ne=ANCESTRAL_NE,
                           random_seed=seed)
    print(f"[{time.strftime('%H:%M:%S')}] recap {time.perf_counter()-t0:.1f}s "
          f"edges={ts.num_edges}", flush=True)

    ksamp = []
    for y in TIMES:
        for i in gi[y]:
            ksamp.extend(ts.samples(population=i, time=TIMES[y]))
    # filter_populations=False is REQUIRED -- the ts.samples(population=i) queries use the
    # ORIGINAL cluster-row index (CLAUDE.md 2, commit c5963ae).
    ts = ts.simplify(samples=ksamp, filter_populations=False)
    gc.collect()
    print(f"[{time.strftime('%H:%M:%S')}] simplified -> {ts.num_samples} samples", flush=True)

    ts = msprime.sim_mutations(
        ts, rate=mu,
        model=msprime.SLiMMutationModel(type=0, next_id=pyslim.next_slim_mutation_id(ts)),
        keep=True, random_seed=seed)
    print(f"[{time.strftime('%H:%M:%S')}] mutated -> {ts.num_sites} sites "
          f"({time.perf_counter()-t0:.1f}s total)", flush=True)

    if save_path:
        ts.dump(Path(save_path))
    return ts, gi


def main(a):
    if ldc.spec_hash() != EMPIRICAL_SPEC:
        raise SystemExit(
            f"ld_common spec {ldc.spec_hash()} != {EMPIRICAL_SPEC}, the spec the empirical run "
            f"(ldCalcOut.txt) was computed under. The two curves are not comparable. Either "
            f"revert the change to ld_common, or re-run CalculateLD.py (~7.7 h) and update "
            f"EMPIRICAL_SPEC.")
    print(f"ld_common spec {ldc.spec_hash()} -- matches the empirical run")

    if a.slim_only:
        if not a.raw_ts:
            raise SystemExit("--slim-only needs --raw-ts to move the tree to")
        run_slim(a.popmult, a.recomb, a.seed, raw_ts=a.raw_ts)
        return

    rec = {"popmult": a.popmult, "mu": a.mu, "recomb": a.recomb, "anc_ne": ANCESTRAL_NE,
           "seed": a.seed, "spec": ldc.spec_hash(), "min_bin": a.min_bin,
           "thin": a.thin,
           "when": time.strftime("%Y-%m-%dT%H:%M:%S")}

    # ---- empirical targets, the real fitted mask, and the ordering assertion 7.5e wants -----
    obs, keep = {}, {}
    for y in TIMES:
        keep[y] = ABC.get_keep_mask(y)
        s, c, labels = read_empirical_ld(y)
        spec = site_names(y)
        if labels != spec:
            raise ValueError(
                f"{y}: averaged_ldDecay column order != specifier-matrix order.\n"
                f"  ld:        {labels[:3]} ...\n  specifier: {spec[:3]} ...\n"
                f"CLAUDE.md 7.5e checked these were identical on 2026-09-06; if the popfiles were "
                f"regenerated they no longer are, and the columns must be remapped.")
        if len(labels) != len(keep[y]):
            raise ValueError(f"{y}: {len(labels)} LD columns vs {len(keep[y])} specifier rows")
        obs[y] = (s, c)

    # ---- forward phase ----------------------------------------------------------------------
    t0 = time.perf_counter()
    if a.load_ts and Path(a.load_ts).exists():
        print(f"[{time.strftime('%H:%M:%S')}] loading {a.load_ts}", flush=True)
        ts = tskit.load(Path(a.load_ts))
        gi = {y: genome_indices(pd.read_csv(Path("../data/cluster_data.csv")), y) for y in TIMES}
    else:
        rec["slim_s"] = (None if a.skip_slim
                         else run_slim(a.popmult, a.recomb, a.seed, raw_ts=a.raw_ts))
        ts, gi = build_mutated_ts(a.recomb, a.mu, a.seed, save_path=a.save_ts,
                                  raw_ts=a.raw_ts)
    rec["raw_ts"] = a.raw_ts
    rec["setup_s"] = time.perf_counter() - t0
    rec["num_sites"] = int(ts.num_sites)
    rec["num_samples"] = int(ts.num_samples)

    # ---- the measurement --------------------------------------------------------------------
    outdir = Path("../out/ld_probe")
    outdir.mkdir(parents=True, exist_ok=True)
    fit = ldc.BIN_EDGES[:-1] >= a.min_bin
    rng = np.random.default_rng(a.seed)

    rec["years"], per_year_loss, ld_seconds = {}, [], 0.0
    for y in TIMES:
        t1 = time.perf_counter()
        # The REAL production function, not a reimplementation (CLAUDE.md 3).
        s_sim, c_sim = ATS.calculate_ld_decay(
            ts, gi[y], time=TIMES[y], year=y,
            # seed in the name: replicates at one POPMULT would otherwise overwrite each
            # other's per-deme CSVs (the jsonl records survive, the curves did not).
            output_path=outdir / f"ld_{y}_pop{a.popmult}_thin{a.thin}_s{a.seed}.csv", rng=rng,
            thin=a.thin)
        dt = time.perf_counter() - t1
        ld_seconds += dt

        m_sim, n_sim = pooled(s_sim, c_sim, keep[y])
        m_obs, n_obs = pooled(*obs[y], keep[y])
        ok = fit & np.isfinite(m_sim) & np.isfinite(m_obs)
        if not ok.any():
            raise ValueError(f"{y}: no fittable bins at or above {a.min_bin} bp")
        loss = float(np.abs(m_sim[ok] - m_obs[ok]).mean())
        per_year_loss.append(loss)

        rec["years"][y] = {
            "ld_s": dt, "loss": loss, "n_fit_bins": int(ok.sum()),
            "sim": [None if not np.isfinite(v) else float(v) for v in m_sim],
            "obs": [float(v) for v in m_obs],
            "sim_pairs": [int(v) for v in n_sim], "obs_pairs": [int(v) for v in n_obs]}

        print(f"\n[{y}]  {dt:.1f}s   fitted bins {int(ok.sum())}   "
              f"|sim-obs| over bins >= {a.min_bin} bp = {loss:.5f}")
        print(f"{'bin':>17} {'sim r2':>9} {'obs r2':>9} {'diff':>9}  {'sim pairs':>12}")
        for i in range(ldc.N_BINS):
            mark = "*" if ok[i] else " "
            sv = "--" if not np.isfinite(m_sim[i]) else f"{m_sim[i]:.4f}"
            dv = "--" if not np.isfinite(m_sim[i]) else f"{m_sim[i]-m_obs[i]:+.4f}"
            print(f"{mark}{ldc.BIN_EDGES[i]:>8}-{ldc.BIN_EDGES[i+1]:<7} {sv:>9} "
                  f"{m_obs[i]:>9.4f} {dv:>9}  {int(n_sim[i]):>12,}")
        ldc.report_halfway(s_sim[:, keep[y]], c_sim[:, keep[y]], prefix=f"    sim {y}: ")

    rec["ld_s_total"] = ld_seconds
    rec["ld_loss"] = float(np.mean(per_year_loss))
    print(f"\n=== POPMULT={a.popmult}: ld_loss (bins >= {a.min_bin} bp) = {rec['ld_loss']:.5f}")
    print(f"    LD cost {ld_seconds:.1f}s at thin={a.thin} for 3 years, "
          f"{len(gi['2015'])} demes "
          f"-- against a ~1.2 h trial that is {ld_seconds/3600/1.2*100:.1f}% overhead")

    JSONL.parent.mkdir(parents=True, exist_ok=True)
    with open(JSONL, "a", encoding="utf-8") as f:
        f.write(json.dumps(rec) + "\n")
    print(f"    appended to {JSONL}")


def summarize(a):
    if not JSONL.exists():
        raise SystemExit(f"{JSONL} not found -- run at least one point first")
    recs = [json.loads(l) for l in open(JSONL, encoding="utf-8") if l.strip()]
    specs = {r["spec"] for r in recs}
    if len(specs) > 1:
        raise SystemExit(f"records span multiple ld_common specs {specs} -- not comparable")
    recs.sort(key=lambda r: r["popmult"])
    print(f"{'POPMULT':>8} {'thin':>5} {'sites':>9} {'ld_s':>8} {'ld_loss':>9}  per-year loss")
    for r in recs:
        py = " ".join(f"{y}:{r['years'][y]['loss']:.5f}" for y in sorted(r["years"]))
        print(f"{r['popmult']:>8} {r.get('thin', 1):>5} {r['num_sites']:>9,} "
              f"{r['ld_s_total']:>8.1f} {r['ld_loss']:>9.5f}  {py}")

    if len(recs) < 2:
        print("\nonly one point -- the question this script asks needs at least two POPMULTs")
        return
    # THE decision number: does the fitted region move with POPMULT, and by how much against the
    # sim-obs gap it would be fitting?
    fit = ldc.BIN_EDGES[:-1] >= recs[0]["min_bin"]
    lo, hi = recs[0], recs[-1]
    print(f"\nmovement in the fitted bins (>= {recs[0]['min_bin']} bp):")
    for y in sorted(lo["years"]):
        a_ = np.array([np.nan if v is None else v for v in lo["years"][y]["sim"]], dtype=float)
        b_ = np.array([np.nan if v is None else v for v in hi["years"][y]["sim"]], dtype=float)
        o_ = np.array(lo["years"][y]["obs"], dtype=float)
        ok = fit & np.isfinite(a_) & np.isfinite(b_)
        move = float(np.abs(a_[ok] - b_[ok]).mean())
        gap = float(np.abs(a_[ok] - o_[ok]).mean())
        ratio = move / gap if gap else float("nan")
        print(f"  {y}: POPMULT {lo['popmult']}->{hi['popmult']} moves sim r2 by {move:.5f}; "
              f"sim-obs gap {gap:.5f}; ratio {ratio:.2f}")
    print("\nA ratio well under 1 means LD cannot discriminate POPMULT at this cut -- the curve "
          "sits at a fixed offset from the target and moving POPMULT barely shifts it.")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--popmult", type=int, default=5000)
    p.add_argument("--mu", type=float, default=ABC.DEFAULT_MUTATION_RATE)
    p.add_argument("--recomb", type=float, default=ABC.DEFAULT_RECOMBINATION_RATE)
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--min-bin", type=int, default=DEFAULT_MIN_BIN,
                   help="lowest bin_lo entering ld_loss (default 562, see module docstring)")
    p.add_argument("--thin", type=int, default=1,
                   help="keep every Nth site before pair enumeration. THE cost knob: cost is "
                        "quadratic in density, so thin=1 is ~141 s/deme (~3.9 h/trial) and "
                        "thin=25 is ~0.22 s/deme (~1.1 min/trial) at ~98 bp resolution -- ample "
                        "for the fitted bin_lo >= 562 cut. Unbiased, and NOT part of spec_hash, "
                        "so it does not invalidate the empirical run")
    p.add_argument("--skip-slim", action="store_true",
                   help="reuse ../out/simTreeSeq.trees as it stands -- you are asserting it was "
                        "produced at --popmult")
    p.add_argument("--raw-ts", help="per-point path for SLiM's RAW output. SLiM writes a fixed "
                   "path so its runs must be serial, but giving each point its own file lets the "
                   "recapitation half run in parallel. --slim-only stops after the forward phase")
    p.add_argument("--slim-only", action="store_true",
                   help="run SLiM, move the tree to --raw-ts, and exit (phase 1 of a parallel sweep)")
    p.add_argument("--save-ts", help="dump the mutated tree so re-runs cost seconds")
    p.add_argument("--load-ts", help="reuse a mutated tree from --save-ts")
    p.add_argument("--summarize", action="store_true")
    args = p.parse_args()
    summarize(args) if args.summarize else main(args)
