"""Re-score a finished batch's LD against ONE chromosome instead of the 17-chromosome pool.

THE QUESTION. LD was retired (CLAUDE.md 7.5) because the empirical target pools 17 chromosomes
whose effective recombination rate spans >=10x, while the simulation has one r. A mixture of
shifted curves is BROADER than any single curve, so no single r can match the pooled target.
One chromosome is (mostly) one r, so it removes most of that mismatch. Does LD then say anything
useful about POPMULT, and does it agree with F_c and F_st?

HOW THE CHROMOSOME IS CHOSEN -- blind to the simulation, and fixed before any scoring.
Every chromosome shares one Ne, so its LD rate rho_k = 4*Ne*r_k is proportional to its own r_k,
and its decay distance d_k is proportional to 1/r_k. The model's r is the linkage-map GENOME
AVERAGE (6.8.1). The chromosome whose decay distance sits at the MEDIAN of the 17 is the one
whose r is closest to typical, so it is the one the model's r is most defensible for. That rule
uses only how the chromosomes compare WITH EACH OTHER. Choosing a chromosome because it fits the
simulation well is the cherry-picking 7.5 rejected; `--all` prints every chromosome as a
sensitivity display, and must never be used to pick one.

Caveats the output cannot show:
  * The r-N trade-off is MOVED, not removed. LD constrains 4*N*r. If the chosen chromosome's
    true r is off the map average by k, the N it prefers is off by k. Read the answer as
    "N conditional on r", like 7.4.3's kernel.
  * The median is unweighted by chromosome length (lengths are not in this repo); the map r is a
    length-weighted average. Linked selection also moves rho between chromosomes, so
    "median rho" is only approximately "median r".
  * Each chromosome is itself a small mixture (r varies along it).

NO SIMULATION. Reads:
  out/batch3/abc_results.csv                                     parameters + the batch's losses
  out/batch3_raw/detailed_sim_results_<job>/run<iteration+1>/ld_<year>.csv   simulated curves
  data/ld_per_chr/chr<k>_<year>_ldDecay.csv                      observed, per chromosome
  out/fc_loss_floor/p<POPMULT>_*/ld_<year>.csv                   fixed-parameter replicates

CROSS-CHECK FIRST: scoring every trial against the genome-wide target must reproduce the batch's
own ld_loss column. If it does not, the trial->folder mapping or the pooling is wrong and nothing
after it is trustworthy, so the script stops.

Run from diagnostics/:  python ld_one_chrom.py [--all] [--chrom K]
"""
import argparse
import csv
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "Python_Code"))
import ABCAnalysisNoRedis as ABC  # noqa: E402  -- the REAL mask, reader, pooling rule, LD_MIN_BIN
import ld_common as ldc  # noqa: E402
from AnalyzeTreeSeq import _real_sample_sizes  # noqa: E402  -- for the 1/n_hap floor
from collect_batch import _rank, _zrank, PRIOR_POP  # noqa: E402

YEARS = ("2015", "2019", "2023")
N_CHR = 17
PER_CHR = Path("../data/ld_per_chr")
GENOME = Path("../data/empiricalStats")
FLOOR_DIR = Path("../out/fc_loss_floor")
LEVELS = (0.20, 0.10, 0.05, 0.02, 0.01)
PARAMS = ("pop", "total_migration", "m", "numClusters")

# decay-distance definition, identical to ld_chrom_check.test_shape (CLAUDE.md 7.5)
REF_LO = 100
FRAC = 0.25

FIT = ldc.BIN_EDGES[:-1] >= ABC.LD_MIN_BIN


# ------------------------------------------------------------------ curves

def pooled_curve(sum_, cnt, keep):
    """Demes pooled by pair count, sum(sum_r2)/sum(cnt) -- ABC._ld_pooled_mean_abs_diff's rule."""
    s = np.asarray(sum_)[:, keep].sum(axis=1)
    c = np.asarray(cnt)[:, keep].sum(axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(c > 0, s / np.where(c > 0, c, 1), np.nan)


def read_obs(path, year):
    """(sum, cnt) for an empirical LD file, with its column order ASSERTED against the specifier
    matrix -- the same guard getObservedData applies to the genome-wide file."""
    s, c, labels = ABC._read_ld_decay(path)
    spec = ABC._specifier_site_names(year)
    if labels != spec:
        raise ValueError(f"{path}: column order != specifier order ({labels[:3]} vs {spec[:3]})")
    return s, c


def obs_target(year, chrom):
    """Pooled observed curve for one year; chrom=None means the genome-wide (17-chr) target."""
    path = (GENOME / f"averaged_ldDecay_{year}.csv" if chrom is None
            else PER_CHR / f"chr{chrom}_{year}_ldDecay.csv")
    s, c = read_obs(path, year)
    return pooled_curve(s, c, ABC.get_keep_mask(year))


def loss(sim, obs):
    """ld_loss for many trials at once. sim: (N, 3, bins), obs: (3, bins) -> (N,).
    Mean |diff| over fitted bins finite on both sides, then mean over years."""
    ok = FIT[None, None, :] & np.isfinite(sim) & np.isfinite(obs)[None]
    d = np.where(ok, np.abs(sim - obs[None]), 0.0)
    with np.errstate(invalid="ignore"):
        per_year = d.sum(axis=2) / ok.sum(axis=2)
    return np.nanmean(per_year, axis=1)


def decay_distance(year, chrom):
    """Where the pooled excess over the 1/n_hap floor falls to FRAC of its REF_LO-bp value.
    Log-linear interpolation. Same definition as ld_chrom_check.test_shape."""
    path = PER_CHR / f"chr{chrom}_{year}_ldDecay.csv"
    s, c = read_obs(path, year)
    keep = ABC.get_keep_mask(year)
    pooled = pooled_curve(s, c, keep)
    floor_deme = 1.0 / (2 * np.asarray(_real_sample_sizes(int(year)), dtype=float))
    w = np.asarray(c)[:, keep]
    with np.errstate(invalid="ignore", divide="ignore"):
        fl = (floor_deme[keep] * w).sum(axis=1) / w.sum(axis=1)
    e = pooled - fl
    lo = ldc.BIN_EDGES[:-1]
    mid = np.sqrt(ldc.BIN_EDGES[:-1] * ldc.BIN_EDGES[1:])
    i0 = int(np.flatnonzero(lo == REF_LO)[0])
    target = FRAC * e[i0]
    for i in range(i0, len(e) - 1):
        if e[i] >= target > e[i + 1]:
            t = (e[i] - target) / (e[i] - e[i + 1])
            return float(np.exp(np.log(mid[i]) + t * (np.log(mid[i + 1]) - np.log(mid[i]))))
    return float("nan")


# ------------------------------------------------------------------ loading

def load_batch(results, raw):
    """Parameters, batch losses and pooled simulated curves for every trial with a raw folder."""
    with open(results, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    keep = {y: ABC.get_keep_mask(y) for y in YEARS}
    curves, kept, missing = [], [], 0
    for r in rows:
        d = raw / f"detailed_sim_results_{int(float(r['job_id']))}" / f"run{int(float(r['iteration'])) + 1}"
        paths = [d / f"ld_{y}.csv" for y in YEARS]
        if not all(p.exists() for p in paths):
            missing += 1
            continue
        yc = []
        for y, p in zip(YEARS, paths):
            s, c, _ = ABC._read_ld_decay(p)
            yc.append(pooled_curve(s, c, keep[y]))
        curves.append(yc)
        kept.append(r)
    A = {k: np.array([float(r[k]) for r in kept]) for k in
         PARAMS + ("ld_loss", "fst_loss", "fc_loss", "pi_loss")}
    return A, np.array(curves), missing, len(rows)


def load_floor():
    """{POPMULT: (reps, 3, bins)} from the fixed-parameter replicates (7.9.10C)."""
    out = {}
    keep = {y: ABC.get_keep_mask(y) for y in YEARS}
    for d in sorted(FLOOR_DIR.glob("p*_*")):
        paths = [d / f"ld_{y}.csv" for y in YEARS]
        if not all(p.exists() for p in paths):
            continue
        pm = int(d.name.split("_")[0][1:])
        yc = []
        for y, p in zip(YEARS, paths):
            s, c, _ = ABC._read_ld_decay(p)
            yc.append(pooled_curve(s, c, keep[y]))
        out.setdefault(pm, []).append(yc)
    return {k: np.array(v) for k, v in out.items()}


# ------------------------------------------------------------------ read-outs

def accept(pop, dist):
    """Median pop and IQR ratio (vs the uniform prior's IQR) at each acceptance level."""
    prior_iqr = 0.5 * (PRIOR_POP[1] - PRIOR_POP[0])
    order = np.argsort(dist, kind="mergesort")
    out = []
    for lv in LEVELS:
        k = max(1, int(round(lv * len(dist))))
        p = pop[order[:k]]
        q1, med, q3 = np.percentile(p, [25, 50, 75])
        out.append((lv, k, float(med), float((q3 - q1) / prior_iqr)))
    return out


def decile_argmin(pop, dist):
    """Median loss per pop decile; index of the lowest, and the medians."""
    edges = np.percentile(pop, np.linspace(0, 100, 11))
    med = []
    for i in range(10):
        m = (pop >= edges[i]) & (pop <= edges[i + 1] if i == 9 else pop < edges[i + 1])
        med.append(float(np.median(dist[m])))
    return int(np.argmin(med)), med, edges


def unique_r2(A, y):
    """Rank-space unique R2 of each parameter -- collect_batch.report_decomposition's method."""
    X = np.column_stack([np.ones(len(y))] + [_zrank(A[p]) for p in PARAMS])
    yz = _zrank(y)
    sst = float(np.sum((yz - yz.mean()) ** 2))
    b, *_ = np.linalg.lstsq(X, yz, rcond=None)
    r2 = 1 - float(np.sum((yz - X @ b) ** 2)) / sst
    uniq = {}
    for k, p in enumerate(PARAMS):
        cols = [0] + [j + 1 for j in range(len(PARAMS)) if j != k]
        b2, *_ = np.linalg.lstsq(X[:, cols], yz, rcond=None)
        uniq[p] = r2 - (1 - float(np.sum((yz - X[:, cols] @ b2) ** 2)) / sst)
    return r2, uniq


def floor_ratio(floor, obs):
    """Replicate sd of the loss at each POPMULT, and pooled sd / |mean(2000) - mean(5000)|."""
    stats = {pm: loss(c, obs) for pm, c in floor.items()}
    if not {2000, 5000} <= set(stats):
        return stats, float("nan")
    a, b = stats[2000], stats[5000]
    sd = np.sqrt(0.5 * (a.var(ddof=1) + b.var(ddof=1)))
    return stats, float(sd / abs(a.mean() - b.mean()))


def report(name, A, dist, floor=None, obs=None):
    print(f"\n--- {name} ---")
    r2, uq = unique_r2(A, dist)
    print(f"  rank R2 total {r2:.3f}   unique: " +
          "  ".join(f"{p} {uq[p]:+.3f}" for p in PARAMS))
    for lv, k, med, iqr in accept(A["pop"], dist):
        print(f"  top {lv:>4.0%} (n={k:4d})  pop median {med:7.0f}   IQR ratio {iqr:.2f}")
    j, med, edges = decile_argmin(A["pop"], dist)
    print(f"  lowest median loss in pop decile {j + 1} ({edges[j]:.0f}-{edges[j + 1]:.0f})"
          f"   [decile medians {med[0]:.5f} ... {med[-1]:.5f}]")
    rec = {"r2": r2, "unique": uq, "accept": accept(A["pop"], dist), "argmin_decile": j + 1}
    if floor is not None and obs is not None:
        stats, ratio = floor_ratio(floor, obs)
        for pm in sorted(stats):
            print(f"  fixed-param replicates, POPMULT {pm}: {stats[pm].mean():.5f} "
                  f"+/- {stats[pm].std(ddof=1):.5f} (n={len(stats[pm])})")
        print(f"  floor/signal (pooled sd / |2000-5000 difference|): {ratio:.2f}")
        rec["floor_signal"] = ratio
    return rec


# ------------------------------------------------------------------ main

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default="../out/batch3/abc_results.csv")
    ap.add_argument("--raw", default="../out/batch3_raw")
    ap.add_argument("--chrom", type=int, default=None,
                    help="override the median-decay rule (sensitivity only; state why)")
    ap.add_argument("--all", action="store_true",
                    help="also score every chromosome -- a sensitivity display, NOT a way to choose")
    ap.add_argument("--log", default="../out/ld_one_chrom.jsonl")
    a = ap.parse_args()

    print(f"ld_common spec {ldc.spec_hash()} (empirical targets: {ABC.LD_EMPIRICAL_SPEC})   "
          f"LD_MIN_BIN {ABC.LD_MIN_BIN} bp ({int(FIT.sum())} fitted bins)")
    if ldc.spec_hash() != ABC.LD_EMPIRICAL_SPEC:
        sys.exit("spec mismatch -- the two sides are binned differently")

    # 1. choose the chromosome, before any trial is read
    print("\n=== 1. CHROMOSOME SELECTION (simulation-blind) ===")
    dd = np.array([[decay_distance(y, k) for y in YEARS] for k in range(1, N_CHR + 1)])
    gm = np.exp(np.nanmean(np.log(dd), axis=1))
    order = np.argsort(gm)
    med_chr = int(order[N_CHR // 2]) + 1
    print("  decay distance (bp), geometric mean over years, fastest -> slowest:")
    for i in order:
        tag = "  <- median" if i + 1 == med_chr else ""
        print(f"    chr{i + 1:<3d} {gm[i]:9,.0f}   ({', '.join(f'{v:,.0f}' for v in dd[i])}){tag}")
    print(f"  2015 range {np.nanmin(dd[:, 0]):,.0f} -> {np.nanmax(dd[:, 0]):,.0f} bp "
          f"(CLAUDE.md 7.5 records 1,985 -> 26,063)")
    chosen = a.chrom or med_chr
    if a.chrom:
        print(f"  OVERRIDE: using chr{chosen} instead of the median chr{med_chr}")

    # 2. load and cross-check
    print("\n=== 2. LOAD + CROSS-CHECK ===")
    t0 = time.time()
    A, sim, missing, total = load_batch(Path(a.results), Path(a.raw))
    print(f"  {len(A['pop'])} of {total} trials have LD curves ({missing} missing), "
          f"{time.time() - t0:.0f} s")
    genome = np.array([obs_target(y, None) for y in YEARS])
    re = loss(sim, genome)
    err = np.nanmax(np.abs(re - A["ld_loss"]))
    print(f"  re-scored genome-wide ld_loss vs the batch's own column: max |diff| {err:.2e}")
    if not err < 1e-9:
        sys.exit("CROSS-CHECK FAILED -- trial->folder mapping or pooling is wrong; stopping")
    floor = load_floor()
    print(f"  fixed-parameter replicates: " +
          ", ".join(f"POPMULT {k}: {len(v)}" for k, v in sorted(floor.items())))

    # 3. read-outs
    print("\n=== 3. READ-OUTS ===")
    print("  IQR ratio = accepted pop IQR / prior IQR (1.0 = learned nothing). Single-statistic"
          "\n  rejection, no weights -- this asks what each statistic says ALONE.")
    log = {"stamp": time.strftime("%Y%m%d-%H%M%S"), "results": a.results,
           "ld_spec": ldc.spec_hash(), "ld_min_bin": ABC.LD_MIN_BIN,
           "median_chr": med_chr, "chosen_chr": chosen,
           "decay_gm": {f"chr{k + 1}": float(gm[k]) for k in range(N_CHR)}, "readouts": {}}
    log["readouts"]["genome"] = report("LD, genome-wide 17-chr pool (the retired target)",
                                       A, re, floor, genome)
    one = np.array([obs_target(y, chosen) for y in YEARS])
    log["readouts"][f"chr{chosen}"] = report(f"LD, chr{chosen} only", A, loss(sim, one), floor, one)
    for s in ("fc_loss", "fst_loss"):
        log["readouts"][s] = report(f"{s} alone, for comparison", A, A[s])

    if a.all:
        print("\n=== 4. EVERY CHROMOSOME -- SENSITIVITY ONLY, DO NOT CHOOSE FROM THIS ===")
        print("  If the preferred pop moves with decay distance, that is the r-N trade-off:")
        print("  slow-decay (low-r) chromosomes should prefer SMALLER pop.")
        print(f"  {'chr':>5s} {'decay bp':>9s} {'top5% pop':>10s} {'top1% pop':>10s} "
              f"{'argmin dec':>10s} {'floor/sig':>9s}")
        per = {}
        for i in order:
            k = i + 1
            obs = np.array([obs_target(y, k) for y in YEARS])
            dist = loss(sim, obs)
            acc = accept(A["pop"], dist)
            j, _, _ = decile_argmin(A["pop"], dist)
            _, ratio = floor_ratio(floor, obs)
            per[f"chr{k}"] = {"accept": acc, "argmin_decile": j + 1, "floor_signal": ratio}
            print(f"  chr{k:<2d} {gm[i]:9,.0f} {acc[2][2]:10.0f} {acc[4][2]:10.0f} "
                  f"{j + 1:10d} {ratio:9.2f}")
        rho = np.corrcoef(_rank(np.log(gm[order])),
                          _rank(np.array([per[f'chr{i + 1}']['accept'][2][2] for i in order])))[0, 1]
        print(f"  spearman(decay distance, top-5% pop median) = {rho:+.2f}")
        log["per_chr"] = per

    with open(a.log, "a", encoding="utf-8") as f:
        f.write(json.dumps(log) + "\n")
    print(f"\nappended one record to {a.log}")


if __name__ == "__main__":
    main()
