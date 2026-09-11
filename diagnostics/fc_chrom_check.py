"""Between-chromosome tests on the empirical temporal F_c, plus where the target lands.

This is 7.8.4's tool applied to F_c instead of LD, and it is the first thing the per-chromosome
files (data/Fc_per_chr/, 17 chr x 19 field pairs) are for. It runs in seconds and simulates
nothing.

WHY IT MATTERS BEFORE ANY PILOT BATCH. F_c passed its selectivity test (7.9.6: F_st's implied N
swings 5.3x across the kernel prior, F_c's swings 1.07-1.22x) and its replicate noise floor
(7.9.8: sd 0.00124 at POPMULT=2000, 11% of the POPMULT 2000->5000 signal). Neither of those says
anything about the EMPIRICAL side, and LD died on exactly that -- ld_loss cleared its own noise
floor at 2.8% and was still worthless, because the target was a 17-chromosome mixture the
single-r model could not be (7.5). Three questions, all answerable here:

  TEST 0 -- integrity. Sum sum_fc over chromosomes, divide by summed n_loci, and check it
  reproduces the averaged file pair by pair. Pooling, not averaging (invariant 4). If this fails
  the two files disagree about what they hold and nothing below means anything.

  TEST 1 -- how precise is the target? Leave-one-chromosome-out jackknife on the GAP-POOLED mean,
  which is the exact quantity _fc_loss differences. Read the SE against two things: the replicate
  noise floor (0.00124), and the POPMULT 2000->5000 signal (0.01119). An empirical SE comparable
  to the signal would end F_c the way the mixture ended LD -- the target would not be known well
  enough to fit.

  TEST 2 -- is F_c homogeneous across chromosomes? It SHOULD be, and that is the interesting part.
  F_c is a frequency change over 8-16 generations; r sets only how many effectively independent
  loci a chromosome contributes, i.e. PRECISION, not expectation. So the >=10x between-chromosome
  spread in effective r that killed LD (7.5) should NOT produce a systematic F_c spread here.
  If it does anyway -- and especially if the ranking repeats the LD one (chr6, chr5, chr15) --
  then F_c has inherited the same mixture problem and the pilot batch is not worth running.

  TEST 4 -- level or slope? F_c ~ pedestal + t/(2N), and we have TWO gaps, so the level and the
  drift rate can be read separately. A t-independent contaminant on the empirical side (genotyping
  error between timepoints, within-field family structure -- 7.2 shows some sites are samples of
  close relatives) shifts the LEVEL of both gaps equally and leaves the CONTRAST alone. This is the
  one thing LD could never do: every LD distance bin was contaminated the same way, so there was no
  internal contrast to fall back on.

  TEST 3 -- where does the target sit relative to the simulated range? Reads the n_real arm of
  out/temporal_fc.jsonl (spec-matched, Q=100 constants) and prints fc_loss per POPMULT. LD's fatal
  sign was a target OUTSIDE the prior, with ld_loss monotone across all of it. Bracketing is what
  we need to see.

Run from diagnostics/:  python fc_chrom_check.py
"""
import argparse
import csv
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "Python_Code"))
import fc_common as fcc  # noqa: E402  -- the SHARED spec, never a reimplementation

PER_CHR = Path("../data/Fc_per_chr")
POOLED = Path("../data/empiricalStats/averaged_temporalFc.csv")
SIM_JSONL = Path("../out/temporal_fc.jsonl")

# 7.9.8: replicate noise floor on the pooled statistic, and the POPMULT 2000->5000 signal,
# both n_real at the Q=100 constants. Quoted here so the ratios below need no lookup.
FLOOR_SD = 0.00124
POPMULT_SIGNAL = 0.01119


def read_per_chr(d):
    """{chr: {(a, b): (t, sum_fc, n_loci)}} -- raw sums, never per-chromosome ratios."""
    out = {}
    for p in sorted(d.glob("chr*_temporalFc.csv")):
        chrom = p.name.split("_")[0]
        rows = {}
        with open(p, newline="", encoding="utf-8") as f:
            for r in csv.DictReader(f):
                rows[(r["site_a"].strip(), r["site_b"].strip())] = (
                    int(r["t"]), float(r["sum_fc"]), int(r["n_loci"]))
        out[chrom] = rows
    if not out:
        raise SystemExit(f"no per-chromosome files under {d}")
    return out


def read_pooled(p):
    """{(a, b): (t, fc, pedestal, n_loci)} from the averaged target."""
    out = {}
    with open(p, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            out[(r["site_a"].strip(), r["site_b"].strip())] = (
                int(r["t"]), float(r["fc"]), float(r["pedestal"]), int(r["n_loci"]))
    return out


def gap_pooled(per_chr, chroms):
    """Pool the named chromosomes, then mean over pairs within each gap.

    This is exactly what _fc_loss does to one side: sum/sum per pair (invariant 4), plain mean
    over the pairs of a gap. Anything else here would be measuring a different statistic than the
    one the fit uses.
    """
    num, den, gap_of = defaultdict(float), defaultdict(int), {}
    for c in chroms:
        for key, (t, s, n) in per_chr[c].items():
            num[key] += s
            den[key] += n
            gap_of[key] = t
    by_gap = defaultdict(list)
    for key in num:
        by_gap[gap_of[key]].append(num[key] / den[key])
    return {t: float(np.mean(v)) for t, v in by_gap.items()}, \
           {key: num[key] / den[key] for key in num}


def per_chrom_gap(per_chr, chrom):
    by_gap = defaultdict(list)
    for key, (t, s, n) in per_chr[chrom].items():
        by_gap[t].append(s / n)
    return {t: float(np.mean(v)) for t, v in by_gap.items()}


def sim_by_popmult(path, spec):
    """(means, skipped, per-seed) from the n_real arm, spec-matched.

    The per-seed dict is what TEST 4 needs: the slope is a DIFFERENCE of two gaps, and its
    replicate noise has to be computed seed by seed, not from the two averaged numbers.
    """
    if not path.exists():
        return {}, [], {}
    recs, skipped = [], []
    for line in path.open(encoding="utf-8"):
        line = line.strip()
        if not line:
            continue
        r = json.loads(line)
        if r.get("arm") != "n_real":
            continue
        if r.get("fc_spec") != spec:
            skipped.append((r.get("popmult"), r.get("seed"), r.get("fc_spec")))
            continue
        recs.append(r)
    out = defaultdict(lambda: defaultdict(list))
    per_seed = defaultdict(dict)                 # {popmult: {(gap, seed): fc}}
    for r in recs:
        for gk, v in r["res"].items():
            if not re.fullmatch(r"t\d+", gk):   # res also carries non-gap keys
                continue
            out[float(r["popmult"])][int(gk[1:])].append(v["fc_mean"])
            per_seed[float(r["popmult"])][(int(gk[1:]), r.get("seed"))] = v["fc_mean"]
    return ({p: {t: (float(np.mean(v)), len(v)) for t, v in g.items()}
             for p, g in out.items()}, skipped, per_seed)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--per-chr", default=str(PER_CHR))
    ap.add_argument("--pooled", default=str(POOLED))
    ap.add_argument("--sim", default=str(SIM_JSONL))
    a = ap.parse_args()

    per_chr = read_per_chr(Path(a.per_chr))
    pooled = read_pooled(Path(a.pooled))
    chroms = sorted(per_chr, key=lambda c: int(c[3:]))
    print(f"fc_common spec {fcc.spec_hash()}")
    print(f"{len(chroms)} chromosomes, {len(pooled)} field pairs\n")

    # ---- TEST 0: pooling integrity -----------------------------------------
    all_gap, all_pair = gap_pooled(per_chr, chroms)
    worst, worst_key = 0.0, None
    for key, (t, fc, _ped, n_loci) in pooled.items():
        if key not in all_pair:
            raise SystemExit(f"{key} is in the pooled file but not in the per-chromosome files")
        d = abs(all_pair[key] - fc)
        if d > worst:
            worst, worst_key = d, key
        n_sum = sum(per_chr[c][key][2] for c in chroms)
        if n_sum != n_loci:
            print(f"  ! n_loci mismatch for {key}: per-chr sum {n_sum} vs pooled {n_loci}")
    print(f"TEST 0  pooling integrity: max |per-chr pooled - averaged file| = {worst:.3e} "
          f"({worst_key[0]} vs {worst_key[1]})")
    print("        " + ("PASS -- the two files agree" if worst < 1e-6 else
                        "FAIL -- the files disagree; stop here"))
    print()

    # ---- TEST 1: jackknife --------------------------------------------------
    gaps = sorted(all_gap)
    print("TEST 1  leave-one-chromosome-out jackknife on the gap-pooled mean")
    print(f"{'gap':>5} {'pooled F_c':>12} {'jack SE':>10} {'SE/floor':>9} {'SE/signal':>10}")
    jack_se = {}
    for t in gaps:
        loo = [gap_pooled(per_chr, [c for c in chroms if c != drop])[0][t] for drop in chroms]
        m = float(np.mean(loo))
        n = len(chroms)
        se = float(np.sqrt((n - 1) / n * np.sum((np.array(loo) - m) ** 2)))
        jack_se[t] = se
        print(f"{t:>5} {all_gap[t]:>12.6f} {se:>10.6f} "
              f"{se / FLOOR_SD:>9.2f} {se / POPMULT_SIGNAL:>10.1%}")
    print(f"        floor sd {FLOOR_SD} and POPMULT 2000->5000 signal {POPMULT_SIGNAL} are 7.9.8's,")
    print("        n_real arm at the Q=100 constants.")
    print()

    # ---- TEST 2: heterogeneity ---------------------------------------------
    print("TEST 2  per-chromosome gap-pooled F_c (is the target a mixture, as LD's was?)")
    head = "  ".join(f"t{t:<9}" for t in gaps)
    print(f"{'chr':>6}  {head}")
    per_c = {c: per_chrom_gap(per_chr, c) for c in chroms}
    for c in chroms:
        print(f"{c:>6}  " + "  ".join(f"{per_c[c][t]:<10.6f}" for t in gaps))
    print()
    for t in gaps:
        vals = np.array([per_c[c][t] for c in chroms])
        # If chromosomes differed only by sampling, their spread would be ~sqrt(n)*jackknife SE.
        expected = jack_se[t] * np.sqrt(len(chroms))
        print(f"   gap {t}: between-chromosome sd {vals.std(ddof=1):.6f}, "
              f"jackknife-implied {expected:.6f}, ratio {vals.std(ddof=1) / expected:.2f}")
        order = [c for _, c in sorted(zip(vals, chroms), reverse=True)]
        print(f"           highest -> lowest: {' '.join(order)}")
    print("   ratio ~1 means the chromosomes are exchangeable, i.e. NOT the LD mixture problem.")
    print("   A ratio >>1 with the SAME ranking in both gaps is the 7.5 failure recurring.")
    print()

    # ---- TEST 3: where the target sits -------------------------------------
    sim, skipped, sim_seeds = sim_by_popmult(Path(a.sim), fcc.spec_hash())
    if skipped:
        print(f"   (skipped {len(skipped)} n_real records at a different fc_spec)")
    if not sim:
        print("TEST 3  no spec-matched n_real records in the simulated jsonl -- skipped")
        return
    print("TEST 3  simulated (n_real) against the empirical target, and the implied fc_loss")
    print(f"{'POPMULT':>8} {'seeds':>6} " + " ".join(f"{'sim t' + str(t):>12}" for t in gaps) +
          f" {'fc_loss':>10}")
    for p in sorted(sim):
        row, per_gap_abs = [], []
        seeds = 0
        for t in gaps:
            if t not in sim[p]:
                row.append("        --  ")
                continue
            v, k = sim[p][t]
            seeds = max(seeds, k)
            row.append(f"{v:>12.6f}")
            per_gap_abs.append(abs(v - all_gap[t]))
        loss = float(np.mean(per_gap_abs)) if per_gap_abs else float("nan")
        print(f"{p:>8.0f} {seeds:>6} " + " ".join(row) + f" {loss:>10.6f}")
    print("   observed  " + " " * 6 + " ".join(f"{all_gap[t]:>12.6f}" for t in gaps))
    print()
    print("   READ THE SIGN. If the simulated values BRACKET the observed one, F_c has an interior")
    print("   optimum inside the prior and the pilot batch is worth running. If they are all on the")
    print("   same side and monotone, that is LD's failure mode (7.5) and the optimum is outside")
    print("   the prior -- find out why before spending a batch.")
    print()

    # ---- TEST 4: split the LEVEL from the SLOPE -----------------------------
    if len(gaps) < 2:
        return
    lo, hi = gaps[0], gaps[-1]
    print(f"TEST 4  the two gaps are two equations: LEVEL offset vs drift SLOPE (t{hi} - t{lo})")
    print("   F_c ~ pedestal + t/(2N). The pedestal is identical on both sides by construction")
    print("   (same n_i), so it cancels in obs-sim -- but ANY t-independent contaminant on the")
    print("   empirical side (genotyping error, within-field family structure) does not. The gap")
    print("   CONTRAST removes every t-independent term from both sides at once, so it reads the")
    print("   drift rate alone. Its cost is noise: it differences two noisy numbers, and the two")
    print("   gaps are largely DISJOINT field sets, so it assumes the contaminant is the same in")
    print("   both groups. Read it as a second opinion on the level, never as a replacement.")
    loo_gaps = [gap_pooled(per_chr, [c for c in chroms if c != d])[0] for d in chroms]
    n = len(chroms)
    sl = np.array([g[hi] - g[lo] for g in loo_gaps])
    obs_slope = all_gap[hi] - all_gap[lo]
    obs_se = float(np.sqrt((n - 1) / n * np.sum((sl - sl.mean()) ** 2)))
    print(f"\n   observed slope {obs_slope:.6f}  (jackknife SE {obs_se:.6f})")
    print(f"{'POPMULT':>8} {'slope':>10} {'sd(seeds)':>11} {'level offset t' + str(lo):>17} "
          f"{'level offset t' + str(hi):>17}")
    fit = []
    for p in sorted(sim):
        if lo not in sim[p] or hi not in sim[p]:
            continue
        s_lo, _ = sim[p][lo]
        s_hi, _ = sim[p][hi]
        paired = [sim_seeds[p][(hi, sd)] - sim_seeds[p][(lo, sd)]
                  for (g, sd) in sim_seeds[p] if g == lo and (hi, sd) in sim_seeds[p]]
        sd_txt = f"{np.std(paired, ddof=1):.6f}" if len(paired) > 1 else "--"
        print(f"{p:>8.0f} {s_hi - s_lo:>10.6f} {sd_txt:>11} "
              f"{all_gap[lo] - s_lo:>+17.6f} {all_gap[hi] - s_hi:>+17.6f}")
        fit.append((p, s_hi - s_lo))
    print("   A level offset that is the SAME in both gaps is a pure t-independent contaminant.")
    print("   One that DIFFERS between gaps is a genuine drift-rate mismatch.")
    if len(fit) >= 2 and all(v > 0 for _, v in fit):
        (p0, s0), (p1, s1) = fit[0], fit[-1]
        b = (np.log(s1) - np.log(s0)) / (np.log(p1) - np.log(p0))
        a = np.log(s0) - b * np.log(p0)
        print(f"\n   simulated slope ~ POPMULT^{b:.3f}   (drift theory says -1)")
        for tgt, lab in ((obs_slope, "observed"), (obs_slope + obs_se, "obs +1 SE"),
                         (max(obs_slope - obs_se, 1e-12), "obs -1 SE")):
            print(f"   POPMULT matching the {lab:<10} slope: {np.exp((np.log(tgt) - a) / b):>10.0f}")
        print("   Two points, so this is a two-point slope extrapolated -- a bearing, not a fix.")


if __name__ == "__main__":
    main()
