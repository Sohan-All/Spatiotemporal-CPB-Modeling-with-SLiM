"""Temporal Fc probe -- CLAUDE.md 7.9. THE GATE on the temporal statistic.

THE QUESTION. F_st identifies only the product N*m (CLAUDE.md 7.4.2), and 7.9's information-budget
argument says why: only 0.5-7% of the genealogy sits inside the POPMULT-dependent forward window,
and that sliver's size is itself set by the SCATTERING phase, whose only parameter is N*m. Every
statistic computed on that tree therefore inherits the confound. Temporal Fc is the one candidate
that escapes it, because it measures a CHANGE:

  * The ancestral phase is shared history between two timepoints, so it cancels EXACTLY. Fc is
    100% forward-phase, against 0.5-7% for pi and F_st. It escapes the budget rather than living
    inside it.
  * It is mu-free and denominator-free (a variance of frequencies), so 5.1's callable-site problem
    and the mu nuisance that dominates pi_loss (7.6.1) both vanish.
  * It is r-free in EXPECTATION -- recombination sets how many effectively independent loci you
    have, i.e. precision, not the value. It does not walk into 7.8's wall.
  * It reads an ABSOLUTE DRIFT RATE over a known number of generations, ~ t/(2N), rather than a
    coalescence-vs-migration probability. Migration enters through reversion to the metapopulation
    mean, ~(1-m)^t, which is a function of m ALONE, not of N*m. Non-parallel contours to F_st's is
    exactly what identification requires.

WHY 11's REJECTION OF THE TEMPORAL METHOD DOES NOT APPLY. That rejection was: at n=5-8 the
sampling correction 1/(2*S_a) + 1/(2*S_b) ~ 0.14 is 11-21x the drift signal. Decisive if you are
CORRECTING it to get a point estimate of Ne. In ABC you do not correct it, you REPRODUCE it --
subsample the simulated demes to the same n_i (invariant 1) and the identical offset appears on
both sides and cancels. What survives is only a precision question, and precision is bought with
loci (Waples 1989: doubling S, t or L buys about the same). This probe measures both arms so the
cost of that pedestal is a number rather than an argument.

WHAT IT MEASURES. 19 field-pairs are resampled across years at IDENTICAL coordinates (checked from
the specifier matrices): 9 at t=16 generations (2015 vs 2023), 8 at t=8 (2019 vs 2023), 2 at t=8
(2015 vs 2019). The simulation already Remembers every individual at times 16/8/0, so nothing about
the pipeline changes -- this runs on a tree sequence that already exists.

Two arms per point:
  n_real : each timepoint cut to the field's REAL n_i (5-8 diploids). What you could actually
           measure. Carries the ~0.14 sampling pedestal.
  n_big  : both timepoints cut to --big-n diploids. The pedestal shrinks to ~1/big_n, so this is
           the clean read of whether Fc tracks N at all. NOT a proposal for production -- the
           empirical side cannot do it -- it is the diagnostic that separates "the statistic has
           no signal" from "the signal is buried under the pedestal".

THE DECISION NUMBER is the contrast between POPMULT points relative to the replicate spread, in
each arm. Theory to check against: whole-deme Fc ~ t/(2*N_deme), so POPMULT 2000 -> 5000 (deme
N 202 -> 505) should roughly halve it in the n_big arm.

USAGE (the trees in out/ from the 7.5.4 LD sweep are ALREADY recapitated+simplified+mutated, so
this costs seconds, not the ~30 min a fresh point would):
    python temporal_fc.py --load-ts ../out/ld_probe_p2000.trees --popmult 2000 --reps 5
    python temporal_fc.py --load-ts ../out/ld_probe_p5000.trees --popmult 5000 --reps 5
    python temporal_fc.py --summarize

From a raw SLiM tree instead: --raw-ts <path> recapitates, simplifies and mutates first.
"""

import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent.parent / 'Python_Code'))
import scale_constants as _sc
import recapitate_util as _ru

import argparse
import json
import math
import sys
import time
from pathlib import Path

import numpy as np
import tskit

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "Python_Code"))
import ABCAnalysisNoRedis as ABC  # noqa: E402  (loss/statistics half only -- no sim stack)
import fc_common as fcc           # noqa: E402  the SHARED spec -- see fc_common.py

DATA = Path(__file__).resolve().parent.parent / "data"
OUT_JSONL = Path(__file__).resolve().parent.parent / "out" / "temporal_fc.jsonl"

# THE SPEC LIVES IN fc_common.py and is shared with the empirical side (ToUseOnBeagles/).
# Re-exported here under the old names. Do NOT redefine either of these locally -- a divergence
# between the two sides is silent and makes the F_c values incomparable, which is exactly the
# failure ld_common.py was created to prevent.
YEAR_TIME = fcc.YEAR_TIME
MIN_MAF = fcc.MIN_MAF


# ---------------------------------------------------------------------------
# The resampled fields
# ---------------------------------------------------------------------------

def _spec_path(year):
    return DATA / "Genetic_Data" / f"specifier_matrix_{year}.csv"


def _pop_path(year):
    return DATA / "Genetic_Data" / f"popFile{year}"


def _cluster_of_specifier_row(year):
    """specifier row -> cluster row, built the same way analyze_tree_sequence builds it.

    cluster_data.csv has one row per CLUSTER; 'Genome Assignment {year}' names the specifier row
    that cluster carries. Inverting it gives specifier row -> cluster row, which is the index
    ts.samples(population=...) wants (and why filter_populations=False is mandatory, CLAUDE.md 2).
    """
    import pandas as pd
    cd = pd.read_csv(DATA / "cluster_data.csv")
    col = cd[f"Genome Assignment {year}"]
    out = [-1] * (int(max(col.dropna())) + 1)
    for i in range(len(col)):
        if not math.isnan(col[i]):
            j = int(col[i])
            if out[j] == -1:
                out[j] = i
    if -1 in out:
        raise ValueError(f"{year}: specifier row {out.index(-1)} has no cluster assignment")
    return out


def _haversine_m(lat1, lon1, lat2, lon2):
    """Metres. Same formula as GenerateClusterData.distance, without its int() truncation."""
    p = math.pi / 180
    a = (0.5 - math.cos((lat2 - lat1) * p) / 2
         + math.cos(lat1 * p) * math.cos(lat2 * p) * (1 - math.cos((lon2 - lon1) * p)) / 2)
    return 2 * 6371000 * math.asin(math.sqrt(a))


def _nearest_cluster(lat, lon):
    import pandas as pd
    cd = pd.read_csv(DATA / "cluster_data.csv")
    d = [_haversine_m(lat, lon, la, lo) for la, lo in zip(cd["Latitude"], cd["Longitude"])]
    k = int(np.argmin(d))
    return k, d[k]


def matched_pairs(deme_choice="nearest", verbose=True):
    """Fields resampled in two years, matched on IDENTICAL coordinates.

    Coordinates, not names: the same field is 'BrilowskiHome-2015' and 'BrikalskiPats-2023', and
    'GarrisonNE-2015' is 'OkrayGarrisonNE-2019'. Matching on names would lose those silently.

    WHICH SIMULATED DEME STANDS FOR THE FIELD -- and this is NOT a detail. The production mapping
    (GenerateClusterData.assign_genomes_to_clusters_idv_year) is greedy in specifier-row order and
    ONE-TO-ONE per year: each site takes its nearest cluster that no earlier site of that year has
    already claimed. So THE SAME PHYSICAL FIELD IS A DIFFERENT DEME IN DIFFERENT YEARS -- measured,
    Arlington is cluster 10 in 2015 and cluster 21 in 2019. Computing Fc across those two would
    measure divergence between two different demes, not drift within one, so it cannot be used
    here. See CLAUDE.md 7.9.2.

      deme_choice="nearest" (default): the cluster whose centroid is nearest the field, identical
        at both timepoints and independent of year. The physically right stand-in, and the only
        choice under which Fc means what it is supposed to mean.
      deme_choice="later": the cluster production assigns in the later year. Kept so the
        sensitivity of the result to this choice can be checked.
    """
    # The field pairs themselves come from the SHARED spec, so the empirical side scores the
    # identical set. Only the deme attribution below is simulation-side.
    pairs, n_masked = fcc.matched_field_pairs(
        {y: _spec_path(y) for y in YEAR_TIME}, {y: _pop_path(y) for y in YEAR_TIME})
    clus = {y: _cluster_of_specifier_row(y) for y in YEAR_TIME}
    spec = {y: fcc.specifier_rows(y, _spec_path(y)) for y in YEAR_TIME}

    n_split = 0
    for p in pairs:
        ca = clus[p["year_a"]][p["row_a"]]
        cb = clus[p["year_b"]][p["row_b"]]
        n_split += (ca != cb)
        _, la, lo = spec[p["year_a"]][p["row_a"]]
        near, near_m = _nearest_cluster(la, lo)
        p.update({"cluster": near if deme_choice == "nearest" else cb,
                  "cluster_prod_a": ca, "cluster_prod_b": cb,
                  "cluster_nearest": near, "nearest_km": round(near_m / 1000, 2)})

    if verbose and n_masked:
        print(f"  {n_masked} field-pair(s) dropped by the n < {fcc.MIN_SUBPOP_N} mask (7.0).")
    if verbose and n_split:
        print(f"  NOTE: {n_split}/{len(pairs)} field-pairs are assigned to DIFFERENT demes in "
              f"their two years by the production mapping (CLAUDE.md 7.9.2).")
        print(f"  Using deme_choice='{deme_choice}'.")
    return pairs


# ---------------------------------------------------------------------------
# Fc
# ---------------------------------------------------------------------------

def _subsample_nodes(ts, pop_idx, tval, n_diploid, rng):
    """n_diploid whole INDIVIDUALS from one deme at one timepoint -> their sample nodes.

    Individuals, not loose nodes: the empirical unit is a diploid. Same argument and the same
    construction as AnalyzeTreeSeq._subsample_nodes / fst_subsample.py (CLAUDE.md 6.6, 7.5a).
    """
    nodes = ts.samples(population=pop_idx, time=tval)
    by_ind = {}
    for nd in nodes:
        by_ind.setdefault(ts.node(nd).individual, []).append(nd)
    inds = np.array(sorted(by_ind), dtype=np.int64)
    if len(inds) > n_diploid:
        inds = rng.choice(inds, size=n_diploid, replace=False)
    return np.array(sorted(nd for i in inds for nd in by_ind[i]), dtype=np.int64), len(inds)


def waples_fc(x, y, sa, sb, t, valid=None):
    """Delegates to fc_common so the empirical side computes the identical thing."""
    return fcc.waples_fc(x, y, sa, sb, t, valid=valid)


def freqs_for_groups(ts, groups, max_bytes=None):
    """Derived-allele frequency per site for each node group, ONE genotype_matrix call.

    Batched for the reason 7.5.2 measured the hard way: genotype_matrix's cost tracks SITES AND
    TREES, not sample count, so calling it once per group pays the whole extraction per group.
    One call over the union costs what one call over four nodes costs.
    """
    # DEDUPLICATE. Groups routinely share nodes: two field-pairs can resolve to the same nearest
    # cluster, and a field in all three years contributes the same (cluster, time) to two pairs.
    # genotype_matrix raises TSK_ERR_DUPLICATE_SAMPLE on a repeated node, so index into a unique
    # list rather than concatenating.
    uniq = sorted({int(n) for g in groups for n in g})

    # CHUNK if the matrix would be too large to hold. sites x nodes as int8: at POPMULT=5000 with
    # --big-n 50 that is ~243k x 3.3k = 800 MB, which was enough to get this killed on a 15 GB box
    # with other work running. Chunking costs one extra genotype_matrix traversal per chunk
    # (extraction tracks SITES AND TREES, not sample count -- CLAUDE.md 7.5.2), so it is a real
    # time-for-memory trade, not free. One chunk reproduces the unchunked result exactly.
    if max_bytes and ts.num_sites * len(uniq) > max_bytes:
        n_chunk = max(1, int(max_bytes // max(1, ts.num_sites)))
        chunks = [uniq[i:i + n_chunk] for i in range(0, len(uniq), n_chunk)]
    else:
        chunks = [uniq]

    if len(chunks) > 1:
        out = [np.zeros(ts.num_sites) for _ in groups]
        multi = [np.zeros(ts.num_sites, dtype=bool) for _ in groups]
        counts = [0] * len(groups)
        for ch in chunks:
            idx = {n: i for i, n in enumerate(ch)}
            Gc = ts.genotype_matrix(samples=np.array(ch, dtype=np.int64))
            for k, g in enumerate(groups):
                cols = [idx[int(n)] for n in g if int(n) in idx]
                if cols:
                    sub = Gc[:, cols]
                    out[k] += (sub == 1).sum(axis=1, dtype=np.float64)
                    multi[k] |= (sub > 1).any(axis=1)
                    counts[k] += len(cols)
            del Gc
        return [(o / c, m) for o, c, m in zip(out, counts, multi)], ts.num_sites

    pos = {n: i for i, n in enumerate(uniq)}
    G = ts.genotype_matrix(samples=np.array(uniq, dtype=np.int64))
    # Take the mean per group WITHOUT materialising a float copy of the whole matrix: at
    # POPMULT=5000 that is 256k sites x ~3.3k nodes, 845 MB as int8 and 3.4 GB as float32, and
    # this runs on the same box as a ~7.6 GB recapitation. The per-group temp is ~1/30 the size.
    # SLiM stacks mutation ids, so "derived" is any non-zero, not == 1.
    # Frequency of ALLELE 1 specifically, plus a flag for any allele index above 1. NOT
    # `G > 0`: that folds every derived allele together, which is what the synthetic round-trip
    # caught disagreeing with the empirical side (fc_common.waples_fc). 7.1% of sites in a real
    # simulated tree carry more than two alleles, so this is not a corner case.
    out = []
    for g in groups:
        cols = [pos[int(n)] for n in g]
        sub = G[:, cols]
        out.append(((sub == 1).mean(axis=1, dtype=np.float64), (sub > 1).any(axis=1)))
    return out, G.shape[0]


def mean_hudson_fst(ts, pairs, arm_n, rng):
    """Mean pairwise Hudson F_st over the matched-field demes at the 2023 timepoint.

    Reported alongside Fc so the two statistics can be read off the SAME tree. This is what makes
    the iso-Nm test possible: F_st tracks 1/(1+4Nm) (CLAUDE.md 6.2), so two points with equal Nm
    should agree on F_st -- and whether Fc still separates them is the whole question.

    Hudson, computed from pi and d_xy, NOT ts.Fst -- which is Nei and ~2x low (CLAUDE.md 6.7,
    invariant 9). Batched indexes= traversal, per CLAUDE.md 2.
    """
    demes = sorted({p["cluster"] for p in pairs})
    sets, sizes = [], []
    for c in demes:
        nodes, n = _subsample_nodes(ts, c, 0, arm_n or 8, rng)
        sets.append(nodes)
        sizes.append(n)
    idx = [(i, j) for i in range(len(sets)) for j in range(i + 1, len(sets))]
    div = ts.diversity(sets, mode="site", span_normalise=True)
    dxy = ts.divergence(sets, indexes=idx, mode="site", span_normalise=True)
    vals = []
    for k, (i, j) in enumerate(idx):
        if dxy[k] > 0:
            vals.append(1.0 - 0.5 * (div[i] + div[j]) / dxy[k])
    return float(np.mean(vals)), len(demes), float(np.mean(sizes))


def run_point(ts, pairs, popmult, arm_n, reps, seed, max_bytes=None,
              simplify_first=False):
    """Fc for every matched pair, `reps` independent subsample draws."""
    per_rep = []
    for rep in range(reps):
        rng = np.random.default_rng(seed * 1000 + rep)
        groups, meta = [], []
        for p in pairs:
            na = p["n_a"] if arm_n is None else arm_n
            nb = p["n_b"] if arm_n is None else arm_n
            ga, sa = _subsample_nodes(ts, p["cluster"], YEAR_TIME[p["year_a"]], na, rng)
            gb, sb = _subsample_nodes(ts, p["cluster"], YEAR_TIME[p["year_b"]], nb, rng)
            groups += [ga, gb]
            meta.append((p, sa, sb))
        # SIMPLIFY to just this rep's nodes before extracting. Measured on the POPMULT=5000
        # kernel tree: 70,078 samples / 5.55M edges / 242,685 sites -> 3,300 samples / 1.02M
        # edges / 57,816 sites, so the genotype matrix drops ~4x AND the traversal gets much
        # faster. Exact for this statistic: allele frequencies among the retained samples are
        # untouched by simplify, and the sites it drops (filter_sites=True) are monomorphic in
        # the retained set, which the MAF filter would have removed anyway.
        # filter_populations=False for the usual reason (CLAUDE.md 2) -- nothing downstream of
        # here queries by population, but the habit is what keeps that bug dead.
        work, gwork = ts, groups
        if simplify_first:
            allnodes = np.unique(np.concatenate(groups))
            work, nmap = ts.simplify(samples=allnodes, filter_populations=False, map_nodes=True)
            gwork = [nmap[g] for g in groups]
            if any((g < 0).any() for g in gwork):
                raise RuntimeError("simplify dropped a node that a group needs")
        fr, n_sites = freqs_for_groups(work, gwork, max_bytes)

        rows = []
        for i, (p, sa, sb) in enumerate(meta):
            x, mx = fr[2 * i]
            y, my = fr[2 * i + 1]
            # biallelic PER PAIR (invariant 7), the same rule the empirical side applies.
            # fcc.waples_fc applies fcc.maf_mask itself -- do NOT pre-filter here as well, or the
            # ascertainment gets applied twice and the two sides stop matching.
            fc, ped, ne, nloc = waples_fc(x, y, sa, sb, p["t"], valid=~mx & ~my)
            rows.append({"site_a": p["site_a"], "site_b": p["site_b"], "t": p["t"],
                         "sa": sa, "sb": sb, "fc": fc, "pedestal": ped,
                         "ne_hat": ne, "n_loci": nloc})
        per_rep.append(rows)

    # Pool the way the ABC would: one number per generation-gap, then the mean over gaps.
    out = {}
    for t in sorted({p["t"] for p in pairs}):
        vals = [np.mean([r["fc"] for r in rows if r["t"] == t]) for rows in per_rep]
        # Mean pedestal, so the drift term can be read off. Do NOT average per-pair ne_hat:
        # at n=7 the per-pair drift term is a small difference of two noisy numbers, so 1/drift
        # has near-zero denominators and its mean is meaningless. Pool Fc first, then invert.
        ped = float(np.mean([r["pedestal"] for r in per_rep[0] if r["t"] == t]))
        fcm = float(np.mean(vals))
        out[f"t{t}"] = {"fc_mean": fcm, "fc_sd": float(np.std(vals, ddof=1)) if reps > 1 else 0.0,
                        "pedestal_mean": ped,
                        "ne_hat": (t / (2 * (fcm - ped))) if fcm > ped else float("inf"),
                        "n_pairs": sum(1 for p in pairs if p["t"] == t)}
    allv = [np.mean([r["fc"] for r in rows]) for rows in per_rep]
    out["all"] = {"fc_mean": float(np.mean(allv)),
                  "fc_sd": float(np.std(allv, ddof=1)) if reps > 1 else 0.0}
    out["per_pair"] = per_rep[0]
    return out


# ---------------------------------------------------------------------------

def prepare(args):
    if args.load_ts:
        ts = tskit.load(args.load_ts)
        if ts.num_sites == 0:
            raise SystemExit(f"{args.load_ts} carries no sites -- that is a RAW SLiM tree. "
                             f"Pass it as --raw-ts so it gets recapitated and mutated.")
        return ts, 0.0
    if not args.raw_ts:
        raise SystemExit("need --load-ts (mutated tree) or --raw-ts (raw SLiM output)")

    import msprime
    import pyslim
    t0 = time.time()
    ts = tskit.load(args.raw_ts)
    ts = _ru.recapitate(ts, recombination_rate=args.recomb, ancestral_Ne=args.anc_ne)
    keep = []
    for y, tv in YEAR_TIME.items():
        for c in _cluster_of_specifier_row(y):
            keep.extend(ts.samples(population=c, time=tv))
    # filter_populations=False is MANDATORY (CLAUDE.md 2) -- the default renumbers populations and
    # every ts.samples(population=cluster_row) below would then read the wrong deme.
    ts = ts.simplify(samples=keep, filter_populations=False)
    ts = msprime.sim_mutations(ts, rate=args.mu,
                               model=msprime.SLiMMutationModel(
                                   type=0, next_id=pyslim.next_slim_mutation_id(ts)),
                               keep=True)
    if args.save_ts:
        ts.dump(args.save_ts)
    return ts, time.time() - t0


def summarize():
    recs = [json.loads(l) for l in open(OUT_JSONL, encoding="utf-8") if l.strip()]
    if not recs:
        raise SystemExit(f"no records in {OUT_JSONL}")
    print("=" * 78)
    print("TEMPORAL Fc -- does it move with POPMULT?  (CLAUDE.md 7.9)")
    print("=" * 78)
    for arm in ("n_real", "n_big"):
        rows = [r for r in recs if r["arm"] == arm]
        if not rows:
            continue
        rows.sort(key=lambda r: r["popmult"])
        print(f"\n  ARM: {arm}" + ("   (real n_i, 5-8 diploids -- what is measurable)"
                                   if arm == "n_real" else
                                   f"   (n={rows[0]['big_n']} diploids -- pedestal suppressed)"))
        print(f"    {'POPMULT':>8} {'tm':>6} {'demeN':>6} {'Nm':>6} {'Fc(t=16)':>10} "
              f"{'+-sd':>8} {'Fc-ped':>8} {'Fc(t=8)':>9} {'Fst':>9}")
        for r in rows:
            deme = 3.33 * r["popmult"] / (33 * r["num_clusters"])
            tm = r.get("total_migration", float("nan"))
            a = r["res"].get("t16", {})
            b = r["res"].get("t8", {})
            ped = a.get("pedestal_mean", float("nan"))
            print(f"    {r['popmult']:8.0f} {tm:6.3f} {deme:6.0f} {deme * tm:6.1f} "
                  f"{a.get('fc_mean', float('nan')):10.5f} {a.get('fc_sd', 0):8.5f} "
                  f"{a.get('fc_mean', float('nan')) - ped:8.5f} "
                  f"{b.get('fc_mean', float('nan')):9.5f} "
                  f"{r['res'].get('fst_hudson', float('nan')):9.5f}")
        if len(rows) >= 2:
            lo, hi = rows[0], rows[-1]
            for key, lab in (("t16", "t=16"), ("t8", "t=8 ")):
                if key not in lo["res"] or key not in hi["res"]:
                    continue
                d = lo["res"][key]["fc_mean"] - hi["res"][key]["fc_mean"]
                sd = math.hypot(lo["res"][key]["fc_sd"], hi["res"][key]["fc_sd"])
                print(f"    {lab}: POPMULT {lo['popmult']:.0f} -> {hi['popmult']:.0f} moves Fc by "
                      f"{d:+.5f} ({d / lo['res'][key]['fc_mean'] * 100:+.1f}%), "
                      f"subsample sd {sd:.5f}" + (f"  -> {abs(d) / sd:.1f} sd" if sd > 0 else ""))
    print("\n  NOTE: sd here is the SUBSAMPLE-draw spread on ONE tree only. It is NOT the 7.3-style")
    print("  noise floor, which needs SLiM, recapitation and the mutation overlay re-rolled.")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--load-ts", help="an already recapitated+simplified+mutated tree")
    p.add_argument("--raw-ts", help="raw SLiM output; recapitate+simplify+mutate first")
    p.add_argument("--save-ts", help="dump the mutated tree so re-runs cost seconds")
    p.add_argument("--popmult", type=float, required=False)
    p.add_argument("--num-clusters", type=int, default=1, help="raw draw; x33 demes")
    p.add_argument("--total-migration", type=float, default=0.05,
                   help="LABEL ONLY -- records the tm the tree was built at, for the iso-Nm test")
    p.add_argument("--mu", type=float, default=ABC.DEFAULT_MUTATION_RATE)
    p.add_argument("--recomb", type=float, default=ABC.DEFAULT_RECOMBINATION_RATE)
    p.add_argument("--anc-ne", type=int, default=_sc.ANCESTRAL_NE)
    p.add_argument("--big-n", type=int, default=50,
                   help="diploids per timepoint in the pedestal-suppressed arm")
    p.add_argument("--reps", type=int, default=5, help="independent subsample draws")
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--deme-choice", choices=("nearest", "later"), default="nearest",
                   help="which simulated deme stands for a resampled field (see matched_pairs)")
    p.add_argument("--max-bytes", type=int, default=0,
                   help="cap the genotype matrix at this many int8 entries; 0 = one call. "
                        "Chunking costs one extra traversal per chunk but bounds peak memory.")
    p.add_argument("--simplify-first", action="store_true",
                   help="simplify to each rep's own nodes before extracting genotypes. Exact, "
                        "and much smaller/faster on a big tree -- see run_point.")
    p.add_argument("--summarize", action="store_true")
    args = p.parse_args()

    if args.summarize:
        return summarize()
    if args.popmult is None:
        raise SystemExit("--popmult is required (it labels the record)")

    pairs = matched_pairs(args.deme_choice)
    print(f"fc_common spec {fcc.spec_hash()} -- the empirical side must print the same hash")
    print(f"{len(pairs)} resampled field-pairs: " +
          ", ".join(f"t={t}:{sum(1 for q in pairs if q['t'] == t)}"
                    for t in sorted({q['t'] for q in pairs})))

    ts, setup_s = prepare(args)
    print(f"tree: {ts.num_samples} samples, {ts.num_sites} sites, {ts.num_populations} pops"
          f"{f' (setup {setup_s:.0f} s)' if setup_s else ''}")

    for arm, arm_n in (("n_real", None), ("n_big", args.big_n)):
        t0 = time.time()
        res = run_point(ts, pairs, args.popmult, arm_n, args.reps, args.seed,
                        args.max_bytes, args.simplify_first)
        fst, n_demes, mean_n = mean_hudson_fst(ts, pairs, arm_n,
                                               np.random.default_rng(args.seed))
        res["fst_hudson"] = fst
        res["fst_n_demes"] = n_demes
        res["fst_mean_n"] = mean_n
        rec = {"arm": arm, "big_n": args.big_n if arm_n else None, "popmult": args.popmult,
               "deme_choice": args.deme_choice, "total_migration": args.total_migration,
               "num_clusters": args.num_clusters, "mu": args.mu, "recomb": args.recomb,
               "anc_ne": args.anc_ne, "seed": args.seed, "reps": args.reps,
               "min_maf": MIN_MAF, "fc_spec": fcc.spec_hash(), "ts": args.load_ts or args.raw_ts,
               "when": time.strftime("%Y-%m-%dT%H:%M:%S"), "wall_s": time.time() - t0,
               "res": res}
        OUT_JSONL.parent.mkdir(parents=True, exist_ok=True)
        with open(OUT_JSONL, "a", encoding="utf-8") as f:
            f.write(json.dumps(rec) + "\n")
        deme = 3.33 * args.popmult / (33 * args.num_clusters)
        print(f"  {arm:7s} Fc(t=16)={res.get('t16', {}).get('fc_mean', float('nan')):.5f} "
              f"Fc(t=8)={res.get('t8', {}).get('fc_mean', float('nan')):.5f}  "
              f"[t/(2N) at deme N={deme:.0f} would be {16 / (2 * deme):.5f}]  "
              f"Fst={fst:.5f}  {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
