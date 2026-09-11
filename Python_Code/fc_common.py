"""Shared temporal-F_c spec -- imported by BOTH sides (CLAUDE.md 7.9.4, 7.9.6).

The two sides run on different machines, so this file must be **COPIED** next to
`ToUseOnBeagles/` exactly as `ld_common.py` is. `spec_hash()` is printed by both and a mismatch
means the two F_c values are not comparable. Never reimplement any of this on one side only --
that is the mistake `ld_common` exists to prevent, and F_c is MORE fragile than LD was, because
its MAF filter is an *ascertainment* rather than a cosmetic cut (see MIN_MAF).

WHAT F_c IS. Waples (1989) standardised variance in allele frequency between two samples of the
same field taken t generations apart:

    F_c = mean_loci  (x - y)^2 / ( (x + y)/2 - x*y )
    E[F_c] ~ t/(2*Ne) + 1/(2*S_a) + 1/(2*S_b)      for S diploids sampled at each timepoint

**The pedestal is NOT subtracted.** In ABC it is reproduced identically on both sides (invariant 1)
and cancels in |F_c,sim - F_c,obs|. Subtracting it would be correct for a point estimate of Ne and
is wrong here -- and at n=7 it is not even reliable (see MIN_MAF).

WHY IT IS WORTH THE TROUBLE. Everything else in this project reads the same 0.5-7% sliver of the
genealogy and reads it as an `Nm` quantity (CLAUDE.md 7.9.1). F_c measures a *change* between two
timepoints, so the shared ancestral phase cancels exactly and it is 100% forward-phase. Measured
(7.9.6): swinging the dispersal kernel 12x moves F_st by +170 to +201% and moves F_c by less than
1.5 sd of its own noise. Across the kernel prior, F_st's implied N swings 5.3-5.8x and F_c's
1.07-1.22x.
"""

import hashlib
import math

import numpy as np

# ---------------------------------------------------------------------------
# The spec. EVERY constant here changes the statistic. Changing one invalidates
# any empirical run made under the old value -- bump nothing casually.
# ---------------------------------------------------------------------------

# Pooled-sample minor-allele floor, applied to (x+y)/2 for the PAIR being scored.
#
# THIS IS AN ASCERTAINMENT, NOT A COSMETIC FILTER, and it is the single most important thing to
# match. At n=7 the sample holds 14 haplotypes, so observable frequencies are multiples of 1/14
# and conditioning on (x+y)/2 -- the very quantity forming F_c's denominator -- biases the result.
# Measured consequence (7.9.4): the `F_c - pedestal` decomposition is fine at n=50 (drift ratio
# 2.49 against a deme-size ratio of 2.50) and COLLAPSES at n=7, going negative at POPMULT=5000.
# So: raw F_c is the statistic; the decomposition is a diagnostic for large-n arms only, and
# **ne_hat must never be read at real n**.
#
# Common variants also happen to be the band Beagle imputation handles best, which matters because
# 7.8.5 records that the empirical files are hard calls with no genotype probabilities, no DR2 and
# no missing-data marker -- nothing records what was imputed.
MIN_MAF = 0.05

# Generations before the end of the forward run at which each year's sample sits. 2 generations
# per year (Cohen et al. 2022; CLAUDE.md 11), and SLiM Remembers at ticks 308/316/324.
YEAR_TIME = {"2015": 16, "2019": 8, "2023": 0}

# Year pairs to match on, in the order they are emitted.
YEAR_PAIRS = (("2015", "2019"), ("2019", "2023"), ("2015", "2023"))

# Fields with fewer than this many sequenced diploids are dropped, matching the cut the FITTED
# statistics already use (CLAUDE.md 7.0): at n<=3 a site scatters around a noise floor rather than
# measuring anything. It also removes the `Arlington2015` typo duplicate (n=2), which would
# otherwise double-count that field.
MIN_SUBPOP_N = 4

# Coordinate rounding used to decide that two years sampled the SAME field.
COORD_DECIMALS = 5

# Individual samples removed from F_c on BOTH sides. {sample_id: (site, reason)}. CLAUDE.md 7.2.2F.
#
# WHY. A sample holding relatives, or one whose genotype calls failed, carries fewer effectively
# independent haplotypes than its diploid count claims. That inflates the sampling variance of a
# site frequency, and F_c with it. The simulation draws unrelated individuals with clean calls, so
# it reproduces neither. Measured 2026-09-11: among the 15 matched pairs with n >= 7, the five
# containing such samples were exactly the five largest F_c excesses (exact p = 1/3003).
#
# THE RULE, re-derived from out/kinship_out by diagnostics/kinship_correct.py, which fails if it
# disagrees with this list:
#   failed   -- per-individual het-miscall rate e > 0.4 (observed het 0.004-0.134 against ~0.27).
#               Dropped outright. The 0.4 cut was chosen with the F_c values in view; the
#               separation holds for any cut in 0.35-0.45.
#   relative -- corrected KING kinship above the FIRST-degree cut (0.177) with another sample of
#               the SAME field. ONE member per pair is dropped -- the one whose e is further from
#               its year's median, i.e. the less typical genotype quality -- and the other is kept.
# Only MATCHED fields are listed: Mortensen9-2015 (a whole family) and Alsum25-2015's pairs sit in
# fields no other year resampled, so F_c never sees them.
#
# WHY FIRST-DEGREE AND NOT SECOND. Relatives in a random sample of a small deme are part of the
# drift signal (Waples & Anderson 2017), so a relative may only be dropped here if the simulation
# could not have produced it. Measured with diagnostics/fc_miscall.py from SLiM's recorded parents:
# random simulated field samples contain a half-sib pair 19% of the time at POPMULT 2000 and 8% at
# 5000, and NEVER a full-sib pair (0 of 18,432 pairs). Those half-sibs carry ~20% of the 2000->5000
# F_c signal, so purging the empirical side's half-sibs (S22 in Alsum18-2023, S306 in
# OkrayGrosheks40-2019) would bias F_c low -- they are KEPT (Sohan's call, 2026-09-11). The
# full-sib pair has no simulated counterpart, so one member goes.
#
# Resulting n: Alsum59-2019 7->6, Alsum18-2023 7->6, H15-2023 7->6, H41-2023 7->4 -- still
# >= MIN_SUBPOP_N, so all 19 field pairs survive.
#
# PART OF THE SPEC HASH: editing this list changes the statistic, and the empirical target must be
# re-run under the new hash.
EXCLUDED_SAMPLES = {
    "S32":  ("Alsum59-2019",         "failed: e=0.456, observed het 0.134"),
    "S27":  ("Alsum18-2023",         "relative: 1st-degree with S26 (kept); S27 e=0.058"),
    "S212": ("H15-2023",             "failed: e=0.631, observed het 0.084"),
    "S221": ("H41-2023",             "failed: e=0.954, observed het 0.004"),
    "S222": ("H41-2023",             "failed: e=0.816, observed het 0.040"),
    "S223": ("H41-2023",             "failed: e=0.558, observed het 0.106"),
}


def spec_hash():
    """Short hash of everything that changes the statistic. Both sides print it; a mismatch
    means the curves are not comparable and the empirical run must be redone."""
    payload = repr((MIN_MAF, sorted(YEAR_TIME.items()), YEAR_PAIRS, MIN_SUBPOP_N,
                    COORD_DECIMALS, "biallelic-per-pair-allele-index<=1",
                    tuple(sorted(EXCLUDED_SAMPLES))))
    return hashlib.sha256(payload.encode()).hexdigest()[:12]


# ---------------------------------------------------------------------------
# The resampled fields
# ---------------------------------------------------------------------------

def specifier_rows(year, path):
    """[(site name, lat, lon)] in specifier-matrix ROW ORDER (CLAUDE.md 4, invariant 5)."""
    rows = []
    with open(path, encoding="utf-8") as f:
        for line in f:
            if line.strip():
                p = line.split(",")
                rows.append((p[0].strip(),
                             round(float(p[1]), COORD_DECIMALS),
                             round(float(p[2]), COORD_DECIMALS)))
    return rows


def read_popfile(path):
    """[(sample_id, site)] in file order, from 'sampleID<TAB>site' lines. NO exclusions applied."""
    rows = []
    with open(path, encoding="utf-8") as f:
        for line in f:
            p = line.split()
            if len(p) >= 2:
                rows.append((p[0].strip(), p[1].strip()))
    return rows


def popfile_members(path):
    """{site: [sample ids]} in file order, with EXCLUDED_SAMPLES removed.

    The empirical side takes each field's VCF columns from here and the simulated side takes each
    field's sample size from here (via popfile_counts), so one exclusion shrinks n identically on
    both. Never re-parse a popfile on one side only.
    """
    out = {}
    for sid, site in read_popfile(path):
        if sid not in EXCLUDED_SAMPLES:
            out.setdefault(site, []).append(sid)
    return out


def popfile_counts(path):
    """{site name: number of sequenced diploids AFTER EXCLUDED_SAMPLES}."""
    return {site: len(ids) for site, ids in popfile_members(path).items()}


def check_exclusions(popfile_paths):
    """Raise unless every EXCLUDED_SAMPLES id appears in a popfile under the site it is listed with.

    A typo in the list would otherwise exclude nothing -- or a different beetle -- in silence, and
    the spec hash would still change, so nothing downstream would look wrong.
    """
    where = {}
    for p in popfile_paths.values():
        for sid, site in read_popfile(p):
            where[sid] = site
    for sid, (site, _) in EXCLUDED_SAMPLES.items():
        if sid not in where:
            raise ValueError(f"EXCLUDED_SAMPLES lists {sid!r} ({site}) but no popfile contains it")
        if where[sid] != site:
            raise ValueError(f"EXCLUDED_SAMPLES lists {sid!r} under {site}, but the popfile puts it "
                             f"in {where[sid]} -- wrong id or wrong site")


def matched_field_pairs(specifier_paths, popfile_paths):
    """Fields sampled in two different years, matched on IDENTICAL coordinates.

    COORDINATES, NOT NAMES. The same physical field is `BrilowskiHome-2015` and
    `BrikalskiPats-2023`, and `GarrisonNE-2015` is `OkrayGarrisonNE-2019`; name matching loses
    both silently. Measured: 19 pairs survive the MIN_SUBPOP_N cut -- 9 at t=16 generations
    (2015 vs 2023), 10 at t=8.

    specifier_paths / popfile_paths: {year: path}. Returns dicts with the site names, the real
    per-timepoint sample sizes, and the generation gap.
    """
    check_exclusions(popfile_paths)
    spec = {y: specifier_rows(y, p) for y, p in specifier_paths.items()}
    cnt = {y: popfile_counts(p) for y, p in popfile_paths.items()}

    pairs, n_masked = [], 0
    for ya, yb in YEAR_PAIRS:
        idx_b = {(la, lo): i for i, (_, la, lo) in enumerate(spec[yb])}
        for i, (name_a, la, lo) in enumerate(spec[ya]):
            j = idx_b.get((la, lo))
            if j is None:
                continue
            name_b = spec[yb][j][0]
            na, nb = cnt[ya].get(name_a), cnt[yb].get(name_b)
            if na is None or nb is None:
                raise ValueError(f"{name_a}/{name_b}: absent from a popfile")
            if min(na, nb) < MIN_SUBPOP_N:
                n_masked += 1
                continue
            pairs.append({"year_a": ya, "year_b": yb, "site_a": name_a, "site_b": name_b,
                          "row_a": i, "row_b": j, "n_a": na, "n_b": nb,
                          "t": YEAR_TIME[ya] - YEAR_TIME[yb]})
    return pairs, n_masked


# ---------------------------------------------------------------------------
# The statistic
# ---------------------------------------------------------------------------

def maf_mask(x, y):
    """Sites kept for a pair: pooled frequency inside [MIN_MAF, 1-MIN_MAF]. See MIN_MAF."""
    pooled = 0.5 * (np.asarray(x) + np.asarray(y))
    return (pooled >= MIN_MAF) & (pooled <= 1.0 - MIN_MAF)


def fc_terms(x, y):
    """Per-locus F_c contributions for already-masked frequency vectors.

    Returns (sum of per-locus values, count). Kept separate from the mean so chromosomes can be
    POOLED rather than averaged: genome-wide F_c is (sum over all loci)/(total loci), never the
    mean of per-chromosome F_c values (CLAUDE.md 5.2, invariant 4).
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    den = 0.5 * (x + y) - x * y
    ok = den > 0
    return float(np.sum((x[ok] - y[ok]) ** 2 / den[ok])), int(ok.sum())


def waples_fc(x, y, sa, sb, t, valid=None):
    """F_c for one field pair, plus the plan-II Ne it implies.

    x, y: frequency of allele 1 at each timepoint. **NOT "frequency of anything non-ancestral".**
    valid: sites this pair may use at all, before the MAF cut. Both sides MUST pass a
      BIALLELIC mask here (invariant 7): a site is usable only if no allele index above 1 appears
      in either timepoint's samples.

    WHY THIS ARGUMENT EXISTS. A synthetic round-trip of the two implementations disagreed by up to
    2.6e-3 in F_c -- small enough to look like rounding and far too large to be. The cause was
    multi-allelic sites, handled differently on each side: the simulated code folded every derived
    allele together as `G > 0`, while the empirical code counted only allele 1 and silently shrank
    the denominator. **7.1% of sites in a real simulated tree carry more than two alleles**
    (recurrent mutation), so this was not a corner case. Dropping them is unbiased for F_c --
    recurrent-mutation sites differ in mutation rate, not in drift -- and, far more importantly,
    it is now the SAME rule on both sides.

    Sites are judged biallelic PER PAIR, on that pair's own samples, exactly as the LD code judges
    it per deme (CLAUDE.md 7.5.2): a site can be biallelic in one field and not another.

    sa, sb are DIPLOID sample sizes. ne_hat is a readability check only -- do NOT read it at the
    real n (see MIN_MAF). The pedestal is returned, never subtracted from the statistic.
    """
    keep = maf_mask(x, y)
    if valid is not None:
        keep = keep & np.asarray(valid, dtype=bool)
    total, n_loci = fc_terms(np.asarray(x)[keep], np.asarray(y)[keep])
    fc = total / n_loci if n_loci else float("nan")
    pedestal = 1.0 / (2.0 * sa) + 1.0 / (2.0 * sb)
    drift = fc - pedestal
    ne = (t / (2.0 * drift)) if drift > 0 else float("inf")
    return fc, pedestal, ne, n_loci


# ---------------------------------------------------------------------------
# Reproducing the empirical genotype miscalls -- SIMULATED SIDE ONLY (CLAUDE.md 7.2.2F)
# ---------------------------------------------------------------------------
#
# Every empirical sample carries a per-individual rate e at which true heterozygotes were called
# homozygous for one of their alleles (median ~0.24; diagnostics/kinship_correct.py). Per
# individual that acts like an inbreeding coefficient, inflating the sampling variance of a site
# frequency. MEASURED on simulated trees (diagnostics/fc_miscall.py): it raises pooled F_c by
# +0.034 at POPMULT 2000 and +0.036 at 5000, i.e. 5.8x the whole 2000->5000 F_c signal. So the
# simulated side puts the same miscalls into its genotypes before taking frequencies; that keeps
# ~79% of the signal.
#
# NOT IN spec_hash(), deliberately. The empirical side computes nothing differently -- its data
# already carry the miscalls -- so hashing this would force a multi-hour empirical re-run for a
# change that cannot move the empirical number. The rates file is the simulated side's calibration
# input, like the mutation rate.

MISCALL_RATES_FILE = "miscall_rates.csv"    # in data/empiricalStats/; kinship_correct.py --write-rates


def read_miscall_rates(path):
    """{sample_id: e} from the rates file."""
    import csv
    with open(path, newline="", encoding="utf-8") as f:
        out = {r["sample_id"]: float(r["e"]) for r in csv.DictReader(f)}
    if not out:
        raise ValueError(f"{path} holds no miscall rates")
    return out


def field_miscall_rates(rates, popfile_path, site):
    """e for one field's RETAINED samples (EXCLUDED_SAMPLES removed), clipped to [0, 1].

    Clipped below at 0 because a sample with MORE heterozygotes than HWE allows cannot be imitated
    by miscalling; clipped above at 1 because a rate cannot exceed it. Raises on a missing sample
    rather than inventing a rate for it.
    """
    ids = popfile_members(popfile_path).get(site, [])
    missing = [s for s in ids if s not in rates]
    if missing:
        raise ValueError(f"{site}: no miscall rate for {missing[:3]} -- regenerate {MISCALL_RATES_FILE}")
    return np.clip([rates[s] for s in ids], 0.0, 1.0)


def apply_miscall(G, e, rng):
    """Imitate the empirical miscalls on simulated genotypes.

    G: (sites, n, 2) allele indices for n diploid individuals. e: n rates, one per individual --
    the caller assigns the field's real rates to simulated individuals in random order.
    Independently at every heterozygous site, with probability e[k], individual k becomes
    homozygous for one of its two alleles chosen at random. Returns a new array.
    """
    g0, g1 = G[:, :, 0], G[:, :, 1]
    flip = (g0 != g1) & (rng.random(g0.shape) < np.asarray(e)[None, :])
    pick = np.where(rng.integers(0, 2, size=g0.shape) == 0, g0, g1)
    out = G.copy()
    out[:, :, 0] = np.where(flip, pick, g0)
    out[:, :, 1] = np.where(flip, pick, g1)
    return out


def pool_by_gap(per_pair):
    """Pool per-pair F_c into one number per generation gap, then the mean over gaps.

    Mirrors how calculate_losses normalises out per-year entry counts (CLAUDE.md 7), so a gap with
    10 field pairs does not outvote one with 9. `per_pair` is an iterable of dicts carrying 't'
    and 'fc'.
    """
    out = {}
    for t in sorted({r["t"] for r in per_pair}):
        vals = [r["fc"] for r in per_pair if r["t"] == t and not math.isnan(r["fc"])]
        out[f"t{t}"] = {"fc_mean": float(np.mean(vals)), "n_pairs": len(vals)}
    allv = [r["fc"] for r in per_pair if not math.isnan(r["fc"])]
    out["all"] = {"fc_mean": float(np.mean(allv)), "n_pairs": len(allv)}
    return out
