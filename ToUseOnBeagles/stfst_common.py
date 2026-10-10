"""Shared spatio-temporal F_st spec -- imported by BOTH sides (CLAUDE.md 7.9.14B).

The two sides run on different machines, so this file must be **COPIED** next to
`ToUseOnBeagles/` exactly as `fc_common.py` is (it imports `sample_exclusions`, which is copied
there already). `spec_hash()` is printed by both; a mismatch means the two sides are not computing
the same statistic. Never reimplement any of this on one side only.

WHAT IT IS. Pairwise F_st between site i sampled in year A and site j sampled in year B, for EVERY
pair of retained sites across all three collections -- same-year pairs (gap 0) included, so one
file holds both halves of the comparison. Under structure that persists through time, cross-year
F_st(i, j) ~ same-year F_st at the same distance plus a small drift term; a model that drifts too
much shows cross-year F_st climbing above same-year. So the DIAGNOSTIC is the contrast

    median F_st(gap g, distance bin b)  -  median F_st(gap 0, distance bin b)

by distance x gap (`contrast`), on each side, compared between sides. The old temporal F_c is
roughly the distance < SAME_POP_KM, gap > 0 cell of this design.

DESIGN RULES (CLAUDE.md 7.9.14B, from the professor's input):
  * ALL pairs, binned -- not "the farthest few". IBD is null in every year, so far pairs are not
    special, and extremes are noisy.
  * MEDIANS, not means: single-site isolates dominate any mean (7.1-7.5).
  * "Same place" is DISTANCE (<= SAME_POP_KM), never site name (7.9.14C/E3).
  * Sample sizes match the real ones exactly (invariant 1 -- F_st's expectation moves with n), and
    one draw per (year, site) is reused in every pair, as the empirical side reuses its beetles.
  * ESTIMATOR: the simulated side is HUDSON, the empirical side pixy's WC, as for the production
    F_st (invariant 9; at this differentiation the two agree, Nei/ts.Fst is ~half). WC also keeps
    the ~24% het deficit (7.2.2) out of F_ST -- it goes to F_IS -- so no miscalls are imitated.

PER-PAIR FILE FORMAT (both sides write it; `read_pairs` reads it):
    year_a, site_a, year_b, site_b, gap, dist_km, n_a, n_b, fst
  year_a <= year_b (collection order); within a year, specifier row order with row_a < row_b.
  n_a / n_b are the diploids ACTUALLY used. The simulated side may also add deme_a, deme_b.

EMPIRICAL SIDE: ToUseOnBeagles/CalcSpatiotemporalFst.py. pixy computes F_st only between
populations of ONE popfile on ONE VCF -- which is already the situation: each chr{i}_cpb.vcf.gz
comes from one Beagle file holding all three years (imputed together), sample IDs are unique across
years and site labels carry the year. So one combined popfile gives every pair, cross-year and
same-year, over identical SNPs; no VCF merge or SNP reconciliation is needed.

Distances are between the REAL site coordinates (specifier matrices) on both sides -- the deme a
site maps to in the simulation does not move its distance.
"""

import csv
import hashlib
import math
from pathlib import Path

import numpy as np

import sample_exclusions

YEARS = ("2015", "2019", "2023")
# Sampling time of each collection in generations before 2023 (SLiM ticks 308/316/324, 2 per year).
YEAR_TIME = {"2015": 16, "2019": 8, "2023": 0}

# A site needs this many diploids (after sample_exclusions) to enter the summary -- the same cut
# as the fitted statistics (ABCAnalysisNoRedis.MIN_SUBPOP_N, CLAUDE.md 7.0). Pairs are still
# WRITTEN for smaller sites; the cut is applied in `summarise`.
MIN_SUBPOP_N = 4

# "Same local population" across collections (CLAUDE.md 7.9.14E3: resampling sits a median
# 0.5-1.8 km from the earlier site, against 2.9-4.8 km to the nearest same-year neighbour).
SAME_POP_KM = 2.0

# Distance bins, km. The first is SAME_POP_KM. Same-year pairs under 2 km barely exist (4 of 602),
# so the 0-2 km cross-year cells are contrasted against the 2-10 km same-year cell (`contrast`).
# 2-10 is one bin because 2-5 km holds only 8/7/8 same-year pairs per year.
DIST_EDGES_KM = (0.0, 2.0, 10.0, 20.0, 50.0, 100.0, math.inf)

# A (gap, bin) cell with fewer pairs than this reports NaN rather than a median of a handful.
MIN_PAIRS_PER_CELL = 5

R_EARTH_KM = 6371.0088

FIELDS = ["year_a", "site_a", "year_b", "site_b", "gap", "dist_km", "n_a", "n_b", "fst"]


def spec_hash():
    """Short hash of everything that changes the statistic or its summary."""
    payload = repr((YEARS, sorted(YEAR_TIME.items()), MIN_SUBPOP_N, SAME_POP_KM, DIST_EDGES_KM,
                    MIN_PAIRS_PER_CELL, "sim-hudson/emp-pixy-wc", "median",
                    sample_exclusions.exclusion_id()))
    return hashlib.sha256(payload.encode()).hexdigest()[:12]


# ---------------------------------------------------------------------------
# Sites and pairs
# ---------------------------------------------------------------------------

def site_table(genetic_data_dir):
    """Every site of every year, in collection order then SPECIFIER ROW ORDER (CLAUDE.md 4):
    [{"year", "row", "site", "lat", "lon", "n"}], n = diploids after sample_exclusions."""
    gd = Path(genetic_data_dir)
    out = []
    for y in YEARS:
        counts = sample_exclusions.popfile_counts(gd / f"popFile{y}")
        with open(gd / f"specifier_matrix_{y}.csv", encoding="utf-8") as f:
            rows = [line.split(",") for line in f if line.strip()]
        for r, p in enumerate(rows):
            site = p[0].strip()
            if site not in counts:
                raise ValueError(f"{y}: specifier site {site} absent from popFile{y}")
            out.append({"year": y, "row": r, "site": site, "lat": float(p[1]),
                        "lon": float(p[2]), "n": counts[site]})
    return out


def haversine_km(la1, lo1, la2, lo2):
    p = math.pi / 180
    h = (math.sin((la2 - la1) * p / 2) ** 2
         + math.cos(la1 * p) * math.cos(la2 * p) * math.sin((lo2 - lo1) * p / 2) ** 2)
    return 2 * R_EARTH_KM * math.asin(math.sqrt(h))


def all_pairs(sites):
    """[(ia, ib, gap, dist_km)] over indexes into `sites` (from site_table): every unordered pair,
    same-year and cross-year, with ia before ib in site_table order (so year_a <= year_b)."""
    out = []
    for ia in range(len(sites)):
        for ib in range(ia + 1, len(sites)):
            a, b = sites[ia], sites[ib]
            gap = YEAR_TIME[a["year"]] - YEAR_TIME[b["year"]]
            out.append((ia, ib, gap, haversine_km(a["lat"], a["lon"], b["lat"], b["lon"])))
    return out


def write_pairs(path, rows, extra=()):
    """rows: dicts with FIELDS (+ any `extra` columns)."""
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(FIELDS + list(extra))
        for r in rows:
            w.writerow([r["year_a"], r["site_a"], r["year_b"], r["site_b"], r["gap"],
                        f"{r['dist_km']:.4f}", r["n_a"], r["n_b"], f"{r['fst']:.10g}"]
                       + [r[k] for k in extra])


def read_pairs(path):
    with open(path, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    for r in rows:
        r["gap"] = int(r["gap"])
        r["dist_km"] = float(r["dist_km"])
        r["n_a"], r["n_b"] = int(r["n_a"]), int(r["n_b"])
        r["fst"] = float(r["fst"])
    return rows


# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

def dist_bin(d):
    for k in range(len(DIST_EDGES_KM) - 1):
        if DIST_EDGES_KM[k] <= d < DIST_EDGES_KM[k + 1]:
            return k
    raise ValueError(f"distance {d} outside DIST_EDGES_KM")


def bin_label(k):
    lo, hi = DIST_EDGES_KM[k], DIST_EDGES_KM[k + 1]
    return f"{lo:g}-{hi:g}km" if math.isfinite(hi) else f">{lo:g}km"


def retained_sites(genetic_data_dir, exclude_sites=()):
    """{(year, site)} entering the summary: n >= MIN_SUBPOP_N and not in exclude_sites. Pass the
    same whole-site exclusions the fitted statistics use (ABCAnalysisNoRedis.EXCLUDED_SITES)."""
    return {(s["year"], s["site"]) for s in site_table(genetic_data_dir)
            if s["n"] >= MIN_SUBPOP_N and s["site"] not in set(exclude_sites)}


def summarise(rows, keep):
    """{(gap, bin): {"n_pairs", "median"}} over pairs whose BOTH sites are in `keep`.
    median is NaN below MIN_PAIRS_PER_CELL. Gaps are 0, 8, 16."""
    cells = {}
    for r in rows:
        if (r["year_a"], r["site_a"]) in keep and (r["year_b"], r["site_b"]) in keep:
            cells.setdefault((r["gap"], dist_bin(r["dist_km"])), []).append(r["fst"])
    out = {}
    for g in sorted({YEAR_TIME[a] - YEAR_TIME[b] for a in YEARS for b in YEARS
                     if YEAR_TIME[a] >= YEAR_TIME[b]}):
        for k in range(len(DIST_EDGES_KM) - 1):
            v = cells.get((g, k), [])
            out[(g, k)] = {"n_pairs": len(v),
                           "median": float(np.median(v)) if len(v) >= MIN_PAIRS_PER_CELL
                           else math.nan}
    return out


def contrast(summary):
    """{(gap, bin): median cross-year F_st - median same-year F_st} for gap > 0. The same-year
    reference for the 0-2 km bin is the 2-10 km cell (same-year pairs under 2 km barely exist)."""
    out = {}
    for (g, k), c in summary.items():
        if g == 0:
            continue
        ref = summary[(0, max(k, 1))]["median"]
        out[(g, k)] = c["median"] - ref
    return out
