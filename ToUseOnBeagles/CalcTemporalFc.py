#!/usr/bin/env python3
"""Empirical temporal F_c -- the observed target for CLAUDE.md 7.9.4 / 7.9.6.

Runs on the machine holding the Beagle VCFs, NOT in the repo. Output goes to
`data/empiricalStats/averaged_temporalFc.csv` in the repo, alongside the other targets.

WHAT IT COMPUTES. Waples (1989) standardised variance in allele frequency for each field that was
sequenced in TWO different years:

    F_c = mean_loci  (x - y)^2 / ( (x + y)/2 - x*y )

x and y are the ALT-allele frequencies of that field's samples at the two timepoints. 19 field
pairs survive -- 9 at t = 16 generations (2015 vs 2023), 10 at t = 8. See fc_common.

WHY IT MATTERS. Every other statistic in this project reads the same thin, `Nm`-shaped sliver of
the genealogy (CLAUDE.md 7.9.1). F_c measures a *change*, so shared history cancels and it reads
an absolute drift rate. Measured on the simulated side (7.9.6): swinging the dispersal kernel 12x
moves F_st by +170 to +201% and moves F_c by less than 1.5 sd of its own noise.

--------------------------------------------------------------------------------------------
FIVE THINGS THAT WILL SILENTLY RUIN THIS IF YOU CHANGE THEM
--------------------------------------------------------------------------------------------
1. **`fc_common.py` must be COPIED next to this file**, exactly as `ld_common.py` is, and the
   `spec_hash` printed here must match the one the simulated side prints. Different spec => the
   two F_c values are not comparable. This script hard-stops if the copy is missing.
2. **The MAF filter is an ASCERTAINMENT, not a cosmetic cut.** It is applied to the POOLED
   frequency (x+y)/2 of the pair being scored, at exactly `fc_common.MIN_MAF`. At n=7 the sample
   has 14 haplotypes, so conditioning on the quantity that forms F_c's denominator biases it --
   which is fine only because the simulated side does the identical thing (invariant 1).
3. **The sampling pedestal is NOT subtracted here.** 1/(2*S_a)+1/(2*S_b) is reproduced on the
   simulated side by subsampling to the same n_i and cancels in |F_c,sim - F_c,obs|. Subtracting
   it would be right for a point estimate of Ne and wrong for this.
4. **Pool, never average** (CLAUDE.md 5.2, invariant 4). Genome-wide F_c is
   (sum of per-locus values over ALL chromosomes) / (total loci), not the mean of 17
   per-chromosome F_c values. This script accumulates sums and counts.
5. **Fields are matched on COORDINATES, not names** -- the same field is `BrilowskiHome-2015` and
   `BrikalskiPats-2023`. fc_common does the matching; do not reimplement it.
6. **Some samples are EXCLUDED, through fc_common.EXCLUDED_SAMPLES** -- relatives and failed
   genotype calls (CLAUDE.md 7.2.2F). The list is in the spec hash, and n_a/n_b in the output are
   the counts AFTER exclusion, which is what the simulated side subsamples to.

Polarisation does NOT matter here, unusually: F_c is invariant under x -> 1-x, y -> 1-y (the
denominator (x+y)/2 - xy maps to itself), so ALT-vs-derived cannot introduce a bias. That is one
whole class of bug this statistic is immune to.

USAGE
    python CalcTemporalFc.py            # ~ one pass over the 17 VCFs
Outputs
    fc_out/chr{i}_temporalFc.csv        per-chromosome sums and counts -- KEEP THESE. They are the
                                        only route to a between-chromosome jackknife or a re-pool
                                        without re-reading every VCF (7.8.4 needed exactly this
                                        for LD, and would have been impossible without them).
    fc_out/averaged_temporalFc.csv      the target: one row per field pair, plus pooled rows.
"""
import csv
import glob
import os
import re
import sys

import numpy as np
import allel

try:
    import fc_common as fcc
except ImportError:
    raise SystemExit(
        "fc_common.py is not importable. COPY it next to this script from Python_Code/ -- both "
        "sides must run the identical spec, exactly as ld_common.py is copied (CLAUDE.md 9).")

VCF_FILES = sorted(glob.glob("chr*_cpb.vcf.gz"))     # the 17 chromosome files
POPFILES = {"2015": "popFile2015", "2019": "popFile2019", "2023": "popFile2023"}
SPECIFIER = "specifier_matrix_{year}.csv"
OUTDIR = "fc_out"


def sample_columns(popfile_path, site, sample_pos):
    """VCF column indices for one site's samples, with fc_common.EXCLUDED_SAMPLES removed.

    Read through fc_common rather than parsed here, so the samples dropped from this side are by
    construction the ones the simulated side's reduced sample sizes reflect (CLAUDE.md 7.2.2F).
    """
    cols = []
    for sid in fcc.popfile_members(popfile_path).get(site, []):
        if sid not in sample_pos:
            raise SystemExit(f"sample {sid!r} ({site}) is not in the VCF")
        cols.append(sample_pos[sid])
    return cols


def chunk_frequencies(gt, cols_by_group):
    """Per group: (frequency of allele 1, whether any allele index above 1 is present).

    max_allele=3 rather than 1 ON PURPOSE. Counting only alleles 0 and 1 does not *exclude* a
    third allele, it silently shrinks the denominator at that site -- which is exactly the bug the
    synthetic round-trip caught (see fc_common.waples_fc). Here the extra alleles are counted so
    the site can be DROPPED, matching the simulated side and invariant 7.

    Denominator is called alleles, matching CalcGenRel.py. Beagle imputes everything, so AN == 2*n
    in practice; the guard makes a change in the upstream data show up as NaN, not as a wrong
    number.
    """
    ac = gt.count_alleles_subpops(cols_by_group, max_allele=3)
    out = {}
    for g in cols_by_group:
        counts = ac[g]
        alt = counts[:, 1].astype(float)
        an = counts.sum(axis=1).astype(float)
        with np.errstate(invalid="ignore", divide="ignore"):
            freq = np.where(an > 0, alt / an, np.nan)
        multi = counts[:, 2:].sum(axis=1) > 0
        out[g] = (freq, multi)
    return out


def main():
    if not VCF_FILES:
        raise SystemExit("no chr*_cpb.vcf.gz found in the working directory")
    os.makedirs(OUTDIR, exist_ok=True)

    print(f"fc_common spec {fcc.spec_hash()} -- the simulated side must print the same hash")
    print(f"{len(fcc.EXCLUDED_SAMPLES)} samples excluded (relatives / failed calls, CLAUDE.md "
          f"7.2.2F): " + ", ".join(f"{s} ({v[0]})" for s, v in sorted(fcc.EXCLUDED_SAMPLES.items())))

    pairs, n_masked = fcc.matched_field_pairs(
        {y: SPECIFIER.format(year=y) for y in POPFILES},
        {y: p for y, p in POPFILES.items()})
    print(f"{len(pairs)} resampled field-pairs ({n_masked} dropped by the "
          f"n < {fcc.MIN_SUBPOP_N} mask): " +
          ", ".join(f"t={t}:{sum(1 for q in pairs if q['t'] == t)}"
                    for t in sorted({q["t"] for q in pairs})))

    # sample order, read once; assumed identical across chromosomes as in CalcGenRel.py
    _, samples, _, _ = allel.iter_vcf_chunks(VCF_FILES[0], fields=["variants/POS"])
    sample_pos = {s: i for i, s in enumerate(samples)}

    # one column group per (year, site); several pairs can share a group
    groups = {}
    for p in pairs:
        for side in ("a", "b"):
            key = (p[f"year_{side}"], p[f"site_{side}"])
            if key not in groups:
                groups[key] = sample_columns(POPFILES[key[0]], key[1], sample_pos)
            n_expected = p[f"n_{side}"]
            if len(groups[key]) != n_expected:
                raise SystemExit(
                    f"{key[1]}: popfile says {n_expected} diploids but {len(groups[key])} VCF "
                    f"columns matched. The popfile and the VCF disagree about this site.")
    cols_by_group = {f"{y}|{s}": c for (y, s), c in groups.items()}

    totals = {i: [0.0, 0] for i in range(len(pairs))}     # pair index -> [sum, n_loci]
    for vcf in VCF_FILES:
        chrom = re.search(r"chr[^0-9]*([0-9]+)", os.path.basename(vcf)).group(1)
        per_chrom = {i: [0.0, 0] for i in range(len(pairs))}
        _, _, _, chunks = allel.iter_vcf_chunks(vcf, fields=["calldata/GT"])
        for chunk in chunks:
            gt = allel.GenotypeArray(chunk[0]["calldata/GT"])
            freqs = chunk_frequencies(gt, cols_by_group)
            for i, p in enumerate(pairs):
                x, mx = freqs[f"{p['year_a']}|{p['site_a']}"]
                y, my = freqs[f"{p['year_b']}|{p['site_b']}"]
                # biallelic PER PAIR (invariant 7) -- a site with a third allele in either
                # timepoint's samples is dropped, exactly as the simulated side drops it.
                ok = np.isfinite(x) & np.isfinite(y) & ~mx & ~my
                keep = ok & fcc.maf_mask(np.where(ok, x, 0.0), np.where(ok, y, 0.0))
                tot, n = fcc.fc_terms(x[keep], y[keep])
                per_chrom[i][0] += tot
                per_chrom[i][1] += n
        with open(f"{OUTDIR}/chr{chrom}_temporalFc.csv", "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["site_a", "site_b", "t", "n_a", "n_b", "sum_fc", "n_loci"])
            for i, p in enumerate(pairs):
                w.writerow([p["site_a"], p["site_b"], p["t"], p["n_a"], p["n_b"],
                            f"{per_chrom[i][0]:.10g}", per_chrom[i][1]])
                totals[i][0] += per_chrom[i][0]
                totals[i][1] += per_chrom[i][1]
        print(f"  chr{chrom}: {sum(v[1] for v in per_chrom.values())} pair-loci", flush=True)

    # POOLED across chromosomes: sum of per-locus values / total loci. NOT a mean of per-chromosome
    # F_c values (CLAUDE.md 5.2, invariant 4).
    rows = []
    for i, p in enumerate(pairs):
        tot, n = totals[i]
        fc = tot / n if n else float("nan")
        ped = 1.0 / (2.0 * p["n_a"]) + 1.0 / (2.0 * p["n_b"])
        rows.append({**p, "fc": fc, "pedestal": ped, "n_loci": n})

    out = f"{OUTDIR}/averaged_temporalFc.csv"
    with open(out, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["site_a", "site_b", "t", "n_a", "n_b", "fc", "pedestal", "n_loci"])
        for r in rows:
            w.writerow([r["site_a"], r["site_b"], r["t"], r["n_a"], r["n_b"],
                        f"{r['fc']:.10g}", f"{r['pedestal']:.10g}", r["n_loci"]])

    pooled = fcc.pool_by_gap(rows)
    print(f"\n{out}")
    for key in sorted(pooled):
        v = pooled[key]
        print(f"  {key:>5}: F_c = {v['fc_mean']:.5f}  over {v['n_pairs']} field pairs")
    print("\nCopy averaged_temporalFc.csv into data/empiricalStats/ in the repo.")
    print("Keep fc_out/chr*_temporalFc.csv -- they are the only route to a between-chromosome")
    print("jackknife or a re-pool without re-reading every VCF.")


if __name__ == "__main__":
    main()
