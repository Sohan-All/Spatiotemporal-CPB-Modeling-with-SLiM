#!/usr/bin/env python3
"""Individual-level kinship within each sampled field -- the check on CLAUDE.md 7.2's isolates.

Runs on the machine holding the Beagle VCFs, NOT in the repo.

--------------------------------------------------------------------------------------------
THE QUESTION, AND WHY IT IS NOT THE ONE IT LOOKS LIKE
--------------------------------------------------------------------------------------------
Two site-samples -- `Mortensen9-2015` and `H41-2023` -- are simultaneously the lowest-pi and
(near-)highest self-relatedness sample of their year, with mean F_st 3-5.7x the next site (7.2).
Removing them moves the F_st-implied POPMULT from ~12,500 to ~26,400, a factor of 2.12, and it
takes the three years from disagreeing 2.88x to agreeing 1.35x -- with 2019, which has nothing to
remove, sitting still in the middle. So something is wrong with those two samples.

**This script is NOT here to decide whether they are siblings or a founder event.** That
distinction changes nothing: the simulation draws its n_i individuals at random from a deme of
several hundred, so it cannot produce either one, and the sample is unmatched either way.

It is here to answer a question nothing else can: **HOW MANY samples are affected.** Per-site pi
and relatedness only reveal OUTLIERS -- a site is flagged by being unlike the others. If every
field sample is somewhat clumped, no site stands out and the whole collection quietly violates the
model's sampling assumption. That is the difference between a correction you apply to two sites
and a caveat you attach to every number in the project, and it is the reason to run this.

The prior evidence points at "only a few", and this is the test of it: pixy's WC F_st scatters
properly through zero (24/24/64 negative pairs), which a large uniform inflation could not do.
But that bounds the AVERAGE effect, not the COUNT of affected samples.

--------------------------------------------------------------------------------------------
FOUR THINGS THAT WILL SILENTLY RUIN THIS
--------------------------------------------------------------------------------------------
1. **BEAGLE IMPUTATION BIASES KINSHIP UPWARD, and there is no correcting it here.** Imputation
   fills genotypes by copying haplotypes from the panel's own samples, which makes individuals
   look more alike than they are. So read the output ASYMMETRICALLY: kinship near zero everywhere
   is trustworthy and settles the question; high kinship is suggestive and cannot on its own prove
   relatedness. What survives the bias is the COMPARISON between sites -- every sample went
   through the identical imputation -- and the comparison is what the decision needs.
2. **The reference is NOT zero.** Individuals from one deme of a few hundred are genuinely related
   at some background level, so "kinship > 0" means nothing on its own. This script therefore
   computes the SAME estimator on BETWEEN-site pairs of the same year as a baseline, on the same
   markers, through the same imputation. Judge a site against that baseline, never against 0.
3. **MAF is computed over the WHOLE YEAR's samples, not per site.** Filtering per site would
   condition on the very quantity under test: a clumped sample has distorted allele frequencies,
   so a per-site MAF cut would preferentially discard exactly the SNPs that reveal it.
4. **Pool, never average** (5.2, invariant 4). KING is a ratio of counts summed over SNPs, so
   genome-wide = (sum of numerators over all 17 chromosomes) / (sum of denominators). It is NOT
   the mean of 17 per-chromosome kinships. This script accumulates counts.

--------------------------------------------------------------------------------------------
THE ESTIMATOR
--------------------------------------------------------------------------------------------
KING-robust (Manichaikul et al. 2010, Bioinformatics 26:2867), the between-family form:

    phi_ij = ( N_AaAa - 2 * N_AAaa ) / ( N_Aa_i + N_Aa_j )

N_AaAa  = SNPs where BOTH individuals are heterozygous
N_AAaa  = SNPs where they are OPPOSITE homozygotes
N_Aa_i  = SNPs where individual i is heterozygous

Chosen over a straight allele-sharing or GRM estimate for one specific reason: **it uses no
allele frequencies at all**, so it does not need a reference population and is robust to the
population structure that is itself under investigation. A frequency-based estimator computed
against the year-wide frequencies would confound "this pair is related" with "this site is
differentiated", which is precisely the confusion to avoid here.

Expected values: duplicate 0.5, parent-offspring or full sib 0.25, half sib 0.125, first cousin
0.0625, unrelated 0. Degree cut points are the paper's own powers of two (see DEGREE_CUTS).

NEGATIVE KINSHIP IS EXPECTED AND IS NOT A BUG. KING-robust uses no allele frequencies, which is
what makes it robust to structure, but the price is that a pair drawn from two DIFFERENT
populations returns a systematically negative value -- the more differentiated the pair, the more
negative. So the between-site baseline will sit slightly below zero and its spread across pairs is
partly real differentiation, not noise. That is exactly why it is reported as a baseline rather
than assumed to be 0.

Note KING assumes SNPs are independent; ours are not. That inflates nothing in the point estimate
-- it only means the standard error is smaller than the SNP count suggests, which does not matter
at 36M markers where the estimate is precise regardless.

--------------------------------------------------------------------------------------------
USAGE
    python CalcKinship.py                 # one pass over the 17 VCFs
Outputs (kinship_out/)
    chr{i}_kinshipCounts.csv    per-chromosome counts for WITHIN-site pairs only -- KEEP THESE,
                                the standing rule (7.8.4 needed exactly this for LD). Between-site
                                pairs are ~50x more numerous and are only ever wanted pooled.
    kinship_pairs_{year}.csv    every pair, pooled genome-wide: counts, kinship, degree call
    kinship_site_summary.csv    THE OUTPUT THAT ANSWERS THE QUESTION -- one row per site, with its
                                within-site kinship distribution against that year's between-site
                                baseline, and how many pairs land in each KING degree class.

If a simulated null is ever wanted -- what within-deme kinship looks like when there is no family
structure by construction -- run the identical estimator on the simulated deme samples. That would
need a shared spec module and a spec_hash on both sides, exactly as ld_common/fc_common do. There
is no such counterpart today and this script deliberately does not pretend otherwise.
"""
import csv
import glob
import os
import re
from collections import defaultdict

import numpy as np
import allel

VCF_FILES = sorted(glob.glob("chr*_cpb.vcf.gz"))
POPFILES = {"2015": "popFile2015", "2019": "popFile2019", "2023": "popFile2023"}
OUTDIR = "kinship_out"

MIN_MAF = 0.05          # matches ld_common.MIN_MAF / fc_common.MIN_MAF
CHUNK_SNPS = 20000      # per-chunk SNPs; the matmuls below are (n_ind x CHUNK) so this is cheap

# Manichaikul et al. 2010 Table 1 -- the powers-of-two boundaries, not eyeballed round numbers.
DEGREE_CUTS = ((0.354, "duplicate"),
               (0.177, "1st"),
               (0.0884, "2nd"),
               (0.0442, "3rd"))


def degree(phi):
    for cut, name in DEGREE_CUTS:
        if phi > cut:
            return name
    return "unrelated"


def read_popfile(path):
    """{sample_id: site} for one year."""
    out = {}
    with open(path, encoding="utf-8") as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 2:
                out[parts[0].strip()] = parts[1].strip()
    if not out:
        raise SystemExit(f"{path} produced no sample->site rows")
    return out


def king_counts(gt, cols):
    """(N_AaAa, N_AAaa, N_Aa per individual, n_snps) for one chunk, over the given VCF columns.

    Done as three boolean matmuls rather than a pair loop: with ~150 individuals and a 20k-SNP
    chunk the products are (150 x 20000) @ (20000 x 150), which is milliseconds, against ~11k
    pairwise passes over the chunk. The counts are additive over chunks and over chromosomes,
    which is what makes the genome-wide pooling in main() exact.

    Multi-allelic and non-polymorphic sites are dropped before counting (invariant 7): a third
    allele makes "opposite homozygote" ambiguous, and monomorphic sites contribute zero to every
    term while still costing time.
    """
    g = gt.take(cols, axis=1)
    ac = g.count_alleles(max_allele=3)
    biallelic = (ac[:, 2:].sum(axis=1) == 0)

    an = ac[:, :2].sum(axis=1).astype(float)
    with np.errstate(invalid="ignore", divide="ignore"):
        maf = np.where(an > 0, np.minimum(ac[:, 0], ac[:, 1]) / an, 0.0)
    keep = biallelic & (maf >= MIN_MAF)
    if not keep.any():
        n = len(cols)
        return (np.zeros((n, n), np.int64), np.zeros((n, n), np.int64),
                np.zeros(n, np.int64), 0)

    g = g.compress(keep, axis=0)
    het = np.ascontiguousarray(g.is_het().T.astype(np.float32))        # (n_ind, n_snp)
    hom_ref = np.ascontiguousarray(g.is_hom_ref().T.astype(np.float32))
    hom_alt = np.ascontiguousarray(g.is_hom_alt().T.astype(np.float32))

    n_hethet = het @ het.T
    opp = hom_ref @ hom_alt.T
    n_opphom = opp + opp.T
    n_het = het.sum(axis=1)
    return (n_hethet.astype(np.int64), n_opphom.astype(np.int64),
            n_het.astype(np.int64), int(keep.sum()))


def main():
    if not VCF_FILES:
        raise SystemExit("no chr*_cpb.vcf.gz found in the working directory")
    os.makedirs(OUTDIR, exist_ok=True)

    _, samples, _, _ = allel.iter_vcf_chunks(VCF_FILES[0], fields=["variants/POS"])
    sample_pos = {s: i for i, s in enumerate(samples)}

    # ---- who is in which year and site ------------------------------------
    year_cols, year_site, year_ids = {}, {}, {}
    for year, pf in POPFILES.items():
        s2site = read_popfile(pf)
        missing = [s for s in s2site if s not in sample_pos]
        if missing:
            raise SystemExit(f"{pf}: {len(missing)} samples are not in the VCF, e.g. {missing[:3]}")
        ids = sorted(s2site, key=lambda s: sample_pos[s])
        year_ids[year] = ids
        year_cols[year] = [sample_pos[s] for s in ids]
        year_site[year] = [s2site[s] for s in ids]
        n_sites = len(set(year_site[year]))
        print(f"{year}: {len(ids)} individuals over {n_sites} sites")

    # per-chromosome files hold WITHIN-site pairs only (see the docstring); everything is
    # accumulated genome-wide regardless.
    tot = {y: [np.zeros((len(year_ids[y]),) * 2, np.int64),
               np.zeros((len(year_ids[y]),) * 2, np.int64),
               np.zeros(len(year_ids[y]), np.int64), 0] for y in POPFILES}

    for vcf in VCF_FILES:
        chrom = re.search(r"chr[^0-9]*([0-9]+)", os.path.basename(vcf)).group(1)
        per = {y: [np.zeros_like(tot[y][0]), np.zeros_like(tot[y][1]),
                   np.zeros_like(tot[y][2]), 0] for y in POPFILES}
        _, _, _, chunks = allel.iter_vcf_chunks(
            vcf, fields=["calldata/GT"], chunk_length=CHUNK_SNPS)
        for chunk in chunks:
            gt = allel.GenotypeArray(chunk[0]["calldata/GT"])
            for y in POPFILES:
                hh, oh, h, n = king_counts(gt, year_cols[y])
                per[y][0] += hh
                per[y][1] += oh
                per[y][2] += h
                per[y][3] += n

        with open(f"{OUTDIR}/chr{chrom}_kinshipCounts.csv", "w", newline="",
                  encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["year", "site", "ind_a", "ind_b", "n_hethet", "n_opphom",
                        "n_het_a", "n_het_b", "n_snps"])
            for y in POPFILES:
                ids, sites = year_ids[y], year_site[y]
                for i in range(len(ids)):
                    for j in range(i + 1, len(ids)):
                        if sites[i] != sites[j]:
                            continue
                        w.writerow([y, sites[i], ids[i], ids[j], per[y][0][i, j],
                                    per[y][1][i, j], per[y][2][i], per[y][2][j], per[y][3]])
                for k in range(4):
                    tot[y][k] += per[y][k]
        print(f"  chr{chrom}: " +
              ", ".join(f"{y} {per[y][3]} SNPs" for y in sorted(POPFILES)), flush=True)

    # ---- pooled kinship, per year -----------------------------------------
    summary = []
    for y in sorted(POPFILES):
        hh, oh, h, n_snps = tot[y]
        ids, sites = year_ids[y], year_site[y]
        rows, within, between = [], defaultdict(list), []
        for i in range(len(ids)):
            for j in range(i + 1, len(ids)):
                den = h[i] + h[j]
                phi = (hh[i, j] - 2.0 * oh[i, j]) / den if den else float("nan")
                same = sites[i] == sites[j]
                rows.append((ids[i], ids[j], sites[i], sites[j], int(same),
                             hh[i, j], oh[i, j], h[i], h[j], phi, degree(phi)))
                (within[sites[i]] if same else between).append(phi)

        with open(f"{OUTDIR}/kinship_pairs_{y}.csv", "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["ind_a", "ind_b", "site_a", "site_b", "within_site", "n_hethet",
                        "n_opphom", "n_het_a", "n_het_b", "kinship", "degree"])
            for r in rows:
                w.writerow(list(r[:9]) + [f"{r[9]:.6f}", r[10]])

        base = np.array(between, float)
        print(f"\n{y}: {n_snps} SNPs, between-site baseline "
              f"median {np.median(base):+.5f}, 99th pct {np.percentile(base, 99):+.5f}, "
              f"max {base.max():+.5f}")
        print(f"{'site':<22} {'n':>3} {'pairs':>6} {'median':>9} {'max':>9} "
              f"{'>=3rd':>6} {'>=2nd':>6} {'>=1st':>6}  vs baseline")
        for site in sorted(within):
            v = np.array(within[site], float)
            if v.size == 0:
                continue
            counts = [int((v > c).sum()) for c in (0.0442, 0.0884, 0.177)]
            # how unusual is this site's median against the same year's between-site pairs?
            pct = float((base < np.median(v)).mean() * 100)
            summary.append(dict(year=y, site=site, n=int(np.isfinite(v).sum()),
                                n_ind=sites.count(site), median=float(np.median(v)),
                                mean=float(v.mean()), max=float(v.max()),
                                ge_3rd=counts[0], ge_2nd=counts[1], ge_1st=counts[2],
                                baseline_pct=pct, n_snps=n_snps))
            print(f"{site:<22} {sites.count(site):>3} {v.size:>6} {np.median(v):>+9.5f} "
                  f"{v.max():>+9.5f} {counts[0]:>6} {counts[1]:>6} {counts[2]:>6}  "
                  f"{pct:>5.1f} pct")

    with open(f"{OUTDIR}/kinship_site_summary.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(summary[0].keys()))
        w.writeheader()
        for r in summary:
            w.writerow(r)

    n_flag = sum(1 for r in summary if r["ge_2nd"] > 0)
    print(f"\n{n_flag} of {len(summary)} site-samples contain at least one 2nd-degree-or-closer "
          f"pair.")
    print("HOW TO READ THIS, and it is the whole point of the run:")
    print("  * a handful of flagged sites  -> the clumping is local. Dropping those sites is a")
    print("    correction, and the 2.12x shift in inferred N is real.")
    print("  * most or all sites flagged   -> the collection systematically violates the model's")
    print("    sampling assumption. Do NOT drop sites; report it as a caveat on every estimate.")
    print("  * nothing flagged anywhere    -> clumping is not the explanation for 7.2's isolates,")
    print("    and a founder/bottleneck event at those fields becomes the leading account.")
    print("Remember imputation biases these UPWARD (docstring pt 1): a null result is strong, a")
    print("positive one is suggestive. Judge every site against its own year's baseline, not 0.")
    print(f"\nCopy {OUTDIR}/kinship_site_summary.csv into the repo under data/ ; keep the")
    print(f"per-chromosome counts in {OUTDIR}/ for a jackknife or a re-pool.")


if __name__ == "__main__":
    main()
