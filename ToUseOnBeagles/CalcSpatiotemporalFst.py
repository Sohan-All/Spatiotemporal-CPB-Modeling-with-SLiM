"""Empirical spatio-temporal F_st (CLAUDE.md 7.9.14B): pixy WC F_st for every pair of sites across
all three collections, same-year pairs included, pooled over the 17 chromosomes and written in
stfst_common's per-pair format -- the empirical counterpart of
AnalyzeTreeSeq.calculate_spatiotemporal_fst.

NO VCF MERGE IS NEEDED. Each chr{i}_cpb.vcf.gz comes from ONE Beagle file holding every sample of
every year (the years were imputed together; PixyTheFiles.py already runs each year's popfile
against the same VCFs). Sample IDs are unique across years and every site label carries its year,
so one combined populations file gives pixy all 1830 pairs over identical SNPs -- cross-year and
same-year pairs on the same site set, as the contrast needs. The script checks all of that with
bcftools before running pixy, rather than assuming it.

Steps:
  1. popFileAll_pixy = popFile2015 + popFile2019 + popFile2023 minus sample_exclusions
     (pixy has no exclude option; an absent sample is simply unused). Raises on a duplicate
     sample ID or a site label shared between years.
  2. `bcftools query -l` on every VCF: each popfile sample must be present, and all 17 VCFs must
     carry the same sample list.
  3. pixy --stats fst per chromosome, one window spanning the chromosome (as PixyTheFiles.py)
     -> statsST_chr{i}/chr{i}_all_fst.txt. Old output is deleted first, so a failed pixy run
     leaves nothing for step 4 to re-read.
  4. Pool per pair as AverageData.py pools the production F_st: SNP-count-weighted mean of
     per-chromosome avg_wc_fst (pixy emits no variance components, CLAUDE.md 5). Every one of the
     1830 pairs must be found on every chromosome, or it raises.
  5. Write stfst_out/averaged_spatiotemporalFst.csv (stfst_common.write_pairs; n = diploids after
     sample_exclusions, distances between the real site coordinates).

Copy into data/empiricalStats/ afterwards. stfst_common.spec_hash() printed here must match the
simulated side's (71d30c708bab as of 2026-10-10).

Usage, from the directory holding chr*_cpb.vcf.gz, popFile{year}, specifier_matrix_{year}.csv,
stfst_common.py and sample_exclusions.py:
    python CalcSpatiotemporalFst.py               # pixy (17 runs) + pool
    python CalcSpatiotemporalFst.py --pool-only   # re-pool existing statsST_chr*/ output
"""
import argparse
import os
import subprocess
import sys

try:
    import stfst_common as stc
    from sample_exclusions import EXCLUDED_SAMPLES, exclusion_id
except ImportError as e:
    raise SystemExit(f"{e}. COPY stfst_common.py and sample_exclusions.py next to this script "
                     f"from Python_Code/ -- both sides must run the identical spec (CLAUDE.md 9).")

# --- CONFIG ---------------------------------------------------------------
NUM_CHRS = 17
VCF = "chr{i}_cpb.vcf.gz"
POPFILE_OUT = "popFileAll_pixy"
STATS_DIR = "statsST_chr{i}"
PREFIX = "chr{i}_all"
OUTDIR = "stfst_out"
OUTFILE = "averaged_spatiotemporalFst.csv"
N_CORES = 8
# pixy fst header, checked against every file rather than trusting fixed offsets
FST_HEADER = ["pop1", "pop2", "chromosome", "window_pos_1", "window_pos_2", "avg_wc_fst", "no_snps"]
# --------------------------------------------------------------------------


def write_combined_popfile():
    """All three years' popfiles, exclusions dropped. Returns {sample: site}."""
    members, site_year, dropped = {}, {}, []
    for y in stc.YEARS:
        with open(f"popFile{y}", encoding="utf-8") as f:
            for line in f:
                p = line.split()
                if len(p) < 2:
                    continue
                sid, site = p[0], p[1]
                if sid in EXCLUDED_SAMPLES:
                    dropped.append(sid)
                    continue
                if sid in members:
                    raise SystemExit(f"sample {sid} appears twice ({members[sid]}, {site})")
                if site_year.setdefault(site, y) != y:
                    raise SystemExit(f"site label {site} is used in {site_year[site]} and {y}; "
                                     f"one combined popfile needs year-unique labels")
                members[sid] = site
    with open(POPFILE_OUT, "w", encoding="utf-8") as f:
        for sid, site in members.items():
            f.write(f"{sid}\t{site}\n")
    print(f"{POPFILE_OUT}: {len(members)} samples, {len(site_year)} sites; "
          f"sample_exclusions [{exclusion_id()}] -> dropped {sorted(dropped) or 'none'}")
    return members


def check_vcfs(members):
    ref = None
    for i in range(1, NUM_CHRS + 1):
        vcf = VCF.format(i=i)
        if not os.path.exists(vcf):
            raise SystemExit(f"{vcf} not found")
        samples = subprocess.run(["bcftools", "query", "-l", vcf], check=True,
                                 capture_output=True, text=True).stdout.split()
        if ref is None:
            ref = samples
            missing = sorted(set(members) - set(samples))
            if missing:
                raise SystemExit(f"{vcf}: popfile samples absent from the VCF: {missing}")
        elif samples != ref:
            raise SystemExit(f"{vcf}: sample list differs from {VCF.format(i=1)}")
    print(f"bcftools: all {NUM_CHRS} VCFs carry the same {len(ref)} samples, "
          f"including every popfile sample")


def chrom_end(vcf):
    """Last POS in the file -- pixy's window spans the whole chromosome (as PixyTheFiles.py)."""
    out = subprocess.run(f"zcat {vcf} | tail -n 1", shell=True, check=True,
                         capture_output=True, text=True).stdout.split()
    return out[1]


def run_pixy():
    for i in range(1, NUM_CHRS + 1):
        vcf, folder, prefix = VCF.format(i=i), STATS_DIR.format(i=i), PREFIX.format(i=i)
        os.makedirs(folder, exist_ok=True)
        stale = os.path.join(folder, f"{prefix}_fst.txt")
        if os.path.exists(stale):
            os.remove(stale)
        print(f"pixy fst, chr {i}", flush=True)
        subprocess.run(["pixy", "--stats", "fst", "--vcf", vcf, "--populations", POPFILE_OUT,
                        "--window_size", chrom_end(vcf), "--n_cores", str(N_CORES),
                        "--output_folder", folder, "--output_prefix", prefix,
                        "--bypass_invariant_check", "--silent"], check=True)


def pool():
    """{frozenset(site_a, site_b): (sum fst*snps, sum snps, chromosomes seen)}."""
    acc = {}
    for i in range(1, NUM_CHRS + 1):
        path = os.path.join(STATS_DIR.format(i=i), f"{PREFIX.format(i=i)}_fst.txt")
        with open(path, encoding="utf-8") as f:
            header = f.readline().split()
            if header != FST_HEADER:
                raise SystemExit(f"{path}: header {header} != expected {FST_HEADER}")
            for line in f:
                p = line.split()
                if not p:
                    continue
                if p[5] in ("NA", "nan") or p[6] in ("NA", "nan"):
                    continue
                key = frozenset((p[0], p[1]))
                s = acc.setdefault(key, [0.0, 0.0, set()])
                s[0] += float(p[5]) * float(p[6])
                s[1] += float(p[6])
                s[2].add(i)
    return acc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pool-only", action="store_true")
    args = ap.parse_args()
    print(f"stfst_common spec {stc.spec_hash()} -- the simulated side must print the same hash")

    members = write_combined_popfile()
    if not args.pool_only:
        check_vcfs(members)
        run_pixy()

    acc = pool()
    sites = stc.site_table(".")
    rows, short = [], []
    for a, b, gap, dist in stc.all_pairs(sites):
        sa, sb = sites[a], sites[b]
        s = acc.get(frozenset((sa["site"], sb["site"])))
        if s is None or s[1] == 0:
            raise SystemExit(f"no pixy F_st for {sa['site']} x {sb['site']}")
        if len(s[2]) != NUM_CHRS:
            short.append((sa["site"], sb["site"], len(s[2])))
        rows.append({"year_a": sa["year"], "site_a": sa["site"], "year_b": sb["year"],
                     "site_b": sb["site"], "gap": gap, "dist_km": dist, "n_a": sa["n"],
                     "n_b": sb["n"], "fst": s[0] / s[1]})
    if short:
        raise SystemExit(f"{len(short)} pairs missing chromosomes, e.g. {short[:5]}")

    os.makedirs(OUTDIR, exist_ok=True)
    out = os.path.join(OUTDIR, OUTFILE)
    stc.write_pairs(out, rows)
    print(f"{out}: {len(rows)} pairs, all pooled over {NUM_CHRS} chromosomes")

    # Preview only -- the written file keeps every pair. Mirrors ABCAnalysisNoRedis.EXCLUDED_SITES
    # (not importable here); the repo side re-summarises with the real constant.
    keep = stc.retained_sites(".", exclude_sites={"Mortensen9-2015"})
    summ = stc.summarise(rows, keep)
    for g in (0, 8, 16):
        print(f"  gap {g:2d}: " + "  ".join(
            f"{stc.bin_label(k)} {summ[(g, k)]['median']:.5f} (n={summ[(g, k)]['n_pairs']})"
            for k in range(len(stc.DIST_EDGES_KM) - 1)))
    con = stc.contrast(summ)
    for g in (8, 16):
        print(f"  contrast gap {g:2d}: " + "  ".join(
            f"{stc.bin_label(k)} {con[(g, k)]:+.5f}" for k in range(len(stc.DIST_EDGES_KM) - 1)))


if __name__ == "__main__":
    sys.exit(main())
