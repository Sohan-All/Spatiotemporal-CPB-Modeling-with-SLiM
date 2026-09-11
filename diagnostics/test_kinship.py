"""Known-truth test for ToUseOnBeagles/CalcKinship.py. Seconds, no VCFs, no SLiM.

CLAUDE.md 10: verify any new estimator against known truth BEFORE running it on real data. That
rule has caught real bugs three times, most recently the multi-allelic mismatch that
test_fc_roundtrip.py found on its first run (7.9.8C) -- a 2.6e-3 disagreement, too small to notice
in a results table and far too large to be rounding.

Kinship is the easy case for known truth, because a PEDIGREE gives exact expected values with no
coalescent theory in the way. This builds one directly -- founders drawn at random allele
frequencies, then Mendelian transmission -- so every pair has a textbook kinship:

    full sibs           0.25          parent-offspring    0.25
    half sibs           0.125         unrelated           0.0

FOUR THINGS IT CHECKS
  1. The estimator recovers all four relationships from a hand-built pedigree.
  2. Parent-offspring carries ZERO opposite homozygotes, by Mendelian necessity. This separates it
     from full sibs, which share the same kinship of 0.25 -- and it is a check no amount of
     agreement on the kinship value alone can give.
  3. MULTI-ALLELIC SITES ARE DROPPED, NOT MISCOUNTED. This is the F_c bug class. Adding tri-allelic
     sites must leave every kinship EXACTLY unchanged, because the retained SNP set is identical --
     so this one is asserted at 1e-12, not at the sampling tolerance.
  4. Per-chromosome counts POOL to the genome-wide answer (invariant 4). Splitting the same SNPs
     across two files must give bit-identical counts to one file holding all of them.

It also exercises the real VCF path: the file written here OMITS the ##FORMAT=<ID=GT> header line,
exactly as ConvertBeagleToVCF.py does (5.4), so scikit-allel takes the same default-GT fallback it
takes on the real data.

Run from diagnostics/:  python test_kinship.py
"""
import os
import shutil
import sys
import tempfile
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
TB = ROOT / "ToUseOnBeagles"

L = 30000               # SNPs; sets the sampling error on each kinship estimate
N_TRI = 400             # tri-allelic sites added for check 3
SEED = 7

# Sampling tolerance, NOT float error. KING's SE goes as ~1/sqrt(#het sites); at L=30000 and
# MAF ~ U(0.05,0.5) that is well under 0.01, so 0.025 is a loose but honest bound. Check 3 is
# asserted at 1e-12 instead, because there the retained SNP set is literally identical.
TOL = 0.025

rng = np.random.default_rng(SEED)


def founders(n):
    """n unrelated individuals as (n, 2, L) haplotypes at frequencies p."""
    return rng.random((n, 2, L)) < P


def child(mother, father):
    """One offspring: a random one of each parent's two haplotypes, per SNP (free recombination).

    Linkage is irrelevant to the point estimate -- KING sums per-SNP counts -- so independent
    transmission per site keeps the pedigree exact without simulating a genetic map.
    """
    pick_m = rng.integers(0, 2, L)
    pick_f = rng.integers(0, 2, L)
    return np.stack([np.where(pick_m == 0, mother[0], mother[1]),
                     np.where(pick_f == 0, father[0], father[1])])


P = rng.uniform(0.05, 0.5, L)


def build():
    """A pedigree laid out as four 'fields', with every expected pairwise kinship known."""
    people, site, expect = {}, {}, {}

    # FieldClean: 6 unrelated individuals -> every pair 0.0
    for i in range(6):
        people[f"clean{i}"] = founders(1)[0]
        site[f"clean{i}"] = "FieldClean-2015"

    # FieldSibs: two unrelated parents, 6 full-sib offspring -> every pair 0.25
    mum, dad = founders(2)
    for i in range(6):
        people[f"sib{i}"] = child(mum, dad)
        site[f"sib{i}"] = "FieldSibs-2015"

    # FieldHalf: one shared father, 6 different mothers -> every pair 0.125
    shared_dad = founders(1)[0]
    for i in range(6):
        people[f"half{i}"] = child(founders(1)[0], shared_dad)
        site[f"half{i}"] = "FieldHalf-2015"

    # FieldMixed: 4 unrelated + one parent-offspring pair. The realistic case -- a small related
    # cluster hiding inside an otherwise ordinary sample, which is exactly what 7.2's isolates
    # would look like if clumping is the explanation.
    for i in range(4):
        people[f"mix{i}"] = founders(1)[0]
        site[f"mix{i}"] = "FieldMixed-2015"
    parent = founders(1)[0]
    people["mixP"] = parent
    people["mixC"] = child(parent, founders(1)[0])
    site["mixP"] = site["mixC"] = "FieldMixed-2015"

    def grp(pre, n):
        return [f"{pre}{i}" for i in range(n)]

    for a in grp("clean", 6):
        for b in grp("clean", 6):
            if a < b:
                expect[(a, b)] = 0.0
    for a in grp("sib", 6):
        for b in grp("sib", 6):
            if a < b:
                expect[(a, b)] = 0.25
    for a in grp("half", 6):
        for b in grp("half", 6):
            if a < b:
                expect[(a, b)] = 0.125
    for a in grp("mix", 4):
        for b in grp("mix", 4):
            if a < b:
                expect[(a, b)] = 0.0
    expect[("mixC", "mixP")] = 0.25
    return people, site, expect


def write_vcf(path, people, order, n_tri=0):
    """Minimal phased VCF. No ##FORMAT=<ID=GT> line, matching ConvertBeagleToVCF.py (5.4).

    Tri-allelic sites carry a genotype of 2 in a few individuals. They are ADVERSARIAL on purpose:
    placed so that, if they were silently treated as biallelic, they would inflate the apparent
    opposite-homozygote count and drag kinship down.
    """
    geno = {nm: people[nm][0].astype(int) + people[nm][1].astype(int) for nm in order}
    with open(path, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                + "\t".join(order) + "\n")
        for k in range(L):
            calls = []
            for nm in order:
                g = geno[nm][k]
                calls.append("0|0" if g == 0 else ("0|1" if g == 1 else "1|1"))
            f.write(f"9\t{k + 1}\t.\tA\tT\t.\tPASS\t.\tGT\t" + "\t".join(calls) + "\n")
        for k in range(n_tri):
            calls = []
            for i, nm in enumerate(order):
                calls.append("2|2" if i % 3 == 0 else ("0|0" if i % 3 == 1 else "1|1"))
            f.write(f"9\t{L + k + 1}\t.\tA\tT,G\t.\tPASS\t.\tGT\t" + "\t".join(calls) + "\n")


def run(tmp, vcf_names):
    """Run the REAL script over a directory and return {(a,b): (kinship, n_opphom)}."""
    cwd = os.getcwd()
    os.chdir(tmp)
    sys.path.insert(0, str(tmp))
    for m in ("CalcKinship",):
        sys.modules.pop(m, None)
    import CalcKinship as K
    K.VCF_FILES = vcf_names
    K.POPFILES = {"2015": "popFile2015"}
    K.main()
    os.chdir(cwd)
    out = {}
    with open(Path(tmp) / "kinship_out" / "kinship_pairs_2015.csv") as f:
        import csv as _csv
        for r in _csv.DictReader(f):
            key = tuple(sorted((r["ind_a"], r["ind_b"])))
            out[key] = (float(r["kinship"]), int(r["n_opphom"]))
    return out


def main():
    people, site, expect = build()
    order = sorted(people)
    tmp = Path(tempfile.mkdtemp(prefix="kintest_"))
    fails = []
    try:
        shutil.copy(TB / "CalcKinship.py", tmp / "CalcKinship.py")
        with open(tmp / "popFile2015", "w") as f:
            for nm in order:
                f.write(f"{nm}\t{site[nm]}\n")

        # ---- checks 1 and 2 -------------------------------------------------
        write_vcf(tmp / "chr1_cpb.vcf", people, order)
        got = run(tmp, ["chr1_cpb.vcf"])

        print("\nCHECK 1  pedigree kinship")
        by_class = {}
        for (a, b), exp in expect.items():
            key = tuple(sorted((a, b)))
            if key not in got:
                fails.append(f"pair {key} missing from the output")
                continue
            by_class.setdefault(exp, []).append(got[key][0])
        for exp in sorted(by_class):
            v = np.array(by_class[exp])
            ok = abs(v - exp).max() <= TOL
            print(f"   expected {exp:<6.3f}  n={len(v):>2}  got mean {v.mean():+.5f} "
                  f"range [{v.min():+.5f}, {v.max():+.5f}]   {'OK' if ok else 'FAIL'}")
            if not ok:
                fails.append(f"kinship class {exp}: worst |err| {abs(v - exp).max():.5f} > {TOL}")

        po = got[tuple(sorted(("mixC", "mixP")))]
        sib = got[tuple(sorted(("sib0", "sib1")))]
        print("\nCHECK 2  opposite homozygotes separate parent-offspring from full sibs")
        print(f"   parent-offspring n_opphom {po[1]}  (Mendelian necessity: 0)")
        print(f"   full sibs        n_opphom {sib[1]}  (must be > 0)")
        if po[1] != 0:
            fails.append(f"parent-offspring has {po[1]} opposite homozygotes, must be 0")
        if sib[1] <= 0:
            fails.append("full sibs show no opposite homozygotes; the counting is wrong")

        # ---- check 3: multi-allelic sites are DROPPED -----------------------
        write_vcf(tmp / "chr1_cpb.vcf", people, order, n_tri=N_TRI)
        got_tri = run(tmp, ["chr1_cpb.vcf"])
        worst = max(abs(got_tri[k][0] - got[k][0]) for k in got)
        print(f"\nCHECK 3  {N_TRI} adversarial tri-allelic sites added")
        print(f"   max |kinship change| = {worst:.3e}   (must be ~0: they are dropped, so the")
        print("   retained SNP set is identical)")
        if worst > 1e-12:
            fails.append(f"tri-allelic sites changed kinship by {worst:.3e}; they are being "
                         f"counted, not dropped")

        # ---- check 4: per-chromosome pooling --------------------------------
        lines = open(tmp / "chr1_cpb.vcf").read().splitlines()
        head, body = lines[:2], lines[2:]
        half = len(body) // 2
        for i, part in enumerate((body[:half], body[half:]), start=1):
            with open(tmp / f"chr{i}_split.vcf", "w") as f:
                f.write("\n".join(head + part) + "\n")
        got_split = run(tmp, ["chr1_split.vcf", "chr2_split.vcf"])
        worst2 = max(abs(got_split[k][0] - got_tri[k][0]) for k in got_tri)
        print(f"\nCHECK 4  same SNPs split across two chromosome files")
        print(f"   max |kinship change| = {worst2:.3e}   (counts pool; must be ~0)")
        if worst2 > 1e-12:
            fails.append(f"splitting across files changed kinship by {worst2:.3e}; the "
                         f"per-chromosome accumulation is not pooling")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    print()
    if fails:
        for m in fails:
            print("FAIL:", m)
        sys.exit(1)
    print("all checks passed")


if __name__ == "__main__":
    main()
