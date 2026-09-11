"""Round-trip regression test for the EMPIRICAL temporal-F_c path. Seconds; no SLiM, no VCFs.

CLAUDE.md 10: verify a new statistic against a known-truth simulation BEFORE running it on real
data. That convention has now caught a real bug three times, and this is the fourth -- see below.

WHAT IT DOES. Simulates a small tree sequence, writes it as a VCF with popfiles and specifier
matrices laid out exactly as the Beagle machine has them, runs `ToUseOnBeagles/CalcTemporalFc.py`'s
real `main()` over it, and compares every field pair against F_c computed straight from the
genotype matrix.

WHAT IT CAUGHT, and why it is worth keeping. The first run disagreed by up to 2.6e-3 in F_c -- too
small to notice by eye in a result table, far too large to be rounding. Cause: the two sides
handled MULTI-ALLELIC sites differently. The simulated side folded every derived allele together
as `G > 0`; the empirical side passed `max_allele=1` to scikit-allel, which does not exclude a
third allele but silently shrinks that site's denominator. **7.1% of sites in a real simulated
tree carry more than two alleles**, so this was not a corner case. Both sides now drop
non-biallelic sites per pair (invariant 7), through `fc_common`.

TOLERANCE. The two paths accumulate in different orders -- the script sums per chunk per
chromosome, the reference sums once -- so float64 summation noise of ~1e-10 is expected and
meaningless. Anything above TOL is a difference in LOGIC, not in arithmetic.

Run:  python test_fc_roundtrip.py
"""
import os, sys, shutil, tempfile, csv
from pathlib import Path
import numpy as np, msprime, allel

ROOT = Path(__file__).resolve().parent.parent
PY = ROOT / "Python_Code"
TB = ROOT / "ToUseOnBeagles"
sys.path.insert(0, str(PY))
import fc_common as fcc

# Exercise the EXCLUSION path (CLAUDE.md 7.2.2F) with synthetic ids. One sample goes from each side
# from FieldB-2019 and one from FieldA-2023; each field sits in two year-pairs, so n drops to 6 in
# FOUR pairs, and the script's columns, its n_a/n_b, and the direct computation below must all
# agree on which beetles are gone.
fcc.EXCLUDED_SAMPLES.clear()
fcc.EXCLUDED_SAMPLES.update({"ind23": ("FieldB-2019", "round-trip test"),
                             "ind30": ("FieldA-2023", "round-trip test")})

YEARS =["2015", "2019", "2023"]
SITES = {"2015": ["FieldA-2015", "FieldB-2015"],
         "2019": ["FieldA-2019", "FieldB-2019"],
         "2023": ["FieldA-2023", "FieldB-2023"]}
COORDS = {"FieldA": (44.10000, -89.50000), "FieldB": (44.20000, -89.60000)}
N = 7

# float64 summation-order noise only; see TOLERANCE above.
TOL = 1e-8

tmp = Path(tempfile.mkdtemp(prefix="fctest_"))
try:
    ts = msprime.sim_ancestry(samples=len(YEARS) * 2 * N, ploidy=2, sequence_length=5e4,
                              recombination_rate=1e-7, population_size=800, random_seed=11)
    ts = msprime.sim_mutations(ts, rate=2e-6, random_seed=11)
    names, member = [], {}
    k = 0
    for y in YEARS:
        for s in SITES[y]:
            member[s] = []
            for _ in range(N):
                nm = f"ind{k}"
                names.append(nm); member[s].append(nm); k += 1
    with open(tmp / "chr1_cpb.vcf", "w") as f:
        ts.write_vcf(f, individual_names=names)
    for y in YEARS:
        with open(tmp / f"specifier_matrix_{y}.csv", "w") as f:
            for s in SITES[y]:
                la, lo = COORDS[s.split("-")[0]]
                f.write(f"{s},{la},{lo}\n")
        with open(tmp / f"popFile{y}", "w") as f:
            for s in SITES[y]:
                for nm in member[s]:
                    f.write(f"{nm}\t{s}\n")

    # ---- run the real script ---------------------------------------------------------------
    shutil.copy(PY / "fc_common.py", tmp / "fc_common.py")
    shutil.copy(TB / "CalcTemporalFc.py", tmp / "CalcTemporalFc.py")
    cwd = os.getcwd(); os.chdir(tmp); sys.path.insert(0, str(tmp))
    import CalcTemporalFc as C
    C.VCF_FILES = ["chr1_cpb.vcf"]
    C.main()
    os.chdir(cwd)

    got = {}
    with open(tmp / "fc_out" / "averaged_temporalFc.csv") as f:
        for r in csv.DictReader(f):
            got[(r["site_a"], r["site_b"])] = float(r["fc"])

    # ---- independent computation straight from the genotype matrix -------------------------
    G = allel.GenotypeArray(ts.genotype_matrix().reshape(ts.num_sites, len(names), 2))
    pos = {nm: i for i, nm in enumerate(names)}
    def freq(site):
        """Frequency of ALLELE 1, and whether any allele index above 1 appears.

        Must mirror the spec exactly (fc_common.waples_fc): NOT to_n_alt(), which folds every
        non-reference allele together and is precisely the bug this round-trip first caught.
        """
        cols = [pos[nm] for nm in member[site] if nm not in fcc.EXCLUDED_SAMPLES]
        sub = np.asarray(G[:, cols, :])            # (sites, n, 2) allele indices
        f = (sub == 1).sum(axis=(1, 2)) / (2.0 * len(cols))
        multi = (sub > 1).any(axis=(1, 2))
        return f, multi
    pairs, _ = fcc.matched_field_pairs(
        {y: str(tmp / f"specifier_matrix_{y}.csv") for y in YEARS},
        {y: str(tmp / f"popFile{y}") for y in YEARS})

    print(f"\nchecking {len(pairs)} pairs against a direct computation")
    worst = 0.0
    for p in pairs:
        x, mx = freq(p["site_a"])
        y, my = freq(p["site_b"])
        fc, _, _, nloc = fcc.waples_fc(x, y, p["n_a"], p["n_b"], p["t"], valid=~mx & ~my)
        d = abs(fc - got[(p["site_a"], p["site_b"])])
        worst = max(worst, d)
        print(f"  {p['site_a']:>12} -> {p['site_b']:<12} t={p['t']:2d}  "
              f"script {got[(p['site_a'], p['site_b'])]:.10f}  direct {fc:.10f}  diff {d:.2e}"
              f"  ({nloc} loci)")
    # The exclusion must shrink n in the script's output exactly as it shrinks the sample sets.
    n_ok = True
    with open(tmp / "fc_out" / "averaged_temporalFc.csv") as f:
        for r in csv.DictReader(f):
            want = tuple(sum(nm not in fcc.EXCLUDED_SAMPLES for nm in member[s])
                         for s in (r["site_a"], r["site_b"]))
            if (int(r["n_a"]), int(r["n_b"])) != want:
                n_ok = False
                print(f"  n MISMATCH {r['site_a']}->{r['site_b']}: script {r['n_a']}/{r['n_b']}, "
                      f"expected {want[0]}/{want[1]}")
    n_excl = sum(1 for p in pairs if 6 in (p["n_a"], p["n_b"]))
    print(f"exclusions: {n_excl} pairs carry n=6, n_a/n_b consistent = {n_ok}")
    n_ok = n_ok and n_excl == 4

    # A typo or a wrong site in the list must RAISE rather than silently exclude nothing.
    guard_ok = True
    saved = dict(fcc.EXCLUDED_SAMPLES)
    for bad in ({"nobody": ("FieldA-2015", "absent id")}, {"ind0": ("FieldB-2015", "wrong site")}):
        fcc.EXCLUDED_SAMPLES.update(bad)
        try:
            fcc.matched_field_pairs({y: str(tmp / f"specifier_matrix_{y}.csv") for y in YEARS},
                                    {y: str(tmp / f"popFile{y}") for y in YEARS})
            guard_ok = False
            print(f"  guard did NOT raise for {bad}")
        except ValueError:
            pass
        fcc.EXCLUDED_SAMPLES.clear()
        fcc.EXCLUDED_SAMPLES.update(saved)
    print(f"exclusion typo guard raises = {guard_ok}")

    ok = worst < TOL and n_ok and guard_ok
    print()
    print(f"max |script - direct| = {worst:.3e}  (tolerance {TOL:.0e})")
    print("PASS" if ok else "*** FAIL ***")
    sys.exit(0 if ok else 1)
finally:
    shutil.rmtree(tmp, ignore_errors=True)
