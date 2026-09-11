"""Per-individual heterozygote-miscall rates and corrected KING kinship -- CLAUDE.md 7.2.2.
Seconds; no VCFs, no simulation. Reads out/kinship_out/ from ToUseOnBeagles/CalcKinship.py.

WHY THE RAW KINSHIP CANNOT BE READ. Unrelated pairs come out at about -0.30, not 0, because the
genotype calls carry a per-individual heterozygote deficit. Model: individual i has a rate e_i at
which a true heterozygote is called homozygous for a random allele. Then, exactly in expectation
under HWE within a year,

    KING_obs  = (phi_true - ebar) / (1 - ebar),     ebar = (e_i + e_j) / 2
    het_obs,i = (1 - e_i) * H

so e_i is fitted by least squares on BETWEEN-site pairs (phi_true ~ 0) and every pair is corrected
with phi_true = phi_obs * (1 - ebar) + ebar. The check the fit does not enforce: raw heterozygote
counts should fall into line with (1 - e_i), i.e. het/(1-e) should be nearly constant.

The counts cannot tell miscalling from genuine inbreeding (identical per individual); miscalling
is the inferred reading. The correction is ill-conditioned as e -> 1 -- never trust a corrected
value involving e > FAILED_E.

ALSO THE SOURCE OF fc_common.EXCLUDED_SAMPLES. derive_exclusions() applies the rule written next to
that list; main() exits non-zero if the two disagree, so the list cannot drift from the data.

AND OF THE SIMULATED SIDE'S MISCALL RATES. --write-rates writes every sample's e to
data/empiricalStats/miscall_rates.csv, which AnalyzeTreeSeq.calculate_temporal_fc reads to imitate
the miscalls in simulated genotypes (fc_common.apply_miscall).

Run from diagnostics/:  python kinship_correct.py [--write-rates]
"""
import argparse
import csv
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
KO = ROOT / "out" / "kinship_out"
GD = ROOT / "data" / "Genetic_Data"
sys.path.insert(0, str(ROOT / "Python_Code"))
import fc_common as fcc  # noqa: E402

YEARS = ("2015", "2019", "2023")
FAILED_E = 0.4           # e above this = failed sample (see fc_common.EXCLUDED_SAMPLES)
THIRD_DEGREE = 0.0442    # KING's third-degree boundary (Manichaikul et al. 2010) -- reporting only
FIRST_DEGREE = 0.177     # KING's first-degree boundary -- the F_c exclusion cut (fc_common: the
                         # simulation makes half-sibs naturally but never full sibs)
DEGREE_CUTS = ((0.354, "dup"), (0.177, "1st"), (0.0884, "2nd"), (0.0442, "3rd"))


def _degree(phi):
    for cut, name in DEGREE_CUTS:
        if phi > cut:
            return name
    return "-"


def _n_snps():
    with open(KO / "kinship_site_summary.csv", encoding="utf-8") as f:
        return {r["year"]: int(r["n_snps"]) for r in csv.DictReader(f)}


def fit_year(year):
    """Fit e_i for one year and correct every pair. Returns a dict of arrays keyed by pair index,
    plus per-individual `names`, `site_of`, `het` (fraction of SNPs) and `e`."""
    with open(KO / f"kinship_pairs_{year}.csv", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    ids, site_of, het_n = {}, {}, {}
    for r in rows:
        for k, s, h in (("ind_a", "site_a", "n_het_a"), ("ind_b", "site_b", "n_het_b")):
            if r[k] not in ids:
                ids[r[k]] = len(ids)
                site_of[r[k]] = r[s]
                het_n[r[k]] = int(r[h])
    names = list(ids)
    a = np.array([ids[r["ind_a"]] for r in rows])
    b = np.array([ids[r["ind_b"]] for r in rows])
    phi = np.array([float(r["kinship"]) for r in rows])
    within = np.array([r["within_site"] == "1" for r in rows])

    bw = ~within
    target = -phi[bw] / (1.0 - phi[bw])
    X = np.zeros((int(bw.sum()), len(names)))
    X[np.arange(X.shape[0]), a[bw]] = 0.5
    X[np.arange(X.shape[0]), b[bw]] = 0.5
    e = np.linalg.lstsq(X, target, rcond=None)[0]
    resid = target - X @ e
    r2 = 1.0 - (resid ** 2).sum() / ((target - target.mean()) ** 2).sum()

    ebar = 0.5 * (e[a] + e[b])
    het = np.array([het_n[s] for s in names], float) / _n_snps()[year]
    return dict(year=year, names=names, site_of=site_of, het=het, e=e, a=a, b=b, phi=phi,
                within=within, corrected=phi * (1.0 - ebar) + ebar, r2=r2)


def miscall_rates():
    """{sample_id: (year, site, e)} over all three years."""
    out = {}
    for y in YEARS:
        fy = fit_year(y)
        for j, s in enumerate(fy["names"]):
            out[s] = (y, fy["site_of"][s], float(fy["e"][j]))
    return out


def _matched_sites():
    """{year: set of site names} that appear in any temporal-F_c field pair."""
    pairs, _ = fcc.matched_field_pairs(
        {y: GD / f"specifier_matrix_{y}.csv" for y in YEARS},
        {y: GD / f"popFile{y}" for y in YEARS})
    out = {y: set() for y in YEARS}
    for p in pairs:
        out[p["year_a"]].add(p["site_a"])
        out[p["year_b"]].add(p["site_b"])
    return out


def derive_exclusions(fits, matched):
    """Apply fc_common's written rule. Returns {sample_id: reason}."""
    out = {}
    for y in YEARS:
        fy = fits[y]
        names, e, site_of = fy["names"], fy["e"], fy["site_of"]
        med = float(np.median(e))
        for j, s in enumerate(names):
            if site_of[s] in matched[y] and e[j] > FAILED_E:
                out[s] = f"failed e={e[j]:.3f}"
        cand = [k for k in np.where(fy["within"])[0]
                if site_of[names[fy["a"][k]]] in matched[y]
                and max(e[fy["a"][k]], e[fy["b"][k]]) <= FAILED_E
                and fy["corrected"][k] > FIRST_DEGREE]
        for k in sorted(cand, key=lambda k: -fy["corrected"][k]):
            i, j = fy["a"][k], fy["b"][k]
            if names[i] in out or names[j] in out:
                continue
            drop = i if abs(e[i] - med) >= abs(e[j] - med) else j
            keep = j if drop == i else i
            out[names[drop]] = (f"relative {_degree(fy['corrected'][k])} "
                                f"({fy['corrected'][k]:+.3f}) with {names[keep]} (kept)")
    return out


def write_rates(fits, path):
    """Every sample's fitted e, for the simulated side (fc_common.read_miscall_rates)."""
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["sample_id", "year", "site", "e"])
        for y in YEARS:
            fy = fits[y]
            for j, s in enumerate(fy["names"]):
                w.writerow([s, y, fy["site_of"][s], f"{fy['e'][j]:.6f}"])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--write-rates", action="store_true",
                    help=f"write data/empiricalStats/{fcc.MISCALL_RATES_FILE}")
    args = ap.parse_args()
    fits = {y: fit_year(y) for y in YEARS}
    print("PER-INDIVIDUAL HET-MISCALL FIT (between-site pairs only)")
    print(f"{'year':>5} {'inds':>5} {'R^2':>7} {'median e':>9} {'e range':>17} "
          f"{'corr(het,1-e)':>14} {'CV het -> het/(1-e)':>20} {'true het':>9}")
    for y in YEARS:
        fy = fits[y]
        e, het = fy["e"], fy["het"]
        H = het / (1.0 - e)
        print(f"{y:>5} {len(e):>5} {fy['r2']:>7.4f} {np.median(e):>9.3f} "
              f"[{e.min():+.3f}, {e.max():+.3f}] {np.corrcoef(het, 1 - e)[0, 1]:>14.4f} "
              f"{het.std() / het.mean():>9.3f} -> {H.std() / H.mean():.3f} {H.mean():>9.4f}")

    print("\nFAILED SAMPLES (e > %.1f) -- corrected kinship involving these is meaningless" % FAILED_E)
    for y in YEARS:
        fy = fits[y]
        for j in np.argsort(-fy["e"]):
            if fy["e"][j] <= FAILED_E:
                break
            s = fy["names"][j]
            print(f"  {s:>6} {fy['site_of'][s]:<24} e={fy['e'][j]:.3f}  observed het {fy['het'][j]:.3f}")

    print(f"\nWITHIN-SITE PAIRS ABOVE THE THIRD-DEGREE CUT (both e <= {FAILED_E})")
    for y in YEARS:
        fy = fits[y]
        n, e = fy["names"], fy["e"]
        for k in sorted(np.where(fy["within"])[0], key=lambda k: -fy["corrected"][k]):
            i, j = fy["a"][k], fy["b"][k]
            if fy["corrected"][k] <= THIRD_DEGREE:
                break
            if max(e[i], e[j]) > FAILED_E:
                continue
            print(f"  {n[i]:>6} x {n[j]:<6} {fy['site_of'][n[i]]:<24} obs {fy['phi'][k]:+.4f} "
                  f"e {e[i]:+.3f}/{e[j]:+.3f} -> {fy['corrected'][k]:+.4f} {_degree(fy['corrected'][k])}")

    derived = derive_exclusions(fits, _matched_sites())
    listed = set(fcc.EXCLUDED_SAMPLES)
    print("\nEXCLUSION RULE applied to matched F_c fields:")
    for s, why in sorted(derived.items()):
        print(f"  {s:>6}  {why}")
    if set(derived) != listed:
        print(f"\n*** fc_common.EXCLUDED_SAMPLES DISAGREES with the rule: only in the rule "
              f"{sorted(set(derived) - listed)}, only in the list {sorted(listed - set(derived))}")
        sys.exit(1)
    print(f"matches fc_common.EXCLUDED_SAMPLES ({len(listed)} samples), spec {fcc.spec_hash()}")
    if args.write_rates:
        path = ROOT / "data" / "empiricalStats" / fcc.MISCALL_RATES_FILE
        write_rates(fits, path)
        print(f"wrote {sum(len(fits[y]['names']) for y in YEARS)} miscall rates to {path}")


if __name__ == "__main__":
    main()
