"""How close was each year's sampling to the previous collection's? (TODO 2026-09-29 step 3)

The professor says collectors deliberately sampled the field nearest the one sampled before, so local
populations should be persistent even though the fields rotate (CLAUDE.md 7.9.14C). This measures
it from the specifier matrices (headerless; col0 site, col1 lat, col2 lon -- CLAUDE.md 4):

  - for every 2019 site, distance to the nearest 2015 site; for every 2023 site, to the nearest 2019
    site and to the nearest 2015 site;
  - for context, each site's nearest SAME-year neighbour, which is what "close" has to beat;
  - which sites the fc_common matched pairs use (0.5 km rule), so the rest are visible as new.

Great-circle (haversine) distance. Runs in under a second.

Usage (from diagnostics/):  python resample_distances.py
"""
import csv
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
YEARS = ["2015", "2019", "2023"]
R_EARTH_KM = 6371.0088
EDGES_KM = [0.5, 1, 2, 5, 10, 25]


def read_spec(year):
    rows = [r for r in csv.reader(open(ROOT / "data" / "Genetic_Data" / f"specifier_matrix_{year}.csv",
                                       newline="", encoding="utf-8")) if r]
    return [r[0].strip() for r in rows], np.array([[float(r[1]), float(r[2])] for r in rows])


def hav(a, b):
    """km between every row of a and every row of b (lat, lon in degrees)."""
    la1, lo1 = np.radians(a[:, 0])[:, None], np.radians(a[:, 1])[:, None]
    la2, lo2 = np.radians(b[:, 0])[None, :], np.radians(b[:, 1])[None, :]
    h = np.sin((la2 - la1) / 2) ** 2 + np.cos(la1) * np.cos(la2) * np.sin((lo2 - lo1) / 2) ** 2
    return 2 * R_EARTH_KM * np.arcsin(np.sqrt(h))


def summarise(label, d):
    cum = "  ".join(f"<={e:g}km:{int(np.sum(d <= e)):2d}" for e in EDGES_KM)
    print(f"  {label:30s} n={len(d):2d}  median {np.median(d):6.2f} km   {cum}")


def main():
    spec = {y: read_spec(y) for y in YEARS}

    print("Nearest SAME-year neighbour (the scale 'close' has to beat):")
    for y in YEARS:
        D = hav(spec[y][1], spec[y][1])
        np.fill_diagonal(D, np.inf)
        summarise(y, D.min(1))

    print("\nNearest site in an EARLIER collection:")
    per_site = []
    for later, earlier in [("2019", "2015"), ("2023", "2019"), ("2023", "2015")]:
        D = hav(spec[later][1], spec[earlier][1])
        j = D.argmin(1)
        dmin = D.min(1)
        summarise(f"{later} -> nearest {earlier}", dmin)
        for i, name in enumerate(spec[later][0]):
            per_site.append((later, earlier, name, spec[earlier][0][j[i]], dmin[i]))

    print("\nPer site (later -> earlier), sorted by distance:")
    for later, earlier, a, b, dist in sorted(per_site, key=lambda r: (r[0], r[1], r[4])):
        print(f"  {later}->{earlier}  {a:26s} nearest {b:26s} {dist:7.2f} km")


if __name__ == "__main__":
    main()
