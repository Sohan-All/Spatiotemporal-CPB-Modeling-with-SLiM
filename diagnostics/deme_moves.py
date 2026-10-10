"""How often does production's site->deme matching split a persistent local population? (TODO 2026-09-29 step 4)

Production matches each year's sites to demes one-to-one by the Hungarian algorithm, each year
independently (GenerateClusterData.assign_genomes_to_clusters_idv_year). Sampling deliberately
followed the population to the nearest field (CLAUDE.md 7.9.14C/E3: "same population" = within
2 km across collections), so a pair of sites that are the same local population in two years should
land in the same deme. This counts how often they do not, at each numClusters the prior can draw,
and compares the nearest-centroid rule temporal F_c already uses (AnalyzeTreeSeq._nearest_cluster_rows).

Rebuilds the clusters exactly as Main.main does (KMeans, seed Main.KMEANS_SEED) but WRITES NOTHING
to data/ -- it stops before create_cluster_distance_matrix / the cluster_data.csv write.

Usage (from diagnostics/):  conda run -n cpb-env python deme_moves.py [--same-km 2]
"""
import argparse
import os
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "Python_Code"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
os.chdir(ROOT / "Python_Code")          # the pipeline's paths are relative to Python_Code/

import CollectData                       # noqa: E402
import GenerateClusterData as GCD        # noqa: E402
from Main import KMEANS_SEED             # noqa: E402
from resample_distances import read_spec, hav, YEARS   # noqa: E402

YEAR_IDX = {"2015": 0, "2019": 1, "2023": 2}
NUM_CLUSTERS = [33, 66, 99]             # the prior's {1,2,3} x 33


def build(nc):
    fields = CollectData.read_csv(Path("../data/final_data_for_modeling.csv"))
    GCD.cluster_coordinates(fields, n_clusters=nc, iters=2000, random_state=KMEANS_SEED)
    clusters = GCD.populate_cluster_objects(fields, estimate_data=True)
    for y in YEARS:
        GCD.assign_genomes_to_clusters_idv_year(
            clusters, int(y), specifier_matrix=Path(f"../data/Genetic_Data/specifier_matrix_{y}.csv"))
    cent = np.array([[c.latitude, c.longitude] for c in clusters])
    hung = {}
    for y in YEARS:
        n = len(read_spec(y)[0])
        d = np.full(n, -1)
        for j, c in enumerate(clusters):
            i = c.genome_assignments[YEAR_IDX[y]]
            if i is not None and i == i:          # skip None / NaN
                d[int(i)] = j
        assert (d >= 0).all(), f"{y}: unassigned site"
        hung[y] = d
    return cent, hung


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--same-km", type=float, default=2.0)
    a = ap.parse_args()
    spec = {y: read_spec(y) for y in YEARS}
    pairs = [("2015", "2019"), ("2019", "2023"), ("2015", "2023")]

    for nc in NUM_CLUSTERS:
        cent, hung = build(nc)
        near = {y: hav(spec[y][1], cent).argmin(axis=1) for y in YEARS}
        off = {y: int((hung[y] != near[y]).sum()) for y in YEARS}
        print(f"\n=== numClusters {nc} ===")
        print("  sites not on their nearest deme (Hungarian): "
              + ", ".join(f"{y} {off[y]}/{len(hung[y])}" for y in YEARS))
        tot_h = tot_n = tot = 0
        for ya, yb in pairs:
            D = hav(spec[ya][1], spec[yb][1])
            ii, jj = np.nonzero(D <= a.same_km)
            sh = int((hung[ya][ii] == hung[yb][jj]).sum())
            sn = int((near[ya][ii] == near[yb][jj]).sum())
            tot += len(ii); tot_h += sh; tot_n += sn
            print(f"  {ya}->{yb}: {len(ii):2d} same-population pairs (<= {a.same_km:g} km) | "
                  f"same deme: Hungarian {sh}/{len(ii)}, nearest {sn}/{len(ii)}")
            for i, j in zip(ii, jj):
                if hung[ya][i] != hung[yb][j]:
                    sep = hav(cent[[hung[ya][i]]], cent[[hung[yb][j]]])[0, 0]
                    print(f"      split: {spec[ya][0][i]} (deme {hung[ya][i]}) / {spec[yb][0][j]} "
                          f"(deme {hung[yb][j]})  sites {D[i, j]:.2f} km apart, demes {sep:.1f} km")
        print(f"  TOTAL same deme: Hungarian {tot_h}/{tot} ({tot_h / tot:.0%}), "
              f"nearest {tot_n}/{tot} ({tot_n / tot:.0%})")


if __name__ == "__main__":
    main()
