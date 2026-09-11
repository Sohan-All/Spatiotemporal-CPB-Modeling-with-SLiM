import numpy as np
import pandas as pd
import random
import matplotlib.pyplot as plt
from pathlib import Path

import DataWrappers
import CollectData
import GenerateSimulationParams
import GenerateClusterData
import AnalyzeTreeSeq
import scale_constants as sc

import subprocess
import math

import warnings
import platform
import sys

#WARNING: don't run this file in VSCode. Run it in the terminal instead.

# Fixed so cluster identity stays stable across ABC iterations.
KMEANS_SEED = 42

def main(num_clusters, migration_rates_modifier, population_modifier, total_migration=0.05,
         mutation_rate=None, recombination_rate=None, ancestral_Ne=sc.ANCESTRAL_NE, silent=False):
    # No default for mutation_rate: it sets the diversity scale, so fail loudly rather than
    # silently reproduce a wrong one (CLAUDE.md 10.1).
    if mutation_rate is None or recombination_rate is None:
        raise ValueError(
            "main() requires explicit mutation_rate AND recombination_rate -- there are "
            "deliberately no defaults, because they set the diversity and linkage scales and the "
            "output files record no scale (CLAUDE.md 10.1). recombination_rate carried a "
            "hardcoded 2.75e-6 until 2026-09-09, which 10.1 claimed had been removed and had not "
            "(6.8.1). Pass ABCAnalysisNoRedis.DEFAULT_MUTATION_RATE / "
            "DEFAULT_RECOMBINATION_RATE, or scale_constants directly.")
    #for cleanliness
    warnings.filterwarnings("ignore")
    
    #Query for mutation rate
    #mutation_rate = 2.1e-9 #float(input("Enter the mutation rate (default 1e-7): ").strip() or 2.1e-9)

    #Query for recombination rate
    #recombination_rate = 2.75e-6 #float(input("Enter the recombination rate (default 1e-8): ").strip() or 2.75e-6)
        
    #Start by reading the data from final_data_for_modeling.csv
    if not silent:
        print("Setting up data for simulations...")
    field_data = CollectData.read_csv(Path('../data/final_data_for_modeling.csv'))
    
    #Cluster the coordinates using KMeans
    GenerateClusterData.cluster_coordinates(field_data, n_clusters=num_clusters, iters=2000, random_state=KMEANS_SEED)
    
    #Put the data for clusters into a list of Cluster objects
    clusters = GenerateClusterData.populate_cluster_objects(field_data, estimate_data=True)
    
    #Assign genomes to clusters based on the specifier matrix for a certain year
    GenerateClusterData.assign_genomes_to_clusters(clusters)
    
    #Generate a distance matrix for the clusters
    distances = GenerateClusterData.create_cluster_distance_matrix(clusters, output_path=Path('../data/cluster_distances.csv'))   
     
    #Save the cluster data to a CSV file
    GenerateClusterData.cluster_data_to_csv(clusters, output_path=Path('../data/cluster_data.csv'))
        
    #Generate migration rates based on the cluster distance matrix
    GenerateSimulationParams.determine_migration_rates(distances, total_migration=total_migration, scale=migration_rates_modifier, output_path=Path('../data/migration_rates.csv'))

    # THE REAL CEILING ON DEME COUNT, and it is NOT recapitation (CLAUDE.md 7.9.5).
    # CPBSampleSim*.slim line 30 does
    #     sim.addSubpop("p"+i, asInteger(Average Count[i] * POPMULT / numSubpops))
    # and SLiM refuses an empty subpopulation:
    #     ERROR (Population::AddSubpopulation): subpopulation p38 empty.
    # asInteger TRUNCATES, so any deme whose Average Count * POPMULT / numSubpops lands under 1.0
    # kills the run -- after the forward sim has already been set up, with an error that names a
    # SLiM subpop id and nothing about clusters. Measured: numClusters=200 at POPMULT=500 puts
    # 2 of 200 demes under 1.0 and dies; the same 200 demes at POPMULT=2000 is fine.
    # Note this couples numClusters and POPMULT: deme size falls as 1/numSubpops, so raising the
    # deme count SHRINKS every deme unless POPMULT rises with it.
    # Checked HERE, before SLiM, because the message SLiM gives is unactionable.
    _ac = [c.data[0] for c in clusters]   # Average Count, as cluster_data_to_csv writes it
    _min_size = min(_ac) * population_modifier / len(clusters)
    if _min_size < 1.0:
        _needed = math.ceil(len(clusters) / min(_ac))
        raise ValueError(
            f"numClusters={len(clusters)} with POPMULT={population_modifier} gives a smallest "
            f"deme of {_min_size:.3f} individuals, which SLiM rejects as an empty subpopulation. "
            f"Deme size is Average Count * POPMULT / numSubpops, so it falls as 1/numSubpops: "
            f"this cluster layout needs POPMULT >= {_needed} at {len(clusters)} clusters. "
            f"Raise POPMULT, or lower numClusters. See CLAUDE.md 7.9.5.")
    
    #Run the SLiM simulation to create the tree sequence
    if not silent:
        print("Running SLiM simulation...")
    
    #Run the appropriate SLiM script based on the operating system, passing in the population modifier as a parameter
    if platform.system() == "Windows":
        slim_script = Path('../SLiM_Code/CPBSampleSimWin.slim')
    elif platform.system() == "Linux":
        slim_script = Path('../SLiM_Code/CPBSampleSimLinux.slim')
    else:
        raise OSError(f"Unsupported operating system: {platform.system()}")

    subprocess.run(['slim', '-l', '0',
                    '-d', f'POPMULT={population_modifier}',
                    '-d', f'RECOMB={recombination_rate!r}',
                    str(slim_script)], check=True)
    
    #Does recapitation and mutation addition, then gets diversity and divergence statistics
    if not silent:
        print("Recapitating tree sequence...")
    AnalyzeTreeSeq.analyze_tree_sequence(mutation_rate=mutation_rate, recombination_rate=recombination_rate, ancestral_Ne=ancestral_Ne)
    if not silent:
        print("Successfully generated diversity and divergence statistics from tree sequence.")
    
    
    

if __name__ == "__main__":
    #Ask for input on number of clusters
    num_clusters = int(input("Enter the number of clusters (default 99): ").strip() or 99)
    
    #Query for a migration rates modifier
    migration_rates_modifier = float(input("Enter the migration rates modifier (default 0.0001): ").strip() or 0.0001)
    
    #Query for population modifier
    population_modifier = float(input("Enter the total population size (default 10000): ").strip() or 10000)
    
    # Imported here, not at module scope, so the simulation path has no import-time dependency
    # on the ABC driver.
    from ABCAnalysisNoRedis import DEFAULT_MUTATION_RATE, DEFAULT_RECOMBINATION_RATE

    #Query for the mutation rate. Not a biological rate -- report theta=4*Ne*mu (CLAUDE.md 6.1).
    mutation_rate = float(input(
        f"Enter the mutation rate (default {DEFAULT_MUTATION_RATE:g}, calibrated): ").strip()
        or DEFAULT_MUTATION_RATE)

    main(num_clusters, migration_rates_modifier, population_modifier,
         mutation_rate=mutation_rate, recombination_rate=DEFAULT_RECOMBINATION_RATE)
