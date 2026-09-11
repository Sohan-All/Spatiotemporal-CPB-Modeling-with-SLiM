from pathlib import Path
import numpy as np
import pandas as pd
import random
from math import pi, cos, asin, sqrt
from sklearn.cluster import KMeans
from scipy.optimize import linear_sum_assignment
import DataWrappers

# Cluster the coordinates using KMeans and add cluster labels to the field data    
def cluster_coordinates(field_data, n_clusters=5, iters=2000, random_state=random.randint(0, 1000)):
    # Prepare data for clustering
    coordinates = np.array([[field.latitude, field.longitude] for field in field_data.values()])
    
    # Perform KMeans clustering
    kmeans = KMeans(n_clusters=n_clusters, max_iter=iters, random_state=random_state, n_init='auto')
    kmeans.fit(coordinates)
    
    # Add cluster labels to the field data
    for i, fvid in enumerate(field_data.keys()):
        field_data[fvid].set_cluster(kmeans.labels_[i])
    
    return field_data


# This function is used to take field data that has been run through cluster_coordinates and create a new list
# that contains averaged data for each cluster such as location and count
def populate_cluster_objects(field_data, estimate_data=True):
    # Start by creating a list of lists of lists which sorts each node into its respective cluster. 
    # Each cluster has a list of latitudes, longitudes, average counts, and field_fvids which will be operated on to create generalized cluster data
    clusters = [None] * (1+max(field_data[fvid].cluster for fvid in field_data.keys()))  # Create a list of None with length equal to the number of clusters
    
    for fvid in field_data.keys():
        cluster_id = field_data[fvid].cluster
        
        #print(len(clusters), cluster_id)
        if clusters[cluster_id] is None:
            clusters[cluster_id] = DataWrappers.Cluster(cluster_id)
        clusters[cluster_id].add_field(field_data[fvid])
    
    # Calculate the yearly data and coordinates for each cluster based on inputted fields
    for cluster in clusters:
        cluster.calculate_coordinates()
        cluster.average_data()
            
    return clusters

def assign_genomes_to_clusters(clusters):
    """
    This function assigns specific groups of genomes to clusters based on how close the
    genetic sample was taken to the cluster. It only assigns each genome once.
    
    parameters:
    clusters (list): List of Cluster objects to which genomes will be assigned.
    genetic_data (str): Path to the genetic data CSV file.
    """
    assign_genomes_to_clusters_idv_year(clusters, 2023, specifier_matrix=Path("../data/Genetic_Data/specifier_matrix_2023.csv"))
    assign_genomes_to_clusters_idv_year(clusters, 2019, specifier_matrix=Path("../data/Genetic_Data/specifier_matrix_2019.csv"))
    assign_genomes_to_clusters_idv_year(clusters, 2015, specifier_matrix=Path("../data/Genetic_Data/specifier_matrix_2015.csv"))
    
    

def assign_genomes_to_clusters_idv_year(clusters, year, specifier_matrix=Path("../data/Genetic_Data/specifier_matrix_2023.csv"), report=False):
    """Assign each sequenced site of `year` to a distinct cluster, minimising TOTAL displacement.

    One site per cluster is deliberate: two sites sharing a deme would have IDENTICAL simulated
    statistics, which is worse for element-wise fitting than a small position error. What is NOT
    deliberate is how that constraint used to be enforced.

    WAS (until 2026-09-09): greedy in specifier-row order -- each site took its nearest cluster
    that no EARLIER site of that year had already claimed. Order-dependent, and badly so
    (CLAUDE.md 7.9.2, measured at 33 clusters):

        year   sites off their nearest cluster   median site->deme   optimal
        2015           16/24  (67%)                  10.74 km        7.41 km
        2019           12/17  (71%)                  11.31 km        8.75 km
        2023           12/20  (60%)                   8.01 km        7.48 km

    Worst single case: Paramount652-2015 was represented by a deme 26.7 km away when a cluster
    sat 4.8 km from it. That is injected noise in exactly the geography the dispersal kernel and
    the IBD regression are trying to read, and it was also year-dependent -- the SAME FIELD landed
    in different clusters in different years (Arlington: cluster 10 in 2015, 21 in 2019), so the
    model had no notion of a persistent field.

    NOW: the optimal one-to-one matching (Hungarian algorithm, scipy.optimize). Same constraint,
    globally minimal total displacement, and no dependence on row order. Cuts total displacement
    12-18% at 33 clusters, and far more as the cluster count rises, since a good matching exists
    once clusters outnumber sites comfortably.

    NOTE this does NOT by itself make a field persistent across years -- each year is still
    matched independently, and a field can still move if its neighbours differ between years.
    Anything comparing the same field across years must still pick its own deme (see
    diagnostics/temporal_fc.py).

    report=True prints the per-year displacement summary.
    """
    # Read the genetic data file. The specifier CSVs do not include a header row,
    # so use header=None to avoid pandas treating the first sample as column names
    specifier = pd.read_csv(specifier_matrix, header=None)

    # Extract latitude and longitude pairs from the specifier matrix
    genome_coords = list(zip(specifier.iloc[:, 1], specifier.iloc[:, 2]))

    yearIdx = {2015: 0, 2019: 1, 2023: 2}[year]

    if len(clusters) < len(genome_coords):
        raise ValueError(
            f"{year}: {len(genome_coords)} sequenced sites but only {len(clusters)} clusters. "
            f"One site per cluster is required, so numClusters must be >= the largest year "
            f"(24, in 2015).")

    # cost[i][j] = metres from site i to cluster j. distance() returns int metres.
    cost = np.array([[distance(lat, lon, c.latitude, c.longitude) for c in clusters]
                     for lat, lon in genome_coords], dtype=float)
    rows, cols = linear_sum_assignment(cost)
    for i, j in zip(rows, cols):
        clusters[j].genome_assignments[yearIdx] = int(i)

    if report:
        d = np.sort(cost[rows, cols]) / 1000.0
        nearest = cost.min(axis=1) / 1000.0
        print(f"  {year}: {len(genome_coords)} sites -> {len(clusters)} clusters | "
              f"site->deme median {np.median(d):.2f} km, max {d.max():.2f} km "
              f"(nearest-cluster floor: median {np.median(nearest):.2f} km) | "
              f"{int((cost[rows, cols] > cost.min(axis=1)).sum())} not on their nearest")
        
    
    
    
    

# Function to calculate distance between two coordinates in km using Haversine formula
def distance(lat1, lon1, lat2, lon2):
    r = 6371  # km
    p = pi / 180

    a = 0.5 - cos((lat2 - lat1) * p) / 2 + cos(lat1 * p) * cos(lat2 * p) * (1 - cos((lon2 - lon1) * p)) / 2
    return int(2 * r * asin(sqrt(a)) * 1000)  # Convert km to meters for more precision



# This function takes a list of clusters and calculates the distance between each pair of clusters, storing the results in a distance matrix 
#The csv file is outputted to the specified path and returned by the function.
def create_cluster_distance_matrix(clusters, output_path=Path("../data/cluster_distances.csv")):

    # Create a distance matrix with cluster IDs as both row and column headers
    grid = np.eye(len(clusters) + 1, dtype=int)
    for i in range(len(clusters)):
        grid[i + 1][0] = clusters[i].cluster_id 
        grid[0][i + 1] = clusters[i].cluster_id 
    grid[0][0] = 0

    # Calculate distances and fill the distance matrix storing data only if dist is <= cutoff
    for i in range(len(clusters)):
        for j in range(i, len(clusters)):
            lat1 = clusters[i].latitude
            lon1 = clusters[i].longitude
            lat2 = clusters[j].latitude
            lon2 = clusters[j].longitude
        
            dist = distance(lat1, lon1, lat2, lon2)  # Convert km to meters for more precision
            grid[i + 1][j + 1] = dist
            grid[j + 1][i + 1] = dist  # Ensure symmetry in the distance matrix

    np.savetxt(output_path, grid, delimiter=",", fmt='%i')
    return grid


# This function takes a list of cluster objects and saves their data to a CSV file
def cluster_data_to_csv(clusters, output_path=Path('../data/cluster_data.csv')):
    # Create a DataFrame to hold the cluster data
    data = {
        'Cluster ID': [],
        'Latitude': [],
        'Longitude': [],
        'Average Count': [],
        'Average GDD': [],
        'Genome Assignment 2015': [],
        'Genome Assignment 2019': [],
        'Genome Assignment 2023': []
    }
    
    for cluster in clusters:
        data['Cluster ID'].append(cluster.cluster_id)
        data['Latitude'].append(cluster.latitude)
        data['Longitude'].append(cluster.longitude)
        data['Average Count'].append(cluster.data[0])  # Assuming first element is average count
        data['Average GDD'].append(cluster.data[1])  # Assuming second element is average GDD
        if cluster.genome_assignments[0] is not None:
            data['Genome Assignment 2015'].append(cluster.genome_assignments[0])
        else:
            data['Genome Assignment 2015'].append('')
        if cluster.genome_assignments[1] is not None:
            data['Genome Assignment 2019'].append(cluster.genome_assignments[1])
        else:
            data['Genome Assignment 2019'].append('')
        if cluster.genome_assignments[2] is not None:
            data['Genome Assignment 2023'].append(cluster.genome_assignments[2])
        else:
            data['Genome Assignment 2023'].append('')
    
    df = pd.DataFrame(data)
    df.to_csv(output_path, index=False)
    
    
# assign_genomes_to_clusters(None, specifier_matrix="../data/Genetic_Data/specifier_matrix_2023.csv")