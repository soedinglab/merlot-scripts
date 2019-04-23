"""
A script that will perform local averaging for a dimensionality-reduced
simulation.
"""

import argparse
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans

def cluster_kmeans(x, n_clust):
    kmeans = KMeans(n_clusters=n_clust).fit(x)
    reduced_kmeans = kmeans.cluster_centers_
    return reduced_kmeans, kmeans.labels_

def main(full_loc, dim):
    dimensionality = np.array(pd.read_csv(full_loc + ".csv", sep=" ", header=0, index_col=0))
    full_coords = dimensionality[:, 0:int(dim)]
    del dimensionality
    n_clust = int(4 * np.sqrt(full_coords.shape[0]))
    reduced_kmeans, kmeans_labels = cluster_kmeans(full_coords, n_clust)
    np.savetxt(full_loc + "_knn.txt", reduced_kmeans)

if __name__ == "__main__":
    PARSER = argparse.ArgumentParser(description='Perform k-means clustering \
                                     on a diffusion map.')
    PARSER.add_argument("-i", "--input", dest="input", help="Input file of \
                        relevant R object", metavar="FILE")
    PARSER.add_argument("-d", "--dim", dest="dim",
                        help="How many dimensions to keep",
                        metavar="INT")

    args = PARSER.parse_args()

    main(args.input, args.dim)
