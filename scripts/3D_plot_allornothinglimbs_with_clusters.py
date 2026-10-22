#!/usr/bin/env python3

import sys
import argparse

import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import AgglomerativeClustering
from vedo import *


# ============================================================
# 1. Load input data
# ============================================================

parser = argparse.ArgumentParser(
    description=(
        "Perform hierarchical clustering of 3D fluorescence "
        "point data, display the clusters and centroids, "
        "and generate the Limb80.mp4 movie."
    ),
    add_help=True
)

parser.add_argument(
    "input_file",
    nargs="?",
    help="Input text file containing exactly three columns of numerical coordinates."
)

args = parser.parse_args()

if args.input_file is None:
    print("Usage: python3 script.py <input_file>")
    sys.exit(1)

file = args.input_file

try:
    pointdata = np.loadtxt(file)
except (OSError, ValueError) as e:
    print(f"Error reading input file '{file}': {e}")
    sys.exit(1)

pointdata = np.atleast_2d(pointdata)

if pointdata.shape[1] != 3:
    raise ValueError("Input file must contain exactly three columns.")

if not np.all(np.isfinite(pointdata)):
    raise ValueError("Input data contains NaN or infinite values.")

totalpoints = len(pointdata)
numof_clusters_wanted = 8

if totalpoints < numof_clusters_wanted:
    raise ValueError("At least 8 points are required for clustering.")

print("Total number of points:", totalpoints)


# ============================================================
# 2. Set up the Matplotlib 3D plot
# ============================================================

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.set_xlabel('6FAM/TAMRA')
ax.set_ylabel('Alexa405/TAMRA')
ax.set_zlabel('Alexa647/TAMRA')


# ============================================================
# 3. Normalise the coordinate dimensions
# ============================================================

# Preserve the original maximum-based normalisation.

print("\n3D point data:")

max_values = np.max(pointdata, axis=0)

maxx, maxy, maxz = max_values

print("\nOriginal coordinate maxima:")
print("X:", maxx)
print("Y:", maxy)
print("Z:", maxz)

if np.any(max_values == 0):
    raise ValueError(
        "At least one coordinate maximum is zero; cannot normalise."
    )

pointdata = pointdata / max_values

print("\nFirst point, normalised Y coordinate:")
print(pointdata[0, 1])


# ============================================================
# 4. Perform hierarchical clustering
# ============================================================

# Ward linkage with Euclidean distance.
# Support both older and newer scikit-learn versions.

try:
    cluster = AgglomerativeClustering(
        n_clusters=numof_clusters_wanted,
        metric='euclidean',
        linkage='ward'
    )
except TypeError:
    cluster = AgglomerativeClustering(
        n_clusters=numof_clusters_wanted,
        affinity='euclidean',
        linkage='ward'
    )

cluster.fit_predict(pointdata)

print("\nClustering complete.")


# ============================================================
# 5. Calculate cluster centroids
# ============================================================

centroids = np.zeros((numof_clusters_wanted, 3))
numin_cluster = np.zeros(numof_clusters_wanted, dtype=int)

for i in range(numof_clusters_wanted):

    cluster_points = pointdata[cluster.labels_ == i]

    numin_cluster[i] = len(cluster_points)

    if numin_cluster[i] > 0:
        centroids[i] = np.mean(cluster_points, axis=0)

print("\nCluster centroids:")

for i in range(numof_clusters_wanted):
    print(f"\nCluster {i}:")
    print("Number of points:", numin_cluster[i])
    print("X:", centroids[i, 0])
    print("Y:", centroids[i, 1])
    print("Z:", centroids[i, 2])


# ============================================================
# 6. Plot points and centroids using Matplotlib
# ============================================================

# Plot the cluster centroids in black.
ax.scatter3D(
    centroids[:, 0],
    centroids[:, 1],
    centroids[:, 2],
    color='#000000'
)

# Plot the individual points, coloured by cluster.
ax.scatter3D(
    pointdata[:, 0],
    pointdata[:, 1],
    pointdata[:, 2],
    c=cluster.labels_,
    cmap='rainbow'
)

plt.show()





