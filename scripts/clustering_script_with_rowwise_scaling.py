"""
Cluster gene-expression profiles and spatial samples.

Input:
    Excel table with genes in rows and samples/spatial regions in columns.
    The first column contains gene names.

Outputs:
    genes_cluster_1.txt
    genes_cluster_2.txt

Plots:
    Row/column clustered heatmap
    Sample dendrogram
    Gene dendrogram
    Gene-clustering silhouette analysis
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.cluster.hierarchy import dendrogram, fcluster, linkage
from scipy.stats import zscore
from sklearn.metrics import silhouette_score


# =============================================================================
# SETTINGS
# =============================================================================

INPUT_FILE = Path("../processed_data_combined_L127L128.xlsx")
SHEET_NAME = 0
INDEX_COLUMN = 0

# Set True only if samples are rows and genes are columns in the Excel file.
TRANSPOSE_INPUT = False

N_GENE_CLUSTERS = 2
LINKAGE_METHOD = "ward"
HEATMAP_CMAP = "viridis"
OUTPUT_DIR = Path(".")
SILHOUETTE_MAX_CLUSTERS = 14

# Colours used for the gene-clustering dendrogram.
GENE_CLUSTER_COLOURS = {1: "green", 2: "red"}


# =============================================================================
# DATA PREPARATION
# =============================================================================

def load_expression_data(path, sheet_name=0, index_col=0, transpose=False):
    """Load numeric expression data with genes as rows and samples as columns."""

    data = pd.read_excel(
        path,
        sheet_name=sheet_name,
        index_col=index_col,
    )

    if transpose:
        data = data.T

    data = data.apply(pd.to_numeric, errors="coerce")
    data = data.replace([np.inf, -np.inf], np.nan)

    # Remove entirely empty rows and columns.
    data = data.dropna(axis=0, how="all").dropna(axis=1, how="all")
    data = data.loc[data.notna().any(axis=1)]

    if data.index.has_duplicates:
        raise ValueError("Gene names must be unique.")

    if data.columns.has_duplicates:
        raise ValueError("Sample/region names must be unique.")

    # Replace missing values with the mean for the corresponding gene.
    data = data.T.fillna(data.mean(axis=1)).T

    # Remove genes with zero variance.
    data = data.loc[data.var(axis=1, ddof=0) > 0]

    if data.shape[0] < 2:
        raise ValueError("At least two variable genes are required.")

    if data.shape[1] < 2:
        raise ValueError("At least two samples/regions are required.")

    return data


def zscore_rows(expression):
    """Z-score each gene across samples/regions."""

    scaled = expression.apply(
        lambda row: zscore(
            row.to_numpy(dtype=float),
            nan_policy="omit",
        ),
        axis=1,
        result_type="expand",
    )

    scaled.index = expression.index
    scaled.columns = expression.columns

    return scaled.replace([np.inf, -np.inf], np.nan).fillna(0.0)


# =============================================================================
# PLOTTING
# =============================================================================

def plot_clustered_heatmap(expression_scaled):
    """Plot a heatmap with hierarchical clustering of genes and samples."""

    row_linkage = linkage(
        expression_scaled.to_numpy(),
        method=LINKAGE_METHOD,
    )

    col_linkage = linkage(
        expression_scaled.T.to_numpy(),
        method=LINKAGE_METHOD,
    )

    grid = sns.clustermap(
        expression_scaled,
        row_linkage=row_linkage,
        col_linkage=col_linkage,
        cmap=HEATMAP_CMAP,
        figsize=(10, 10),
        dendrogram_ratio=(0.18, 0.18),
        tree_kws={"linewidths": 1.2},
    )

    grid.fig.suptitle(
        "Gene-expression heatmap (row z-scores)",
        y=1.02,
    )

    plt.show()


def plot_sample_dendrogram(linkage_matrix, sample_names):
    """Plot hierarchical clustering of spatial samples/regions."""

    plt.figure(figsize=(10, 5))

    dendrogram(
        linkage_matrix,
        labels=sample_names,
        leaf_rotation=90,
    )

    plt.title("Sample clustering (gene-wise z-scored expression)")
    plt.xlabel("Sample / spatial region")
    plt.ylabel("Ward linkage distance")
    plt.tight_layout()

    plt.show()


def make_gene_dendrogram_colour_func(linkage_matrix, leaf_cluster_labels):
    """Colour gene dendrogram branches according to cluster membership."""

    n_leaves = len(leaf_cluster_labels)

    descendants = {
        leaf: {int(leaf_cluster_labels[leaf])}
        for leaf in range(n_leaves)
    }

    colours = {}

    for row_idx, row in enumerate(linkage_matrix):
        left, right = int(row[0]), int(row[1])
        node_id = n_leaves + row_idx

        descendants[node_id] = descendants[left] | descendants[right]
        node_clusters = descendants[node_id]

        if len(node_clusters) == 1:
            cluster_id = next(iter(node_clusters))
            colours[node_id] = GENE_CLUSTER_COLOURS.get(
                cluster_id,
                "grey",
            )
        else:
            colours[node_id] = "black"

    return lambda node_id: colours.get(node_id, "black")


def plot_gene_dendrogram(linkage_matrix, genes, cluster_labels):
    """Plot the gene dendrogram, coloured by gene-cluster membership."""

    colour_func = make_gene_dendrogram_colour_func(
        linkage_matrix,
        cluster_labels,
    )

    plt.figure(figsize=(max(12, len(genes) * 0.12), 6))

    dendrogram(
        linkage_matrix,
        labels=list(genes),
        leaf_rotation=90,
        leaf_font_size=6,
        color_threshold=0,
        link_color_func=colour_func,
    )

    plt.title("Gene clustering (row z-scores)")
    plt.xlabel("Gene")
    plt.ylabel("Ward linkage distance")
    plt.tight_layout()

    plt.show()


def plot_silhouette_analysis(expression_scaled, gene_linkage):
    """Plot silhouette scores for alternative gene-cluster numbers."""

    n_genes = expression_scaled.shape[0]
    max_k = min(SILHOUETTE_MAX_CLUSTERS, n_genes - 1)

    if max_k < 2:
        print("Silhouette analysis skipped: too few genes.")
        return

    candidate_ks = []
    scores = []

    for k in range(2, max_k + 1):
        labels = fcluster(
            gene_linkage,
            k,
            criterion="maxclust",
        )

        if len(np.unique(labels)) < 2 or len(np.unique(labels)) >= n_genes:
            continue

        score = silhouette_score(
            expression_scaled.to_numpy(),
            labels,
        )

        candidate_ks.append(k)
        scores.append(score)

    if not scores:
        print("Silhouette analysis skipped: no valid cluster solutions.")
        return

    plt.figure(figsize=(7, 4))
    plt.plot(candidate_ks, scores, marker="o")
    plt.xlabel("Requested number of clusters")
    plt.ylabel("Silhouette score")
    plt.title("Silhouette analysis (row z-scores)")
    plt.tight_layout()

    plt.show()


# =============================================================================
# MAIN ANALYSIS
# =============================================================================

def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    expression = load_expression_data(
        INPUT_FILE,
        sheet_name=SHEET_NAME,
        index_col=INDEX_COLUMN,
        transpose=TRANSPOSE_INPUT,
    )

    print(
        f"Loaded {expression.shape[0]} genes across "
        f"{expression.shape[1]} samples/regions."
    )

    # Row-wise z-scores are used for clustering.
    expression_scaled = zscore_rows(expression)

    # Plot 1: clustered heatmap.
    plot_clustered_heatmap(expression_scaled)

    # Cluster spatial samples/regions and plot their dendrogram.
    sample_linkage = linkage(
        expression_scaled.T.to_numpy(),
        method=LINKAGE_METHOD,
        metric="euclidean",
    )

    plot_sample_dendrogram(
        sample_linkage,
        expression_scaled.columns.tolist(),
    )

    # Hierarchically cluster genes.
    gene_linkage = linkage(
        expression_scaled.to_numpy(),
        method=LINKAGE_METHOD,
        metric="euclidean",
    )

    gene_cluster_labels = fcluster(
        gene_linkage,
        N_GENE_CLUSTERS,
        criterion="maxclust",
    )

    # Save only the gene lists required for downstream PANTHER analysis.
    gene_clusters = {}

    for gene, cluster_id in zip(
        expression_scaled.index,
        gene_cluster_labels,
    ):
        gene_clusters.setdefault(int(cluster_id), []).append(str(gene))

    for cluster_id, genes in sorted(gene_clusters.items()):
        output_file = OUTPUT_DIR / f"genes_cluster_{cluster_id}.txt"

        with output_file.open("w") as f:
            f.write("\n".join(genes) + "\n")

        print(
            f"Saved {len(genes)} genes for cluster {cluster_id}: "
            f"{output_file}"
        )

    # Plot 3: gene dendrogram coloured by cluster.
    plot_gene_dendrogram(
        gene_linkage,
        expression_scaled.index,
        gene_cluster_labels,
    )

    # Plot 4: silhouette analysis of alternative gene-cluster numbers.
    plot_silhouette_analysis(
        expression_scaled,
        gene_linkage,
    )


if __name__ == "__main__":
    main()
