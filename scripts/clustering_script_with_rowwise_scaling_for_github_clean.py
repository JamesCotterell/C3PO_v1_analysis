"""
Cluster gene-expression profiles and spatial samples.

Expected input:
    Excel table with genes in rows and samples/spatial regions in columns.
    The first column contains gene names (index_column=0).

Outputs:
    distance_matrix.csv
    marker_genes_cluster_<ID>.csv
    gene_clusters.csv
    gene_module_scores.csv
    gene_module_enrichment.csv

Plots:
    row/column clustered heatmap, sample dendrogram, gene dendrogram,
    gene-clustering silhouette analysis
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.cluster.hierarchy import dendrogram, fcluster, linkage
from scipy.spatial.distance import pdist, squareform
from scipy.stats import ttest_ind, zscore
from sklearn.metrics import silhouette_score
from statsmodels.stats.multitest import multipletests


# =============================================================================
# SETTINGS
# =============================================================================

INPUT_FILE = Path(
    "./processed_data_combined_L127L128.xlsx"
)
SHEET_NAME = 0
INDEX_COLUMN = 0

# Input is expected to have genes as rows and samples/regions as columns.
# Set True only if the Excel table instead has samples as rows and genes as columns.
TRANSPOSE_INPUT = False

N_SAMPLE_CLUSTERS = 2
N_GENE_CLUSTERS = 2
LINKAGE_METHOD = "ward"
HEATMAP_CMAP = "viridis"
OUTPUT_DIR = Path(".")
SILHOUETTE_MAX_CLUSTERS = 14

# Explicit colours for the two-cluster dendrogram/plots.
GENE_CLUSTER_COLOURS = {1: "green", 2: "red"}


# =============================================================================
# DATA PREPARATION
# =============================================================================

def load_expression_data(path, sheet_name=0, index_col=0, transpose=False):
    """Load numeric expression data with genes as rows and samples as columns."""
    data = pd.read_excel(path, sheet_name=sheet_name, index_col=index_col)

    if transpose:
        data = data.T

    data = data.apply(pd.to_numeric, errors="coerce")
    data = data.replace([np.inf, -np.inf], np.nan)

    # Drop entirely empty rows/columns, then remove genes with no usable values.
    data = data.dropna(axis=0, how="all").dropna(axis=1, how="all")
    data = data.loc[data.notna().any(axis=1)]

    if data.index.has_duplicates:
        raise ValueError("Gene names (row index) must be unique.")
    if data.columns.has_duplicates:
        raise ValueError("Sample/region names (column labels) must be unique.")

    # Fill remaining missing values per gene with that gene's mean. This avoids
    # silently converting missing measurements to zero.
    data = data.T.fillna(data.mean(axis=1)).T

    # Remove genes that still have zero variance or invalid values.
    data = data.loc[data.var(axis=1, ddof=0) > 0]
    if data.shape[0] < 2:
        raise ValueError("At least two variable genes are required.")
    if data.shape[1] < 2:
        raise ValueError("At least two samples/regions are required.")

    return data


def zscore_rows(expression):
    """Z-score each gene across samples/regions."""
    scaled = expression.apply(
        lambda row: zscore(row.to_numpy(dtype=float), nan_policy="omit"),
        axis=1,
        result_type="expand",
    )
    scaled.index = expression.index
    scaled.columns = expression.columns
    return scaled.replace([np.inf, -np.inf], np.nan).fillna(0.0)


# =============================================================================
# CLUSTERING AND STATISTICS
# =============================================================================

def cluster_genes_and_score_modules(expression_scaled, n_clusters):
    """Cluster genes and calculate each module's mean scaled profile."""
    linkage_matrix = linkage(
        expression_scaled.to_numpy(), method=LINKAGE_METHOD, metric="euclidean"
    )
    labels = fcluster(linkage_matrix, n_clusters, criterion="maxclust")
    gene_to_cluster = dict(zip(expression_scaled.index, labels))

    module_scores = expression_scaled.groupby(labels, sort=True).mean()
    module_scores.index.name = "gene_cluster"

    return module_scores, gene_to_cluster, linkage_matrix, labels


def find_enriched_modules(module_scores, sample_cluster_map):
    """Welch-test each module's score in each sample cluster versus all others."""
    scores_by_sample = module_scores.T
    scores_by_sample["sample_cluster"] = scores_by_sample.index.map(sample_cluster_map)

    results = []
    module_ids = module_scores.index.tolist()
    cluster_ids = sorted(scores_by_sample["sample_cluster"].dropna().unique())

    for cluster_id in cluster_ids:
        in_cluster = scores_by_sample[
            scores_by_sample["sample_cluster"] == cluster_id
        ]
        out_cluster = scores_by_sample[
            scores_by_sample["sample_cluster"] != cluster_id
        ]

        for module_id in module_ids:
            inside = in_cluster[module_id].dropna().to_numpy(dtype=float)
            outside = out_cluster[module_id].dropna().to_numpy(dtype=float)

            if len(inside) < 2 or len(outside) < 2:
                p_value = np.nan
            else:
                _, p_value = ttest_ind(inside, outside, equal_var=False)

            results.append(
                {
                    "sample_cluster": cluster_id,
                    "gene_module": module_id,
                    "mean_in_cluster": np.mean(inside) if len(inside) else np.nan,
                    "mean_outside_cluster": np.mean(outside) if len(outside) else np.nan,
                    "mean_difference": (
                        np.mean(inside) - np.mean(outside)
                        if len(inside) and len(outside)
                        else np.nan
                    ),
                    "p_value": p_value,
                }
            )

    result_df = pd.DataFrame(results)
    if not result_df.empty:
        result_df["p_value_adj"] = np.nan
        valid = result_df["p_value"].notna()
        if valid.any():
            result_df.loc[valid, "p_value_adj"] = multipletests(
                result_df.loc[valid, "p_value"], method="fdr_bh"
            )[1]
        result_df = result_df.sort_values(
            ["sample_cluster", "p_value_adj", "mean_difference"],
            ascending=[True, True, False],
            na_position="last",
        )

    return result_df


def find_marker_genes(expression, sample_cluster_map):
    """
    Compare each sample cluster with all remaining samples for every gene.

    `mean_difference` is mean expression in the cluster minus mean expression
    outside it. It is not a log2 fold change.
    """
    results = {}

    for cluster_id in sorted(set(sample_cluster_map.values())):
        inside_samples = [
            sample for sample in expression.columns
            if sample_cluster_map.get(sample) == cluster_id
        ]
        outside_samples = [
            sample for sample in expression.columns
            if sample_cluster_map.get(sample) != cluster_id
        ]

        rows = []
        for gene in expression.index:
            inside = expression.loc[gene, inside_samples].dropna().to_numpy(float)
            outside = expression.loc[gene, outside_samples].dropna().to_numpy(float)

            if len(inside) < 2 or len(outside) < 2:
                p_value = np.nan
            elif np.var(inside) == 0 and np.var(outside) == 0:
                p_value = np.nan
            else:
                _, p_value = ttest_ind(inside, outside, equal_var=False)

            rows.append(
                {
                    "gene": gene,
                    "p_value": p_value,
                    "mean_difference": (
                        np.mean(inside) - np.mean(outside)
                        if len(inside) and len(outside)
                        else np.nan
                    ),
                    "mean_in_cluster": np.mean(inside) if len(inside) else np.nan,
                    "mean_outside_cluster": (
                        np.mean(outside) if len(outside) else np.nan
                    ),
                }
            )

        marker_df = pd.DataFrame(rows)
        marker_df["p_value_adj"] = np.nan
        valid = marker_df["p_value"].notna()
        if valid.any():
            marker_df.loc[valid, "p_value_adj"] = multipletests(
                marker_df.loc[valid, "p_value"], method="fdr_bh"
            )[1]

        results[cluster_id] = marker_df.sort_values(
            ["mean_difference", "p_value_adj"],
            ascending=[False, True],
            na_position="last",
        ).reset_index(drop=True)

    return results


# =============================================================================
# PLOTTING HELPERS
# =============================================================================

def plot_clustered_heatmap(expression_scaled):
    """Plot a heatmap with hierarchical clustering on genes and samples."""
    row_linkage = linkage(expression_scaled.to_numpy(), method=LINKAGE_METHOD)
    col_linkage = linkage(expression_scaled.T.to_numpy(), method=LINKAGE_METHOD)

    grid = sns.clustermap(
        expression_scaled,
        row_linkage=row_linkage,
        col_linkage=col_linkage,
        cmap=HEATMAP_CMAP,
        figsize=(10, 10),
        dendrogram_ratio=(0.18, 0.18),
        tree_kws={"linewidths": 1.2},
    )
    grid.fig.suptitle("Gene-expression heatmap (row z-scores)", y=1.02)
    plt.show()


def plot_sample_dendrogram(linkage_matrix, sample_names):
    plt.figure(figsize=(10, 5))
    dendrogram(linkage_matrix, labels=sample_names, leaf_rotation=90)
    plt.title("Sample clustering (gene-wise z-scored expression)")
    plt.xlabel("Sample / spatial region")
    plt.ylabel("Ward linkage distance")
    plt.tight_layout()
    plt.show()


def make_gene_dendrogram_colour_func(linkage_matrix, leaf_cluster_labels):
    """
    Colour branches green/red when all descendant leaves belong to one cluster;
    colour mixed-cluster branches black.
    """
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
            colours[node_id] = GENE_CLUSTER_COLOURS.get(cluster_id, "grey")
        else:
            colours[node_id] = "black"

    return lambda node_id: colours.get(node_id, "black")


def plot_gene_dendrogram(linkage_matrix, genes, cluster_labels):
    colour_func = make_gene_dendrogram_colour_func(linkage_matrix, cluster_labels)

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
    n_genes = expression_scaled.shape[0]
    max_k = min(SILHOUETTE_MAX_CLUSTERS, n_genes - 1)

    if max_k < 2:
        print("Silhouette analysis skipped: too few genes.")
        return None

    candidate_ks = []
    scores = []

    for k in range(2, max_k + 1):
        labels = fcluster(gene_linkage, k, criterion="maxclust")
        if len(np.unique(labels)) < 2 or len(np.unique(labels)) >= n_genes:
            continue

        scores.append(silhouette_score(expression_scaled.to_numpy(), labels))
        candidate_ks.append(k)

    if not scores:
        print("Silhouette analysis skipped: no valid cluster solutions.")
        return None

    plt.figure(figsize=(7, 4))
    plt.plot(candidate_ks, scores, marker="o")
    plt.xlabel("Requested number of clusters")
    plt.ylabel("Silhouette score")
    plt.title("Silhouette analysis (row z-scores)")
    plt.tight_layout()
    plt.show()

    best_k = candidate_ks[int(np.argmax(scores))]
    print(f"Best silhouette result: requested k={best_k}, score={max(scores):.3f}")
    return best_k



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

    # Row z-scores are used for clustering and module profiles.
    expression_scaled = zscore_rows(expression)

    # Heatmap with hierarchical gene and sample clustering.
    plot_clustered_heatmap(expression_scaled)

    # Sample/spatial-region clustering.
    sample_linkage = linkage(
        expression_scaled.T.to_numpy(), method=LINKAGE_METHOD, metric="euclidean"
    )
    sample_labels = expression_scaled.columns.tolist()
    sample_cluster_labels = fcluster(
        sample_linkage, N_SAMPLE_CLUSTERS, criterion="maxclust"
    )
    sample_cluster_map = dict(zip(sample_labels, sample_cluster_labels))

    distance_matrix = squareform(
        pdist(expression_scaled.T.to_numpy(), metric="euclidean")
    )
    pd.DataFrame(
        distance_matrix, index=sample_labels, columns=sample_labels
    ).to_csv(OUTPUT_DIR / "distance_matrix.csv")

    plot_sample_dendrogram(sample_linkage, sample_labels)

    print("Sample clusters:", {
        cluster_id: [
            sample for sample, assigned in sample_cluster_map.items()
            if assigned == cluster_id
        ]
        for cluster_id in sorted(set(sample_cluster_labels))
    })

    # Marker genes are tested on the original (unscaled) expression values.
    marker_results = find_marker_genes(expression, sample_cluster_map)
    for cluster_id, marker_df in marker_results.items():
        output_path = OUTPUT_DIR / f"marker_genes_cluster_{cluster_id}.csv"
        marker_df.to_csv(output_path, index=False)
        print(f"Saved marker genes for sample cluster {cluster_id}: {output_path}")

    # Gene modules and module scores.
    print("Clustering genes into expression modules...")
    module_scores, gene_to_module, gene_linkage, gene_cluster_labels = (
        cluster_genes_and_score_modules(expression_scaled, N_GENE_CLUSTERS)
    )

    module_scores.to_csv(OUTPUT_DIR / "gene_module_scores.csv")
    module_enrichment = find_enriched_modules(module_scores, sample_cluster_map)
    module_enrichment.to_csv(
        OUTPUT_DIR / "gene_module_enrichment.csv", index=False
    )

    for sample_cluster in sorted(set(sample_cluster_labels)):
        subset = module_enrichment[
            module_enrichment["sample_cluster"] == sample_cluster
        ]
        print(f"\nModules associated with sample cluster {sample_cluster}:")
        if subset.empty:
            print("  No testable module results.")
        else:
            print(subset.head(10).to_string(index=False))

    # Save gene cluster assignments.
    gene_cluster_map = pd.Series(
        gene_cluster_labels, index=expression_scaled.index, name="gene_cluster"
    )
    gene_cluster_map.index.name = "gene"
    gene_cluster_map.reset_index().to_csv(
        OUTPUT_DIR / "gene_clusters.csv", index=False
    )

    plot_gene_dendrogram(
        gene_linkage, expression_scaled.index, gene_cluster_labels
    )
    
    plot_silhouette_analysis(expression_scaled, gene_linkage)

    sample_cluster_series = pd.Series(
        sample_cluster_labels,
        index=expression_scaled.columns,
        name="sample_cluster",
    )
    


if __name__ == "__main__":
    main()


