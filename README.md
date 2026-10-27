# C3PO v1 — Analysis Code and Data

This repository contains analysis code and associated data for **Cell 3D Positioning by Optical encoding (C3PO)**, as described in the manuscript:

*Cell 3D Positioning by Optical encoding (C3PO) and its application to spatial transcriptomics.*

C3PO is a method for preserving three-dimensional positional information in biological tissues before dissociation, enabling spatially resolved molecular profiling of cells.

This repository provides selected data and analysis scripts associated with the study, supporting transparency, reproducibility and reuse of the results.

## Repository contents

The repository is organised as follows:

| File or directory                       | Description                                                                                      |
| --------------------------------------- | ------------------------------------------------------------------------------------------------ |
| `scripts/`                              | Scripts used for data analysis, processing and visualisation.                                    |
| `Digitizations/`                        | Digitised tissue geometries and associated data.                                                 |
| `limb127_data/`                         | Data associated with experimental sample Limb 127, also referred to as Limb 5 in the manuscript. |
| `limb128_data/`                         | Data associated with experimental sample Limb 128, also referred to as Limb 4 in the manuscript. |
| `Normalized_FACS_data/`                 | Normalised fluorescence-activated cell sorting (FACS) data.                                      |
| `processed_data_combined_L127L128.xlsx` | Combined processed data for experimental samples Limb 127 and Limb 128.                          |

## Dataset naming conventions and experimental groups

The repository contains eight principal datasets, referred to by different names in the repository and manuscript. The datasets fall into three experimental groups, each supporting a different aspect of the C3PO method.

### 1. Transcriptomics datasets

These datasets combine C3PO spatial encoding with gene-expression measurements and are used for the spatial transcriptomics analyses.

| Experimental dataset | Manuscript name | Description                                         |
| -------------------- | --------------- | --------------------------------------------------- |
| Limb 128             | Limb 4          | C3PO spatial encoding with transcriptomic profiling |
| Limb 127             | Limb 5          | C3PO spatial encoding with transcriptomic profiling |

These datasets are used by the following scripts:

* `clustering_script_with_rowwise_scaling.py` — hierarchical clustering of spatial gene-expression profiles.
* `2D_plot_gene_expression.py` — visualisation of gene expression across FACS gates.
* `create_C3PO_limbnet_json.py` — mapping gene-expression measurements onto digitised 3D limb geometries.
* `digitization_visualization.py` — visualisation of values mapped onto digitised 3D tissue geometries.

### 2. Ectoderm and mesenchyme separation datasets

These datasets were generated using C3PO without transcriptomic profiling. The ectoderm and mesenchyme were separated before fluorescence-activated cell sorting (FACS), allowing the fluorescence distributions of the two tissue compartments to be analysed independently.

| Experimental dataset | Manuscript name | Description                                                |
| -------------------- | --------------- | ---------------------------------------------------------- |
| Limb 76              | Limb 1          | C3PO encoding; ectoderm and mesenchyme analysed separately |
| Limb 89              | Limb 2          | C3PO encoding; ectoderm and mesenchyme analysed separately |
| Limb 87              | Limb 3          | C3PO encoding; ectoderm and mesenchyme analysed separately |

These datasets are used by:

* `3D_plot_Mesenchyme_and_Ectoderm.py` — visualisation and comparison of the 3D fluorescence distributions of ectoderm and mesenchyme.

### 3. Optical encoding validation datasets

These datasets were generated using C3PO without transcriptomic profiling. They were used to validate the ability to resolve the eight optical addresses encoded into the tissue by photobleaching.

| Experimental dataset | Manuscript reference   | Description                               |
| -------------------- | ---------------------- | ----------------------------------------- |
| Limb 74              | Figure S5, top limb    | Validation of the eight optical addresses |
| Limb 80              | Figure S5, middle limb | Validation of the eight optical addresses |
| Limb 69              | Figure S5, bottom limb | Validation of the eight optical addresses |

These datasets are used by:

* `3D_plot_allornothinglimbs_with_clusters.py` — hierarchical clustering and visualisation of fluorescence measurements to assess whether the eight intended optical encoding zones can be distinguished.

### Naming conventions

The experimental dataset numbers identify the original experiments, whereas the manuscript limb numbers indicate the order in which the transcriptomics and ectoderm–mesenchyme datasets are presented in the manuscript. The three optical encoding validation datasets are identified by their positions in Figure S5.

Consequently, the same dataset may be referred to by different names depending on the context. For example, files associated with `limb127_data/` correspond to Limb 5 in the manuscript, while files associated with `limb128_data/` correspond to Limb 4.

Transcriptomic analyses are restricted to Limbs 127 and 128. The remaining six datasets support the evaluation of C3PO optical encoding and the analysis of tissue-compartment-specific fluorescence distributions.

Not all data generated during the study were included in the manuscript. This repository contains selected data and analysis resources associated with the study; therefore, not every file necessarily corresponds to a result or figure presented in the paper.

## Using the code

The `scripts/` directory contains the analysis and visualisation scripts used in the study.

Individual scripts may require specific input files, Python packages or software environments. Please refer to the individual scripts for their input requirements and usage instructions.

Where supported, command-line help can be displayed using the `-h` option. For example:

```bash
python3 script_name.py -h
```

Replace `script_name.py` with the name of the script you wish to run.

Some analyses may require input data from the directories included in this repository. Ensure that the relevant data files are available and that any required software dependencies are installed before running a script.

### `3D_plot_allornothinglimbs_with_clusters.py`

**Purpose:** Identify and visualise eight clusters in three-dimensional fluorescence data, corresponding to the eight intended optical encoding zones generated by C3PO.

This script is associated with the optical encoding validation experiments using Limbs 74, 80 and 69, shown in Figure S5. These datasets do not include transcriptomic profiling.

In this experiment, all-or-nothing photobleaching along three orthogonal axes was used to generate eight spatial zones, each representing a distinct combination of the three binary encoding dimensions. The script analyses three fluorescence measurements corresponding to these dimensions to assess the distribution of the encoded points in three-dimensional fluorescence space.

**Analysis workflow:**

1. **Load the data:** Read an input text file containing three numerical fluorescence coordinates per point.
2. **Normalise the coordinates:** Divide each fluorescence coordinate by the maximum value observed in its respective dimension.
3. **Perform hierarchical clustering:** Apply agglomerative clustering with Ward linkage and Euclidean distance, specifying eight clusters.
4. **Calculate cluster properties:** Determine the centroid and number of points assigned to each cluster.
5. **Visualise the results:** Generate a three-dimensional scatter plot showing individual points coloured by cluster assignment and the eight cluster centroids.

The clustering is intended to evaluate whether the fluorescence measurements resolve into eight groups corresponding to the eight intended C3PO optical encoding zones.

**Usage:**

```bash
python3 3D_plot_allornothinglimbs_with_clusters.py <input_file>
```

The input file must contain exactly three columns of finite numerical values, with one point per row. At least eight points are required.

**Output:** The script prints the number of input points, the normalisation maxima, and the size and centroid coordinates of each cluster. It also displays a three-dimensional plot of the clustered points and their centroids.

### `3D_plot_Mesenchyme_and_Ectoderm.py`

**Purpose:** Visualise and compare the three-dimensional fluorescence point distributions of ectoderm and mesenchyme isolated from a limb bud and analysed separately by fluorescence-activated cell sorting (FACS).

This script is associated with the ectoderm–mesenchyme separation experiments using Limbs 76, 89 and 87, corresponding to manuscript Limbs 1, 2 and 3. These datasets do not include transcriptomic profiling.

The script displays the two tissue compartments together in a common three-dimensional coordinate system, allowing their fluorescence point distributions to be inspected and compared visually.

**Analysis workflow:**

1. **Load the data:** Read the ectoderm and mesenchyme point clouds from two separate input files.
2. **Distinguish the tissue compartments:** Display ectoderm points in blue and mesenchyme points in red.
3. **Remove outliers:** Apply the same cutting plane to both point clouds to exclude points on one side of the specified plane.
4. **Adjust visual scaling:** Apply the original axis-specific scaling factors to each point cloud.
5. **Reduce point density:** Subsample each point cloud to reduce the number of closely spaced points displayed.
6. **Visualise the results:** Display both point clouds in the same interactive three-dimensional plot, with axes and a white background.

**Usage:**

```bash
python3 3D_plot_Mesenchyme_and_Ectoderm.py <ectoderm_file> <mesenchyme_file>
```

Provide the paths to the ectoderm and mesenchyme input files as the first and second arguments, respectively.

**Output:** An interactive three-dimensional visualisation showing ectoderm in blue and mesenchyme in red.

**Note:** The script applies different scaling factors to the two point clouds and performs geometric filtering and subsampling for visualisation. Consequently, the displayed coordinates are transformed for plotting and should not be interpreted as an unmodified quantitative comparison of the original fluorescence measurements.

### `2D_plot_gene_expression.py`

**Purpose:** Visualise the spatial distribution of a selected gene's normalised expression across the fluorescence-activated cell sorting (FACS) gates for Limb 127 (Limb 5) or Limb 128 (Limb 4).

This script uses the transcriptomics datasets from Limbs 127 and 128, which combine C3PO spatial encoding with gene-expression measurements.

The script combines normalised two-dimensional cell-position data with gene expression counts obtained from the individual FACS gates. It generates a spatial visualisation in which the cells associated with each gate are coloured according to the normalised expression level of the selected gene in that gate.

**Analysis workflow:**

1. **Select the dataset and gene:** Specify Limb 127 or Limb 128 and the gene to visualise.
2. **Load cell positions:** Read the normalised two-dimensional cell coordinates, centroids and gate-specific cell-position data.
3. **Load gene counts:** Read the featureCounts output files corresponding to the individual FACS gates.
4. **Calculate normalised expression:** For each gate, divide the selected gene's read count by the total read counts across the genes included in the denominator. The first two gene records are excluded from the denominator, following the original analysis.
5. **Map expression to cell positions:** Assign each gate's normalised expression value to its corresponding two-dimensional cell-position point cloud.
6. **Visualise and export:** Display the gate-specific point clouds using a grey colour scale, with darker and lighter shades representing different normalised expression levels, and save a screenshot as a PNG file.

**Dataset and gate selection:**

* **Limb 127 (Limb 5):** Uses FACS gates P5–P19.
* **Limb 128 (Limb 4):** Uses FACS gates P21–P34.

These gate ranges correspond to the datasets configured in the script.

**Usage:**

To plot a gene in Limb 127:

```bash
python3 2D_plot_gene_expression.py 127 --gene Sox9
```

To plot the same gene in Limb 128:

```bash
python3 2D_plot_gene_expression.py 128 --gene Sox9
```

Replace `Sox9` with the name of the gene of interest.

**Input files:** The script requires the normalised cell-position data, gate-specific cell-position files and corresponding featureCounts gene-count files for the selected limb dataset. These files are expected in the appropriate `limb127_data/` or `limb128_data/` directory.

**Output:** An interactive two-dimensional visualisation of the gate-specific cell-position point clouds, together with a PNG screenshot named after the selected gene. Screenshots are saved in the dataset's configured screenshot directory.

**Interpretation:** The visualisation shows how the selected gene's normalised read-count abundance varies across the FACS gates and maps those differences onto their associated cell-position distributions. It provides a way to examine gene expression patterns in the reconstructed spatial context of the limb bud.

### `clustering_script_with_rowwise_scaling.py`

**Purpose:** Perform hierarchical clustering of the combined spatial gene-expression data from Limb 127 (Limb 5) and Limb 128 (Limb 4), identifying groups of genes with similar expression patterns and visualising relationships between spatial samples or regions.

This analysis uses the two transcriptomics datasets. The other six datasets do not contain transcriptomic measurements and are not included in this gene-expression clustering analysis.

The script analyses the combined gene-expression matrix, in which genes are represented by rows and spatial samples or regions by columns. It applies row-wise z-score normalisation so that gene clustering emphasises the relative expression pattern of each gene across spatial regions, rather than differences in absolute expression levels between genes.

**Analysis workflow:**

1. **Load the combined expression matrix:** Read `../processed_data_combined_L127L128.xlsx`, containing gene-expression measurements for spatial samples or regions from both limbs.
2. **Prepare the data:** Convert expression values to numeric format, remove entirely empty rows and columns, impute remaining missing values using the corresponding gene's mean, and remove genes with zero variance.
3. **Normalise expression profiles:** Calculate row-wise z-scores for each gene across all spatial samples or regions.
4. **Cluster spatial samples or regions:** Apply hierarchical clustering using Ward linkage and Euclidean distance to generate the sample dendrogram.
5. **Cluster genes:** Apply hierarchical clustering to the normalised gene-expression profiles, requesting two gene clusters.
6. **Export gene lists:** Save the genes assigned to each cluster in separate plain-text files for downstream pathway-enrichment analysis using PANTHER.
7. **Visualise the results:** Generate a hierarchically clustered gene-expression heatmap, a sample-clustering dendrogram, a gene-clustering dendrogram and a silhouette analysis plot to evaluate alternative numbers of gene clusters.

**Input file:**

`../processed_data_combined_L127L128.xlsx`

The input is expected to contain genes in rows and spatial samples or regions in columns, with gene names in the first column.

**Pre-processing note:** Low-complexity voxels have already been excluded from this input dataset before the downstream clustering analysis. The script therefore operates on the pre-filtered dataset and does not perform the initial voxel-filtering step itself.

**Clustering settings:**

* **Gene clusters:** 2 requested
* **Linkage method:** Ward
* **Distance metric:** Euclidean
* **Expression scaling:** Row-wise z-score normalisation across spatial samples or regions

**Output files:**

| File                  | Description                                               |
| --------------------- | --------------------------------------------------------- |
| `genes_cluster_1.txt` | Gene names assigned to gene cluster 1, one gene per line. |
| `genes_cluster_2.txt` | Gene names assigned to gene cluster 2, one gene per line. |

The script writes one file for each cluster returned by the clustering algorithm. The output files contain gene names suitable for use as input gene lists for PANTHER pathway-enrichment analysis.

**Visualisations:**

* **Clustered heatmap:** Displays normalised gene-expression patterns with hierarchical clustering of both genes and spatial samples or regions.
* **Sample dendrogram:** Shows hierarchical relationships between spatial samples or regions based on their normalised gene-expression profiles.
* **Gene dendrogram:** Shows hierarchical relationships between genes, with branches coloured to distinguish the two gene clusters where possible.
* **Silhouette analysis:** Evaluates clustering solutions across alternative numbers of gene clusters and displays their silhouette scores. This analysis does not change the primary gene-clustering solution.

**Usage:**

Run the script from the `scripts` directory, assuming the Excel file is located in the repository's parent directory:

```bash
python3 clustering_script_with_rowwise_scaling.py
```

The script displays the plots interactively and saves the gene-list text files in the current working directory. The output directory is configured by `OUTPUT_DIR` in the script.

**Interpretation:** Gene clustering identifies groups of genes with similar relative expression profiles across the combined spatial dataset. The resulting gene lists can be used for downstream pathway-enrichment analysis. The sample dendrogram and clustered heatmap provide complementary visualisations of similarities between spatial regions.

### `compute_FDR.py`

**Purpose:** Apply the Benjamini–Hochberg false discovery rate (FDR) correction to pathway-enrichment p-values for the red and green gene clusters identified in the C3PO analysis.

The script uses pathway-enrichment p-values obtained from external analyses to calculate FDR-adjusted p-values across five signalling pathways: BMP, TGF-β, FGF, Wnt and Hedgehog (HH). The p-values are hard-coded in the script, with BMP pathway values obtained from Reactome and the other pathway values obtained from PANTHER.

**Analysis workflow:**

1. **Load pathway-enrichment results:** Use the predefined p-values for the five signalling pathways in each of the red and green gene clusters.
2. **Correct for multiple testing:** Apply the Benjamini–Hochberg procedure separately to the five pathway p-values for each cluster.
3. **Display the results:** Print the pathway names, original p-values and FDR-adjusted p-values, sorted by adjusted p-value.

**Pathways analysed:**

* BMP
* TGF-β
* FGF
* Wnt
* Hedgehog (HH)

**Usage:**

Run the script using Python:

```bash
python3 compute_FDR.py
```

No input files or command-line arguments are required, as the pathway names and p-values are defined within the script.

**Output:** The script prints separate tables for the red and green clusters, listing each pathway's original p-value and Benjamini–Hochberg-adjusted FDR value.

**Interpretation:** The adjusted p-values account for multiple pathway tests within each cluster, helping assess the statistical significance of the pathway-enrichment results after correction for multiple comparisons.

### `create_C3PO_limbnet_json.py`

**Purpose:** Generate gene-specific JSON files containing normalised gene-expression values mapped onto the digitised 3D limb geometries used by LimbNet.

This script uses the transcriptomics datasets from Limb 127 (Limb 5) and Limb 128 (Limb 4). It is not intended for the six datasets without transcriptomic profiling.

The script combines gate-mapped JSON files with gene-expression counts from Limb 127 and Limb 128. For a specified gene, it extracts the expression count for each FACS gate, normalises the values within each limb by that limb's maximum expression value, and adds the normalised values to the corresponding JSON structure. It also generates an average expression JSON by averaging the normalised values across the limbs in which expression was detected.

**Inputs:**

* **Gene name:** The gene to process, specified as a positional command-line argument.
* **Limb 127 data directory:** Contains the L127 gene-count files.
* **Limb 128 data directory:** Contains the L128 gene-count files.
* **Digitisations directory:** Contains `limb4_colour_mapped.json` and `limb5_colour_mapped.json`, which provide the gate-mapped limb geometries.
* **Output directory:** The destination for the generated JSON files.

The script expects gene-count files to follow the naming convention:

```text
counts_L127P<gate>_doubletrimmed_sorted_markdup.txt
counts_L128P<gate>_doubletrimmed_sorted_markdup.txt
```

For each gate, the script reads the specified gene's count from the first column of the gene-count data, skipping the first two lines. Expression values are normalised independently for each limb by dividing each gate's value by the maximum value for that gene in that limb. If the gene has no detected expression in a limb, the script skips that limb's output file. If expression is detected in only one limb, the average output is based on that limb alone.

**Usage:**

Run the script from the directory containing `create_C3PO_limbnet_json.py`, specifying the gene and all four directories:

```bash
python3 create_C3PO_limbnet_json.py Sox9 \
    --limb127-dir ./limb127_data \
    --limb128-dir ./limb128_data \
    --gates-dir ./Digitizations \
    --output-dir ./C3PO_limbnet_digitizations
```

Replace `Sox9` with the gene of interest and adjust the directory paths to match your file structure. The directories can be supplied as relative or absolute paths.

For a list of available arguments, run:

```bash
python3 create_C3PO_limbnet_json.py -h
```

**Outputs:**

The script writes the following files to the specified output directory, provided expression is detected in at least one limb:

* `limb4_expression_<gene>.json` — normalised expression mapped onto Limb 4 (Limb 128).
* `limb5_expression_<gene>.json` — normalised expression mapped onto Limb 5 (Limb 127).
* `limb_avg_expression_<gene>.json` — average normalised expression across the limbs in which expression was detected.

The output directory is created automatically if it does not already exist.

**Important notes:**

* Normalisation is performed separately for each limb; values therefore represent relative expression within each limb, rather than directly comparable absolute expression levels.
* The average JSON is calculated by averaging normalised values at corresponding positions in the two input sample arrays. The script uses the first available limb's JSON as the structural template for this output.
* Samples that do not correspond to a valid FACS gate label receive a value of zero.
* The script maps gate-level expression measurements onto the digitised geometry; the values are not single-cell expression measurements.

### `digitization_visualization.py`

**Purpose:** Visualise gene-expression patterns or other numerical sample values mapped onto digitised 3D tissue geometries stored in JSON files. The script uses Vedo to display the digitised centroid positions as coloured 3D points, with configurable colour maps and normalisation options.

This is a general-purpose visualisation script. It can be used to inspect JSON files containing digitised geometries and associated sample values, including gene-expression JSON files generated by `create_C3PO_limbnet_json.py`.

**Inputs:**

* **JSON file:** A JSON file containing `centroids` (3D coordinates) and `samples` (numerical values associated with each centroid).
* **Colour map:** Specifies the colour gradient used to represent the sample values.
* **Normalisation:** Determines how the sample values are transformed before visualisation.

**Available colour maps:**

| Colour map    | Gradient       | Background |
| ------------- | -------------- | ---------- |
| `white-green` | White to green | White      |
| `white-red`   | White to red   | White      |
| `white-blue`  | White to blue  | White      |
| `white-black` | White to black | White      |
| `black-green` | Black to green | Black      |
| `black-red`   | Black to red   | Black      |
| `black-blue`  | Black to blue  | Black      |

The background colour is selected automatically based on the chosen colour map: maps beginning with `white-` use a white background, while maps beginning with `black-` use a black background.

**Normalisation options:**

* **`minmax`:** Applies min–max normalisation to scale the input values to the range 0–1. The colour map is explicitly scaled to this range. If all input values are identical, they are set to zero.
* **`none`:** Preserves the original values up to a square-root transformation, which is applied before visualisation. Vedo automatically scales the colour map to the transformed data range. This option requires non-negative input values.

**Usage:**

Run the script with the JSON file, colour map and normalisation method as positional arguments:

```bash
python3 digitization_visualization.py <json_file> <colormap> <normalization>
```

For example, to display the Sox9 expression pattern using a white-to-green colour map and min–max normalisation:

```bash
python3 digitization_visualization.py \
    C3PO_limbnet_digitizations/limb_avg_expression_Sox9.json \
    white-green \
    minmax
```

To use the square-root transformation instead:

```bash
python3 digitization_visualization.py \
    C3PO_limbnet_digitizations/limb_avg_expression_Sox9.json \
    white-green \
    none
```

**Output:**

The script opens an interactive 3D visualisation window displaying the centroid positions coloured according to their associated sample values. It also prints the input JSON file, selected colour map, normalisation method and background colour to the terminal. The script does not save an image or modify the input JSON file.

**Validation:**

The script checks that the JSON contains both `centroids` and `samples`, that their lengths match, that the coordinates and values are finite, and that the file contains at least one sample. Negative values are rejected when using the `none` normalisation option because the square-root transformation requires non-negative values.

## Data and reproducibility

The data provided here support the analysis and visualisation of the C3PO experimental results.

The repository includes processed datasets, normalised FACS measurements, digitised tissue data and analysis scripts. These resources are intended to facilitate inspection of the data, reproduction of the analyses and further exploration of the results.

The eight principal datasets support three complementary aspects of the study: spatial transcriptomics, comparison of ectoderm and mesenchyme fluorescence distributions, and validation of the eight intended optical encoding zones.

Because the repository contains selected data generated during the study, it should not be interpreted as an exhaustive archive of every experiment or intermediate result produced during the project.

## Data processing and filtering

The combined gene-expression dataset, `processed_data_combined_L127L128.xlsx`, contains data from Limb 127 (Limb 5) and Limb 128 (Limb 4). Low-complexity voxels have already been excluded from this dataset prior to the downstream clustering analysis.

Consequently, the clustering script operates on the pre-filtered dataset and does not perform this voxel-filtering step itself.

## Citation

If you use the code or data in this repository, please cite the associated C3PO manuscript:

*Cell 3D Positioning by Optical encoding (C3PO) and its application to spatial transcriptomics.*

