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

## Dataset naming conventions

Two principal datasets are referred to by different names in the repository and manuscript. These names reflect two different numbering systems:

* **Limb 127 = Limb 5**
* **Limb 128 = Limb 4**

The numbers **127 and 128** identify the original experiments, whereas **Limb 5 and Limb 4** refer to the order in which the respective datasets were presented in the manuscript.

Consequently, the same dataset may be referred to by either name, depending on the context. For example, files associated with `limb127_data/` correspond to Limb 5 in the manuscript, while files associated with `limb128_data/` correspond to Limb 4.

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

## Data and reproducibility

The data provided here support the analysis and visualisation of the C3PO experimental results.

The repository includes processed datasets, normalised FACS measurements, digitised tissue data and analysis scripts. These resources are intended to facilitate inspection of the data, reproduction of the analyses and further exploration of the results.

Because the repository contains selected data generated during the study, it should not be interpreted as an exhaustive archive of every experiment or intermediate result produced during the project.

## Citation

If you use the code or data in this reposito

