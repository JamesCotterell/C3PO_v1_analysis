import pandas as pd
from statsmodels.stats.multitest import multipletests

# Select dataset:
# "Red cluster" or "Green cluster"
# BMP values are from Reactome; other pathways are from PANTHER.

datasets = {
    "Red cluster": {
        "pathway": ["BMP", "TGF_beta", "FGF", "Wnt", "HH"],
        "p_value": [0.229, 0.542, 0.111, 0.396, 0.0109],
    },
    "Green cluster": {
        "pathway": ["BMP", "TGF_beta", "FGF", "Wnt", "HH"],
        "p_value": [0.0182, 0.270, 0.592, 0.0278, 0.212],
    },
}

selected_cluster = "Red cluster"

df = pd.DataFrame(datasets[selected_cluster])

# Benjamini–Hochberg correction across the pathways in this dataset
df["FDR"] = multipletests(df["p_value"], method="fdr_bh")[1]

print(f"\nPathway enrichment: {selected_cluster}")
print(df.sort_values("FDR").to_string(index=False))

selected_cluster = "Green cluster"

df = pd.DataFrame(datasets[selected_cluster])

# Benjamini–Hochberg correction across the pathways in this dataset
df["FDR"] = multipletests(df["p_value"], method="fdr_bh")[1]

print(f"\nPathway enrichment: {selected_cluster}")
print(df.sort_values("FDR").to_string(index=False))
