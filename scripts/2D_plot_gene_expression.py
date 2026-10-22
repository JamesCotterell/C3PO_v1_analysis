#!/usr/bin/env python3

import argparse
from pathlib import Path

import numpy as np
from vedo import Plotter, Points, color_map


# ============================================================
# Command-line arguments
# ============================================================

parser = argparse.ArgumentParser(
    description="Visualise one gene's normalised expression across FACS gates."
)

parser.add_argument(
    "limb",
    type=int,
    choices=[127, 128],
    help="Limb dataset: 127 or 128.",
)

parser.add_argument(
    "--gene",
    required=True,
    help="Name of the single gene to display, e.g. Sox9.",
)

args = parser.parse_args()
limb = args.limb
gene = args.gene


# ============================================================
# Settings
# ============================================================

yscale = 2.5
my_vmax = 0.0005
point_radius = 8
point_opacity = 0.5

if limb == 127:
    data_dir = Path("./limb127_data")
    start_gate = 5
    end_gate = 20
    gate_prefix = "L127"
    screenshot_dir = data_dir / "screenshots"
else:
    data_dir = Path("./limb128_data")
    start_gate = 21
    end_gate = 35
    gate_prefix = "L128"
    screenshot_dir = data_dir / "screenshots_whitebackground"

gates = list(range(start_gate, end_gate))

print(f"\nLimb: {limb}")
print(f"Gene: {gene}")
print(f"Gates: {gates}")


# ============================================================
# Load normalised cell positions and centroids
# ============================================================

cells_file = data_dir / f"limb{limb}_normalized_data.txt"
centroids_file = data_dir / f"limb{limb}_normalized_data_centroids.txt"

cells = np.atleast_2d(np.genfromtxt(cells_file, delimiter="\t"))[:, :2]
centroids = np.atleast_2d(
    np.genfromtxt(centroids_file, delimiter="\t")
)[:, :2]

if not np.isfinite(cells).all():
    raise ValueError(f"Invalid coordinates in {cells_file}")
if not np.isfinite(centroids).all():
    raise ValueError(f"Invalid coordinates in {centroids_file}")

cells[:, 1] *= yscale
centroids[:, 1] *= yscale

print("Number of normalised cell positions:", len(cells))
print("Number of centroids:", len(centroids))


# ============================================================
# Load cell positions for each FACS gate
# ============================================================

gate_positions = {}

for gate in gates:
    gate_file = data_dir / f"{gate_prefix}_allsorts_P{gate}.fcs.csv"
    gate_data = np.atleast_2d(np.genfromtxt(gate_file, delimiter=","))

    if gate_data.shape[1] < 2:
        raise ValueError(f"Expected at least two columns in {gate_file}")

    gate_data = gate_data[:, :2]

    if not np.isfinite(gate_data).all():
        raise ValueError(f"Non-finite coordinates found in {gate_file}")

    # Scale y coordinates once, as in the original script.
    gate_data[:, 1] *= yscale
    gate_positions[gate] = gate_data

    print(f"Loaded P{gate}: {len(gate_data)} cells")


# ============================================================
# Load counts and calculate total expression for each gate
# ============================================================

gene_names = None
counts_by_gate = []
total_expression_by_gate = []

for gate in gates:
    count_file = data_dir / (
        f"counts_{gate_prefix}P{gate}_doubletrimmed_sorted_markdup.txt"
    )

    current_gene_names = []
    current_counts = []

    with count_file.open("r") as f:
        for line in f:
            # Skip the featureCounts command comment and column-header row.
            if line.startswith("#") or line.startswith("Geneid\t"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 7:
                continue

            current_gene_names.append(fields[0])
            current_counts.append(int(fields[-1]))

    current_counts = np.asarray(current_counts, dtype=np.int64)

    if gene_names is None:
        gene_names = current_gene_names
    elif current_gene_names != gene_names:
        raise ValueError(f"Gene names or order differ in {count_file}")

    if len(current_counts) != len(gene_names):
        raise ValueError(f"Gene count mismatch in {count_file}")

    # Preserve the original denominator calculation: exclude the first
    # two gene records from the total, but retain them in the gene list.
    total_expression = int(np.sum(current_counts[2:]))

    if total_expression <= 0:
        raise ValueError(f"Total expression is zero for P{gate}")

    counts_by_gate.append(current_counts)
    total_expression_by_gate.append(total_expression)

    print(f"P{gate}: total expression = {total_expression}")

counts_by_gate = np.asarray(counts_by_gate, dtype=np.float64)
total_expression_by_gate = np.asarray(
    total_expression_by_gate, dtype=np.float64
)

gene_to_index = {name: index for index, name in enumerate(gene_names)}

if gene not in gene_to_index:
    raise ValueError(f"Gene '{gene}' was not found in the count files.")

gene_index = gene_to_index[gene]

# Normalised expression for this gene in each gate:
# gene counts / total counts for that gate.
expression = (
    counts_by_gate[:, gene_index] / total_expression_by_gate
)

print(f"\nNormalised expression for {gene}:")
for gate, value in zip(gates, expression):
    print(f"  P{gate}: {value:.8g}")


# ============================================================
# Plot the gene across the FACS gates
# ============================================================

gate_actors = []

for gate_index, gate in enumerate(gates):
    colour = color_map(
        expression[gate_index],
        "gray",
        vmin=0,
        vmax=my_vmax,
    )

    gate_points = Points(
        gate_positions[gate],
        r=point_radius,
        c=colour,
    ).opacity(point_opacity)

    gate_actors.append(gate_points)

plotter = Plotter(axes=0, bg="black")
plotter.show(gate_actors)


# ============================================================
# Save screenshot
# ============================================================

screenshot_dir.mkdir(parents=True, exist_ok=True)
screenshot_file = screenshot_dir / f"{gene}.png"

plotter.screenshot(str(screenshot_file))

print("\nScreenshot saved to:")
print(screenshot_file)

plotter.interactive().close()


