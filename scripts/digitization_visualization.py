
#!/usr/bin/env python3

import argparse
import json
import sys

import numpy as np
import vedo
import matplotlib.colors as mcolors


def main():
    parser = argparse.ArgumentParser(
        description="Display a single C3PO gene or gene-set pattern."
    )

    parser.add_argument(
        "json_file",
        help="Path to the JSON file containing centroids and samples."
    )

    parser.add_argument(
        "colormap",
        choices=[
            "white-green",
            "white-red",
            "white-blue",
            "white-black",
            "black-green",
            "black-red",
            "black-blue",
        ],
        help="Colour map to use."
    )

    parser.add_argument(
        "normalization",
        choices=["minmax", "none"],
        help=(
            "minmax = explicitly normalize values to 0-1; "
            "none = preserve the original square-root transformation "
            "and let Vedo automatically scale the colour map."
        )
    )

    args = parser.parse_args()

    # ---------------------------------------------------------
    # LOAD JSON
    # ---------------------------------------------------------
    with open(args.json_file, "r") as f:
        data = json.load(f)

    if "centroids" not in data or "samples" not in data:
        raise ValueError(
            "JSON must contain 'centroids' and 'samples'."
        )

    centroids = np.asarray(data["centroids"], dtype=float)
    values = np.asarray(data["samples"], dtype=float)

    if len(centroids) != len(values):
        raise ValueError(
            "The number of centroids must equal the number of samples."
        )

    if not np.all(np.isfinite(centroids)):
        raise ValueError("Centroid coordinates contain NaN or infinity.")

    if not np.all(np.isfinite(values)):
        raise ValueError("Expression values contain NaN or infinity.")

    if len(values) == 0:
        raise ValueError("The JSON file contains no samples.")

    # ---------------------------------------------------------
    # NORMALIZATION
    # ---------------------------------------------------------
    if args.normalization == "minmax":

        # Explicitly normalize original values to 0-1.
        vmin = values.min()
        vmax = values.max()

        if vmax == vmin:
            values = np.zeros_like(values)
        else:
            values = (values - vmin) / (vmax - vmin)

        colour_limits = (0, 1)

    else:

        # Preserve the square-root transformation from the original script.
        if np.any(values < 0):
            raise ValueError(
                "Square-root transformation requires non-negative values."
            )

        values = np.sqrt(values)

        # Vedo automatically scales the transformed values.
        colour_limits = None

    # ---------------------------------------------------------
    # COLOUR MAPS
    # ---------------------------------------------------------
    colour_maps = {
        "white-green": ["white", "green"],
        "white-red":   ["white", "red"],
        "white-blue":  ["white", "blue"],
        "white-black": ["white", "black"],
        "black-green": ["black", "green"],
        "black-red":   ["black", "red"],
        "black-blue":  ["black", "blue"],
    }

    cmap = mcolors.LinearSegmentedColormap.from_list(
        args.colormap,
        colour_maps[args.colormap]
    )

    # ---------------------------------------------------------
    # CREATE POINTS
    # ---------------------------------------------------------
    pts = vedo.Points(centroids).point_size(10)
    pts.pointdata["samples"] = values

    if colour_limits is not None:
        pts.cmap(
            cmap,
            on="points",
            vmin=0,
            vmax=1
        )
    else:
        pts.cmap(cmap, on="points")

    # ---------------------------------------------------------
    # BACKGROUND
    # ---------------------------------------------------------
    if args.colormap.startswith("black-"):
        background = "black"
    else:
        background = "white"

    # ---------------------------------------------------------
    # DISPLAY
    # ---------------------------------------------------------
    plt = vedo.Plotter(bg=background)
    plt.add(pts)

    print(f"JSON:          {args.json_file}")
    print(f"Colour map:    {args.colormap}")
    print(f"Normalization: {args.normalization}")
    print(f"Background:    {background}")

    plt.show()


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(1)

