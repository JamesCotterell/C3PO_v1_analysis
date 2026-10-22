#!/usr/bin/env python3

import argparse
from vedo import *


# ============================================================
# 1. Parse command-line arguments
# ============================================================

parser = argparse.ArgumentParser(
    description=(
        "Load and visualise Ectoderm (E) and Mesenchyme (M) "
        "3D point clouds using Vedo."
    )
)

parser.add_argument(
    "ectoderm_file",
    help="Path to the Ectoderm (E) input text file."
)

parser.add_argument(
    "mesenchyme_file",
    help="Path to the Mesenchyme (M) input text file."
)

args = parser.parse_args()


# ============================================================
# 2. Load the two point clouds
# ============================================================

try:
    points_E = Points(args.ectoderm_file)
    points_M = Points(args.mesenchyme_file)
except Exception as e:
    parser.error(f"Could not load input files: {e}")


# ============================================================
# 3. Set colours and point sizes
# ============================================================

points_E.color("blue5").point_size(7)
points_M.color("red5").point_size(7)


# ============================================================
# 4. Remove outliers using the original cutting plane
# ============================================================

points_E.cut_with_plane(origin=(0, 0.4, 0), normal="-y")
points_M.cut_with_plane(origin=(0, 0.4, 0), normal="-y")


# ============================================================
# 5. Preserve the original scaling
# ============================================================

points_E.scale([1, 2, 0.2])
points_M.scale([1, 2, 0.15])


# ============================================================
# 6. Remove points that are very close together
# ============================================================

points_E.subsample(0.01).print()
points_M.subsample(0.01).print()


# ============================================================
# 7. Display both point clouds in the same 3D plot
# ============================================================

show(
    points_E,
    points_M,
    axes=1,
    viewup="z",
    bg="white",
    title="Limb 76: E (blue) and M (red)",
).close()


