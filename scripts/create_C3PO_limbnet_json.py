
#!/usr/bin/env python3

import argparse
import json
import sys
from pathlib import Path


def get_gene_expression(count_file, gene):
    """Return the expression value for a gene, or zero if not found."""
    with count_file.open() as f:
        # Skip the first two header lines, as in the original script.
        next(f, None)
        next(f, None)

        for line in f:
            parts = line.split()
            if parts and parts[0] == gene:
                return float(parts[-1])

    return 0.0


def convert_json(gates_file, data_dir, count_prefix, gene):
    """Add normalized gene expression values to a gate-mapped JSON."""
    with gates_file.open() as f:
        data = json.load(f)

    samples = data["samples"]

    # Extract gate numbers from labels such as P5, P19, P34.
    gates = {
        int(sample[1:])
        for sample in samples
        if isinstance(sample, str)
        and sample.startswith("P")
        and sample[1:].isdigit()
    }

    # Read expression once for each unique gate.
    expr_by_gate = {}

    for gate in gates:
        count_file = data_dir / (
            f"counts_{count_prefix}P{gate}_doubletrimmed_sorted_markdup.txt"
        )

        if not count_file.exists():
            raise FileNotFoundError(
                f"Missing expression file: {count_file}"
            )

        expr_by_gate[gate] = get_gene_expression(count_file, gene)

    # Match expression values to the original sample ordering.
    raw_values = []

    for sample in samples:
        if (
            isinstance(sample, str)
            and sample.startswith("P")
            and sample[1:].isdigit()
        ):
            raw_values.append(expr_by_gate.get(int(sample[1:]), 0.0))
        else:
            raw_values.append(0.0)

    # Normalize by the maximum expression for this gene in this limb.
    max_val = max(raw_values, default=0.0)

    if max_val > 0:
        norm_values = [min(value / max_val, 1.0) for value in raw_values]
    else:
        norm_values = raw_values

    data["samples"] = norm_values
    data["name"] = gene

    return data


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Generate normalized gene-expression JSON files for Limb4, "
            "Limb5, and their average."
        ),
        epilog=(
            "Example: python generate_limb_expression.py Sox9\n"
            "Input count files and gate-mapped JSON files are expected "
            "in the current directory's data folders."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "gene",
        help="Gene name as it appears in the first column of the count files."
    )

    parser.add_argument(
        "--limb127-dir",
        type=Path,
        default=Path("./limb127_data"),
        help="Directory containing L127 count files (default: ./limb127_data)."
    )

    parser.add_argument(
        "--limb128-dir",
        type=Path,
        default=Path("./limb128_data"),
        help="Directory containing L128 count files (default: ./limb128_data)."
    )

    parser.add_argument(
        "--gates-dir",
        type=Path,
        default=Path("."),
        help="Directory containing limb4_colour_mapped.json and "
             "limb5_colour_mapped.json (default: current directory)."
    )

    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("./C3PO_limbnet_digitizations"),
        help="Directory for output JSON files."
    )

    args = parser.parse_args()

    gene = args.gene
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Limb4 uses L128 data; Limb5 uses L127 data.
    limb_configs = {
        "Limb4": {
            "gates_file": args.gates_dir / "limb4_colour_mapped.json",
            "data_dir": args.limb128_dir,
            "count_prefix": "L128",
            "output": args.output_dir / f"limb4_expression_{gene}.json",
        },
        "Limb5": {
            "gates_file": args.gates_dir / "limb5_colour_mapped.json",
            "data_dir": args.limb127_dir,
            "count_prefix": "L127",
            "output": args.output_dir / f"limb5_expression_{gene}.json",
        },
    }

    limb_results = {}

    # Generate limb-specific JSON files.
    for limb, cfg in limb_configs.items():
        if not cfg["gates_file"].exists():
            raise FileNotFoundError(
                f"Missing gate-mapped JSON: {cfg['gates_file']}"
            )

        new_json = convert_json(
            cfg["gates_file"],
            cfg["data_dir"],
            cfg["count_prefix"],
            gene,
        )

        if all(value == 0 for value in new_json["samples"]):
            print(f"Skipping {limb} for {gene}: no expression detected.")
            continue

        with cfg["output"].open("w") as f:
            json.dump(new_json, f, indent=2)

        print(f"Saved: {cfg['output']}")
        limb_results[limb] = new_json

    # Generate average across limbs that contain expression.
    if not limb_results:
        print(f"No expression detected in either limb for {gene}.")
        return

    sample_lists = [
        result["samples"] for result in limb_results.values()
    ]

    avg_samples = [
        sum(values) / len(values)
        for values in zip(*sample_lists)
    ]

    # Use the first available limb's JSON as the structural template.
    avg_json = next(iter(limb_results.values())).copy()
    avg_json["samples"] = avg_samples
    avg_json["name"] = f"{gene}_avg_L4_L5"

    avg_path = args.output_dir / f"limb_avg_expression_{gene}.json"

    with avg_path.open("w") as f:
        json.dump(avg_json, f, indent=2)

    print(f"Saved: {avg_path}")


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, KeyError, json.JSONDecodeError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(1)

