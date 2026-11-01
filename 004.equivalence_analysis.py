#!/usr/bin/env python3

from pathlib import Path
import csv
import math
import re
import sys

import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parent

BASELINE = ROOT / "baseline" / "system"
MODULARIZED = ROOT / "modularized" / "system"
RESULTS = ROOT / "comparison_results"

STATISTICS = RESULTS / "equivalence_statistics.csv"

def input_output_prefix():
    input_path = ROOT / "dn_generic.in"
    if not input_path.is_file():
        raise SystemExit(f"ERROR: Missing canonical input: {input_path}")
    values = []
    for line in input_path.read_text().splitlines():
        fields = line.split()
        if len(fields) >= 2 and fields[0] == "dn_output_prefix":
            values.append(fields[1])
    if len(values) != 1:
        raise SystemExit("ERROR: Expected one dn_output_prefix in dn_generic.in")
    return values[0]


OUTPUT_PREFIX = input_output_prefix()
FINAL_MOL2 = f"{OUTPUT_PREFIX}.completed.denovo_build.mol2"


def split_mol2(path):
    text = path.read_text(errors="replace")
    marker = "@<TRIPOS>MOLECULE"

    return [
        marker + chunk
        for chunk in text.split(marker)[1:]
        if chunk.strip()
    ]


def molecule_count(path):
    if not path.exists():
        return 0

    return len(split_mol2(path))


def molecule_name(block):
    lines = block.splitlines()

    for i, line in enumerate(lines):
        if line.strip() == "@<TRIPOS>MOLECULE":
            if i + 1 < len(lines):
                return lines[i + 1].strip()

    return ""


def parse_properties(block):
    properties = {}

    for line in block.splitlines():
        stripped = line.strip()

        if not stripped.startswith("##########"):
            continue

        payload = stripped.lstrip("#").strip()

        if ":" not in payload:
            continue

        key, value = payload.split(":", 1)
        properties.setdefault(key.strip(), []).append(value.strip())

    return properties


def numeric_value(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def parse_atoms(block):
    atoms = []

    if "@<TRIPOS>ATOM" not in block:
        return atoms

    section = block.split("@<TRIPOS>ATOM", 1)[1]

    match = re.search(r"\n@<TRIPOS>", section)

    if match:
        section = section[:match.start()]

    for line in section.splitlines():
        fields = line.split()

        if len(fields) < 6:
            continue

        try:
            atoms.append(
                (
                    int(fields[0]),
                    float(fields[2]),
                    float(fields[3]),
                    float(fields[4]),
                )
            )
        except ValueError:
            continue

    return atoms


def discover_layers(directory):
    pattern = re.compile(
        rf"{re.escape(OUTPUT_PREFIX)}\.anchor_\d+\.root_layer_(\d+)\.mol2$"
    )

    layers = {}

    for path in directory.glob(f"{OUTPUT_PREFIX}.anchor_*.root_layer_*.mol2"):
        match = pattern.match(path.name)

        if match:
            layers.setdefault(int(match.group(1)), []).append(path)

    return layers


def discover_growth_trees(directory):
    pattern = re.compile(rf"{re.escape(OUTPUT_PREFIX)}\.growth_tree_(\d+)\.mol2$")

    trees = {}

    for path in directory.glob(f"{OUTPUT_PREFIX}.growth_tree_*.mol2"):
        match = pattern.match(path.name)

        if match:
            trees[int(match.group(1))] = path

    return trees


def layer_statistics():
    baseline_layers = discover_layers(BASELINE)
    modularized_layers = discover_layers(MODULARIZED)
    rows = []

    for layer in sorted(set(baseline_layers) | set(modularized_layers)):
        baseline_count = sum(molecule_count(path) for path in baseline_layers.get(layer, []))
        modularized_count = sum(molecule_count(path) for path in modularized_layers.get(layer, []))

        rows.append({
            "layer": layer,
            "baseline": baseline_count,
            "modularized": modularized_count,
            "difference": modularized_count - baseline_count,
        })

    return rows


def growth_tree_statistics():
    baseline_trees = discover_growth_trees(BASELINE)
    modularized_trees = discover_growth_trees(MODULARIZED)
    rows = []

    for tree_id in sorted(set(baseline_trees) | set(modularized_trees)):
        baseline_count = molecule_count(baseline_trees[tree_id]) if tree_id in baseline_trees else 0
        modularized_count = molecule_count(modularized_trees[tree_id]) if tree_id in modularized_trees else 0

        rows.append({
            "tree": tree_id,
            "baseline": baseline_count,
            "modularized": modularized_count,
            "difference": modularized_count - baseline_count,
        })

    return rows


def coordinate_statistics():
    baseline_file = BASELINE / FINAL_MOL2
    modularized_file = MODULARIZED / FINAL_MOL2

    baseline_molecules = split_mol2(baseline_file)
    modularized_molecules = split_mol2(modularized_file)

    if len(baseline_molecules) != len(modularized_molecules):
        raise ValueError("Final molecule counts differ; coordinate analysis is not aligned.")

    squared = []
    distances = []

    for baseline_mol, modularized_mol in zip(
        baseline_molecules,
        modularized_molecules,
    ):
        baseline_atoms = parse_atoms(baseline_mol)
        modularized_atoms = parse_atoms(modularized_mol)

        if len(baseline_atoms) != len(modularized_atoms):
            raise ValueError("Final molecule atom counts differ; coordinate analysis is not aligned.")
        for atom_a, atom_b in zip(
            baseline_atoms,
            modularized_atoms,
        ):
            if atom_a[0] != atom_b[0]:
                raise ValueError("Final molecule atom IDs differ; coordinate analysis is not aligned.")

            dx = atom_b[1] - atom_a[1]
            dy = atom_b[2] - atom_a[2]
            dz = atom_b[3] - atom_a[3]

            d2 = dx * dx + dy * dy + dz * dz

            squared.append(d2)
            distances.append(math.sqrt(d2))

    if not distances:
        return {
            "atoms": 0,
            "mean": math.nan,
            "maximum": math.nan,
            "rmsd": math.nan,
        }

    return {
        "atoms": len(distances),
        "mean": sum(distances) / len(distances),
        "maximum": max(distances),
        "rmsd": math.sqrt(sum(squared) / len(squared)),
    }


def common_numeric_properties():
    baseline_file = BASELINE / FINAL_MOL2
    modularized_file = MODULARIZED / FINAL_MOL2

    baseline_molecules = split_mol2(baseline_file)
    modularized_molecules = split_mol2(modularized_file)

    if len(baseline_molecules) != len(modularized_molecules):
        raise ValueError("Final molecule counts differ; score analysis is not aligned.")

    property_pairs = {}

    for baseline_mol, modularized_mol in zip(
        baseline_molecules,
        modularized_molecules,
    ):
        props_a = parse_properties(baseline_mol)
        props_b = parse_properties(modularized_mol)

        for key in sorted(set(props_a) & set(props_b)):
            values_a = props_a[key]
            values_b = props_b[key]

            if len(values_a) != len(values_b):
                continue

            for value_a, value_b in zip(values_a, values_b):
                number_a = numeric_value(value_a)
                number_b = numeric_value(value_b)

                if number_a is None or number_b is None:
                    continue

                property_pairs.setdefault(key, []).append(
                    (number_a, number_b)
                )

    return property_pairs


def select_score_property(property_pairs):
    if not property_pairs:
        return None, []

    preferred = [
        "Descriptor_Score",
        "descriptor_score",
        "Grid_Score",
        "grid_score",
        "HMS_Score",
        "hms_score",
        "Score",
        "score",
    ]

    for preferred_key in preferred:
        for key, values in property_pairs.items():
            if key.lower() == preferred_key.lower() and values:
                return key, values

    score_candidates = [
        (key, values)
        for key, values in property_pairs.items()
        if "score" in key.lower() and values
    ]

    if score_candidates:
        score_candidates.sort(
            key=lambda item: len(item[1]),
            reverse=True,
        )

        return score_candidates[0]

    candidates = [
        (key, values)
        for key, values in property_pairs.items()
        if values
    ]

    if not candidates:
        return None, []

    candidates.sort(
        key=lambda item: len(item[1]),
        reverse=True,
    )

    return candidates[0]


def score_statistics(values):
    if not values:
        return {
            "n": 0,
            "mae": math.nan,
            "rmse": math.nan,
            "maximum": math.nan,
        }

    differences = [
        modularized - baseline
        for baseline, modularized in values
    ]

    absolute = [abs(value) for value in differences]

    return {
        "n": len(values),
        "mae": sum(absolute) / len(absolute),
        "rmse": math.sqrt(
            sum(value * value for value in differences)
            / len(differences)
        ),
        "maximum": max(absolute),
    }


def save_statistics(
    layer_rows,
    tree_rows,
    coordinate_stats,
    score_key,
    score_stats,
):
    max_layer_difference = max(
        (
            abs(row["difference"])
            for row in layer_rows
        ),
        default=0,
    )

    max_tree_difference = max(
        (
            abs(row["difference"])
            for row in tree_rows
        ),
        default=0,
    )

    final_baseline = molecule_count(
        BASELINE / FINAL_MOL2
    )

    final_modularized = molecule_count(
        MODULARIZED / FINAL_MOL2
    )

    rows = [
        ("final_molecules_baseline", final_baseline),
        ("final_molecules_modularized", final_modularized),
        ("growth_layers_compared", len(layer_rows)),
        ("max_growth_layer_population_difference",
         max_layer_difference),
        ("growth_trees_compared", len(tree_rows)),
        ("max_growth_tree_size_difference",
         max_tree_difference),
        ("atoms_compared", coordinate_stats["atoms"]),
        ("mean_coordinate_displacement_A",
         coordinate_stats["mean"]),
        ("maximum_coordinate_displacement_A",
         coordinate_stats["maximum"]),
        ("coordinate_RMSD_A",
         coordinate_stats["rmsd"]),
        ("score_property", score_key or "not_found"),
        ("scores_compared", score_stats["n"]),
        ("score_MAE", score_stats["mae"]),
        ("score_RMSE", score_stats["rmse"]),
        ("maximum_absolute_score_difference",
         score_stats["maximum"]),
    ]

    with STATISTICS.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["metric", "value"])
        writer.writerows(rows)


def plot_growth_layers(rows):
    layers = [row["layer"] for row in rows]
    baseline = [row["baseline"] for row in rows]
    modularized = [row["modularized"] for row in rows]

    plt.figure(figsize=(8, 5))

    plt.plot(
        layers,
        baseline,
        marker="o",
        linewidth=2,
        label="Baseline DRef",
    )

    plt.plot(
        layers,
        modularized,
        marker="x",
        linestyle="--",
        linewidth=1.5,
        label="Modularized DRef",
    )

    plt.xlabel("Growth Layer")
    plt.ylabel("Number of Molecules")
    plt.title("DRef Growth-Layer Population Equivalence")
    plt.xticks(layers)
    plt.legend()
    plt.tight_layout()

    plt.savefig(
        RESULTS / "growth_layer_population.png",
        dpi=300,
    )

    plt.close()


def plot_growth_tree_identity(rows):
    baseline = [row["baseline"] for row in rows]
    modularized = [row["modularized"] for row in rows]

    plt.figure(figsize=(6, 6))

    plt.scatter(
        baseline,
        modularized,
        s=28,
        alpha=0.65,
    )

    if baseline and modularized:
        minimum = min(baseline + modularized)
        maximum = max(baseline + modularized)

        plt.plot(
            [minimum, maximum],
            [minimum, maximum],
            linestyle="--",
            linewidth=1.5,
            label="Identity line",
        )

    plt.xlabel("Baseline Growth-Tree Size")
    plt.ylabel("Modularized Growth-Tree Size")
    plt.title("Growth-Tree Identity")
    plt.legend()
    plt.tight_layout()

    plt.savefig(
        RESULTS / "growth_tree_identity.png",
        dpi=300,
    )

    plt.close()


def plot_growth_tree_distribution(rows):
    baseline = [row["baseline"] for row in rows]
    modularized = [row["modularized"] for row in rows]

    unique_sizes = sorted(
        set(baseline) | set(modularized)
    )

    baseline_frequency = [
        baseline.count(size)
        for size in unique_sizes
    ]

    modularized_frequency = [
        modularized.count(size)
        for size in unique_sizes
    ]

    positions = list(range(len(unique_sizes)))
    width = 0.38

    plt.figure(figsize=(9, 5))

    plt.bar(
        [position - width / 2 for position in positions],
        baseline_frequency,
        width=width,
        label="Baseline DRef",
    )

    plt.bar(
        [position + width / 2 for position in positions],
        modularized_frequency,
        width=width,
        label="Modularized DRef",
    )

    plt.xlabel("Growth-Tree Size")
    plt.ylabel("Number of Growth Trees")
    plt.title("Growth-Tree Size Distribution")
    plt.xticks(positions, unique_sizes)
    plt.legend()
    plt.tight_layout()

    plt.savefig(
        RESULTS / "growth_tree_distribution.png",
        dpi=300,
    )

    plt.close()


def plot_score_identity(score_key, values):
    baseline = [value[0] for value in values]
    modularized = [value[1] for value in values]

    plt.figure(figsize=(6, 6))

    plt.scatter(
        baseline,
        modularized,
        s=25,
        alpha=0.65,
    )

    if baseline and modularized:
        minimum = min(baseline + modularized)
        maximum = max(baseline + modularized)

        plt.plot(
            [minimum, maximum],
            [minimum, maximum],
            linestyle="--",
            linewidth=1.5,
            label="Identity line",
        )

    plt.xlabel(f"Baseline {score_key}")
    plt.ylabel(f"Modularized {score_key}")
    plt.title("Score Identity")
    plt.legend()
    plt.tight_layout()

    plt.savefig(
        RESULTS / "score_identity.png",
        dpi=300,
    )

    plt.close()


def plot_score_difference(score_key, values):
    differences = [
        modularized - baseline
        for baseline, modularized in values
    ]

    molecule_ids = list(range(1, len(differences) + 1))

    plt.figure(figsize=(9, 5))

    plt.scatter(
        molecule_ids,
        differences,
        s=20,
        alpha=0.65,
    )

    plt.axhline(
        0.0,
        linestyle="--",
        linewidth=1.5,
    )

    plt.xlabel("Score Observation")
    plt.ylabel("Modularized - Baseline")
    plt.title(f"{score_key} Difference")
    plt.tight_layout()

    plt.savefig(
        RESULTS / "score_difference.png",
        dpi=300,
    )

    plt.close()


def main():
    RESULTS.mkdir(parents=True, exist_ok=True)

    final_baseline = BASELINE / FINAL_MOL2
    final_modularized = MODULARIZED / FINAL_MOL2

    if not final_baseline.exists():
        print(
            f"ERROR: Missing baseline {FINAL_MOL2}",
            file=sys.stderr,
        )
        return 1

    if not final_modularized.exists():
        print(
            f"ERROR: Missing modularized {FINAL_MOL2}",
            file=sys.stderr,
        )
        return 1

    layer_rows = layer_statistics()
    tree_rows = growth_tree_statistics()
    coordinate_stats = coordinate_statistics()

    property_pairs = common_numeric_properties()
    score_key, score_values = select_score_property(
        property_pairs
    )

    score_stats = score_statistics(score_values)

    save_statistics(
        layer_rows,
        tree_rows,
        coordinate_stats,
        score_key,
        score_stats,
    )

    plot_growth_layers(layer_rows)
    plot_growth_tree_identity(tree_rows)
    plot_growth_tree_distribution(tree_rows)

    if score_key and score_values:
        plot_score_identity(
            score_key,
            score_values,
        )

        plot_score_difference(
            score_key,
            score_values,
        )
    else:
        print(
            "WARNING: No common numeric score property found."
        )

    print()
    print("DRef Equivalence Analysis")
    print()
    print(
        f"Growth layers compared:       {len(layer_rows)}"
    )
    print(
        f"Growth trees compared:        {len(tree_rows)}"
    )
    print(
        f"Atoms compared:               "
        f"{coordinate_stats['atoms']}"
    )
    print(
        f"Coordinate RMSD:              "
        f"{coordinate_stats['rmsd']:.10f} A"
    )
    print(
        f"Maximum coordinate difference:"
        f" {coordinate_stats['maximum']:.10f} A"
    )
    print(
        f"Score property:               "
        f"{score_key or 'not found'}"
    )
    print(
        f"Scores compared:              "
        f"{score_stats['n']}"
    )
    print(
        f"Score MAE:                    "
        f"{score_stats['mae']:.10f}"
    )
    print(
        f"Score RMSE:                   "
        f"{score_stats['rmse']:.10f}"
    )
    print(
        f"Maximum score difference:     "
        f"{score_stats['maximum']:.10f}"
    )
    print()
    print(f"Statistics: {STATISTICS}")
    print()
    print("Plots:")
    print(
        RESULTS / "growth_layer_population.png"
    )
    print(
        RESULTS / "growth_tree_identity.png"
    )
    print(
        RESULTS / "growth_tree_distribution.png"
    )

    if score_key and score_values:
        print(
            RESULTS / "score_identity.png"
        )
        print(
            RESULTS / "score_difference.png"
        )

    print()

    return 0


if __name__ == "__main__":
    sys.exit(main())
