#!/usr/bin/env python3

from pathlib import Path
import csv
import math
import re
import sys

import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parent

INHERITED = ROOT / "inherited" / "system"
ISOLATED = ROOT / "isolated" / "system"
RESULTS = ROOT / "comparison_results"

STATISTICS = RESULTS / "equivalence_statistics.csv"

FINAL_MOL2 = "output.denovo_build.mol2"


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
        r"output\.anchor_\d+\.root_layer_(\d+)\.mol2$"
    )

    layers = {}

    for path in directory.glob("output.anchor_*.root_layer_*.mol2"):
        match = pattern.match(path.name)

        if match:
            layers[int(match.group(1))] = path

    return layers


def discover_growth_trees(directory):
    pattern = re.compile(r"output\.growth_tree_(\d+)\.mol2$")

    trees = {}

    for path in directory.glob("output.growth_tree_*.mol2"):
        match = pattern.match(path.name)

        if match:
            trees[int(match.group(1))] = path

    return trees


def layer_statistics():
    inherited_layers = discover_layers(INHERITED)
    isolated_layers = discover_layers(ISOLATED)

    all_layers = sorted(
        set(inherited_layers) | set(isolated_layers)
    )

    rows = []

    for layer in all_layers:
        inherited_file = inherited_layers.get(layer)
        isolated_file = isolated_layers.get(layer)

        inherited_count = (
            molecule_count(inherited_file)
            if inherited_file
            else 0
        )

        isolated_count = (
            molecule_count(isolated_file)
            if isolated_file
            else 0
        )

        rows.append(
            {
                "layer": layer,
                "inherited": inherited_count,
                "isolated": isolated_count,
                "difference": isolated_count - inherited_count,
            }
        )

    return rows


def growth_tree_statistics():
    inherited_trees = discover_growth_trees(INHERITED)
    isolated_trees = discover_growth_trees(ISOLATED)

    common_ids = sorted(
        set(inherited_trees) & set(isolated_trees)
    )

    rows = []

    for tree_id in common_ids:
        inherited_count = molecule_count(
            inherited_trees[tree_id]
        )

        isolated_count = molecule_count(
            isolated_trees[tree_id]
        )

        rows.append(
            {
                "tree": tree_id,
                "inherited": inherited_count,
                "isolated": isolated_count,
                "difference": isolated_count - inherited_count,
            }
        )

    return rows


def coordinate_statistics():
    inherited_file = INHERITED / FINAL_MOL2
    isolated_file = ISOLATED / FINAL_MOL2

    inherited_molecules = split_mol2(inherited_file)
    isolated_molecules = split_mol2(isolated_file)

    squared = []
    distances = []

    for inherited_mol, isolated_mol in zip(
        inherited_molecules,
        isolated_molecules,
    ):
        inherited_atoms = parse_atoms(inherited_mol)
        isolated_atoms = parse_atoms(isolated_mol)

        for atom_a, atom_b in zip(
            inherited_atoms,
            isolated_atoms,
        ):
            if atom_a[0] != atom_b[0]:
                continue

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
    inherited_file = INHERITED / FINAL_MOL2
    isolated_file = ISOLATED / FINAL_MOL2

    inherited_molecules = split_mol2(inherited_file)
    isolated_molecules = split_mol2(isolated_file)

    property_pairs = {}

    for inherited_mol, isolated_mol in zip(
        inherited_molecules,
        isolated_molecules,
    ):
        props_a = parse_properties(inherited_mol)
        props_b = parse_properties(isolated_mol)

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
        isolated - inherited
        for inherited, isolated in values
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

    final_inherited = molecule_count(
        INHERITED / FINAL_MOL2
    )

    final_isolated = molecule_count(
        ISOLATED / FINAL_MOL2
    )

    rows = [
        ("final_molecules_inherited", final_inherited),
        ("final_molecules_isolated", final_isolated),
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
    inherited = [row["inherited"] for row in rows]
    isolated = [row["isolated"] for row in rows]

    plt.figure(figsize=(8, 5))

    plt.plot(
        layers,
        inherited,
        marker="o",
        linewidth=2,
        label="Inherited DRef",
    )

    plt.plot(
        layers,
        isolated,
        marker="x",
        linestyle="--",
        linewidth=1.5,
        label="Isolated DRef",
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
    inherited = [row["inherited"] for row in rows]
    isolated = [row["isolated"] for row in rows]

    plt.figure(figsize=(6, 6))

    plt.scatter(
        inherited,
        isolated,
        s=28,
        alpha=0.65,
    )

    if inherited and isolated:
        minimum = min(inherited + isolated)
        maximum = max(inherited + isolated)

        plt.plot(
            [minimum, maximum],
            [minimum, maximum],
            linestyle="--",
            linewidth=1.5,
            label="Identity line",
        )

    plt.xlabel("Inherited Growth-Tree Size")
    plt.ylabel("Isolated Growth-Tree Size")
    plt.title("Growth-Tree Identity")
    plt.legend()
    plt.tight_layout()

    plt.savefig(
        RESULTS / "growth_tree_identity.png",
        dpi=300,
    )

    plt.close()


def plot_growth_tree_distribution(rows):
    inherited = [row["inherited"] for row in rows]
    isolated = [row["isolated"] for row in rows]

    unique_sizes = sorted(
        set(inherited) | set(isolated)
    )

    inherited_frequency = [
        inherited.count(size)
        for size in unique_sizes
    ]

    isolated_frequency = [
        isolated.count(size)
        for size in unique_sizes
    ]

    positions = list(range(len(unique_sizes)))
    width = 0.38

    plt.figure(figsize=(9, 5))

    plt.bar(
        [position - width / 2 for position in positions],
        inherited_frequency,
        width=width,
        label="Inherited DRef",
    )

    plt.bar(
        [position + width / 2 for position in positions],
        isolated_frequency,
        width=width,
        label="Isolated DRef",
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
    inherited = [value[0] for value in values]
    isolated = [value[1] for value in values]

    plt.figure(figsize=(6, 6))

    plt.scatter(
        inherited,
        isolated,
        s=25,
        alpha=0.65,
    )

    if inherited and isolated:
        minimum = min(inherited + isolated)
        maximum = max(inherited + isolated)

        plt.plot(
            [minimum, maximum],
            [minimum, maximum],
            linestyle="--",
            linewidth=1.5,
            label="Identity line",
        )

    plt.xlabel(f"Inherited {score_key}")
    plt.ylabel(f"Isolated {score_key}")
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
        isolated - inherited
        for inherited, isolated in values
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
    plt.ylabel("Isolated - Inherited")
    plt.title(f"{score_key} Difference")
    plt.tight_layout()

    plt.savefig(
        RESULTS / "score_difference.png",
        dpi=300,
    )

    plt.close()


def main():
    RESULTS.mkdir(parents=True, exist_ok=True)

    final_inherited = INHERITED / FINAL_MOL2
    final_isolated = ISOLATED / FINAL_MOL2

    if not final_inherited.exists():
        print(
            f"ERROR: Missing inherited {FINAL_MOL2}",
            file=sys.stderr,
        )
        return 1

    if not final_isolated.exists():
        print(
            f"ERROR: Missing isolated {FINAL_MOL2}",
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
