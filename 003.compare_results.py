#!/usr/bin/env python3

from pathlib import Path
import hashlib
import math
import os
import re
import sys


ROOT = Path(__file__).resolve().parent

INHERITED = ROOT / "inherited" / "system"
ISOLATED = ROOT / "isolated" / "system"
RESULTS = ROOT / "comparison_results"

REPORT = RESULTS / "equivalence_report.txt"

COORD_TOL = 1.0e-6
CHARGE_TOL = 1.0e-6
NUMERIC_TOL = 1.0e-6

FINAL_MOL2 = "output.denovo_build.mol2"


def sha256(path):
    h = hashlib.sha256()

    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)

    return h.hexdigest()


def split_mol2(path):
    text = path.read_text(errors="replace")

    marker = "@<TRIPOS>MOLECULE"

    chunks = text.split(marker)

    return [
        marker + chunk
        for chunk in chunks[1:]
        if chunk.strip()
    ]


def section_lines(block, section):
    marker = f"@<TRIPOS>{section}"

    if marker not in block:
        return []

    part = block.split(marker, 1)[1]

    next_section = re.search(r"\n@<TRIPOS>", part)

    if next_section:
        part = part[:next_section.start()]

    return [
        line.strip()
        for line in part.splitlines()
        if line.strip()
    ]


def molecule_name(block):
    lines = block.splitlines()

    for i, line in enumerate(lines):
        if line.strip() == "@<TRIPOS>MOLECULE":
            if i + 1 < len(lines):
                return lines[i + 1].strip()

    return ""


def parse_atoms(block):
    atoms = []

    for line in section_lines(block, "ATOM"):
        fields = line.split()

        if len(fields) < 6:
            continue

        atom = {
            "id": int(fields[0]),
            "name": fields[1],
            "x": float(fields[2]),
            "y": float(fields[3]),
            "z": float(fields[4]),
            "type": fields[5],
            "charge": float(fields[8]) if len(fields) > 8 else None,
        }

        atoms.append(atom)

    return atoms


def parse_bonds(block):
    bonds = []

    for line in section_lines(block, "BOND"):
        fields = line.split()

        if len(fields) < 4:
            continue

        bonds.append(
            (
                int(fields[0]),
                int(fields[1]),
                int(fields[2]),
                fields[3],
            )
        )

    return bonds


def parse_properties(block):
    properties = {}

    for line in block.splitlines():
        stripped = line.strip()

        if not stripped.startswith("##########"):
            continue

        payload = stripped.lstrip("#").strip()

        if ":" in payload:
            key, value = payload.split(":", 1)
            properties.setdefault(key.strip(), []).append(value.strip())

    return properties


def compare_numeric_text(a, b, tol):
    try:
        av = float(a)
        bv = float(b)
    except ValueError:
        return a == b

    return math.isclose(av, bv, rel_tol=tol, abs_tol=tol)


def compare_properties(a, b):
    differences = []

    keys = sorted(set(a) | set(b))

    for key in keys:
        va = a.get(key)
        vb = b.get(key)

        if va is None or vb is None:
            differences.append(f"property missing: {key}")
            continue

        if len(va) != len(vb):
            differences.append(
                f"property count differs: {key} "
                f"({len(va)} vs {len(vb)})"
            )
            continue

        for index, (x, y) in enumerate(zip(va, vb), start=1):
            if not compare_numeric_text(x, y, NUMERIC_TOL):
                differences.append(
                    f"property differs: {key}[{index}] "
                    f"({x} vs {y})"
                )

    return differences


def compare_molecule(block_a, block_b, index):
    differences = []

    name_a = molecule_name(block_a)
    name_b = molecule_name(block_b)

    if name_a != name_b:
        differences.append(
            f"molecule {index}: name differs "
            f"({name_a!r} vs {name_b!r})"
        )

    atoms_a = parse_atoms(block_a)
    atoms_b = parse_atoms(block_b)

    if len(atoms_a) != len(atoms_b):
        differences.append(
            f"molecule {index}: atom count differs "
            f"({len(atoms_a)} vs {len(atoms_b)})"
        )
    else:
        for atom_a, atom_b in zip(atoms_a, atoms_b):

            if atom_a["id"] != atom_b["id"]:
                differences.append(
                    f"molecule {index}: atom ID differs "
                    f"({atom_a['id']} vs {atom_b['id']})"
                )

            if atom_a["name"] != atom_b["name"]:
                differences.append(
                    f"molecule {index} atom {atom_a['id']}: "
                    f"name differs"
                )

            if atom_a["type"] != atom_b["type"]:
                differences.append(
                    f"molecule {index} atom {atom_a['id']}: "
                    f"type differs"
                )

            for axis in ("x", "y", "z"):
                if not math.isclose(
                    atom_a[axis],
                    atom_b[axis],
                    rel_tol=COORD_TOL,
                    abs_tol=COORD_TOL,
                ):
                    differences.append(
                        f"molecule {index} atom {atom_a['id']}: "
                        f"{axis} differs "
                        f"({atom_a[axis]} vs {atom_b[axis]})"
                    )

            qa = atom_a["charge"]
            qb = atom_b["charge"]

            if qa is None or qb is None:
                if qa != qb:
                    differences.append(
                        f"molecule {index} atom {atom_a['id']}: "
                        "charge field differs"
                    )
            elif not math.isclose(
                qa,
                qb,
                rel_tol=CHARGE_TOL,
                abs_tol=CHARGE_TOL,
            ):
                differences.append(
                    f"molecule {index} atom {atom_a['id']}: "
                    f"charge differs ({qa} vs {qb})"
                )

    bonds_a = parse_bonds(block_a)
    bonds_b = parse_bonds(block_b)

    if bonds_a != bonds_b:
        differences.append(
            f"molecule {index}: bond definitions differ"
        )

    properties_a = parse_properties(block_a)
    properties_b = parse_properties(block_b)

    for item in compare_properties(properties_a, properties_b):
        differences.append(f"molecule {index}: {item}")

    return differences


def compare_mol2(path_a, path_b):
    molecules_a = split_mol2(path_a)
    molecules_b = split_mol2(path_b)

    differences = []

    if len(molecules_a) != len(molecules_b):
        differences.append(
            f"molecule count differs: "
            f"{len(molecules_a)} vs {len(molecules_b)}"
        )
        return differences, len(molecules_a), len(molecules_b)

    for index, (mol_a, mol_b) in enumerate(
        zip(molecules_a, molecules_b),
        start=1,
    ):
        differences.extend(
            compare_molecule(mol_a, mol_b, index)
        )

    return differences, len(molecules_a), len(molecules_b)


def generated_mol2_files(directory):
    return {
        path.name: path
        for path in directory.glob("output*.mol2")
        if path.is_file()
    }


def write(line="", report=None):
    print(line)

    if report is not None:
        report.write(line + "\n")


def main():
    RESULTS.mkdir(parents=True, exist_ok=True)

    failures = []
    warnings = []

    with REPORT.open("w") as report:

        write("DRef Equivalence Test", report)
        write("", report)

        inherited_success = (INHERITED / ".success").exists()
        isolated_success = (ISOLATED / ".success").exists()

        write(
            f"Inherited execution: "
            f"{'PASS' if inherited_success else 'FAIL'}",
            report,
        )

        write(
            f"Isolated execution:  "
            f"{'PASS' if isolated_success else 'FAIL'}",
            report,
        )

        if not inherited_success:
            failures.append("Inherited calculation did not complete.")

        if not isolated_success:
            failures.append("Isolated calculation did not complete.")

        input_a = INHERITED / "dn_generic.in"
        input_b = ISOLATED / "dn_generic.in"

        if input_a.exists() and input_b.exists():

            hash_a = sha256(input_a)
            hash_b = sha256(input_b)

            input_equal = hash_a == hash_b

            write("", report)
            write(
                f"Input files:         "
                f"{'PASS' if input_equal else 'FAIL'}",
                report,
            )
            write(f"Inherited SHA256:   {hash_a}", report)
            write(f"Isolated SHA256:    {hash_b}", report)

            if not input_equal:
                failures.append(
                    "Inherited and isolated input files differ."
                )

        else:
            failures.append("One or both input files are missing.")

        inherited_files = generated_mol2_files(INHERITED)
        isolated_files = generated_mol2_files(ISOLATED)

        names_a = set(inherited_files)
        names_b = set(isolated_files)

        write("", report)
        write(
            f"Generated MOL2 files: "
            f"{len(names_a)} inherited, "
            f"{len(names_b)} isolated",
            report,
        )

        missing_isolated = sorted(names_a - names_b)
        missing_inherited = sorted(names_b - names_a)

        if missing_isolated:
            failures.append(
                "Missing from isolated: "
                + ", ".join(missing_isolated)
            )

        if missing_inherited:
            failures.append(
                "Missing from inherited: "
                + ", ".join(missing_inherited)
            )

        common_files = sorted(names_a & names_b)

        structural_failures = 0
        exact_matches = 0

        write("", report)
        write("MOL2 comparison", report)

        for name in common_files:

            path_a = inherited_files[name]
            path_b = isolated_files[name]

            exact = sha256(path_a) == sha256(path_b)

            differences, count_a, count_b = compare_mol2(
                path_a,
                path_b,
            )

            structural_equal = len(differences) == 0

            if exact:
                exact_matches += 1

            if not structural_equal:
                structural_failures += 1

            write(
                f"{name}: "
                f"molecules={count_a}/{count_b}, "
                f"structural="
                f"{'PASS' if structural_equal else 'FAIL'}, "
                f"byte="
                f"{'PASS' if exact else 'DIFFERENT'}",
                report,
            )

            if differences:
                for difference in differences[:25]:
                    write(f"  {difference}", report)

                if len(differences) > 25:
                    write(
                        f"  ... {len(differences) - 25} "
                        f"additional differences",
                        report,
                    )

        final_a = INHERITED / FINAL_MOL2
        final_b = ISOLATED / FINAL_MOL2

        if not final_a.exists():
            failures.append(
                f"Inherited final MOL2 missing: {FINAL_MOL2}"
            )

        if not final_b.exists():
            failures.append(
                f"Isolated final MOL2 missing: {FINAL_MOL2}"
            )

        if structural_failures:
            failures.append(
                f"{structural_failures} common MOL2 file(s) "
                "failed structural equivalence."
            )

        if common_files and exact_matches != len(common_files):
            warnings.append(
                f"{len(common_files) - exact_matches} common MOL2 "
                "file(s) were structurally equivalent but not "
                "byte-identical."
            )

        write("", report)
        write("Summary", report)
        write(
            f"Common MOL2 files:    {len(common_files)}",
            report,
        )
        write(
            f"Byte-identical MOL2:  {exact_matches}",
            report,
        )
        write(
            f"Structural failures:  {structural_failures}",
            report,
        )

        if warnings:
            write("", report)
            write("Warnings", report)

            for warning in warnings:
                write(f"- {warning}", report)

        write("", report)

        if failures:
            write("OVERALL RESULT: FAIL", report)
            write("", report)

            for failure in failures:
                write(f"- {failure}", report)

            result = 1

        else:
            write("OVERALL RESULT: PASS", report)
            write(
                "The isolated implementation is structurally "
                "equivalent to the inherited implementation "
                "for this test case.",
                report,
            )

            result = 0

        write("", report)
        write(f"Report: {REPORT}", report)

    return result


if __name__ == "__main__":
    sys.exit(main())
