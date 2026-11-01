#!/usr/bin/env python3
"""Compare each inherited/isolated system-anchor pair and summarize equivalence."""
from pathlib import Path
import csv
import hashlib
import math
import re
import sys

ROOT = Path(__file__).resolve().parent
RUNS = ROOT / 'runs'
RESULTS = ROOT / 'comparison_results'
TOL = 1e-6


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parameter(path, key):
    values = [f[1] for line in path.read_text().splitlines()
              if (f := line.split()) and len(f) > 1 and f[0] == key]
    if len(values) != 1:
        raise ValueError(f'Expected one {key}: {path}')
    return values[0]


def split_molecules(path):
    text = path.read_text(errors='replace')
    return ['@<TRIPOS>MOLECULE' + x for x in text.split('@<TRIPOS>MOLECULE')[1:] if x.strip()]


def section(molecule, name):
    marker = f'@<TRIPOS>{name}'
    if marker not in molecule:
        return []
    piece = molecule.split(marker, 1)[1]
    piece = re.split(r'\n@<TRIPOS>', piece, maxsplit=1)[0]
    return [line.split() for line in piece.splitlines() if line.strip()]


def properties(molecule):
    props = {}
    for line in molecule.splitlines():
        if line.lstrip().startswith('##########') and ':' in line:
            k, v = line.strip().lstrip('#').strip().split(':', 1)
            props.setdefault(k.strip(), []).append(v.strip())
    return props


def near(a, b):
    try:
        return math.isclose(float(a), float(b), abs_tol=TOL, rel_tol=TOL)
    except (TypeError, ValueError):
        return a == b


def compare_molecule(a, b):
    errors, distances, score_differences = [], [], []
    name_a = a.splitlines()[1].strip() if len(a.splitlines()) > 1 else ''
    name_b = b.splitlines()[1].strip() if len(b.splitlines()) > 1 else ''
    if name_a != name_b:
        errors.append('molecule name')
    aa, ab = section(a, 'ATOM'), section(b, 'ATOM')
    if len(aa) != len(ab):
        errors.append('atom count')
    for x, y in zip(aa, ab):
        if len(x) < 6 or len(y) < 6:
            errors.append('malformed atom record')
            continue
        if x[0:2] != y[0:2] or x[5] != y[5] or (x[6:8] != y[6:8]):
            errors.append('atom identity/type')
        for i in (2, 3, 4):
            if not near(x[i], y[i]):
                errors.append('coordinate')
        d2 = sum((float(x[i]) - float(y[i])) ** 2 for i in (2, 3, 4))
        distances.append(d2)
        if len(x) > 8 or len(y) > 8:
            if len(x) <= 8 or len(y) <= 8 or not near(x[8], y[8]):
                errors.append('charge')
    if section(a, 'BOND') != section(b, 'BOND'):
        errors.append('bond definitions')
    pa, pb = properties(a), properties(b)
    if pa.keys() != pb.keys():
        errors.append('property keys')
    for key in pa.keys() & pb.keys():
        if len(pa[key]) != len(pb[key]):
            errors.append('property count: ' + key)
            continue
        for x, y in zip(pa[key], pb[key]):
            if not near(x, y):
                errors.append('property: ' + key)
            if 'score' in key.lower():
                try:
                    score_differences.append(float(y) - float(x))
                except ValueError:
                    pass
    return errors, distances, score_differences


def compare_pair(system, anchor):
    left = RUNS / 'inherited' / system / f'anc_{anchor}'
    right = RUNS / 'isolated' / system / f'anc_{anchor}'
    outdir = RESULTS / system / f'anc_{anchor}'
    outdir.mkdir(parents=True, exist_ok=True)
    errors, distances, score_differences = [], [], []
    if not (left / '.success').is_file() or not (right / '.success').is_file():
        errors.append('One or both runs lack a success marker')
    if not (left / 'dn_generic.in').is_file() or not (right / 'dn_generic.in').is_file():
        errors.append('Input missing')
        prefix = None
    else:
        if sha(left / 'dn_generic.in') != sha(right / 'dn_generic.in'):
            errors.append('Inputs have different hashes')
        prefix = parameter(left / 'dn_generic.in', 'dn_output_prefix')
    files_a = {p.name: p for p in left.glob('*.mol2') if p.is_file()}
    files_b = {p.name: p for p in right.glob('*.mol2') if p.is_file()}
    if set(files_a) != set(files_b):
        errors.append('MOL2 filename sets differ: inherited-only=' +
                      str(sorted(set(files_a) - set(files_b))) + '; isolated-only=' +
                      str(sorted(set(files_b) - set(files_a))))
    common = sorted(set(files_a) & set(files_b))
    if not common:
        errors.append('No common MOL2 results')
    final_names = [n for n in common if 'denovo_build' in n and n.endswith('.mol2')]
    if not final_names:
        errors.append('No common final denovo_build MOL2 file')
    if prefix and not any(n.startswith(prefix) for n in common):
        errors.append('No MOL2 output uses configured dn_output_prefix')
    growth = []
    tree = []
    final_a = final_b = 0
    for name in common:
        molecules_a = split_molecules(files_a[name])
        molecules_b = split_molecules(files_b[name])
        if len(molecules_a) != len(molecules_b):
            errors.append(f'{name}: molecule count {len(molecules_a)} != {len(molecules_b)}')
        if 'root_layer_' in name:
            growth.append((name, len(molecules_a), len(molecules_b)))
        if 'growth_tree_' in name:
            tree.append((name, len(molecules_a), len(molecules_b)))
        if name in final_names:
            final_a += len(molecules_a)
            final_b += len(molecules_b)
        for index, (ma, mb) in enumerate(zip(molecules_a, molecules_b), 1):
            faults, ds, ss = compare_molecule(ma, mb)
            if faults:
                errors.append(f'{name}: molecule {index}: {", ".join(sorted(set(faults)))}')
            if name in final_names:
                distances.extend(ds)
                score_differences.extend(ss)
    if not growth:
        errors.append('No growth-layer files found')
    if not tree:
        errors.append('No growth-tree files found')
    if not score_differences:
        errors.append('No comparable numerical score properties found')
    rmsd = math.sqrt(sum(distances) / len(distances)) if distances else math.nan
    mae = sum(abs(d) for d in score_differences) / len(score_differences) if score_differences else math.nan
    rmse = math.sqrt(sum(d * d for d in score_differences) / len(score_differences)) if score_differences else math.nan
    result = 'PASS' if not errors else 'FAIL'
    report = [f'System: {system}', f'Anchor: {anchor}', f'Overall: {result}',
              f'Inherited final molecules: {final_a}', f'Isolated final molecules: {final_b}',
              f'Growth layer files: {len(growth)}', f'Growth tree files: {len(tree)}',
              f'Coordinate RMSD (A): {rmsd}', f'Score MAE: {mae}', f'Score RMSE: {rmse}',
              'Growth populations: ' + repr(growth), 'Growth trees: ' + repr(tree),
              'Differences:'] + (errors if errors else ['None'])
    (outdir / 'equivalence_report.txt').write_text('\n'.join(report) + '\n')
    return [system, anchor, result, final_a, final_b, len(growth), len(tree), rmsd, mae, rmse, len(errors)]


def main():
    tasks = RUNS / 'tasks.tsv'
    if not tasks.is_file():
        raise SystemExit('ERROR: Run 001.setup_test.sh first')
    RESULTS.mkdir(parents=True, exist_ok=True)
    rows = []
    for line in tasks.read_text().splitlines():
        system, anchor = line.split('\t')
        try:
            row = compare_pair(system, anchor)
        except Exception as exc:
            folder = RESULTS / system / f'anc_{anchor}'
            folder.mkdir(parents=True, exist_ok=True)
            (folder / 'equivalence_report.txt').write_text(f'FAIL: {exc}\n')
            row = [system, anchor, 'FAIL', '', '', '', '', '', '', '', str(exc)]
        rows.append(row)
        print(system, 'anc_' + anchor, row[2])
    with (RESULTS / 'equivalence_summary.csv').open('w', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(['system', 'anchor', 'equivalence', 'inherited_final_molecules',
                         'isolated_final_molecules', 'growth_layer_files', 'growth_tree_files',
                         'coordinate_rmsd_A', 'score_mae', 'score_rmse', 'differences'])
        writer.writerows(rows)
    print('Report:', RESULTS / 'equivalence_summary.csv')
    return int(any(row[2] != 'PASS' for row in rows) or not rows)


if __name__ == '__main__':
    sys.exit(main())
