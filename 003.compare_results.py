#!/usr/bin/env python3
"""Compare inherited/isolated DOCK de novo outputs for each matched system-anchor."""
from pathlib import Path
from itertools import zip_longest
import csv
import hashlib
import math
import re
import sys

ROOT = Path(__file__).resolve().parent
RUNS = ROOT / 'runs'
RESULTS = ROOT / 'comparison_results'
TOL = 1e-6


def parameter(path, key):
    values = [parts[1] for line in path.read_text().splitlines()
              if (parts := line.split()) and len(parts) > 1 and parts[0] == key]
    if len(values) != 1:
        raise ValueError(f'Expected exactly one {key} in {path}')
    return values[0]


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def mol_blocks(path):
    """Stream molecules without loading the complete MOL2 output into memory."""
    block = []
    with path.open(errors='replace') as handle:
        for line in handle:
            if line.startswith('@<TRIPOS>MOLECULE'):
                if block:
                    yield ''.join(block)
                block = [line]
            elif block:
                block.append(line)
        if block:
            yield ''.join(block)


def section(block, kind):
    marker = '@<TRIPOS>' + kind
    if marker not in block:
        return []
    remainder = block.split(marker, 1)[1]
    remainder = re.split(r'\n@<TRIPOS>', remainder, maxsplit=1)[0]
    return [line.split() for line in remainder.splitlines() if line.strip() and not line.startswith('#')]


def properties(block):
    data = {}
    for line in block.splitlines():
        if line.lstrip().startswith('##########') and ':' in line:
            key, value = line.strip().lstrip('#').strip().split(':', 1)
            data.setdefault(key.strip(), []).append(value.strip())
    return data


def equal(a, b):
    try:
        x, y = float(a), float(b)
        return math.isfinite(x) and math.isfinite(y) and math.isclose(x, y, rel_tol=TOL, abs_tol=TOL)
    except ValueError:
        return a == b


def inspect_molecule(a, b, score_pairs=None, coord_squares=None):
    faults = set()
    al, bl = a.splitlines(), b.splitlines()
    if (al[1].strip() if len(al) > 1 else '') != (bl[1].strip() if len(bl) > 1 else ''):
        faults.add('name')
    aa, ba = section(a, 'ATOM'), section(b, 'ATOM')
    if len(aa) != len(ba):
        faults.add('atom count')
    for x, y in zip(aa, ba):
        if len(x) < 6 or len(y) < 6:
            faults.add('malformed atom')
            continue
        if x[:2] != y[:2] or x[5:8] != y[5:8]:
            faults.add('atom identity/type')
        try:
            delta = [float(x[i]) - float(y[i]) for i in (2, 3, 4)]
            if coord_squares is not None:
                coord_squares.append(sum(v*v for v in delta))
            if any(not equal(x[i], y[i]) for i in (2, 3, 4)):
                faults.add('coordinates')
        except (ValueError, IndexError):
            faults.add('malformed coordinates')
        if (x[8] if len(x) > 8 else None) != (y[8] if len(y) > 8 else None):
            if len(x) <= 8 or len(y) <= 8 or not equal(x[8], y[8]):
                faults.add('charge')
    if section(a, 'BOND') != section(b, 'BOND'):
        faults.add('bonds')
    pa, pb = properties(a), properties(b)
    if set(pa) != set(pb):
        faults.add('property keys')
    for key in set(pa) & set(pb):
        if len(pa[key]) != len(pb[key]):
            faults.add('property count')
            continue
        for av, bv in zip(pa[key], pb[key]):
            if not equal(av, bv):
                faults.add('property values')
            if score_pairs is not None and 'score' in key.lower():
                try:
                    x, y = float(av), float(bv)
                    if math.isfinite(x) and math.isfinite(y):
                        score_pairs.setdefault(key, []).append((x, y))
                except ValueError:
                    pass
    return faults


def compare_final(a, b):
    total_a = total_b = 0
    mismatches = 0
    examples = []
    scores = {}
    coord_squares = []
    for index, (ma, mb) in enumerate(zip_longest(mol_blocks(a), mol_blocks(b)), 1):
        if ma is not None:
            total_a += 1
        if mb is not None:
            total_b += 1
        if ma is None or mb is None:
            mismatches += 1
            continue
        faults = inspect_molecule(ma, mb, scores, coord_squares)
        if faults:
            mismatches += 1
            if len(examples) < 20:
                examples.append(f'Molecule {index}: {", ".join(sorted(faults))}')
    rmsd = math.sqrt(sum(coord_squares) / len(coord_squares)) if coord_squares else math.nan
    return total_a, total_b, mismatches, examples, scores, rmsd


def compare_tree(a, b):
    """Accept numeric-tolerance equivalence when files differ byte-for-byte."""
    if digest(a) == digest(b):
        return True
    for ma, mb in zip_longest(mol_blocks(a), mol_blocks(b)):
        if ma is None or mb is None or inspect_molecule(ma, mb):
            return False
    return True


def molecule_count(path):
    count = 0
    with path.open('rb') as handle:
        for line in handle:
            if line.startswith(b'@<TRIPOS>MOLECULE'):
                count += 1
    return count


def scatter(path, x, y, title, xlabel, ylabel):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(x, y, s=8, alpha=0.45)
    lo, hi = min(x+y), max(x+y)
    if lo == hi:
        lo, hi = lo-1, hi+1
    ax.plot([lo, hi], [lo, hi], '--', linewidth=1)
    ax.set(xlabel=xlabel, ylabel=ylabel, title=title)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def compare_pair(system, anchor):
    a = RUNS / 'inherited' / system / f'anc_{anchor}'
    b = RUNS / 'isolated' / system / f'anc_{anchor}'
    out = RESULTS / system / f'anc_{anchor}'
    out.mkdir(parents=True, exist_ok=True)
    errors, notes = [], []
    ia, ib = a / 'dn_generic.in', b / 'dn_generic.in'
    if not ia.is_file() or not ib.is_file():
        raise RuntimeError('Missing one or both calculation inputs')
    if digest(ia) != digest(ib):
        errors.append('Input files differ')
    prefix = parameter(ia, 'dn_output_prefix')
    fa = a / f'{prefix}.completed.denovo_build.mol2'
    fb = b / f'{prefix}.completed.denovo_build.mol2'
    if not fa.is_file() or not fb.is_file():
        raise RuntimeError(f'Missing final output: {fa if not fa.is_file() else fb}')
    if not (a / '.success').exists() or not (b / '.success').exists():
        notes.append('One or both .success markers missing; inspect SLURM exit status')

    na, nb, bad, samples, scores, rmsd = compare_final(fa, fb)
    if na == 0 or nb == 0:
        errors.append('Final MOL2 contains no molecules')
    if na != nb:
        errors.append(f'Final molecule counts differ: {na} vs {nb}')
    if bad:
        errors.append(f'{bad} molecule records differ (coordinate/property tolerance {TOL})')
    if not scores:
        notes.append('No numeric score properties found in final molecules')

    # Compare growth-tree outputs by ID; a changed byte hash is flagged separately.
    pattern = re.compile(re.escape(prefix) + r'\.growth_tree_(\d+)\.mol2$')
    def trees(folder):
        ans = {}
        for p in folder.glob(f'{prefix}.growth_tree_*.mol2'):
            m = pattern.fullmatch(p.name)
            if m:
                ans[int(m.group(1))] = p
        return ans
    ta, tb = trees(a), trees(b)
    only_a, only_b = set(ta) - set(tb), set(tb) - set(ta)
    if only_a or only_b:
        errors.append(f'Growth-tree IDs differ: inherited-only={len(only_a)}, isolated-only={len(only_b)}')
    ids = sorted(set(ta) & set(tb))
    sizes_a, sizes_b = [], []
    differing_trees = 0
    for tree_id in ids:
        pa, pb = ta[tree_id], tb[tree_id]
        ca, cb = molecule_count(pa), molecule_count(pb)
        sizes_a.append(ca)
        sizes_b.append(cb)
        if ca != cb or not compare_tree(pa, pb):
            differing_trees += 1
    if differing_trees:
        errors.append(f'{differing_trees} growth-tree files differ structurally or in population')
    if not ta and not tb:
        notes.append('No growth-tree files; growth-tree equivalence not assessed')

    # Growth-layer files are optional; their absence must not create false failures.
    layer_a = {p.name for p in a.glob(f'{prefix}.*root_layer_*.mol2')}
    layer_b = {p.name for p in b.glob(f'{prefix}.*root_layer_*.mol2')}
    if layer_a != layer_b:
        errors.append('Growth-layer output filename sets differ')
    elif not layer_a:
        notes.append('No separate growth-layer MOL2 files were written')
    else:
        for name in sorted(layer_a):
            if not compare_tree(a / name, b / name):
                errors.append(f'Growth-layer file differs: {name}')

    if ids:
        scatter(out / 'growth_tree_identity.png', sizes_a, sizes_b,
                'Growth-tree population identity', 'Inherited molecules', 'Isolated molecules')
    for key, values in sorted(scores.items()):
        safe = re.sub(r'[^A-Za-z0-9_-]+', '_', key)
        if values:
            scatter(out / f'score_identity_{safe}.png', [x for x, _ in values],
                    [y for _, y in values], f'{key} identity', 'Inherited score', 'Isolated score')
    result = 'PASS' if not errors else 'FAIL'
    report = [f'System: {system}', f'Anchor: {anchor}', f'Equivalence: {result}',
              f'Final molecules: inherited={na}, isolated={nb}',
              f'Final molecule mismatches: {bad}', f'Coordinate RMSD (A): {rmsd}',
              f'Growth trees: inherited={len(ta)}, isolated={len(tb)}',
              f'Growth-tree content mismatches: {differing_trees}',
              'Score observations: ' + ', '.join(f'{k}={len(v)}' for k, v in scores.items()),
              'Differences:'] + (errors or ['None']) + ['Examples:'] + (samples or ['None']) + ['Notes:'] + (notes or ['None'])
    (out / 'equivalence_report.txt').write_text('\n'.join(report) + '\n')
    return [system, anchor, result, na, nb, bad, len(ta), len(tb), differing_trees, rmsd, '; '.join(errors)]


def main():
    tasks = RUNS / 'tasks.tsv'
    if not tasks.is_file():
        sys.exit('ERROR: Missing runs/tasks.tsv')
    RESULTS.mkdir(parents=True, exist_ok=True)
    rows = []
    for line in tasks.read_text().splitlines():
        if not line.strip():
            continue
        system, anchor = line.split('\t')
        try:
            row = compare_pair(system, anchor)
        except Exception as exc:
            folder = RESULTS / system / f'anc_{anchor}'
            folder.mkdir(parents=True, exist_ok=True)
            (folder / 'equivalence_report.txt').write_text(f'Comparison error: {exc}\n')
            row = [system, anchor, 'ERROR', '', '', '', '', '', '', '', str(exc)]
        rows.append(row)
        print(system, f'anc_{anchor}', row[2], flush=True)
    with (RESULTS / 'equivalence_summary.csv').open('w', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(['system', 'anchor', 'equivalence', 'inherited_final', 'isolated_final',
                         'final_molecule_mismatches', 'inherited_growth_trees', 'isolated_growth_trees',
                         'growth_tree_mismatches', 'coordinate_rmsd_A', 'details'])
        writer.writerows(rows)
    print('Summary:', RESULTS / 'equivalence_summary.csv')
    return int(not rows or any(row[2] != 'PASS' for row in rows))


if __name__ == '__main__':
    sys.exit(main())
