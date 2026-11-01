#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/000.config.sh"

# Generate one identical input pair for each system and anchor.
python3 - "$SYSTEM_LIST" "$SYSTEM_ROOT" "$ANCHOR_ROOT" "$N_anc" \
    "$TEMPLATE_INPUT" "$RUN_ROOT" "$INHERITED_DOCK_BIN" "$ISOLATED_DOCK_BIN" \
    "$SCRIPT_DIR" <<'PY'
from pathlib import Path
import hashlib
import re
import sys

list_file, system_root, anchor_root, n_str, template_file, run_root, inherited_bin, isolated_bin, work = sys.argv[1:]
list_file, system_root = Path(list_file), Path(system_root)
anchor_root, template_file = Path(anchor_root), Path(template_file)
run_root, work = Path(run_root), Path(work)
n = int(n_str)
if n < 1:
    raise SystemExit('ERROR: N_anc must be positive')
if run_root.exists():
    raise SystemExit(f'ERROR: Runs already exist: {run_root}')
for name, path in [('system list', list_file), ('validated template', template_file),
                   ('inherited binary', Path(inherited_bin)), ('isolated binary', Path(isolated_bin))]:
    if not path.is_file():
        raise SystemExit(f'ERROR: Missing {name}: {path}')
if not Path(inherited_bin).stat().st_mode & 0o111 or not Path(isolated_bin).stat().st_mode & 0o111:
    raise SystemExit('ERROR: Both binaries must be executable')

systems = []
for line in list_file.read_text().splitlines():
    name = line.split('#', 1)[0].strip()
    if not name:
        continue
    if not re.fullmatch(r'[A-Za-z0-9_-]+', name) or name in systems:
        raise SystemExit(f'ERROR: Invalid or duplicate system: {name}')
    systems.append(name)
if not systems:
    raise SystemExit('ERROR: No systems listed')

original = template_file.read_text().splitlines()
required = {
    'conformer_search_type': 'denovo',
    'dn_mode': 'sample',
    'dn_sampling_method': 'graph',
    'dn_user_specified_anchor': 'yes',
    'descriptor_score_primary': 'yes',
    'descriptor_use_hungarian': 'yes',
    'descriptor_use_dynamic_reference': 'yes',
    'dn_write_growth_trees': 'yes',
    'simplex_random_seed': '0',
}
keys = [line.split()[0] for line in original if line.split() and not line.lstrip().startswith('#')]
for key, value in required.items():
    matches = [line.split()[1] for line in original if line.split() and line.split()[0] == key and len(line.split()) > 1]
    if matches != [value]:
        raise SystemExit(f'ERROR: Template must contain exactly {key} {value}')
change_keys = ('receptor_site_file', 'descriptor_grid_score_grid_prefix',
               'descriptor_hms_score_ref_filename', 'dn_fraglib_anchor_file')
for key in change_keys + ('dn_output_prefix',):
    if keys.count(key) != 1:
        raise SystemExit(f'ERROR: Expected exactly one {key} in template')

# Validate the entire experiment before writing any run directories.
tasks = []
for system in systems:
    sysdir = system_root / system
    replacements = {
        'receptor_site_file': str(sysdir / f'{system}.rec.clust.close.sph'),
        'descriptor_grid_score_grid_prefix': str(sysdir / f'{system}.rec'),
        'descriptor_hms_score_ref_filename': str(sysdir / f'{system}.lig.am1bcc.mol2'),
    }
    for path in (Path(replacements['receptor_site_file']),
                 Path(replacements['descriptor_grid_score_grid_prefix'] + '.nrg'),
                 Path(replacements['descriptor_grid_score_grid_prefix'] + '.bmp'),
                 Path(replacements['descriptor_hms_score_ref_filename'])):
        if not path.is_file():
            raise SystemExit(f'ERROR: Missing system input: {path}')
    manifest = work / 'anchors' / f'{system}.txt'
    if manifest.exists():
        anchors = []
        for line in manifest.read_text().splitlines():
            entry = line.split('#', 1)[0].strip()
            if entry:
                path = Path(entry)
                anchors.append(path if path.is_absolute() else anchor_root / entry)
    else:
        # Supports the shared anchor_1.mol2, anchor_2.mol2, ... convention.
        anchors = [anchor_root / f'anchor_{i}.mol2' for i in range(1, n + 1)]
    if len(anchors) < n:
        raise SystemExit(f'ERROR: Only {len(anchors)} anchors specified for {system}')
    for i, anchor in enumerate(anchors[:n], 1):
        if not anchor.is_file():
            raise SystemExit(f'ERROR: Missing anchor: {anchor}')
        parameters = {**replacements, 'dn_fraglib_anchor_file': str(anchor)}
        lines = []
        for line in original:
            parts = line.split()
            if parts and parts[0] in parameters:
                line = f'{parts[0]:<65} {parameters[parts[0]]}'
            lines.append(line)
        generated = '\n'.join(lines) + '\n'
        # A copied input must not retain receptor/reference paths for another system.
        found_systems = set(re.findall(r'/002_DRef_Systems/([^/\s]+)/', generated))
        if found_systems and found_systems != {system}:
            raise SystemExit(f'ERROR: Stale system paths for {system}: {found_systems}')
        tasks.append((system, i, generated))

run_root.mkdir(parents=True)
for system, i, generated in tasks:
    for mode in ('inherited', 'isolated'):
        target = run_root / mode / system / f'anc_{i}'
        target.mkdir(parents=True)
        (target / 'dn_generic.in').write_text(generated)
    left = (run_root / 'inherited' / system / f'anc_{i}' / 'dn_generic.in').read_bytes()
    right = (run_root / 'isolated' / system / f'anc_{i}' / 'dn_generic.in').read_bytes()
    if left != right:
        raise SystemExit(f'ERROR: Nonidentical inputs for {system} anchor {i}')
(run_root / 'tasks.tsv').write_text(''.join(f'{system}\t{i}\n' for system, i, _ in tasks))
with (run_root / 'binary_sha256.txt').open('w') as handle:
    for mode, binary in (('inherited', inherited_bin), ('isolated', isolated_bin)):
        h = hashlib.sha256(Path(binary).read_bytes()).hexdigest()
        handle.write(f'{mode}\t{h}\t{binary}\n')
print(f'Prepared {len(tasks)} matched pairs across {len(systems)} system(s).')
print(f'Task list: {run_root / "tasks.tsv"}')
PY
