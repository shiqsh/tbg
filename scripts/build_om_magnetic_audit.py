#!/usr/bin/env python3
"""Build a source snapshot with its own tr; never modify production constants.

Requires an existing fpm build of stdlib, gfortran, and macOS Accelerate.
Homebrew OpenBLAS is an optional alternative (--blas openblas). Example:
    python3 scripts/build_om_magnetic_audit.py --tr 4
"""
import argparse
from pathlib import Path
import re
import subprocess
import hashlib
import json

root = Path(__file__).resolve().parents[1]
p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--tr', type=int, default=4)
p.add_argument('--blas', choices=['accelerate','openblas'], default='accelerate')
args = p.parse_args()
if args.tr < 1:
    p.error('tr must be at least 1')
build = root / 'fortran/build'
dest = build / f'om_magnetic_audit_tr{args.tr}_{args.blas}'
dest.mkdir(parents=True, exist_ok=True)
mods = sorted(build.glob('*/stdlib_linalg.mod'), key=lambda x: x.stat().st_mtime, reverse=True)
libs = sorted(build.glob('*/tbg/libtbg.a'), key=lambda x: x.stat().st_mtime, reverse=True)
if not mods or not libs:
    raise SystemExit('First build the project with fpm to provide stdlib.')
sources = ['src/constants.f90', 'src/model_interface.f90', 'src/tbg/tbg.f90',
           'src/heterostrain/strain.f90', 'src/linalg.f90',
           'src/orbital-magnetization/om.f90', 'src/orbital-magnetization/om_cutoff.f90',
           'src/orbital-magnetization/om_magnetic.f90',
           'app/orbital-magnetization/tbg_strain/compare_om_magnetic.f90']
hashes = {}
for rel in sources:
    source = root / 'fortran' / rel
    content = source.read_text()
    if source.name == 'constants.f90':
        content, count = re.subn(r'(integer, parameter :: tr\s*=\s*)\d+',
                                lambda m: m[1] + str(args.tr), content)
        assert count == 1
    (dest / source.name).write_text(content)
    hashes[rel] = hashlib.sha256(content.encode()).hexdigest()
(dest / 'source_sha256.json').write_text(json.dumps(hashes, indent=2) + '\n')
cmd = ['gfortran', '-O3', '-fopenmp', '-ffree-line-length-none',
       '-J.', '-I' + str(mods[0].parent), *[Path(s).name for s in sources], str(libs[0]),
       *(['-framework', 'Accelerate'] if args.blas == 'accelerate' else
         ['-L/opt/homebrew/opt/openblas/lib', '-lopenblas']), '-o', 'compare_om_magnetic']
subprocess.run(cmd, cwd=dest, check=True)
print(dest / 'compare_om_magnetic')
