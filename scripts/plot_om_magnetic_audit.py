#!/usr/bin/env python3
"""Summarize only completed OM magnetic-audit runs; no recalculation."""
from pathlib import Path
import csv
import os
root = Path(__file__).resolve().parents[1]
folder = root / 'data/tbg/fortran/om/magnetic_audit'
os.environ.setdefault('MPLCONFIGDIR', '/private/tmp/tbg-om-matplotlib')
os.environ.setdefault('XDG_CACHE_HOME', '/private/tmp/tbg-om-cache')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

def read(name):
    lines = (folder / name).read_text().splitlines()
    if not any(x.startswith('wall_seconds') for x in lines):
        raise ValueError(f'incomplete run: {name}')
    out = {'formula': {}, 'name': name}
    for line in lines:
        fields = line.split()
        if not fields or fields[0] == '#':
            continue
        if fields[0] == 'formula':
            out['formula'][int(fields[1])] = np.array(fields[2:], float)
        else:
            out[fields[0]] = np.array(fields[1:], float)
    return out

runs = [read(p.name) for p in sorted(folder.glob('tr*.txt'))]
with (folder / 'summary.csv').open('w') as f:
    w = csv.writer(f)
    w.writerow(['run', 'method', 'ncut', 'Mx_muB_cell', 'My_muB_cell'])
    for r in runs:
        for cut, val in r['formula'].items():
            w.writerow([r['name'], 'formula', cut, *val])
        for method in ['Omega_FD', 'HF_raw', 'HF_active_only', 'BZ_boundary', 'HF_minus_boundary', 'closure_error']:
            if method in r:
                w.writerow([r['name'], method, '', *r[method]])

plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False})
fig, ax = plt.subplots(2, 2, figsize=(11, 7.6), constrained_layout=True)
basis = [read('tr3_nk12_B0p01_phi0.txt'), read('tr4_nk12_B0p001_phi0.txt'), read('tr5_nk12_B0p01_phi0.txt')]
ax[0,0].plot([3,4,5], [r['HF_raw'][0] for r in basis], 'o-', label='Raw grand-potential derivative (HF)')
ax[0,0].plot([3,4,5], [r['BZ_boundary'][0] for r in basis], 'x--', label='Independent BZ boundary')
ax[0,0].set(xlabel='Plane-wave cutoff tr', ylabel=r'$M_x$ ($\mu_B$/cell)', xticks=[3,4,5], title='Large background grows with basis cutoff')
ax[0,0].legend(fontsize=8)

meshes = [read('tr4_nk8_B0p01_phi0_accelerate.txt'), read('tr4_nk12_B0p001_phi0.txt'), read('tr4_nk40_B0p01_phi0.txt')]
ax[0,1].semilogy([8,12,40], [max(abs(r['closure_error'])) for r in meshes], 'o-')
ax[0,1].set(xlabel='k-mesh nk', ylabel=r'Max component error ($\mu_B$/cell)', title='HF minus boundary vs full-band formula', xticks=[8,12,40])

main = meshes[-1]
cuts = np.array([c for c in sorted(main['formula']) if c >= 20])
for i, name in enumerate(['x', 'y']):
    values = np.array([main['formula'][c][i] for c in cuts])
    ax[1,i].plot(cuts, values, 'o-', color=['#2166ac','#b35806'][i])
    ax[1,i].axvline(100, ls='--', lw=1, color='gray', label='Production ncut=100')
    ax[1,i].set(xlabel='Band-window cutoff ncut', ylabel=rf'$M_{name}$ ($\mu_B$/cell)', title=f'Formula: component {name}, tr=4, nk=40')
    ax[1,i].legend(fontsize=8)
for a in ax.flat:
    a.grid(alpha=0.18)
fig.suptitle(r'TBG: $\theta=1.4^\circ$, $\epsilon=0.1\%$, $\phi=0^\circ$, $\mu=0.5$ meV' + '\nSingle spin and valley; raw continuum grand potential is unregularized', fontsize=12)
fig.savefig(folder / 'comparison.png', dpi=180)
fig.savefig(folder / 'comparison.pdf')
print(folder / 'summary.csv')
print(folder / 'comparison.png')
