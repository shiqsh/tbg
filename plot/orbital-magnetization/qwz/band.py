#!/usr/bin/env python3  

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[3]

# ============================================================
# Parameters
# ============================================================

numk = 100

m = 1.0
lambda_ = 0.2
t_perp = 0.2

# ============================================================
# File name
# ============================================================

name = f"qwz_band_lambda_{lambda_:.3f}_m_{m:.3f}_tperp_{t_perp:.3f}"

data_file = (
    ROOT
    / "data"
    / "orbital-magnetization"
    / "qwz"
    / "band"
    / f"{name}.dat"
)

figure_file = (
    ROOT
    / "figure"
    / "orbital-magnetization"
    / "qwz"
    / "band"
    / f"{name}.pdf"
)

# ============================================================
# Read data
# ============================================================

data = np.loadtxt(data_file)

dk = data[:, 0]
k = data[:, 1:3]
energy = data[:, 3:7]

# ============================================================
# Plot
# ============================================================
fig, ax = plt.subplots(figsize=(8, 8))

for iband in range(4):
    ax.plot(dk, energy[:, iband])

ticks = dk[[0, numk, 2*numk, 3*numk]]

ax.set_xticks(
    ticks,
    [r"$\Gamma$", "X", "M", r"$\Gamma$"]
)

for x in ticks:
    ax.axvline(x, linewidth=0.8)

ax.axhline(0.0, color="k", linewidth=0.8)
# plt.axhline(-1.0, color="k", linewidth=0.8)

ax.set_ylabel("Energy")
ax.set_xlim(dk[0], dk[-1])
ax.set_title(
    rf"Bilayer QWZ band structure "
    rf"($m={m:.1f}$, $t_\perp={t_perp:.1f}$, $\lambda={lambda_:.3f}$)"
)

ax.legend()
fig.tight_layout()

figure_file.parent.mkdir(parents=True, exist_ok=True)

fig.savefig(figure_file, dpi=300)
plt.show()