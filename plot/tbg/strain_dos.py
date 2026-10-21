#!/usr/bin/env python3
"""Plot saved heterostrain triangle DOS, without recalculation or smoothing."""

import argparse
from pathlib import Path

import numpy as np


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_NAME = "dos_theta1p4_epsilon0p1_phi90_nk100_ne2000_eminm50_emax50.dat"


def read_dos(path):
    """Read the two numeric columns and the parameter header saved by strain_dos."""
    metadata = {}
    with path.open(encoding="utf-8") as source:
        for line in source:
            if line.startswith("#") and "=" in line:
                key, value = line[1:].strip().split("=", 1)
                metadata[key.strip()] = value.strip()
    required = ("theta_deg", "epsilon_percent", "phi_deg", "nkx", "nky", "n_intervals",
                "E_min_meV", "E_max_meV", "method", "normalization")
    missing = [key for key in required if key not in metadata]
    if missing:
        raise ValueError(f"Missing metadata: {', '.join(missing)}; use strain_dos output.")
    data = np.loadtxt(path, comments="#", ndmin=2)
    if data.shape[1] != 2 or data.shape[0] < 2 or not np.all(np.isfinite(data)):
        raise ValueError("Expected two finite columns with at least two energy points.")
    if not np.all(np.diff(data[:, 0]) > 0):
        raise ValueError("Energy points must increase strictly.")
    if data.shape[0] != int(metadata["n_intervals"]) + 1:
        raise ValueError("The energy-point count does not match the file header.")
    if metadata["method"] != "bin_averaged_analytic_triangle":
        raise ValueError("Expected method=bin_averaged_analytic_triangle.")
    if metadata["normalization"] != "states_per_eV_per_moire_cell_single_spin_single_valley":
        raise ValueError("Unrecognized DOS normalization.")
    return data[:, 0]*1000.0, data[:, 1], metadata


def plot_dos(energy_mev, density, metadata):
    import matplotlib.pyplot as plt

    plt.rcParams.update({
        "font.family": "serif", "font.size": 13, "mathtext.fontset": "dejavuserif",
        "text.usetex": False, "axes.linewidth": 1.2,
    })
    theta = float(metadata["theta_deg"])
    epsilon = float(metadata["epsilon_percent"])
    phi = float(metadata["phi_deg"])
    nkx, nky = int(metadata["nkx"]), int(metadata["nky"])
    fig, ax = plt.subplots(figsize=(6.4, 6))
    ax.plot(energy_mev, density, color="#0072B2", linewidth=1.7)
    ax.axvline(0.0, color="#CC0000", linewidth=0.8, linestyle="--", zorder=0)
    ax.set_xlim(float(metadata["E_min_meV"]), float(metadata["E_max_meV"]))
    ax.set_ylim(bottom=0.0)
    ax.set_xlabel(r"$E$ (meV)", fontsize=17)
    ax.set_ylabel(r"$D(E)$ [states / (eV moiré cell)]", fontsize=16)
    ax.set_title(rf"$\theta={theta:g}^\circ,\ \epsilon={epsilon:g}\%,\ \phi={phi:g}^\circ$"
                 f"\nAnalytic triangle DOS, {nkx} × {nky} grid cells", fontsize=14, pad=12)
    ax.tick_params(direction="in", which="both", top=True, right=True, length=5)
    ax.ticklabel_format(axis="y", style="plain", useOffset=False)
    ax.grid(False)
    fig.text(0.5, 0.025, "Single spin, single valley · no broadening or smoothing",
             ha="center", fontsize=10)
    fig.tight_layout(rect=(0, 0.055, 1, 1))
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path,
                        default=PROJECT_ROOT / "data" / "tbg" / "fortran" / "dos" / DEFAULT_NAME)
    parser.add_argument("--output", type=Path, help="PNG or PDF path; saves both formats using this stem")
    parser.add_argument("--no-show", action="store_true", help="Save files without opening a window")
    args = parser.parse_args()
    if not args.input.is_file():
        parser.error(f"Missing DOS data: {args.input}; run fpm run --target strain_dos first.")
    if args.no_show:
        import matplotlib
        matplotlib.use("Agg")
    try:
        energy, density, metadata = read_dos(args.input)
    except (OSError, ValueError) as exc:
        parser.error(str(exc))
    fig = plot_dos(energy, density, metadata)
    output = args.output or (PROJECT_ROOT / "figure" / "tbg" / f"{args.input.stem}.png")
    output = output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    for suffix in (".png", ".pdf"):
        destination = output.with_suffix(suffix)
        fig.savefig(destination, dpi=300, bbox_inches="tight")
        print(f"Saved {destination}")
    if not args.no_show:
        import matplotlib.pyplot as plt
        plt.show()


if __name__ == "__main__":
    main()
