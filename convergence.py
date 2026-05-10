"""Convergence test for ecutwfc and k-point mesh for bulk Si with Quantum ESPRESSO."""

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from ase.build import bulk
from ase.calculators.espresso import Espresso, EspressoProfile

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
LATTICE_PARAM = 5.431           # Å, experimental Si
PSEUDO_DIR = Path("pseudo")
PSEUDO_FILE = "Si.pbe-n-rrkjus_psl.1.0.0.UPF"
DUAL = 8                        # ecutrho = DUAL * ecutwfc (ultrasoft pseudo)
CONV_THR = 1.0e-9               # Ry, SCF convergence threshold

KPTS_FIXED = (6, 6, 6)
ECUTWFC_VALUES = [30, 40, 50, 60, 70]   # Ry

KPTS_VALUES = [(4, 4, 4), (6, 6, 6), (8, 8, 8), (10, 10, 10)]

CONVERGENCE_THR_EV = 1.0e-3    # eV/atom (1 meV/atom)

QE_RUNS_DIR = Path("qe_runs")
RESULTS_DIR = Path("results")

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def make_profile() -> EspressoProfile:
    return EspressoProfile(
        command="pw.x",
        pseudo_dir=str(PSEUDO_DIR.resolve()),
    )


def make_si_cell() -> object:
    return bulk("Si", "diamond", a=LATTICE_PARAM)


def run_scf(atoms, ecutwfc: int, kpts: tuple, run_dir: Path) -> float:
    """Run a single SCF calculation; return energy per atom in eV. Exits on failure."""
    run_dir.mkdir(parents=True, exist_ok=True)

    calc = Espresso(
        profile=make_profile(),
        pseudopotentials={"Si": PSEUDO_FILE},
        input_data={
            "control": {
                "calculation": "scf",
                "restart_mode": "from_scratch",
                "prefix": "si",
                "verbosity": "low",
            },
            "system": {
                "ecutwfc": ecutwfc,
                "ecutrho": DUAL * ecutwfc,
            },
            "electrons": {
                "conv_thr": CONV_THR,
            },
        },
        kpts=kpts,
        directory=str(run_dir),
    )

    atoms_copy = atoms.copy()
    atoms_copy.calc = calc

    try:
        energy = atoms_copy.get_potential_energy()
    except Exception as exc:
        print(f"\nERROR: QE run failed in {run_dir}")
        print(f"  {exc}")
        print("Inspect the QE output files in that directory for details.")
        sys.exit(1)

    n_atoms = len(atoms_copy)
    return energy / n_atoms


# ---------------------------------------------------------------------------
# Ecutwfc sweep
# ---------------------------------------------------------------------------

def sweep_ecutwfc() -> tuple[list[int], list[float]]:
    print("=" * 60)
    print("Sweep 1: ecutwfc (k-mesh fixed at 6×6×6)")
    print("=" * 60)

    si = make_si_cell()
    energies = []

    for ecut in ECUTWFC_VALUES:
        run_dir = QE_RUNS_DIR / f"conv_ecutwfc_{ecut}"
        print(f"  ecutwfc = {ecut:3d} Ry  ->  {run_dir}", end="", flush=True)
        e_per_atom = run_scf(si, ecut, KPTS_FIXED, run_dir)
        energies.append(e_per_atom)
        print(f"  E = {e_per_atom:.6f} eV/atom")

    return ECUTWFC_VALUES, energies


def find_converged_ecutwfc(ecuts: list[int], energies: list[float]) -> int:
    ref = energies[-1]
    for i, (ecut, e) in enumerate(zip(ecuts, energies)):
        if abs(e - ref) < CONVERGENCE_THR_EV:
            return ecut
    return ecuts[-1]


# ---------------------------------------------------------------------------
# K-mesh sweep
# ---------------------------------------------------------------------------

def sweep_kmesh(ecutwfc_fixed: int) -> tuple[list[int], list[float]]:
    print()
    print("=" * 60)
    print(f"Sweep 2: k-mesh (ecutwfc fixed at {ecutwfc_fixed} Ry)")
    print("=" * 60)

    si = make_si_cell()
    k_sizes = []
    energies = []

    for kpts in KPTS_VALUES:
        n = kpts[0]
        run_dir = QE_RUNS_DIR / f"conv_kmesh_{n}x{n}x{n}"
        print(f"  k-mesh = {n}×{n}×{n}  ->  {run_dir}", end="", flush=True)
        e_per_atom = run_scf(si, ecutwfc_fixed, kpts, run_dir)
        k_sizes.append(n)
        energies.append(e_per_atom)
        print(f"  E = {e_per_atom:.6f} eV/atom")

    return k_sizes, energies


def find_converged_kmesh(k_sizes: list[int], energies: list[float]) -> int:
    ref = energies[-1]
    for n, e in zip(k_sizes, energies):
        if abs(e - ref) < CONVERGENCE_THR_EV:
            return n
    return k_sizes[-1]


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_ecutwfc(ecuts: list[int], energies: list[float]) -> None:
    RESULTS_DIR.mkdir(exist_ok=True)
    ref = energies[-1]
    delta = [abs(e - ref) * 1000 for e in energies]   # meV/atom

    fig, ax = plt.subplots()
    ax.plot(ecuts, delta, "o-")
    ax.axhline(1.0, color="red", linestyle="--", label="1 meV/atom threshold")
    ax.set_xlabel("ecutwfc (Ry)")
    ax.set_ylabel("|ΔE| (meV/atom) relative to 70 Ry")
    ax.set_title("Ecutwfc convergence — bulk Si")
    ax.legend()
    ax.set_yscale("log")
    out = RESULTS_DIR / "conv_ecutwfc.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"Saved: {out}")


def plot_kmesh(k_sizes: list[int], energies: list[float]) -> None:
    RESULTS_DIR.mkdir(exist_ok=True)
    ref = energies[-1]
    delta = [abs(e - ref) * 1000 for e in energies]   # meV/atom

    fig, ax = plt.subplots()
    ax.plot(k_sizes, delta, "s-")
    ax.axhline(1.0, color="red", linestyle="--", label="1 meV/atom threshold")
    ax.set_xlabel("k-grid N (N×N×N Monkhorst-Pack)")
    ax.set_ylabel("|ΔE| (meV/atom) relative to 10×10×10")
    ax.set_title("K-mesh convergence — bulk Si")
    ax.legend()
    ax.set_yscale("log")
    out = RESULTS_DIR / "conv_kmesh.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"Saved: {out}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    QE_RUNS_DIR.mkdir(exist_ok=True)
    RESULTS_DIR.mkdir(exist_ok=True)

    ecuts, e_ecut = sweep_ecutwfc()
    plot_ecutwfc(ecuts, e_ecut)

    rec_ecut = find_converged_ecutwfc(ecuts, e_ecut)
    print(f"\nRecommended ecutwfc: {rec_ecut} Ry")

    k_sizes, e_kmesh = sweep_kmesh(rec_ecut)
    plot_kmesh(k_sizes, e_kmesh)

    rec_k = find_converged_kmesh(k_sizes, e_kmesh)
    print(f"Recommended k-mesh : {rec_k}×{rec_k}×{rec_k}")

    params = {
        "ecutwfc": rec_ecut,
        "ecutrho": DUAL * rec_ecut,
        "kpts": [rec_k, rec_k, rec_k],
    }
    out = RESULTS_DIR / "converged_params.json"
    out.write_text(json.dumps(params, indent=2))
    print(f"\nSaved converged params to: {out}")
    print(json.dumps(params, indent=2))


if __name__ == "__main__":
    main()
