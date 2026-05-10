"""Volume sweep production run for bulk Si EOS."""

import csv
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from ase.build import bulk
from ase.calculators.espresso import Espresso, EspressoProfile

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
LATTICE_PARAM_EXP = 5.431           # Å, experimental Si
SCALE_FACTORS = [0.97, 0.98, 0.99, 1.00, 1.01, 1.02, 1.03]

PSEUDO_DIR = Path("pseudo")
PSEUDO_FILE = "Si.pbe-n-rrkjus_psl.1.0.0.UPF"

PARAMS_FILE = Path("results/converged_params.json")
QE_RUNS_DIR = Path("qe_runs")
RESULTS_DIR = Path("results")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def load_params() -> dict:
    if not PARAMS_FILE.exists():
        print(f"ERROR: {PARAMS_FILE} not found. Run convergence.py first.")
        sys.exit(1)
    params = json.loads(PARAMS_FILE.read_text())
    print(f"Loaded converged params: ecutwfc={params['ecutwfc']} Ry, "
          f"kpts={params['kpts']}")
    return params


def make_profile() -> EspressoProfile:
    return EspressoProfile(
        command="pw.x",
        pseudo_dir=str(PSEUDO_DIR.resolve()),
    )


def run_scf(a: float, params: dict) -> tuple[float, float, float]:
    """Build Si cell at lattice param a, run SCF, return (a, V/atom, E/atom)."""
    atoms = bulk("Si", "diamond", a=a)
    n_atoms = len(atoms)
    vol_per_atom = atoms.get_volume() / n_atoms

    label = f"{a:.4f}"
    run_dir = QE_RUNS_DIR / f"vol_{label}_AA"
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
                "ecutwfc": params["ecutwfc"],
                "ecutrho": params["ecutrho"],
            },
            "electrons": {
                "conv_thr": 1.0e-9,
            },
        },
        kpts=tuple(params["kpts"]),
        directory=str(run_dir),
    )

    atoms.calc = calc

    try:
        energy = atoms.get_potential_energy()
    except Exception as exc:
        print(f"\nERROR: QE run failed in {run_dir}")
        print(f"  {exc}")
        print("Inspect the QE output files in that directory for details.")
        sys.exit(1)

    e_per_atom = energy / n_atoms
    return a, vol_per_atom, e_per_atom


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    QE_RUNS_DIR.mkdir(exist_ok=True)
    RESULTS_DIR.mkdir(exist_ok=True)

    params = load_params()

    print()
    print("=" * 60)
    print("Volume sweep — bulk Si")
    print("=" * 60)

    rows = []
    for s in SCALE_FACTORS:
        a = LATTICE_PARAM_EXP * s
        print(f"  a = {a:.4f} Å (scale {s:.2f}) ...", end="", flush=True)
        a_val, vol, energy = run_scf(a, params)
        rows.append((a_val, vol, energy))
        print(f"  V = {vol:.4f} Å³/atom  E = {energy:.6f} eV/atom")

    # Save CSV
    csv_path = RESULTS_DIR / "eos_data.csv"
    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["lattice_param_AA", "volume_per_atom_AA3", "energy_per_atom_eV"])
        writer.writerows(rows)
    print(f"\nSaved: {csv_path}")

    # Plot raw E vs V
    vols = [r[1] for r in rows]
    energies = [r[2] for r in rows]

    fig, ax = plt.subplots()
    ax.plot(vols, energies, "o-")
    ax.set_xlabel("Volume per atom (Å³)")
    ax.set_ylabel("Energy per atom (eV)")
    ax.set_title("E vs V — bulk Si (raw DFT data)")
    out = RESULTS_DIR / "eos_raw.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"Saved: {out}")


if __name__ == "__main__":
    main()
