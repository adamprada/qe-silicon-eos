"""Birch-Murnaghan EOS fit and final report for bulk Si."""

import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

RESULTS_DIR = Path("results")
EOS_CSV = RESULTS_DIR / "eos_data.csv"

# 1 eV/Å³ = 160.21766208 GPa
EV_ANG3_TO_GPA = 160.21766208

# Expected PBE ranges for Si
A0_MIN, A0_MAX = 5.43, 5.48   # Å
B0_MIN, B0_MAX = 85.0, 95.0   # GPa


# ---------------------------------------------------------------------------
# Birch-Murnaghan EOS (3rd order)
# ---------------------------------------------------------------------------

def birch_murnaghan(V, E0, V0, B0, B0p):
    """E(V) — B0 in eV/Å³, V in Å³, E in eV."""
    x = (V0 / V) ** (2.0 / 3.0)
    return E0 + (9.0 * V0 * B0 / 16.0) * (
        (x - 1.0) ** 3 * B0p + (x - 1.0) ** 2 * (6.0 - 4.0 * x)
    )


# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------

def load_data() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    a_vals, vols, energies = [], [], []
    with open(EOS_CSV) as f:
        reader = csv.DictReader(f)
        for row in reader:
            a_vals.append(float(row["lattice_param_AA"]))
            vols.append(float(row["volume_per_atom_AA3"]))
            energies.append(float(row["energy_per_atom_eV"]))
    return np.array(a_vals), np.array(vols), np.array(energies)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    a_vals, vols, energies = load_data()

    # Initial guesses
    i_min = np.argmin(energies)
    p0 = [energies[i_min], vols[i_min], 0.56, 4.0]   # E0, V0, B0(eV/Å³), B0p

    popt, pcov = curve_fit(birch_murnaghan, vols, energies, p0=p0, maxfev=10000)
    E0, V0, B0_raw, B0p = popt
    perr = np.sqrt(np.diag(pcov))

    B0_gpa = B0_raw * EV_ANG3_TO_GPA

    # Equilibrium cubic lattice parameter: V_per_atom = a³/8
    a0 = (8.0 * V0) ** (1.0 / 3.0)

    # Fit residuals (RMS, meV/atom)
    e_fit = birch_murnaghan(vols, *popt)
    residuals_meV = (energies - e_fit) * 1000.0
    rms_meV = np.sqrt(np.mean(residuals_meV ** 2))

    # -----------------------------------------------------------------------
    # Report
    # -----------------------------------------------------------------------
    print("=" * 60)
    print("Birch-Murnaghan EOS fit — bulk Si (PBE)")
    print("=" * 60)
    print(f"  a0   = {a0:.4f} Å")
    print(f"  V0   = {V0:.4f} Å³/atom")
    print(f"  E0   = {E0:.6f} eV/atom")
    print(f"  B0   = {B0_gpa:.2f} GPa")
    print(f"  B0'  = {B0p:.3f}")
    print(f"  RMS residual = {rms_meV:.4f} meV/atom")
    print()

    a0_ok = A0_MIN <= a0 <= A0_MAX
    b0_ok = B0_MIN <= B0_gpa <= B0_MAX
    print(f"  Expected PBE a0: {A0_MIN}–{A0_MAX} Å  →  "
          f"{'PASS' if a0_ok else 'OUTSIDE EXPECTED RANGE'}")
    print(f"  Expected PBE B0: {B0_MIN}–{B0_MAX} GPa  →  "
          f"{'PASS' if b0_ok else 'OUTSIDE EXPECTED RANGE'}")
    print("=" * 60)

    # -----------------------------------------------------------------------
    # Plot
    # -----------------------------------------------------------------------
    v_dense = np.linspace(vols.min() * 0.995, vols.max() * 1.005, 300)
    e_dense = birch_murnaghan(v_dense, *popt)

    fig, ax = plt.subplots()
    ax.plot(vols, energies, "o", label="DFT data", zorder=3)
    ax.plot(v_dense, e_dense, "-", label="BM fit")
    ax.axvline(V0, color="gray", linestyle=":", linewidth=0.8, label=f"V0 = {V0:.3f} Å³")
    ax.set_xlabel("Volume per atom (Å³)")
    ax.set_ylabel("Energy per atom (eV)")
    ax.set_title(f"BM EOS fit — bulk Si\na0 = {a0:.4f} Å,  B0 = {B0_gpa:.1f} GPa")
    ax.legend()
    out_png = RESULTS_DIR / "eos_fit.png"
    fig.savefig(out_png, dpi=150)
    plt.close(fig)
    print(f"\nSaved: {out_png}")

    # -----------------------------------------------------------------------
    # JSON summary
    # -----------------------------------------------------------------------
    summary = {
        "a0_AA": round(a0, 6),
        "V0_AA3_per_atom": round(V0, 6),
        "E0_eV_per_atom": round(E0, 8),
        "B0_GPa": round(B0_gpa, 4),
        "B0_prime": round(B0p, 6),
        "rms_residual_meV_per_atom": round(rms_meV, 6),
        "a0_in_PBE_range": bool(a0_ok),
        "B0_in_PBE_range": bool(b0_ok),
    }
    out_json = RESULTS_DIR / "eos_summary.json"
    out_json.write_text(json.dumps(summary, indent=2))
    print(f"Saved: {out_json}")


if __name__ == "__main__":
    main()
