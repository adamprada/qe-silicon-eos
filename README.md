# qe-silicon-eos

Compute the equilibrium lattice parameter and bulk modulus of bulk silicon from first-principles density functional theory by running a Quantum ESPRESSO volume sweep and fitting a Birch-Murnaghan equation of state.

## Overview

The pipeline drives Quantum ESPRESSO from Python through ASE. A two-atom diamond-cubic primitive cell of Si is generated, total energies are computed at a small set of unit-cell volumes spanning the equilibrium, and a Birch-Murnaghan equation of state is fit to recover the equilibrium lattice parameter, the equilibrium energy, and the bulk modulus.

A separate convergence script verifies that the chosen plane-wave cutoff and Brillouin-zone sampling are sufficient before the production sweep is run.

## Design choices and assumptions

The system is bulk silicon in the diamond-cubic structure, treated in its 2-atom primitive cell. Silicon is chosen because it has a small primitive cell, a well-known equation of state from both experiment and DFT references, and a closed-shell electronic structure that does not require spin polarization or Hubbard corrections.

The DFT package is Quantum ESPRESSO, an open-source plane-wave pseudopotential code. The plane-wave formulation is appropriate for periodic solids and the Apple Silicon conda-forge build is robust. The exchange-correlation functional is PBE, with a pseudopotential from the Standard Solid State Pseudopotentials (SSSP) efficiency library. PBE is known to systematically overestimate the lattice parameter of Si by a fraction of a percent relative to experiment, and to underestimate the bulk modulus by several percent. These are accepted reference deviations of the functional, not bugs in the workflow.

The plane-wave cutoff and k-point mesh are selected by an explicit convergence test rather than picked from defaults. The production sweep covers seven volumes ranging from roughly 0.93 to 1.07 of the experimental equilibrium volume. The equation of state is fit with the third-order Birch-Murnaghan form, which is appropriate over this volume range.

No phonon, finite-temperature, or zero-point contributions are included. Reported lattice parameter and bulk modulus are static, electronic-energy values at 0 K, suitable for comparison to other static DFT references but not directly to room-temperature experiment.

## Requirements

* Python 3.10 or later
* Quantum ESPRESSO 7.x (`qe` from conda-forge recommended)
* ASE
* NumPy, SciPy, Matplotlib

Conda is the most reliable way to install Quantum ESPRESSO on macOS, including Apple Silicon. Alternatively, QE can be compiled from source as described below.

## Compiling Quantum ESPRESSO on Apple Silicon (M1/M2/M3)

Install dependencies via Homebrew:

```
brew install gcc open-mpi cmake fftw
```

Download the desired QE release from https://www.quantum-espresso.org and unzip it. Then build with CMake:

```
cd quantum-espresso-x.y.z
mkdir build
cd build
cmake -DCMAKE_C_COMPILER=mpicc -DCMAKE_Fortran_COMPILER=mpif90 -DCMAKE_INSTALL_PREFIX=/path/to/install ..
make -jN all
```

where `N` is the number of CPU cores to use for parallel compilation (e.g. `make -j8 all`). Optionally run `make install` to copy the binaries to the prefix directory. Add the `bin/` subdirectory of your install prefix (or `build/bin/` if no prefix was set) to your `PATH` so that `pw.x` is accessible.

## Installation

```
git clone <repo-url> qe-silicon-eos
cd qe-silicon-eos
conda env create -f environment.yml
conda activate qe-silicon-eos
```

For Python dependencies only (if `pw.x` is already on PATH):

```
pip install -r requirements.txt
```

Pseudopotentials are downloaded by `setup_pseudos.py` on first run.

## Usage

```
python setup_pseudos.py     # download Si pseudopotential into pseudo/
python convergence.py       # converge cutoff and k-mesh
python volume_sweep.py      # production volume sweep
python fit_eos.py           # Birch-Murnaghan fit and report
```

Quantum ESPRESSO output files for each calculation are written to `qe_runs/`. Final summary, energy-volume curve, and EOS fit are written to `results/`.

## Project structure

```
qe-silicon-eos/
├── README.md
├── report.md               # full technical report with figures
├── environment.yml         # conda environment spec
├── requirements.txt        # pip-installable Python dependencies
├── setup_pseudos.py        # download SSSP Si pseudopotential
├── convergence.py          # cutoff and k-mesh convergence
├── volume_sweep.py         # E(V) production runs
├── fit_eos.py              # Birch-Murnaghan fit, plot, report
├── pseudo/                 # SSSP pseudopotential files (gitignored)
├── qe_runs/                # QE input/output (gitignored)
└── results/                # output figures, CSV data, JSON summaries
    ├── conv_ecutwfc_energy.png
    ├── conv_ecutwfc.png
    ├── conv_kmesh_energy.png
    ├── conv_kmesh.png
    ├── converged_params.json
    ├── eos_data.csv
    ├── eos_raw.png
    ├── eos_fit.png
    └── eos_summary.json
```

## Verification

The fitted equilibrium lattice parameter for Si in PBE should be in the range of about 5.43 to 5.48 Å, and the fitted bulk modulus in the range of about 85 to 95 GPa. These are standard PBE reference values reproduced widely in the literature. Falling outside these ranges indicates a setup problem (insufficient cutoff, too coarse a k-mesh, mismatched pseudopotential) rather than a successful calculation.

This workflow produced a₀ = 5.4696 Å and B₀ = 88.71 GPa (ecutwfc = 40 Ry, 10×10×10 k-mesh, RMS fit residual 0.008 meV/atom), both within the expected PBE ranges.

## References

* P. Giannozzi et al., "Quantum ESPRESSO: a modular and open-source software project for quantum simulations of materials" (2009)
* F. Birch, "Finite Elastic Strain of Cubic Crystals" (1947)
* Standard Solid State Pseudopotentials: https://www.materialscloud.org/discover/sssp/
* ASE Quantum ESPRESSO calculator: https://wiki.fysik.dtu.dk/ase/ase/calculators/espresso.html

## License

MIT
