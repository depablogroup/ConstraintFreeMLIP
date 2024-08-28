# Evaluate Free Energy Surface with Allegro

This directory contains the code and saved models used to perform free energy surface (FES) estimation with **Allegro**.

### Scripts

The scripts in this folder utilize **PySAGES**, **ASE**, and **Allegro** to learn the FES of an ADP molecule using Spectral-ABF sampling.

- **simulationpysages.py**: This is the main script that performs MD simulations with Spectral-ABF and analyzes the FES. It has three arguments:
  - `--savedir`: Directory where results and models are saved.
  - `--expname`: Experiment name for files.
  - `--restart`: Flag to restart the simulation from the last saved state.
  - `--fixbond`: Flag to add bond constraints
- **run.sh**: A bash script for running the simulation.
- **simulation-pes**: A python script that analyzes the PES by computing histograms of dihedral angles and their corresponding potential energies.
- **simulationpysages_water**: A python script that performs MD simulations with Spectral-ABF with water in the simulation box.
- **convert_units.py**: A python script for unit conversion.

### Water + adp

Models trained on ADP with flexible water data. The directory is organized into two subfolders, each varying by sampling method:

- **unbiased**
  - an MLIP model trained with 10,000 snapshots generated from unbiased simulations at 600 K.
- **uniform**
  - an MLIP model trained with 10,000 snapshots generated from uniform (enhanced sampling) simulations at 600 K.

### Explicit

Models trained on ADP with explicit solvent, using force-only data. The directory is organized into six subfolders, each varying by dataset size, sampling methods, or temperatures:

- **nequip-forces_are_not_enough**
  - models trained with corrected units and uncorrected unit conversion reproduced from paper.
  - a sample nequip/allegro configuration file for reproducing the model.
- **unbiased-1k**
  - an MLIP model trained with 1,000 snapshots generated from unbiased simulations at 600 K.
- **unbiased-40k**
  - an MLIP model trained with 40,000 snapshots generated from unbiased simulations at 600 K.
- **unbiased-40k-300K**
  - an MLIP model trained with 40,000 snapshots generated from unbiased simulations at 300 K.
- **uniform-1k**
  - an MLIP model trained with 1,000 snapshots generated from uniform (enhanced sampling) simulations at 600 K.
- **uniform-40k**
  - an MLIP model trained with 40,000 snapshots generated from uniform (enhanced sampling) simulations at 600 K.

### Implicit

ADP in implicit solvent. The directory is organized into two subfolders, each varying by sampling methods:

- **unbiased**
  - an MLIP model trained with 1,000 snapshots generated from unbiased simulations at 600 K.
- **uniform**
  - an MLIP model trained with 1,000 snapshots generated from uniform (enhanced sampling) simulations at 600 K.

### Vacuum

ADP in vacuum. The directory is organized into four subfolders, each varying by dataset size, sampling methods, or temperatures:

- **unbiased**
  - an MLIP model trained with 1,000 snapshots generated from unbiased simulations at 600 K.
- **unbiased-dft**
  - an MLIP model trained with 1,000 snapshots generated from unbiased simulations at 600 K, calculated at the DFT level.
- **uniform**
  - an MLIP model trained with 1,000 snapshots generated from uniform (enhanced sampling) simulations at 600 K.
- **uniform-dft**
  - an MLIP model trained with 1,000 snapshots generated from uniform (enhanced sampling) simulations at 600 K, calculated at the DFT level.
