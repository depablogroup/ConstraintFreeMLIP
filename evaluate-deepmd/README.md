# Evaluate Free Energy Surface with DeePMD

This directory stores the code and saved models we used to perform free energy
surface estimation with **DeePMD**.

## Scripts

This scripts under this folder uses **PySAGES**, **ASE**, and **DeePMD** to learn the FES of a ADP molecule using Spectral-ABF sampling.

- **pysages-deepmd.py**: This is the main script that performs MD simulations with Spectral-ABF and analyzes the FES. It has three arguments:
  - `--savedir`: Directory where results and models are saved.
  - `--expname`: Experiment name for files.
  - `--restart`: Flag to restart the simulation from the last saved state.
- **run.sh**: A bash script for running the simulation
- **pysages-deepmd-water**: A python script for running the simulation with water in the simulation box

## Water + adp

ADP with flexible water data.

- **unbiased**
  - an MLIP model trained with 10,000 snapshots generated from unbiased simulations at 600 K.

- **uniform**
  - an MLIP model trained with 10,000 snapshots generated from uniform simulations at 600 K.

## Explicit

ADP with explicit solvent, using force-only data.

- **unbiased40k**
  - an MLIP model trained with 40,000 snapshots generated from unbiased simulations at 600 K.

- **uniform40k**
  - an MLIP model trained with 40,000 snapshots generated from uniform simulations at 600 K.

## Vacuum

ADP in vacuum 

  - an MLIP model trained with 1,000 snapshots generated from unbiased simulations at 600 K.





