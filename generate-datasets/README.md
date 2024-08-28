# Training Data Generation

This directory contains the datasets and scripts used for data generation.

### Water + adp

This directory contains subfolders and a python script for generating the training data for ADP with flexible water.

- **unbiased**
  - restart and parameter files for unbiased simulations
- **uniform**
  - restart and parameter files for enhanced sampling simulations
  - **runcorr.py** corrects the biasing forces generated from the enhanced sampling simulation by performing energy/force calcultions using the **OpenMM** library on a given a trajectory.
- **run.py**
  - runs a MD simulation using the **OpenMM** library, combined with enhanced sampling methods provided by **PySAGES**. It supports both unbiased and biased simulations. It has two arguments:
  - `--unbiased`: run an unbiased simulation if is set `True`.
  - `--restart`: restart from a previously saved state.

### Explicit

This directory contains subfolders and a python script for generating the training data for ADP in explicit solvent.

- **unbiased**
  - restart and parameter files for unbiased simulations
- **uniform**
  - restart and parameter files for enhanced sampling simulations
  - **runcorr.py** corrects the biasing forces generated from the enhanced sampling simulation by performing energy/force calcultions using the **OpenMM** library on a given a trajectory.
- **run.py**
  - runs a MD simulation using the **OpenMM** library, combined with enhanced sampling methods provided by **PySAGES**. It supports both unbiased and biased simulations. It has two arguments:
  - `--unbiased`: run an unbiased simulation if is set `True`.
  - `--restart`: restart from a previously saved state.


### Implicit

This directory contains subfolders and two python scripts for generating the training data for ADP in implicit solvent.

- **unbiased**
  - restart and parameter files for unbiased simulations.
- **uniform**
  - restart and parameter files for enhanced sampling simulations
- **pysages-ala.py**
  - runs a MD simulation using the **ASE** library, combined with enhanced sampling methods provided by **PySAGES**. It supports both unbiased and biased simulations. It has one argument:
  - `--unbiased`: run an unbiased simulation if is set `True`.
- **readtrajtoexyz.py**
  - converts trajectory from **ASE** `traj` format to `extxyz`.

### Vacuum

This directory contains subfolders for generating the training data for ADP in vacuum.

- **spc-dft**
  - a python script for single point energy calculations using the **CP2K** calculator in **ASE**.
- **unbiased**
  - **LAMMPS** input file for unbiased simulations
  - **run.py** runs an unbiased MD simulation using **LAMMPS** and **PySAGES**.

- **uniform**
  - restart and parameter files for enhanced sampling simulations
  - **pysages-ala.py** runs a MD simulation using **ASE** and **PySAGES**.

