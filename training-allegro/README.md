# Training

This directory contains the datasets and a sample configuration file for training **Allegro** model for ADP.

### Scripts

We provide a sample configuration for the allegro training. To replicate our training process, you should replace the dataset specified in the configuration file with the dataset you would like to use, as well as the `n_train` and `n_val` parameters.

- **test.yaml**: This is a sample configuration file for training **Allegro** model for ADP.

### Water + adp

This directory contains datasets with flexible water data. The directory is organized into 2 subfolders, each varying by sampling methods:

- **unbiased**
  - 10,000 snapshots generated from unbiased simulations at 600 K.
- **uniform**
  - 10,000 snapshots generated from uniform simulations at 600 K.

### Explicit

This directory contains datasets with force-only data. The directory is organized into 3 subfolders, each varying by sampling methods, or temperatures:

- **unbiased**
  - 40,000 snapshots generated from unbiased simulations at 600 K.
- **unbiased-300K**
  - 40,000 snapshots generated from uniform simulations at 300 K.
- **uniform**
  - 40,000 snapshots generated from uniform simulations at 600 K.

### Implicit

This directory contains models trained on ADP in implicit solvent. The directory is organized into 2 subfolders, each varying by sampling methods:

- **unbiased**
  - 1,000 snapshots generated from unbiased simulations at 600 K.
- **uniform**
  -  1,000 snapshots generated from uniform simulations at 600 K.

### Vacuum

This directory contains models trained on ADP in vacuum. The directory is organized into four subfolders, each varying by sampling methods or calculation levels:

- **unbiased**
  - 1,000 snapshots generated from unbiased simulations at 600 K.
- **unbiased-dft**
  - 1,000 snapshots generated from unbiased simulations at 600 K, calculated at the DFT level.
- **uniform**
  - 1,000 snapshots generated from uniform sampling at 600 K.
- **uniform-dft**
  - 1,000 snapshots generated from uniform sampling at 600 K, calculated at the DFT level.
