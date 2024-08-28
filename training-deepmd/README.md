# Training

This directory contains the datasets and a sample configuration file for training **DeePMD** model for ADP.

### Scripts

We provide a sample configuration for training. To replicate our training process, you should replace the dataset specified in the configuration file with the dataset you would like to use, as well as the `system` parameter.

- **input.json**: This is a sample configuration file for training **DeePMD** model for ADP.

### Water + adp

This directory contains datasets with flexible water data. The directory is organized into 2 subfolders, each varying by sampling methods:

- **unbiased**
  - 10,000 snapshots generated from unbiased simulations at 600 K.
- **uniform**
  - 10,000 snapshots generated from enhanced sampling at 600 K.

### Explicit

This directory contains datasets with force-only data. The directory is organized into 2 subfolders, each varying by sampling methods:

- **unbiased40k**
  - 40,000 snapshots generated from unbiased simulations at 600 K.
- **uniform40k**
  - 40,000 snapshots generated from enhanced sampling at 600 K.
