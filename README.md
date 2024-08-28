# Re-evaluating Benchmarks for MLIPs

This repository contains files related to the data generation, training, and evaluation for MLIP models. The repository is organized as follows:

## Directory Structure

- [**Dataset Generation/**](#dataset-generation)
- [**MLIP Training/**](#mlip-training)
  - [**Allegro MLIP Models**](#allegro-mlip-models)
  - [**DeePMD MLIP Models**](#deepmd-mlip-models)
- [**MLIP Evaluation/**](#mlip-evaluation)
  - [**Allegro MLIP Models**](#allegro-mlip-evaluation)
  - [**DeePMD MLIP Models**](#deepmd-mlip-evaluation)

Go to each subfolder for a detailed description.

## Dataset Generation

The `generate-datasets/` directory contains scripts to generate datasets for MLIP model training.

## MLIP Training

### Allegro MLIP Models

The `training-allegro/` directory provides datasets and sample configuration files for the training process of **Allegro** MLIP models.

### DeePMD MLIP Models

The `training-deepmd/` directory provides datasets and sample configuration files for the training process of **DeePMD** MLIP models.

## MLIP Evaluation

### Allegro MLIP Models

The `evaluate-allegro/` directory contains scripts to perform FES estimation with **Allegro**.

### DeePMD MLIP Models

The `evaluate-deepmd/` directory contains scripts to perform FES estimation with **DeePMD**.
