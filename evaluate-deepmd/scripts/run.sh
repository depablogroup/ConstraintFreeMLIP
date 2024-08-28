#!/bin/bash

# To run the simulation, designate a directory and run the following command
SAVEDIR="../explicit/unbiased40k/"  # Replace with your experiment of interest
python pysages-deepmd.py --savedir=$SAVEDIR
