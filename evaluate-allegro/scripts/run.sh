#!/bin/bash

# To run the simulation, designate a directory and run the following command
SAVEDIR="../explicit/unbiased-1k/"  # Replace with your experiment of interest
python simulationpysages.py --savedir=$SAVEDIR
