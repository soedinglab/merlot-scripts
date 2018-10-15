#!/usr/bin/bash

mscripts=$1 # location of merlot-scripts
job=$2      # name of the simulation, e.g. "test34"
out=$3      # output folder for the particular simulation
n=$4        # number of bifurcations

source activate py36
python "$mscripts/generate_simN.py" -j "$job" -o "$out" -n "$n"
source deactivate py36