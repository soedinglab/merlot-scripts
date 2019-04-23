#!/usr/bin/bash

mscripts=$1 # location of merlot-scripts
job=$2      # name of the simulation, e.g. "test34"
out=$3      # output folder for the particular simulation
n=$4        # number of bifurcations
input=$5    # name of the simulation we are resampling

source activate py36
python "$mscripts"/scripts/resample_simN.py -j "$job" -o "$out" -n "$n" -i "$input"
source deactivate py36