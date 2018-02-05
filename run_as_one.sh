#!/usr/bin/bash -eux

mscripts=$1 # location of merlot-scripts
job=$2      # name of the simulation, e.g. "test34"
out=$3      # output folder for the particular simulation
n=$4        # number of bifurcations

dim=$(expr $n + 1) # number of allowed dimensions
dim=$(($dim>2?$dim:2))


# # run diffusion maps, which we need for destiny and MERLoT+destiny
# Rscript "$mscripts/run_destiny.R" -o "$out/" -j "$job" -d "$dim" -n -l
# Rscript "$mscripts/run_destiny.R" -o "$out/" -j "$job" -d "$dim" -l

# # # run monocle2, needed for itself and MERLoT+monocle2
# Rscript "$mscripts/run_monocle.R" -o "$out/" -j "$job" -d "$dim"
# Rscript "$mscripts/run_monocle.R" -o "$out/" -j "$job" --unconstrained

# predict trees and evaluate them against labels
bash "$mscripts/timed_benchmarksN.sh" "$mscripts" "$out/" "$job" "$dim"