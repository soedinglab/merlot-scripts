#!/usr/bin/bash -eux

mscripts=$1 # location of merlot-scripts repository
batch=$2    # name of each simulation, e.g. "test", "sim"
out=$3      # output folder
cores=$4    # how many cores to use
n=$5        # number of bifurcations

# my python stuff is in conda
# if you have R in conda you need to make sure that the version works for
# the libraries we need
# if you are not using conda, remove the offending lines

# benchmark plan:
# 1. simulate a number of datasets with $n bifurcations.
# 2. run predictions on them, save files
# 3. benchmark predictions

# how many datasets do we produce
many=100

# create the input file for parallel
: > "$out/par_params"
for ((i=0; i<"$many"; i++))
do
  name="$batch$i"
  mkdir "$out/$name/"
  echo "$name" "$out/$name" >> "$out/par_params"
done

# unnecessary if simulations are present
# parallel -j "$cores" --colsep " " --results "$out/" -a "$out/par_params" bash "$mscripts/run_simN.sh" "$mscripts" "{1}" "{2}" "$n"
parallel -j "$cores" --colsep " " --results "$out/" -a "$out/par_params" bash "$mscripts/run_as_one.sh" "$mscripts" "{1}" "{2}" "$n"

# parse all scaffold trees and record the number of endpoints, branch points and branches
# Rscript ~/Documents/repos/mscripts/parse_scaffolds.R -o "$out" -b "$out"