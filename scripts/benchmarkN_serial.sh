#!/usr/bin/bash -eux

mscripts=$1   # location of merlot-scripts repository
batch=$2    # name of each simulation folder
out=$3      # name of folder for parallel output
many=$4     # how many simulations to perform
cores=$5    # how many cores to use
n=$6        # number of bifurcations
sim=$7      # which simulation suite to use (if empty only runs evaluation)

# general benchmark plan:
# 1. simulate $many scRNA-seq datasets for $n bifurcations. Save each
#    in its own folder.
# 2. run predictions and evaluate them

# make the folders
mkdir "${out}"
: > "${out}"/par_params

# make the parameter file we need for parallel
for ((i=0; i<"$many"; i++))
do
  name=${batch}${i}
  mkdir "$out"/"$name"/
  echo "$name" "${out}/${name}" >> "${out}"/par_params
done

# make the simulations
if [ "$sim" = "splatter" ]; then
  # echo "splatter"
  parallel -j "$cores" --colsep " " --results "$out"/ -a "$out"/par_params bash "$mscripts"/scripts/run_splatN.sh "$mscripts" {1} {2} "$n"
elif [ "$sim" = "PROSSTT" ]; then
  # echo "prosstt"
  parallel -j "$cores" --colsep " " --results "$out"/ -a "$out"/par_params bash "$mscripts"/scripts/run_simN.sh "$mscripts" {1} {2} "$n"
fi

# run predictions and evaluate them
parallel -j "$cores" --colsep " " --results "$out"/ -a "$out"/par_params bash "$mscripts"/scripts/run_as_one.sh "$mscripts" {1} {2} "$n"

# parse all scaffold trees and record the number of endpoints, branch points and branches.
# This is optional - only run if you want to see the performance of the scaffold trees.
# Rscript "$mscripts"/scripts/parse_scaffolds.R -o "$out" -b "$out"