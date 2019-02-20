#!/usr/bin/bash -eux

mscripts=$1   # location of merlot-scripts repository
batch=$2      # name of each simulation folder
out=$3        # name of folder for parallel output
many=$4       # how many simulations to perform
cores=$5      # how many cores to use
n=$6          # number of bifurcations
resample=$7   # in which folder to look if we are resampling

# general benchmark plan:
# 1. simulate $many scRNA-seq datasets for $n bifurcations. Save each
#    in its own folder.
# 2. run predictions and evaluate them

# make the folders
mkdir "${out}"
mkdir "${resample}"
: > "${resample}"/par_params

# make the parameter file we need for parallel
for ((i=0; i<"$many"; i++))
do
  name=${batch}${i}
  mkdir "$out"/"$name"/
  mkdir "$resample"/"$name"/
  echo "$name" "${out}/${name}/${name}" "${resample}/${name}" >> "${resample}"/par_params
done

# parallel -j "$cores" --colsep " " --results "$out"/ -a "$resample"/par_params bash "$mscripts"/scripts/resample.sh "$mscripts" {1} {3} "$n" {2}

# run predictions and evaluate them
parallel -j "$cores" --colsep " " --results "$out"/ -a "$resample"/par_params bash "$mscripts"/scripts/run_as_one.sh "$mscripts" {1} {3} "$n"

# parse all scaffold trees and record the number of endpoints, branch points and branches.
# This is optional - only run if you want to see the performance of the scaffold trees.
# Rscript "$mscripts"/scripts/parse_scaffolds.R -o "$out" -b "$out"