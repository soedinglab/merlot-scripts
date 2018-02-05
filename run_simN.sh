#!/usr/bin/bash

mscripts=$1
job=$2
out=$3
n=$4

source activate py36
python "$mscripts/generate_simN.py" -j "$job" -o "$out" -n "$n"
source deactivate py36