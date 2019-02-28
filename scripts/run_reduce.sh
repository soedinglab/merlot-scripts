#!/usr/bin/bash

mscripts=$1 # location of merlot-scripts
job=$2      # name of the input file
dim=$3      # number of dimensions

source activate py36
python "$mscripts/reduce.py" -i "$job" -d "${dim}"
source deactivate py36
