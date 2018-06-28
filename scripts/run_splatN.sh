#!/usr/bin/bash

mscripts=$1
job=$2
out=$3
matH=$4
n=$5

Rscript $mscripts/scripts/generate_splatN.R -j $job -o $out -n $n
