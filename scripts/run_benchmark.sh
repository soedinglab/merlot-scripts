#!/usr/bin/bash -eux

# this script will run a benchmark with PROSSTT 100 simulations per bifurcation
# level on 25 cores. Each simulation will have the prefix "sim" and all simulations
# for one bifurcation level will be saved in some/dir/benchmark$n
mscripts="path/to/mscripts"

for n in 1 2 3 4 5 6 7 8 9 10;
do
    bash $mscripts/scripts/benchmarkN_serial.sh $mscripts sim some/dir/benchmark$n 100 25 $n PROSSTT
done