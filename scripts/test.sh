#!/usr/bin/bash -eux

#BSUB -q mpi-long+
#BSUB -o log.%J
#BSUB -e log.%J
#BSUB -W 48:00
#BSUB -n 16
#BSUB -m hh
#BSUB -R cbscratch
#BSUB -R "span[ptile=16]"

time Rscript /home/users/npapado/bin/merlot-scripts/scripts/benchmark_slingshot.R -o /cbscratch/npapado/resample/benchmark1 -j sim3 -d 2 -i /cbscratch/npapado/resample/benchmark1/sim3/sim3_destiny_log_k -t /home/users/npapado/bin/merlot-scripts; 2>> /cbscratch/npapado/resample/benchmark1/sim3/sim3.exec