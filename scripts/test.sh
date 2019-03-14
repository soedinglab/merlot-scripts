#!/usr/bin/bash -eux

#BSUB -q mpi-long+
#BSUB -o log.%J
#BSUB -e log.%J
#BSUB -W 48:00
#BSUB -n 16
#BSUB -m hh
#BSUB -R cbscratch
#BSUB -R "span[ptile=16]"

module load R/3.5.0
time Rscript /home/users/npapado/bin/merlot-scripts/scripts/benchmark_slingshot.R -o /cbscratch/npapado/resample/benchmark5/sim0/ -j sim0 -d 6 -i /cbscratch/npapado/resample/benchmark5/sim0/sim0_destiny_log_k -t /home/users/npapado/bin/merlot-scripts 2>> /cbscratch/npapado/resample/benchmark5/sim0/sim0.exec

