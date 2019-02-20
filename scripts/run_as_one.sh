#!/usr/bin/bash -eux

mscripts=$1             # location of merlot-scripts
job=$2                  # name of the simulation, e.g. "test34"
out=$3                  # output folder for the particular simulation
n=$4                    # number of bifurcations
dim=$((n + 1))          # number of allowed dimensions. n+2 for splatter sims might yield better results.
dim=$((dim>2?dim:2))    # minimum number of dimensions if users somehow set this to 1? d=3 for splatter
                        # might yield better results.

scripts=$mscripts/scripts

# # run Destiny, needed for all methods that use diffusion maps as input.
# Rscript "${scripts}"/run_destiny.R -o "${out}"/ -j "${job}" -d "${dim}" -n -l -s "none"
# Rscript "${scripts}"/run_destiny.R -o "${out}"/ -j "${job}" -d "${dim}" -l -s "none"

# perform local averaging for the diffusion maps
source activate py36
python "${scripts}"/reduce.py -i "${out}"/"${job}"_destiny_log -d "${dim}"
python "${scripts}"/reduce.py -i "${out}"/"${job}"_destiny_log_k -d "${dim}"
source deactivate py36

# # run monocle2, needed for all methods that use DDRTree as input.
# Rscript "${scripts}"/run_monocle.R -o "${out}"/ -j "${job}" --unconstrained -s "none"
# Rscript "${scripts}"/run_monocle.R -o "${out}"/ -j "${job}" -d "${dim}" -s "none"

# predict trees and evaluate them against labels
# bash "${scripts}"/timed_benchmarksN.sh "${mscripts}" "${out}"/ "${job}" "${dim}"

# read the results and plot them:
# Rscript "${scripts}"/plot_predictions.R $mscripts "${out}"/ "${job}"