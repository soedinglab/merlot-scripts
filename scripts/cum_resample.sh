#!/usr/bin/bash -eux

module load R/3.5.0
mscripts="$1"   # location of merlot-scripts repository
job="$2"      # name of each simulation folder
input="$3"    # name of input folder (sim to resample)
out="$4"      # name of output folder
n="$5"        # number of bifurcations


dim=$((n + 1)) # number of allowed dimensions
dim=$((dim>2?dim:2))

scripts=$mscripts/scripts

# resample
bash "${scripts}/resample.sh" "$scripts" "$job" "$out" "$n" "$input"

# run diffusion maps but check if they exist first
if [[ ! -f "${out}"/"${job}"_destiny_log_k ]]; then
  Rscript "${scripts}/run_destiny.R" -o "${out}"/ -j "${job}" -d "${dim}" -n -l -s none
else
  echo "${out}/${job}_destiny_log_k already present!"
fi

if [[ ! -f "${out}/${job}_destiny_log" ]]; then
  Rscript "${scripts}/run_destiny.R" -o "${out}/" -j "${job}" -d "${dim}" -l -s none
else
  echo "${out}"/"${job}"_destiny_log already present!
fi

# perform local averaging for the diffusion maps
# but check if they exist first!
if [[ ! -f "${out}"/"${job}"_destiny_log_knn.txt ]]; then
  bash "${scripts}"/run_reduce.sh "${scripts}" "${out}"/"${job}"_destiny_log "${dim}"
else
  echo "${out}"/"${job}"_destiny_log_knn.txt already present!
fi

if [[ ! -f "${out}"/"${job}"_destiny_log_k_knn.txt ]]; then
  bash "${scripts}"/run_reduce.sh "${scripts}" "${out}"/"${job}"_destiny_log_k "${dim}"
else
  echo "${out}"/"${job}"_destiny_log_k_knn.txt already present!
fi

# run monocle
# Rscript "${scripts}"/run_monocle.R -o "${out}"/ -j "${job}" --unconstrained -s none
# Rscript "${scripts}"/run_monocle.R -o "${out}"/ -j "${job}" -d "${dim}" -s none

# run predictions and benchmark them
bash "${scripts}"/timed_resamplesN.sh "${mscripts}" "${out}"/ "${job}" "${dim}"

