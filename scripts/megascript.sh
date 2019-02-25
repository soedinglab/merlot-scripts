#!/usr/bin/bash -eux

module load R/3.5.0
# super mega script that does  E V E R Y T H I N G

hhtree=$1   # location of merlot-scripts repository
job=$2    # name of each simulation folder
out=$3      # name of output folder
n=$4        # number of bifurcations
sim=$5      # which simulation suite to use (if empty only runs evaluation)

scripts=$hhtree/scripts

# determine the simulations
if [ "$sim" = "splatter" ]; then
  # echo "splatter"
  simscript="$hhtree"/scripts/run_splatN.sh
  dim=$((n + 2)) # number of allowed dimensions
  dim=$((dim>2?dim:3))
elif [ "$sim" = "PROSSTT" ]; then
  simscript="$hhtree"/scripts/run_simN.sh
  dim=$((n + 1)) # number of allowed dimensions
  dim=$((dim>2?dim:2))
fi

# simulate
bash "$simscript" "$hhtree" "$job" "$out" bla "$n"

# run diffusion maps
Rscript "${scripts}"/run_destiny.R -o "${out}"/ -j "${job}" -d "${dim}" -n -l -s none
Rscript "${scripts}"/run_destiny.R -o "${out}"/ -j "${job}" -d "${dim}" -l -s none

# run monocle
Rscript "${scripts}"/run_monocle.R -o "${out}"/ -j "${job}" --unconstrained -s none
Rscript "${scripts}"/run_monocle.R -o "${out}"/ -j "${job}" -d "${dim}" -s none

# run predictions and benchmark them
bash "${scripts}"/timed_benchmarksN.sh "${hhtree}" "${out}"/ "${job}" "${dim}"

