#!/usr/bin/bash -eux

outfile=$1
touch "$outfile"
mscripts="/home/users/npapado/bin/merlot-scripts"

###############################
###############################
###############################
bench="benchmark8"
dim=9
sims=(
    sim12
    sim62
    sim65
    sim71
    sim92
)

out="/cbscratch/npapado/resample/$bench"
for sim in "${sims[@]}";
do
    timefile="${out}/${sim}/${sim}.exec"
    # echo "Rscript ${mscripts}/scripts/run_monocle.R -o ${out}/${sim}/ -j ${sim} -d ${dim} -s none" >> "$outfile"
    echo "time Rscript $mscripts/scripts/benchmark_monocl2.R -o ${out}/${sim}/ -j ${sim} -d ${dim} -t ${mscripts} 2>> ${timefile} &" >> "$outfile"
done

###############################
###############################
###############################
bench="benchmark9"
dim=10
sims=(
    sim13
    sim2
    sim35
    sim75
)

out="/cbscratch/npapado/resample/$bench"
for sim in "${sims[@]}";
do
    timefile="${out}/${sim}/${sim}.exec"
    # echo "Rscript ${mscripts}/scripts/run_monocle.R -o ${out}/${sim}/ -j ${sim} -d ${dim} -s none" >> "$outfile"
    echo "time Rscript $mscripts/scripts/benchmark_monocl2.R -o ${out}/${sim}/ -j ${sim} -d ${dim} -t ${mscripts} 2>> ${timefile} &" >> "$outfile"
done

###############################
###############################
###############################
bench="benchmark10"
dim=11
sims=(
    sim17
    sim19
    sim34
    sim45
    sim46
    sim62
    sim91
)

out="/cbscratch/npapado/resample/$bench"
for sim in "${sims[@]}";
do
    timefile="${out}/${sim}/${sim}.exec"
    # echo "Rscript ${mscripts}/scripts/run_monocle.R -o ${out}/${sim}/ -j ${sim} -d ${dim} -s none" >> "$outfile"
    echo "time Rscript $mscripts/scripts/benchmark_monocl2.R -o ${out}/${sim}/ -j ${sim} -d ${dim} -t ${mscripts} 2>> ${timefile} &" >> "$outfile"
done


split -l 16 --additional-suffix .sh -d "$outfile" sub

for sub in ./sub*.sh; do
  code=$(cat "$sub")
  echo "#!/usr/bin/bash -eux" > "$sub"
  {
    echo ""
    echo "#BSUB -q mpi-long+"
    echo "#BSUB -o log.%J"
    echo "#BSUB -e log.%J"
    echo "#BSUB -W 100:00"
    echo "#BSUB -n 16"
    echo "#BSUB -m hh"
    echo "#BSUB -R cbscratch"
    echo "#BSUB -R \"span[ptile=16]\""
    echo ""
    echo "module load R/3.5.0"
    echo "$code"
    echo ""
    echo "wait"
  } >> "$sub"
done
