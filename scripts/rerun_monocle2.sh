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
    sim18
    sim5
    sim62
    sim64
    sim65
    sim66
    sim71
    sim7
    sim81
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
    sim16
    sim17
    sim25
    sim29
    sim2
    sim35
    sim37
    sim39
    sim44
    sim46
    sim50
    sim61
    sim63
    sim65
    sim69
    sim75
    sim76
    sim81
    sim91
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
    sim0
    sim10
    sim17
    sim18
    sim19
    sim1
    sim22
    sim24
    sim28
    sim30
    sim34
    sim35
    sim3
    sim42
    sim44
    sim45
    sim46
    sim50
    sim53
    sim57
    sim58
    sim61
    sim62
    sim63
    sim67
    sim68
    sim73
    sim79
    sim86
    sim87
    sim90
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
