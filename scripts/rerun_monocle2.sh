#!/usr/bin/bash -eux

outfile=$1
touch "$outfile"
mscripts="/home/users/npapado/bin/merlot-scripts"

###############################
###############################
###############################
bench="benchmark7"
dim=8
sims=(
    sim73
    sim23
    sim37
    sim61
    sim96
    sim90
    sim11
    sim17
    sim33
    sim77
    sim46
    sim64
    sim48
    sim56
    sim82
    sim32
    sim53
    sim3
    sim75
    sim95
    sim24
    sim65
    sim93
    sim14
    sim36
)

out="/cbscratch/npapado/resample/$bench"
for sim in "${sims[@]}";
do
    input="/cbscratch/npapado/resample/$bench/$sim/${sim}_destiny_log_k"
    timefile="${out}/${sim}/${sim}.exec"
    echo "time Rscript $mscripts/scripts/benchmark_monocl2.R -o ${out}/${sim}/ -j ${sim} -d ${dim} -t ${mscripts}; 2>> ${timefile} &" >> "$outfile"
done

###############################
###############################
###############################
bench="benchmark8"
dim=9
sims=(
    sim12
    sim15
    sim18
    sim21
    sim22
    sim27
    sim28
    sim31
    sim34
    sim35
    sim38
    sim3
    sim40
    sim41
    sim42
    sim44
    sim45
    sim47
    sim48
    sim49
    sim56
    sim5
    sim62
    sim63
    sim64
    sim65
    sim66
    sim67
    sim71
    sim74
    sim75
    sim78
    sim79
    sim7
    sim81
    sim82
    sim83
    sim84
    sim85
    sim86
    sim8
    sim92
    sim93
    sim94
)

out="/cbscratch/npapado/resample/$bench"
for sim in "${sims[@]}";
do
    input="/cbscratch/npapado/resample/$bench/$sim/${sim}_destiny_log_k"
    timefile="${out}/${sim}/${sim}.exec"
    echo "time Rscript $mscripts/scripts/benchmark_monocl2.R -o ${out}/${sim}/ -j ${sim} -d ${dim} -t ${mscripts}; 2>> ${timefile} &" >> "$outfile"
done

###############################
###############################
###############################
bench="benchmark9"
dim=10
sims=(
    sim10
    sim11
    sim12
    sim13
    sim15
    sim16
    sim17
    sim23
    sim24
    sim25
    sim27
    sim29
    sim2
    sim30
    sim32
    sim33
    sim35
    sim37
    sim39
    sim42
    sim43
    sim44
    sim46
    sim48
    sim50
    sim52
    sim54
    sim55
    sim5
    sim61
    sim63
    sim64
    sim65
    sim66
    sim67
    sim69
    sim75
    sim76
    sim79
    sim7
    sim80
    sim81
    sim84
    sim86
    sim88
    sim91
    sim92
    sim93
    sim9
)

out="/cbscratch/npapado/resample/$bench"
for sim in "${sims[@]}";
do
    input="/cbscratch/npapado/resample/$bench/$sim/${sim}_destiny_log_k"
    timefile="${out}/${sim}/${sim}.exec"
    echo "time Rscript $mscripts/scripts/benchmark_monocl2.R -o ${out}/${sim}/ -j ${sim} -d ${dim} -t ${mscripts}; 2>> ${timefile} &" >> "$outfile"
done

###############################
###############################
###############################
bench="benchmark10"
dim=11
sims=(
    sim0
    sim10
    sim13
    sim14
    sim15
    sim17
    sim18
    sim19
    sim1
    sim20
    sim22
    sim23
    sim24
    sim26
    sim28
    sim29
    sim2
    sim30
    sim32
    sim34
    sim35
    sim36
    sim37
    sim3
    sim42
    sim44
    sim45
    sim46
    sim47
    sim49
    sim50
    sim53
    sim55
    sim57
    sim58
    sim59
    sim5
    sim60
    sim61
    sim62
    sim63
    sim65
    sim66
    sim67
    sim68
    sim6
    sim72
    sim73
    sim75
    sim79
    sim80
    sim83
    sim85
    sim86
    sim87
    sim89
    sim90
    sim91
    sim92
    sim93
    sim94
    sim9
)

out="/cbscratch/npapado/resample/$bench"
for sim in "${sims[@]}";
do
    input="/cbscratch/npapado/resample/$bench/$sim/${sim}_destiny_log_k"
    timefile="${out}/${sim}/${sim}.exec"
    echo "time Rscript $mscripts/scripts/benchmark_monocl2.R -o ${out}/${sim}/ -j ${sim} -d ${dim} -t ${mscripts}; 2>> ${timefile} &" >> "$outfile"
done


split -l 16 --additional-suffix .sh -d "$outfile" sub

for sub in sub*.sh; do
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
