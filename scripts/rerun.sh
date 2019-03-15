#!/usr/bin/bash -eux

outfile=$1
touch "$outfile"
mscripts="/home/users/npapado/bin/merlot-scripts"

bench="benchmark1"
dim=2
sims=(
    sim3
    sim56
)

out="/cbscratch/npapado/resample/$bench"
for sim in "${sims[@]}";
do
    input="/cbscratch/npapado/resample/$bench/$sim/${sim}_destiny_log_k"
    timefile="${out}/${sim}/${sim}.exec"
    echo "time Rscript $mscripts/scripts/benchmark_slingshot.R -o ${out}/${sim}/ -j ${sim} -d ${dim} -i $input -t ${mscripts}; 2>> ${timefile } &" >> "$outfile"
done

###############################
###############################
###############################
bench="benchmark5"
dim=6
sims=(
    sim11
    sim2
    sim34
    sim35
    sim51
    sim5
    sim66
    sim68
    sim73
    sim85
    sim92
)

out="/cbscratch/npapado/resample/$bench"
for sim in "${sims[@]}";
do
    input="/cbscratch/npapado/resample/$bench/$sim/${sim}_destiny_log_k"
    timefile="${out}/${sim}/${sim}.exec"
    echo "time Rscript $mscripts/scripts/benchmark_slingshot.R -o ${out}/${sim}/ -j ${sim} -d ${dim} -i $input -t ${mscripts}; 2>> ${timefile } &" >> "$outfile"
done

###############################
###############################
###############################
bench="benchmark6"
dim=7
sims=(
    sim0
    sim11
    sim12
    sim13
    sim17
    sim19
    sim22
    sim28
    sim32
    sim35
    sim36
    sim37
    sim39
    sim3
    sim40
    sim44
    sim46
    sim49
    sim50
    sim52
    sim53
    sim57
    sim59
    sim69
    sim6
    sim74
    sim75
    sim78
    sim80
    sim81
    sim84
    sim90
    sim93
)

out="/cbscratch/npapado/resample/$bench"
for sim in "${sims[@]}";
do
    input="/cbscratch/npapado/resample/$bench/$sim/${sim}_destiny_log_k"
    timefile="${out}/${sim}/${sim}.exec"
    echo "time Rscript $mscripts/scripts/benchmark_slingshot.R -o ${out}/${sim}/ -j ${sim} -d ${dim} -i $input -t ${mscripts}; 2>> ${timefile } &" >> "$outfile"
done

###############################
###############################
###############################
bench="benchmark7"
dim=8
sims=(
    sim10
    sim13
    sim16
    sim18
    sim19
    sim1
    sim21
    sim22
    sim23
    sim25
    sim27
    sim28
    sim29
    sim2
    sim31
    sim32
    sim35
    sim37
    sim40
    sim42
    sim46
    sim47
    sim48
    sim49
    sim4
    sim51
    sim5
    sim60
    sim63
    sim65
    sim67
    sim68
    sim69
    sim6
    sim72
    sim73
    sim74
    sim79
    sim7
    sim81
    sim82
    sim83
    sim85
    sim87
    sim88
    sim8
    sim91
    sim95
    sim96
    sim97
    sim99
)

out="/cbscratch/npapado/resample/$bench"
for sim in "${sims[@]}";
do
    input="/cbscratch/npapado/resample/$bench/$sim/${sim}_destiny_log_k"
    timefile="${out}/${sim}/${sim}.exec"
    echo "time Rscript $mscripts/scripts/benchmark_slingshot.R -o ${out}/${sim}/ -j ${sim} -d ${dim} -i $input -t ${mscripts}; 2>> ${timefile } &" >> "$outfile"
done

###############################
###############################
###############################
bench="benchmark8"
dim=9
sims=(
    sim10
    sim11
    sim12
    sim13
    sim14
    sim15
    sim16
    sim18
    sim19
    sim1
    sim22
    sim23
    sim24
    sim25
    sim27
    sim28
    sim29
    sim30
    sim31
    sim32
    sim34
    sim35
    sim36
    sim37
    sim39
    sim3
    sim40
    sim41
    sim44
    sim45
    sim4
    sim52
    sim53
    sim55
    sim56
    sim57
    sim58
    sim59
    sim5
    sim60
    sim61
    sim62
    sim64
    sim65
    sim66
    sim67
    sim68
    sim69
    sim70
    sim74
    sim75
    sim77
    sim78
    sim79
    sim80
    sim81
    sim83
    sim85
    sim86
    sim87
    sim88
    sim8
    sim90
    sim92
    sim93
    sim94
    sim95
    sim96
    sim97
    sim98
)

out="/cbscratch/npapado/resample/$bench"
for sim in "${sims[@]}";
do
    input="/cbscratch/npapado/resample/$bench/$sim/${sim}_destiny_log_k"
    timefile="${out}/${sim}/${sim}.exec"
    echo "time Rscript $mscripts/scripts/benchmark_slingshot.R -o ${out}/${sim}/ -j ${sim} -d ${dim} -i $input -t ${mscripts}; 2>> ${timefile } &" >> "$outfile"
done

###############################
###############################
###############################
bench="benchmark9"
dim=10
sims=(
    sim0
    sim10
    sim12
    sim13
    sim14
    sim15
    sim16
    sim17
    sim18
    sim21
    sim22
    sim24
    sim26
    sim27
    sim28
    sim29
    sim30
    sim31
    sim32
    sim33
    sim34
    sim35
    sim36
    sim37
    sim38
    sim39
    sim40
    sim41
    sim42
    sim43
    sim46
    sim47
    sim48
    sim49
    sim4
    sim50
    sim55
    sim56
    sim57
    sim58
    sim59
    sim5
    sim60
    sim61
    sim63
    sim64
    sim65
    sim66
    sim67
    sim68
    sim69
    sim6
    sim70
    sim71
    sim72
    sim73
    sim74
    sim75
    sim76
    sim77
    sim78
    sim79
    sim7
    sim80
    sim82
    sim83
    sim84
    sim85
    sim87
    sim89
    sim8
    sim91
    sim92
    sim93
    sim96
    sim97
    sim98
)

out="/cbscratch/npapado/resample/$bench"
for sim in "${sims[@]}";
do
    input="/cbscratch/npapado/resample/$bench/$sim/${sim}_destiny_log_k"
    timefile="${out}/${sim}/${sim}.exec"
    echo "time Rscript $mscripts/scripts/benchmark_slingshot.R -o ${out}/${sim}/ -j ${sim} -d ${dim} -i $input -t ${mscripts}; 2>> ${timefile } &" >> "$outfile"
done

###############################
###############################
###############################
bench="benchmark10"
dim=11
sims=(
    sim10
    sim11
    sim12
    sim13
    sim15
    sim17
    sim18
    sim20
    sim21
    sim22
    sim23
    sim24
    sim26
    sim28
    sim29
    sim30
    sim32
    sim33
    sim34
    sim35
    sim36
    sim37
    sim38
    sim39
    sim3
    sim40
    sim41
    sim42
    sim43
    sim44
    sim45
    sim46
    sim47
    sim48
    sim4
    sim50
    sim51
    sim52
    sim53
    sim54
    sim55
    sim56
    sim57
    sim58
    sim60
    sim61
    sim62
    sim63
    sim64
    sim65
    sim66
    sim68
    sim6
    sim70
    sim71
    sim73
    sim75
    sim76
    sim77
    sim79
    sim7
    sim80
    sim81
    sim82
    sim84
    sim85
    sim86
    sim87
    sim88
    sim89
    sim8
    sim90
    sim91
    sim93
    sim94
    sim95
    sim96
    sim97
    sim98
    sim99
    sim9
)

out="/cbscratch/npapado/resample/$bench"
for sim in "${sims[@]}";
do
    input="/cbscratch/npapado/resample/$bench/$sim/${sim}_destiny_log_k"
    timefile="${out}/${sim}/${sim}.exec"
    echo "time Rscript $mscripts/scripts/benchmark_slingshot.R -o ${out}/${sim}/ -j ${sim} -d ${dim} -i $input -t ${mscripts}; 2>> ${timefile } &" >> "$outfile"
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
