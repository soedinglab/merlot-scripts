#!/usr/bin/bash -eux

mscripts=$1       # location of merlot-scripts repository
batch=$2        # name of each simulation folder
out=$3          # name of overall output folder
sub_file=$4     # the big file that holds all of the jobs to be submitted
res_source=$5   # the source for the resampling

mkdir -p "$out"
:> "$sub_file"

for b in {1..10}
do
  benchmark="$out"/benchmark$b
  mkdir -p "$benchmark"
  for i in {0..29}
  do
    name=${batch}${b}
    output="${benchmark}/${name}/"
    mkdir -p "${output}"
    input="${res_source}/benchmark${i}/${name}/${name}"
    echo "$mscripts/scripts/cum_resample.sh $mscripts $name $input $output $b &" >> "$sub_file"
  done
done

split -l 16 --additional-suffix .sh -d "$sub_file" sub

for sub in sub*.sh; do
  code=$(cat "$sub")
  echo "#!/usr/bin/bash -eux" > "$sub"
  {
    echo ""
    echo "#BSUB -q mpi-long+"
    echo "#BSUB -o log.%J"
    echo "#BSUB -e log.%J"
    echo "#BSUB -W 48:00"
    echo "#BSUB -n 16"
    echo "#BSUB -m hh"
    echo "#BSUB -R cbscratch"
    echo "#BSUB -R \"span[ptile=16]\""
    echo ""
    echo "$code"
    echo ""
    echo "wait"
  } >> "$sub"
done
