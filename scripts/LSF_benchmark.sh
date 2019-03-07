#!/usr/bin/bash -eux

hhtree=$1     # location of merlot-scripts repository
batch=$2      # name of each simulation folder
out=$3        # name of overall output folder
sim=$4        # which simulation suite to use (if empty only runs evaluation)
howmany=$5    # the number of simulations to create per benchmark
sub_file=$6   # the big file that holds all of the jobs to be submitted

mkdir -p "$out"
:> "$sub_file"

for b in 1 2 3 4 5 6 7 8 9 10;
do
  benchmark="$out"/benchmark$b
  mkdir -p "$benchmark"
  for ((i=0; i<"$howmany"; i++))
  do
    name=${batch}${i}
    mkdir -p "$benchmark"/"$name"
    echo "$hhtree"/scripts/megascript.sh "$hhtree" "$name" "$benchmark"/"$name"/ "$b" "$sim" "&" >> "$sub_file"
  done
done

split -l 16 --additional-suffix .sh -d "$sub_file" sub

for sub in sub*.sh; do
  test=$(cat "$sub")
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
    echo "$test"
    echo ""
    echo "wait"
  } >> "$sub"
done
