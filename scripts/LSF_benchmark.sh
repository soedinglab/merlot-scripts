#!/usr/bin/bash -eux

hhtree=$1     # location of merlot-scripts repository
batch=$2      # name of each simulation folder
out=$3        # name of overall output folder
sim=$4        # which simulation suite to use (if empty only runs evaluation)
howmany=$5    # the number of simulations to create per benchmark
benchmarks=$

mkdir -p $out
j=0

for ((i=0; i<"$howmany"; i++))
do
  if !(($i % 16 != 0)); then
    j=$((j + 1))
    sub_file="/usr/users/npapado/sub"${j}".sh"
    :> "$sub_file"
  fi
  name=${batch}${i}
  mkdir -p "$out"/"$name"/
  echo "$hhtree"/scripts/megascript.sh "$hhtree" "$name" "$out"/"$name"/ "$n" "$sim" >> "$sub_file"
done


