#!/usr/bin/bash -eux

mscripts=$1
scripts=$1/scripts
out=$2
job=$3
dim=$4

timefile="${out}/${job}".exec
################################################################################
# final benchmark that keeps changing                                          #
################################################################################

# this seems to help, even though I deactivated the environment previously too
source deactivate py36
touch "${timefile}"

# destiny in log space with and without k optimization
name="destiny_log_k"
echo "${name}" > "${timefile}"
{ time timeout 60m Rscript "${scripts}"/benchmark_destiny.R -o "${out}"/ -j "${job}" -d "${dim}" -l -n -t "$mscripts"; } 2>> "${timefile}"
rc=$?
if [[ $rc == 124 ]]; then
    Rscript "${scripts}"/benchmark_stopped.R "${out}"/ "${job}" "${name}"
fi

n="-n"
for f in "" "-f";
do
    for b in "" "-b";
    do
        for sel in "none";# "merlot" "monocle" "seurat"
        do
            for reduced in "none" "knn";
            do
                # echo "benchmark_MERLoT_dest log k fixed interp=none" $b
                if [ "$b" ]; then emb="emb"; else emb="el"; fi
                if [ "$f" ]; then fixed="fixed"; else fixed="free"; fi
                name="LPGraph_log_k_"$fixed"_"$emb"_"$sel"_"$reduced
                echo "${name}" >> "${timefile}"
                { time timeout 60m Rscript "${scripts}"/benchmark_MERLoT_dest.R -o "${out}"/ -j "${job}" -d "${dim}" --log $b $f $n --sens --select "${sel}" -t "$mscripts" -r "$reduced"; } 2>> "${timefile}"
                rc=$?
                if [[ $rc == 124 ]]; then
                    Rscript "${scripts}"/benchmark_stopped.R "${out}"/ "${job}" "${name}"
                fi
            done
        done
    done
done