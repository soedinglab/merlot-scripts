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

name="destiny_log"
echo "${name}" >> "${timefile}"
{ time timeout 60m Rscript "${scripts}"/benchmark_destiny.R -o "${out}"/ -j "${job}" -d "${dim}" -l -t "$mscripts"; } 2>> "${timefile}"
rc=$?
if [[ $rc == 124 ]]; then
    Rscript "${scripts}"/benchmark_stopped.R "${out}"/ "${job}" "${name}"
fi

# SLICER, as-is
name="SLICER"
echo "${name}" >> "${timefile}"
{ time timeout 60m Rscript "${scripts}"/benchmark_SLICER.R -o "${out}"/ -j "${job}" -d "${dim}" -t "$mscripts"; } 2>> "${timefile}"
rc=$?
if [[ $rc == 124 ]]; then
    Rscript "${scripts}"/benchmark_stopped.R "${out}"/ "${job}" "${name}"
fi

# monocle2, as-is
name="monocl2"
echo "${name}" >> "${timefile}"
{ time timeout 60m Rscript "${scripts}"/benchmark_monocl2.R -o "${out}"/ -j "${job}" -d "${dim}" -t "$mscripts"; } 2>> "${timefile}"
rc=$?
if [[ $rc == 124 ]]; then
    Rscript "${scripts}"/benchmark_stopped.R "${out}"/ "${job}" "${name}"
fi

# monocle2 with only 2 dimensions
name="monocl2_unc"
echo "${name}" >> "${timefile}"
{ time timeout 60m Rscript "${scripts}"/benchmark_monocl2.R -o "${out}"/ -j "${job}" --unconstrained -t "$mscripts"; } 2>> "${timefile}"
rc=$?
if [[ $rc == 124 ]]; then
    Rscript "${scripts}"/benchmark_stopped.R "${out}"/ "${job}" "${name}"
fi

# MERLoT+monocle2 in the following combinations:
# free, elastic
# free, embedded
# fixed, elastic
# fixed, embedded
# always in sensitive mode
for f in "" "-f";
do
    for b in "" "-b";
    do
        for u in "" "-u";
        do
            for sel in "none"; #"merlot" "monocle" "seurat";
            do
                if [ "$b" ]; then emb="emb"; else emb="el"; fi
                if [ "$f" ]; then fixed="fixed"; else fixed="free"; fi
                if [ "$u" ]; then unc="unc_"; else unc=""; fi
                name="monTree_"$unc""$fixed"_"$emb"_"$sel
                echo "${name}" >> "${timefile}"
                echo "monTree" $f $b $u
                { time timeout 60m Rscript "${scripts}"/benchmark_MERLoT_mon.R -o "${out}"/ -j "${job}" -d "${dim}" $b $f $u --sens --select "${sel}" -t "$mscripts"; } 2>> "${timefile}"
                rc=$?
                if [[ $rc == 124 ]]; then
                    Rscript "${scripts}"/benchmark_stopped.R "${out}"/ "${job}" "${name}"
                fi
            done
        done
    done
done

# MERLoT+destiny (log, k-opt), in the following combinations:
# free elastic
# free embedded
# fixed elastic
# fixed embedded
n="-n"
for f in "" "-f";
do
    for b in "" "-b";
    do
        for sel in "none";# "merlot" "monocle" "seurat"
        do
            echo "benchmark_MERLoT_dest log k fixed interp=none" $b
            if [ "$b" ]; then emb="emb"; else emb="el"; fi
            if [ "$f" ]; then fixed="fixed"; else fixed="free"; fi
            name="LPGraph_log_k_"$fixed"_"$emb"_"$sel
            echo "${name}" >> "${timefile}"
            { time timeout 60m Rscript "${scripts}"/benchmark_MERLoT_dest.R -o "${out}"/ -j "${job}" -d "${dim}" --log $b $f $n --sens --select "${sel}" -t "$mscripts"; } 2>> "${timefile}"
            rc=$?
            if [[ $rc == 124 ]]; then
                Rscript "${scripts}"/benchmark_stopped.R "${out}"/ "${job}" "${name}"
            fi
        done
    done
done

# TSCAN
name="TSCAN"
echo "${name}" >> "${timefile}"
{ time timeout 60m Rscript "${scripts}"/benchmark_TSCAN.R -o "${out}"/ -j "${job}" -t "$mscripts"; } 2>> "${timefile}"
rc=$?
if [[ $rc == 124 ]]; then
    Rscript "${scripts}"/benchmark_stopped.R "${out}"/ "${job}" "${name}"
fi

# slingshot
possible=$(find "${out}" -maxdepth 1 -type f | grep "\." -v | grep "${job}")
for line in $possible; do
    name="slingshot "$line
    echo "${name}" >> "${timefile}"
    { time timeout 60m Rscript "${scripts}"/benchmark_slingshot.R -o "${out}"/ -j "${job}" -d "${dim}" -i "${line}" -t "$mscripts"; } 2>> "${timefile}"
    rc=$?
    if [[ $rc == 124 ]]; then
        Rscript "${scripts}"/benchmark_stopped.R "${out}"/ "${job}" "${name}"
    fi
done
