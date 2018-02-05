#!/usr/bin/bash -eux

mscripts=$1 # location of merlot-scripts
out=$2      # output folder for the particular simulation
job=$3      # name of the simulation, e.g. "test34"
dim=$4      # number of dimensions to be used

################################################################################
# final benchmark (?)                                                          #
################################################################################

# this seems to help, even though I deactivated the environment previously too
source deactivate py36

# destiny in log space with and without k optimization
name="destiny_log_k"
timeout 60m Rscript "$mscripts/benchmark_destiny.R" -o "$out/" -j "$job" -d "$dim" -l -n
rc=$?
if [[ $rc == 124 ]]; then
    Rscript "$mscripts/benchmark_stopped.R" "$out/" "$job" "$name"
fi

name="destiny_log"
timeout 60m Rscript "$mscripts/benchmark_destiny.R" -o "$out/" -j "$job" -d "$dim" -l
rc=$?
if [[ $rc == 124 ]]; then
    Rscript "$mscripts/benchmark_stopped.R" "$out/" "$job" "$name"
fi

# # SLICER, as-is
name="SLICER"
timeout 60m Rscript "$mscripts/benchmark_SLICER.R" -o "$out/" -j "$job" -d "$dim"
rc=$?
if [[ $rc == 124 ]]; then
    Rscript "$mscripts/benchmark_stopped.R" "$out/" "$job" "$name"
fi

# monocle2, as-is
name="monocl2"
timeout 60m Rscript "$mscripts/benchmark_monocl2.R" -o "$out/" -j "$job" -d "$dim"
rc=$?
if [[ $rc == 124 ]]; then
    Rscript "$mscripts/benchmark_stopped.R" "$out/" "$job" "$name"
fi

# monocle2 with only 2 dimensions
name="monocl2_unc"
timeout 60m Rscript "$mscripts/benchmark_monocl2.R" -o "$out/" -j "$job" --unconstrained
rc=$?
if [[ $rc == 124 ]]; then
    Rscript "$mscripts/benchmark_stopped.R" "$out/" "$job" "$name"
fi

# MERLoT+destiny (log, k-opt), in the following combinations:
# auto elastic
# auto embedded
# fixed elastic
# fixed embedded

source activate py36

n="-n"
for f in "" "-f";
do
    for b in "" "-b";
    do
        echo "benchmark_MERLoT_dest log k fixed interp=none" $b
        if [ "$b" ]; then emb="emb"; else emb="el"; fi
        if [ "$f" ]; then fixed="fixed"; else fixed="auto"; fi
        name="MERLoT_dest_log_k_${fixed}_${emb}"
        timeout 60m Rscript "$mscripts/benchmark_MERLoT_dest.R" -o "$out/" -j "$job" -d "$dim" --log $b $f $n --sens
        rc=$?
        if [[ $rc == 124 ]]; then
            Rscript "$mscripts/benchmark_stopped.R" "$out/" "$job" "$name"
        fi
    done
done

# MERLoT+monocle2 in the following combinations:
# auto, elastic
# auto, embedded
# fixed, elastic
# fixed, embedded
# always in sensitive mode
for f in "" "-f";
do
    for b in "" "-b";
    do
        for u in "" "-u";
        do
            if [ "$b" ]; then emb="emb"; else emb="el"; fi
            if [ "$f" ]; then fixed="fixed"; else fixed="auto"; fi
            if [ "$u" ]; then unc="_unc"; else unc=""; fi
            name="MERLoT_mon_${unc}${fixed}_${emb}"
            echo "MERLoT_mon" $f $b $u
            timeout 60m Rscript "$mscripts/benchmark_MERLoT_mon.R" -o "$out/" -j "$job" -d "$dim" $b $f $u --sens
            rc=$?
            if [[ $rc == 124 ]]; then
                Rscript "$mscripts/benchmark_stopped.R" "$out/" "$job" "$name"
            fi
        done
    done
done

source deactivate py36