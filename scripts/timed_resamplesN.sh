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
# name="destiny_log_k"
# echo "${name}"
# echo "${name}" > "${timefile}"
# { time Rscript "${scripts}"/benchmark_destiny.R -o "${out}"/ -j "${job}" -d "${dim}" -l -n -t "$mscripts"; } 2>> "${timefile}"
# rc=$?
# if [[ $rc == 124 ]]; then
#     Rscript "${scripts}"/benchmark_stopped.R "${out}"/ "${job}" "${name}"
# fi

# for f in "" "-f";
# do
#     for b in "" "-b";
#     do
#         for n in "" "-n";
#         do
#             for elpi in "" "--elpi";
#             do
#                 # echo "benchmark_MERLoT_dest log k fixed interp=none" $b
#                 if [ -n "$b" ]; then emb="emb"; else emb="el"; fi
#                 if [ -n "$f" ]; then fixed="fixed"; else fixed="free"; fi
#                 if [ -n "$elpi" ]; then par="elpi"; else par="merlot"; fi
#                 name="MERLoT_log_k_${fixed}_${emb}_none_knn_${par}"
#                 echo "${name}" >> "${timefile}"
#                 echo "${name}"
#                 PARAMS="${b} ${f} ${n} ${elpi}"
#                 echo Rscript "${scripts}"/benchmark_MERLoT_dest.R -t "${mscripts}" -o "${out}"/ -j "${job}" -d "${dim}" --log $PARAMS --sens --select none -r knn
#                 { time Rscript "${scripts}"/benchmark_MERLoT_dest.R -t "${mscripts}" -o "${out}"/ -j "${job}" -d "${dim}" --log $PARAMS --sens --select none -r knn; } 2>> "${timefile}"
#                 rc=$?
#                 if [[ $rc == 124 ]]; then
#                     Rscript "${scripts}"/benchmark_stopped.R "${out}"/ "${job}" "${name}"
#                 fi
#             done
#         done
#     done
# done

# # slingshot
# possible=$(find "${out}" -maxdepth 1 -type f | grep "\." -v | grep "${job}")
# for line in $possible; do
#     name="slingshot "$line
#     echo "${name}" >> "${timefile}"
#     { time Rscript "${scripts}"/benchmark_slingshot.R -o "${out}"/ -j "${job}" -d "${dim}" -i "${line}" -t "$mscripts"; } 2>> "${timefile}"
#     rc=$?
#     if [[ $rc == 124 ]]; then
#         Rscript "${scripts}"/benchmark_stopped.R "${out}"/ "${job}" "${name}"
#     fi
# done

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
