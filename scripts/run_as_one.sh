#!/usr/bin/bash -eux

mscripts=$1
job=$2
out=$3
matH=$4 # relic from previous versions, left to not break script calls
n=$5 # number of branch points
dim=$((n + 1)) # number of allowed dimensions
dim=$((dim>2?dim:2))

scripts=$mscripts/scripts

# first perform gene selection the monocle and MERLoT way
# Rscript "${scripts}"/select_genes_merlot.R  -o "${out}"/ -j "${job}"
# Rscript "${scripts}"/select_genes_monocle.R -o "${out}"/ -j "${job}"
Rscript "${scripts}"/select_genes_seurat.R -o "${out}"/ -j "${job}" -t "${mscripts}"

# run diffusion maps, which we need for destiny and MERLoT+destiny
for gene_selection in "none" "seurat";
do
    Rscript "${scripts}"/run_destiny.R -o "${out}"/ -j "${job}" -d "${dim}" -n -l -s "${gene_selection}"
    Rscript "${scripts}"/run_destiny.R -o "${out}"/ -j "${job}" -d "${dim}" -l -s "${gene_selection}"
done

# run monocle2, needed for itself and MERLoT+monocle2. No need to run with "none" for gene selection,
# since monocle selection is the default anyway
for gene_selection in "none"; # "merlot" "monocle" "seurat";
do
    Rscript "${scripts}"/run_monocle.R -o "${out}"/ -j "${job}" --unconstrained -s "${gene_selection}"
    Rscript "${scripts}"/run_monocle.R -o "${out}"/ -j "${job}" -d "${dim}" -s "${gene_selection}"
done

# # last but not least PCA for slingshot
# Rscript "${scripts}"/run_PCA.R -o "${out}"/ -j "${job}"

# predict trees and evaluate them against labels
bash "${scripts}"/timed_benchmarksN.sh "${mscripts}" "${out}"/ "${job}" "${dim}"

# read the results and plot them:
# Rscript "${scripts}"/plot_predictions.R $mscripts "${out}"/ "${job}"