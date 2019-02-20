"""
A script that will read a simulation with N bifurcations and resample it
in much higher depth.
"""

import ast
import argparse

import numpy as np
from numpy import random
import pandas as pd

from prosstt import simulation as sim
from prosstt import tree

def save_files(job_id, save_dir, X, labs, brns, scalings):
    # make the data more presentable by adding gene and cell names
    cell_names = ["cell_" + str(i) for i in range(X.shape[0])]
    gene_names = ["gene_" + str(i) for i in range(X.shape[1])]

    expr_mat = pd.DataFrame(X, columns=gene_names,
                            index=cell_names).astype(int)
    cell_params = pd.DataFrame({"pseudotime": labs,
                                "branches": brns,
                                "scalings": scalings},
                               index=cell_names,
                               columns=["pseudotime", "branches", "scalings"])

    expr_mat.to_csv(save_dir + "/" + job_id + "_simulation.txt", sep="\t")
    cell_params.to_csv(save_dir + "/" + job_id + "_cellparams.txt", sep="\t")

def save_params(job_id, save_dir, lineage_tree, rseed):
    paramfile = save_dir + "/" + job_id + "_params.txt"
    with open(paramfile, 'w') as out:
        out.write("Genes: " + str(lineage_tree.G) + "\n")
        out.write("pseudotimes: " + str(list(lineage_tree.time.values)) + "\n")
        out.write("topology: " + str(lineage_tree.topology) + "\n")
        out.write("#modules: " + str(lineage_tree.modules) + "\n")
        out.write("random seed: " + str(rseed))

def read_paramfile(job):
    num_genes = 0
    pseudotimes = []
    topology = []
    modules = 0
    with open(job + "_params.txt", 'r') as paramfile:
        for line in paramfile:
            if line.startswith("Genes"):
                num_genes = int(line.rstrip().split()[-1])
            elif line.startswith("pseudotimes"):
                pt_string = line.rstrip().split(":")[-1].lstrip()
                pseudotimes = ast.literal_eval(pt_string)
            elif line.startswith("topology"):
                top_string = line.rstrip().split(":")[-1].lstrip()
                topology = ast.literal_eval(top_string)
            elif line.startswith("#modules"):
                modules = int(line.rstrip().split(" ")[-1])

    times = {}
    for branch, time in zip(np.sort(np.unique(topology)), pseudotimes):
        times[branch] = time
    return num_genes, times, topology, modules

def main(job_id, save_dir, input_dir, bif):
    rseed = np.random.randint(100000)
    random.seed(rseed)

    # read dataset and general parameters
    G, times, topology, modules = read_paramfile(input_dir)

    lineage = tree.Tree(topology=topology,
                        time=times,
                        num_branches=2*bif+1,
                        branch_points=bif,
                        modules=modules,
                        G=G,
                        root=0)

    # read gene expression
    geneparams = pd.read_csv(
        input_dir + "_geneparams.txt", sep="\t", index_col=0)
    gene_scale = np.array(geneparams["genescale"])
    Ms = {}
    for branch in lineage.branches:
        ums = np.loadtxt(input_dir + "_ums" + str(branch) + ".txt")
        Ms[branch] = np.exp(ums) * gene_scale

    lineage.add_genes(Ms)

    # draw cells and normalize
    X, pseudotime, brns, scalings = sim.sample_density(lineage,
                                                    (2*bif+1) * 500,
                                                    alpha=geneparams["alpha"],
                                                    beta=geneparams["beta"])
    X = (X.transpose() / scalings).transpose()

    # save new cell_params file
    save_params(job_id, save_dir, lineage, rseed)
    save_files(job_id, save_dir, X, pseudotime, brns, scalings)


if __name__ == "__main__":
    PARSER = argparse.ArgumentParser(description='Generate a simulated scRNAseq \
                                     dataset with an n-fold bifurcation.')
    PARSER.add_argument("-j", "--job", dest="job", help="Job ID (prepended to \
                        all generated files)", metavar="FILE")
    PARSER.add_argument("-i", "--in", dest="input", help="Input d \
                        all generated files)", metavar="FILE")
    PARSER.add_argument("-o", "--out", dest="outdir", help="Directory where  \
                        output files are saved", metavar="DIR")
    PARSER.add_argument("-n", "--num_brpoints", dest="n",
                        help="How many branching points the simulation contains",
                        metavar="INT")

    args = PARSER.parse_args()

    main(args.job, args.outdir, args.input, int(args.n))
