"""
A script that will produce a simulation with N bifurcations. Every branch
will have the same pseudotime length (50), and hyperparameters for all sanity
checks and count sampling take default values.
"""

import warnings
import argparse

import numpy as np
from numpy import random
import pandas as pd
import scipy as sp

import matplotlib.pyplot as plt
import pylab

from prosstt import simulation as sim
from prosstt import tree
from prosstt import sim_utils as sut

def maxes(x):
    res = np.zeros(len(x.keys()))
    for i, b in enumerate(x.keys()):
        res[i] = np.max(x[b])
    return res

def save_files(job_id, save_dir, X, labs, brns, scalings, uMs, H, gene_scale, alpha, beta):
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
    gene_params = pd.DataFrame({"alpha": alpha,
                                "beta": beta,
                                "genescale": gene_scale},
                               index=gene_names,
                               columns=["alpha", "beta", "genescale"])

    expr_mat.to_csv(save_dir + "/" + job_id + "_simulation.txt", sep="\t")
    cell_params.to_csv(save_dir + "/" + job_id + "_cellparams.txt", sep="\t")
    gene_params.to_csv(save_dir + "/" + job_id + "_geneparams.txt", sep="\t")

    np.savetxt(fname=save_dir + "/" + job_id + "_h.txt", X=H)

    for branch in uMs.keys():
        np.savetxt(fname=save_dir + "/" + job_id + "_ums" + str(branch) + ".txt",
                   X=uMs[branch])

def save_params(job_id, save_dir, lineage_tree, rseed):
    paramfile = save_dir + "/" + job_id + "_params.txt"
    with open(paramfile, 'w') as out:
        out.write("Genes: " + str(lineage_tree.G) + "\n")
        out.write("pseudotimes: " + str(list(lineage_tree.time.values)) + "\n")
        out.write("topology: " + str(lineage_tree.topology) + "\n")
        out.write("#modules: " + str(lineage_tree.modules) + "\n")
        out.write("random seed: " + str(rseed))

def main(job_id, save_dir, num_brpoints):
    rseed = np.random.randint(100000)
    random.seed(rseed)

    # sample the parameters randomly:
    G = random.randint(100, 1001)

    alpha = np.exp(random.normal(loc=np.log(0.2), scale=np.log(1.5), size=G))
    beta = np.exp(random.normal(loc=np.log(1), scale=np.log(1.5), size=G)) + 1

    num_branches = 2 * num_brpoints + 1
    modules = 5 * num_brpoints + np.random.randint(3, 20)
    top = tree.Tree.gen_random_topology(num_brpoints)

    branches = np.unique(np.array(top).flatten())
    time = {b: 50 for b in branches}

    # define tree
    t = tree.Tree(topology=top, time=time, num_branches=num_branches,
                  G=G, branch_points=num_brpoints, modules=modules)

    # pick a for parameters according to the number of modules, since
    # more modules mean more regulatory choices --> you can't have all
    # programs firing from the start
    mya = np.min([0.05, 1 / t.modules])

    uMs, Ws, H = sim.simulate_lineage(t, a=mya, intra_branch_tol=0, inter_branch_tol=0)
    # test = maxes(uMs)
    # print(test)

    gene_scale = sut.simulate_base_gene_exp(t, uMs, abs_max=10000)
    Ms = {}
    for branch in t.branches:
        Ms[branch] = np.exp(uMs[branch]) * gene_scale
    t.add_genes(Ms)

    X, pseudotime, brns, scalings = sim.sample_density(
        t, t.num_branches * 50, alpha=alpha, beta=beta)
    # X, pseudotime, brns, scalings = sim.sample_whole_tree(t, 1, alpha=alpha, beta=beta)

    save_params(job_id, save_dir, t, rseed)
    save_files(job_id, save_dir, X, pseudotime, brns, scalings, uMs, H,
               gene_scale, alpha, beta)
    # print(rseed, ":", np.max(X))


if __name__ == "__main__":
    PARSER = argparse.ArgumentParser(description='Generate a simulated scRNAseq \
                                     dataset with an n-fold bifurcation.')
    PARSER.add_argument("-j", "--job", dest="job", help="Job ID (prepended to \
                        all generated files)", metavar="FILE")
    PARSER.add_argument("-o", "--out", dest="outdir", help="Directory where  \
                        output files are saved", metavar="DIR")
    PARSER.add_argument("-n", "--num_brpoints", dest="n",
                        help="How many branching points the simulation contains",
                        metavar="INT")

    args = PARSER.parse_args()

    main(args.job, args.outdir, int(args.n))
