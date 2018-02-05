#!/usr/bin/env python
# coding: utf-8

import argparse

import numpy as np
from numpy import random
import pandas as pd
import scipy as sp

from prosstt import simulation as sim
from prosstt import tree


def save_files(job_id, save_dir, X, labs, brns, scalings, uMs, Hs, gene_scale, alpha, beta):
    # make the data more presentable by adding gene and cell names
    cell_names = ["cell_" + str(i) for i in range(X.shape[0])]
    gene_names = ["gene_" + str(i) for i in range(X.shape[1])]

    pdX = pd.DataFrame(X, columns=gene_names, index=cell_names).astype(int)
    pdCells = pd.DataFrame({"pseudotime": labs, "branches": brns, "scalings": scalings}, index=cell_names,
                           columns=["pseudotime", "branches", "scalings"])
    pdGenes = pd.DataFrame({"alpha": alpha, "beta": beta, "genescale": gene_scale}, index=gene_names,
                           columns=["alpha", "beta", "genescale"])

    pdX.to_csv(save_dir + "/" + job_id + "_simulation.txt", sep="\t")
    pdCells.to_csv(save_dir + "/" + job_id + "_cellparams.txt", sep="\t")
    pdGenes.to_csv(save_dir + "/" + job_id + "_geneparams.txt", sep="\t")

    num_branches = len(uMs)
    for i in range(num_branches):
        np.savetxt(fname=save_dir + "/" + job_id + "_ums" + str(i) + ".txt", X=uMs[i])
        np.savetxt(fname=save_dir + "/" + job_id + "_hs" + str(i) + ".txt", X=Hs[i])


def save_params(job_id, save_dir, G, br_lengths, br_compl, rseed, topology):
    paramfile = save_dir + "/" + job_id + "_params.txt"
    with open(paramfile, 'w') as out:
        out.write("Genes: " + str(G) + "\n")
        out.write("pseudotimes: " + str(br_lengths) + "\n")
        out.write("topology: " + str(topology) + "\n")
        out.write("#modules: " + str(br_compl) + "\n")
        out.write("random seed: " + str(rseed))


def main(job_id, save_dir, num_brpoints):
    rseed = np.random.randint(1000)
    random.seed(rseed)

    # sample the parameters randomly:
    G = random.random_integers(100, 1000)
    
    alpha = np.exp(random.normal(loc=np.log(0.2), scale=np.log(1.5), size=G))
    beta = np.exp(random.normal(loc=np.log(1), scale=np.log(1.5), size=G)) + 1

    num_branches = 2*num_brpoints + 1
    top = tree.Tree.gen_random_topology(num_brpoints)

    br_lengths = {}
    for branch in np.unique(top):
        br_lengths[branch] = 50

    br_compl = 10
    
    t = tree.Tree(topology=top, time=br_lengths, num_branches=num_branches,
                  branch_points=num_brpoints, modules=br_compl, G=G)

    sample_time = np.arange(0, t.get_max_time())

    Ms = {}
    while not sim.are_lengths_ok(Ms, abs_max=1000, rel_dif=0.05):
        uMs, Ws, Hs = sim.simulate_branching_data(t, tol=0.2)
        gene_scale = np.exp(sp.stats.norm.rvs(loc=0.8, scale=1, size=G))
        for branch in t.time.keys():
            Ms[branch] = np.exp(uMs[branch]) * gene_scale
    
    print(len(Ms))
    t.add_genes(Ms)

    X, labs, brns, scalings = sim.sample_data_balanced(1, G, t, sample_time, alpha, beta, scale_v=0.8)

    # job_id = "test"
    # save_dir = "/home/npapado/Desktop"
    save_params(job_id, save_dir, G, br_lengths, br_compl, rseed, top)
    save_files(job_id, save_dir, X, labs, brns, scalings, uMs, Hs, gene_scale, alpha, beta)
    scalefile = save_dir + "/" + job_id + "_scalings.txt"
    np.savetxt(scalefile, scalings)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate a simulated scRNAseq dataset with n bifurcation.')
    parser.add_argument("-j", "--job", dest="job",
                  help="Job ID (prepended to all generated files)", metavar="FILE")
    parser.add_argument("-o", "--out", dest="outdir",
                  help="Directory where output files are saved", metavar="FILE")
    parser.add_argument("-n", "--num_brpoints", dest="n",
                  help="How many branching points the simulation contains", metavar="FILE")

    args = parser.parse_args()

    main(args.job, args.outdir, int(args.n))
