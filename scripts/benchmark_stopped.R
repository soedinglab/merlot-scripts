#!/usr/bin/Rscript

args <- commandArgs(trailingOnly = TRUE)
mscripts <- args[1]
JobFolder <- args[2]
JobName <- paste(args[3], "_", sep="")
method <- args[4]

various <- paste(mscripts, "/scripts/various.R", sep="")
evaluat <- paste(mscripts, "/scripts/evaluate_method.R", sep="")

suppressPackageStartupMessages(source(various))
suppressPackageStartupMessages(source(evaluat))

res <- evaluate_method(method, NA, NA, NA, NA, returnNA=TRUE)
read_write_output(res, paste(JobFolder, JobName, "_eval.txt", sep=""))