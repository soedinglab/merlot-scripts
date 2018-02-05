#!/usr/bin/Rscript

args <- commandArgs(trailingOnly = TRUE)
hhtree <- args[1]
JobFolder <- args[2]
JobName <- paste(args[3], "_", sep="")
method <- args[4]
# hhtree <- "~/Documents/repos/hhtree"
# JobFolder <- "~/Desktop/benchmark/test0/"
# JobName <- paste("test0", "_", sep="")
# dimensions <- 2

various <- paste(hhtree, "/scripts/various.R", sep="")
evaluat <- paste(hhtree, "/scripts/evaluate_method.R", sep="")

suppressPackageStartupMessages(source(various))
suppressPackageStartupMessages(source(evaluat))

res <- evaluate_method(method, NA, NA, NA, NA, returnNA=TRUE)
read_write_output(res, paste(JobFolder, JobName, "_eval.txt", sep=""))