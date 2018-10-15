#!/usr/bin/Rscript

suppressPackageStartupMessages(library("optparse"))

parse_scaffolds <- function(subfolders, input_dir, scaf_names) {
  branches <- matrix(0, ncol = length(scaf_names), nrow = length(subfolders))
  bpoints <- matrix(0, ncol = length(scaf_names), nrow = length(subfolders))
  epoints <- matrix(0, ncol = length(scaf_names), nrow = length(subfolders))
  rownames(branches) <- subfolders
  rownames(bpoints) <- subfolders
  rownames(epoints) <- subfolders
  colnames(branches) <- scaf_names
  colnames(bpoints) <- scaf_names
  colnames(epoints) <- scaf_names
  
  pb <- txtProgressBar(min = 5, max = 100, style = 3)
  for (i in 1:length(subfolders)) {
    s = subfolders[i]
    setTxtProgressBar(pb, i)
    for (c in scaf_names) {
      to_read <- paste(input_folder, "/", s, "/", s, "_", c, ".scaf", sep="")
      scaf <- readRDS(to_read)
      if (all(scaf$Branchpoints == 0)) {
        bpoints[s, c] <- 0
        branches[s, c] <- 1
        epoints[s, c] <- 2
      } else {
        bpoints[s, c] <- length(scaf$Branchpoints)
        epoints[s, c] <- length(scaf$Endpoints)
        branches[s, c] <- dim(scaf$Branches)[1]
      }
    }
  }
  res <- list(bpoints, branches, epoints)
  return (res)
}

option_list = list(
  make_option(c("-o", "--out"), type = "character", default = "",
              help = "Location where the output is stored", metavar = "/output/folder"),
  make_option(c("-b", "--benchmark"), type = "character", 
              help = "Location of the benchmark.", metavar = "/benchmark/folder")
);


opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

output_folder <- opt$out
input_folder <- opt$benchmark

if (output_folder == "") {
  output_folder <- input_folder
}

subfolders <- dir(input_folder)
remove <- c()

scaffolds <- c("hhTree_log_k_fixed_el",
               "hhTree_log_k_free_el_sens",
               "LPGraph_log_k_fixed_el",
               "LPGraph_log_k_free_el_sens",
               "monGraph_fixed_el",
               "monGraph_free_el_sens",
               "monGraph_unc_fixed_el",
               "monGraph_unc_free_el_sens",
               "monTree_fixed_el",
               "monTree_free_el_sens",
               "monTree_unc_fixed_el",
               "monTree_unc_free_el_sens")

# remove missing files
for (s in subfolders) {
  scaf1 <- paste(input_folder, "/", s, "/", s, "_", scaffolds[1], ".scaf", sep="")
  scaf2 <- paste(input_folder, "/", s, "/", s, "_", scaffolds[2], ".scaf", sep="")
  scaf3 <- paste(input_folder, "/", s, "/", s, "_", scaffolds[3], ".scaf", sep="")
  scaf4 <- paste(input_folder, "/", s, "/", s, "_", scaffolds[4], ".scaf", sep="")
  
  scaf_locs <- c(file.exists(scaf1), file.exists(scaf2), file.exists(scaf3), file.exists(scaf4))
  if (!all(scaf_locs)) {
    remove <- cbind(remove, s)
  }
}

subfolders <- subfolders[!is.element(subfolders, remove)]
print(paste(length(subfolders), "/100", sep=""))

gonzalo_info <- parse_scaffolds(subfolders, input_dir, scaffolds)
bpoints <- gonzalo_info[[1]]
branches <- gonzalo_info[[2]]
epoints <- gonzalo_info[[3]]

location = strsplit(input_folder, split = "/")[[1]]

write.csv(bpoints, paste(output_folder, "/", location[length(location)], "_branchpoints.csv", sep = ""))
write.csv(branches, paste(output_folder, "/", location[length(location)], "_branches.csv", sep = ""))
