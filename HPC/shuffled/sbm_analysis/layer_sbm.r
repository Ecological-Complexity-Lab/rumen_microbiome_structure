#! /gpfs0/shai/projects/R4/R-4.0.3/bin/Rscript
# ------ layer_sbm.r ----------------------------------------
# This script calculated the probability based Transitivity in a single layer
# they will later be moved to more logically relevant places in the repo
#-----------------------------------------------------------------------
.libPaths("/gpfs0/shai/projects/R4/R-4.0.3/lib64/R/library")
print(.libPaths())
print(sessionInfo())


#------ includes ----------
library(tidyverse)
library(magrittr)
library(reshape2)
library(igraph)
library(blockmodels)

#------ Arguments ---------
JOB_ID <- Sys.getenv("JOB_ID")
if (length(commandArgs(trailingOnly=TRUE))==0) {
  stop('No arguments were found!')
} else {
  args <- commandArgs(trailingOnly=TRUE)
  input_folder <- as.character(args[1])
}

#----- Run -------
# Read current shuffled data
# add column for counting shuffled networks with empiric link
files <- list.files(path = input_folder , pattern = '_edge_list.csv', recursive = T,full.names = T)

gps <- NULL
for (l in files) {
  print(l)
  l_net <- read_csv(l) %>% filter(edge_type=="pos")
  farm <- pull(l_net[1, "level_name"])
  id <- str_split(basename(l), "_")[[1]][1]

  g <- graph.data.frame(l_net[,1:3])
  adj <- get.adjacency(g,sparse=FALSE, attr='weight')
  
  sbm_model <- BM_bernoulli$new("SBM", adj)
  sbm_model$estimate()
  max_group <- which.max(sbm_model$ICL)
  
  line <- tibble(id=id, farm=farm, shuf_max_ICL=max_group)
  gps <- rbind(gps, line)
}

# save sbm results
write_csv(gps, "shuff_layer_SBM_results.csv", append = T)


