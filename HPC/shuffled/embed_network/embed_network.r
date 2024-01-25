#! /gpfs0/shai/projects/R4/R-4.0.3/bin/Rscript
# ------ embed_network.r ----------------------------------------
# -----------------------------------------------------------
.libPaths("/gpfs0/shai/projects/R4/R-4.0.3/lib64/R/library")
print(.libPaths())
print(sessionInfo())

#------ includes ----------
library(tidyverse)
library(magrittr)
library(reshape2)
library(data.table)
library(node2vec)
library(Rtsne)

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

teb_all <- NULL
for (l in files[1:2]) {
  print(l)
  l_net <- read_csv(l) %>% filter(edge_type=="pos")
  farm <- pull(l_net[1, "level_name"])
  id <- str_split(basename(l), "_")[[1]][1]
  
  node2vec_model <- node2vecR(as.data.frame(l_net[,1:3]), walk_length = 10)
  
  # Perform t-SNE for 2D visualization
  tsne_res <- Rtsne(node2vec_model, dims = 2)$Y
  
  teb <- data.frame(tsne_res) %>% select(x=1, y=2) %>% add_column(farm=l)
  teb_all <- rbind(teb_all, teb)
}

# save sbm results
write_csv(teb_all, paste("results/", id, "_shuff_net_embeding.csv", sep = ""))

avrg <- teb_all %>% group_by(farm) %>% summarise(xx=mean(x), yy=mean(y)) %>%
          add_column(run=id)
med <- teb_all %>% group_by(farm) %>% summarise(xx=median(x), yy=median(y)) %>%
  add_column(run=id)

write_csv(avrg, "mean_shuff_net_embeding.csv", append = T)
write_csv(med, "median_shuff_net_embeding.csv", append = T)

