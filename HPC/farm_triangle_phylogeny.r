#----- farm_triangle_phylogeny.r ----------
# this script saves the pylogenetic distances of every 
# triangle in a specific farm
#------

#! /gpfs0/shai/projects/R4/R-4.0.3/bin/Rscript
.libPaths("/gpfs0/shai/projects/R4/R-4.0.3/lib64/R/library")
print(.libPaths())
print(sessionInfo())

library(tidyverse)
library(magrittr)
library(vegan)
library(dendextend)
library(ape)
library(igraph)
source('../functions.R')
getwd()

# --- functions ------
# the function returns the longest distance in the triangle
max_tri_dist <- function(x, distances) {
  asvs <- pull(x[,1])
  prs <- combn(asvs, 2) # All pair combinations
  
  dists <- c()
  for (pair in 1:ncol(prs)) {
    dists[pair] <- PatristicDistMatrix[prs[1,pair], prs[2,pair]]
  }
  out <- data.frame(max_dist=max(dists), 
                    mean_dist=mean(dists), 
                    min_dist=min(dists),
                    sum_dis=sum(dists),
                    asvs=paste(asvs[1], asvs[2], asvs[3], sep = " + "))
  return(out)
}

# --- argument handling ------

# args <- c('UK2')
JOB_ID <- Sys.getenv("JOB_ID")
if (length(commandArgs(trailingOnly=TRUE))==0) {
  stop('No arguments were found!')
} else {
  args <- commandArgs(trailingOnly=TRUE)
  curr_farm <- as.character(args[1])
}


# --- run the analysis -----

# filter file name with "Farm_"+ curr_farm +"_edge_list.csv" in the end
farm_net_file <- list.files(path = "." , 
                            pattern = paste("Farm", curr_farm, "edge_list.csv",sep="_"), 
                            recursive = T,full.names = T)


# read it and filter only the positive co-occarances
farm_multilayer_pos_30 <- read_csv(farm_net_file) %>% filter(edge_type=="pos")
 
# create the triangles
g <- graph_from_data_frame(farm_multilayer_pos_30, directed = FALSE, vertices = NULL)
tri <- names(triangles(g))
tri <- tibble(tri_id=rep(1:(length(tri)/3),each=3), ASV_ID=tri)


# read distance matrix made from the phylogenetic tree
PatristicDistMatrix <- readRDS("../distance_mat_from_tree.rds")


# calculate per triangle the max distance, (also the min and sum mean) 
res <- tri %>% 
  filter(tri_id<=500) %>% # TODO this is temporary to work with a small amount first.
  group_by(tri_id) %>% 
  group_modify(~max_tri_dist(.x, PatristicDistMatrix))

write_csv(res, paste('triangles_phylo_distances_',curr_farm, '.csv', sep = ""))


