#-----
# runs partner fidelity on one of the shuffled networks. 
# this is meant to be ran on the HPC
# input argument: experiment id number from the experiments.csv file
# input files: _edge_list.csv - file with intralayer edges
#              fitted_asvs_phylo_tree.rds - file with the phylogenetic tree of all the ASVs
# other requirements: functions.R
# all need to be in the working directory.
#------

#! /gpfs0/shai/projects/R4/R-4.0.3/bin/Rscript
.libPaths("/gpfs0/shai/projects/R4/R-4.0.3/lib64/R/library")
print(.libPaths())
print(sessionInfo())

library(tidyverse)
library(magrittr)
library(vegan)
library(dendextend)
library(GUniFrac)
source('../functions.R')
getwd()

# parse arguments
JOB_ID <- Sys.getenv("JOB_ID")
if (length(commandArgs(trailingOnly=TRUE))==0) {
  stop('No arguments were found!')
} else {
  args <- commandArgs(trailingOnly=TRUE)
  net_id <- as.numeric(args[1])
}

#net_id=1 
print(paste("Suffled network number:", net_id ))


# read the shuffled networks and analyze -----
# read shuffled network layers
print("Reading network files." )
farm_net_shuff <- list.files(path = "." , pattern = paste('_edge_list.csv',sep=""), recursive = T,full.names = T)
farm_net_shuff <- sapply(farm_net_shuff, read_csv, simplify=FALSE) %>% 
  bind_rows(.id = "id") 

# fidelity analysis on the shuffled network
print("Calculate fidelity for network." )

otherway <- farm_net_shuff %>% 
            relocate(to, from) %>%
            rename(to=from, from=to)
both <- bind_rows(farm_net_shuff, otherway) %>% 
      # keep only the ones that appear in 2 or more farms
      group_by(from) %>%
      mutate(num_farms_from=n_distinct(level_name)) %>%
      filter(num_farms_from>=2)
fidelity_shuff <- both %>% 
                  group_by(from) %>% 
                  group_modify(~calculate_PF_J(.x))
fidelity_shuff$run <- net_id
write_csv(fidelity_shuff, 'fidelity_shuff_farm_30.csv')


print("Calculate UniFrec fidelity for network." )


# read tree from phylo data
phylo_tree <- readRDS("fitted_asvs_phylo_tree.rds")

unifreq_shuff <- both %>% 
                 group_by(from) %>% 
                 group_modify(~calculate_PF_U(.x, phylo_tree$tree))
unifreq_shuff$run <- net_id
write_csv(unifreq_shuff, 'uniFrec_shuff_farm_30.csv')

print("Done analysing the network." )

