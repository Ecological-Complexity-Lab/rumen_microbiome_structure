#! /gpfs0/shai/projects/R4/R-4.0.3/bin/Rscript
# ------ validate_link.r ----------------------------------------
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

#------ Arguments ---------
JOB_ID <- Sys.getenv("JOB_ID")
if (length(commandArgs(trailingOnly=TRUE))==0) {
  stop('No arguments were found!')
} else {
  args <- commandArgs(trailingOnly=TRUE)
  l <- as.character(args[1])
}

#----- Run -------
net_shuffled <- read_csv("all_shuff_networks_r0_30_500_jac_intra.csv")

# Read ampiric data
# add column for counting shuffled networks with empiric link
with_comp <- read_csv('farm_multilayer_pos_30.csv') %>% 
              select(layer=level_name, node_from=from, node_to=to, weight) %>% 
              add_column(count=0)

head(net_shuffled)
nets <- net_shuffled %>% filter(farm == l)
with_comp %<>% filter(layer == l)

# count links that appear in the empiric network
for (r in 1:nrow(with_comp)) {
  if(r %% 200 == 0) {
    print(100*r/nrow(with_comp))
  }
  row1 <- nets %>% filter(farm == pull(with_comp[r, "layer"]),
                          from == pull(with_comp[r, "node_from"]),
                          to == pull(with_comp[r, "node_to"]))
  # need to check both because its not directed and random
  row2 <- nets %>% filter(farm == pull(with_comp[r, "layer"]),
                          from == pull(with_comp[r, "node_to"]),
                          to == pull(with_comp[r, "node_from"]))
  # record number of times they appear in shuffled network
  with_comp[r, "count"] <- nrow(row1) + nrow(row2) # number of times the link appeared across the networks
}
print(100*r/nrow(with_comp))

# calculate p value for each link
with_comp %<>% mutate(p_val=count/500)

# Save link validations
write_csv(with_comp, paste(l,'link_validation_r0_50_500.csv', sep = "_"))

