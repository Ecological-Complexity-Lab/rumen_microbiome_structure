#! /gpfs0/shai/projects/R4/R-4.2.0/bin/Rscript
# ------ test_interlayer_amount.r ----------------------------------------
# This script calculated the probability based Transitivity in a single layer
# they will later be moved to more logically relevant places in the repo
#-----------------------------------------------------------------------
.libPaths("/gpfs0/shai/projects/R4/R-4.2.0/lib64/R/library")
print(.libPaths())
print(sessionInfo())


#------ includes ----------
library(tidyverse)
library(magrittr)
library(infomapecology)

# handle 
JOB_ID <- Sys.getenv("JOB_ID")
if (length(commandArgs(trailingOnly=TRUE))==0) {
  stop('No arguments were found!')
} else {
  args <- commandArgs(trailingOnly=TRUE)
  itr_id <- as.numeric(args[1])
}

print(paste("This iteration:", itr_id, ". job ID:", JOB_ID, "."))


# Run -----
set.seed(itr_id)

mln <- read_csv("../mln_multi_interlayer_edges_30.csv")
layers <- tibble(layer_id=1:7, 
                 short_name=c('UK1', 'UK2', 'IT1', 'IT2', 'IT3', 'FI1', 'SE1'))

intra <- mln %>% filter(layer_from==layer_to)
inter <- mln %>% filter(layer_from!=layer_to)

mods <- NULL
modules_details <- NULL
for (p in seq(from = 0, to = 1, by = 0.1)) {
  n_new_rows <- p*nrow(inter)
  print(paste("Running percentage:", p*100, "%. edge amount:", n_new_rows, "."))
  
  if (n_new_rows == 0) {
    n_new_rows <- 1
  }
  
  sampled <- inter %>% sample_n(n_new_rows)
  new_nlm <- rbind(intra, sampled)
  
  # make all nodes etc
  all_nodes <- sort(unique(c(new_nlm$node_from, new_nlm$node_to)))
  all_nodes <- tibble(node_id=1:length(all_nodes), node_name=all_nodes)
  
  # make network be represented by indexes
  mln_inds <- 
    new_nlm %>%
    left_join(layers, by = c('layer_from' = 'short_name')) %>% 
    left_join(layers, by = c('layer_to' = 'short_name')) %>% 
    left_join(all_nodes, by = c('node_from' = 'node_name')) %>% 
    left_join(all_nodes, by = c('node_to' = 'node_name')) %>% 
    select(layer_from=layer_id.x, node_from=node_id.x, layer_to=layer_id.y, node_to=node_id.y, weight)
  
  # run the iteration of the infomap
  inf_obj <- create_multilayer_object(extended = mln_inds,
                                      nodes = all_nodes,
                                      layers = layers)
  m <- infomapecology::run_infomap_multilayer(inf_obj, silent = T,
                                              flow_model = 'undirected',  
                                              trials = 200, relax = F, seed=111)
  
  # save results
  lne <- tibble(percent=p, modules=m$m, itr=itr_id)
  mods <- rbind(mods, lne)
  
  to_add <- m$modules %>% add_column(itr = itr_id, percent=p)
  modules_details <- rbind(modules_details, to_add)
}

#save results
write_csv(mods, '../sampled_multi_edge_30_n_modules.csv', append = T)
write_csv(modules_details, '../sampled_multi_edge_30_modules.csv', append = T)
