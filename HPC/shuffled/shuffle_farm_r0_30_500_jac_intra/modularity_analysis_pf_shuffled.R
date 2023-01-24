#-----
# runs infomap on a multilayer network with partner fidelity based interlayer edges.
# input argumemnst: experiment id number from the experiments.csv file
# input files: _edge_list.csv file with intralayer edges
# other requirements: working infomap exe (make sure it is the ubuntu version, and has running permission)
#                     functions.R
# all need to be in the working directory.
#------


#! /gpfs0/shai/projects/R4/R-4.0.3/bin/Rscript
.libPaths("/gpfs0/shai/projects/R4/R-4.0.3/lib64/R/library")
print(.libPaths())
print(sessionInfo())

library(tidyverse)
library(magrittr)
library(reshape2) 
library(igraph)
library(vegan)
#library(infomapecology) # to be uncommented once R version is updated on the HPC
library(extRC)
library(GUniFrac)

# parse arguments
JOB_ID <- Sys.getenv("JOB_ID")
if (length(commandArgs(trailingOnly=TRUE))==0) {
  stop('No arguments were found!')
} else {
  args <- commandArgs(trailingOnly=TRUE)
  net_id <- as.numeric(args[1])
}

#net_id=10

print(getwd())
source('functions.R')
source('../temp_infomapecology.r') # to be removed once R version is updated on the HPC
check_infomap()

# reading co-occurrance intralayer edges ----
one_shuff_farm_f <- list.files(pattern = paste("^",net_id,"_",sep=''), full.names = T)
print(one_shuff_farm_f)
one_shuff_farm_f <- str_subset(one_shuff_farm_f, "edge_list.csv$")
print(one_shuff_farm_f)
#bind all farms to a single multilayer
farm_multilayer <- sapply(one_shuff_farm_f, read_csv, simplify=FALSE) %>% 
  bind_rows(.id = "id") %>% select(-id)

# For the analysis separate the positive and negative networks
farm_multilayer_pos <- farm_multilayer %>% filter(edge_type=='pos')
farm_multilayer_neg <- farm_multilayer %>% filter(edge_type=='neg')

all_nodes <- sort(unique(c(farm_multilayer_pos$from, farm_multilayer_pos$to)))
all_nodes <- tibble(node_id=1:length(all_nodes), node_name=all_nodes)

layers <- tibble(layer_id=1:7, layer_name=unique(farm_multilayer_pos$level_name))

# building partner fidelity interlayer edges ----
# Define interlayer edges as Jaccard similarity in neighbors between farms----

# bind 'from' and 'to'
farm_multilayer_pos_ASV2 <-
  farm_multilayer_pos %>%
  relocate(to, from) %>%
  rename(to=from, from=to)
farm_multilayer_pos_final <- bind_rows(farm_multilayer_pos, farm_multilayer_pos_ASV2)

# keep only the ones that appear in 2 or more farms
farm_multilayer_pos_final %<>%
  group_by(from) %>%
  mutate(num_farms_from=n_distinct(level_name)) %>%
  filter(num_farms_from>=2)

# calculate intelayer edges
fidelity_inter_jaccard <- NULL
for (i in unique(farm_multilayer_pos_final$from)) {
  print(i)
  partners_mat_ASV <- 
    farm_multilayer_pos_final %>%
    filter(from==i) %>%
    group_by(to) %>%
    select(c(to,level_name)) %>%
    mutate(present=1) %>%
    spread(to, present, fill = 0) %>%
    column_to_rownames("level_name")
  beta_farms_ASV <- 1-as.matrix(vegdist(partners_mat_ASV, "jaccard"))
  beta_farms_ASV_m <- melt(as.matrix(extRC::tril(beta_farms_ASV)))
  inter_fid <- beta_farms_ASV_m %>%
    tibble() %>%
    filter(value!=0) %>%
    subset(Var1 != Var2) %>%
    mutate(ASV_ID=i) %>%
    select(c(ASV_ID,layer_from=Var1, layer_to=Var2, weight=value))
  fidelity_inter_jaccard <- rbind(fidelity_inter_jaccard,inter_fid)
}

# Remove nodes from the interlayer edges that do not have intralayer edges in
# any layer and are therefore not in the all_nodes tibble
nodes_to_remove_j <- setdiff(unique(fidelity_inter_jaccard$ASV_ID),unique(all_nodes$node_name))
fidelity_inter_jaccard %<>% filter(!ASV_ID%in%nodes_to_remove_j)

# Connect the intra and interlayer edges and change names to IDs
intra_j <- 
  farm_multilayer_pos %>% 
  select(layer_from=level_name, node_from=from, layer_to=level_name, node_to=to, weight) %>% 
  left_join(layers, by = c('layer_from' = 'layer_name')) %>% 
  left_join(layers, by = c('layer_to' = 'layer_name')) %>% 
  left_join(all_nodes, by = c('node_from' = 'node_name')) %>% 
  left_join(all_nodes, by = c('node_to' = 'node_name')) %>% 
  select(layer_from=layer_id.x, node_from=node_id.x, layer_to=layer_id.y, node_to=node_id.y, weight)

# mold the df to the right format
inter_j <- 
  fidelity_inter_jaccard %>% 
  ungroup() %>% 
  left_join(layers, by = c('layer_from' = 'layer_name')) %>% 
  left_join(layers, by = c('layer_to' = 'layer_name')) %>% 
  left_join(all_nodes, by = c('ASV_ID' = 'node_name')) %>% 
  select(layer_from=layer_id.x, node_from=node_id, layer_to=layer_id.y, node_to=node_id, weight)

multilayer_jaccard <- rbind(intra_j %>% mutate(type='intra'),
                    inter_j %>% mutate(type='inter'))
table(multilayer_jaccard$type)

# Run Infomap ----
x_j <- create_multilayer_object(extended = multilayer_jaccard[,1:5], nodes = all_nodes, layers = layers)
m_j <- run_infomap_multilayer(x_j, flow_model = 'undirected', silent = F, trials = 100, relax = F)
farm_modules_pos_j <- m_j$modules %>% left_join(layers)

# In the file name include net_id, multilayer
write_csv(multilayer_jaccard, paste(net_id,'_multilayer_pf_jaccard.csv',sep=""))

# In the file name include net_id, farm
write_csv(farm_modules_pos_j, paste(net_id,'_farm_modules_pf_jaccard.csv',sep=""))

# Run Infomap multilevel----
m_j_ls <- run_infomap_multilayer_multilevel(x_j, two_level = F, 
                                            flow_model = 'undirected', silent = F, 
                                            trials = 100, relax = F, seed=NULL)
modules <- m_j_ls$modules %>% left_join(layers)
# See that the file name matches the network before writing
write_csv(modules, paste(net_id,'_farm_modules_pos_30_J_multilevel.csv',sep=""))

# Write summary when done ------------------
m_file_j <- tibble(e_id=net_id, JOB_ID=JOB_ID, call=c(m_j$call, m_j_ls$call), 
                                               L=c(m_j$L, m_j_ls$L), 
                                               modules_num=c(m_j$m, m_j_ls$m))

summary_file_j <- '../farm_modulation_summary_pf_jaccard.csv'
write_csv(m_file_j, summary_file_j, append = T)

# Define interlayer edges with unifrac----
phylo_tree <- readRDS("fitted_asvs_phylo_tree.rds")
tree <- phylo_tree$tree

# calculate intelayer edges
fidelity_unif_inter <- NULL
for (i in unique(farm_multilayer_pos_final$from)) {
  print(i)
  ASV_net <- farm_multilayer_pos_final %>%
    filter(from==i)
  
  # prune the tree
  included_asvs <- unique(ASV_net$to)
  unincluded <- tree$tip.label[!tree$tip.label %in% included_asvs]
  pruned <- dendextend::prune(tree, unincluded)
  
  mat_farm_ASV <- 
    farm_multilayer_pos_final %>%
    filter(from==i) %>%
    group_by(to) %>%
    select(c(to,level_name)) %>%
    mutate(present=1) %>%
    spread(to, present, fill = 0) %>%
    column_to_rownames("level_name")
  # run unifrec
  unifracs <- GUniFrac(mat_farm_ASV, pruned, alpha=c(0, 0.5, 1))$unifracs
  d_UW_ASV_mat <- 1-(unifracs[, , "d_UW"])
  d_UW_ASV_mat_m <- melt(as.matrix(extRC::tril(d_UW_ASV_mat)))
  inter_fid_unif <- d_UW_ASV_mat_m %>%
    tibble() %>%
    filter(value!=0) %>%
    subset(Var1 != Var2) %>%
    mutate(ASV_ID=i) %>%
    select(c(ASV_ID,layer_from=Var1, layer_to=Var2, weight=value))
  fidelity_unif_inter <- rbind(fidelity_unif_inter,inter_fid_unif)
}

# Remove nodes from the interlayer edges that do not have intralayer edges in
# any layer and are therefore not in the all_nodes tibble
nodes_to_remove_u <- setdiff(unique(fidelity_unif_inter$ASV_ID),unique(all_nodes$node_name))
fidelity_unif_inter %<>% filter(!ASV_ID%in%nodes_to_remove_u)

# Connect the intra and interlayer edges and change names to IDs
intra_u <- 
  farm_multilayer_pos %>% 
  select(layer_from=level_name, node_from=from, layer_to=level_name, node_to=to, weight) %>% 
  left_join(layers, by = c('layer_from' = 'layer_name')) %>% 
  left_join(layers, by = c('layer_to' = 'layer_name')) %>% 
  left_join(all_nodes, by = c('node_from' = 'node_name')) %>% 
  left_join(all_nodes, by = c('node_to' = 'node_name')) %>% 
  select(layer_from=layer_id.x, node_from=node_id.x, layer_to=layer_id.y, node_to=node_id.y, weight)

# mold the df to the right format
inter_u <- 
  fidelity_unif_inter %>% 
  ungroup() %>% 
  left_join(layers, by = c('layer_from' = 'layer_name')) %>% 
  left_join(layers, by = c('layer_to' = 'layer_name')) %>% 
  left_join(all_nodes, by = c('ASV_ID' = 'node_name')) %>% 
  select(layer_from=layer_id.x, node_from=node_id, layer_to=layer_id.y, node_to=node_id, weight)

multilayer_unif <- rbind(intra_u %>% mutate(type='intra'),
                         inter_u %>% mutate(type='inter'))
table(multilayer_unif$type)

# Run Infomap ----
x_u <- create_multilayer_object(extended = multilayer_unif[,1:5], nodes = all_nodes, layers = layers)
m_u <- run_infomap_multilayer(x_u, flow_model = 'undirected', silent = F, trials = 100, relax = F)
farm_modules_pos_unif <- m_u$modules %>% left_join(layers)

# In the file name include net_id, multilayer
write_csv(multilayer_unif, paste(net_id,'_multilayer_pf_unif.csv',sep=""))

# In the file name include net_id, farm
write_csv(farm_modules_pos_unif, paste(net_id,'_farm_modules_pf_unif.csv',sep=""))


# Run Infomap multilevel----
m_u_ls <- run_infomap_multilayer_multilevel(x_u, two_level = F, 
                                            flow_model = 'undirected', silent = F, 
                                            trials = 100, relax = F, seed=NULL)
modules <- m_u_ls$modules %>% left_join(layers)
# See that the file name matches the network before writing
write_csv(modules, paste(net_id,'_farm_modules_pos_30_U_multilevel.csv',sep=""))


# Write summary when done ------------------
m_file_u <- tibble(e_id=net_id, JOB_ID=JOB_ID, 
                   call=c(m_u$call, m_u_ls$call), 
                   L=c(m_u$L, m_u_ls$L), 
                   modules_num=c(m_u$m, m_u_ls$m), 
                   time_stamp=Sys.time())

summary_file_u <- '../farm_modulation_summary_pf_unif.csv'
write_csv(m_file_u, summary_file_u, append = T)
