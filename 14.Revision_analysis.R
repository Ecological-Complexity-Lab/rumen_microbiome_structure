# ------ 14.Revision_analysis.r ----------------------------------------
# This script temporarily holds the analysis requested in the revision
# they will later be moved to more logically relevant places in the repo
#-----------------------------------------------------------------------

#------ includes ----------
library(tidyverse)
library(magrittr)
library(reshape2)
library(data.table)
library(igraph)
library(blockmodels)
library(vegan)

#------ consts ---------
single_prob_file <- 'local_output/single_asv_occur_prob_80.csv'
combo_prob_file <- 'combo_asv_occur_prob.csv'

#------ run ------------
## Read data to be used -----
# read network:
farm_multilayer_pos_30 <- read_csv('local_output/farm_multilayer_pos_30.csv')
all_nodes <- sort(unique(c(farm_multilayer_pos_30$from, farm_multilayer_pos_30$to)))
all_nodes <- tibble(node_id=1:length(all_nodes), node_name=all_nodes)
layers <- tibble(layer_id=1:7, layer_name=c('NUDC', 'Park', 'Bian', 'Fran','Gand','Mink','Raab'),
                 short_name=c('UK1', 'UK2', 'IT1', 'IT2', 'IT3', 'FI1', 'SE1'))
intras <- farm_multilayer_pos_30 %>% 
            select(layer=level_name, node_from=from, node_to=to, weight)


# read asv data
ASV_data_80 <- read_csv("local_output/core_ASV_80.csv")
ASV_data_30 <- read_csv("local_output/core_ASV_30.csv")


# ------ Cow level: ---------
## link vs shuffled -----
# Compare link from empiric to all in 500 shuffled, 
# to see if they exist in less then 5% of 500 shuffled networks so they are significant

# read shuffled networks
parent.folder <- "HPC/shuffled/shuffle_farm_r0_30_500_jac_intra"
files <- list.files(path = parent.folder , pattern = '_multilayer_pf_unif.csv', recursive = T,full.names = T)

# read all shuffled networks, and saving as full names
net_shuffled <- NULL
for (s in files) {
  print(s)
  folder <- dirname(s)
  node_file <- list.files(path = folder , pattern = '_farm_modules_pf_unif.csv', recursive = T,full.names = T)
  nodes <- fread(node_file) %>% select(node_id, node_name) %>% distinct()
  lyrs <- fread(node_file) %>% select(layer_id, layer_name) %>% distinct()
  
  shuf_net <- fread(s) %>% filter(layer_from == layer_to) %>% select(layer = layer_from, node_from, node_to)
  shuf_net$id <- str_split_fixed(s, pattern = '/', n = 5)[4]
  
  # add asv ids
  shuf_net %<>% left_join(lyrs, by=c("layer"="layer_id")) %>%
    left_join(nodes, by=c("node_from"="node_id")) %>% rename(from=node_name) %>%
    left_join(nodes, by=c("node_to"="node_id"))   %>% rename(to=node_name) %>%
    select(id, farm=layer_name, from, to)
  
  net_shuffled <- rbind(net_shuffled, shuf_net)
}

net_shuffled <- as_tibble(net_shuffled)
# this is used in the HPC analysis
write_csv(net_shuffled, 'local_output/all_shuff_networks_r0_30_500_jac_intra.csv')

# --- HPC side step 

# read validation results produced on the HPC
all_ <- NULL
for (l in layers$short_name) {
  s <- paste('HPC/validate_link/', l,'_link_validation_r0_50_500.csv', sep = "")
  res <- fread(s)
  
  all_ <- rbind(all_, res)
  
  print(sum(res$p_val <0.05))
}

# plot the layers
all_ %>% ggplot(aes(x=p_val, fill = as.factor(layer)))+
  geom_histogram(aes(y = after_stat(density)), alpha=0.4, position='identity')


# ------ Farm level: ---------
## SBM on layer -----
# find group number per layer in empiric network
gps <- NULL
mems_table <- NULL
for (l in layers$short_name) {
  print(l)
  lay <- intras %>% filter(layer == l)
  
  g <- graph.data.frame(lay[,2:4])
  adj <- get.adjacency(g,sparse=FALSE, attr='weight')
  
  sbm_model <- BM_bernoulli$new("SBM", adj)
  sbm_model$estimate()
  max_group <- which.max(sbm_model$ICL)
  mem <- sbm_model$memberships[[max_group]]$Z
  
  nds <- row.names(adj)
  row.names(mem) <- nds

  # membership documenting - find grouping
  grp_mem <- apply(mem, 1, function(x) match(max(x), x)) # this will take the first group with the max value as the node's group.
  mems_tbl <- tibble(farm=l, asv_id=nds, membership=grp_mem)

  gps <- c(gps, max_group)
  mems_table <- rbind(mems_table, mems_tbl)
}
groups <- layers %>% select(short_name) %>% add_column(emp_max_ICL=gps)

write_csv(groups, "local_output/layer_SBM_results.csv")
write_csv(mems_table, "local_output/layer_SBM_membership_results.csv")

groups <- read_csv("local_output/layer_SBM_results.csv")
mems_table <- read_csv("local_output/layer_SBM_membership_results.csv")


# --- shuffled data from the HPC

# process results from HPC
res_file <- read_csv(gps, "HPC/shuffled/sbm_analysis/shuff_layer_SBM_results.csv")

# read empiric data
emp <- read_csv("local_output/layer_SBM_results.csv") %>% add_column(id="000") %>%
              select(id, farm=short_name, max_ICL=emp_max_ICL)

all_groups <- rbind(emp, res_file %>% rename(max_ICL=shuf_max_ICL))

all_groups %>% ggplot(aes(x=max_ICL, fill = as.factor(farm))) +
  geom_histogram(aes(y = after_stat(density)), alpha=0.4, position='identity')


## network embedding -----
# Process results from HPC run:
# read embedding data that was ran on the HPC
mean_embd <- read_csv("HPC/shuffled/embed_network/mean_shuff_net_embeding.csv")
med_embd <- read_csv("HPC/shuffled/embed_network/median_shuff_net_embeding.csv")

# Plot mean embedding data
ggplot(mean_embd, aes(mean_embd, yy, color = farm)) +
  geom_point() +
  theme_minimal() +
  ggtitle("Mean Graph Embeddings")

# Plot mean embedding data
ggplot(med_embd, aes(xx, yy, color = farm)) +
  geom_point() +
  theme_minimal() +
  ggtitle("Median Graph Embeddings")



# ------ Interfarm level: ---------

## taxonomic beta-diversity ------ 
# or phylogeny partner fidelity

# read phylogenetic data
ASV_taxa <- read_csv('local_output/ASV_full_taxa.csv') %>% 
  select(ASV_ID, everything(), -seq16S)

# filter only taxa that exist in the networks
asvs <- sort(unique(c(intras$node_from, intras$node_to)))
all_taxa <- ASV_taxa %>% filter(ASV_ID %in% asvs)

get_taxa_pf <- function(taxa_level="Order") {
  print(sum(is.na(all_taxa[,taxa_level]))) #FYI
  
  # Handle order taxa 
  ord <- intras %>% 
    left_join(all_taxa, by = c('node_from' = 'ASV_ID')) %>%
    select(layer, node_from, node_to, weight, taxa_from = all_of(taxa_level)) %>%
    left_join(all_taxa, by = c('node_to' = 'ASV_ID')) %>%
    select(layer, node_from, node_to, weight, taxa_from, taxa_to = all_of(taxa_level))
  
  # get the number of connections between each taxa pair in a layer
  taxa_pairs <- ord %>% group_by(layer, taxa_from, taxa_to) %>% 
    summarise(count = n()) %>% drop_na() # we remove lined with unknown taxas
  otherway <- taxa_pairs %>% 
    rename(taxa_to=taxa_from, taxa_from=taxa_to)
  both <- bind_rows(taxa_pairs, otherway)
  
  # make it unique
  both %<>% group_by(layer, taxa_from, taxa_to) %>% 
    summarise(count = sum(count))
  
  # keep only taxas that occur in 2 or more farms
  both %<>%
    group_by(taxa_from) %>%
    mutate(num_farms_from=n_distinct(layer)) %>%
    filter(num_farms_from>=2)
  
  ## PF_T observed network:
  PF_T_obs <-
    both %>%
    group_by(taxa_from) %>%
    group_modify(~calculate_PF_T(.x)) %>% as_tibble()
  
  return(PF_T_obs)
}


hist(get_taxa_pf("Class")$PF_T)
hist(get_taxa_pf()$PF_T)
hist(get_taxa_pf("Family")$PF_T)
hist(get_taxa_pf("Genus")$PF_T)


## NMI between farms ----
# TODO implement - between sbm results


## modularity - phylogenetic composition --------
# includes randomality check
# TODO implement


## NMI of clusters and hypothesis ---------
# TODO implement



# --- not sure will be done ----

## Transitivity ---------

### single probabilities ----
# calculate single probability of each ASV for every layer
probabilities <- NULL
for (layer in layers$short_name) {
  # filter out the layer
  curr_layer <- ASV_data_80 %>% filter(Farm == layer)
  cows <- as.numeric(curr_layer %>% summarise(cows=n_distinct(Cow_Code)))
  
  print(cows)

  for (asv in all_nodes$node_name) {
    cows_with_asv <- as.numeric(curr_layer %>% filter(ASV_ID == asv)  %>% 
                          summarise(cows=n_distinct(Cow_Code)))
  
    probabilities <- rbind(tibble(ASV=asv, Farm=layer, rand_prob=cows_with_asv/cows),
                         probabilities)
  }
}

# save probabilities
write_csv(probabilities, single_prob_file)

### combo probabilities ----
# The combo properties is run on the HPC, here the results are processed
HPC_output_folder <- 'HPC/Transitivity 80/'

# for every layer, get the the phylogenetic groups in it, for every level
for (layer in layers$short_name) {
  # filter out the layer
  layer_asv <- ASV_data_80 %>% filter(Farm == layer)
  layer_net <- intras %>% filter(layer == layer)
  
  # run over results - 
  # TODO some analysis should be done here, this hist is not good enough - talk to shai
  rrr <- read_csv(paste(HPC_output_folder, layer, '_', combo_prob_file, sep = ""), 
                  col_names = c("ASV_1", "ASV_2", "ASV_3", "common_cows", "random", "emp_prob"))
  ns <- nrow(rrr)
  p1 <- hist(rrr$emp_prob)                    
  p2 <- hist(rrr$random)                     
  plot( p2, col=rgb(1,0,0,1/4), xlim=c(0,1), ylim = c(0,ns))
  plot( p1, col=rgb(0,0,1,1/4), xlim=c(0,1), add=T) 
  
}



