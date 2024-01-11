# ------ 14.Revision_analysis.r ----------------------------------------
# This script temporarily holds the analysis requested in the revision
# they will later be moved to more logically relevant places in the repo
#-----------------------------------------------------------------------

#------ includes ----------
library(tidyverse)
library(magrittr)
library(reshape2)
library(data.table)
library(sbm)

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

# TODO implement - also compare SBM to shuffled counterpart layer to have the histogram

## network embedding -----
# TODO to compare a layer to shuffled in another analysis, in a different scale

# ------ Interfarm level: ---------

# taxonomic beta-diversity ------ 
# or phylogeny partner fidelity
# read phylogenetic data
ASV_taxa <- read_csv('local_output/ASV_full_taxa.csv') %>% 
  select(ASV_ID, everything(), -seq16S)


# for every layer, get the the phylogenetic groups in it, for every level
for (l in layers$short_name) {
  # filter out the layer
  layer_asv <- ASV_data_30 %>% filter(Farm == l)
  layer_net <- intras %>% filter(layer == l)
  
  # add treatment of a single layer here
}

# start with one layer:
l <- "SE1"

#layer_asv <- ASV_data_30 %>% filter(Farm == layer)
layer_net <- intras %>% filter(layer == layerr)
asvs <- sort(unique(c(layer_net$node_from, layer_net$node_to)))

# filter taxa by layer
layer_taxa <- ASV_taxa %>% filter(ASV_ID %in% asvs)

sum(is.na(layer_taxa$Order))

# Handle order taxa 
ord <- layer_net %>% 
  left_join(layer_taxa, by = c('node_from' = 'ASV_ID')) %>%
  select(layer, node_from, node_to, weight, Order_from = Order) %>%
  left_join(layer_taxa, by = c('node_to' = 'ASV_ID')) %>%
  select(layer, node_from, node_to, weight, Order_from, Order_to = Order)

# use the data above to find partner fidelity?
phylo_pf <- ord %>% group_by(Order_from, Order_to) %>% summarise(count = n())


otherway <- ord %>% 
  relocate(Order_from, Order_to) %>%
  rename(Order_to=Order_from, Order_from=Order_to)

calculate_PF_T <- function(variables) {
  mat_ASV=
    x %>%
    group_by(to) %>%
    select(c(to,level_name)) %>%
    mutate(present=1) %>%
    spread(to, present, fill = 0) %>%
    column_to_rownames("level_name")
  beta_ASV <- 1-vegdist(mat_ASV, "jaccard") # Convert to similarity
  PF_J <- mean(beta_ASV)
  PF_J_sd <- sd(beta_ASV)
  num_layers <- nrow(as.matrix(beta_ASV))
  out <- data.frame(PF_J, PF_J_sd, num_layers=num_layers)
  return(out)
}

both <- bind_rows(ord, otherway) %>% 
  # keep only the ones that appear in 2 or more farms
  group_by(Order_from) %>% group_by(Order_from, Order_to) %>% summarise(count = n()) #?
  mutate(num_farms_from=n_distinct(level_name)) %>%
  filter(num_farms_from>=2)
fidelity_shuff <- both %>% 
  group_by(from) %>% 
  group_modify(~calculate_PF_T(.x))




unique(bind_rows(ord, otherway)$Order_from)


# NMI between farms ----
# TODO implement - between sbm results


# modularity - phylogenetic composition --------
# includes randomality check
# TODO implement


# NMI of clusters and hypothesis ---------
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



