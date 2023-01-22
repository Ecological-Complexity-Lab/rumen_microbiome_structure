setwd("/Users/dafnaar/GitHub/microbiome_structure_v2")
library(tidyverse)
library(magrittr)
library(reshape2) 
library(vegan) 
library(cooccur)
library(igraph)
library(cowplot)
library(infomapecology)
check_infomap()
source('functions.R')
# co-occurrence phylogeny network analysis
# Genus----
core_ASV_30 <- read_csv('output/core_ASV_30.csv')
ASV_full_taxa <- read_csv('output/ASV_full_taxa.csv')
ASV_30_genus <- core_ASV_30 %>%
  left_join(ASV_full_taxa) %>%
  select(Farm,Cow_Code,Genus) %>%
  drop_na()

ASV_cow <- ASV_30_genus %>%
  filter(Farm=='UK1') %>% 
  distinct(Cow_Code,Genus) %>% 
  mutate(weight=1)
ASV_occurrence_mat_cow <- acast(ASV_cow, Genus~Cow_Code, fill = 0)
# Find significant pairwise co-occurrences.
ASV_co <- cooccur(ASV_occurrence_mat_cow, spp_names = TRUE)
ASV_co <- ASV_co$results
ASV_co_final <- ASV_co %>%
  as_tibble() %>% 
  mutate(level='Farm') %>% 
  mutate(level_name='UK1') %>% 
  select(sp1,sp2,weight=prob_cooccur, everything()) %>% 
  mutate(edge_type=case_when(
    p_lt < 0.05 & p_gt >= 0.05 ~ "neg",
    p_lt >= 0.05 & p_gt < 0.05 ~ "pos",
    p_lt >= 0.05 & p_gt >= 0.05 ~ "not_significant"
  )) %>%
  select(from=sp1_name, to=sp2_name, weight, edge_type, level, level_name)


# a for loop to create network for each farm
edge_list_genus <- NULL
for (i in unique(ASV_30_genus$Farm)) {
  ASV_cow <- ASV_30_genus %>%
    filter(Farm==i) %>% 
    distinct(Cow_Code,Genus) %>% 
    mutate(weight=1)
  # Create the bipartite matrix
  ASV_occurrence_mat_cow <- acast(ASV_cow, Genus~Cow_Code, fill = 0)
  # Find significant pairwise co-occurrences.
  ASV_co <- cooccur(ASV_occurrence_mat_cow, spp_names = TRUE)
  ASV_co <- ASV_co$results
  ASV_co_final <- ASV_co %>%
    as_tibble() %>% 
    mutate(level='Farm') %>% 
    mutate(level_name=i) %>% 
    select(sp1,sp2,weight=prob_cooccur, everything()) %>% 
    mutate(edge_type=case_when(
      p_lt < 0.05 & p_gt >= 0.05 ~ "neg",
      p_lt >= 0.05 & p_gt < 0.05 ~ "pos",
      p_lt >= 0.05 & p_gt >= 0.05 ~ "not_significant"
    )) %>%
    select(from=sp1_name, to=sp2_name, weight, edge_type, level, level_name)
  edge_list_genus <- rbind(edge_list_genus,ASV_co_final)
}

edge_list_genus_pos <- edge_list_genus %>%
  filter(edge_type=='pos')

# Edge weight distributions
png(filename = 'output/figures/edge_weight_dist_genus.png', width = 1600, height = 1500, res = 300)
ggplot(edge_list_genus_pos, aes(weight))+
  geom_density(alpha=0.5, color='blue', fill='blue')+
  labs(x='Edge weight', y='Density', title='Edge weight distributions')+
  theme_bw()+
  theme(panel.grid=element_blank(),
        axis.text = element_text(size=9, color='black'),
        axis.title = element_text(size=14, color='black'),
        legend.position = c(0.9,0.9))+
  facet_wrap(~level_name)
dev.off()

#network density
farm_networks_density <- NULL
for (i in unique(edge_list_genus_pos$level_name)) { 
  print(i)
  farm_multilayer_pos_dense <- edge_list_genus_pos %>%
    filter(level_name==i) %>%
    graph_from_data_frame(directed = FALSE, vertices = NULL) %>%
    edge_density()
  farm_multilayer_density <- edge_list_genus_pos %>%
    filter(level_name==i) %>%
    mutate(network_density=farm_multilayer_pos_dense) 
  farm_networks_density <- rbind(farm_networks_density,farm_multilayer_density)
}

farm_networks_density %>%
  group_by(level_name) %>%
  summarise(network_density=unique(network_density))

physical_nodes <- edge_list_genus_pos %>%
  group_by(node_id) %>%
  summarise(Modules=n_distinct(module))

state_nodes <- farm_modules_pos %>%
  select(node_id,module)

# physical and state nodes
# bind 'from' and 'to'- to include all nodes:
edge_list_genus_pos2 <-
  edge_list_genus_pos %>%
  relocate(to, from) %>%
  rename(to=from, from=to)
intra_final <- bind_rows(edge_list_genus_pos, edge_list_genus_pos2)
n_distinct(intra_final$to) # number of physical nodes
intra_final %>%
  group_by(level_name) %>%
  summarise(state_nodes=n_distinct(from)) # number of state nodes per network

# shuffled networks
# first we need to read the shuffled data:
farm_shuffled <- list.files(path = "/Users/dafnaar/GitHub/microbiome_structure_v2/HPC/shuffled/shuffle_farm_curveball_30", pattern = paste('shuff_farm',sep=""), full.names = T)
farm_shuffled_f <- sapply(farm_shuffled, read_csv, simplify=FALSE) %>%
  bind_rows(.id = "id") %>% select(-id) 
# join between the shuffled data and the taxa data
farm_shuffled_genus <- farm_shuffled_f %>%
  left_join(ASV_full_taxa) %>%
  select(c(shuff_id,Farm,Cow_Code,Genus)) %>%
  drop_na()

shuffled_co <- NULL
for(j in unique(farm_shuffled_genus$shuff_id)) {
  ASV_cow_shuff_filtered <- farm_shuffled_genus %>%
    filter(shuff_id==j) 
  # a for loop to create network for each farm
  edge_list_genus <- NULL
  for (i in unique(farm_shuffled_genus$Farm)) {
    ASV_cow_shuff <- ASV_cow_shuff_filtered %>%
      filter(Farm==i) %>% 
      distinct(Cow_Code,Genus) %>% 
      mutate(weight=1)
    # Create the bipartite matrix
    ASV_occurrence_mat_cow <- acast(ASV_cow_shuff, Genus~Cow_Code, fill = 0)
    # Find significant pairwise co-occurrences.
    ASV_co <- cooccur(ASV_occurrence_mat_cow, spp_names = TRUE)
    ASV_co <- ASV_co$results
    ASV_co_final <- ASV_co %>%
      as_tibble() %>% 
      mutate(level='Farm') %>% 
      mutate(level_name=i) %>% 
      select(sp1,sp2,weight=prob_cooccur, everything()) %>% 
      mutate(edge_type=case_when(
        p_lt < 0.05 & p_gt >= 0.05 ~ "neg",
        p_lt >= 0.05 & p_gt < 0.05 ~ "pos",
        p_lt >= 0.05 & p_gt >= 0.05 ~ "not_significant"
      )) %>%
      select(from=sp1_name, to=sp2_name, weight, edge_type, level, level_name) %>%
      mutate(id=j, .before = from) 
    edge_list_genus <- rbind(edge_list_genus,ASV_co_final)
  }
  shuffled_co <- rbind(shuffled_co,edge_list_genus)
}

shuffled_co_pos <- shuffled_co %>%
  filter(edge_type=='pos')

# Edge weight distributions for one shuffled network
png(filename = 'output/figures/edge_weight_dist_genus_shuff.png', width = 1600, height = 1500, res = 300)
shuffled_co_pos %>%
  filter(id==55) %>%
  ggplot(aes(weight))+
  geom_density(alpha=0.5, color='blue', fill='blue')+
  labs(x='Edge weight', y='Density', title='Edge weight distributions- shuffled')+
  theme_bw()+
  theme(panel.grid=element_blank(),
        axis.text = element_text(size=9, color='black'),
        axis.title = element_text(size=14, color='black'),
        legend.position = c(0.9,0.9))+
  facet_wrap(~level_name)
dev.off()

# clustering coefficient
# local
calc_CC_local <- function(x){
  g <- graph_from_data_frame(x, directed = FALSE, vertices = NULL)
  CC=transitivity(g, type = 'local')
  tibble(ASV_ID=V(g)$name, CC=CC, k=degree(g))
}

cc_degree <- edge_list_genus_pos %>%
  group_by(level_name) %>% 
  group_modify(~calc_CC_local(.x)) %>%
  drop_na()
cc_degree %>% ggplot(aes(CC,k))+geom_point()+facet_wrap(~level_name)

png(filename = 'output/figures/local_clustering_coefficient_genus_network.png', width = 1600, height = 1500, res = 300)
cc_degree %>% 
  ggplot(aes(CC))+
  geom_histogram(color='white')+
  facet_wrap(~level_name)+
  theme_bw() +
  labs(x='Clustering Coefficient', y='Count' ,title='Local Clustering Coefficient') +
  theme(axis.text = element_text(size=10, color='black'), axis.title = element_text(size=14, color='black'), title = element_text(size=15, color='black'))
dev.off()

cc_mean <- cc_degree %>%
  drop_na() %>%
  group_by(level_name) %>%
  summarise(CC_mean=mean(CC))

# clustering coefficient on one shuffled network
cc_degree_shuffled <- shuffled_co_pos %>%
  filter(id=='1') %>%
  group_by(level_name) %>% 
  select(c(from,to)) %>%
  group_modify(~calc_CC_local(.x)) %>%
  drop_na()
cc_degree_shuffled %>% ggplot(aes(CC,k))+geom_point()+facet_wrap(~level_name)

png(filename = 'output/figures/local_clustering_coefficient_one_shuffled_genus.png', width = 1600, height = 1500, res = 300)
cc_degree_shuffled %>% 
  ggplot(aes(CC))+
  geom_histogram(color='white')+
  facet_wrap(~level_name)+
  theme_bw() +
  labs(x='Clustering Coefficient', y='Count' ,title='Local Clustering Coefficient (shuffled network)') +
  theme(axis.text = element_text(size=10, color='black'), axis.title = element_text(size=14, color='black'), title = element_text(size=13, color='black'))
dev.off()

# all shuffled networks
cc_shuffled_mean <- NULL
for (i in unique(shuffled_co_pos$id)) {
  print(i)
  farm_shuffled_trans <- shuffled_co_pos %>%
    filter(id==i) %>%
    group_by(level_name) %>% 
    select(c(from,to)) %>%
    group_modify(~calc_CC_local(.x)) %>%
    mutate(id=i, .before = level_name) %>%
    drop_na() %>%
    mutate(CC_mean=mean(CC)) %>%
    select(-CC) %>%
    group_by(id, level_name) %>%
    distinct(CC_mean)
  cc_shuffled_mean <- rbind(cc_shuffled_mean,farm_shuffled_trans)
}

png(filename = 'output/figures/local_clustering_coefficient_obs_vs_shff_genus.png', width = 1600, height = 1500, res = 300)
cc_shuffled_mean %>% 
  ggplot(aes(CC_mean))+
  geom_histogram(color='white')+
  facet_wrap(~level_name)+
  theme_bw() +
  geom_vline(aes(xintercept=CC_mean), cc_mean, color= 'red') +
  labs(x='Clustering Coefficient', y='Count' ,title='Local Clustering Coefficient shuffled vs observed') +
  theme(axis.text = element_text(size=10, color='black'), axis.title = element_text(size=14, color='black'), title = element_text(size=12, color='black'))
dev.off()

# for nodes that occur in all layers
png(filename = 'output/figures/local_clustering_coefficient_sd_genus.png', width = 1600, height = 1500, res = 300)
cc_degree %>%
  group_by(ASV_ID) %>%
  mutate(number_of_layers=n_distinct(level_name)) %>%
  filter(number_of_layers==7) %>%
  filter(k>5) %>%
  group_by(ASV_ID) %>%
  summarise(cc_sd=sd(CC)) %>%
  ggplot(aes(cc_sd)) +
  geom_histogram(color='white')+
  theme_bw()+
  labs(x='CC SD', y='Count' ,title='Local Clustering Coefficient SD of nodes in layers') +
  theme(axis.text = element_text(size=14, color='black'), axis.title = element_text(size=14, color='black'), title = element_text(size=12, color='black'))
dev.off()

png(filename = 'output/figures/local_clustering_coefficient_in_farms_genus.png', width = 1600, height = 900, res = 300)
cc_degree %>%
  group_by(ASV_ID) %>%
  mutate(number_of_layers=n_distinct(level_name)) %>%
  filter(number_of_layers==7) %>%
  #filter(k>5) %>%
  ggplot(aes(x = level_name, y = ASV_ID, fill = CC)) +
  geom_tile() + 
  #geom_text()+
  theme_bw()+
  scale_fill_gradient(high = "blue", low = "light blue") +
  labs(x = '',y='', title = 'Local Clustering Coefficient in different layers')+
  theme(axis.text = element_text(size=12, color='black'), title = element_text(size = 10, color = 'black'))
dev.off()

# for density plot
cc_shuffled <- NULL
for (i in unique(shuffled_co_pos$id)) {
  print(i)
  farm_shuffled_trans <- shuffled_co_pos %>%
    filter(id==i) %>%
    group_by(level_name) %>%
    select(c(from,to)) %>%
    group_modify(~calc_CC_local(.x)) %>%
    mutate(id=i, .before = level_name) %>%
    drop_na()
  cc_shuffled <- rbind(cc_shuffled,farm_shuffled_trans)
}

cc_degree$group <- 'obs'
cc_shuffled$group <- 'shuff'
cc_combined <- rbind(cc_degree,cc_shuffled)

png(filename = 'output/figures/local_clust_coeff_density_obs_vs_shff_genus.png', width = 1600, height = 900, res = 150)
cc_combined %>% 
  filter(k>=10) %>%
  ggplot(aes(CC, fill=group))+
  geom_density(alpha=0.5)+
  labs(x='Clustering Coefficient', y='Density', title = 'Local Clustering Coefficient observed vs shuffled')+
  scale_fill_manual(values = c('purple','dark green'))+
  theme_bw()+
  facet_wrap(~level_name)+
  theme(panel.grid=element_blank(),
        axis.text = element_text(size=13, color='black'),
        axis.title = element_text(size=15, color='black'),
        legend.title = element_text(size=12, color='black'),
        legend.text = element_text(size=12, color='black'),
        title = element_text(size=17, color='black'))
dev.off()

# Eigenvector Centrality

edge_list_genus_pos_g <- edge_list_genus_pos %>%
  graph_from_data_frame(directed = FALSE, vertices = NULL) %>%
  eigen_centrality()
central_e_genus <- edge_list_genus_pos_g$vector
genus_centrality <- tibble(genus=names(central_e_genus),
                              centrality = central_e_genus)


genus_centrality_farm <- NULL
for (i in unique(edge_list_genus_pos$level_name)) { 
  print(i)
  edge_list_genus_pos_g <- edge_list_genus_pos %>%
    filter(level_name==i) %>%
    graph_from_data_frame(directed = FALSE, vertices = NULL) %>%
    eigen_centrality()
  central_e_genus <- edge_list_genus_pos_g$vector
  genus_centrality <- tibble(genus=names(central_e_genus),
                             centrality = central_e_genus, layer_name=i)
  genus_centrality_farm <- rbind(genus_centrality_farm,genus_centrality)
  
}

# now we will look only at genus occur at least in 3 layers:
genus_centrality_farm_filtered <- genus_centrality_farm %>%
  group_by(genus) %>%
  mutate(layers=n_distinct(layer_name)) %>%
  filter(layers>=3) %>%
  select(-layers)
# The farms are taken into account in ranking
png(filename = 'output/figures/centrality_in_layers_genus.png', width = 1600, height = 900, res = 150)
genus_centrality_farm_filtered %>%
  group_by(layer_name) %>%
  arrange(desc(centrality)) %>%
  mutate(rank = dense_rank(desc(centrality))) %>%
  ggplot(aes(x = layer_name, y = genus, fill = rank)) +
  geom_tile() + 
  #geom_text()+
  theme_bw()+
  scale_fill_gradient(high = "light blue", low = "blue") +
  labs(x = '',y='', title = 'Centrality of each genus in different layers')+
  theme(axis.text = element_text(size=6, color='black'), title = element_text(size = 10, color = 'black'))
dev.off()

# The farms are not taken into account in ranking
png(filename = 'output/figures/centrality_in_layers_genus_relative.png', width = 1600, height = 900, res = 150)
genus_centrality_farm_filtered %>%
  arrange(desc(centrality)) %>%
  mutate(rank = dense_rank(desc(centrality))) %>%
  ggplot(aes(x = layer_name, y = genus, fill = rank)) +
  geom_tile() + 
  #geom_text()+
  theme_bw()+
  scale_fill_gradient(high = "light blue", low = "blue") +
  labs(x = '',y='', title = 'Centrality of each genus in different layers')+
  theme(axis.text = element_text(size=6, color='black'), title = element_text(size = 10, color = 'black'))
dev.off()

# Family----
core_ASV_30 <- read_csv('output/core_ASV_30.csv')
ASV_full_taxa <- read_csv('output/ASV_full_taxa.csv')
ASV_30_family <- core_ASV_30 %>%
  left_join(ASV_full_taxa) %>%
  select(Farm,Cow_Code,Family) %>%
  drop_na()

ASV_cow <- ASV_30_family %>%
  filter(Farm=='UK1') %>% 
  distinct(Cow_Code,Family) %>% 
  mutate(weight=1)
ASV_occurrence_mat_cow <- acast(ASV_cow, Family~Cow_Code, fill = 0)
# Find significant pairwise co-occurrences.
ASV_co <- cooccur(ASV_occurrence_mat_cow, spp_names = TRUE)
ASV_co <- ASV_co$results
ASV_co_final <- ASV_co %>%
  as_tibble() %>% 
  mutate(level='Farm') %>% 
  mutate(level_name='UK1') %>% 
  select(sp1,sp2,weight=prob_cooccur, everything()) %>% 
  mutate(edge_type=case_when(
    p_lt < 0.05 & p_gt >= 0.05 ~ "neg",
    p_lt >= 0.05 & p_gt < 0.05 ~ "pos",
    p_lt >= 0.05 & p_gt >= 0.05 ~ "not_significant"
  )) %>%
  select(from=sp1_name, to=sp2_name, weight, edge_type, level, level_name)


# a for loop to create network for each farm
edge_list_family <- NULL
for (i in unique(ASV_30_family$Farm)) {
  ASV_cow <- ASV_30_family %>%
    filter(Farm==i) %>% 
    distinct(Cow_Code,Family) %>% 
    mutate(weight=1)
  # Create the bipartite matrix
  ASV_occurrence_mat_cow <- acast(ASV_cow, Family~Cow_Code, fill = 0)
  # Find significant pairwise co-occurrences.
  ASV_co <- cooccur(ASV_occurrence_mat_cow, spp_names = TRUE)
  ASV_co <- ASV_co$results
  ASV_co_final <- ASV_co %>%
    as_tibble() %>% 
    mutate(level='Farm') %>% 
    mutate(level_name=i) %>% 
    select(sp1,sp2,weight=prob_cooccur, everything()) %>% 
    mutate(edge_type=case_when(
      p_lt < 0.05 & p_gt >= 0.05 ~ "neg",
      p_lt >= 0.05 & p_gt < 0.05 ~ "pos",
      p_lt >= 0.05 & p_gt >= 0.05 ~ "not_significant"
    )) %>%
    select(from=sp1_name, to=sp2_name, weight, edge_type, level, level_name)
  edge_list_family <- rbind(edge_list_family,ASV_co_final)
}

edge_list_family_pos <- edge_list_family %>%
  filter(edge_type=='pos')

# Edge weight distributions
png(filename = 'output/figures/edge_weight_dist_family.png', width = 1600, height = 1500, res = 300)
ggplot(edge_list_family_pos, aes(weight))+
  geom_density(alpha=0.5, color='blue', fill='blue')+
  labs(x='Edge weight', y='Density', title='Edge weight distributions- family')+
  theme_bw()+
  theme(panel.grid=element_blank(),
        axis.text = element_text(size=9, color='black'),
        axis.title = element_text(size=14, color='black'),
        legend.position = c(0.9,0.9))+
  facet_wrap(~level_name)
dev.off()

#network density
farm_networks_density <- NULL
for (i in unique(edge_list_family_pos$level_name)) { 
  print(i)
  farm_multilayer_pos_dense <- edge_list_family_pos %>%
    filter(level_name==i) %>%
    graph_from_data_frame(directed = FALSE, vertices = NULL) %>%
    edge_density()
  farm_multilayer_density <- edge_list_family_pos %>%
    filter(level_name==i) %>%
    mutate(network_density=farm_multilayer_pos_dense) 
  farm_networks_density <- rbind(farm_networks_density,farm_multilayer_density)
}

farm_networks_density %>%
  group_by(level_name) %>%
  summarise(network_density=unique(network_density))

# physical and state nodes
# bind 'from' and 'to'- to include all nodes:
edge_list_family_pos2 <-
  edge_list_family_pos %>%
  relocate(to, from) %>%
  rename(to=from, from=to)
intra_final <- bind_rows(edge_list_family_pos, edge_list_family_pos2)
n_distinct(intra_final$to) # number of physical nodes
intra_final %>%
  group_by(level_name) %>%
  summarise(state_nodes=n_distinct(from)) # number of state nodes per network

# shuffled networks- family 
# first we need to read the shuffled data:
farm_shuffled <- list.files(path = "/Users/dafnaar/GitHub/microbiome_structure_v2/HPC/shuffled/shuffle_farm_curveball_30", pattern = paste('shuff_farm',sep=""), full.names = T)
farm_shuffled_f <- sapply(farm_shuffled, read_csv, simplify=FALSE) %>%
  bind_rows(.id = "id") %>% select(-id) 
# join between the shuffled data and the taxa data
farm_shuffled_family <- farm_shuffled_f %>%
  left_join(ASV_full_taxa) %>%
  select(c(shuff_id,Farm,Cow_Code,Family)) %>%
  drop_na()

shuffled_co <- NULL
for(j in unique(farm_shuffled_family$shuff_id)) {
  ASV_cow_shuff_filtered <- farm_shuffled_family %>%
    filter(shuff_id==j) 
  # a for loop to create network for each farm
  edge_list_family <- NULL
  for (i in unique(farm_shuffled_family$Farm)) {
    ASV_cow_shuff <- ASV_cow_shuff_filtered %>%
      filter(Farm==i) %>% 
      distinct(Cow_Code,Family) %>% 
      mutate(weight=1)
    # Create the bipartite matrix
    ASV_occurrence_mat_cow <- acast(ASV_cow_shuff, Family~Cow_Code, fill = 0)
    # Find significant pairwise co-occurrences.
    ASV_co <- cooccur(ASV_occurrence_mat_cow, spp_names = TRUE)
    ASV_co <- ASV_co$results
    ASV_co_final <- ASV_co %>%
      as_tibble() %>% 
      mutate(level='Farm') %>% 
      mutate(level_name=i) %>% 
      select(sp1,sp2,weight=prob_cooccur, everything()) %>% 
      mutate(edge_type=case_when(
        p_lt < 0.05 & p_gt >= 0.05 ~ "neg",
        p_lt >= 0.05 & p_gt < 0.05 ~ "pos",
        p_lt >= 0.05 & p_gt >= 0.05 ~ "not_significant"
      )) %>%
      select(from=sp1_name, to=sp2_name, weight, edge_type, level, level_name) %>%
      mutate(id=j, .before = from) 
    edge_list_family <- rbind(edge_list_family,ASV_co_final)
  }
  shuffled_co <- rbind(shuffled_co,edge_list_family)
}

shuffled_family_co_pos <- shuffled_co %>%
  filter(edge_type=='pos')

# Edge weight distributions for one shuffled network
png(filename = 'output/figures/edge_weight_dist_family_shuff.png', width = 1600, height = 1500, res = 300)
shuffled_family_co_pos %>%
  filter(id==1) %>%
  ggplot(aes(weight))+
  geom_density(alpha=0.5, color='blue', fill='blue')+
  labs(x='Edge weight', y='Density', title='Edge weight distributions- shuffled')+
  theme_bw()+
  theme(panel.grid=element_blank(),
        axis.text = element_text(size=9, color='black'),
        axis.title = element_text(size=14, color='black'),
        legend.position = c(0.9,0.9))+
  facet_wrap(~level_name)
dev.off()

# clustering coefficient
# local
calc_CC_local <- function(x){
  g <- graph_from_data_frame(x, directed = FALSE, vertices = NULL)
  CC=transitivity(g, type = 'local')
  tibble(ASV_ID=V(g)$name, CC=CC, k=degree(g))
}

cc_degree <- edge_list_family_pos %>%
  group_by(level_name) %>% 
  group_modify(~calc_CC_local(.x)) %>%
  drop_na()
cc_degree %>% ggplot(aes(CC,k))+geom_point()+facet_wrap(~level_name)

png(filename = 'output/figures/local_clustering_coefficient_family_network.png', width = 1600, height = 1500, res = 300)
cc_degree %>% 
  ggplot(aes(CC))+
  geom_histogram(color='white')+
  facet_wrap(~level_name)+
  theme_bw() +
  labs(x='Clustering Coefficient', y='Count' ,title='Local Clustering Coefficient- family network') +
  theme(axis.text = element_text(size=10, color='black'), axis.title = element_text(size=14, color='black'), title = element_text(size=14, color='black'))
dev.off()

cc_mean <- cc_degree %>%
  drop_na() %>%
  group_by(level_name) %>%
  summarise(CC_mean=mean(CC))

# clustering coefficient on one shuffled network
cc_degree_shuffled <- shuffled_family_co_pos %>%
  filter(id=='8') %>%
  group_by(level_name) %>% 
  select(c(from,to)) %>%
  group_modify(~calc_CC_local(.x)) %>%
  drop_na()
cc_degree_shuffled %>% ggplot(aes(CC,k))+geom_point()+facet_wrap(~level_name)

png(filename = 'output/figures/local_clustering_coefficient_one_shuffled_family.png', width = 1600, height = 1500, res = 300)
cc_degree_shuffled %>% 
  ggplot(aes(CC))+
  geom_histogram(color='white')+
  facet_wrap(~level_name)+
  theme_bw() +
  labs(x='Clustering Coefficient', y='Count' ,title='Local Clustering Coefficient (shuffled network)') +
  theme(axis.text = element_text(size=10, color='black'), axis.title = element_text(size=14, color='black'), title = element_text(size=13, color='black'))
dev.off()

# all shuffled networks
cc_shuffled_mean <- NULL
for (i in unique(shuffled_family_co_pos$id)) {
  print(i)
  farm_shuffled_trans <- shuffled_family_co_pos %>%
    filter(id==i) %>%
    group_by(level_name) %>% 
    select(c(from,to)) %>%
    group_modify(~calc_CC_local(.x)) %>%
    mutate(id=i, .before = level_name) %>%
    drop_na() %>%
    mutate(CC_mean=mean(CC)) %>%
    select(-CC) %>%
    group_by(id, level_name) %>%
    distinct(CC_mean)
  cc_shuffled_mean <- rbind(cc_shuffled_mean,farm_shuffled_trans)
}

png(filename = 'output/figures/local_clustering_coefficient_obs_vs_shff_family.png', width = 1600, height = 1500, res = 300)
cc_shuffled_mean %>% 
  ggplot(aes(CC_mean))+
  geom_histogram(color='white')+
  facet_wrap(~level_name)+
  theme_bw() +
  geom_vline(aes(xintercept=CC_mean), cc_mean, color= 'red') +
  labs(x='Clustering Coefficient', y='Count' ,title='Local Clustering Coefficient shuffled vs observed') +
  theme(axis.text = element_text(size=10, color='black'), axis.title = element_text(size=14, color='black'), title = element_text(size=12, color='black'))
dev.off()

# for nodes that occur in all layers
png(filename = 'output/figures/local_clustering_coefficient_sd_family.png', width = 1600, height = 1500, res = 300)
cc_degree %>%
  group_by(ASV_ID) %>%
  mutate(number_of_layers=n_distinct(level_name)) %>%
  filter(number_of_layers==7) %>%
  filter(k>5) %>%
  group_by(ASV_ID) %>%
  summarise(cc_sd=sd(CC)) %>%
  ggplot(aes(cc_sd)) +
  geom_histogram(color='white')+
  theme_bw()+
  labs(x='CC SD', y='Count' ,title='Local Clustering Coefficient SD of nodes in layers') +
  theme(axis.text = element_text(size=10, color='black'), axis.title = element_text(size=14, color='black'), title = element_text(size=12, color='black'))
dev.off()

png(filename = 'output/figures/local_clustering_coefficient_in_farms_family.png', width = 1600, height = 900, res = 300)
cc_degree %>%
  group_by(ASV_ID) %>%
  mutate(number_of_layers=n_distinct(level_name)) %>%
  filter(number_of_layers==7) %>%
  #filter(k>5) %>%
  ggplot(aes(x = level_name, y = ASV_ID, fill = CC)) +
  geom_tile() + 
  #geom_text()+
  theme_bw()+
  scale_fill_gradient(high = "blue", low = "light blue") +
  labs(x = '',y='', title = 'Local Clustering Coefficient in different layers')+
  theme(axis.text = element_text(size=10, color='black'), title = element_text(size = 9, color = 'black'))
dev.off()

# for density plot
cc_shuffled <- NULL
for (i in unique(shuffled_family_co_pos$id)) {
  print(i)
  farm_shuffled_trans <- shuffled_family_co_pos %>%
    filter(id==i) %>%
    group_by(level_name) %>%
    select(c(from,to)) %>%
    group_modify(~calc_CC_local(.x)) %>%
    mutate(id=i, .before = level_name) %>%
    drop_na()
  cc_shuffled <- rbind(cc_shuffled,farm_shuffled_trans)
}

cc_degree$group <- 'obs'
cc_shuffled$group <- 'shuff'
cc_combined <- rbind(cc_degree,cc_shuffled)

png(filename = 'output/figures/local_clust_coeff_density_obs_vs_shff_family.png', width = 1600, height = 900, res = 150)
cc_combined %>% 
  filter(k>=10) %>%
  ggplot(aes(CC, fill=group))+
  geom_density(alpha=0.5)+
  labs(x='Clustering Coefficient', y='Density', title = 'Local Clustering Coefficient observed vs shuffled')+
  scale_fill_manual(values = c('purple','dark green'))+
  theme_bw()+
  facet_wrap(~level_name)+
  theme(panel.grid=element_blank(),
        axis.text = element_text(size=13, color='black'),
        axis.title = element_text(size=15, color='black'),
        legend.title = element_text(size=12, color='black'),
        legend.text = element_text(size=12, color='black'),
        title = element_text(size=17, color='black'))
dev.off()


# Eigenvector Centrality

edge_list_family_pos_g <- edge_list_family_pos %>%
  graph_from_data_frame(directed = FALSE, vertices = NULL) %>%
  eigen_centrality()
central_e_family <- edge_list_family_pos_g$vector
family_centrality <- tibble(family=names(central_e_family),
                           centrality = central_e_family)


family_centrality_farm <- NULL
for (i in unique(edge_list_family_pos$level_name)) { 
  print(i)
  edge_list_family_pos_g <- edge_list_family_pos %>%
    filter(level_name==i) %>%
    graph_from_data_frame(directed = FALSE, vertices = NULL) %>%
    eigen_centrality()
  central_e_family <- edge_list_family_pos_g$vector
  family_centrality <- tibble(family=names(central_e_family),
                             centrality = central_e_family, layer_name=i)
  family_centrality_farm <- rbind(family_centrality_farm,family_centrality)
  
}

png(filename = 'output/figures/centrality_in_layers_family.png', width = 1600, height = 900, res = 150)
family_centrality_farm %>%
  #filter(centrality>0.1) %>%
  group_by(layer_name) %>%
  arrange(desc(centrality)) %>%
  mutate(rank = dense_rank(desc(centrality))) %>%
  #filter(rank<15) %>%
  ggplot(aes(x = layer_name, y = family, fill = rank)) +
  geom_tile() + 
  #geom_text()+
  theme_bw()+
  scale_fill_gradient(high = "light blue", low = "blue") +
  labs(x = '',y='', title = 'Centrality of each family in different layers')+
  theme(axis.text = element_text(size=6, color='black'), title = element_text(size = 10, color = 'black'))
dev.off()


