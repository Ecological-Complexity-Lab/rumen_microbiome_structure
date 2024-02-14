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
library(NMI)
library(pheatmap)
library(ape)
library(dendextend)

source('functions.R')

#------ consts ---------
single_prob_file <- 'local_output/single_asv_occur_prob_80.csv'
combo_prob_file <- 'combo_asv_occur_prob.csv'

#------ functions ------
# function to calculate traits:
produce_network_traits <- function(net, grph) {
  
  # 1. Node number
  n_c <- length(unique(c(net$from, net$to)))
  
  # 2. edge number
  n_e <- nrow(net)
  
  # 3. connectivity
  den <- n_e/n_c
  
  # 4. density
  potential_edges <- n_c*(n_c-1)/2
  conn <- n_e / potential_edges
  
  # 5. Diameter (longest path)
  diam <- diameter(grph, directed = FALSE)
  
  # 6. average path length
  avg_path <- mean_distance(grph, directed = FALSE, unconnected = FALSE)
  
  # 7. betweenness
  bet <- betweenness(grph, directed = FALSE)[[1]]
  
  # 8. Mean Clustering Coefficient
  cc <- transitivity(grph, type = 'local')
  cc_mean <- mean(cc[!is.na(cc)])
  
  # 9. mean degree
  degs <- igraph::degree(grph, loops = FALSE)
  deg_mean <- mean(degs)
  
  
  # prepate output
  df_line <- data.frame(n_nodes=n_c,
                        n_edges=n_e,
                        density=den,
                        connectanse=conn,
                        diameter=diam,
                        avg_path=avg_path, 
                        between = bet, 
                        mean_cc = cc_mean, 
                        mean_degree = deg_mean)
  return(df_line)
}

# function to process dataframe for pca:
remove_constants <- function(net_to_clean) {
  # remove columns that are constant 
  to_remove <- c()
  for (col in colnames(net_to_clean)) {
    if (length(unique(net_to_clean[,col])) == 1L) {
      to_remove <- c(to_remove, col)
    }
  }
  
  return(net_to_clean %>% select(-all_of(to_remove)))
}

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


# ------ Cow level: ---------------------------------
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
  geom_histogram(aes(y = after_stat(density)), alpha=0.4, position='identity') +
  geom_vline(xintercept = 0.05, linetype="dashed", color = "red", size=0.5) + 
  labs(fill="Farm")


# ------ Farm level: --------------------------------
## SBM on layer -----
# find group number per layer in empiric network
gps <- NULL
mems_table <- NULL
for (l in layers$short_name) {
  print(l)
  lay <- intras %>% filter(layer == l)
  
  g <- graph.data.frame(lay[1:200,2:4])
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

### get SBM for multilayer - not on individual layers -----
# how: node label will be made up from farm_asv (so its a state node)
# prepare multilayer - turn physical nodes to state nodes
mln <- intras %>% mutate(from=paste(layer, node_from, sep = "_"), 
                         to=paste(layer, node_to, sep = "_"))      %>% 
                  select(from, to, weight)

g <- graph.data.frame(mln)
adj <- get.adjacency(g,sparse=FALSE, attr='weight')

# run SBM on mln
sbm_model <- BM_bernoulli$new("SBM", adj)
sbm_model$estimate() # really long step (days)
max_group <- which.max(sbm_model$ICL)
mem <- sbm_model$memberships[[max_group]]$Z

nds <- row.names(adj)
row.names(mem) <- nds

# membership documenting - find grouping
grp_mem <- apply(mem, 1, function(x) match(max(x), x)) # this will take the first group with the max value as the node's group.
mems_tbl <- tibble(asv_id=nds, membership=grp_mem)

# number of groups in the mln: 
max_group # - 42

saveRDS(sbm_model, "local_output/mln_sbm_estimate.rds")
write_csv(mems_tbl, "local_output/mln_SBM_membership_results.csv")


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
### compare the different layers ----
# observed:
# get network traits
tifs <- NULL
all_net_traits <- NULL
for (l in layers$short_name) {
  print(l)
  nett <- intras %>% filter(layer == l) %>% select(from=node_from, to=node_to, weight)
  
  # motifs:
  g <- graph.data.frame(nett)
  m <- motifs(g)
  tifs <- rbind(tifs ,m) 
  
  # other traits:
  traits <- produce_network_traits(nett, g)
  rownames(traits) <- l
  all_net_traits <- rbind(all_net_traits, traits)
}

farm_motifs <- as.data.frame(tifs) %>% select_if(~ !any(is.na(.)))
rownames(farm_motifs) <- layers$short_name

# add cow number
farm_traits <- read_csv("local_output/paper_table1.csv")
all_net_traits$cows <- farm_traits$Cows

for_pca <- merge(all_net_traits, farm_motifs, by = 'row.names', all = FALSE) 

for_pca %<>% column_to_rownames("Row.names")

for_pca <- remove_constants(for_pca)

write.csv(for_pca, "local_output/network_traits_for_pca_obs.csv", row.names=TRUE)
for_pca <- read.csv("local_output/network_traits_for_pca_obs.csv", row.names = 1)

# run the PCA analysis
res.pca <- prcomp(for_pca, scale = TRUE)
fviz_eig(res.pca)
fviz_pca_ind(res.pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE)     # Avoid text overlapping


### compare the observed layer to its shuffled ----
# read shuffled networks
parent.folder <- "HPC/shuffled/shuffle_farm_r0_30_500_jac_intra"
files <- list.files(path = parent.folder , pattern = '_edge_list.csv', recursive = T,full.names = T)

tfs_shuff <- NULL
shuff_net_traits <- NULL
for (f in files) {
  print(f)
  shuf_net <- fread(f) %>% filter(edge_type=="pos")
  s_id <- str_split_fixed(f, pattern = '/', n = 5)[4]
  l <- str_split_fixed(basename(f), pattern = '_', n = 5)[4]
  rname <- paste(s_id, l, sep = "_")
  
  # motifs:
  g <- graph.data.frame(shuf_net[,1:3])
  m <- motifs(g)
  s_motifs <- t(as.data.frame(m))
  rownames(s_motifs) <- rname
  
  tfs_shuff <- rbind(tfs_shuff ,s_motifs) 
  
  # other traits:
  traits <- produce_network_traits(shuf_net, g)
  rownames(traits) <- rname
  shuff_net_traits <- rbind(shuff_net_traits, traits)
}

# merge
to_process <- merge(shuff_net_traits, tfs_shuff, by = 'row.names') 

# save traits data to a file
write_csv(to_process, "local_output/network_traits_for_pca_shuff.csv")

# cleanup motives table - no NA or constant columns
pca_input <- remove_constants(to_process) %>% column_to_rownames("Row.names")

# run the PCA analysis per layer:
# add observed
obs_pca <- read.csv("local_output/network_traits_for_pca_obs.csv", row.names = 1) %>%
            select(-cows)
pca_input <- rbind(obs_pca, pca_input)

all_plot_data <- NULL
for (frm in layers$short_name) {
  print(frm)
  curr <- pca_input %>% filter(grepl(frm, row.names(pca_input))) %>% 
    select(-n_nodes) # networks of the same farm have the same number of nodes
  
  res.pca <- prcomp(curr, scale = TRUE)
  res.ind <- get_pca_ind(res.pca)

  to_plot <- as.data.frame(res.ind$coord) %>% 
    select(Dim.1, Dim.2) %>% add_column(type="shuff", farm=frm) 
  to_plot[1, "type"] <- "obs"
  
  all_plot_data <- rbind(all_plot_data, to_plot)
}

# plot the observed vs shuffled PCA results, per farm
pdf("local_output/figures/PCA_per_layer_output.pdf")
ggplot(all_plot_data, aes(x=Dim.1, y=Dim.2, color=type)) + 
        geom_point(aes(size=type)) +
        scale_color_manual(values=c('red', '#999999'))+
        scale_size_manual(values=c(2,1))+
        paper_figs_theme + facet_wrap(~ farm, ncol=3)
dev.off()

#plot all in one: - results from 7 different PCA run are plotted together
ggplot(all_plot_data, aes(x=Dim.1, y=Dim.2, color=type)) + 
  geom_point(aes(size=type)) +
  scale_color_manual(values=c('red', '#999999'))+
  scale_size_manual(values=c(2,1))+
  paper_figs_theme + ggtitle("All farms PCA results (7 differend PCA analysis)")


# run PCA on all the networks we have, on one scale:
res.pca_all <- prcomp(pca_input, scale = TRUE)
res.ind_all <- get_pca_ind(res.pca_all)

to_plot_all <- as.data.frame(res.ind_all$coord) %>% 
  select(Dim.1, Dim.2) %>% add_column(type="shuff") 
to_plot_all[1:7, "type"] <- "obs"

# plot it: - results from one PCA run for 507 individuals are plotted together
ggplot(to_plot_all, aes(x=Dim.1, y=Dim.2, color=type)) + 
  geom_point(aes(size=type)) +
  scale_color_manual(values=c('red', '#999999'))+
  scale_size_manual(values=c(2,1))+
  paper_figs_theme + ggtitle("PCA result for all farms")


# ------ Inter-farm level: --------------------------
## taxonomic beta-diversity ------ 
# or phylogeny partner fidelity

# read phylogenetic data
ASV_taxa <- read_csv('local_output/ASV_full_taxa.csv') %>% 
  select(ASV_ID, everything(), -seq16S)

# filter only taxa that exist in the networks
asvs <- sort(unique(c(intras$node_from, intras$node_to)))
all_taxa <- ASV_taxa %>% filter(ASV_ID %in% asvs)

# collect data:
c <- get_taxa_pf(intras, all_taxa, "Class") %>% add_column(taxa="Class")
o <- get_taxa_pf(intras, all_taxa, "Order") %>% add_column(taxa="Order")
f <- get_taxa_pf(intras, all_taxa, "Family") %>% add_column(taxa="Family")
g <- get_taxa_pf(intras, all_taxa, "Genus") %>% add_column(taxa="Genus")

PF_T_obs <- rbind(c, o, f, g) %>% add_column(run=0)

ggplot(PF_T_obs, aes(PF_T, fill=taxa)) +
  geom_histogram(alpha=0.5, color='white', position="identity")+
  labs(x='Taxa Partner fidelity score', y='Count')+
  paper_figs_theme +
  theme(panel.grid=element_blank(),
        axis.text = element_text(size=10, color='black'),
        axis.title = element_text(size=10, color='black'))

# read shuffled taxa PF results: 
# Folder from HPC containing the 001-500 sub-folders
parent.folder <- "HPC/shuffled/shuffle_farm_r0_30_500_jac_intra"
sub.folders <- list.dirs(parent.folder, recursive=TRUE)[-1]

# read all shuffled networks
PF_T_shuff <- NULL
for (dir in sub.folders) {
  print(dir)
  shuff_fid <- fread(paste(dir,"/taxa_pf_shuff_farm_30.csv", sep="")) # Faster to read with this
  PF_T_shuff <- rbind(PF_T_shuff, shuff_fid)
}
PF_T_shuff <- as_tibble(PF_T_shuff)

write_csv(PF_T_shuff, 'local_output/PF_T_pos_30_shuffled_r0.csv')
PF_T_shuff <- read_csv('local_output/PF_T_pos_30_shuffled_r0.csv')

# calculate Z score: includes all taxa t
PF_T_z_score <- 
  PF_T_shuff %>%
  group_by(taxa_from, taxa) %>%
  summarise(PF_T_shuff_mean=mean(PF_T), PF_T_shuff_sd=sd(PF_T)) %>% 
  inner_join(PF_T_obs %>% select(taxa_from, PF_T)) %>%
  mutate(z=(PF_T-PF_T_shuff_mean)/PF_T_shuff_sd) %>% 
  mutate(signif=case_when(z>1.96 ~ 'above', # Obs is more than the shuffled
                          z< -1.96 ~ 'below', # Obs is lower than the shuffled
                          z<=1.96 | z>=-1.96 ~ 'not signif'))

# What proportion of ASVs have a statistical significant PF_J?
PF_T_z_score %>% 
  group_by(signif) %>% 
  summarise(n=n(),prop=n/nrow(PF_T_z_score))


# Pie figure:
PF_T_z_score %>% group_by(taxa, signif) %>% 
  summarise(n=n()) %>% 
  mutate(prop = n / sum(n)) %>% 
  mutate(N = sum(n)) %>% 
  mutate(signif=factor(signif, levels=c('not signif','below','above'))) %>% 
  mutate(ypos = cumsum(prop)- 0.5*prop ) %>% 
  ggplot(aes(x="", y=prop, fill=signif))+
  facet_wrap(~taxa)+
  geom_bar(stat="identity", width=1) +
  scale_fill_manual(values = c('blue','orange','#32a852'))+
  coord_polar("y", start=0)+
  paper_figs_theme_no_legend+theme_void()


## NMI between farms ----
# read groups per farm (created up in this script)
mems_table <- read_csv("local_output/layer_SBM_membership_results.csv")

infomap_table <- read_csv("local_output/farm_modules_pos_30_U.csv") %>%
  select(farm=short_name, node_name, module) # infomap with unifrec interlayer edges
  # note that here the is a multilayer network. 
  # meaning module 6 in one layer is the same module as module 6 in another layer.

NMI_from_membership <- function(table, names) {
  NMIs <- matrix(0, nrow = length(names), ncol = length(names))
  colnames(NMIs) <- rownames(NMIs) <- names
  
  for (i in rownames(NMIs)) {
    farm1 <- table %>% filter(farm==i) %>% select(-farm)
      for (j in colnames(NMIs)) {
        farm2 <- table %>% filter(farm==j) %>% select(-farm)
        
        NMIs[i,j] <- NMI(farm1, farm2)$value
      }
  }
  return(NMIs)
}


sbms     <- NMI_from_membership(mems_table, layers$short_name)
infomaps <- NMI_from_membership(infomap_table, layers$short_name) 
            # Note: NMI = 0 if all the labels of the layer are in the same group
# save results
write_csv(as.data.frame(sbms), "local_output/NMI_SBM_layers_30")
write_csv(as.data.frame(infomaps), "local_output/NMI_Infomap_layers_30")

# read results to present
sbms <- read_csv("local_output/NMI_SBM_layers_30")
infomaps <- read_csv("local_output/NMI_Infomap_layers_30")

pheatmap(sbms)
pheatmap(infomaps) # maybe do this without removing small modules?


## modularity - phylogenetic composition --------
# get membership data, merge with edgelist - only infomap
infomap_table <- read_csv("local_output/farm_modules_pos_30_U.csv") %>%
  select(farm=short_name, asv_id=node_name, membership=module)

nets <- intras %>% 
  left_join(infomap_table, by = c('layer'='farm', 'node_from'='asv_id')) %>%
  select(layer, node_from, node_to, weight, from_module = membership) %>%
  left_join(infomap_table, by = c('layer'='farm', 'node_to'='asv_id')) %>%
  select(layer, node_from, node_to, weight, from_module, to_module = membership) %>% 
  mutate(same_module=from_module==to_module)

# filter only taxa that exist in the networks
asvs <- sort(unique(c(intras$node_from, intras$node_to)))

# read and process tree:
phylo_tree <- readRDS("local_output/fitted_asvs_phylo_tree.rds")
tree <- phylo_tree$tree
# prune the tree
unincluded <- tree$tip.label[!tree$tip.label %in% asvs]
pruned <- dendextend::prune(tree, unincluded)

# calculate asv distances
distances <- cophenetic.phylo(pruned)

# get edges distance:
nets$phylo_dist <- apply(intras, 1, function(x) distances[x[2], x[3]])

write_csv(nets, "local_output/modules_phylogenetic_composition.csv")

# read results for plotting
nets <- read_csv("local_output/modules_phylogenetic_composition.csv")


# plot distance distribution between in module and out module links 
nets %>%
  ggplot( aes(x=phylo_dist, fill=same_module)) +
  geom_histogram(color="#e9ecef", alpha=0.6, position = 'stack') +
  scale_fill_manual(values=c("#69b3a2", "#404080")) +
  labs(fill="") + ggtitle("Edges between and within modules")

# plot distance distribution between modules
nets %>% filter(same_module == TRUE) %>% select(from_module, phylo_dist) %>%
  transform(from_module=as.character(from_module)) %>%
  ggplot(aes(x=phylo_dist, fill=from_module)) +
  geom_histogram(color="#e9ecef", alpha=0.6, position = 'identity') +
  labs(fill="Modul number") + ggtitle("Edges within each module")


## NMI of clusters and hypothesis ---------
infomap_table <- read_csv("local_output/farm_modules_pos_30_U.csv") %>%
  select(farm=short_name, asv_id=node_name, membership=module)

# built hypothesis table:
# H1 + H2 + H3:
Hs <- infomap_table %>% 
               mutate(H1=farm) %>% 
               mutate(H2=case_when(farm %in% c("FI1", "SE1") ~ "north",
                                  !farm %in% c("FI1", "SE1") ~ "south")) %>%
               add_column(H3=1) %>% 
               mutate(label=paste(farm, asv_id, sep = "_"))

NMI(Hs %>% select(label, membership, -farm),
    Hs %>% select(label, H1, -farm)) # 0.8566478 - higher value
NMI(Hs %>% select(label, membership, -farm),
    Hs %>% select(label, H2, -farm)) # 0.5038392
NMI(Hs %>% select(label, membership, -farm),
    Hs %>% select(label, H3, -farm)) # 0 - because its one big group

# using sbm results:
# use multilayer SBM method for this
membs <- read_csv("local_output/mln_SBM_membership_results.csv") # SBM with labels

# built hypothesis table:
# H1 + H2 + H3:
Hs <- membs %>% rename(label=asv_id) %>%
  separate(label, c('farm'), remove = FALSE, extra = 'drop') %>% 
  mutate(H1=farm) %>% 
  mutate(H2=case_when(farm %in% c("FI1", "SE1") ~ "north",
                      !farm %in% c("FI1", "SE1") ~ "south")) %>%
  add_column(H3=1)

NMI(Hs %>% select(label, membership, -farm),
    Hs %>% select(label, H1, -farm)) # 0.6659741 - higher value
NMI(Hs %>% select(label, membership, -farm),
    Hs %>% select(label, H2, -farm)) # 0.2618777
NMI(Hs %>% select(label, membership, -farm),
    Hs %>% select(label, H3, -farm)) # 0 - because its one big group

# --- not sure will be done --------------------------

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



