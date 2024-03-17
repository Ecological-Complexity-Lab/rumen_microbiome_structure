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
library(factoextra)
library(infomapecology)

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
  conn <- n_e/n_c
  
  # 4. density
  potential_edges <- n_c*(n_c-1)/2
  den <- n_e / potential_edges
  
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
  net_to_clean <- as.data.frame(net_to_clean)
  # remove columns that are constant 
  to_remove <- c()
  for (col in colnames(net_to_clean)) {
    if (length(unique(net_to_clean[,col])) == 1L) {
      to_remove <- c(to_remove, col)
    }
  }
  
  return(net_to_clean %>% select(-all_of(to_remove)))
}

# produce a pairwise NMI values between the farms
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

#------ run ------------
## Read data to be used -----
# read network:
farm_multilayer_pos_30 <- read_csv('local_output/farm_multilayer_pos_30.csv')
all_nodes <- sort(unique(c(farm_multilayer_pos_30$from, farm_multilayer_pos_30$to)))
all_nodes <- tibble(node_id=1:length(all_nodes), node_name=all_nodes)
layers <- tibble(layer_id=1:7,
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
# calculate density for observed farms: ----
obs_den <- NULL
for (f in layers$short_name) {
  print(f)
  fa <- intras %>% filter(layer == f)
  
  # density of one farm:
  n_c <- length(unique(c(fa$node_from, fa$node_to))) # Node number
  n_e <- nrow(fa) # edge number
  
  # density
  potential_edges <- n_c*(n_c-1)/2
  den <- n_e / potential_edges
  
  new_line <- tibble(farm=f, run="000", density=den)
  obs_den <- rbind(obs_den, new_line)
}

# calculate density for each shuffled farm
# read shuffled networks
parent.folder <- "HPC/shuffled/shuffle_farm_r0_30_500_jac_intra"
files <- list.files(path = parent.folder , pattern = '_edge_list.csv', recursive = T,full.names = T)

shuff_net_density <- NULL
for (f in files) {
  print(f)
  shuf_net <- fread(f) %>% filter(edge_type=="pos")
  s_id <- str_split_fixed(f, pattern = '/', n = 5)[4]
  l <- str_split_fixed(basename(f), pattern = '_', n = 5)[4]

  # density of one farm:
  n_c <- length(unique(c(shuf_net$from, shuf_net$to))) # Node number
  n_e <- nrow(shuf_net) # edge number
  
  # density
  potential_edges <- n_c*(n_c-1)/2
  den <- n_e / potential_edges
  
  new_line <- tibble(farm=l, run=s_id, density=den)
  shuff_net_density <- rbind(shuff_net_density, new_line)
}

# plot density
ggplot(data = shuff_net_density, aes(x=density)) +
  geom_histogram(fill = "steelblue") +
  labs(y = "Count", x = "Network Density") +
  geom_vline(data = obs_den, mapping = aes(xintercept = density), 
             colour="#BB0000", linetype="dashed") +
  facet_grid(farm ~ .)


## SBM on layer -----
### find group number per layer in empiric network ----
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

pdf("local_output/figures/SI_SBM_group_size.pdf", 6, 4)
mems_table %>% group_by(farm, membership) %>% summarise(n=n()) %>%
  ggplot(aes(x=membership, y=n)) + 
  geom_bar(stat="identity") + paper_figs_theme + 
  facet_wrap(. ~ farm, ncol=3) + labs(x="Group ID", y="Number of ASVs") +
  geom_text(groups %>% rename(farm=short_name),
            mapping = aes(x=0, y=75, label=paste("k=", emp_max_ICL, sep = "")),
            hjust = 0)
dev.off()

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

### get SBM for multilayer - with unifrac intelayer edges -----
multilayer_unif <- read_csv('local_output/multilayer_unif.csv')
mln <- multilayer_unif %>% # get layers names
  left_join(layers, by = c('layer_from'='layer_id')) %>%
  select(node_from, node_to, weight, from_farm = short_name, layer_to, type) %>%
  left_join(layers, by = c('layer_to'='layer_id')) %>%
  select(node_from, node_to, weight, from_farm, to_farm=short_name, type) %>%
  # prepare format for sbm
  mutate(from=paste(from_farm, node_from, sep = "_"), 
         to=paste(to_farm, node_to, sep = "_"))         %>% 
  select(from, to, weight) 

g <- graph.data.frame(mln)
adj <- get.adjacency(g, sparse=FALSE, attr='weight')

# run SBM on mln with interlayer edges
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
max_group # - 44

saveRDS(sbm_model, "local_output/mln_interlayer_sbm_estimate.rds")
write_csv(mems_tbl, "local_output/mln_interlayer_SBM_membership_results.csv")


# --- shuffled data from the HPC ? ----

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

# Visualize PCA
fviz_eig(res.pca, addlabels = TRUE)
fviz_pca_var(res.pca, col.var = "cos2",
             gradient.cols = c("black", "orange", "green"),
             repel = TRUE)
fviz_pca_ind(res.pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE)     # Avoid text overlapping
fviz_cos2(res.pca, choice = "var", axes = 1:2)
summary(res.pca)


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
to_process <- read_csv("local_output/network_traits_for_pca_shuff.csv")

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


# run PCA on all the networks we have, on one scale:
res.pca_all <- prcomp(pca_input, scale = TRUE)
res.ind_all <- get_pca_ind(res.pca_all)

to_plot_all <- as.data.frame(res.ind_all$coord) %>% 
  select(Dim.1, Dim.2) %>% add_column(Type="shuff") 
to_plot_all[,"Farm"] <- str_split_fixed(rownames(to_plot_all), pattern="_", n=2)[,2]
to_plot_all[1:7,"Farm"] <- rownames(to_plot_all)[1:7]
to_plot_all[1:7, "Type"] <- "obs"

# plot it: - results from one PCA run for 507 individuals are plotted together
farm_pca <- ggplot(to_plot_all, aes(x=Dim.1, y=Dim.2, color=Farm)) + 
  geom_point(aes(shape=Type, size=Type)) +
  scale_shape_manual(values=c(16, 3)) +
  scale_size_manual(values=c(3, 1)) +
  paper_figs_theme #+ ggtitle("PCA result for all farms")

pdf("local_output/figures/SI_farm_pca.pdf",5,4)
farm_pca
dev.off()

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

PF_T_obs <- rbind(c, o, f, g) %>% add_column(run=0) %>% 
                mutate(taxa = fct_relevel(taxa, 
                      "Class", "Order", "Family", "Genus"))

ggplot(PF_T_obs, aes(PF_T, fill=taxa)) +
  geom_histogram(alpha=0.5, color='white', position="identity")+
  labs(x='Taxa Partner fidelity score', y='Count')+
  paper_figs_theme +
  theme(panel.grid=element_blank(),
        axis.text = element_text(size=10, color='black'),
        axis.title = element_text(size=10, color='black'))

p <- ggplot(PF_T_obs, aes(PF_T, fill=taxa)) +
  geom_histogram(alpha=0.5, color='white', position="identity")+
  labs(x='Bray-Curtis score', y='Count')+
  paper_figs_theme +
  theme(panel.grid=element_blank(),
        axis.text = element_text(size=10),
        axis.title = element_text(size=10), 
        legend.position = "none")+
  facet_grid(taxa ~ .)

pdf("local_output/figures/SI_taxa_pf.pdf",4,5)
p
dev.off()

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

sbms     <- NMI_from_membership(mems_table, layers$short_name)
infomaps <- NMI_from_membership(infomap_table, layers$short_name) 
            # Note: NMI = 0 if all the labels of the layer are in the same group
# save results
write.csv(as.data.frame(sbms), "local_output/NMI_SBM_layers_30", row.names = TRUE)
write.csv(as.data.frame(infomaps), "local_output/NMI_Infomap_layers_30", row.names = TRUE)

# read results to present
sbms <- read.csv("local_output/NMI_SBM_layers_30", row.names = 1)
pheatmap(sbms, treeheight_row=0 , treeheight_col=0, 
         filename="local_output/figures/farm_sbm_nmi_paper.pdf", width = 3.3,height = 2.8)

infomaps <- read.csv("local_output/NMI_Infomap_layers_30", row.names = 1)
pheatmap(infomaps) # maybe do this without removing small modules?


# Read NMI for sbm grouping + Jaccard between farms
sbms <- read.csv("local_output/NMI_SBM_layers_30", row.names = 1)
jac <- read.csv("local_output/jaccard_beta_diversity_30", row.names = 1)

# perform mantel test: (spearman)
m_res1 <- mantel(sbms, jac, method = "spearman", permutations = 9999)
m_res1 # Mantel statistic r: 0.587 , p-val: 0.005754 

# perform mantel test: (pearson)
m_res2 <- mantel(sbms, jac, permutations = 9999)
m_res2 # Mantel statistic r: 0.6009 , p-val: 0.003373 


# compare to SBM MNI to shuffled memberships ----
mems_table <- read_csv("local_output/layer_SBM_membership_results.csv")
sbms <- read.csv("local_output/NMI_SBM_layers_30", row.names = 1)

# prepare the matrix with shuffled NMI values:
result <- array(0, dim = c(7,7,1000))
rownames(result) <- colnames(result) <- layers$short_name

# shuffle membership then run nmi 1000 times
for (i in 1:1000) {
  curr_shuff <- NULL
  #reesembel shuffled farms memberships
  for (f in layers$short_name) {
    l <- mems_table %>% filter(farm==f)
    shuffled <- sample(l$membership)
    
    sh_farm <- cbind(l[,1:2], shuffled)
    
    curr_shuff <- rbind(curr_shuff, sh_farm)
  }

  shuff_mnis <- NMI_from_membership(curr_shuff, layers$short_name)
  result[,,i] <- shuff_mnis
  
  print(i)
}

# compare shuffled to observed:

#check p-value per farm pair
pvals <- array(1, dim = c(7,7))
rownames(pvals) <- colnames(pvals) <- layers$short_name

for (f1 in layers$short_name) {
  for (f2 in layers$short_name) {
    vals <- result[f1,f2,]
    obs <- sbms[f1,f2]
    one_pval <- sum(vals>obs)/1000
    pvals[f1,f2] <- one_pval
  }
}

#turn tensore into long format (to visualize)
shuff_nmis_long <- NULL
for (i in 1:1000) {
  mat <- result[,,i]
  mat[!lower.tri(mat)] <- 0
  
  long <- reshape2::melt(mat) %>% filter(value > 0) %>% 
    select(farm1=Var1, farm2=Var2, nmi_val=value) %>% add_column(iter=i)
  shuff_nmis_long <- rbind(shuff_nmis_long, long)
}

write_csv(shuff_nmis_long, "local_output/NMI_SBM_shuffled_membership_30.csv")

# plot the NMI comparosones
f_ls <- c("UK1", "UK2", "IT1", "IT2", "IT3", "FI1", "SE1")
shuff_nmis_long <- as.data.frame(read_csv("local_output/NMI_SBM_shuffled_membership_30.csv"))
shuff_nmis_long$farm1 <- factor(shuff_nmis_long$farm1, levels = f_ls)
shuff_nmis_long$farm2 <- factor(shuff_nmis_long$farm2, levels = f_ls)
sbms <- read.csv("local_output/NMI_SBM_layers_30", row.names = 1)

# plot the histograms
# prepare observed data be long format for vline
mat_obs <- as.matrix(sbms)
mat_obs[!lower.tri(mat_obs)] <- 0
long_obs <- reshape2::melt(mat_obs) %>% filter(value > 0) %>% 
  select(farm1=Var1, farm2=Var2, nmi_val=value) %>% add_column(iter=0)
# plot
p <- ggplot(data = shuff_nmis_long, aes(x=nmi_val)) +
  geom_histogram(fill = "steelblue") +
  labs(y = "Count", x = "NMI value") +
  geom_vline(data = long_obs, mapping = aes(xintercept = nmi_val), 
             colour="#BB0000", linetype="dashed") +
  facet_grid(farm1 ~ farm2)

pdf("local_output/figures/SI_farm_nmi_vs_shuff.pdf", 5, 4)
p
dev.off()

# extract NMI range:
by_val <- long_obs %>% arrange(nmi_val)
lowest <- by_val[1,]   # FI1, IT3: 0.2165693
highest <- by_val[21,] # IT2, UK1: 0.5793949

## modularity - phylogenetic composition --------
# prepare phylogenetic data to calculate distances
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


# get membership data, merge with edgelist - only infomap
infomap_table <- read_csv("local_output/farm_modules_pos_30_U.csv") %>%
  select(farm=short_name, asv_id=node_name, membership=module)

nets <- intras %>% 
  left_join(infomap_table, by = c('layer'='farm', 'node_from'='asv_id')) %>%
  select(layer, node_from, node_to, weight, from_module = membership) %>%
  left_join(infomap_table, by = c('layer'='farm', 'node_to'='asv_id')) %>%
  select(layer, node_from, node_to, weight, from_module, to_module = membership) %>% 
  mutate(same_module=from_module==to_module)

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
  labs(fill="") + ggtitle("Edges between and within modules - Infomap")

# plot distance distribution between modules
nets %>% filter(same_module == TRUE) %>% select(from_module, phylo_dist) %>%
  transform(from_module=as.character(from_module)) %>%
  ggplot(aes(x=phylo_dist, fill=from_module)) +
  geom_histogram(color="#e9ecef", alpha=0.6, position = 'stack') +
  labs(fill="Modul number") + ggtitle("Edges within each module - Infomap")


# Do the same but with SBM
# read SBM mln results:
membs <- read_csv("local_output/mln_SBM_membership_results.csv") %>% # SBM with labels
          separate(asv_id, c('farm','asv', 'id')) %>% 
          unite(asv_id, c(asv, id), sep = "_", remove = TRUE)
nets_sbm <- intras %>% 
  left_join(membs, by = c('layer'='farm', 'node_from'='asv_id')) %>%
  select(layer, node_from, node_to, weight, from_module = membership) %>%
  left_join(membs, by = c('layer'='farm', 'node_to'='asv_id')) %>%
  select(layer, node_from, node_to, weight, from_module, to_module = membership) %>% 
  mutate(same_module=from_module==to_module)

# get edges distance:
nets_sbm$phylo_dist <- apply(intras, 1, function(x) distances[x[2], x[3]])

# save point
write_csv(nets_sbm, "local_output/modules_phylogenetic_composition_SBM.csv")
nets_sbm <- read_csv("local_output/modules_phylogenetic_composition_SBM.csv")

# plotting
# plot distance distribution between in module and out module links 
nets_sbm %>%
  ggplot( aes(x=phylo_dist, fill=same_module)) +
  geom_histogram(color="#e9ecef", alpha=0.6, position = 'stack') +
  scale_fill_manual(values=c("#69b3a2", "#404080")) +
  labs(fill="") + ggtitle("Edges between and within modules - SBM")

# plot distance distribution between modules
nets_sbm %>% filter(same_module == TRUE) %>% select(from_module, phylo_dist) %>%
  transform(from_module=as.character(from_module)) %>%
  ggplot(aes(x=phylo_dist, fill=from_module)) +
  geom_histogram(color="#e9ecef", alpha=0.6, position = 'stack') +
  labs(fill="Modul number") + ggtitle("Edges within each module - SBM")

# plot distance distribution between farms 
nets_sbm %>%
  ggplot( aes(x=phylo_dist, fill=layer)) +
  geom_histogram(color="#e9ecef", alpha=0.6, position = 'stack') +
  #scale_fill_manual(values=c("#69b3a2", "#404080")) +
  labs(fill="") + ggtitle("Edge distance in each farm")

## NMI of clusters and hypothesis ---------
### infomap -----
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

### SBM (no interlayer edges) -----
# use multilayer SBM method for this
membss <- read_csv("local_output/mln_SBM_membership_results.csv") # SBM with labels

# built hypothesis table:
# H1 + H2 + H3:
Hs <- membss %>% rename(label=asv_id) %>%
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

# check asvs per farm per module
data <- membss %>% 
  separate(asv_id, c('farm','asv', 'id')) %>% 
  unite(asv_id, c(asv, id), sep = "_", remove = TRUE) %>%
  mutate(farm=factor(farm, levels = c("UK1","UK2","IT1","IT2","IT3","FI1",'SE1'))) %>%
  group_by(farm) %>%
  mutate(nodes_in_layers=n_distinct(asv_id)) %>%
  group_by(farm, membership) %>%
  mutate(nodes_in_modules=n_distinct(asv_id)) %>%
  mutate(nodes_percent=nodes_in_modules/nodes_in_layers) %>%
  distinct(farm, membership, nodes_percent) %>% 
  arrange(membership, farm)
# plot
data %>% filter(nodes_percent >0.03) %>%
  ggplot(aes(x = membership, y = farm, fill=nodes_percent))+
  geom_tile(color='white')+
  scale_x_continuous(breaks = seq(1, max(membs$membership), 1))+
  scale_fill_viridis_c(limits = c(0, 1))+
  theme_bw()+
  labs(x='Module ID', y='', title = "with threshold")+
  theme(panel.grid=element_blank(),
        axis.text = element_text(size=10, color='black'),
        axis.title = element_text(size=10, color='black'),
        title = element_text(size=10, color='black'),
        plot.tag = element_text(face = "bold")) +
  paper_figs_theme_no_legend


### SBM (WITH interlayer edges) -----
memstbl <- read_csv("local_output/mln_interlayer_SBM_membership_results.csv")

# built hypothesis table:
# H1 + H2 + H3:
Hs <- memstbl %>% rename(label=asv_id) %>%
  separate(label, c('farm'), remove = FALSE, extra = 'drop') %>% 
  mutate(H1=farm) %>% 
  mutate(H2=case_when(farm %in% c("FI1", "SE1") ~ "north",
                      !farm %in% c("FI1", "SE1") ~ "south")) %>%
  add_column(H3=1)

NMI(Hs %>% select(label, membership, -farm),
    Hs %>% select(label, H1, -farm)) # 0.6549356 - higher value
NMI(Hs %>% select(label, membership, -farm),
    Hs %>% select(label, H2, -farm)) # 0.2584942
NMI(Hs %>% select(label, membership, -farm),
    Hs %>% select(label, H3, -farm)) # 0 - because its one big group

# check asvs per farm per module
data2 <- memstbl %>% 
  separate(asv_id, c('farm', 'id')) %>% 
  mutate(farm=factor(farm, levels = c("UK1","UK2","IT1","IT2","IT3","FI1",'SE1'))) %>%
  group_by(farm) %>%
  mutate(nodes_in_layers=n_distinct(id)) %>%
  group_by(farm, membership) %>%
  mutate(nodes_in_modules=n_distinct(id)) %>%
  mutate(nodes_percent=nodes_in_modules/nodes_in_layers) %>%
  distinct(farm, membership, nodes_percent) %>% 
  arrange(membership, farm)
# plot
data2 %>% filter(nodes_percent >0.03) %>%
  ggplot(aes(x = membership, y = farm, fill=nodes_percent))+
  geom_tile(color='white')+
  scale_x_continuous(breaks = seq(1, max(membs$membership), 1))+
  scale_fill_viridis_c(limits = c(0, 1))+
  theme_bw()+
  labs(x='Module ID', y='', title = "With interlayer and with threshold")+
  theme(panel.grid=element_blank(),
        axis.text = element_text(size=10, color='black'),
        axis.title = element_text(size=10, color='black'),
        title = element_text(size=10, color='black'),
        plot.tag = element_text(face = "bold")) +
  paper_figs_theme_no_legend



## Infomap of MLN with interlayer edges TO COMMON partners -----
# this *should* negate the effect that multiple intra-edges have in 
# overwhelming the infomap towards inside the farm.

# function generating interlayer edges for a pair:
get_inter_for_farm_pair <- function(intra, f1, f2) {
  # get farms
  Farm1 <- intra %>% filter(layer == f1)
  Farm2 <- intra %>% filter(layer == f2)
  
  # filter only ASVs that appear in both farms
  f1_asvs <- unique(c(Farm1$node_from, Farm1$node_to))
  f2_asvs <- unique(c(Farm2$node_from, Farm2$node_to))
  intrsc <- intersect(f1_asvs, f2_asvs)
  
  # filter edges with nodes that appear in both farms 
  relevant1 <- Farm1 %>% filter((node_from %in% intrsc) & (node_to %in% intrsc))
  relevant2 <- Farm2 %>% filter((node_from %in% intrsc) & (node_to %in% intrsc))
  
  # add flipped edges to allow accessing the nodes using only "from"
  dbl_1 <- rbind(relevant1 %>% select(node_from, node_to),
                 relevant1 %>% select(node_from=node_to, node_to=node_from)) 
  dbl_2 <- rbind(relevant2 %>% select(node_from, node_to),
                 relevant2 %>% select(node_from=node_to, node_to=node_from))
  
  output <- NULL
  #go over only the filtered node that appear in one of the farms:
  for (node in intrsc) {
    # filter only edges of the current node
    node_edgesin_1 <- dbl_1 %>% filter(node_from==node)
    node_edgesin_2 <- dbl_2 %>% filter(node_from==node)
    
    # get only edges that appear in both farms
    common_edges_from_1 <- node_edgesin_1 %>% 
      filter(node_to %in% node_edgesin_2$node_to)
    # note: no need to do the other direction as well as it will be the same.
    
    # produce the interlayer edges:
    node_inter <- tibble(layer_from=f1, node_from=node, 
                         layer_to=f2,   node_to=common_edges_from_1$node_to)
    
    output <- rbind(output, node_inter)
  }
  
  return(output)
}

get_common_partners_mln <- function(intras, layers) {
  inter_all <- NULL
  for (f1 in layers$short_name) {
    farm_ind <- 1
    f2 = layers$short_name[farm_ind]
    while (f2 != f1) { # to make sure we go through a pair of 
      print(paste("f1:", f1, ", f2:", f2))
      # get the set of new interlayer edges:
      inter_pair <- get_inter_for_farm_pair(intras, f1, f2)
      inter_all <- rbind(inter_all, inter_pair)
      
      # for next iteration
      farm_ind <- farm_ind + 1
      f2 = layers$short_name[farm_ind]
    }
  }
  
  # run infomap on this new mln
  intra_big <- intras %>% select(layer_from=layer, node_from,
                                 layer_to=layer, node_to)
  mln <- rbind(intra_big, inter_all) %>% add_column(weight=1)
  
  return(mln)
}

infomap_from_mln <- function(filepath){
  # prepare data for infomap run (30%)
  mln_mile <- read_csv(filepath) 
  all_nodes <- sort(unique(c(mln_mile$node_from, mln_mile$node_to)))
  all_nodes <- tibble(node_id=1:length(all_nodes), node_name=all_nodes)
  layers <- tibble(layer_id=1:7, 
                   short_name=c('UK1', 'UK2', 'IT1', 'IT2', 'IT3', 'FI1', 'SE1'))
  
  nrow(mln_mile %>% filter(layer_from==layer_to)) # intra: 269,285
  nrow(mln_mile %>% filter(layer_from!=layer_to)) # inter: 287,900
  nrow(mln_mile) # all: 557,185
  
  # make network be represented by indexes
  mln_inds <- 
    mln_mile %>%
    left_join(layers, by = c('layer_from' = 'short_name')) %>% 
    left_join(layers, by = c('layer_to' = 'short_name')) %>% 
    left_join(all_nodes, by = c('node_from' = 'node_name')) %>% 
    left_join(all_nodes, by = c('node_to' = 'node_name')) %>% 
    select(layer_from=layer_id.x, node_from=node_id.x, layer_to=layer_id.y, node_to=node_id.y, weight)
  
  # Run Infomap
  multilayer_for_infomap <- create_multilayer_object(extended = mln_inds, 
                                                     nodes = all_nodes, 
                                                     layers = layers)
  m <- infomapecology::run_infomap_multilayer(multilayer_for_infomap, silent = F,
                                              flow_model = 'undirected',  
                                              trials = 200, relax = F, seed=111)
  return(m)
}

# check asvs per farm per module
modules_to_plot <- function(modules_df) {
  data <- modules_df %>% left_join(layers) %>% 
    mutate(farm=factor(short_name, levels = c("UK1","UK2","IT1","IT2","IT3","FI1",'SE1'))) %>%
    group_by(farm) %>%
    mutate(nodes_in_layers=n_distinct(node_name)) %>%
    group_by(farm, module) %>%
    mutate(nodes_in_modules=n_distinct(node_name)) %>%
    mutate(nodes_percent=nodes_in_modules/nodes_in_layers) %>%
    distinct(farm, module, nodes_percent) %>% 
    arrange(module, farm)
  # plot
  p <- data %>% filter(nodes_percent > 0.03) %>%
    ggplot(aes(x = module, y = farm, fill=nodes_percent))+
    geom_tile(color='white')+
    scale_x_continuous(breaks = seq(1, max(data$module), 1))+
    scale_fill_viridis_c(limits = c(0, 1))+
    theme_bw()+
    labs(x='Module ID', y='')+
    theme(panel.grid=element_blank(),
          axis.text = element_text(size=10, color='black'),
          axis.title = element_text(size=10, color='black'),
          title = element_text(size=10, color='black'),
          plot.tag = element_text(face = "bold")) +
    paper_figs_theme_no_legend
  
  return(p)
}


# get the relevant intra layers (30 / 20 / 10 / 5)
# get multilayer from filtering: 
# 30%
intra30 <- read_csv('local_output/farm_multilayer_pos_30.csv') %>% 
             select(layer=level_name, node_from=from, node_to=to, weight)
mln <- get_common_partners_mln(intra30, layers)
write_csv(mln, "local_output/mln_multi_interlayer_edges_30.csv")

# 20%
intra20 <- read_csv('local_output/farm_multilayer_pos_20.csv') %>% 
  select(layer=level_name, node_from=from, node_to=to, weight)
mln <- get_common_partners_mln(intra20, layers)
write_csv(mln, "local_output/mln_multi_interlayer_edges_20.csv")

# 10%
intra10 <- read_csv('local_output/farm_multilayer_pos_10.csv') %>% 
  select(layer=level_name, node_from=from, node_to=to, weight)
mln <- get_common_partners_mln(intra10, layers)
write_csv(mln, "local_output/mln_multi_interlayer_edges_10.csv")

# 5%
intra05 <- read_csv('local_output/farm_multilayer_pos_05.csv') %>% 
  select(layer=level_name, node_from=from, node_to=to, weight)
mln <- get_common_partners_mln(intra05, layers)
write_csv(mln, "local_output/mln_multi_interlayer_edges_05.csv")

# run the infomap
m_30 <- infomap_from_mln("local_output/mln_multi_interlayer_edges_30.csv")
m_30$m # 1

m_20 <- infomap_from_mln("local_output/mln_multi_interlayer_edges_20.csv")
m_20$m # 2

m_10 <- infomap_from_mln("local_output/mln_multi_interlayer_edges_10.csv")
m_10$m # 4

m_05 <- infomap_from_mln("local_output/mln_multi_interlayer_edges_05.csv")
m_05$m # 7


write_csv(m_30$modules %>% left_join(layers), 'local_output/multi_inter_edge_30_modules.csv')
write_csv(m_20$modules %>% left_join(layers), 'local_output/multi_inter_edge_20_modules.csv')
write_csv(m_10$modules %>% left_join(layers), 'local_output/multi_inter_edge_10_modules.csv')
write_csv(m_05$modules %>% left_join(layers), 'local_output/multi_inter_edge_05_modules.csv')

# find division - distribution across farms plot
modules_to_plot(m_30$modules) + labs(title = "30%")
modules_to_plot(m_20$modules) + labs(title = "20%")
modules_to_plot(m_10$modules) + labs(title = "10%")
modules_to_plot(m_05$modules) + labs(title = "5%")


# Run Infomap - Multi-level
ml <- run_infomap_multilayer_multilevel(multilayer_for_infomap, 
                                        two_level = F, 
                                        flow_model = 'undirected', silent = F, 
                                        trials = 200, relax = F, seed=123)
modules <- ml$modules %>% left_join(layers)
# See that the file name matches the network before writing
write_csv(modules, 'local_output/figures/multi_edge_intra_30_modules_multilevel.csv')


# run infomap of multi-interlayer-edges mln for one of the shuffled networks
# shuff - example
intra_sh <- read_csv('HPC/shuffled/shuffle_farm_r0_30_500_jac_intra/001/1_multilayer_pf_unif.csv') %>%
  filter(type=="intra") %>%  
  left_join(layers, by=c("layer_from"="layer_id")) %>% 
  select(layer=short_name, node_from, node_to, weight)
mln_sh <- get_common_partners_mln(intra_sh, layers)

all_nodes <- sort(unique(c(mln_sh$node_from, mln_sh$node_to)))
all_nodes <- tibble(node_id=1:length(all_nodes), node_name=all_nodes)

mln_sh_inds <- 
  mln_sh %>%
  left_join(layers, by = c('layer_from' = 'short_name')) %>% 
  left_join(layers, by = c('layer_to' = 'short_name')) %>% 
  select(layer_from=layer_id.x, node_from, layer_to=layer_id.y, node_to, weight)

# Run Infomap
multilayer_for_infomap <- create_multilayer_object(extended = mln_sh_inds,
                                                   nodes = all_nodes,
                                                   layers = layers)
m <- infomapecology::run_infomap_multilayer(multilayer_for_infomap, silent = F,
                                            flow_model = 'undirected',  
                                            trials = 200, relax = F, seed=111)
# it is all one big module.


