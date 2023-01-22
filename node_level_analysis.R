library(tidyverse)
library(magrittr)
library(reshape2) 
library(vegan) 
library(igraph)
library(cowplot)
library(stringr)
library(BSDA)
library(dendextend)
library(GUniFrac)
library(data.table) # Read CSVs fast

source('functions.R')


# Get data ----------------------------------------------------------------

farm_multilayer_pos_30 <- read_csv('local_output/farm_multilayer_pos_30.csv') # for 30% filter- positive co-occurrence


# Motifs ------------------------------------------------------------------
farm_multilayer_pos_30 <- read_csv('fixed_data/local_output/farm_multilayer_pos_30.csv') # for 30% filter- positive co-occurrence
farm_multilayer_neg_30 <- read_csv('fixed_data/local_output/farm_multilayer_neg_30.csv') # for 30% filter- positive co-occurrence
d <- bind_rows(
  farm_multilayer_pos_30 %>% filter(level_name=='SE1'),
  farm_multilayer_neg_30 %>% filter(level_name=='SE1'))

d %<>% mutate(edge_type=ifelse(edge_type=='pos','+','-'))

# Get all 3-node combinations


# d %<>%
#   left_join(all_nodes, by=c('from'='node_name')) %>%
#   left_join(all_nodes, by=c('to'='node_name')) %>%
#   select(-from, -to) %>%
#   select(from=node_id.x, to=node_id.y, everything())


d <- data.frame(from=c('a','a','b','b','c','a'), to=c('b','c','c','d','e','e'), 
                edge_type=c('+','+','+','+','+','-'))

G <- graph_from_data_frame(d[,-3], directed = F)
G
plot(G, edge.label=E(G)$edge_type)
tri <- names(triangles(G))
length(tri)/3
tri <- split(tri,ceiling(seq_along(tri) / 3))
tri

openTriList <- unique(do.call(c, lapply(as_ids(V(G)), function(v) {
  do.call(c, lapply(as_ids(neighbors(G, v)), function(v1) {
    v2 <- as_ids(neighbors(G, v1))
    v2 <- v2[shortest.paths(G, v, v2) == 2]
    
    if(length(v2) != 0) {
      lapply(v2, function(vv2) { c(v, v1, vv2)[order(c(v, v1, vv2))] })
    } else { list() }
  }))
})))

get_motif <- function(n){
  g <- subgraph(G, n)
  x <- E(g)$edge_type
  paste(length(which(x=='+')),'+',length(which(x=='-')),'-',sep = '')
}

tri_open <- do.call(rbind, openTriList)
tri_closed <- do.call(rbind, tri)

motif_counts_closed <- NULL
for (i in 1:nrow(tri_closed)){
  motif_counts_closed <- c(motif_counts_closed, get_motif(tri_closed[i,]))
}
table(motif_counts_closed)

motif_counts_open <- NULL
for (i in 1:nrow(tri_open)){
  motif_counts_open <- c(motif_counts_open, get_motif(tri_open[i,]))
}
table(motif_counts_open)


# For a specific node
rows_with_node <- which(tri_open=='c', arr.ind = T)[,1]


# Clustering coefficient---------------------

# test <- data.frame(from=c('a','b','c','a'), to=c('d','d','b','b'))
# g <- graph_from_data_frame(test, directed = FALSE, vertices = NULL)
# g
# plot(g)
# transitivity(g, type = 'local')
# test2=
# bind_rows(test, 
#           test %>%
#             relocate(to, from) %>%
#             rename(to=from, from=to))
# g=graph_from_data_frame(test2, directed = FALSE, vertices = NULL)
# g
# plot(g)
# transitivity(g, type = 'local')

CC_obs <- 
  farm_multilayer_pos_30 %>%
  group_by(level_name) %>% 
  group_modify(~calc_CC_local(.x)) %>%
  drop_na()

# How does CC correlate with degree?
CC_obs %>% 
  ggplot(aes(k, CC))+
  geom_point()+facet_wrap(~level_name, scales = 'free_x')

# Distributions of observed CC
CC_obs %>% 
  filter(k>=10) %>% 
  ggplot(aes(CC))+
  geom_histogram(color='white')+
  facet_wrap(~level_name, scales='free_y')+
  theme_bw() +
  labs(x='Clustering Coefficient', y='Count' ,title='Local Clustering Coefficient') +
  theme(axis.text = element_text(size=10, color='black'), axis.title = element_text(size=14, color='black'), title = element_text(size=15, color='black'))

CC_obs_mean <- 
  CC_obs %>%
  mutate(level_name=factor(level_name, levels = c("UK1","UK2","IT1","IT2","IT3","FI1",'SE1'))) %>%
  drop_na() %>%
  group_by(level_name) %>%
  summarise(CC_mean=mean(CC))


## Compare to shuffled networks -------------------------------------------

# Folder from HPC containing the 001-500 sub-folders
parent.folder <- "fixed_data/HPC/shuffled/shuffle_farm_r0_30_500/"
sub.folders <- list.dirs(parent.folder, recursive=TRUE)[-1]

# read all shuffled networks and join to one object
shuffled_networks <- NULL
for (s in sub.folders) {
  print(s)
  farm_net_shuff <- list.files(path = s , pattern = paste('_edge_list.csv',sep=""), recursive = T, full.names = T)
  farm_net_shuff <- 
    sapply(farm_net_shuff, fread, simplify=FALSE) %>% 
    bind_rows(.id = "id") %>% 
    mutate(id=which(sub.folders==s))
  shuffled_networks <- rbind(shuffled_networks,farm_net_shuff)
}
shuffled_networks <- as_tibble(shuffled_networks)
write_csv(shuffled_networks, 'fixed_data/local_output/all_shuffled_networks_farm_r0_30_shuff_500.csv')

# Or read the networks if previous steps have been concluded
shuffled_networks <- fread('fixed_data/local_output/all_shuffled_networks_farm_r0_30_shuff_500.csv')

# Calculate CC on shuffled networks
CC_shuff <- NULL
for (i in unique(shuffled_networks$id)) {
  print(i)
  CC_shuff <- 
    bind_rows(CC_shuff,
              shuffled_networks %>%
                filter(id==i) %>%
                select(from,to,level_name) %>%
                group_by(level_name) %>% 
                group_modify(~calc_CC_local(.x)) %>%
                drop_na() %>% 
                mutate(id=i, .before = level_name)
                 )
}

CC_combined <- 
  bind_rows(
    CC_obs %>% mutate(id=0, group='Obs') %>% select(id, everything()),
    CC_shuff %>% mutate(group='Shuff'))
write_csv(CC_combined, 'fixed_data/local_output/CC_pos_30_combined_r0.csv')
# Can also get CC_obs and CC_shuff from the CC_combined
CC_combined <- read_csv('local_output/CC_pos_30_combined_r0.csv')
CC_obs <- CC_combined %>% filter(group=='Obs') %>% select(-id,-group)
CC_shuff <- CC_combined %>% filter(group!='Obs') %>% select(-id,-group)

# Calculate z-scores. 
CC_z_score <- 
  inner_join(CC_obs %>% filter(k>=10),
             CC_shuff %>% filter(k>=10) %>% 
                          group_by(level_name,ASV_ID) %>% 
                          summarise(cc_shuff_mean=mean(CC), cc_shuff_sd=sd(CC), n=n())
  ) %>% 
  drop_na() %>% 
  mutate(z=(CC-cc_shuff_mean)/cc_shuff_sd)

CC_z_score %<>%
  mutate(signif=case_when(z>1.96 ~ 'above', # Obs is more than the shuffled
                          z< -1.96 ~ 'below', # Obs is lower than the shuffled
                          z<=1.96 | z>=-1.96 ~ 'not signif'
  ))

## Plots for paper --------

# pdf('~/Dropbox (BGU)/Apps/Overleaf/Rumen microbiome coocurrence/CC.pdf', 10, 6)
# pdf('fixed_data/local_output/figures/clustering_coefficient.pdf',8,8)
png(filename = 'local_output/figures/clustering_coefficient.png', width = 1300, height = 900, res = 300)
CC_obs_plot <- 
CC_obs %>% 
  filter(k>=10) %>% 
  mutate(level_name=factor(level_name, levels = c("UK1","UK2","IT1","IT2","IT3","FI1",'SE1'))) %>%
  # ggplot(aes(CC, fill=level_name))+
  # geom_histogram(color='white', alpha=0.8)+
  # facet_wrap(~level_name, scales='free_y')+
  # geom_vline(data=CC_obs_mean, aes(color=level_name, xintercept=CC_mean), size=1)+
  ggplot(aes(x=level_name, y=CC))+
  geom_boxplot()+
  theme_bw() +
  labs(y='Clustering coefficient', x='Farm', tag = "(A)") +
  theme(panel.grid=element_blank(),
        axis.text = element_text(size=10, color='black'),
        axis.title = element_text(size=10, color='black'),
        plot.tag = element_text(face = "bold"))
  paper_figs_theme_no_legend
dev.off()

# Set colors for significance categories
signif_colors <- c('blue','orange','#32a852')
names(signif_colors) <- c('not signif','below','above')
# pdf('~/Dropbox (BGU)/Apps/Overleaf/Rumen microbiome coocurrence/CC_prop_signif_r0.pdf', 10, 6)
# pdf('fixed_data/local_output/figures/CC_prop_signif_r0.pdf',10, 6)
png(filename = 'local_output/figures/clustering_coefficient.png', width = 1300, height = 900, res = 300)
CC_z_score %>%
  filter(k>=10) %>% 
  mutate(level_name=factor(level_name, levels = c("UK1","UK2","IT1","IT2","IT3","FI1",'SE1'))) %>%
  group_by(level_name,signif) %>% 
  summarise(n=n()) %>% 
  mutate(prop = n / sum(n)) %>% 
  mutate(N = sum(n)) %>% 
  mutate(signif=factor(signif, levels=c('not signif','below','above'))) %>%
  # left_join(signif_colors) %>% 
  mutate(ypos = cumsum(prop)- 0.5*prop ) %>% 
  ggplot(aes(x="", y=prop, fill=signif))+
  facet_grid(~level_name)+
  geom_bar(stat="identity", width=1) +
  scale_fill_manual(values = signif_colors)+
  # geom_text(aes(y = ypos, label = round(prop,2)), color = "white", size=3) +
  coord_polar("y", start=0)+
  paper_figs_theme_no_legend+
  theme(strip.text = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.border = element_blank())
dev.off()

CC_z_score %>% 
  filter(k>=10) %>% 
  mutate(level_name=factor(level_name, levels = c("UK1","UK2","IT1","IT2","IT3","FI1",'SE1'))) %>%
  ggplot(aes(x=level_name, y=z, fill=level_name))+
  geom_boxplot()+
  labs(y='Z-score', x='Farm')+
  geom_hline(yintercept = c(-1.96, 1.96), color = 'black', linetype='dashed')+
  paper_figs_theme


# Taxonomy in triangles ---------------------------------------------------
ASV_taxa <- read_csv('local_output/ASV_full_taxa.csv') %>% select(ASV_ID, everything(), -seq16S)

x <- farm_multilayer_pos_30 %>% filter(level_name=='UK2')
g <- graph_from_data_frame(x, directed = FALSE, vertices = NULL)
tri <- names(triangles(g))
tri <- tibble(tri_id=rep(1:(length(tri)/3),each=3), ASV_ID=tri)
max(tri$tri_id)

CC_taxonomy <- 
  tri %>% 
  left_join(ASV_taxa) %>% select(-Kingdom)

triangles_taxonomy <- 
CC_taxonomy %>% 
  # group_by(tri_id) %>% dplyr::sample_n(tri_id, size=50000) %>% 
  filter(tri_id<=50000) %>%
  gather(key=tax_level, value, -tri_id, -ASV_ID) %>% 
  group_by(tri_id, tax_level) %>% filter(all(!is.na(value))) %>% 
  summarise(num_taxa=n_distinct(value)) %>% 
  mutate(tax_level=factor(tax_level, levels = c('Kingdom','Phylum','Class','Order','Family','Genus','Species')))

triangles_taxonomy %>% 
  group_by(tax_level, num_taxa) %>% 
  summarise(n=n_distinct(tri_id)) %>% 
  ggplot(aes(fill=as.factor(num_taxa), y=n, x=tax_level)) + 
  geom_bar(position="fill", stat="identity")

  
# ggplot(triangles_taxonomy, aes(fill=num_taxa, x=tax_level, y=))+geom_bar()+
  



# Partner Fidelity with Jaccard ---------------------------
# Because this is an undirected network, not all ASVs are in the from column,
# which we use for the analysis. So duplicate the links to ensure that all ASVs
# are in the from column.
setdiff(farm_multilayer_pos_30$from,farm_multilayer_pos_30$to) # These ASVs were missing
setdiff(farm_multilayer_pos_30$to,farm_multilayer_pos_30$from) # These ASVs were missing


farm_multilayer_pos_final <- 
  bind_rows(farm_multilayer_pos_30, 
            farm_multilayer_pos_30 %>%
              relocate(to, from) %>%
              rename(to=from, from=to))

n_distinct(farm_multilayer_pos_30$from)
n_distinct(farm_multilayer_pos_final$from)
setdiff(farm_multilayer_pos_final$from,farm_multilayer_pos_30$to) # These ASVs were missing

# keep only ASVs that occur in 2 or more farms
farm_multilayer_pos_final %<>%
  group_by(from) %>%
  mutate(num_farms_from=n_distinct(level_name)) %>%
  filter(num_farms_from>=2)

## PF_J observed network ---------------------------

# # Demo - calculate PF_J for a single ASV
# mat_farm_ASV_00001 <-
#   farm_multilayer_pos_final %>%
#   filter(from=="ASV_00001") %>%
#   group_by(to) %>%
#   select(c(to,level_name)) %>%
#   mutate(present=1) %>%
#   spread(to, present, fill = 0) %>%
#   column_to_rownames("level_name")
# beta_farms_ASV_00001 <- vegdist(mat_farm_ASV_00001, "jaccard")
# beta_farms_ASV_00001 <- 1-beta_farms_ASV_00001
# PF_J <- mean(beta_farms_ASV_00001)
# PF_J_sd <- sd(beta_farms_ASV_00001)
# num_layers <- nrow(as.matrix(beta_farms_ASV_00001))
# tibble(ASV="ASV_00001", PF_J, PF_J_sd, num_layers=num_layers)

PF_J_obs <- 
  farm_multilayer_pos_final %>% 
  group_by(from) %>% 
  group_modify(~calculate_PF_J(.x))

PF_J_obs <- as_tibble(PF_J_obs)
write_csv(PF_J_obs, 'fixed_data/local_output/PF_J_pos_30_obs.csv')

## PF_J shuffled networks ---------------------------

# First, perform analysis on HPC (see Wiki for pipeline description)

# Folder from HPC containing the 001-500 sub-folders
parent.folder <- "fixed_data/HPC/shuffled/shuffle_farm_r0_30_500/"
sub.folders <- list.dirs(parent.folder, recursive=TRUE)[-1]

# read all shuffled networks
PF_J_shuff <- NULL
for (dir in sub.folders) {
  print(dir)
  # suppressMessages(shuff_fid <- read_csv(paste(dir,"/fidelity_shuff_farm_30.csv", sep="")))
  shuff_fid <- fread(paste(dir,"/fidelity_shuff_farm_30.csv", sep="")) # Faster to read with this
  PF_J_shuff <- rbind(PF_J_shuff, shuff_fid)
}
PF_J_shuff <- as_tibble(PF_J_shuff)
write_csv(PF_J_shuff, 'fixed_data/local_output/PF_J_pos_30_shuffled_r0.csv')


## START HERE FOR READING ALREADY WRITTEN RESULTS ------------

PF_J_obs <- read_csv('local_output/PF_J_pos_30_obs.csv')

# Observed PF
mean(PF_J_obs$PF_J)
sd(PF_J_obs$PF_J)
median(PF_J_obs$PF_J)
PF_J_obs %>%
  ggplot(aes(PF_J)) +
  geom_histogram(fill='purple',color='white')+
  paper_figs_theme_no_legend
# Correlation with number of farms ASV occur in
PF_J_obs %>%  
  ggplot(aes(x = num_layers, y = PF_J, group=num_layers))+
  scale_x_continuous(breaks = 1:7)+
  geom_boxplot()+
  labs(title = 'fidelity as function of the number of farms')+
  theme_bw()+theme(panel.grid.minor = element_blank())


PF_J_shuff <- read_csv('local_output/PF_J_pos_30_shuffled_r0.csv')
names(PF_J_shuff)[1:4] <- names(PF_J_obs)
mean(PF_J_shuff$PF_J)
median(PF_J_shuff$PF_J)
PF_J_shuff %>%
  ggplot(aes(PF_J)) +
  geom_histogram(fill='purple',color='white')+
  paper_figs_theme_no_legend
# Correlation with number of farms ASV occur in
PF_J_shuff %>%  
  ggplot(aes(x = num_layers, y = PF_J, group=num_layers))+
  scale_x_continuous(breaks = 1:7)+
  geom_boxplot()+
  labs(title = 'fidelity as function of the number of farms')+
  theme_bw()+theme(panel.grid.minor = element_blank())

## PF_J: Observed vs shuffled --------------------------------
PF_J_obs <- read_csv('local_output/PF_J_pos_30_obs.csv')
PF_J_shuff <- read_csv('local_output/PF_J_pos_30_shuffled_r0.csv')

PF_J_z_score <- 
  PF_J_shuff %>%
  group_by(from) %>%
  summarise(PF_J_shuff_mean=mean(PF_J), PF_J_shuff_sd=sd(PF_J)) %>% 
  inner_join(PF_J_obs) %>%
  mutate(z=(PF_J-PF_J_shuff_mean)/PF_J_shuff_sd) %>% 
  mutate(signif=case_when(z>1.96 ~ 'above', # Obs is more than the shuffled
                          z< -1.96 ~ 'below', # Obs is lower than the shuffled
                          z<=1.96 | z>=-1.96 ~ 'not signif'))

# What proportion of ASVs have a statistical significant PF_J?
PF_J_z_score %>% 
  group_by(signif) %>% 
  summarise(n=n(),prop=n/nrow(PF_J_z_score))


# Partner Fidelity with UniFrac ---------------------------


## NOTE THAT UniFrac is calculated as dissimilarity in the function
## calculate_PF_U. 

## PF_U observed network ---------------------------

# read tree from phylo data
# set working directory
phylo_tree <- readRDS("local_output/fitted_asvs_phylo_tree.rds")

PF_U_obs <- 
  farm_multilayer_pos_final %>% 
  group_by(from) %>% 
  group_modify(~calculate_PF_U(.x, phylo_tree$tree))
names(PF_U_obs) <- c("from", "PF_U", "PF_U_sd", "num_layers","UniFrac_type")
PF_U_obs %<>% filter(UniFrac_type=='d_UW')
PF_U_obs$PF_U=1-PF_U_obs$PF_U # Work with similarity instead of dissimilarity
PF_U_obs <- as_tibble(PF_U_obs)
write_csv(PF_U_obs, 'fixed_data/local_output/PF_U_pos_30_obs.csv')

## PF_U shuffled networks ---------------------------

# First, perform analysis on HPC (see Wiki for pipeline description)

# Folder from HPC containing the 001-500 sub-folders
parent.folder <- "fixed_data/HPC/shuffled/shuffle_farm_r0_30_500/"
sub.folders <- list.dirs(parent.folder, recursive=TRUE)[-1]

# read all shuffled networks
PF_U_shuff <- NULL
for (dir in sub.folders) {
  print(dir)
  # suppressMessages(shuff_fid <- read_csv(paste(dir,"/uniFrec_shuff_farm_30.csv", sep="")))
  shuff_fid <- fread(paste(dir,"/uniFrec_shuff_farm_30.csv", sep="")) # Faster to read with this
  shuff_fid <- subset(shuff_fid, unifrec_type=='d_UW')
  PF_U_shuff <- rbind(PF_U_shuff, shuff_fid)
}
names(PF_U_shuff) <- c("from", "PF_U", "PF_U_sd", "num_layers","UniFrac_type", "run" )
PF_U_shuff$PF_U=1-PF_U_shuff$PF_U # Work with similarity instead of dissimilarity
PF_U_shuff <- as_tibble(PF_U_shuff)
write_csv(PF_U_shuff, 'fixed_data/local_output/PF_U_pos_30_shuffled_r0.csv')

## START HERE FOR READING ALREADY WRITTEN RESULTS ------------

PF_U_obs <- read_csv('fixed_data/local_output/PF_U_pos_30_obs.csv')

# Observed PF
mean(PF_U_obs$PF_U)
median(PF_U_obs$PF_U)
sd(PF_U_obs$PF_U)
PF_U_obs %>%
  ggplot(aes(PF_U)) +
  geom_histogram(fill='purple',color='white')+
  paper_figs_theme_no_legend
# Correlation with number of farms ASVs occur in
PF_U_obs %>%  
  ggplot(aes(x = num_layers, y = PF_U, group=num_layers))+
  scale_x_continuous(breaks = 1:7)+
  geom_boxplot()+
  labs(title = 'fidelity as function of the number of farms')+
  theme_bw()+theme(panel.grid.minor = element_blank())


PF_U_shuff <- read_csv('local_output/PF_U_pos_30_shuffled_r0.csv')
mean(PF_U_shuff$PF_U)
sd(PF_U_shuff$PF_U)
median(PF_U_shuff$PF_U)
PF_U_shuff %>%
  ggplot(aes(PF_U)) +
  geom_histogram(fill='purple',color='white')+
  paper_figs_theme_no_legend
# Correlation with number of farms ASV occur in
PF_U_shuff %>%  
  ggplot(aes(x = num_layers, y = PF_U, group=num_layers))+
  scale_x_continuous(breaks = 1:7)+
  geom_boxplot()+
  labs(title = 'fidelity as function of the number of farms')+
  theme_bw()+theme(panel.grid.minor = element_blank())

## PF_U: Observed vs shuffled --------------------------------
PF_U_obs <- read_csv('local_output/PF_U_pos_30_obs.csv')
PF_U_shuff <- read_csv('local_output/PF_U_pos_30_shuffled_r0.csv')

PF_U_z_score <- 
  PF_U_shuff %>%
  group_by(from) %>%
  summarise(PF_U_shuff_mean=mean(PF_U), PF_U_shuff_sd=sd(PF_U)) %>% 
  inner_join(PF_U_obs) %>%
  mutate(z=(PF_U-PF_U_shuff_mean)/PF_U_shuff_sd) %>% 
  mutate(signif=case_when(z>1.96 ~ 'above', # Obs is more than the shuffled
                          z< -1.96 ~ 'below', # Obs is lower than the shuffled
                          z<=1.96 | z>=-1.96 ~ 'not signif'))

# What proportion of ASVs have a statistical significant UF?
PF_U_z_score %>% 
  group_by(signif) %>% 
  summarise(n=n(),prop=n/nrow(PF_U_z_score))


# Plot Jaccard and UniFrac for paper ------------------------------------------
PF_J_obs <- read_csv('fixed_data/local_output/PF_J_pos_30_obs.csv')
PF_U_obs <- read_csv('fixed_data/local_output/PF_U_pos_30_obs.csv')

PF_score_plot <- 
bind_rows(
  PF_J_obs %>% select(from, PF=PF_J) %>% mutate(type='Jaccard'),
  PF_U_obs %>% select(from, PF=PF_U) %>% mutate(type='UniFrac')
) %>% 
  ggplot(aes(PF, fill=type)) +
<<<<<<< Updated upstream
  geom_histogram(alpha=1, color='white')+
  labs(x='Partner fidelity score', y='Count', tag = "(B)")+
  # facet_grid(~type)+
  # geom_vline(xintercept = c(-1.96, 1.96), color = 'red')+
  # geom_vline(xintercept = c(-2.5, 2.5), color = 'red', linetype='dashed')+
  paper_figs_theme_no_legend+
  theme(panel.grid=element_blank(),
        axis.text = element_text(size=10, color='black'),
        axis.title = element_text(size=10, color='black'),
        plot.tag = element_text(face = "bold"))+
  guides(fill=guide_legend(title=NULL))
=======
  geom_histogram(alpha=0.6, color='white', position="identity")+
  labs(x='Partner fidelity score', y='Count')+
  # facet_grid(~type)+
  # geom_vline(xintercept = c(-1.96, 1.96), color = 'red')+
  # geom_vline(xintercept = c(-2.5, 2.5), color = 'red', linetype='dashed')+
  paper_figs_theme +
  theme(panel.grid=element_blank(),
        axis.text = element_text(size=10, color='black'),
        axis.title = element_text(size=10, color='black'),
        legend.position = c(0.6, 70))
>>>>>>> Stashed changes
  
PF_score_plot

PF_z_score_plot <- 
  bind_rows(
    PF_U_z_score %>% select(from, z) %>% mutate(type='UniFrac'),
    PF_J_z_score %>% select(from, z) %>% mutate(type='Jaccard')
  ) %>% 
  ggplot(aes(z, fill=type)) +
  geom_histogram(alpha=1, color='white')+
  labs(x='z-score', y='Count', fill='PF index')+
  geom_vline(xintercept = c(-1.96, 1.96), color = 'black', linetype='dashed')+
  # geom_vline(xintercept = c(-2.5, 2.5), color = 'red', linetype='dashed')+
  paper_figs_theme
PF_z_score_plot

# pdf('~/Dropbox (BGU)/Apps/Overleaf/Rumen microbiome coocurrence/partner_fidelity.pdf', 10, 6)
# cowplot::plot_grid(PF_score_plot+theme(legend.position = 'none'),PF_z_score_plot+theme(legend.position = c(0.8,0.9)))
# dev.off()

# pdf('~/Dropbox (BGU)/Apps/Overleaf/Rumen microbiome coocurrence/partner_fidelity.pdf', 10, 6)
# pdf('fixed_data/local_output/figures/partner_fidelity.pdf',10,6)
png(filename = 'local_output/figures/partner_fidelity.png', width = 1300, height = 900, res = 300)
PF_score_plot
dev.off()

# pdf('~/Dropbox (BGU)/Apps/Overleaf/Rumen microbiome coocurrence/partner_fidelity_signif.pdf', 10, 6)
# pdf('fixed_data/local_output/figures/partner_fidelity_signif.pdf',10,6)
png(filename = 'local_output/figures/partner_fidelity_signif.png', width = 1300, height = 900, res = 300)
bind_rows(
  PF_U_z_score %>% select(from, z,signif) %>% mutate(type='UniFrac'),
  PF_J_z_score %>% select(from, z,signif) %>% mutate(type='Jaccard')
) %>% 
  group_by(type,signif) %>% 
  summarise(n=n()) %>% 
  mutate(prop = n / sum(n)) %>% 
  mutate(N = sum(n)) %>% 
  mutate(signif=factor(signif, levels=c('not signif','below','above'))) %>% 
  mutate(ypos = cumsum(prop)- 0.5*prop ) %>% 
  ggplot(aes(x="", y=prop, fill=signif))+
  facet_wrap(~type)+
  geom_bar(stat="identity", width=1) +
  scale_fill_manual(values = c('blue','orange','#32a852'))+
  # geom_text(aes(y = ypos, label = round(prop,2)), color = "white", size=3) +
  coord_polar("y", start=0)+
  paper_figs_theme_no_legend+theme_void()
dev.off()

pdf('local_output/figures/clustering_coefficient_pf.pdf',10,6)
pdf('/Users/Geut/Dropbox/for_processing/rumen/clustering_coefficient_pf_no_pies.pdf',10,4)
# png(filename = 'local_output/figures/clustering_coefficient_pf.png', width = 1500, height = 500, res = 300)
plot_grid(CC_obs_plot,PF_score_plot, labels = c('(A)','(B)'))
dev.off()
