# analysis on different filters (20%,10%,5%)
# Includes ------------
library(tidyverse)
library(magrittr)
library(cowplot)
library(reshape2) 
library(vegan)
library(igraph)
library(GUniFrac)
library(infomapecology)
check_infomap()
source('functions.R')
# Load data --------------------------------
# change the number of the file to the wanted filtration level
ASV_Core <- read_csv('local_output/core_ASV_fixed_20.csv') %>% 
  mutate(Farm=factor(Farm, levels = c("UK1","UK2","IT1","IT2","IT3","FI1",'SE1')))
e_id <- 17

## @knitr END
# Supplementary
##1 Richness per cow  --------------------------------
pdf('local_output/exploratory_analysis_20.pdf',8,8)
# plt_richness_per_cow <- 
  ASV_Core %>%
  group_by(Cow_Code) %>%
  summarise(Richness=n_distinct(ASV_ID)) %>%
  arrange(desc(Richness)) %>% 
  ggplot(aes(Richness))+
  geom_histogram(fill='dark green', color='white')+
  labs(x='ASVs per cow', y='Count') +
  html_figs_theme_no_legend

##2 Richness per farm  --------------------------------
richness_per_farm <- ASV_Core %>%
  group_by(Farm) %>%
  summarise(ASV_Richness=n_distinct(ASV_ID)) 

richness_per_farm %>%
  arrange(desc(ASV_Richness)) %>% 
  ggplot(aes(Farm,ASV_Richness))+
  geom_col(fill='dark red', color='white')+
  theme_classic() +
  labs(y='Number of ASVs') +
  html_figs_theme_no_legend



##3 Richness per cow per farm   --------------------------------
Richness_per_cow_farm <- 
  ASV_Core %>%
  group_by(Farm,Cow_Code) %>%
  summarise(S_cow=n_distinct(ASV_ID)) %>% 
  group_by(Farm) %>% 
  summarise(S_mean=round(mean(S_cow),1),
            S_median=round(median(S_cow),1),
            S_sd=round(sd(S_cow),1),
            S_min=min(S_cow),
            S_max=max(S_cow)) %>% 
  mutate(ASV_summary=paste(S_mean,' [',S_min,'-',S_max,']',sep=''))


ASV_Core %>%
  group_by(Country,Farm,Cow_Code) %>%
  summarise(richness=n_distinct(ASV_ID)) %>%
  arrange(desc(richness)) %>% 
  ggplot(aes(x=Farm, y=richness, Country))+
  geom_boxplot(aes(color=Country))+
  theme_bw() +
  labs(y = 'ASV richness per cow per farm') +
  html_figs_theme_no_legend

## @knitr END

#4 ASV Beta diversity between farms ----------------
ASV_occurrence_farm <- 
  ASV_Core %>%
  select(-c(Cow_Code)) %>%
  distinct(Farm, ASV_ID) %>%
  mutate(present=1) %>% 
  spread(ASV_ID, present, fill = 0) %>%
  column_to_rownames("Farm")

## Using Jaccard ----
beta_diver_farms <- as.matrix(1-vegdist(ASV_occurrence_farm, "jaccard"))
diag(beta_diver_farms) <- 1
write.csv(beta_diver_farms, "local_output/jaccard_beta_diversity_30", row.names = TRUE)

# Heatmap
beta_diver_farms[lower.tri(beta_diver_farms, diag = F)] <- NA
beta_diver_farms_m <- reshape2::melt(beta_diver_farms) 

plt_J <- 
  ggplot(beta_diver_farms_m, aes(x = Var1, y = Var2, fill = value, label=round(value,2))) +
  geom_tile() + 
  geom_text(color = "black", size = 4)+
  scale_y_discrete(limits=rev)+
  scale_fill_gradient(high = "blue", low = "light blue", na.value = 'white') +
  html_figs_theme_no_legend+theme(axis.title = element_blank())

## Using UniFrac ----
# tree <- readRDS("local_output/rooted_phylo_tree.rds") # only for 5%
phylo_tree <- readRDS("local_output/fitted_asvs_phylo_tree.rds")
tree <- phylo_tree$tree
# prune the tree
included_asvs <- unique(ASV_Core$ASV_ID)
unincluded <- tree$tip.label[!tree$tip.label %in% included_asvs]
pruned <- dendextend::prune(tree, unincluded)
unifracs <- GUniFrac(ASV_occurrence_farm, pruned, alpha=c(0, 0.5, 1))$unifracs
d_UW <- 1-(unifracs[, , "d_UW"])
# Heatmap 
d_UW[lower.tri(d_UW)] <- NA
beta_diver_farms_UF <- reshape2::melt(d_UW)
plt_U <- 
  ggplot(beta_diver_farms_UF, aes(x = Var1, y = Var2, fill = value, label=round(value,2))) +
  geom_tile() + 
  geom_text(color = "black", size = 4)+
  scale_y_discrete(limits=rev)+
  scale_fill_gradient(high = "#ff8c00", low = "#fff494", na.value = 'white') +
  html_figs_theme_no_legend+theme(axis.title = element_blank())

# combine the two plots of beta diversity
plot_grid(plt_J,plt_U, labels = c('(A)','(B)'))

## @knitr END


## @knitr betadiv

#5 Taxonomic analysis ------------------
ASV_taxa <- read_csv('local_output/ASV_full_taxa.csv') %>% 
  select(ASV_ID, everything(), -seq16S)
n_distinct(ASV_Core$ASV_ID)
## Phylum-level composition --------------------------------------------

ASV_Core %>% 
  left_join(ASV_taxa) %>%
  group_by(Farm, Phylum) %>% 
  summarise(ASV_num=n_distinct(ASV_ID)) %>%
  drop_na() %>% 
  mutate(relative_richness=ASV_num/sum(ASV_num)) %>% 
  mutate(ypos = cumsum(relative_richness)- 0.5*relative_richness) %>% 
  ggplot(aes(x="", y=relative_richness, fill=Phylum))+
  facet_wrap(~Farm)+
  geom_bar(stat="identity", width=1) +
  # scale_fill_manual(values = signif_colors)+
  # geom_text(aes(y = ypos, label = round(prop,2)), color = "white", size=3) +
  coord_polar("y", start=0)+
  paper_figs_theme+
  theme(axis.text = element_blank(),
        axis.title = element_blank())

#6 Number of cows in which microbes occur -----------
ASV_Core %>%
  group_by(ASV_ID) %>%
  summarise(habitats=n_distinct(Cow_Code)) %>%
  ggplot(aes(habitats))+
  geom_histogram(fill='brown', color='white')+
  labs(x='Number of cows', y='Count')+
  html_figs_theme_no_legend

#7 Number of farms in which microbes occur -----------
ASV_Core %>%
  group_by(ASV_ID) %>%
  summarise(habitats=n_distinct(Farm)) %>%
  arrange(desc(habitats)) %>% 
  group_by(habitats) %>% 
  summarise(n=n_distinct(ASV_ID)) %>% 
  ggplot(aes(habitats, n))+
  geom_col(fill='dark green')+
  labs(x='Number of farms', y='Count')+
  scale_x_continuous(breaks = seq(0,7,1))+
  paper_figs_theme_no_legend


dev.off() # close the handle that collects all the network exploratory analysis


#------ Modularity per filtration level ()-------------------------------
ASV_Core <- read_csv('local_output/core_ASV_fixed_05.csv') %>% 
  mutate(Farm=factor(Farm, levels = c("UK1","UK2","IT1","IT2","IT3","FI1",'SE1')))
e_id <- 19 #change the number of the e_id to the wanted experiment
percent_str <- "05"

run_summary <- read_csv('HPC/run_summary.csv', 
                        col_names = c('exp_id','level','level_name','JOB_ID','data_file','time_stamp')) %>%
                arrange(exp_id,level)
run_summary %>% filter(exp_id==e_id)
layers <- tibble(layer_id=1:7, layer_name=c('NUDC', 'Park', 'Bian', 'Fran','Gand','Mink','Raab'),
                 short_name=c('UK1', 'UK2', 'IT1', 'IT2', 'IT3', 'FI1', 'SE1'))

## run infomap--------------------------------
# Create a multilayer network for 7 farms with intralayer edges ----
farm_multilayer <- NULL
setwd(paste('HPC/exp_',e_id, sep=''))
lyrs_list <- layers$short_name
if (e_id>4) { # this is compatible with short farm name (the new ones)
  lyrs_list <- layers$short_name
} else{ # this is compatible with full farm name (the old ones)
  lyrs_list <- layers$layer_name
}
# build each layer separately first
for (l in lyrs_list){ 
  print('------------------')
  print(l)
  print('------------------')
  x <- parse_networks(e_id = e_id, Level = 'Farm', Level_name = l)
  farm_multilayer <- rbind(farm_multilayer,x$edge_list)
  assign(paste('net_farm',l,sep = '_'), x)
}
setwd('../../')

# For the analysis separate the positive and negative networks
farm_multilayer_pos <- farm_multilayer %>% filter(edge_type=='pos')

# for CC:
write_csv(farm_multilayer_pos, 
          paste('local_output/core_range/farm_multilayer_pos_', percent_str,'.csv', sep = ''))

all_nodes <- sort(unique(c(farm_multilayer_pos$from, farm_multilayer_pos$to)))
all_nodes <- tibble(node_id=1:length(all_nodes), node_name=all_nodes)

intra <- 
  farm_multilayer_pos %>% 
  select(layer_from=level_name, node_from=from, layer_to=level_name, node_to=to, weight) %>% 
  left_join(layers, by = c('layer_from' = 'short_name')) %>% 
  left_join(layers, by = c('layer_to' = 'short_name')) %>% 
  left_join(all_nodes, by = c('node_from' = 'node_name')) %>% 
  left_join(all_nodes, by = c('node_to' = 'node_name')) %>% 
  select(layer_from=layer_id.x, node_from=node_id.x, layer_to=layer_id.y, node_to=node_id.y, weight)

## Interlayer links with UniFrac -------------------------------------------
# Because this is an undirected network, not all ASVs are in the from column,
# which we use for the analysis. So duplicate the links to ensure that all ASVs
# are in the from column.
setdiff(farm_multilayer_pos$from,farm_multilayer_pos$to) # These ASVs were missing
setdiff(farm_multilayer_pos$to,farm_multilayer_pos$from) # These ASVs were missing

farm_multilayer_pos_final <- 
  bind_rows(farm_multilayer_pos, 
            farm_multilayer_pos %>%
              relocate(to, from) %>%
              rename(to=from, from=to))

n_distinct(farm_multilayer_pos$from)
n_distinct(farm_multilayer_pos_final$from)
setdiff(farm_multilayer_pos_final$from,farm_multilayer_pos$to) # These ASVs were missing

# keep only ASVs that occur in 2 or more farms
farm_multilayer_pos_final %<>%
  group_by(from) %>%
  mutate(num_farms_from=n_distinct(level_name)) %>%
  filter(num_farms_from>=2)
phylo_tree <- readRDS("local_output/rooted_phylo_tree.rds")
# a for loop that calculates all the interlayer edges based on unifrac
inter_PF_U <- NULL
for (i in unique(farm_multilayer_pos_final$from)) {
  print(i)
  ASV_net <- farm_multilayer_pos_final %>%
    filter(from==i)
  tree <- phylo_tree
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
  inter_PF_U <- rbind(inter_PF_U,inter_fid_unif)
}

inter <- 
  inter_PF_U %>% 
  ungroup() %>% 
  left_join(layers, by = c('layer_from' = 'short_name')) %>% 
  left_join(layers, by = c('layer_to' = 'short_name')) %>% 
  left_join(all_nodes, by = c('ASV_ID' = 'node_name')) %>% 
  select(layer_from=layer_id.x, node_from=node_id, layer_to=layer_id.y, node_to=node_id, weight)

multilayer_unif <- rbind(intra %>% mutate(type='intra'),
                         inter %>% mutate(type='inter'))
table(multilayer_unif$type)

write_csv(multilayer_unif, 
          paste('local_output/core_range/multilayer_net_', percent_str,'.csv', sep = ''))

# open pdf handle
pdf(paste('local_output/core_range/modularity_analysis_', percent_str,'.pdf', sep = ''),8,8)

#1 Edge weight distributions
ggplot(multilayer_unif, aes(weight, fill=type))+
  geom_density(alpha=0.5)+
  labs(x='Edge weight', y='Density', title='Edge weight distributions')+
  scale_fill_manual(values = c('blue','orange'))+
  theme_bw()+
  theme(panel.grid=element_blank(),
        axis.text = element_text(size=22, color='black'),
        axis.title = element_text(size=22, color='black'),
        legend.position = c(0.9,0.9))

## Run Infomap ------------------------------------------------------
multilayer_unif <- read_csv(paste('local_output/core_range/multilayer_net_', percent_str,'.csv', sep = ''))
net <- multilayer_unif[,1:5]

# Run Infomap
multilayer_for_infomap <- create_multilayer_object(extended = net, nodes = all_nodes, layers = layers)
m <- infomapecology::run_infomap_multilayer(multilayer_for_infomap, silent = F,
                                            flow_model = 'undirected',  
                                            trials = 200, relax = F, seed=123)

modules <- m$modules %>% left_join(layers)

write_csv(modules, 
          paste('local_output/core_range/multilayer_modules_', percent_str,'.csv', sep = ''))
modules %<>% rename(level1=module)

## Analyze observed modularity results -------------------------------
modules <- read_csv(paste('local_output/core_range/multilayer_modules_', percent_str,'.csv', sep = ''))

# no threshold
modules %>%
  mutate(short_name=factor(short_name, levels = c("UK1","UK2","IT1","IT2","IT3","FI1",'SE1'))) %>%
  group_by(short_name) %>%
  mutate(nodes_in_layers=n_distinct(node_id)) %>%
  group_by(short_name,module) %>%
  mutate(nodes_in_modules=n_distinct(node_id)) %>%
  mutate(nodes_percent=nodes_in_modules/nodes_in_layers) %>%
  distinct(short_name, module, nodes_percent) %>% 
  arrange(module, short_name) %>%
  # Only include modules that contain at least 3% of the ASVs in the layer
  # Plot
  ggplot(aes(x = module, y = short_name, fill=nodes_percent))+
  geom_tile(color='white')+
  scale_x_continuous(breaks = seq(1, max(modules$module), 1))+
  scale_fill_viridis_c(limits = c(0, 1))+
  labs(x='Module ID', y='')+ ggtitle('Observed Modules')+
  paper_figs_theme

# with threshold
modules %>%
  mutate(short_name=factor(short_name, levels = c("UK1","UK2","IT1","IT2","IT3","FI1",'SE1'))) %>%
  group_by(short_name) %>%
  mutate(nodes_in_layers=n_distinct(node_id)) %>%
  group_by(short_name,module) %>%
  mutate(nodes_in_modules=n_distinct(node_id)) %>%
  mutate(nodes_percent=nodes_in_modules/nodes_in_layers) %>%
  distinct(short_name, module, nodes_percent) %>% 
  arrange(module, short_name) %>%
  # Only include modules that contain at least 3% of the ASVs in the layer
  filter(nodes_percent>=0.03) %>%
  # Plot
  ggplot(aes(x = module, y = short_name, fill=nodes_percent))+
  geom_tile(color='white')+
  scale_x_continuous(breaks = seq(1, max(modules$module), 1))+
  scale_fill_viridis_c(limits = c(0, 1))+
  labs(x='Module ID', y='')+ ggtitle('Observed Modules')+
  paper_figs_theme
dev.off() # this closes the handle that prints modularity plots

# For paper -----

# Clustering coefficient---------------------

# Distributions of observed CC - only the 30% threshold:
farm_multilayer_pos <- read_csv('local_output/farm_multilayer_pos_30.csv')

CC_obs <- 
  farm_multilayer_pos %>%
  group_by(level_name) %>% 
  group_modify(~calc_CC_local(.x)) %>%
  drop_na()

CC_obs %>% 
  filter(k>=10) %>% 
  mutate(level_name=factor(level_name, levels = c("UK1","UK2","IT1","IT2","IT3","FI1",'SE1'))) %>%
  ggplot(aes(x=level_name, y=CC))+
  geom_boxplot()+
  theme_bw() +
  labs(y='Clustering coefficient', x='Farm') +
  paper_figs_theme_no_legend


# Distributions of observed CC - Comparison between the thresholds:
farm_multilayer_pos_30 <-read_csv('local_output/farm_multilayer_pos_30.csv')
farm_multilayer_pos_20 <-read_csv('local_output/farm_multilayer_pos_20.csv') 
farm_multilayer_pos_10 <-read_csv('local_output/farm_multilayer_pos_10.csv') 
farm_multilayer_pos_05 <-read_csv('local_output/farm_multilayer_pos_05.csv')

CC_obs_30 <- farm_multilayer_pos_30 %>%
             group_by(level_name) %>% 
             group_modify(~calc_CC_local(.x)) %>%
             drop_na() %>%
             mutate(Threshold='30')

CC_obs_20 <- farm_multilayer_pos_20 %>%
             group_by(level_name) %>% 
             group_modify(~calc_CC_local(.x)) %>%
             drop_na() %>%
             mutate(Threshold='20')

CC_obs_10 <- farm_multilayer_pos_10 %>%
             group_by(level_name) %>% 
             group_modify(~calc_CC_local(.x)) %>%
             drop_na() %>%
             mutate(Threshold='10')

CC_obs_05 <- farm_multilayer_pos_05 %>%
             group_by(level_name) %>% 
             group_modify(~calc_CC_local(.x)) %>%
             drop_na() %>%
             mutate(Threshold='05')

CC_obs_all <- bind_rows(CC_obs_30, CC_obs_20, CC_obs_10, CC_obs_05)

# Distributions of observed CC
pdf('local_output/figures/clustering_coefficient_thresholds.pdf',10,6)
CC_obs_all %>% 
  filter(k>=10) %>% 
  mutate(level_name=factor(level_name, levels = c("UK1","UK2","IT1","IT2","IT3","FI1",'SE1'))) %>%
  ggplot(aes(x=level_name, y=CC))+
  geom_boxplot(aes(color=Threshold))+
  theme_bw() +
  labs(y='Clustering coefficient', x='Farm') +
  paper_figs_theme
dev.off()

# Partner Fidelity with Jaccard ---------------------------
# Because this is an undirected network, not all ASVs are in the from column,
# which we use for the analysis. So duplicate the links to ensure that all ASVs
# are in the from column.
setdiff(farm_multilayer_pos$from,farm_multilayer_pos$to) # These ASVs were missing
setdiff(farm_multilayer_pos$to,farm_multilayer_pos$from) # These ASVs were missing

farm_multilayer_pos_final <- 
  bind_rows(farm_multilayer_pos, 
            farm_multilayer_pos %>%
              relocate(to, from) %>%
              rename(to=from, from=to))

n_distinct(farm_multilayer_pos$from)
n_distinct(farm_multilayer_pos_final$from)
setdiff(farm_multilayer_pos_final$from,farm_multilayer_pos$to) # These ASVs were missing

# keep only ASVs that occur in 2 or more farms
farm_multilayer_pos_final %<>%
  group_by(from) %>%
  mutate(num_farms_from=n_distinct(level_name)) %>%
  filter(num_farms_from>=2)

## PF_J observed network ---------------------------
PF_J_obs <-
  farm_multilayer_pos_final %>%
  group_by(from) %>%
  group_modify(~calculate_PF_J(.x))

PF_J_obs <- as_tibble(PF_J_obs)

# Partner Fidelity with UniFrac ---------------------------

## PF_U observed network ---------------------------

# read tree from phylo data
# set working directory
phylo_tree <- readRDS("local_output/fitted_asvs_phylo_tree.rds")
tree <- phylo_tree$tree
# tree <- readRDS("local_output/rooted_phylo_tree.rds") # only for 5%

PF_U_obs <-
  farm_multilayer_pos_final %>%
  group_by(from) %>%
  group_modify(~calculate_PF_U(.x, tree))
names(PF_U_obs) <- c("from", "PF_U", "PF_U_sd", "num_layers","UniFrac_type")
PF_U_obs %<>% filter(UniFrac_type=='d_UW')
PF_U_obs$PF_U=1-PF_U_obs$PF_U # Work with similarity instead of dissimilarity
PF_U_obs <- as_tibble(PF_U_obs)

# Plot Jaccard and UniFrac for paper ------------------------------------------
PF_score_plot <- 
  bind_rows(
    PF_J_obs %>% select(from, PF=PF_J) %>% mutate(type='Jaccard'),
    PF_U_obs %>% select(from, PF=PF_U) %>% mutate(type='UniFrac')
  ) %>% 
  ggplot(aes(PF, fill=type)) +
  geom_histogram(alpha=1, color='white')+
  labs(x='Partner fidelity score', y='Count')+
  paper_figs_theme
PF_score_plot




