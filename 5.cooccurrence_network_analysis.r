#------ cooccurrence_network_analysis.r -------------------------------------
# This script analyses the structure of the co-occurrence network built from 
# ASV abundance in the farm level.
#------------------------------------------------------------

#------ includes ------------
library(tidyverse)
library(magrittr)
library(igraph)
library(cowplot)#?
source('functions.R')

#------ network creation -------------------------------

# Co-occurrence networks at regional, country and farm levels are all generated
# on the HPC using the core microbe data set created above.

#------ run --------------------------------
# Co-occurrence network for farm scale 30%
e_id <- 16
run_summary <- read_csv('HPC/run_summary.csv', 
                        col_names = c('exp_id','level','level_name','JOB_ID','data_file','time_stamp')) %>%
               arrange(exp_id,level)
run_summary %>% filter(exp_id==e_id)
layers <- tibble(layer_id=1:7, layer_name=c('NUDC', 'Park', 'Bian', 'Fran','Gand','Mink','Raab'),
                 short_name=c('UK1', 'UK2', 'IT1', 'IT2', 'IT3', 'FI1', 'SE1'))

# Create a multilayer network for 7 farms with intralayer edges ----
farm_multilayer <- NULL
lyrs_list <- NULL
if (e_id>4) { # this is compatible with short farm name (the new ones)
  lyrs_list <- layers$short_name
} else{ # this is compatible with full farm name (the old ones)
  lyrs_list <- layers$layer_name
}

# Build network based on p value instead of the one used so far ------
setwd(paste('HPC/exp_',e_id, sep=''))
for (l in lyrs_list){ 
  print('------------------')
  print(l)
  print('------------------')
  x <- parse_networks_from_cooc(e_id = e_id, Level = 'Farm', Level_name = l)
  farm_multilayer <- rbind(farm_multilayer,x$edge_list)
}
setwd('../../')

# save the positive edges
write_csv(farm_multilayer, 'local_output/farm_multilayer_pos_30.csv')


# Analyse the structure ------------------------- 

# General network attributes -----------------
## Number of co-occurrence links per farm -----
farm_multilayer_pos <- read_csv('local_output/farm_multilayer_pos_30.csv')
links_per_farm <- farm_multilayer_pos %>%
  group_by(level_name) %>%
  mutate(link=1) %>%
  summarise(links=sum(link))

links_per_farm %>%
  ggplot(aes(level_name, links)) +
  geom_col(fill='dark green', color='white') +
  theme_bw() +
  labs(y = 'Number of links') +
  theme(axis.text = element_text(size=10, color='black'), title = element_text(size=14), axis.text.x = element_text(angle = 90))

## Density (connectance per farm) -------
farm_net_density <- function(d){
  g <- graph_from_data_frame(d, directed = FALSE, vertices = NULL)
  out <- tibble(d=round(igraph::edge_density(g),2))
  return(out)
}

farm_density <- 
  farm_multilayer_pos %>% 
  group_by(level_name) %>% 
  group_modify(~farm_net_density(.x))

# complete the table started on exploratory_analysis.r to make paper table 1
incomplete_table <- read_csv('local_output/summary_table_pre_network.csv')
incomplete_table %>% 
  left_join(links_per_farm, by = c('Farm'='level_name')) %>%
  left_join(farm_density, by = c('Farm'='level_name')) %>%
  as_tibble() %>% 
  rename(Cows=cow_num, ASVs=ASV_Richness, `ASVs per cow`=ASV_summary,Links=links,Density=d) %>% 
  write_csv('local_output/paper_table1.csv')


# Clustering coefficient---------------------
farm_multilayer_pos <- read_csv('local_output/farm_multilayer_pos_30.csv')

CC_obs <- 
  farm_multilayer_pos %>%
  group_by(level_name) %>% 
  group_modify(~calc_CC_local(.x)) %>%
  drop_na()


# Distributions of observed CC
pdf('/Users/Geut/Dropbox/for_processing/rumen/cc_farm_boxplot.pdf', 5, 5)
CC_obs %>% 
  filter(k>=10) %>% 
  mutate(level_name=factor(level_name, levels = c("UK1","UK2","IT1","IT2","IT3","FI1",'SE1'))) %>%
  ggplot(aes(x=level_name, y=CC))+
  geom_boxplot()+
  theme_bw() +
  labs(y='Clustering coefficient', x='Farm') +
  paper_figs_theme_no_legend
dev.off()

# CC statistics in each farm
CC_obs %>%
  filter(k>=10) %>% 
  summarise(cc_mean=mean(CC),
            cc_median=median(CC),
            cc_sd=sd(CC),
            cc_min=min(CC),
            cc_max=max(CC)) %>% 
  mutate(ASV_summary=paste(cc_mean,' [',cc_min,'-',cc_max,']',sep=''))

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

## PF_J observed network:
PF_J_obs <-
  farm_multilayer_pos_final %>%
  group_by(from) %>%
  group_modify(~calculate_PF_J(.x))

PF_J_obs <- as_tibble(PF_J_obs)

# Partner Fidelity with UniFrac ---------------------------

# read tree from phylo data
# set working directory
phylo_tree <- readRDS("local_output/fitted_asvs_phylo_tree.rds")
tree <- phylo_tree$tree

## PF_U observed network:
PF_U_obs <-
  farm_multilayer_pos_final %>%
  group_by(from) %>%
  group_modify(~calculate_PF_U(.x, tree))
names(PF_U_obs) <- c("from", "PF_U", "PF_U_sd", "num_layers","UniFrac_type")
PF_U_obs %<>% filter(UniFrac_type=='d_UW')
PF_U_obs$PF_U=1-PF_U_obs$PF_U # Work with similarity instead of dissimilarity
PF_U_obs <- as_tibble(PF_U_obs)

# Plot Jaccard and UniFrac for paper ------------------------------------------
PF_J_obs <- read_csv('local_output/PF_J_pos_30_obs.csv') # change the number according to the filter level
PF_U_obs <- read_csv('local_output/PF_U_pos_30_obs.csv')

pdf('/Users/Geut/Dropbox/for_processing/rumen/cc_pf_hist.pdf', 5, 6)
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
dev.off()

bind_rows(
  PF_U_obs %>% select(from, PF=PF_U) %>% mutate(type='UniFrac')
) %>% 
  ggplot(aes(PF)) +
  geom_histogram(alpha=1, color='white', fill='lightblue')+
  labs(x='UniFrac score', y='Count')+
  # facet_grid(~type)+
  # geom_vline(xintercept = c(-1.96, 1.96), color = 'red')+
  # geom_vline(xintercept = c(-2.5, 2.5), color = 'red', linetype='dashed')+
  paper_figs_theme



