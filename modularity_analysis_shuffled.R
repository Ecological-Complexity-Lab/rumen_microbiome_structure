# includes -----------
library(tidyverse)
library(dplyr)
library(magrittr)
library(igraph)
library(reshape2)
library(cowplot)#?
library(vegan)
source('functions.R')

# shuffled microbes within farms----
# read one modularity shuffled network
farm_modules <- list.files(path = "HPC/shuffled/shuffle_farm_r0_30_500_jac_intra/001", 
                           pattern = paste('_farm_modules.csv',sep=""), full.names = T)
farm_modules_f <- sapply(farm_modules, read_csv, simplify=FALSE) %>%
  bind_rows(.id = "id") %>% select(-id) %>%
  mutate(id='1')

# Folder containing sub-folders
parent.folder <- "HPC/shuffled/shuffle_farm_r0_30_500_jac_intra"

# Sub-folders
sub.folders <- list.dirs(parent.folder, recursive=TRUE)[-1]

# read all 
farm_modules_shuffled <- NULL
for (script in sub.folders) {
  print(script)
  farm_modules_shuff <- list.files(path = script , pattern = paste('_farm_modules_pf_unif.csv',sep=""), recursive = T,full.names = T)
  farm_modules_shuff <- sapply(farm_modules_shuff, read_csv, simplify=FALSE) %>% 
    bind_rows(.id = "id") 
  farm_modules_shuff_f <- farm_modules_shuff %>% 
    mutate(id= (str_split_fixed(farm_modules_shuff$id[1], pattern = '/', n = 10)[4]))
  farm_modules_shuffled <- rbind(farm_modules_shuffled,farm_modules_shuff_f)
}

write_csv(farm_modules_shuffled, 'local_output/farm_modules_shuffled_100.csv')

# number of modules in each permutation
pdf('local_output/figures/modules_per_perm_shuff_within#2.pdf',10,6)
farm_modules_shuffled %>%
  group_by(id) %>%
  summarise(module_sum=n_distinct(module)) %>%
  ggplot(aes(id,module_sum)) +
  geom_col() +
  theme_bw() +
  labs(x='', y='', title='Number of modules in each permutation (shuffled within farms)')+
  theme(axis.text.x = element_text(angle = 90, size=6, color='black'))
dev.off()

number_of_modules <- farm_modules_shuffled %>%
  group_by(id) %>%
  summarise(module_sum=n_distinct(module))

# comparison between observed and shuffled
png(filename = 'output/figures/observed_shuffled_within_hist_curveball.png', width = 1600, height = 900, res = 300)
number_of_modules %>%
  ggplot(aes(module_sum)) +
  geom_histogram(color='white') +
  theme_bw() +
  geom_vline(aes(xintercept=11), color= 'red') + 
  geom_text(aes(x=10.5, label='Observed', y=10, angle=90)) +
  theme(axis.text = element_text(size=14, color='black'),
        axis.title = element_text(size = 14, color='black'))+
  labs(x='Number of modules' ,y='count',title = 'Observed vs Shuffled (within farms)')
dev.off()

n_distinct(farm_modules_shuffled$module)

# number of layers in modules
pdf('/Users/dafnaar/GitHub/microbiome_structure_v2/output/figures/farms_in_modules_within_hist#2.pdf',10,6)
farm_modules_shuffled %>%
  group_by(id, module) %>%
  summarise(number_of_farms=n_distinct(layer_name)) %>%
  ggplot(aes(number_of_farms)) +
  geom_histogram(color='white')+
  labs(x='Number of farms', title='Number of farms in modules (shuffled within farms)')+
  theme(axis.text.x = element_text(size=10, color='black'))+
  theme_bw()
dev.off()

# comparison between the L value in observed and shuffled
farm_modulation_summary <- read_csv('HPC/shuffled/shuffle_farm_r0_30_500_jac_intra/farm_modulation_summary_pf_unif.csv', col_names = FALSE)
png(filename = 'output/figures/observed_shuffled_L_value_curveball.png', width = 1600, height = 900, res = 300)
farm_modulation_summary %>% 
  select(L=X4) %>%
  ggplot(aes(L)) +
  geom_histogram(color='white') +
  theme_bw() +
  geom_vline(aes(xintercept=9.10248), color= 'red') + 
  #geom_text(aes(x=10.5, label='Observed', y=10, angle=90)) +
  theme(axis.text = element_text(size=15, color='black'),
        axis.title = element_text(size = 15, color='black'),
        title = element_text(size = 15, color='black'))+
  labs(y='count',title = 'Observed vs Shuffled- L value')
dev.off()

# nodes in modules
png(filename = 'output/figures/nodes_in_modules_shuffled_curveball.png', width = 1600, height = 900, res = 300)
farm_modules_shuffled %>%
  group_by(module) %>%
  summarise(physical_nodes=n_distinct(node_name)) %>%
  arrange(desc(physical_nodes)) %>%
  ggplot(aes(physical_nodes)) +
  geom_histogram(color='white') +
  theme_bw()+
  labs(x='Number of physical nodes', title='Nodes in modules (shuffled networks)')+
  theme(axis.text = element_text(size=15, color='black'), axis.title = element_text(size=15, color='black'), title = element_text(size = 16))
dev.off()

# shuffled microbes between farms----
# read one modularity shuffled network
farm_modules <- list.files(path = "HPC/shuffled/shuffle_farm_r0_30_500_jac_intra/001", 
                           pattern = paste('_farm_modules_pf_unif.csv',sep=""), full.names = T)
farm_modules_f <- sapply(farm_modules, read_csv, simplify=FALSE) %>%
  bind_rows(.id = "id") %>% select(-id) %>%
  mutate(id='1')

# Folder containing sub-folders
parent.folder <- "HPC/shuffled/shuffle_farm_r0_30_500_jac_intra"

# Sub-folders
sub.folders <- list.dirs(parent.folder, recursive=TRUE)[-1]

# read all 
farm_modules_shuffall <- NULL
for (script in sub.folders) {
  print(script)
  farm_modules_shuff <- list.files(path = script , pattern = paste('_farm_modules_pf_unif.csv',sep=""), recursive = T,full.names = T)
  farm_modules_shuff <- sapply(farm_modules_shuff, read_csv, simplify=FALSE) %>% 
    bind_rows(.id = "id") 
  farm_modules_shuffall_f <- farm_modules_shuff %>% 
    mutate(id= (str_split_fixed(farm_modules_shuff$id[1], pattern = '/', n = 10)[4]))
  farm_modules_shuffall <- rbind(farm_modules_shuffall,farm_modules_shuffall_f)
}

# number of modules in each permutation
pdf('/Users/dafnaar/GitHub/microbiome_structure_v2/output/figures/modules_per_perm.pdf',10,6)
farm_modules_shuffall %>%
  group_by(id) %>%
  summarise(module_sum=n_distinct(module)) %>%
  ggplot(aes(id,module_sum)) +
  geom_col() +
  labs(x='', y='', title='Number of modules in each permutation')+
  theme(axis.text.x = element_text(angle = 90, size=6, color='black'))
dev.off()

# number of layers in modules
pdf('/Users/dafnaar/GitHub/microbiome_structure_v2/output/figures/farms_in_modules_between_hist.pdf',10,6)
farm_modules_shuffall %>%
  group_by(id, module) %>%
  summarise(number_of_farms=n_distinct(layer_name)) %>%
  ggplot(aes(number_of_farms)) +
  geom_histogram(color='white')+
  labs(x='Number of farms', title='Number of farms in modules (shuffled between farms)')+
  theme(axis.text.x = element_text(size=10, color='black'))+
  theme_bw()
dev.off()

# number of physical nodes in modules
pdf('output/figures/physical_nodes_modules_between_hist.pdf',10,6)
farm_modules_shuffall %>%
  group_by(id, module) %>%
  summarise(physical_nodes=n_distinct(node_name)) %>%
  arrange(desc(physical_nodes)) %>%
  ggplot(aes(physical_nodes)) +
  geom_histogram(color='white')+ scale_y_log10() +
  theme_bw()+
  labs(x='Number of physical nodes', title='Number of physical nodes in modules (shuffled between)')+
  theme(axis.text = element_text(size=12, color='black'), axis.title = element_text(size=15, color='black'), title = element_text(size = 16))
dev.off()

# modules in each permutation
number_of_modules <- farm_modules_shuffall %>%
  group_by(id) %>%
  summarise(module_sum=n_distinct(module))

png(filename = 'output/figures/shuffled_all_hist.png', width = 1600, height = 900, res = 300)
number_of_modules %>%
  ggplot(aes(module_sum)) +
  geom_histogram() +
  theme_bw() +
  labs(y='count',title = 'Modules in each permutation')
dev.off()

# comparison between observed and shuffled
png(filename = 'output/figures/observed_shuffled_all_hist.png', width = 1600, height = 900, res = 300)
number_of_modules %>%
  ggplot(aes(module_sum)) +
  geom_histogram(color='white',binwidth = 9) +
  theme_bw() +
  scale_x_continuous(breaks = seq(0,1000,100)) +
  geom_vline(aes(xintercept=11), color= 'red') + 
  geom_text(aes(x=100, label='Observed', y=20)) +
  theme(axis.text.x = element_text(size=10, color='black'),
        axis.text.y = element_text(size = 10, color='black'))+
  labs(y='count',title = 'Observed vs Shuffled')
dev.off()

n_distinct(farm_modules_shuffled$module)

# looking at the biggest module----
# observed
farm_modules_pos <- read_csv('local_output/farm_modules_pos_30_U.csv')
observed_biggest_module <- farm_modules_pos %>%
  group_by(module) %>%
  mutate(layers_in_modules=n_distinct(layer_name)) %>%
  arrange(desc(layers_in_modules)) %>%
  filter(layers_in_modules==max(layers_in_modules)) %>%
  mutate(nodes_in_layers=n_distinct(node_name)) %>%
  arrange(desc(nodes_in_layers)) %>%
  filter(module == 1) %>%
  mutate(id=0) %>%
  select(c(id,module, node_name, layer_name=short_name, layers_in_modules, nodes_in_layers)) %>%
  transform(id = as.character(id)) %>%
  tibble()
# shuffled within
shuffled_biggest_module <- farm_modules_shuffled %>%
  select(-c(node_id,layer_id, flow)) %>%
  group_by(id,module) %>%
  mutate(layers_in_modules=n_distinct(layer_name)) %>%
  arrange(desc(layers_in_modules)) %>%
  filter(layers_in_modules==7) %>%
  group_by(id,module,layer_name) %>%
  mutate(nodes_in_layers=n_distinct(node_name)) %>%
  arrange(desc(nodes_in_layers)) %>%
  filter(module == 1) %>%
  transform(module = as.integer(module)) %>%
  tibble()

# combine observed and shuffled
observed_biggest_module$group <- 'obs'
shuffled_biggest_module$group <- 'shuff'
biggest_module_combined <- rbind(observed_biggest_module,shuffled_biggest_module)

obs <- biggest_module_combined %>% 
  group_by(group,id,layer_name) %>% 
  summarise(n=n_distinct(node_name)) %>% 
  filter(group=='obs')

png(filename = 'output/figures/observed_shuffled_biggest_module_curveball.png', width = 1600, height = 900, res = 300)
biggest_module_combined %>% 
  group_by(group,id,layer_name) %>% 
  summarise(n=n_distinct(node_name)) %>% 
  filter(group=='shuff') %>% 
  ggplot(aes(n))+geom_histogram()+
  facet_wrap(~layer_name)+ 
  theme_bw()+
  geom_vline(data=obs, aes(xintercept = n), color= 'red')+
  theme(axis.text = element_text(size=8 ,color='black'), axis.title = element_text(size=18 ,color='black'),title = element_text(size=16 ,color='black'))+
  labs(x = 'number of nodes',y='count', title = 'Observed vs Shuffled (large module)')
dev.off()

# histogram 
png(filename = 'output/figures/observed_shuffled_biggest_module.png', width = 1600, height = 900, res = 300)
biggest_module_combined %>%
  group_by(layer_name,group,id) %>%
  summarise(number_of_nodes=n_distinct(node_name)) %>%
  ggplot(aes(number_of_nodes, fill=group)) +
  geom_histogram(binwidth = 10, position = 'dodge') +
  theme_bw() +
  facet_wrap(~layer_name)+
  labs(x = 'number of nodes',y='count', title = 'Observed vs Shuffled')
dev.off()

# for each farm 
png(filename = 'output/figures/observed_shuffled_biggest_module_UK2.png', width = 1600, height = 1500, res = 300)
shuffled_biggest_module %>%
  filter(layer_name=="UK1") %>%
  group_by(group,id) %>%
  summarise(number_of_nodes=n_distinct(node_name)) %>%
  ggplot(aes(number_of_nodes)) +
  geom_histogram(color='white',binwidth = 14, position = 'dodge') +
  theme_bw() +
  geom_vline(aes(xintercept=482), color= 'red') +
  labs(x = 'Number of nodes',y='Count', title = 'Observed vs Shuffled-farm UK2')+
  theme(axis.text = element_text(size=18 ,color='black'), axis.title = element_text(size=18 ,color='black'),title = element_text(size=18 ,color='black'))
dev.off()

# intercept
biggest_module_combined %>%
  group_by(layer_name,group,id) %>%
  summarise(number_of_nodes=n_distinct(node_name)) %>%
  filter(layer_name=='UK1')

# combine all modules
farm_modules_observed <- farm_modules_pos %>%
  mutate(id='000') %>%
  select(c(id,module, node_name, layer_name=short_name)) %>%
  transform(id = as.character(id)) %>%
  tibble()

farm_modules_shuffled_within <- farm_modules_shuffled %>%
  select(-c(node_id,layer_id)) %>%
  transform(module = as.integer(module)) %>%
  tibble()

# combine observed and shuffled
farm_modules_observed$group <- 'obs'
farm_modules_shuffled_within$group <- 'shuff within'
farm_modules_combined <- rbind(farm_modules_observed,
                               farm_modules_shuffled_within %>% select(id, module, node_name, layer_name, group))
  
# comparison between observed and shuffled within only- percentage in each farm
farm_modules_filtered <- farm_modules_combined %>%
  group_by(group,id,node_name) %>%
  mutate(number_of_farms=n_distinct(layer_name)) %>%
  filter(number_of_farms==7)

biggest_module_filtered <- farm_modules_filtered %>%
  filter(module==1)

farm_modules_filtered %<>%
  group_by(group,id,layer_name) %>%
  summarise(number_of_nodes=n_distinct(node_name))

biggest_module_filtered %<>%
  group_by(group,id,layer_name) %>%
  summarise(number_of_nodes_big_module=n_distinct(node_name))

ASVs_big_module_farm <- left_join(farm_modules_filtered,biggest_module_filtered) 

ASVs_big_module_farm %<>%
  mutate(percent = 100 * number_of_nodes_big_module/number_of_nodes ) %>%
  select(-c(number_of_nodes,number_of_nodes_big_module)) %>%
  group_by(layer_name,group) %>%
  summarise(percent_mean=mean(percent))

png(filename = 'output/figures/nodes_occur_in_all_farms_curveball.png', width = 1600, height = 900, res = 300)
ASVs_big_module_farm %>%
  ggplot(aes(x=layer_name, y=percent_mean, fill=group)) +
  geom_bar(stat="identity", position="dodge")+
  theme_bw()+
  labs(x = 'Layer',y='percent of nodes in the largest module', title = 'Percent of nodes in the largest module')+
  theme(axis.text = element_text(size=10 ,color='black'), axis.title = element_text(size=12 ,color='black'), title = element_text(size=10 ,color='black'), ) 
dev.off()

# percent of nodes in the big module
farm_modules_pos %>%
  mutate(number_of_nodes=n_distinct(node_id)) %>%
  filter(module==1) %>%
  mutate(number_of_nodes_big_m=n_distinct(node_id)) %>%
  mutate(percent = 100 * number_of_nodes_big_m/number_of_nodes) %>%
  distinct(percent)

# filter by rank-----
# Distribution of module sizes within layers
mod_layer_size_100 <- 
  farm_modules_shuffled %>%
  filter(id=='100') %>%
  group_by(layer_name,module) %>% 
  summarise(n=n_distinct(node_id)) %>%
  arrange(desc(n))
# A more vivid figure:
mod_layer_size_100$rank <- 1:nrow(mod_layer_size_100)
threshold_rank <- 10

pdf('local_output/figures/module_size_rank#2.pdf',10,6)
mod_layer_size_100 %>% 
  ggplot(aes(rank,n))+
  geom_point(color='#336BFF', size=3)+geom_line(color='#336BFF', size=2)+
  geom_vline(xintercept = threshold_rank, linetype='dashed')+
  labs(x='Module-layer size rank', y='Module-layer size')+
  theme_bw()+
  theme(panel.grid=element_blank(),
        axis.text = element_text(size=22, color='black'),
        axis.title = element_text(size=22, color='black'),
        legend.position = c(0.9,0.9))
dev.off()

# Put the module rank in the module data table
farm_modules_shuffled_100 <- farm_modules_shuffled %>%
  filter(id=='100')
farm_modules_shuffled_100 %<>% left_join(mod_layer_size_100) 

# Modules in farms shuff_100
png(filename = 'local_output/figures/modules_in_layers_all_shuff_100_curveball.png', width = 1300, height = 900, res = 300)
farm_modules_shuffled_100 %>%
  group_by(layer_name) %>%
  mutate(nodes_in_layers=n_distinct(node_id)) %>%
  group_by(layer_name,module) %>%
  mutate(nodes_in_modules=n_distinct(node_id)) %>%
  mutate(nodes_percent=nodes_in_modules/nodes_in_layers) %>%
  group_by(layer_name, module,rank,nodes_percent) %>%
  summarise(nodes=n_distinct(node_id)) %>%
  ggplot(aes(x = module, y = layer_name, fill=nodes_percent))+
  scale_x_continuous(breaks = seq(1,16,2))+
  geom_tile(color='white')+
  scale_fill_viridis_c()+
  labs(x='Module ID', y='')+ ggtitle('shuffled Modules')+
  theme_bw()+
  #labs(tag = "b") +
  theme(panel.grid.minor = element_blank(),
        axis.text = element_text(size=10, color='black'),
        axis.title = element_text(size=10, color='black'), title = element_text(size = 10, color = 'black'),
        legend.title = element_text(size=7, color='black'), legend.text = element_text(size=7, color='black'))
dev.off()

# Filtered by rank size
pdf('output/figures/modules_in_layers_highly_ranked_100.pdf',10,8)
farm_modules_shuffled_100 %>%
  group_by(layer_name, module,rank) %>%
  summarise(nodes=n_distinct(node_id)) %>%
  filter(rank<=threshold_rank) %>%
  ggplot(aes(x = module, y = layer_name, fill=nodes))+
  geom_tile(color='white')+
  scale_fill_viridis_c()+
  labs(x='Module ID', y='Farm', title = 'Modules in layers highly ranked (shuffled 100)')+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        axis.text = element_text(size=20, color='black'),
        axis.title = element_text(size=20, color='black'), title = element_text(size = 20, color = 'black'))
dev.off()

# filter the smallest modules----
# for the observed network
farm_modules_pos_03 <- read_csv('local_output/farm_modules_pos_30_U.csv')
farm_modules_pos_03 %>%
  group_by(module) %>%
  summarise(n=n_distinct(node_id)) %>%
  filter(n>1) 

# for the shuffled networks
farm_modules_shuff_filtered <- farm_modules_shuffled %>%
  group_by(id,module) %>% 
  summarise(n=n_distinct(node_id)) %>%
  filter(n>1) 

png(filename = 'local_output/figures/observed_shuffled_hist_filter.png', width = 1600, height = 900, res = 300)
farm_modules_shuff_filtered %>%
  group_by(id) %>%
  summarise(module_sum=n_distinct(module)) %>%
  ggplot(aes(module_sum)) +
  geom_histogram() +
  theme_bw() +
  geom_vline(aes(xintercept=8), color= 'red') + 
  geom_text(aes(x=7.6,label='Observed', y=450, angle=90)) +
  theme(axis.text.x = element_text(size=10, color='black'),
        axis.text.y = element_text(size = 10, color='black'))+
  labs(y='count', title = 'Observed vs Shuffled (filter)')
dev.off()

# filter by rank-----
# Distribution of module sizes within layers
mod_layer_size_100 <- 
  farm_modules_shuffled %>%
  filter(id=='001') %>%
  group_by(layer_name,module) %>% 
  summarise(n=n_distinct(node_id)) %>%
  arrange(desc(n))
# A more vivid figure:
mod_layer_size_100$rank <- 1:nrow(mod_layer_size_100)
threshold_rank <- 10

pdf('local_output/figures/module_size_rank.pdf',10,6)
mod_layer_size_100 %>%
  ggplot(aes(rank,n))+
  geom_point(color='#336BFF', size=3)+geom_line(color='#336BFF', size=2)+
  geom_vline(xintercept = threshold_rank, linetype='dashed')+
  labs(x='Module-layer size rank', y='Module-layer size')+
  theme_bw()+
  theme(panel.grid=element_blank(),
        axis.text = element_text(size=22, color='black'),
        axis.title = element_text(size=22, color='black'),
        legend.position = c(0.9,0.9))
dev.off()

# Put the module rank in the module data table
farm_modules_shuffled_100 <- farm_modules_shuffled %>%
  filter(id=='001')
farm_modules_shuffled_100 %<>% left_join(mod_layer_size_100) 

# Modules in farms shuff_100
pdf('local_output/figures/modules_in_layers_shuffled_within.pdf',10,8)
farm_modules_shuffled_100 %>%
  group_by(layer_name, module,rank) %>%
  summarise(nodes=n_distinct(node_id)) %>%
  ggplot(aes(x = module, y = layer_name, fill=nodes))+
  scale_x_continuous(breaks = seq(1, number_of_modules$module_sum[100], 2), limits = c(0,number_of_modules$module_sum[001]))+
  geom_tile(color='white')+
  scale_fill_viridis_c()+
  labs(x='Module ID', y='Farm')+ ggtitle('Modules Shuffled within layers')+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        axis.text = element_text(size=20, color='black'),
        axis.title = element_text(size=20, color='black'), title = element_text(size = 20, color = 'black'))
dev.off()

# Filtered by rank size
pdf('output/figures/modules_in_layers_highly_ranked_100.pdf',10,8)
farm_modules_shuffled_100 %>%
  group_by(layer_name, module,rank) %>%
  summarise(nodes=n_distinct(node_id)) %>%
  filter(rank<=threshold_rank) %>%
  ggplot(aes(x = module, y = layer_name, fill=nodes))+
  #scale_x_continuous(breaks = seq(1, m$m, 2), limits = c(0,m$m))+
  geom_tile(color='white')+
  scale_fill_viridis_c()+
  labs(x='Module ID', y='Farm', title = 'Modules in layers highly ranked (shuffled 100)')+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        axis.text = element_text(size=20, color='black'),
        axis.title = element_text(size=20, color='black'), title = element_text(size = 20, color = 'black'))
dev.off()


# pairwise similarity big module-----
# observed:
presence_absence_big_module <- observed_biggest_module %>%
  distinct(layer_name, node_name) %>%
  mutate(present=1) %>% 
  spread(node_name, present, fill = 0)  %>%
  column_to_rownames('layer_name')
presence_absence_big_module <- as.matrix(presence_absence_big_module)
dim(presence_absence_big_module)
beta_div_big_module <- 1-as.matrix(vegdist(presence_absence_big_module, "jaccard"))
# Heatmap 
beta_div_big_module_m <- melt(beta_div_big_module)
png(filename = 'output/figures/shared_nodes_observed_big_module.png', width = 1600, height = 900, res = 300)
ggplot(beta_div_big_module_m, aes(x = Var1, y = Var2, fill = value, label=round(value,2))) +
  geom_tile() + 
  geom_text()+
  theme_bw()+
  scale_fill_gradient(high = "blue", low = "light blue") +
  labs(x = '',y='', title = 'Shared nodes between farms (big module)')+
  theme(axis.text = element_text(size=14, color='black'), title = element_text(size = 14, color = 'black'))
dev.off()

beta_div_big_module_m %<>%
  tibble() %>%
  filter(value!=1) %>%
  distinct(value, .keep_all = TRUE) %>%
  mutate(id='000') %>%
  unite(farm_pairs,Var1:Var2, sep = '-')

# shuffled within farms:
shuffled_biggest_module_sumry <- shuffled_biggest_module %>%
  distinct(id,layer_name, node_name) %>%
  mutate(present=1) %>%
  arrange(id)

# first we create a tibble includes all the farm pairs
# demo for one id
one_id_mat <- shuffled_biggest_module_sumry %>%
filter(id=='001') %>%
  spread(node_name, present, fill = 0)  %>%
  select(-id) %>%
  column_to_rownames('layer_name')
one_id_mat <- as.matrix(one_id_mat)
beta_div_big_module_shuffled <- 1-as.matrix(vegdist(one_id_mat, "jaccard"))
beta_div_big_module_shuffled_m <- melt(beta_div_big_module_shuffled)
beta_div_big_module_shuffled_m %<>%
  tibble() %>%
  filter(value!=1) %>%
  distinct(value, .keep_all = TRUE) %>%
  mutate(id='001') %>%
  unite(farm_pairs,Var1:Var2, sep = '-')

# a loop for all ids
beta_div_big_module_shuffled_final <- NULL
for (i in unique(shuffled_biggest_module_sumry$id)) { 
  print(i)
  presence_absence_big_module_shuffled <- shuffled_biggest_module_sumry %>%
    filter(id==i) %>%
    spread(node_name, present, fill = 0)  %>%
    select(-id) %>%
    column_to_rownames('layer_name')
  presence_absence_big_module_shuffled <- as.matrix(presence_absence_big_module_shuffled)
  beta_div_big_module_shuffled <- 1-as.matrix(vegdist(presence_absence_big_module_shuffled, "jaccard"))
  beta_div_big_module_shuffled_m <- melt(beta_div_big_module_shuffled)
  beta_div_big_module_shuffled_m %<>%
    tibble() %>%
    filter(value!=1) %>%
    distinct(value, .keep_all = TRUE) %>%
    mutate(id=i) %>%
    unite(farm_pairs,Var1:Var2, sep = '-')
  beta_div_big_module_shuffled_final <- rbind(beta_div_big_module_shuffled_final,beta_div_big_module_shuffled_m)
}

beta_div_big_module_m$group <- 'obs'
beta_div_big_module_shuffled_final$group <- 'shuff'
beta_div_big_module_combined <- rbind(beta_div_big_module_m,beta_div_big_module_shuffled_final)

beta_div_big_module_combined %>%
  ggplot(aes(value, fill=group)) +
  geom_histogram(binwidth = 2, position = 'dodge') +
  theme_bw() +
  facet_wrap(~farm_pairs)+
  labs(x = 'nodes similarity',y='count', title = 'Observed vs Shuffled nodes similarity')

# for each pair of farms
beta_div_big_module_shuffled_final %>%
  filter(farm_pairs=='SE1-FI1') %>%
  ggplot(aes(value)) +
  geom_histogram(color='white',binwidth = 0.008, position = 'dodge') +
  theme_bw() +
  geom_vline(aes(xintercept=0.0326), color= 'red')+
  labs(x = 'nodes similarity',y='count', title = 'Observed vs Shuffled nodes similarity SE1-FI1')


# analysis with shuffled between farms-----
# combine for biggest module with the between shuffled
between_shuffled_biggest_module <- farm_modules_shuffall %>%
  select(-c(node_id,layer_id, flow)) %>%
  group_by(id,module) %>%
  mutate(layers_in_modules=n_distinct(layer_name)) %>%
  arrange(desc(layers_in_modules)) %>%
  filter(layers_in_modules==7) %>%
  mutate(nodes_in_layers=n_distinct(node_name)) %>%
  arrange(desc(nodes_in_layers)) %>%
  filter(module == 1) %>%
  transform(module = as.integer(module)) %>%
  tibble()

observed_biggest_module$group <- 'obs'
shuffled_biggest_module$group <- 'shuff within'
between_shuffled_biggest_module$group <- 'shuff between'
biggest_module_combined_all <- rbind(observed_biggest_module,shuffled_biggest_module,between_shuffled_biggest_module)

# how many ASVs in each layer
biggest_module_combined_all %<>%
  group_by(layer_name,group) %>%
  summarise(number_of_ASVs_biggest_module=n_distinct(node_name))

# combine all
farm_modules_observed <- farm_modules_pos %>%
  mutate(id='000') %>%
  select(c(id,module, node_name, layer_name=short_name)) %>%
  transform(id = as.character(id)) %>%
  tibble()

farm_modules_shuffled_within <- farm_modules_shuffled %>%
  select(-c(node_id, layer_id, flow)) %>%
  transform(module = as.integer(module)) %>%
  tibble()

farm_modules_shuffled_between <- farm_modules_shuffall %>%
  select(-c(node_id, layer_id, flow)) %>%
  transform(module = as.integer(module)) %>%
  tibble()

# combine observed and shuffled
farm_modules_observed$group <- 'obs'
farm_modules_shuffled_within$group <- 'shuff within'
farm_modules_shuffled_between$group <- 'shuff between'
farm_modules_combined <- rbind(farm_modules_observed,farm_modules_shuffled_within,farm_modules_shuffled_between)

# how many ASVs in each layer
farm_modules_combined_ASVs <- farm_modules_combined %>%
  group_by(layer_name,group) %>%
  summarise(number_of_ASVs=n_distinct(node_name))  

# join farm_modules_combined and biggest_module_combined_all
farm_modules_comb <- left_join(biggest_module_combined_all,farm_modules_combined_ASVs)
farm_modules_comb %<>%
  mutate(relative_ASVs=number_of_ASVs_biggest_module/number_of_ASVs)

farm_modules_comb %>%
  ggplot(aes(x=layer_name,y=relative_ASVs, fill=group)) +
  geom_bar(stat="identity", position="dodge")+
  theme_bw()+
  theme(axis.text = element_text(size=10, color='black')) 

# how many ASVs that occurr in all 7 farms, are included in the big module 
# comparison between observed, shuffled within and shuffled between
farm_modules_combined_filtered <- farm_modules_combined %>%
  group_by(group,id,node_name) %>%
  mutate(number_of_farms=n_distinct(layer_name)) %>%
  filter(number_of_farms==7)

biggest_module_combined_filterd <- farm_modules_combined_filtered %>%
  filter(module==1)

farm_modules_combined_filtered %<>%
  group_by(group,id) %>%
  summarise(number_of_nodes=n_distinct(node_name))

biggest_module_combined_filterd %<>%
  group_by(group,id) %>%
  summarise(number_of_nodes_big_module=n_distinct(node_name))

ASVs_big_module <- left_join(farm_modules_combined_filtered,biggest_module_combined_filterd) 

ASVs_big_module %<>% melt(id.vars= c("group", "id"))   # converting wide data to long

png(filename = 'output/figures/ASVs_occur_in_all_farms.png', width = 1600, height = 900, res = 300)
ASVs_big_module %>%
  ggplot(aes(x=group, y=value, fill=variable)) +
  geom_bar(stat="identity", position="dodge")+
  theme_bw() +
  labs(x = '',y='number of nodes', title = 'ASVs that occur in all 7 farms, largest module vs all modules')+
  theme(axis.text = element_text(size=9 ,color='black'), axis.title = element_text(size=12 ,color='black'), title = element_text(size=10 ,color='black'), )
dev.off()