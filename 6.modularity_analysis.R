#------ modularity_analysis.r -------------------------------------


# Includes ------------
library(tidyverse)
library(magrittr)
library(igraph)
library(cowplot)#?
library(infomapecology)
check_infomap()
library(reshape2)
library(extRC)
library(GUniFrac)
library(data.table)
source('functions.R')

# Read data ---------------------------------------------------------------
farm_multilayer_pos_30 <-read_csv('local_output/farm_multilayer_pos_30.csv')
all_nodes <- sort(unique(c(farm_multilayer_pos_30$from, farm_multilayer_pos_30$to)))
all_nodes <- tibble(node_id=1:length(all_nodes), node_name=all_nodes)
layers <- tibble(layer_id=1:7, layer_name=c('NUDC', 'Park', 'Bian', 'Fran','Gand','Mink','Raab'),
                 short_name=c('UK1', 'UK2', 'IT1', 'IT2', 'IT3', 'FI1', 'SE1'))


# Build intralayer edges --------------------------------------------------
intra <- 
  farm_multilayer_pos_30 %>% 
  select(layer_from=level_name, node_from=from, layer_to=level_name, node_to=to, weight) %>% 
  left_join(layers, by = c('layer_from' = 'short_name')) %>% 
  left_join(layers, by = c('layer_to' = 'short_name')) %>% 
  left_join(all_nodes, by = c('node_from' = 'node_name')) %>% 
  left_join(all_nodes, by = c('node_to' = 'node_name')) %>% 
  select(layer_from=layer_id.x, node_from=node_id.x, layer_to=layer_id.y, node_to=node_id.y, weight)

# Observed modularity ------------------------------------------------------


## Interlayer links with Jaccard--------------------------------------------


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

# a for loop that calculates all the interlayer edges based on jaccard index
inter_PF_J <- NULL
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
  inter_PF_J <- rbind(inter_PF_J,inter_fid)
}

# Connect the intra and interlayer edges and change names to IDs
inter <- 
  inter_PF_J %>% 
  ungroup() %>% 
  left_join(layers, by = c('layer_from' = 'short_name')) %>% 
  left_join(layers, by = c('layer_to' = 'short_name')) %>% 
  left_join(all_nodes, by = c('ASV_ID' = 'node_name')) %>% 
  select(layer_from=layer_id.x, node_from=node_id, layer_to=layer_id.y, node_to=node_id, weight)

multilayer_jaccard <- rbind(intra %>% mutate(type='intra'),
                            inter %>% mutate(type='inter'))
table(multilayer_jaccard$type)
write_csv(multilayer_jaccard,'local_output/multilayer_jaccard.csv')

# Edge weight distributions
pdf('local_output/figures/edge_weights_jaccard_inter.pdf',8,8)
ggplot(multilayer_jaccard, aes(weight, fill=type))+
  geom_density(alpha=0.5)+
  labs(x='Edge weight', y='Density', title='Edge weight distributions')+
  scale_fill_manual(values = c('blue','orange'))+
  theme_bw()+
  theme(panel.grid=element_blank(),
        axis.text = element_text(size=22, color='black'),
        axis.title = element_text(size=22, color='black'),
        legend.position = c(0.9,0.9))
dev.off()


## Interlayer links with UniFrac -------------------------------------------
farm_multilayer_pos_final <- 
  bind_rows(farm_multilayer_pos_30, 
            farm_multilayer_pos_30 %>%
              relocate(to, from) %>%
              rename(to=from, from=to))
farm_multilayer_pos_final %<>%
  group_by(from) %>%
  mutate(num_farms_from=n_distinct(level_name)) %>%
  filter(num_farms_from>=2)
phylo_tree <- readRDS("local_output/fitted_asvs_phylo_tree.rds")

# a for loop that calculates all the interlayer edges based on unifrac
inter_PF_U <- NULL
for (i in unique(farm_multilayer_pos_final$from)) {
  print(i)
  ASV_net <- farm_multilayer_pos_final %>%
    filter(from==i)
  tree <- phylo_tree$tree
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
write_csv(multilayer_unif,'local_output/multilayer_unif.csv')

multilayer_unif <- read_csv('local_output/multilayer_unif.csv')
 
# Edge weight distributions
pdf('local_output/figures/edge_weights_unifrac_inter.pdf',8,8)
ggplot(multilayer_unif, aes(weight, fill=type))+
  geom_density(alpha=0.5)+
  labs(x='Edge weight', y='Density', tag = "(C)")+
  scale_fill_manual(values = c('blue','orange'))+
  theme_bw()+
  theme(panel.grid=element_blank(),
        axis.text = element_text(size=10, color='black'),
        axis.title = element_text(size=10, color='black'),
        title = element_text(size=10, color='black'),
        plot.tag = element_text(face = "bold"),
        legend.position = c(0.9,0.7))+
  paper_figs_theme
dev.off()

# Plot link distributions together
bind_rows(
  multilayer_jaccard %>% mutate(index='J'),
  multilayer_unif %>% mutate(index='U')
) %>% 
  mutate(grp=case_when(type=='intra' ~ 'intra',
                       type=='inter' & index=='J' ~ 'inter J',
                       type=='inter' & index=='U' ~ 'inter U')
  ) %>% 
  ggplot(aes(weight, fill=grp))+
  geom_density(alpha=0.5)+
  labs(x='Edge weight', y='Density', title='Edge weight distributions')+
  scale_fill_manual(values = c('red','blue','orange'))+
  theme_bw()+
  theme(panel.grid=element_blank(),
        axis.text = element_text(size=22, color='black'),
        axis.title = element_text(size=22, color='black'),
        legend.position = c(0.9,0.9))

## Run Infomap ------------------------------------------------------
multilayer_unif <- read_csv('local_output/multilayer_unif.csv')
net <- multilayer_unif[,1:5]

# Run Infomap
multilayer_for_infomap <- create_multilayer_object(extended = net, nodes = all_nodes, layers = layers)
m <- infomapecology::run_infomap_multilayer(multilayer_for_infomap, silent = F,
                                            flow_model = 'undirected',  
                                            trials = 200, relax = F, seed=123)



# Write summary to later compare with the shuffled networks
tibble(net='multilayer_unif', call=m$call, L=m$L, top_modules=m$m,time_stamp=Sys.time()) %>% 
write_csv('local_output/farm_modules_pos_30_summary.csv', append=T)
modules_obs <- m$modules %>% left_join(layers)
write_csv(modules_obs, 'local_output/farm_modules_pos_30_U.csv')
  
# Read from files if already run
modules_obs <- read_csv('local_output/farm_modules_pos_30_U.csv')
mod_summary_obs <- read_csv('local_output/farm_modules_pos_30_summary.csv', 
                            col_names = c('net', 'call', 'L', 'top_modules', 'time_stamp'))
# get the latest run
mod_summary_obs <- mod_summary_obs[which.max(mod_summary_obs$time_stamp),] 
num_modules_obs <- mod_summary_obs$top_modules
L_obs <- mod_summary_obs$L


## Analyze observed modularity results -------------------------------
# Distribution of module sizes
a <- modules_obs %>%
group_by(module) %>% 
summarise(n=n_distinct(node_id)) %>%
arrange(desc(n))

plot(a)

# And within layers
tot <- modules_obs %>%
  group_by(short_name) %>% 
  summarise(n_farm=n_distinct(node_id)) %>%
  arrange(desc(n_farm))

# percentages of modules in each the farm
modules_obs %>%
group_by(short_name,module) %>% 
summarise(n=n_distinct(node_id)) %>%
left_join(tot, by="short_name") %>%
mutate(per_of_farm=100*n/n_farm) %>% arrange(desc(per_of_farm))


modules_obs %<>% rename(level1=module)

png(filename = 'local_output/figures/modules_unif_no_thresh.png', width = 1200, height = 900, res = 300)
a <- modules_obs %>%
  mutate(short_name=factor(short_name, levels = c("UK1","UK2","IT1","IT2","IT3","FI1",'SE1'))) %>%
  group_by(short_name) %>%
  mutate(nodes_in_layers=n_distinct(node_id)) %>%
  group_by(short_name,level1) %>%
  mutate(nodes_in_modules=n_distinct(node_id)) %>%
  mutate(nodes_percent=nodes_in_modules/nodes_in_layers) %>%
  distinct(short_name, level1, nodes_percent) %>% 
  arrange(level1, short_name) %>%
  # Plot
  ggplot(aes(x = level1, y = short_name, fill=nodes_percent))+
  geom_tile(color='white')+
  scale_x_continuous(breaks = seq(1, max(modules_obs$level1), 1))+
  scale_fill_viridis_c(limits = c(0, 1))+
  theme_bw()+
  labs(x='Module ID', y='', title = "No threshold")+
  theme(panel.grid=element_blank(),
        axis.text = element_text(size=10, color='black'),
        axis.title = element_text(size=10, color='black'),
        title = element_text(size=10, color='black'),
        plot.tag = element_text(face = "bold")) +
  paper_figs_theme_no_legend
a
dev.off()

png(filename = 'local_output/figures/modules_unif_with_thresh.png', width = 1200, height = 900, res = 300)
b <- modules_obs %>%
  mutate(short_name=factor(short_name, levels = c("UK1","UK2","IT1","IT2","IT3","FI1",'SE1'))) %>%
  group_by(short_name) %>%
  mutate(nodes_in_layers=n_distinct(node_id)) %>%
  group_by(short_name,level1) %>%
  mutate(nodes_in_modules=n_distinct(node_id)) %>%
  mutate(nodes_percent=nodes_in_modules/nodes_in_layers) %>%
  distinct(short_name, level1, nodes_percent) %>% 
  arrange(level1, short_name) %>%
  # Only include modules that contain at least 3% of the ASVs in the layer
  filter(nodes_percent>=0.03) %>%
  # Plot
  ggplot(aes(x = level1, y = short_name, fill=nodes_percent))+
  geom_tile(color='white')+
  scale_x_continuous(breaks = seq(1, max(modules_obs$level1), 1))+
  scale_fill_viridis_c(limits = c(0, 1))+
  theme_bw()+
  labs(x='Module ID', y='', title = "With threshold (0.03)")+
  theme(panel.grid=element_blank(),
        axis.text = element_text(size=10, color='black'),
        axis.title = element_text(size=10, color='black'),
        title = element_text(size=10, color='black'),
        plot.tag = element_text(face = "bold"), 
        legend.title=element_blank())
b
dev.off()

pdf('local_output/figures/SI_modules_heir.pdf',10,6)
plot_grid(a, b, nrow = 1, ncol = 2, labels = c('(A)','(B)'),vjust = 1.1)
dev.off()


# Compared to shuffled networks -------------------------------------------
# read observed for comparing:
modules_obs <- read_csv('local_output/farm_modules_pos_30_U.csv')
mod_summary_obs <- read_csv('local_output/farm_modules_pos_30_summary.csv', 
                            col_names = c('net', 'call', 'L', 'top_modules', 'time_stamp'))
# get the latest run
mod_summary_obs <- mod_summary_obs[which.max(mod_summary_obs$time_stamp),] 
num_modules_obs <- mod_summary_obs$top_modules
L_obs <- mod_summary_obs$L

parent.folder <- "HPC/shuffled/shuffle_farm_r0_30_500_jac_intra"

# General stats
summ <- read_csv(paste(parent.folder,'/farm_modulation_summary_pf_unif.csv',sep=''), 
                 col_names = c('e_id','JOB','call','L','num_modules'))
summ %<>% slice(tail(row_number(), 1000)) %>% filter(!str_detect(call, '-2'))
ggplot(summ, aes(L))+
  geom_histogram()+
  geom_vline(xintercept = L_obs)
pvalue <- sum(summ$L<L_obs)/500

ggplot(summ, aes(num_modules))+
  geom_histogram()+
  geom_vline(xintercept = num_modules_obs)

# percent of shuffled networks with 1 big module:
100*sum(summ$num_modules == 1)/nrow(summ)

## Read shuffled files -------------------------------

# Folder containing sub-folders
# Sub-folders
sub.folders <- list.dirs(parent.folder, recursive=TRUE)[-1]

# get model nums (because there was NA in the summ for some reason)
module_nums <- NULL

# read all 
modules_shuffled <- NULL
for (s in sub.folders) {
  print(s)
  modules_run <- list.files(path = s , pattern = '_farm_modules_pf_unif.csv', recursive = T,full.names = T)
  modules_run <- fread(modules_run)
  modules_run$id <- str_split_fixed(s, pattern = '/', n = 4)[4]
  modules_shuffled <- rbind(modules_shuffled,modules_run)
  module_nums <- rbind(module_nums, data.table(origin=s, modules=max(modules_run$module)))
}
modules_shuffled <- as_tibble(modules_shuffled)

# percent of shuffled networks with 1 big module as read from the modules files:
100* (module_nums$modules == 1)/nrow(module_nums)


a2 <- ggplot(summ, aes(L))+
  geom_histogram()+
  geom_vline(xintercept = L_obs, linetype = "dashed") + paper_figs_theme

b2 <- ggplot(module_nums, aes(modules))+
  geom_histogram()+ xlab("Number of modules") +
  geom_vline(xintercept = num_modules_obs, linetype = "dashed")+ 
  paper_figs_theme

pdf('local_output/figures/SI_modularity_obs_shuff.pdf',10,6)
plot_grid(a2, b2, nrow = 1, ncol = 2, labels = c('(A)','(B)'),vjust = 1.1, hjust = 0)
dev.off()


# Flow dynamics -----------------------------------------------------------
modules_obs <- read_csv('local_output/farm_modules_pos_30_U.csv')

ASV_taxa <- read_csv('local_output/ASV_full_taxa.csv') %>% 
  select(node_name=ASV_ID, everything(), -seq16S)
modules_obs %<>% rename(Farm=short_name)

ggplot(modules_obs, aes(flow, fill=Farm))+geom_histogram()+
  facet_wrap(~Farm,scales = 'free_y')

sum(modules_obs$flow)

# This correlates with the number of ASVs in the farm, but not completely
farm_flow <- modules_obs %>% 
  group_by(Farm) %>% 
  summarise(farm_flow=sum(flow), ASV_num=n_distinct(node_name)) %>% arrange(desc(farm_flow))
ggplot(farm_flow, aes(ASV_num, farm_flow))+
  geom_point(size=3, color='orange')+
  paper_figs_theme_no_legend

# Flow reflects the taxonomy
modules_obs %>% 
  left_join(ASV_taxa) %>% 
  group_by(Farm, Phylum) %>% 
  summarise(taxa_flow=sum(flow), ASV_num=n_distinct(node_name)) %>% 
  mutate(relative_flow=taxa_flow/sum(taxa_flow)) %>% 
  mutate(relative_richness=ASV_num/sum(ASV_num)) %>% 
  mutate(ypos = cumsum(relative_flow)- 0.5*relative_flow) %>% 
  ggplot(aes(x="", y=relative_flow, fill=Phylum))+
  facet_wrap(~Farm)+
  geom_bar(stat="identity", width=1) +
  # scale_fill_manual(values = signif_colors)+
  # geom_text(aes(y = ypos, label = round(prop,2)), color = "white", size=3) +
  coord_polar("y", start=0)+
  paper_figs_theme+
  theme(axis.text = element_blank(),
        axis.title = element_blank())

# Flow distribution - node level
modules_obs %>% 
  group_by(node_name) %>% 
  summarise(node_flow=sum(flow)) %>% 
  mutate(node_flow_norm=node_flow/max(node_flow)) %>% 
  arrange(desc(node_flow)) %>% 
  mutate(flow_cum=cumsum(node_flow)) %>% 
  left_join(ASV_taxa) %>% 
  ggplot(aes(node_flow))+
  geom_histogram()+
  paper_figs_theme

top_flow_nodes <- 
modules_obs %>% 
  group_by(node_name) %>% 
  summarise(node_flow=sum(flow)) %>% 
  mutate(node_flow_norm=node_flow/max(node_flow)) %>% 
  arrange(desc(node_flow)) %>% 
  mutate(flow_cum=cumsum(node_flow)) %>% 
  left_join(ASV_taxa) %>% 
  first(93) # 10% most influential nodes
  # filter(flow_cum<=0.8)

sum(top_flow_nodes$node_flow) # the sum of the flow of the top flowing 10% ASVs

top_flow_nodes %>% 
  group_by(Phylum) %>% 
  summarise(R=n_distinct(node_name)) %>% 
  arrange(desc(R))

modules_obs %>% 
  group_by(node_name) %>% 
  summarise(node_flow=sum(flow)) %>% 
  mutate(node_flow_norm=node_flow/max(node_flow)) %>% 
  arrange(desc(node_flow)) %>% 
  mutate(flow_cum=cumsum(node_flow)) %>% 
  left_join(ASV_taxa) %>% 
  ggplot(aes(x=(1:946)/946, y=flow_cum))+geom_point()+
  geom_hline(yintercept = 0.8, color='red')+
  # geom_segment(aes(x = 0, y = 0.8, xend = nrow(top_flow_nodes), yend = 0.8 , colour = "red"))+
  # geom_segment(aes(x = nrow(top80), y = 0.8, xend =nrow(top_flow_nodes), yend = 0, colour = "red"))+
  paper_figs_theme

# Compare shuffled to observed flow distributions
flow_obs <- 
modules_obs %>% 
  group_by(node_name) %>% 
  summarise(node_flow=sum(flow)) %>% 
  arrange(desc(node_flow)) %>% 
  mutate(grp='obs')
flow_shuff <- 
  modules_shuffled %>% 
  group_by(id,node_name) %>% 
  summarise(node_flow=sum(flow)) %>% 
  # group_by(node_name) %>% 
  # summarise(node_flow=mean(node_flow)) %>% 
  mutate(grp='shuff') %>% 
  ungroup() %>%
  select(-id)

bind_rows(flow_obs, flow_shuff) %>% 
  ggplot(aes(node_flow, fill=grp))+geom_density(alpha=0.5)+
  paper_figs_theme


# Cumulative curves
cum_flow_obs <- 
modules_obs %>% 
  group_by(node_name) %>% 
  summarise(node_flow=sum(flow)) %>% 
  arrange(desc(node_flow)) %>% 
  mutate(flow_cum=cumsum(node_flow))

cum_flow_shuff <- 
  modules_shuffled %>% 
  group_by(id,node_name) %>% 
  summarise(node_flow=sum(flow)) %>% 
  group_by(id) %>% 
  mutate(flow_cum=cumsum(node_flow),
         row_id=1:n())


bind_rows(
cum_flow_obs %>%  mutate(grp='obs'),
cum_flow_shuff %>%  
  group_by(row_id) %>%
  summarise(flow_cum_shuff_mean=mean(flow_cum),
            flow_cum_shuff_sd=mean(sd)) %>% mutate(grp='shuff')
) %>% 
  ggplot(aes(x=row_id, y=flow_cum, group=row_id, fill=grp))+
  geom_boxplot(outlier.colour = NA)+
  paper_figs_theme


# A test with a shuffled network
modules_shuffled %>% 
  group_by(id,node_name) %>% 
  summarise(node_flow=sum(flow)) %>% 
  mutate(flow_cum=cumsum(node_flow)) %>% 
  left_join(ASV_taxa) %>% 
  ggplot(aes(node_flow))+geom_histogram()+
  paper_figs_theme



z_score <- 
  modules_shuffled %>%
  group_by(node_name) %>% 
  summarise(node_flow_shuff_mean=mean(flow),
            node_flow_shuff_sd=sd(flow)) %>% 
  inner_join(
    modules_obs %>% 
      group_by(node_name) %>% 
      summarise(node_flow=sum(flow))
  ) %>%
  mutate(z=(node_flow-node_flow_shuff_mean)/node_flow_shuff_sd) %>% 
  mutate(signif=case_when(z>1.96 ~ 'above', # Obs is more than the shuffled
                          z< -1.96 ~ 'below', # Obs is lower than the shuffled
                          z<=1.96 | z>=-1.96 ~ 'not signif')) %>% 
  filter(node_name %in% top_flow_nodes$node_name)
# What proportion of ASVs have a statistical significant change in flow?
z_score %>% 
  group_by(signif) %>% 
  summarise(n=n(),prop=n/nrow(z_score))


# Compared to the regional network ----------------------------------------

# Co-occurrence network for region scale 30%----
region_monolayer <- list.files(path = "HPC/exp_20", pattern = paste('_Region_edge_list.csv',sep=""), full.names = T)
region_monolayer_pos <- sapply(region_monolayer, read_csv, simplify=FALSE) %>%
  bind_rows(.id = "id") %>% 
  select(-id) %>%
  filter (edge_type=='pos')

all_nodes <- sort(unique(c(region_monolayer_pos$from, region_monolayer_pos$to)))
all_nodes <- tibble(node_id=1:length(all_nodes), node_name=all_nodes)

monolayer <- 
  region_monolayer_pos %>% 
  select(from, to, weight) %>% 
  left_join(all_nodes, by = c('from' = 'node_name')) %>% 
  left_join(all_nodes, by = c('to' = 'node_name')) %>% 
  select(from=node_id.x, to=node_id.y, weight) %>%
  transform(from = as.character(from), to = as.character(to)) %>%
  tibble()

# Prepare network objects
network_object <- create_monolayer_object(x=monolayer[,1:3], 
                                          directed = F, bipartite = F, 
                                          node_metadata = all_nodes)

# This is needed because there was a bug in the function in the infomapecology package
run_infomap_monolayer <- function(x, infomap_executable='Infomap', flow_model=NULL, silent=F, trials=100, two_level=T, seed=NULL, signif=F, shuff_method=NULL, nsim=1000, verbose=T, remove_auxillary_files=T, ...){
  
  # Check stuff and prepare
  if(check_infomap(infomap_executable)==F){stop('Error in Infomap standalone file.')}
  if(class(x)!='monolayer'){stop('x must be of class monolayer')}
  if (verbose){ print('Creating a link list...') }
  obs <- create_infomap_linklist(x) # obs = Observed network
  nodes <- x$nodes
  
  # Run Infomap for observed

  # Infomap arguments
  arguments <- paste(' -N ',trials,' -f ',flow_model,sep='')
  arguments <- ifelse(!is.null(seed), paste(arguments, '--seed',seed), arguments)
  arguments <- ifelse(silent, paste(arguments, '--silent'), arguments)
  arguments <- ifelse(two_level, paste(arguments, '--two-level'), arguments)
  arguments <- paste(arguments,...)
  call <- paste('./',infomap_executable, ' infomap.txt . --tree ',arguments,sep='')
  
  # Write temporary file for Infomap
  write_delim(obs$edge_list_infomap, 'infomap.txt', delim = ' ', col_names = F)
  # Run Infomap
  if (verbose){ cat('running: ');cat(call);cat('\n') }
  system(call)
  tmp <- str_split(read_lines('infomap.tree')[5], 'levels')
  num_levels <- parse_number(tmp[[1]][1]) # Get number of levels
  m <- parse_number(tmp[[1]][2]) # Get number of modules
  L <- parse_number(read_lines('infomap.tree')[6]) # Get the map equation value, L
  # Read module results from Infomap
  modules <- suppressMessages(read_delim('infomap.tree', delim = ' ', skip = 8, col_names = c('path', 'flow', 'name', 'node_id')))
  # Get the modules
  suppressWarnings( # Need to suppress warnings because with hierarchical modules there are NA generated.
    modules %<>%
      select(node_id, path, flow) %>%
      mutate(levels=str_count(path, pattern  = ':') ) %>%
      select(node_id, levels, path, flow) %>%
      separate(path, into=paste('module_level',1:num_levels,sep=''), sep = ':') %>%
      mutate_all(as.double) %>%
      full_join(obs$nodes, 'node_id') %>%
      select(node_id, node_name, flow, levels, starts_with('module_level'), everything()) %>%
      # rename(leaf_id=paste('module_level',num_levels,sep='')) %>% 
      arrange(node_id)
  )
  # Prepare first output, before significance testing
  out <- list(call=call, L=L, m=m, modules=modules, edge_list=x$edge_list, L_sim=NULL, m_sim=NULL, pvalue=NULL)
  
  if (signif){
    if (class(shuff_method)=='list'){ # If list with shuflled networks is provided.
      shuffled_linklist <- shuff_method 
      nsim <- length(shuffled_linklist)
    } else { 
      shuffled_linklist <- shuffle_infomap(x, shuff_method=shuff_method, nsim=nsim) # Otherwise, shuffle according to arguments.
    }
    # Run Infomap for shuffled
    L_sim <- m_sim <- NULL
    for (i in 1:nsim){
      write_delim(shuffled_linklist[[i]], 'infomap.txt', delim = ' ', col_names = F) # Write temporary file for Infomap
      if (verbose){ print(paste('Running Infomap on shuffled network ',i,'/',nsim,sep=''))}
      system(call) # Run Infomap
      L_sim <- c(L_sim, parse_number(read_lines('infomap.tree')[6])) # Get the map equation value, L
      m_sim <- c(m_sim, parse_number(str_split(read_lines('infomap.tree')[5], 'levels')[[1]][2])) # Get number of modules
    }
    
    out$L_sim <- L_sim
    out$m_sim <- m_sim
    out$pvalue <- sum(L_sim<L)/nsim 
    if(out$pvalue==0){warning(paste('pvalue is not really 0, it is <',1/nsim,sep=''))}
  }
  
  if (remove_auxillary_files){
    if (verbose){ print('Removing auxilary files...') }
    file.remove('infomap.txt')
    file.remove('infomap.tree')
  }
  class(out) <- 'infomap_monolayer'
  return(out)
}



# Run infomap with hieararchy
infomap_object <- run_infomap_monolayer(network_object, 
                                        infomap_executable='Infomap',
                                        flow_model = 'undirected',
                                        silent=F,trials=200, two_level=F, 
                                        seed = NULL, remove_auxillary_files=F)
x=infomap_object$modules
unique(x$levels)
n_distinct(x$module_level2)
max(x$module_level2)

write_csv(x,'local_output/regional_multilevel_modules.csv')

# Run infomap without hieararchy
infomap_object <- run_infomap_monolayer(network_object, 
                                        infomap_executable='Infomap',
                                        flow_model = 'undirected',
                                        silent=F,trials=200, two_level=T, 
                                        seed = 123, remove_auxillary_files=F)
x=infomap_object$modules
unique(x$module_level1)
table(x$module_level1)

write_csv(x,'local_output/regional_monolevel_modules.csv')
x <- read_csv('local_output/regional_monolevel_modules.csv')


# Draw a multilayer network -----------------------------------------------
multilayer_unif <- read_csv('local_output/multilayer_unif.csv')

layer_network <- 
  multilayer_unif %>%
  filter(layer_from!=layer_to) %>% 
  group_by(layer_from,layer_to) %>% 
  summarise(weight=length(type)) %>% 
  left_join(layers, by=c('layer_from' = 'layer_id')) %>%
  left_join(layers, by=c('layer_to' = 'layer_id')) %>%
    select(layer_from=short_name.x, layer_to=short_name.y, weight) %>% 
  graph_from_data_frame(directed = F)
plot(layer_network)


E(layer_network)$weight
V(layer_network)$name
# V(layer_network)$color <- c('#7E1F30', '#B95466', 'red', 'red', 'yellow', 'gray','gray')
ASV_data <- read_csv('local_output/core_ASV_fixed_30.csv')
num_ASV <- ASV_data %>%
  group_by(Farm) %>%
  summarise(ASV_num=n_distinct(ASV_ID))
V(layer_network)$num_ASV <- num_ASV$ASV_num[match(V(layer_network)$name, num_ASV$Farm)]


pdf('fixed_data/local_output/figures/layer_network.pdf',6,6)
plot(layer_network, edge.width=E(layer_network)$weight/50, edge.color='black', layout=layout.circle,
     vertex.size=V(layer_network)$num_ASV/40+20, vertex.color=NA)
dev.off()

# Draw modules for each farm pie charts -----
modules_obs <- read_csv('local_output/farm_modules_pos_30_U.csv')
modules_obs %<>% select(node=node_id, module, farm=short_name)

pdf('local_output/figures/layer_modules_composition.pdf',6,6)
modules_obs %>% 
  group_by(farm,module) %>% 
  summarise(n=n()) %>% 
  mutate(prop = n / sum(n)) %>% 
  mutate(N = sum(n)) %>% 
  mutate(Module=factor(module, levels=1:max(module))) %>% 
  #mutate(ypos = cumsum(prop)- 0.5*prop ) %>% 
  ggplot(aes(x="", y=prop, fill=Module))+
  facet_wrap(~farm)+
  geom_bar(stat="identity", width=1) +
  #scale_fill_manual(values = c('blue','orange','#32a852'))+
  # geom_text(aes(y = ypos, label = round(prop,2)), color = "white", size=3) +
  coord_polar("y", start=0)+
  paper_figs_theme_no_legend+theme_void()
dev.off()


# compare a farm's density to its connectivity to other farms ------
multilayer_unif <- read_csv('local_output/multilayer_unif.csv')

# count the intralayer first 
intra_nets <- multilayer_unif %>% filter(type=="intra") %>% 
        group_by(layer_from) %>% 
        summarise(n_intra=n()) %>% rename(from=layer_from)

# make the edgelist double and reversed to count each interlayer twice
intr <- multilayer_unif %>% filter(type=="inter") %>% 
  select(from=layer_from, to=layer_to)
otherway <- intr %>% 
  relocate(to, from) %>%
  rename(to=from, from=to)
res <- bind_rows(intr, otherway) %>% 
  group_by(from) %>%
  summarise(n_inter=n()) %>% left_join(intra_nets, by="from") %>%
  # this calculates the denciety
  mutate(outwardness=n_inter/(n_inter+n_intra)) %>% 
  mutate(relative=n_inter/n_intra) %>% 
  left_join(layers, by= c('from' = 'layer_id')) %>%
  select(farm=short_name, outwardness, relative)

#number of modules in each layer:
a <- cbind(res, n_m=c(1, 3, 1, 1, 5, 3, 1))

# is there a correlations? - density vs n_modules
shapiro.test(a$outwardness) #X
shapiro.test(a$n_m) #X

cor(a$outwardness, a$n_m) # calculates correlation coefficient
cor.test(a$outwardness, a$n_m, method="pearson")  

# is there a correlations? - density vs n_modules
shapiro.test(a$relative) #X

cor(a$relative, a$n_m) # calculates correlation coefficient
cor.test(a$relative, a$n_m, method="pearson") 
    