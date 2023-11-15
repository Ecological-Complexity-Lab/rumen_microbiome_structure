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
multilayer_jaccard <- read_csv('local_output/multilayer_jaccard.csv')

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
summ %<>% slice(tail(row_number(), 1000)) %>% filter(str_detect(call, '-2')) # We want only runs with *no* multi-level analysis
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

# statistical tests:

# is there a correlations? - density vs n_modules
shapiro.test(a$outwardness) #X
shapiro.test(a$n_m) #X

cor(a$outwardness, a$n_m) # calculates correlation coefficient
cor.test(a$outwardness, a$n_m, method="pearson")  

# is there a correlations? - density vs n_modules
shapiro.test(a$relative) #X

cor(a$relative, a$n_m) # calculates correlation coefficient
cor.test(a$relative, a$n_m, method="pearson") 
    