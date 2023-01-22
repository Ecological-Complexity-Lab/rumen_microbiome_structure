#-----
# runs infomap on a multilayer network with abundance based interlayer edges.
# input argumemnst: experiment id number from the experiments.csv file
# input files: _edge_list.csv file with intralayer edges
#              experiments.csv file with the names of abundance data files by exp_id
#              <abundance_file> file with abundance data (name read from experiments.csv) - usually the core ASV file.
# other requirements: working infomap exe (make sure it is the ubuntu version, and has running permission)
#                     functions.R
# all need to be in the working directory.
#------


#! /gpfs0/shai/projects/R4/R-4.0.3/bin/Rscript
.libPaths("/gpfs0/shai/projects/R4/R-4.0.3/lib64/R/library")
print(.libPaths())
print(sessionInfo())

library(tidyverse)
library(magrittr)
library(reshape2) 
library(igraph)
library(infomapecology)

# parse arguments
JOB_ID <- Sys.getenv("JOB_ID")
if (length(commandArgs(trailingOnly=TRUE))==0) {
  stop('No arguments were found!')
} else {
  args <- commandArgs(trailingOnly=TRUE)
  net_id <- as.numeric(args[1])
}

#net_id=100

print(getwd())
source('functions.R')
check_infomap()

# building co-occurrance intralayer edges ----
one_shuff_farm_f <- list.files(pattern = paste("^",net_id,"_",sep=''), full.names = T)
print(one_shuff_farm_f)
one_shuff_farm_f <- str_subset(one_shuff_farm_f, "edge_list.csv$")
print(one_shuff_farm_f)
#bind all farms to a single multilayer
farm_multilayer <- sapply(one_shuff_farm_f, read_csv, simplify=FALSE) %>% 
  bind_rows(.id = "id") %>% select(-id)

# For the analysis separate the positive and negative networks
farm_multilayer_pos <- farm_multilayer %>% filter(edge_type=='pos')
farm_multilayer_neg <- farm_multilayer %>% filter(edge_type=='neg')

all_nodes <- sort(unique(c(farm_multilayer_pos$from, farm_multilayer_pos$to)))
all_nodes <- tibble(node_id=1:length(all_nodes), node_name=all_nodes)

layers <- tibble(layer_id=1:7, layer_name=unique(farm_multilayer_pos$level_name))
farm_multilayer_pos %<>% 
  left_join(all_nodes, by = c('from' = 'node_name')) %>% 
  left_join(all_nodes, by = c('to' = 'node_name')) %>% 
  left_join(layers, by = c('level_name' = 'layer_name')) %>% 
  select(layer=layer_id, node_from=node_id.x, node_to=node_id.y, weight)

# building relative abundance interlayer edges ----
# read abundance file using name from experiments.csv
experiment_info <- read_csv('experiments.csv') %>% filter(e_id==net_id)
abundance_file <- experiment_info$Abundance_file
abundance_data <- read_csv(abundance_file)

# Total abundance in each farm
abundances <-
  abundance_data %>%
  group_by(Farm, ASV_ID) %>%
  summarise(tot_abund=sum(Abundance)) # Total abundance of the ASV in the specific farm, across all cows
# Relative abundance in each farm
farm_abundances <-
  abundance_data %>% group_by(Farm) %>% summarise(A=sum(Abundance))
abundances %<>%
  left_join(farm_abundances) %>%
  mutate(rel_abund=tot_abund/A) # Calculate the relative abundance of each ASV in a farm

# Check relative abundances are calculated correctly
abundances %>% group_by(Farm) %>% summarise(p=sum(rel_abund)) 

# remove ASVs that appear in only 1 farm because they do not have interlayer edges
one_farm_asvs <- names(which(table(abundances$ASV_ID)==1))
abundances %<>% filter(!ASV_ID %in% one_farm_asvs)

# create interlayer edges
interlayer_edges <- 
  abundances %>%
  group_by(ASV_ID) %>% 
  group_modify(~calc_abund_inter_ASV(.x))

# Remove nodes from the interlayer edges that do not have intralayer edges in
# any layer and are therefore not in the all_nodes tibble
nodes_to_remove <- setdiff(unique(interlayer_edges$ASV_ID),unique(all_nodes$node_name))
interlayer_edges %<>% filter(!ASV_ID%in%nodes_to_remove)

# mold the df to the right format
inter <- 
  interlayer_edges %>% 
  ungroup() %>% 
  left_join(layers, by = c('layer_from' = 'layer_name')) %>% 
  left_join(layers, by = c('layer_to' = 'layer_name')) %>% 
  left_join(all_nodes, by = c('ASV_ID' = 'node_name')) %>% 
  select(layer_from=layer_id.x, node_from=node_id, layer_to=layer_id.y, node_to=node_id, weight)

nrow(farm_multilayer_pos)
nrow(inter)

# Run Infomap ----
x <- create_multilayer_object(intra = farm_multilayer_pos, inter = inter, 
                              nodes = all_nodes, layers = layers)

m <- run_infomap_multilayer(x, flow_model = 'undirected', silent = F, trials = 10)
farm_modules_pos <- m$modules %>% left_join(layers)

# In the file name include net_id, farm
write_csv(farm_modules_pos, paste(net_id,'_farm_modules.csv',sep=""))

# Write summary when done ------------------
m_file <- tibble(e_id=net_id, JOB_ID=JOB_ID, call=m$call, L=m$L, modules_num=m$m)

summary_file <- 'farm_modulation_summary.csv'
# if data is shuffled save summery one directory up
if (grepl("shuff", experiment_info$data_file, ignore.case = TRUE)) { 
  summary_file <- '../farm_modulation_summary.csv' 
} 
write_csv(m_file, summary_file, append = T)
