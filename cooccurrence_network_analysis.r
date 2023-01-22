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
e_id <- 12
run_summary <- read_csv('HPC/run_summary.csv', 
                        col_names = c('exp_id','level','level_name','JOB_ID','data_file','time_stamp')) %>%
               arrange(exp_id,level)
run_summary %>% filter(exp_id==e_id)
layers <- tibble(layer_id=1:7, layer_name=c('NUDC', 'Park', 'Bian', 'Fran','Gand','Mink','Raab'),
                 short_name=c('UK1', 'UK2', 'IT1', 'IT2', 'IT3', 'FI1', 'SE1'))

# Create a multilayer network for 7 farms with intralayer edges ----

farm_multilayer <- NULL
setwd(paste('HPC/exp_',e_id, sep=''))
lyrs_list <- layers$short_name
if (e_id>4) { # this is compatible with short farm name (the new ones)
  lyrs_list <- layers$short_name
} else{ # this is compatible with full farm name (the old ones)
  lyrs_list <- layers$layer_name
}

# build network layers using old cooccure values -------
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
write_csv(farm_multilayer_pos, 'local_output/farm_multilayer_pos_30.csv')
farm_multilayer_neg <- farm_multilayer %>% filter(edge_type=='neg')
write_csv(farm_multilayer_neg, 'local_output/farm_multilayer_neg_30.csv')


# Build network based on p value instead of the one used so far ------
e_id <- 12
layers <- tibble(layer_id=1:7, layer_name=c('NUDC', 'Park', 'Bian', 'Fran','Gand','Mink','Raab'),
                 short_name=c('UK1', 'UK2', 'IT1', 'IT2', 'IT3', 'FI1', 'SE1'))
farm_multilayer <- NULL
lyrs_list <- layers$short_name
run_summary <- read_csv('run_summary.csv', 
                        col_names = c('exp_id','level','level_name','JOB_ID','data_file','time_stamp')) %>%
                arrange(exp_id,level)

for (l in lyrs_list){ 
  print('------------------')
  print(l)
  print('------------------')
  x <- parse_networks_from_cooc(e_id = e_id, Level = 'Farm', Level_name = l)
  farm_multilayer <- rbind(farm_multilayer,x$edge_list)
}

# add jaccard values, and using them as weight:
farm_multilayer %<>% rename(size_effect=weight) %>% 
  mutate(weight=obs_cooccur/(sp1_inc+sp2_inc-obs_cooccur))

# save the positive edges
write_csv(farm_multilayer, 'local_output/farm_multilayer_pos_30_cooc.csv')

