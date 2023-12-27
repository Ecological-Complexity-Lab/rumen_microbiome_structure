# ------ 14.Revision_analysis.r ----------------------------------------
# This script temporarily holds the analysis requested in the revision
# they will later be moved to more logically relevant places in the repo
#-----------------------------------------------------------------------

#------ includes ----------
library(tidyverse)
library(magrittr)
library(reshape2)

#------ consts ---------
single_prob_file <- 'local_output/single_asv_occur_prob_80.csv'
combo_prob_file <- 'combo_asv_occur_prob.csv'

#------ run ------------
## Read data to be used -----
# read network:
farm_multilayer_pos_30 <-read_csv('local_output/farm_multilayer_pos_30.csv')
all_nodes <- sort(unique(c(farm_multilayer_pos_30$from, farm_multilayer_pos_30$to)))
all_nodes <- tibble(node_id=1:length(all_nodes), node_name=all_nodes)
layers <- tibble(layer_id=1:7, layer_name=c('NUDC', 'Park', 'Bian', 'Fran','Gand','Mink','Raab'),
                 short_name=c('UK1', 'UK2', 'IT1', 'IT2', 'IT3', 'FI1', 'SE1'))


# read asv data
ASV_data_final <- read_csv("local_output/core_ASV_80.csv")

## Transitivity ---------

### single probabilities ----
# calculate single probability of each ASV for every layer
probabilities <- NULL
for (layer in layers$short_name) {
  # filter out the layer
  curr_layer <- ASV_data_final %>% filter(Farm == layer)
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
  layer_asv <- ASV_data_final %>% filter(Farm == layer)
  layer_net <- intras %>% filter(layer == layer)
  
  # run over results - 
  # TODO some analysis should be done here, this hist is not good enough
  rrr <- read_csv(paste(HPC_output_folder, layer, '_', combo_prob_file, sep = ""), 
                  col_names = c("ASV_1", "ASV_2", "ASV_3", "common_cows", "random", "emp_prob"))
  ns <- nrow(rrr)
  p1 <- hist(rrr$emp_prob)                    
  p2 <- hist(rrr$random)                     
  plot( p2, col=rgb(1,0,0,1/4), xlim=c(0,1), ylim = c(0,ns))
  plot( p1, col=rgb(0,0,1,1/4), xlim=c(0,1), add=T) 
  
}



# phylogeny partner fidelity ----------



# modularity - phylogenetic composition --------
# includes randomality check



# NMI of clusters and hypothesis ---------

