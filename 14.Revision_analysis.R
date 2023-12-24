# ------ 14.Revision_analysis.r ----------------------------------------
# This script temporarily holds the analysis requested in the revision
# they will later be moved to more logically relevant places in the repo
#-----------------------------------------------------------------------

#------ includes ----------
library(tidyverse)
library(magrittr)
library(reshape2)

#------ consts ---------
single_prob_file <- 'local_output/single_asv_occur_prob.csv'
combo_prob_file <- 'combo_asv_occur_prob.csv'

#------ functions ---------
# iterate combinations function
get_next_cbn <- function(cbn, n){
  ## Generates the combination that follows the one provided as input, 
  # slower but this way i don't need to save all the combinations
  cbn.bin      <- rep(0, n)
  cbn.bin[cbn] <- 1
  if (tail(cbn.bin, 1) == 0){
    ind <- tail(which(cbn.bin == 1), 1)
    cbn.bin[c(ind, ind+1)] <- c(0, 1)
  }else{
    ind <- 1 + tail(which(diff(cbn.bin) == -1), 1)
    nb  <- sum(cbn.bin[-c(1:ind)] == 1)
    cbn.bin[c(ind-1, (n-nb+1):n)] <- 0
    cbn.bin[ind:(ind+nb)]         <- 1
  }
  cbn <- which(cbn.bin == 1)
}

process_combination <- function(cbn, asv_data, asv_list, single_probs, output) {
  # here:
  # get 3 asvs
  first <- asv_list[cbn[1]]
  second <- asv_list[cbn[2]]
  third <- asv_list[cbn[3]]
  
  # get single probs
  prob1 <- single_probs %>% filter(ASV==first) %>% pull(rand_prob)
  prob2 <- single_probs %>% filter(ASV==second) %>% pull(rand_prob)
  prob3 <- single_probs %>% filter(ASV==third) %>% pull(rand_prob)
  
  # count number of common hosting cows
  host_cows <- count_comon_cows(first, second, third, asv_data)
  # calculate probability to coocur randomly
  rand_prob <- prob1 * prob2 * prob3
  
  # calculate probability
  # save to a line in a tibble 
  new_line <- tibble(ASV_1=first, ASV_2=second, ASV_3=third, 
                     common_cows=host_cows, random=rand_prob)
  output <- rbind(new_line, output)
}

count_comon_cows <- function(ASV_1, ASV_2, ASV_3, asv_data){
  # in how many cows all three exist
  all_three <- asv_data %>% filter((ASV_ID == ASV_1) | 
                                       (ASV_ID == ASV_2) | 
                                       (ASV_ID == ASV_3))
  temp <- table(all_three$Cow_Code)
  output <- sum(temp == 3)
  
  return(output)
}


#------ run ------------
## Read data to be used -----
# read network:
farm_multilayer_pos_30 <-read_csv('local_output/farm_multilayer_pos_30.csv')
all_nodes <- sort(unique(c(farm_multilayer_pos_30$from, farm_multilayer_pos_30$to)))
all_nodes <- tibble(node_id=1:length(all_nodes), node_name=all_nodes)
layers <- tibble(layer_id=1:7, layer_name=c('NUDC', 'Park', 'Bian', 'Fran','Gand','Mink','Raab'),
                 short_name=c('UK1', 'UK2', 'IT1', 'IT2', 'IT3', 'FI1', 'SE1'))

# read cow data
ASV_data_final <- read_csv("local_output/ASV_processed_data.csv")
ASV_data_final %>% group_by(Country,Farm) %>% summarise(cows=n_distinct(Cow_Code))
ASV_data_final %>% group_by(Country,Farm) %>% summarise(ASVs=n_distinct(ASV_ID))

# Match farm names to those in the original paper
ASV_data_final %<>% 
  mutate(Farm=replace(Farm, Farm=='NUDC', 'UK1')) %>% 
  mutate(Farm=replace(Farm, Farm=='Park', 'UK2')) %>% 
  mutate(Farm=replace(Farm, Farm=='Bianchini', 'IT1')) %>% 
  mutate(Farm=replace(Farm, Farm=='Franciosi', 'IT2')) %>% 
  mutate(Farm=replace(Farm, Farm=='Gandolfi', 'IT3')) %>%
  mutate(Farm=replace(Farm, Farm=='MinkiÃ¶', 'FI1')) %>% 
  mutate(Farm=replace(Farm, Farm=='RÃ¶bÃ¤cksdalen', 'SE1'))

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
probabilities <- read_csv(single_prob_file)

## calculate transitivity for each combination in each layer
for (layer in layers$short_name) {
  layer_probs <- probabilities %>% filter(Farm == layer)
  curr_layer <- ASV_data_final %>% filter(Farm == layer)
  cows <- as.numeric(curr_layer %>% summarise(cows=n_distinct(Cow_Code)))
  v <- unique(curr_layer$ASV_ID)

  n   <- length(v)
  k   <- 3
  
  # handle first combination
  ## Iteration example
  cbn <- 1:k
  output <- process_combination(cbn, curr_layer, v, layer_probs, NULL) 
  print(cbn)
  for (i in 2:choose(n, k)){
  #for (i in 2:100){
    cbn <- get_next_cbn(cbn, n)
    print(cbn)
    output <- process_combination(cbn, curr_layer, v, layer_probs, output)
    
    
    if (i %% 1000 == 0){ # make this occasional as i/o actions are costly
      output <-  output %>% mutate(emp_prob=common_cows/cows)
      write_csv(output, combo_prob_file, append = T)
      
      output <- NULL
    }
  }
  
  if (nrow(output) > 0) {
    output <-  output %>% mutate(emp_prob=common_cows/cows)
    write_csv(output, combo_prob_file, append = T)
  }
}



# phylogeny partner fidelity ----------



# modularity - phylogenetic composition --------
# includes randomality check



# NMI of clusters and hypothesis ---------

