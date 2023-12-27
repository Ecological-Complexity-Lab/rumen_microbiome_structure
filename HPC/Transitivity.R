#! /gpfs0/shai/projects/R4/R-4.0.3/bin/Rscript
# ------ Transitivity ----------------------------------------
# This script calculated the probability based Transitivity in a single layer
# they will later be moved to more logically relevant places in the repo
#-----------------------------------------------------------------------
.libPaths("/gpfs0/shai/projects/R4/R-4.0.3/lib64/R/library")
print(.libPaths())
print(sessionInfo())


#------ includes ----------
library(tidyverse)
library(magrittr)
library(reshape2)


#------ Arguments ---------
JOB_ID <- Sys.getenv("JOB_ID")
if (length(commandArgs(trailingOnly=TRUE))==0) {
  stop('No arguments were found!')
} else {
  args <- commandArgs(trailingOnly=TRUE)
  layer <- as.character(args[1])
}

#------ consts ---------
asv_data_file    <- 'core_ASV_80.csv'
single_prob_file <- 'single_asv_occur_prob_80.csv'
combo_prob_file  <- 'combo_asv_occur_prob.csv'
output_file <- paste(layer, combo_prob_file, sep = "_")

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
# read asv data
ASV_data_final <- read_csv(asv_data_file)

### combo probabilities ----
probabilities <- read_csv(single_prob_file)

## calculate transitivity for each combination in each layer
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
    write_csv(output, output_file, append = T)
    
    output <- NULL
  }
}

if (nrow(output) > 0) {
  output <-  output %>% mutate(emp_prob=common_cows/cows)
  write_csv(output, output_file, append = T)
}


