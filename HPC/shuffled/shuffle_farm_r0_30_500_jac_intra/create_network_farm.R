#! /gpfs0/shai/projects/R4/R-4.0.3/bin/Rscript
.libPaths("/gpfs0/shai/projects/R4/R-4.0.3/lib64/R/library")
print(.libPaths())
print(sessionInfo())

library(tidyverse)
library(magrittr)
library(reshape2) 
library(vegan) 
library(cooccur)
library(igraph)
library(dplyr)

JOB_ID <- Sys.getenv("JOB_ID")
if (length(commandArgs(trailingOnly=TRUE))==0) {
  stop('No arguments were found!')
} else {
  args <- commandArgs(trailingOnly=TRUE)
  exp_id <- as.numeric(args[1])
  curr_farm <- as.character(args[2])
}

source('functions.R')

ASV_Core_file <- read_csv('experiments.csv') %>% filter(e_id==exp_id)
ASV_Core_file <- ASV_Core_file$data_file
is_shuff <- ifelse(grepl("shuff", ASV_Core_file, ignore.case = TRUE), TRUE, FALSE)

write_to_log(paste0('START EVENT LOG: EXPERIMENT ',exp_id), exp_id = exp_id, JOB_ID = JOB_ID, level = 'Farm', level_name = curr_farm, append = F)
write_to_log(paste('Data file: ',ASV_Core_file), exp_id = exp_id, JOB_ID = JOB_ID, level = 'Farm', level_name = curr_farm)


# Coocurrence network for country scale -------------------------------------------

write_to_log(paste('========= FARM LEVEL: ',curr_farm,'  ========='), exp_id = exp_id, JOB_ID = JOB_ID, level = 'Farm', level_name = curr_farm)

ASV_data <- read_csv(ASV_Core_file) %>% filter(Farm==curr_farm)

# Create the network. The distinct is not needed per se because cows are
# individuals so ASV-cow pairs must be unique anyways. But it helps making the data binary
ASV_cow <- ASV_data %>% distinct(Cow_Code,ASV_ID) %>% mutate(weight=1)

# Create the bipartite matrix
ASV_occurrence_mat_cow <- acast(ASV_cow, ASV_ID~Cow_Code, fill = 0)

write_to_log('--- ASV_occurrence_mat_cow ---', exp_id = exp_id, JOB_ID = JOB_ID, level = 'Farm', level_name = curr_farm)
write_to_log(dim(ASV_occurrence_mat_cow), exp_id = exp_id, JOB_ID = JOB_ID, level = 'Farm', level_name = curr_farm)
write_to_log(density_bip(ASV_occurrence_mat_cow), exp_id = exp_id, JOB_ID = JOB_ID, level = 'Farm', level_name = curr_farm)

# Are all ASVs in the matrix?
x <-  sort(unique(ASV_cow$ASV_ID))
y <- rownames(ASV_occurrence_mat_cow)
write_to_log(paste('Are all ASVs in the matrix?',setequal(x,y)), exp_id = exp_id, JOB_ID = JOB_ID, level = 'Farm', level_name = curr_farm)
# Are all cows in the matrix?
x <-  sort(unique(ASV_cow$Cow_Code))
y <- colnames(ASV_occurrence_mat_cow)
write_to_log(paste('Are all cows in the matrix?',setequal(x,y)), exp_id = exp_id, JOB_ID = JOB_ID, level = 'Farm', level_name = curr_farm)

if (!is_shuff) {
  # Write the matrix
  write.csv(ASV_occurrence_mat_cow, paste(exp_id,JOB_ID,'Farm',curr_farm,"ASV_cow_mat.csv",sep = "_"))
}

# Find significant pairwise co-occurrences.
ASV_co <- cooccur(ASV_occurrence_mat_cow, spp_names = TRUE)
#ez <- effect.sizes(mod=ASV_co, standardized = TRUE)

ASV_co <- ASV_co$results
# add jaccard values, and using them as weight:
ASV_co %<>% mutate(weight=obs_cooccur/(sp1_inc+sp2_inc-obs_cooccur))
#ASV_co <- left_join(ASV_co, ez, by=c("sp1_name" = "sp1", "sp2_name" = "sp2"))

# Do all the ASVs appear?
x <- sort(unique(ASV_cow$ASV_ID))
y <- sort(unique(c(ASV_co$sp1_name,ASV_co$sp2_name)))
write_to_log(paste('Do all the ASVs appear?',setequal(x,y)), exp_id = exp_id, JOB_ID = JOB_ID, level = 'Farm', level_name = curr_farm)

# Create the co-occurrence network
ASV_co %<>%
  as_tibble() %>% 
  mutate(level='Farm') %>% 
  mutate(level_name=curr_farm) %>% 
  select(sp1,sp2, weight, everything()) %>% 
  mutate(edge_type=case_when(
    p_lt < 0.05 & p_gt >= 0.05 ~ "neg",
    p_lt >= 0.05 & p_gt < 0.05 ~ "pos",
    p_lt >= 0.05 & p_gt >= 0.05 ~ "not_significant"
  ))

write.csv(ASV_co, paste(exp_id,JOB_ID,'Farm',curr_farm,"COOC.csv",sep = "_"))


# Parse the coocurrence to a network --------------------------------------

write_to_log('Parsing networks', exp_id = exp_id, JOB_ID = JOB_ID, level = 'Farm', level_name = curr_farm)

net <- parse_networks_HPC()
write_to_log('Writing files', exp_id = exp_id, JOB_ID = JOB_ID, level = 'Farm', level_name = curr_farm)
if (!is_shuff) {
  write_csv(net$nodes, paste(exp_id,JOB_ID,'Farm',curr_farm,"nodes.csv",sep = "_"))
  write_csv(as.data.frame(net$mat_pos), paste(exp_id,JOB_ID,'Farm',curr_farm,"mat_pos.csv",sep = "_"))
  write_csv(as.data.frame(net$mat_neg), paste(exp_id,JOB_ID,'Farm',curr_farm,"mat_neg.csv",sep = "_"))
  write_csv(net$singletons, paste(exp_id,JOB_ID,'Farm',curr_farm,"singletons.csv",sep = "_"))
}
net$edge_list$level='Farm'
net$edge_list$level_name=curr_farm
write_csv(net$edge_list, paste(exp_id,JOB_ID,'Farm',curr_farm,'edge_list.csv',sep = "_"))

# Write summary when done ------------------
summary_file <- 'run_summary.csv'
# if data is shuffled save summery one directory up
if (is_shuff) { 
  summary_file <- '../run_summary.csv' 
} 

write_csv(tibble(exp_id=exp_id, Level='Farm', Level_name=curr_farm, 
                 JOB_ID=JOB_ID, data_file=ASV_Core_file, time_stamp=Sys.time()),
          summary_file, append = T)
