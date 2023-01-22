#------ prepare_raw_data.r -------------------------------------
# This script prepares the data to be analysed from raw data 
# It converts sequences count to an occurrence table
#------------------------------------------------------------

#------ includes ----------
library(tidyverse)
library(magrittr)
library(reshape2)

#------ run ---------------
# Read the phylogenetic classification
phylo <- read_delim('raw_data/full_taxa.tsv', delim='\t')
phylo %<>%
  # Add a unique ID for each sequence
  mutate(ASV_ID=paste('ASV',str_pad(1:nrow(phylo),5,'left','0'),sep='_')) %>% 
  select(ASV_ID, everything()) %>%
  write_csv("local_output/ASV_full_taxa.csv")
# How many distinct sequences do we have?
phylo %>% distinct(seq16S) %>% count() # Check the number of distinct (unique) sequences. Must be the same as the number of rows
# Any duplicates?
any(duplicated(phylo$seq16S))

# Get the ASV per cow data. 1st column is the cow, the rest are the ASVs. Cells
# are ASV abundance in a cow
ASV_data <- read_delim('raw_data/full_table.nochim.txt', delim = ' ')
colnames(ASV_data)[1:3] # Look at the names of the 3 first columns

# Each sample has two reads. Example with cows UK271 and UK272:
x <- ASV_data[c(801,802,1813,1814),c(1,100:105)]
x$sample

# We need to sum the reads. Here is an example with a subset of the data 
x %>%
  separate(sample, sep = c("_"), c("Prok","seq", "sample")) %>%
  select(-seq, -Prok) %>%
  group_by(sample) %>%
  summarise_all(list(sum=sum))

# A loop that writes a section of 400 columns into a file at a time
for (i in seq(2,12002,400)) {     # run this loop on 400 columns every time
  d=399   
  if (i==12002) { d= length(ASV_data)-i}     # because the number of columns is not divisible by 400
  ASV_data_sum_reads <-  ASV_data[,c(1,i:(i+d))] %>%    # always include the 1st column
    mutate(sample=str_remove(sample,'.fastq')) %>%    # remove the fastq from sample name
    separate(sample, sep = c("_"), c("Prok","seq", "sample")) %>% # Separate name parts to sum across reads
    select(-seq, -Prok) %>% 
    # Sumamrise the samples across reads
    group_by(sample) %>%
    summarise_all(list(sum_reads=sum))
  
  write_csv(ASV_data_sum_reads, paste("raw_data/ASV", i, "_", i+d, ".csv", sep = ""))
  rm(ASV_data_sum_reads)    # remove
}

# A loop that reads the files we wrote and binds them all together to one dataset
ASV_data_sum_reads <- read_csv(paste("raw_data/ASV", 2, "_", 2+399, ".csv", sep = ""))[,1]   # read the 1st column (sample)
for (i in seq(2,12002,400)) {   
  d=399
  if (i==12002) { d= length(ASV_data)-i }
  ASV_temp <- read_csv(paste("raw_data/ASV", i, "_", i+d, ".csv", sep = ""))[,-1]  # read each file without the 1st column
  ASV_data_sum_reads <- cbind(ASV_data_sum_reads, ASV_temp)   # bind each file 
  rm(ASV_temp)
}

write_csv(ASV_data_sum_reads, file = "raw_data/all_ASV_reads.csv")
ASV_data_sum_reads <- read_csv("raw_data/all_ASV_reads.csv")

# Look at a subuset of the data
ASV_data_sum_reads[200:201,801:804]  
# select a single sequence to see how the names changed and for further tests below
seq_to_test <- colnames(ASV_data_sum_reads)[3] 

# add E to the Sweden cows that has only S in their name:
ASV_data_sum_reads$sample <- ifelse(str_starts(ASV_data_sum_reads$sample, 'S'), paste("SE", parse_number(ASV_data_sum_reads$sample)), ASV_data_sum_reads$sample)
# 135 cows in SWEDEN (should be 100)

# remove the spaces from the cow codes
ASV_data_sum_reads %<>%    
  mutate(sample=str_replace(sample, " ", "")) %>%
  write_csv("raw_data/all_ASV_sum_reads.csv")
ASV_data_sum_reads <- read_csv("raw_data/all_ASV_sum_reads.csv")

# Now for the real data:
# Transpose to be able to join the ASV data with the phylogenetic table that has unique ASV ids
ASV_data_trans <- ASV_data_sum_reads %>% 
  gather(seq16S, val, 2:ncol(ASV_data_sum_reads)) %>%
  spread(sample, val) %>% 
  # Remove all the suffix "_sum_reads"
  mutate(seq16S=str_remove(seq16S, '_sum_reads'))

# Check to see the transpose worked
subset(ASV_data_trans, seq16S==str_remove(seq_to_test, '_sum_reads'))[,c(1,801:804)]

# Join the long-format ASV data with the tibble that has the unique IDs.

# First take a look at a single sequence to see this will work
ASV_data_trans %>% inner_join(phylo %>% select(ASV_ID, seq16S), by='seq16S') %>% 
  select(ASV_ID, everything()) %>% 
  # Check to see the join worked
  filter(seq16S==str_remove(seq_to_test, '_sum_reads'))

# Now join everything
ASV_data_matrix <- 
  ASV_data_trans %>% inner_join(phylo %>% select(ASV_ID, seq16S), by='seq16S') %>% 
  select(ASV_ID, everything()) %>% 
  # We don't need the sequence itself
  select(-seq16S) %>% 
  arrange(ASV_ID)

cowdata <- readxl::read_excel('raw_data/RuminOmics_Animal_Phenotypes_for_Mizrahi_v2_plus_rt_quantification_with_total_20170921_and_depth.xlsx', sheet = 3)
cowdata[406:411,1:3] # There are row spaces between countries

cowdata_new <- cowdata %>%
  select(Country, `Farm/Research site code`, `Cow Code`) %>%   # Select only relevant columns
  drop_na()    # Remove all rows with NAs
colnames(cowdata_new) <- c('Country', 'Farm', 'Cow_Code')  
n_distinct(cowdata_new$Cow_Code, na.rm = F)  # Check there are no duplicates in Cow Code
# 100 cows in SWEDEN

# convert the ASV_data_matrix to a long format:
# first I need to transpose again to be able to join the ASV data with the cow data that has the cow codes
ASV_data_matrix_trans <- ASV_data_matrix %>%
  gather(Cow_Code, val, 2:ncol(ASV_data_matrix)) %>%
  spread(ASV_ID, val) %>%
  arrange(Cow_Code)
ASV_data_matrix_long <- ASV_data_matrix_trans %>% 
  melt(id.vars = "Cow_Code") 
colnames(ASV_data_matrix_long) <- c("Cow_Code", "ASV_ID", "Abundance")

# join ASV_data_matrix_long with cowdata_new
# ASV_data_final <- inner_join(cowdata_new, ASV_data_matrix_long, by = "Cow_Code") %>%
#   filter(Abundance > 0)   # if the microbe does not appear in this cow, the line is irrelevant to us
# 
# write_csv(ASV_data_final, "output/ASV_processed_data.csv")


# check the cows which are in cowdata_new but not in ASV_data_matrix_long and vice versa
setdiff(cowdata_new$Cow_Code,ASV_data_matrix_long$Cow_Code)
setdiff(ASV_data_matrix_long$Cow_Code,cowdata_new$Cow_Code)

# fix the labels of the cows- from SE to FI in ASV_data_matrix_long
ASV_data_matrix_long$Cow_Code <- ifelse(str_starts(ASV_data_matrix_long$Cow_Code, 'SE9'), paste("FI", parse_number(ASV_data_matrix_long$Cow_Code)), ASV_data_matrix_long$Cow_Code)

# remove the spaces from the cow codes
ASV_data_matrix_long_fixed <- ASV_data_matrix_long %>%
  mutate(Cow_Code=str_replace(Cow_Code, " ", "")) 

# check the cows which are in cowdata_new but not in ASV_data_matrix_long_fixed and vice versa
setdiff(cowdata_new$Cow_Code,ASV_data_matrix_long_fixed$Cow_Code)
setdiff(ASV_data_matrix_long_fixed$Cow_Code,cowdata_new$Cow_Code)

# join ASV_data_matrix_long_fixed with cowdata_new
ASV_data_final_fixed <- inner_join(cowdata_new, ASV_data_matrix_long_fixed, by = "Cow_Code") %>%
 filter(Abundance > 0)   # if the microbe does not appear in this cow, the line is irrelevant to us

n_distinct(ASV_data_final_fixed$ASV_ID)
write_csv(ASV_data_final_fixed, "local_output/ASV_processed_data_fixed.csv")
