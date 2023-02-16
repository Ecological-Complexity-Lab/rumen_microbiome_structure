# ----- cow_genetic_net.r -----
# Here a genetic network will be build using 
# SNPs data taken from the original paper

# includes ----
library(snpStats)
library(readr)


# consts and functions ----
# Southern breed
bed_h <- "raw_data/Cows_SNPs/Ruminomics_Holstein.bed"
bim_h <- "raw_data/Cows_SNPs/Ruminomics_Holstein.bim"
fam_h <- "raw_data/Cows_SNPs/Ruminomics_Holstein.fam"
# Northern breed
bed_r <- "raw_data/Cows_SNPs/Ruminomics_NordicRed.bed"
bim_r <- "raw_data/Cows_SNPs/Ruminomics_NordicRed.bim"
fam_r <- "raw_data/Cows_SNPs/Ruminomics_NordicRed.fam"

# run -----
raw_h <- read.plink(bed_h, bim_h, fam_h) # 796 cows, 138,892 snps
raw_r <- read.plink(bed_r, bim_r, fam_r) # 199 cows,  76,883 snps

# make sure all the cows have the same list of SNPs 
# (Holstein where sampled with a chip that can detect more types of SNPs. 
#  we need to remove those extra SNPs to compare the breeds)
intr <- intersect(colnames(raw_h$genotypes), 
                  colnames(raw_r$genotypes))

# filtering only the SNPs that exist in both breed data sets
SNPs_h <- raw_h$genotypes[, intr]
SNPs_r <- raw_r$genotypes[, intr]

# as string matrixes:
SNPs_h <- as(SNPs_h, 'character')
SNPs_r <- as(SNPs_r, 'character')

all_cow <- rbind.data.frame(SNPs_h, SNPs_r)


saveRDS(all_cow, "local_output/cows_SNPs_filtered.rds")

all_cow <- readRDS("local_output/cows_SNPs_filtered.rds")
write_csv(all_cow, 'local_output/cows_SNPs_filtered.csv') # this file is big so will not be pushed to github.

# TODO remove, so this is temp.
write.csv(as.data.frame(SNPs_r[1:100, 1:1000]), '~/Desktop/nordic_cows.csv')



# TODO now we have the matrix in a genotype manner. From here calculate the similarity between cows

# hopefully use Noa's code here to construct a similarity matrix, or an edgelist

# TODO Here we have an edgelist using cow IDs


# Get Cow-farm correct labeling
cowdata <- readxl::read_excel('raw_data/RuminOmics_Animal_Phenotypes_for_Mizrahi_v2_plus_rt_quantification_with_total_20170921_and_depth.xlsx', sheet = 3)
cowdata <- cowdata %>%
  select(`Cow ID`, `Farm/Research site code`, `Cow Code`) %>%   # Select only relevant columns
  drop_na()    # Remove all rows with NAs
colnames(cowdata) <- c('ID', 'Farm', 'Cow_Code')
cowdata %<>% 
  mutate(Farm=replace(Farm, Farm=='NUDC', 'UK1')) %>% 
  mutate(Farm=replace(Farm, Farm=='Park', 'UK2')) %>% 
  mutate(Farm=replace(Farm, Farm=='Bianchini', 'IT1')) %>% 
  mutate(Farm=replace(Farm, Farm=='Franciosi', 'IT2')) %>% 
  mutate(Farm=replace(Farm, Farm=='Gandolfi', 'IT3')) %>%
  mutate(Farm=replace(Farm, Farm=='MinkiÃ¶', 'FI1')) %>% 
  mutate(Farm=replace(Farm, Farm=='RÃ¶bÃ¤cksdalen', 'SE1'))

# Join cow labeling to the node or edges so we can connect cow to farm

