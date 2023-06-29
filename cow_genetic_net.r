# ----- cow_genetic_net.r -----
# Here a genetic network will be build using 
# SNPs data taken from the original paper

# includes ----
library(snpStats)
library(readr)
library(tidyverse)
library(dendextend)


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


# Get list of cows that appear in both data sets -----
# compare the set of cows with SNPs to the set with microbiome data

## SNP ----
# read SNP cows ids and combine to one df
nord_snp <- read_delim("raw_data/Cows_SNPs/Ruminomics_NordicRed.fam", 
                       skip_empty_rows = TRUE,col_names = FALSE)
hols_snp <- read_delim("raw_data/Cows_SNPs/Ruminomics_Holstein.fam", 
                     skip_empty_rows = TRUE,col_names = FALSE)
nord_snp$ID <- as.numeric(nord_snp$X1)
hols_snp$ID <- as.numeric(hols_snp$X1)

ambas <- rbind(nord_snp, hols_snp) 
nrow(ambas)

min(ambas$ID)
max(ambas$ID)

## Microbiome cows -----
# read microbiome cows ids - all cow in our data
# read only id of cows used in the analysis:
ASV_Core_30 <- read_csv('local_output/core_ASV_30.csv') %>% 
  mutate(Farm=factor(Farm, levels = c("UK1","UK2","IT1","IT2","IT3","FI1",'SE1')))

microbs <- ASV_Core_30 %>%
  group_by(Cow_Code) %>%
  summarise(microbs=n()) %>%
  arrange(Cow_Code) %>%
  separate(Cow_Code, into = c("country", "ID"), c(2))%>%
  mutate(breed=case_when(country %in% c("FI", "SE") ~ "Nordic",
                         !country %in% c("FI", "SE") ~ "Holstein"))
microbs$ID <- as.numeric(microbs$ID)
nrow(microbs)

min(microbs$ID)
max(microbs$ID)

table(microbs$breed)
# Holstein   Nordic 
# 816      196 

# compare all cows in SNP vs microbiome
micr <- microbs$ID

length(intersect(snp, micr))
setdiff(snp, micr)
setdiff(micr, snp)

## compare datasets for each breed separately -----
# Nordics:
nord_micr <- microbs[microbs$breed == "Nordic",]

# differences
nrow(nord_snp) # SNPs
nrow(nord_micr) # Microbiome

setdiff(nord_snp$ID, nord_micr$ID)
setdiff(nord_micr$ID, nord_snp$ID)
length(intersect(nord_snp$ID, nord_micr$ID))

# Holstein:
hols_micr <- microbs[microbs$breed == "Holstein",]

# differences
nrow(hols_snp) # SNPs
nrow(hols_micr) # microbiome

setdiff(hols_snp$ID, hols_micr$ID)
setdiff(hols_micr$ID, hols_snp$ID)
length(intersect(hols_snp$ID, hols_micr$ID))

# Save intersection cows -----
# output intersections cow ids between SNPs and microbiome data
nord_intr_df <- data.frame(cow_id=intersect(nord_snp$ID, nord_micr$ID),
                           breed="NordicRed")
hols_intr_df <- data.frame(cow_id=intersect(hols_snp$ID, hols_micr$ID),
                           breed="Holstein")
length(intersect(nord_snp$ID, nord_micr$ID))
length(intersect(hols_snp$ID, hols_micr$ID))

# save cow intersection
write.csv(rbind(hols_intr_df, nord_intr_df), 'local_output/SNP_micro_intersect_cows.csv')

# produce the combined cow list with country and farm data -----
a <- read.delim('list cows in combined dataset.txt', header = FALSE, sep = " ")

# Get Cow-farm correct labeling
cowdata <- readxl::read_excel('raw_data/RuminOmics_Animal_Phenotypes_for_Mizrahi_v2_plus_rt_quantification_with_total_20170921_and_depth.xlsx', sheet = 3)
cowdata <- cowdata %>%
  select(`Cow ID`, `Farm/Research site code`, `Cow Code`, Country) %>%   # Select only relevant columns
  drop_na()    # Remove all rows with NAs
colnames(cowdata) <- c('ID', 'Farm', 'Cow_Code', 'Country')
cowdata %<>% 
  mutate(Farm=replace(Farm, Farm=='NUDC', 'UK1')) %>% 
  mutate(Farm=replace(Farm, Farm=='Park', 'UK2')) %>% 
  mutate(Farm=replace(Farm, Farm=='Bianchini', 'IT1')) %>% 
  mutate(Farm=replace(Farm, Farm=='Franciosi', 'IT2')) %>% 
  mutate(Farm=replace(Farm, Farm=='Gandolfi', 'IT3')) %>%
  mutate(Farm=replace(Farm, Farm=='MinkiÃ¶', 'FI1')) %>% 
  mutate(Farm=replace(Farm, Farm=='RÃ¶bÃ¤cksdalen', 'SE1'))


both <- a %>% left_join(cowdata, by=c('V1' = 'ID'))

write.csv(both, 'local_output/cow_list_locations.csv')


# visualize the genetic similarity results produced by Keren ------
res <- read.csv('./cows_genetic_results/genmb_similarity_matrix_weighted.csv', header = TRUE, row.names = 1)
lbls <- read.csv('local_output/cow_list_locations.csv', row.names = 1) %>% 
          select(cow_id=V1, Farm, Country) %>% 
          mutate(breed=case_when(Country %in% c("FI", "SE") ~ "Nordic",
                                 !Country %in% c("FI", "SE") ~ "Holstein"))
lbls <- tibble(lbls)

# make sure it is symmetric
res1 <- res
rr <- paste("X", rownames(res1), sep = "")
rownames(res1) <- rr
isSymmetric(as.matrix(res1))

get_id <- function(strr) {
  both <- strsplit(strr,"_")
  # sanity check
  if (both[[1]][1] != both[[1]][2]) stop('The ID string is now a double cow ID.')
  return(as.integer(both[[1]][1]))
}

labels <- unlist(lapply(rownames(res), get_id))
dnd_lbls <- tibble(cow_id=labels) %>% 
          left_join(lbls, by="cow_id") 
# prepare to make colors:
dnd_lbls2 <- dnd_lbls %>% 
          mutate_at(2, as_factor)%>%  mutate_at(2, as.numeric) %>%
          mutate_at(3, as_factor)%>%  mutate_at(3, as.numeric) %>%
          mutate_at(4, as_factor) %>%  mutate_at(4, as.numeric)

# create a dendagram
dendi <- as.dist(1-res) %>% # convert *similarity* matrix to a *distance* object 
  hclust(method = "ward.D2") %>% # Hierarchical clustering 
  as.dendrogram() # Turn the object into a dendrogram


pdf("./cow_genetics.pdf", 10, 6)
dendi %>% set("labels", "") %>% plot
colored_bars(colors = cbind(dnd_lbls2$breed, dnd_lbls2$Country, dnd_lbls2$Farm),
             rowLabels = c("Breeds", "Country" , "Farm"))
dev.off()

