# ----- cow_genetic_net.r -----
# Here a genetic network will be build using 
# SNPs data taken from the original paper

# includes ----
library(snpStats)

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
write_csv(all_cow, 'local_output/cows_SNPs_filtered.csv') # this file is big so will not be pushed to github.

# TODO remove, so this is temp.
write.csv(as.data.frame(SNPs_r[1:100, 1:1000]), '~/Desktop/nordic_cows.csv')



# TODO now we have the matrix in a genotype manner. From here calculate the similarity between cows





