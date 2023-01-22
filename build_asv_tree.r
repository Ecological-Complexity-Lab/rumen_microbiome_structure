#------ build_ASV_tree.r -------------------------------------
# This script analyses dna sequance of the ASVs to create 
# a phylogenetic tree, and save it to be used in the future

# takes around 2 days to run
#------------------------------------------------------------

#------ includes ------------
library(tidyverse)
library(ape)
library(adegenet)
library(phangorn)

#------ run -------------------------------
#read dna sequences, and ready the data to build the tree
sequences <- read_csv("local_output/ASV_full_taxa.csv")
asv_names <- c(sequences[,1])$ASV_ID

to_convert <- t(sequences[,1:2] %>% column_to_rownames(var = "ASV_ID"))
rownames(to_convert) <- NULL
to_convert <- as.list(to_convert)
names(to_convert) <- asv_names

to_convert <- lapply(to_convert, strsplit, split = "" )
to_convert <- lapply(to_convert, unlist)

# this function receives a list of vectors. each vector contains chars of the sequence of 
binary_for_tree <- as.DNAbin(to_convert)
distences <- dist.dna(binary_for_tree, model="TN93") # long. some are so different that they return a NaN

sum(is.nan(distences)) # how many turned out not to be NaN (distance couldnt be calculated)

# build the tree
tre.ini <- njs(distences) # long

# save it to a file
saveRDS(tre.ini, file = "local_output/njs_distances_for_all_asvs_phylo.rds")

# optimize the initial tree built --------
# read initial tree (so this section will be independent)
tre.ini <- readRDS("local_output/njs_distances_for_all_asvs_phylo.rds")

# what is the likelihood of the data given the model
dna2 <- as.phyDat(binary_for_tree)
fit.ini <- pml(tre.ini, dna2, k=4)

# optimize the tree 
fit <- optim.pml(fit.ini, optNni=TRUE, optBf=TRUE, optQ=TRUE, optGamma=TRUE) # long

# check if the change is significant using anova - it is (has "***" in the results)
anova(fit.ini, fit)

# save the optimized tree. this is a pml object. to get the tree call: fit$tree
saveRDS(fit, file = "local_output/fitted_asvs_phylo_tree.rds")

# lower AIC value is better
AIC(fit.ini) # 822904.5
AIC(fit) # 737123


# Build distance matrix from the build tree ------

# get distance matrix
fit <- readRDS("local_output/fitted_asvs_phylo_tree.rds")
PatristicDistMatrix <- cophenetic.phylo(fit$tree) # long step
saveRDS(PatristicDistMatrix, "local_output/distance_mat_from_tree.rds")

# root the tree using a distinct phylom ----
ASV_taxa <- read_csv('local_output/ASV_full_taxa.csv') %>% select(ASV_ID, everything(), -seq16S)
fit <- readRDS("local_output/fitted_asvs_phylo_tree.rds")

phy <- ASV_taxa %>% filter(Phylum == "Parabasalia")
rooted_tree <- ape::root(fit$tree, outgroup=phy$ASV_ID, resolve.root=TRUE)
is.rooted(rooted_tree)

saveRDS(rooted_tree, file = "local_output/rooted_phylo_tree.rds")

# ideal outgroup is:
# phy <- ASV_taxa %>% filter(Kingdom == "Eukaryota" | Kingdom == "Archaea")
# but it is not monophyletic in our tree..
