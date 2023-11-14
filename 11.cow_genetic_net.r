# ----- cow_genetic_net.r -----
# Here a genetic cow network will be build using 
# SNPs data taken from the original paper

# includes ----
library(snpStats)
library(readr)
library(tidyverse)
library(dendextend)
library(ape)

# run -----

# Get list of cows that appear in both data sets -----
# compare the set of cows with SNPs to the set with microbiome data

## SNP cows ----
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
a <- read.delim('cows_genetic_results/list cows in combined dataset.txt', header = FALSE, sep = " ")

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

# ^ This list was taken to be the basis for calculating the genetic similarity by keren's code.

# visualize the genetic similarity results produced by Keren ------
library(ggtreeExtra)
library(ggtree)
library(ggplot2)
library(ggnewscale)
library(treeio)
library(tidytree)
library(dplyr)
library(ggstar)

# read data
res <- read.csv('cows_genetic_results/genmb_similarity_matrix_weighted.csv', header = TRUE, row.names = 1)
phy <- as.dist(1-res) %>% # convert *similarity* matrix to a *distance* object 
  hclust(method = "ward.D2") %>% as.phylo()
tree <- treeio::as.treedata(phy)

# read metadata
lbls <- read.csv('local_output/cow_list_locations.csv', row.names = 1) %>% 
  select(cow_id=V1, Farm, Country) %>% 
  mutate(breed=case_when(Country %in% c("FI", "SE") ~ "Nordic",
                         !Country %in% c("FI", "SE") ~ "Holstein"))
lbls <- tibble(lbls)

# this is used to parse the cow_id to be able to join with the metadata
get_id <- function(strr) {
  both <- strsplit(strr,"_")
  # sanity check
  if (both[[1]][1] != both[[1]][2]) stop('The ID string is now a double cow ID.')
  return(as.integer(both[[1]][1]))
}

labels <- unlist(lapply(rownames(res), get_id))
dnd_lbls <- tibble(cow_id=labels) %>% left_join(lbls, by="cow_id") 
dnd_lbls$id <- paste(dnd_lbls$cow_id, "_", dnd_lbls$cow_id, sep = "")

breedcolors <- dnd_lbls %>%
  select(c("breed")) %>%
  distinct()
breedcolors$colors <- c("navyblue", "#00B6EB")

#doing this manually because we need a specific order in the colors and legend
countrycolors <- tibble(1:4) 
countrycolors$Country <- c("FI", "SE", "UK", "IT")
countrycolors$colors <- c("#f0a4ff", "#7d37be", "#07ba1b", "darkgreen")

farmcolors <- tibble(1:7)
farmcolors$Farm <- c("FI1", "SE1", "UK1", "UK2", "IT1", "IT2", "IT3")
farmcolors$colors <- c("#fff45d", "#fbcd4f", "#f6a541", "#f27e33", "#ed5724", "#e92f16", "#C80707")

dnd_lbls$breed   <- factor(dnd_lbls$breed, levels=breedcolors$breed)
dnd_lbls$Country <- factor(dnd_lbls$Country, levels=countrycolors$Country)
dnd_lbls$Farm    <-  factor(dnd_lbls$Farm, levels=farmcolors$Farm)

meta <- dnd_lbls %>% select(id, Breed = breed, Country, Farm)

# start building the tree
p <- ggtree(tree, layout="fan", open.angle=5, size=0.2, branch.length = "none")
p <- p %<+% meta

# adding the breed
p1 <-p +
  geom_fruit(geom=geom_tile,
             mapping=aes(fill=Breed),
             width=1.8,
             offset=0.05) +
  scale_fill_manual(name="Breed",
                    values=breedcolors$colors,
                    guide=guide_legend(keywidth=0.3, keyheight=0.5, ncol=1, order=2)) +
  theme(legend.title=element_text(size=10), 
        legend.text=element_text(size=8),
        legend.spacing.y = unit(0.05, "cm"))
p2 <-p1 +
  new_scale_fill() +
  geom_fruit(geom=geom_tile,
             mapping=aes(fill=Country),
             width=1.8,
             offset=0.08) +
  scale_fill_manual(name="Country",
                    values=countrycolors$colors,
                    guide=guide_legend(keywidth=0.3, keyheight=0.5, ncol=1, order=2)) 
p3 <-p2 +
  new_scale_fill() +
  geom_fruit(geom=geom_tile,
             mapping=aes(fill=Farm),
             width=1.8,
             offset=0.08) +
  scale_fill_manual(name="Farm",
                    values=farmcolors$colors,
                    guide=guide_legend(keywidth=0.3, keyheight=0.5, ncol=1, order=2)) 
p3

pdf("local_output/figures/cow_genetics.pdf", 6, 6)
p3
dev.off()

# Analyse genetic tree structure ----
# explore the tree structure in the first two levels, regarding annotations
x <- as_tibble(tree)

y <- full_join(x, dnd_lbls, by = c('label' = 'id'))
new_tree <- as.treedata(y)

y # node number 1 is a leaf with a label (an actual cow)
ancestor(y, 1) # from this we learn node *936* is the top node
child(y, 936) # level 1 division is between nodes 937 and 938


## Level 1 analysis: ----
cows_north <- offspring(y, 937) %>% filter(!is.na(label)) # Northern cluster
cows_south <- offspring(y, 938) %>% filter(!is.na(label)) # Southern cluster

nrow(cows_north) # 194 cows in the northern cluster
nrow(cows_south) # 741 cows in the south

table(cows_north$breed) # only Nordic Red cows in the south
table(cows_south$breed) # only Holstein cows in the south

table(cows_north$Country) # FI=100, SE=94
table(cows_south$Country) # UK=368, IT=373


## Level 2 analysis: ----

# Northern cluster:
child(y, 937) # level 2 - north division is between nodes 949 and 950

nord_1 <- offspring(y, 949) %>% filter(!is.na(label)) # Nord-1
nord_2 <- offspring(y, 950) %>% filter(!is.na(label)) # Nord-2

nrow(nord_1) # 141 cows
nrow(nord_2) # 53  cows

table(nord_1$Country) # FI=95, SE=46
table(nord_2$Country) # FI=5,  SE=48


# Southern cluster:
child(y, 938) # level 2 - south division is between nodes 939 and 940

south_1 <- offspring(y, 939) %>% filter(!is.na(label)) # south-1
south_2 <- offspring(y, 940) %>% filter(!is.na(label)) # south-2

nrow(south_1) # 28 cows
nrow(south_2) # 713  cows

table(south_1$Country) # UK=28,  IT=0
table(south_2$Country) # UK=368, IT=345


