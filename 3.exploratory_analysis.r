#------ exploratory_analysis.r -------------------------------------
# This script explores the filtered data in broad strokes - 
# It finds kingdom distribution, quantifies abundance and richness, 
# and calculates alpha and beta diversity.
#------------------------------------------------------------

## @knitr include
# Includes ------------
library(tidyverse)
library(magrittr)
library(cowplot)
library(reshape2) 
library(vegan)
library(igraph)
library(GUniFrac)
library(infomapecology)
check_infomap()
source('functions.R')

# Load data --------------------------------
ASV_Core_30 <- read_csv('local_output/core_ASV_30.csv') %>% 
  mutate(Farm=factor(Farm, levels = c("UK1","UK2","IT1","IT2","IT3","FI1",'SE1')))

## @knitr END

# ASV analysis   --------------------------------

## @knitr ASV_analysis

## Richness per cow  --------------------------------
plt_richness_per_cow <- 
ASV_Core_30 %>%
  group_by(Cow_Code) %>%
  summarise(Richness=n_distinct(ASV_ID)) %>%
  arrange(desc(Richness)) %>% 
  ggplot(aes(Richness))+
  geom_histogram(fill='dark green', color='white')+
  labs(x='ASVs per cow', y='Count') +
  html_figs_theme_no_legend

## Richness per farm  --------------------------------
richness_per_farm <- ASV_Core_30 %>%
  group_by(Farm) %>%
  summarise(ASV_Richness=n_distinct(ASV_ID)) 

plt_richness_per_farm <- 
  richness_per_farm %>%
  arrange(desc(ASV_Richness)) %>% 
  ggplot(aes(Farm,ASV_Richness))+
  geom_col(fill='dark red', color='white')+
  theme_classic() +
  labs(y='Number of ASVs') +
  html_figs_theme_no_legend

png(filename = 'local_output/figures/ASV_richness.png', width = 2300, height = 1000, res = 300)
plot_grid(plt_richness_per_cow, plt_richness_per_farm, 
          nrow = 1, ncol = 2, labels = c('(A)','(B)'),vjust = 1.1)
dev.off()


## Richness per cow per farm   --------------------------------
Richness_per_cow_farm <- 
ASV_Core_30 %>%
  group_by(Farm,Cow_Code) %>%
  summarise(S_cow=n_distinct(ASV_ID)) %>% 
  group_by(Farm) %>% 
  summarise(S_mean=round(mean(S_cow),1),
            S_median=round(median(S_cow),1),
            S_sd=round(sd(S_cow),1),
            S_min=min(S_cow),
            S_max=max(S_cow)) %>% 
  mutate(ASV_summary=paste(S_mean,' [',S_min,'-',S_max,']',sep=''))

# cow stats per all cows
ASV_Core_30 %>%
  group_by(Farm,Cow_Code) %>%
  summarise(S_cow=n_distinct(ASV_ID)) %>% 
  group_by(Farm) %>% 
  summarise(S_mean=round(mean(S_cow),1),
            S_median=round(median(S_cow),1),
            S_sd=round(sd(S_cow),1),
            S_min=min(S_cow),
            S_max=max(S_cow)) %>% 
  mutate(ASV_summary=paste(S_mean,' [',S_min,'-',S_max,']',sep=''))

plt_richness_per_cow_farm <- 
ASV_Core_30 %>%
  group_by(Country,Farm,Cow_Code) %>%
  summarise(richness=n_distinct(ASV_ID)) %>%
  arrange(desc(richness)) %>% 
  ggplot(aes(x=Farm, y=richness, Country))+
  geom_boxplot(aes(color=Country))+
  theme_bw() +
  labs(y = 'ASV richness per cow per farm') +
  paper_figs_theme_no_legend
## @knitr END

## Number of cows in which microbes occur -----------
p1=ASV_Core_30 %>%
  group_by(ASV_ID) %>%
  summarise(habitats=n_distinct(Cow_Code)) %>%
  ggplot(aes(habitats))+
  geom_histogram(fill='brown', color='white')+
  labs(x='Number of cows', y='Count')+
  paper_figs_theme_no_legend

## Number of farms in which microbes occur -----------
p2=ASV_Core_30 %>%
  group_by(ASV_ID) %>%
  summarise(habitats=n_distinct(Farm)) %>%
  arrange(desc(habitats)) %>% 
  group_by(habitats) %>% 
  summarise(n=n_distinct(ASV_ID)) %>% 
  ggplot(aes(habitats, n))+
  geom_col(fill='dark green')+
  labs(x='Number of farms', y='Count')+
  scale_x_continuous(breaks = seq(0,7,1))+
  paper_figs_theme_no_legend

pdf(paste(paper_output_path, "microbe_occurrence.pdf", sep = ""), 10, 6)
plot_grid(p1,p2,nrow = 1, ncol = 2, labels = c('(A)','(B)'))
dev.off()

# Per-farm network stats --------------

## Number of cows per farm -------
cows_per_farm <- ASV_Core_30 %>%
  group_by(Farm) %>%
  summarise(cow_num=n_distinct(Cow_Code))

cows_per_farm %>%
  ggplot(aes(Farm,cow_num))+
  geom_col(fill='dark green', color='white')+
  theme_bw() +
  labs(y = 'Number of cows per farm') +
  theme(axis.text = element_text(size=10, color='black'), 
        title = element_text(size=14), axis.text.x = element_text(angle = 90))


# Save ASV Summary table for the farms ----
# This it part of Table 1 in the paper
# the rest will be completed after constructing the networks on the HPC
# on scripts "cooccurrence_network_analysis.r"
cows_per_farm %>%
  left_join(richness_per_farm) %>%
  left_join(Richness_per_cow_farm %>% select(Farm, ASV_summary)) %>% 
  write_csv('local_output/summary_table_pre_network.csv')


# Cows summary across farms ------
ASV_Core_30 %>%
  group_by(Cow_Code) %>%
  summarise(S_cow=n_distinct(ASV_ID)) %>% 
  summarise(S_mean=mean(S_cow),
            S_median=median(S_cow),
            S_sd=sd(S_cow),
            S_min=min(S_cow),
            S_max=max(S_cow)) %>% 
  mutate(ASV_summary=paste(S_mean,' [',S_min,'-',S_max,']',sep=''))

## @knitr betadiv
# ASV Beta diversity between farms ----------------
ASV_occurrence_farm <- 
  ASV_Core_30 %>%
  select(-c(Cow_Code)) %>%
  distinct(Farm, ASV_ID) %>%
  mutate(present=1) %>% 
  spread(ASV_ID, present, fill = 0) %>%
  column_to_rownames("Farm")

## Using Jaccard ----
beta_diver_farms <- as.matrix(1-vegdist(ASV_occurrence_farm, "jaccard"))
diag(beta_diver_farms) <- 1
# Heatmap
beta_diver_farms[lower.tri(beta_diver_farms, diag = F)] <- NA
beta_diver_farms_m <- reshape2::melt(beta_diver_farms) 

plt_J <- ggplot(beta_diver_farms_m, aes(x = Var1, y = Var2, fill = value, label=round(value,2))) +
  geom_tile() + 
  geom_text(color = "black", size = 4)+
  scale_y_discrete(limits=rev)+
  scale_fill_gradient(high = "blue", low = "light blue", na.value = 'white') +
  html_figs_theme_no_legend+theme(axis.title = element_blank())

## Using UniFrac ----
phylo_tree <- readRDS("local_output/fitted_asvs_phylo_tree.rds")
tree <- phylo_tree$tree
# prune the tree
included_asvs <- unique(ASV_Core_30$ASV_ID)
unincluded <- tree$tip.label[!tree$tip.label %in% included_asvs]
pruned <- dendextend::prune(tree, unincluded)
unifracs <- GUniFrac(ASV_occurrence_farm, pruned, alpha=c(0, 0.5, 1))$unifracs
d_UW <- 1-(unifracs[, , "d_UW"])
# Heatmap 
d_UW[lower.tri(d_UW)] <- NA
beta_diver_farms_UF <- reshape2::melt(d_UW)
plt_U <- 
ggplot(beta_diver_farms_UF, aes(x = Var1, y = Var2, fill = value, label=round(value,2))) +
  geom_tile() + 
  geom_text(color = "black", size = 4)+
  scale_y_discrete(limits=rev)+
  scale_fill_gradient(high = "#ff8c00", low = "#fff494", na.value = 'white') +
  html_figs_theme_no_legend+theme(axis.title = element_blank())

plt_beta_div <- plot_grid(plt_J,plt_U, labels = c('(A)','(B)'))

## @knitr END

# Taxonomic analysis ------------------
ASV_taxa <- read_csv('local_output/ASV_full_taxa.csv') %>% 
  select(ASV_ID, everything(), -seq16S)

n_distinct(ASV_Core_30$ASV_ID)


## Phylum-level composition --------------------------------------------

ASV_Core_30 %>% 
  left_join(ASV_taxa) %>%
  group_by(Farm, Phylum) %>% 
  summarise(ASV_num=n_distinct(ASV_ID)) %>%
  drop_na() %>% 
  mutate(relative_richness=ASV_num/sum(ASV_num)) %>% 
  mutate(ypos = cumsum(relative_richness)- 0.5*relative_richness) %>% 
  ggplot(aes(x="", y=relative_richness, fill=Phylum))+
  facet_wrap(~Farm)+
  geom_bar(stat="identity", width=1) +
  # scale_fill_manual(values = signif_colors)+
  # geom_text(aes(y = ypos, label = round(prop,2)), color = "white", size=3) +
  coord_polar("y", start=0)+
  paper_figs_theme+
  theme(axis.text = element_blank(),
        axis.title = element_blank())

## composition in all ASVs --------------------------------------------
### Phylum-level
percents_ph <- ASV_Core_30 %>% 
  left_join(ASV_taxa) %>%
  group_by(Phylum) %>% 
  summarise(ASV_num=n_distinct(ASV_ID)) %>%
  drop_na() %>% 
  mutate(relative_richness=ASV_num/sum(ASV_num)) %>% 
  mutate(percnt = 100*ASV_num/sum(ASV_num))

### Family-level
percents_fa <- ASV_Core_30 %>% 
  left_join(ASV_taxa) %>%
  group_by(Family) %>% 
  summarise(ASV_num=n_distinct(ASV_ID)) %>%
  drop_na() %>% 
  mutate(relative_richness=ASV_num/sum(ASV_num)) %>% 
  mutate(percnt = 100*ASV_num/sum(ASV_num))


# For paper----
## Number of cows in which microbes occur -----------
ASV_Core_30 %>%
  group_by(ASV_ID) %>%
  summarise(habitats=n_distinct(Cow_Code)) %>%
  ggplot(aes(habitats))+
  geom_histogram(fill='brown', color='white')+
  labs(x='Number of cows', y='Count')+
  html_figs_theme_no_legend

## Number of farms in which microbes occur -----------
ASV_Core_30 %>%
  group_by(ASV_ID) %>%
  summarise(habitats=n_distinct(Farm)) %>%
  arrange(desc(habitats)) %>% 
  group_by(habitats) %>% 
  summarise(n=n_distinct(ASV_ID)) %>% 
  ggplot(aes(habitats, n))+
  geom_col(fill='dark green')+
  labs(x='Number of farms', y='Count')+
  scale_x_continuous(breaks = seq(0,7,1))+
  paper_figs_theme_no_legend

