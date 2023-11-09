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

## Number of co-occurrence links per farm -----
farm_multilayer_pos <- read_csv('local_output/farm_multilayer_pos_30.csv')
links_per_farm <- farm_multilayer_pos %>%
  group_by(level_name) %>%
  mutate(link=1) %>%
  summarise(links=sum(link))

links_per_farm %>%
  ggplot(aes(level_name, links)) +
  geom_col(fill='dark green', color='white') +
  theme_bw() +
  labs(y = 'Number of links') +
  theme(axis.text = element_text(size=10, color='black'), title = element_text(size=14), axis.text.x = element_text(angle = 90))

## Density (connectance per farm) -------
farm_net_density <- function(d){
  g <- graph_from_data_frame(d, directed = FALSE, vertices = NULL)
  out <- tibble(d=round(igraph::edge_density(g),2))
  return(out)
}

farm_density <- 
  farm_multilayer_pos %>% 
  group_by(level_name) %>% 
  group_modify(~farm_net_density(.x))
  
# Summary table for farms --------
cows_per_farm %>%
  left_join(richness_per_farm) %>%
  left_join(Richness_per_cow_farm %>% select(Farm, ASV_summary)) %>% 
  left_join(links_per_farm, by = c('Farm'='level_name')) %>%
  left_join(farm_density, by = c('Farm'='level_name')) %>%
  as_tibble() %>% 
  rename(Cows=cow_num, ASVs=ASV_Richness, `ASVs per cow`=ASV_summary,Links=links,Density=d) %>% 
  write_csv('local_output/paper_table1.csv')

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


# Clustering coefficient---------------------
farm_multilayer_pos <- read_csv('local_output/farm_multilayer_pos_30.csv')

CC_obs <- 
  farm_multilayer_pos %>%
  group_by(level_name) %>% 
  group_modify(~calc_CC_local(.x)) %>%
  drop_na()


# Distributions of observed CC
pdf('/Users/Geut/Dropbox/for_processing/rumen/cc_farm_boxplot.pdf', 5, 5)
CC_obs %>% 
  filter(k>=10) %>% 
  mutate(level_name=factor(level_name, levels = c("UK1","UK2","IT1","IT2","IT3","FI1",'SE1'))) %>%
  ggplot(aes(x=level_name, y=CC))+
  geom_boxplot()+
  theme_bw() +
  labs(y='Clustering coefficient', x='Farm') +
  paper_figs_theme_no_legend
dev.off()

# CC statistics in each farm
CC_obs %>%
  filter(k>=10) %>% 
  summarise(cc_mean=mean(CC),
            cc_median=median(CC),
            cc_sd=sd(CC),
            cc_min=min(CC),
            cc_max=max(CC)) %>% 
  mutate(ASV_summary=paste(cc_mean,' [',cc_min,'-',cc_max,']',sep=''))

# Partner Fidelity with Jaccard ---------------------------
# Because this is an undirected network, not all ASVs are in the from column,
# which we use for the analysis. So duplicate the links to ensure that all ASVs
# are in the from column.
setdiff(farm_multilayer_pos$from,farm_multilayer_pos$to) # These ASVs were missing
setdiff(farm_multilayer_pos$to,farm_multilayer_pos$from) # These ASVs were missing

farm_multilayer_pos_final <- 
  bind_rows(farm_multilayer_pos, 
            farm_multilayer_pos %>%
              relocate(to, from) %>%
              rename(to=from, from=to))

n_distinct(farm_multilayer_pos$from)
n_distinct(farm_multilayer_pos_final$from)
setdiff(farm_multilayer_pos_final$from,farm_multilayer_pos$to) # These ASVs were missing

# keep only ASVs that occur in 2 or more farms
farm_multilayer_pos_final %<>%
  group_by(from) %>%
  mutate(num_farms_from=n_distinct(level_name)) %>%
  filter(num_farms_from>=2)

## PF_J observed network ---------------------------
PF_J_obs <-
  farm_multilayer_pos_final %>%
  group_by(from) %>%
  group_modify(~calculate_PF_J(.x))

PF_J_obs <- as_tibble(PF_J_obs)

# Partner Fidelity with UniFrac ---------------------------

## PF_U observed network ---------------------------

# read tree from phylo data
# set working directory
phylo_tree <- readRDS("local_output/fitted_asvs_phylo_tree.rds")
tree <- phylo_tree$tree
# tree <- readRDS("local_output/rooted_phylo_tree.rds") # only for 5%

PF_U_obs <-
  farm_multilayer_pos_final %>%
  group_by(from) %>%
  group_modify(~calculate_PF_U(.x, tree))
names(PF_U_obs) <- c("from", "PF_U", "PF_U_sd", "num_layers","UniFrac_type")
PF_U_obs %<>% filter(UniFrac_type=='d_UW')
PF_U_obs$PF_U=1-PF_U_obs$PF_U # Work with similarity instead of dissimilarity
PF_U_obs <- as_tibble(PF_U_obs)

# Plot Jaccard and UniFrac for paper ------------------------------------------
PF_J_obs <- read_csv('local_output/PF_J_pos_30_obs.csv') # change the number according to the filter level
PF_U_obs <- read_csv('local_output/PF_U_pos_30_obs.csv')

pdf('/Users/Geut/Dropbox/for_processing/rumen/cc_pf_hist.pdf', 5, 6)
PF_score_plot <- 
  bind_rows(
    PF_J_obs %>% select(from, PF=PF_J) %>% mutate(type='Jaccard'),
    PF_U_obs %>% select(from, PF=PF_U) %>% mutate(type='UniFrac')
  ) %>% 
  ggplot(aes(PF, fill=type)) +
  geom_histogram(alpha=1, color='white')+
  labs(x='Partner fidelity score', y='Count')+
  paper_figs_theme
PF_score_plot
dev.off()

bind_rows(
  PF_U_obs %>% select(from, PF=PF_U) %>% mutate(type='UniFrac')
) %>% 
  ggplot(aes(PF)) +
  geom_histogram(alpha=1, color='white', fill='lightblue')+
  labs(x='UniFrac score', y='Count')+
  # facet_grid(~type)+
  # geom_vline(xintercept = c(-1.96, 1.96), color = 'red')+
  # geom_vline(xintercept = c(-2.5, 2.5), color = 'red', linetype='dashed')+
  paper_figs_theme


