#------ shuffle_microbe_distribution.r -------------------------------------
# This script performs a shuffling for the networks to comaper the 
# resultd and in that way to validate them. shuffling is done in the farm,
# region, country and global level
#------------------------------------------------------------

#------ includes ------------
library(tidyverse)
library(magrittr)
library(vegan)
library(reshape2)
source('functions.R')

#------ consts ----------
data_file <- "local_output/core_ASV_30.csv"
nsim <- 500

#------ run --------------------------------
# What if bacteria were distributed at random?

# Correlation between the number of times a microbe appears in the data set
# (=# cows in which a microbe is present) and the # of farms
ASV_data <- read_csv(data_file)

inner_join(
  ASV_data %>% group_by(ASV_ID) %>% summarise(n_cows=n_distinct(Cow_Code)),
  ASV_data %>% group_by(ASV_ID) %>% summarise(n_farms=n_distinct(Farm))
) %>% ggplot(aes(n_cows, n_farms))+geom_point()


inner_join(
  ASV_data %>% group_by(ASV_ID) %>% summarise(n_cows=n_distinct(Cow_Code)),
  ASV_data %>% group_by(ASV_ID) %>% summarise(n_countries=n_distinct(Country))
) %>% ggplot(aes(n_cows, n_countries))+geom_point()

# Expected (random) distribution of microbes in farms when microbes are randomly
# distributed across cows in all the region

permuted <- NULL
for (i in 1:500){
  x <- ASV_data %>% 
    mutate(ASV_ID_perm=sample(ASV_ID, replace = F)) %>% # Permute
    group_by(ASV_ID_perm) %>%
    summarise(habitats=n_distinct(Farm)) %>%
    arrange(desc(habitats)) %>% 
    group_by(habitats) %>% 
    summarise(n=n()) %>% mutate(run=i)
  permuted <- rbind(permuted, x)
}

permuted %>% group_by(habitats) %>% summarise(n_mean=mean(n), sd=sd(n)) %>% 
  ggplot(aes(habitats, n_mean))+
  geom_col(fill='brown', color='black')+
  geom_errorbar(aes(xmin=habitats,xmax=habitats, ymin=n_mean-sd, ymax=n_mean+sd), width=0.3)+
  labs(title='Random distribution of microbes in farms')+
  scale_x_continuous(breaks = seq(0,7,1))

# Beta-div between farms
ASV_occurrence_farm <- 
  ASV_data %>%
  mutate(ASV_ID_perm=sample(ASV_ID, replace = F)) %>% # Permute
  select(-c(Cow_Code)) %>%
  distinct(Farm, ASV_ID_perm) %>%
  mutate(present=1) %>% 
  spread(ASV_ID_perm, present, fill = 0) %>%
  column_to_rownames("Farm")
dim(ASV_occurrence_farm)
vegdist(ASV_occurrence_farm, "jaccard")

# save distribution shuffling ---------------------------------------------------------------

# This function permutes the microbes between cows WITHIN a farm. Microbes are not
# shuffled between farms.
output <- "local_output/shuffle_farm_30"
shuffle_farm_microbe(data_file, nsim, output)

# This permutes the microbes between cows between farms. Microbes are 
# shuffled between farms.
output <- "local_output/shuffle_all_30"
shuffle_region_microbe(data_file, nsim, output)

# This permutes the microbes between cows between farms. BUT it fixes the degree
# of the cows as well as the degree of the Microbes. 
# So no change in species number per cow. and no changes in cows per species.
output <- "local_output/shuffle_farm_curveball_30"
suffle_cow_microb_vegan(data_file, nsim, output)

# This permutes the microbes between cows within farms. BUT it fixes only the degree
# of the cows. So no change in species number per cow.
output <- "local_output/shuffle_farm_r0_30_500"
suffle_cow_microb_vegan(data_file, nsim, output, shuff_method="r0")

# This permutes the microbes between cows between farms. it fixes only the fill
# (the number of occurrences in the matrix.
output <- "local_output/shuffle_farm_r00_30"
suffle_cow_microb_vegan(data_file, nsim, output, shuff_method="r00")

# Vaildate shuffling ------------------------------------------------------
ASV_data <- read_csv(data_file)
ASVs_obs <- unique(ASV_data$ASV_ID)
files <- list.files(path = "local_output/shuffle_farm_r0_30_500/", pattern = "shuff_farm", full.names = T)
for (f in files){
  f <- files[1]
  x <- suppressMessages(read_csv(f))
  ASVs_shuff <- unique(x$ASV_ID)
  print(setequal(ASVs_obs, ASVs_shuff))
}
