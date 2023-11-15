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

# ------ consts ----------
data_file <- "local_output/core_ASV_30.csv"
nsim <- 500

# ------ functions ----------
suffle_cow_microb_vegan <- function(data_file_name, nsim=500, output_folder="shuff_vegan_30", shuff_method='curveball'){
  asv_data <- read_csv(data_file_name)
  
  # get layers names
  farm_names <- asv_data %>% select(Farm) %>% distinct() %>% pull(Farm)
  
  layers <- list()
  # make shuffled matrices -
  # have to separate to farms in order to not shuffle between them
  for (frm in farm_names) {
    lyr <- asv_data %>% filter(Farm==frm)
    
    # convert edgelist to matrix
    lyr_mat <- lyr %>%
      select(Cow_Code, ASV_ID) %>% add_column(weight=1) %>%
      dcast(Cow_Code ~ ASV_ID, value.var = "weight", fill = 0) %>%
      column_to_rownames(var="Cow_Code")
    
    # run curveball shuffling
    null <- vegan::nullmodel(lyr_mat, method = shuff_method)
    suff <- simulate(null, nsim = nsim, burnin = 5000, seed = 1234)
    
    layers[[frm]] <- suff
  }
  
  # make files from shuffled
  for (i in 1:nsim) {
    all_farms <- tibble(Country = character(),
                        Farm = character(),
                        Cow_Code = character(),
                        ASV_ID = character(),
                        Abundance = numeric(),
                        shuff_id = numeric())
    # convert multi-mats to "edge lists"
    for (frm in farm_names) {
      matt <- layers[[frm]][,,i]
      # convert
      edge <- melt(matt) %>% filter(value == 1) %>%
        select(Cow_Code=Var1, ASV_ID=Var2, Abundance=value) %>%
        add_column(shuff_id=i) %>%
        add_column(Farm=frm, .before = 1) %>%
        add_column(Country=substr(frm, 1, 2), .before = 1)
      all_farms <- rbind(all_farms, edge)
    }
    # save to file
    write_csv(all_farms, paste(output_folder,'/shuff_farm_',str_pad(i, 3, '0',side='left'),'.csv', sep=''))
  }
  
  write_csv(tibble(e_id=1:nsim, 
                   data_file=paste('shuff_farm_',str_pad(1:nsim, 3, '0',side='left'),'.csv', sep=''),
                   Abundance_file=basename(data_file_name)),
            paste(output_folder,'/experiments.csv', sep=''))
}


# ------ run --------------------------------
# What if bacteria were distributed at random?

# save distribution shuffling ---------------------------------------------------------------

# This permutes the microbes between cows within farms. BUT it fixes the degree
# of the cows. So no change in ASV number per cow.
output <- "local_output/shuffle_farm_r0_30_500"
suffle_cow_microb_vegan(data_file, nsim, output, shuff_method="r0")


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
