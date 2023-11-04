#------ filter_data.r -------------------------------------
# This script filter the data to remove scarce microbes.
# 
#------------------------------------------------------------

#------ includes ------------
library(tidyverse)
library(magrittr)
source('functions.R')

#------ functions ---------
filter_data <- function(x) {
  Core_microbes <- 
    ASV_cow_presence %>%
    group_by(Farm) %>% 
    filter(cows >= x*Total_Cows) # Can change the proportion of cows in which a microbe occurs (p)
  ASV_Core <- inner_join(ASV_data_final, Core_microbes, by = c("Farm", "ASV_ID")) %>%
    select(-c(cows, Total_Cows))
  
  # Rerun sensitivity for abundance with the core microbes
  sensitivity_abundance <- NULL
  for (t in seq(0,0.4,0.05)){
    print(t)
    thresh <-  ASV_Core %>% group_by(Farm) %>% summarise(thresh=quantile(Abundance, t))
    ASV_data_filtered_abund <- 
      ASV_Core %>% 
      left_join(thresh) %>% 
      filter(Abundance>thresh)
    
    sensitivity_abundance <- 
      bind_rows(sensitivity_abundance,
                inner_join(ASV_Core %>% group_by(Farm) %>% summarise(N=length(unique(ASV_ID))),
                           ASV_data_filtered_abund %>% group_by(Farm) %>% summarise(N_filtered=length(unique(ASV_ID)))
                ) %>% 
                  mutate(prop_left=N_filtered/N,
                         t=t)
      )
  }
  ggplot(sensitivity_abundance, aes(t, prop_left, color=Farm))+
    geom_line()+geom_point()+theme_bw()+
    facet_wrap(~Farm)+scale_y_continuous(limits=c(0,1))
  
  # Check the abundance distributions of core microbes
  ASV_Core %>% 
    ggplot(aes(Abundance))+geom_histogram()+facet_wrap(~Farm, scales = 'free')+theme_bw()
  
  ASV_Core %>% group_by(Farm) %>% summarise(N=n_distinct(ASV_ID))
  
  # Write the final data set of filtered microbes
  sx <- sprintf("%02d", x*100)
  write_csv(ASV_Core, paste("local_output/core_ASV_", sx, ".csv", sep = ""))
  
}


# Run --------------------------------
ASV_data_final <- read_csv("local_output/ASV_processed_data.csv")
ASV_data_final %>% group_by(Country,Farm) %>% summarise(cows=n_distinct(Cow_Code))
ASV_data_final %>% group_by(Country,Farm) %>% summarise(ASVs=n_distinct(ASV_ID))

# Match farm names to those in the original paper
ASV_data_final %<>% 
  mutate(Farm=replace(Farm, Farm=='NUDC', 'UK1')) %>% 
  mutate(Farm=replace(Farm, Farm=='Park', 'UK2')) %>% 
  mutate(Farm=replace(Farm, Farm=='Bianchini', 'IT1')) %>% 
  mutate(Farm=replace(Farm, Farm=='Franciosi', 'IT2')) %>% 
  mutate(Farm=replace(Farm, Farm=='Gandolfi', 'IT3')) %>%
  mutate(Farm=replace(Farm, Farm=='MinkiÃ¶', 'FI1')) %>% 
  mutate(Farm=replace(Farm, Farm=='RÃ¶bÃ¤cksdalen', 'SE1'))

ASV_data_final %>% distinct(Country, Farm)

# Total number of ASVs
length(unique(ASV_data_final$ASV_ID))

# Number of ASVs in each farm
ASV_data_final %>%
  group_by(Country, Farm) %>% summarise(n=n_distinct(ASV_ID))


## ASV abundance in farm ---------------------------------------------------

# ASV abundance - filter out microbes that are very rare within the farm
ASV_data_final %>% arrange(Farm, Abundance)

sensitivity_abundance <- NULL
for (t in seq(0,0.4,0.05)){
  print(t)
  thresh <- ASV_data_final %>%
    group_by(Farm) %>% 
    summarise(thresh=quantile(Abundance, t))
  
  ASV_data_filtered_abund <- 
    ASV_data_final %>% 
    left_join(thresh) %>% 
    filter(Abundance>thresh)
  
  sensitivity_abundance <- bind_rows(sensitivity_abundance,
                                     inner_join(
                                       ASV_data_final %>% group_by(Farm) %>% summarise(N=length(unique(ASV_ID))),
                                       ASV_data_filtered_abund %>% group_by(Farm) %>%
                                         summarise(N_filtered=length(unique(ASV_ID)))) %>% 
                                       mutate(prop_left=N_filtered/N,t=t))
}

ggplot(sensitivity_abundance, aes(t, prop_left, color=Farm))+
  geom_line()+geom_point()+theme_bw()+facet_wrap(~Farm)


## ASV relative abundance in cows ------------------------------------------
# Filter out microbes whose relative abundance within the sample (cow) is low

ASV_data_final %>% arrange(Farm, Abundance)

tot_abund <- 
  ASV_data_final %>% group_by(Country,Farm,Cow_Code) %>% 
  summarise(tot=sum(Abundance))
rel_abund <- 
left_join(ASV_data_final,tot_abund) %>% 
  mutate(rel_abund=Abundance/tot)

rel_abund %>% 
  group_by(Farm) %>% summarise(max(rel_abund))

sensitivity_rel_abundance <- NULL
for (t in seq(0,0.05,0.001)){
  print(t)
  ASV_data_filtered_abund <- 
    rel_abund %>% 
    group_by(Cow_Code) %>% 
    filter(rel_abund>t)
  
  sensitivity_rel_abundance <- bind_rows(sensitivity_rel_abundance,
                                     inner_join(
                                       ASV_data_final %>% group_by(Farm) %>% summarise(N=length(unique(ASV_ID))),
                                       ASV_data_filtered_abund %>% group_by(Farm) %>%
                                         summarise(N_filtered=length(unique(ASV_ID)))) %>% 
                                       mutate(prop_left=N_filtered/N,t=t))
}
ggplot(sensitivity_rel_abundance, aes(t, prop_left, color=Farm))+
  geom_line()+geom_point()+theme_bw()+facet_wrap(~Farm)


## Core microbes -----------------------------------------------------------

# Core microbes are defined as those occurring in a certain proportion of cows
# within each farm
ASV_cow_presence <- ASV_data_final %>% 
  group_by(Farm,ASV_ID) %>%
  summarise(cows=n_distinct(Cow_Code)) %>%
  arrange(desc(cows))

# how many cows in each farm:
cows_in_farms <- 
  ASV_data_final %>%
  group_by(Farm) %>% 
  summarise(Total_Cows=n_distinct(Cow_Code))
sum(cows_in_farms$Total_Cows) # Total number of cows in the region

ASV_cow_presence %<>% left_join(cows_in_farms)
# An example for a single farm
ASV_cow_presence %>%
  filter(Farm=='UK1') %>% 
  group_by(cows) %>% 
  count() %>% 
  ggplot(aes(x=cows, y=n))+geom_col()

sensitivity_core <- NULL
for (p in seq(0,0.6,0.05)){
  print(p)
  Core_microbes <- 
    ASV_cow_presence %>%
    group_by(Farm) %>% 
    filter(cows >= p*Total_Cows) # Can change the proportion of cows in which a microbe occurs (p)
  sensitivity_core <- bind_rows(sensitivity_core,
                                inner_join(
                                  ASV_data_final %>% group_by(Farm) %>% summarise(N=n_distinct(ASV_ID)),
                                  Core_microbes %>% group_by(Farm) %>% summarise(N_filtered=n_distinct(ASV_ID))
                                ) %>% 
                                  mutate(prop_left=N_filtered/N,p=p)
  )
}

# Note that from the onset, even with p=0, not all microbes occur in all cows, so the starting value for p=0 is not 1.
pdf(paste(paper_output_path, "sensitivity_core.pdf", sep=""), 10, 6)
ggplot(sensitivity_core, aes(p, prop_left, color=Farm))+
  geom_line()+
  geom_point()+
  facet_wrap(~Farm)+
  geom_vline(xintercept = c(0.3), linetype='dashed')+
  scale_x_continuous(breaks = seq(0, 0.6, 0.15))+
  labs(x='% of cows in which microbes occur', y='Proportion of microbes defined as core')+
  paper_figs_theme_no_legend
dev.off()

pdf('local_output/figures/core_microbes_sensitivity_example.pdf', 10, 6)
sensitivity_core %>% 
  filter(Farm=='UK1') %>% 
  ggplot(aes(p, prop_left, color=Farm))+
  geom_line(size=1.5)+geom_point()+
  # facet_wrap(~Farm)+
  geom_vline(xintercept = c(0.05,0.3,0.5), linetype='dashed')+
  labs(x='Proportion of cows in which microbes occur', y='Proportion of microbes defined as core')+
  theme_bw()+
  theme(axis.text = element_text(size=22),
        axis.title = element_text(size=22),
        legend.position = 'none')
dev.off()


# Compare core micrbes with relative abundance filtering ------------------
p=0.3
Core_microbes <- ASV_cow_presence %>%
  group_by(Farm) %>% 
  filter(cows >= p*Total_Cows) 

t=0.005
ASV_data_filtered_abund <- 
  rel_abund %>% 
  group_by(Cow_Code) %>% 
  filter(rel_abund>t)

n <- intersect(ASV_data_filtered_abund$ASV_ID,Core_microbes$ASV_ID)
d <- union(ASV_data_filtered_abund$ASV_ID,Core_microbes$ASV_ID)
length(n)/length(d)


inner_join(ASV_cow_presence, rel_abund) %>% 
  mutate(prop_cows=cows/Total_Cows) %>% 
  distinct(ASV_ID,prop_cows,rel_abund) %>% 
  filter(prop_cows>0.3) %>% 
  ggplot(aes(prop_cows,rel_abund))+
  geom_point()
  
# Comparing the two filtering methods, the core microbes is better because it
# removes also the scarce microbes, which are not present in many cows by
# definition. So we will subset the data by core microbes. The huge drop in the
# plot is at p=0.05




# Output data -------------------------------------------------------------

# filtering the data by percents
# core microbes in 50% of cows
filter_data(0.5)
# core microbes in 5% of cows
filter_data(0.05)
# core microbes in 30% of cows
filter_data(0.3)
# core microbes in 20% of cows
filter_data(0.2)
# core microbes in 10% of cows
filter_data(0.1)

