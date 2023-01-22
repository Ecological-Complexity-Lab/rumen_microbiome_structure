library(dplyr)
library(magrittr)
library(bipartite)
library(tidyverse)
library(reshape2)
library(ggplot2)
setwd("~/GitHub/microbiome_structure_v2")
Farm_IT1_edge_list <- read_csv('fixed_data/HPC/exp_12/12_3945654_Farm_IT1_edge_list.csv')

setdiff(Farm_IT1_edge_list$from,Farm_IT1_edge_list$to) # These ASVs were missing
setdiff(Farm_IT1_edge_list$to,Farm_IT1_edge_list$from) # These ASVs were missing
Farm_IT1_edge_list_final <- 
  bind_rows(Farm_IT1_edge_list, 
            Farm_IT1_edge_list %>%
              relocate(to, from) %>%
              rename(to=from, from=to))
n_distinct(Farm_IT1_edge_list$from)
n_distinct(Farm_IT1_edge_list$to)
n_distinct(Farm_IT1_edge_list_final$from)
n_distinct(Farm_IT1_edge_list_final$from)
setdiff(Farm_IT1_edge_list_final$from,Farm_IT1_edge_list_final$to) 

Farm_IT1_edge_list_final %<>% mutate(links=1) 
Farm_IT1_net <- acast(Farm_IT1_edge_list_final, from~to, value.var = 'links', fill = 0)

# demo----
small_net <- Farm_IT1_net[1:10,1:10]

# This function removes a single species from both rows and columns
single_extinct_unipartite <- function(m, x){
  {m[x,] <- 0}
  {m[,x] <- 0}
  return(m)
}

# This function counts the number of remaining species in either rows or columns
record_2nd_extinctions_unipartite <- function(m){
  # Because if removal is from rows we need to quantify columns and vice versa
  {x <- sum(rowSums(m)==0)-1}
  #{x <- sum(colSums(m)==0)}
  return(x)
}

small_net <- Farm_IT1_net[1:10,1:10]
num_extinct <- NULL
extinction_sequence <- order(colSums(small_net),decreasing = TRUE)
# Loop through the extinction sequence
for (e in extinction_sequence){
  # print(e)
  if (sum(small_net[,e])==0) {
    num_extinct <- c(num_extinct, 0)
    next
  }
  prev_extinctions <- sum(rowSums(small_net)==0)
  small_net <- single_extinct_unipartite(small_net,e)
  num_extinct <- c(num_extinct, record_2nd_extinctions_unipartite(small_net)+1-prev_extinctions)
  # if (e==2) {
  #   break
  # }
}
# Produce the results
results <- data.frame(num_removed=1:nrow(small_net), num_extinct=num_extinct)
results$prop_removed <- results$num_removed/nrow(small_net)
results %<>%
  mutate(prop_remain = 1-(cumsum(num_extinct)/nrow(small_net)))

# This is the main function
small_net <- Farm_IT1_net[1:10,1:10]
extinction_sequence <- order(colSums(small_net),decreasing = TRUE)
extinct_unipartite <- function(m, extinction_sequence){
  # m is the matrix, margin is from where to remove (1 rows 2 columns)
  
  num_extinct <- NULL # Initialize the number of species going extinct
  # Loop through the extinction sequence
  for (e in extinction_sequence){
    # print(e)
    if (sum(m[,e])==0) {
      num_extinct <- c(num_extinct, 0)
      next
    }
    prev_extinctions <- sum(rowSums(m)==0)
    m <- single_extinct_unipartite(m,e)
    num_extinct <- c(num_extinct, record_2nd_extinctions_unipartite(m)+1-prev_extinctions)
  }
  return(num_extinct)
  # Produce the results
  results <- data.frame(num_removed=0:nrow(m), num_extinct=num_extinct)
  results$prop_removed <- results$num_removed/nrow(m)
  results %<>%
    mutate(prop_remain = 1-(cumsum(num_extinct)/nrow(m)))
  print(results)
}
extinct_unipartite(small_net,extinction_sequence)

play_game_demo <- function(){
  # Read file with matrix
  A <- Farm_IT1_net[1:10,1:10]
  print(A) # Show the matrix
  df <- extinct_unipartite(A, extinction_sequence = order(colSums(A),decreasing = TRUE)) # Run extinction function
  
  # Calculate area under the curve
  y <- df[,2]
  y <- (sum(y) - cumsum(y))/sum(y)
  x <- df$prop_removed
  # x <- (object[, "no"]/max(object[, "no"]))
  ext.curve <- splinefun(x, y)
  ext.area <- integrate(ext.curve, 0, 1)
  R <- as.numeric(ext.area[[1]])
  # Plot
  p <- ggplot(df, aes(prop_removed, prop_remain))+
    geom_point(size=4, color='red')+
    geom_line(size=1, color='red')+
    labs(title=paste('Your score: ',round(100*R,1)), x='Proportion removed', y='Proportion remained') +
    theme_bw() # Make a plot
  print(p)
}






# one farm----

# This function removes a single species from both rows and columns
single_extinct_unipartite <- function(m, x){
  {m[x,] <- 0}
  {m[,x] <- 0}
  return(m)
}

# This function counts the number of remaining species in either rows or columns
record_2nd_extinctions_unipartite <- function(m){
  # Because if removal is from rows we need to quantify columns and vice versa
  {x <- sum(rowSums(m)==0)-1}
  #{x <- sum(colSums(m)==0)}
  return(x)
}

Farm_IT1_net <- acast(Farm_IT1_edge_list_final, from~to, value.var = 'links', fill = 0)
num_extinct <- NULL
extinction_sequence <- order(colSums(Farm_IT1_net),decreasing = TRUE)
# Loop through the extinction sequence
for (e in extinction_sequence){
  # print(e)
  if (sum(Farm_IT1_net[,e])==0) {
    num_extinct <- c(num_extinct, 0)
    next
  }
  prev_extinctions <- sum(rowSums(Farm_IT1_net)==0)
  Farm_IT1_net <- single_extinct_unipartite(Farm_IT1_net,e)
  num_extinct <- c(num_extinct, record_2nd_extinctions_unipartite(Farm_IT1_net)+1-prev_extinctions)
  # if (e==2) {
  #   break
  # }
}
# Produce the results
results <- data.frame(num_removed=1:nrow(Farm_IT1_net), num_extinct=num_extinct)
results$prop_removed <- results$num_removed/nrow(Farm_IT1_net)
results %<>%
  mutate(prop_remain = 1-(cumsum(num_extinct)/nrow(Farm_IT1_net)))

# The main function
Farm_IT1_net <- acast(Farm_IT1_edge_list_final, from~to, value.var = 'links', fill = 0)
extinction_sequence <- order(colSums(Farm_IT1_net),decreasing = TRUE)
extinct_unipartite <- function(m, extinction_sequence){
  # m is the matrix, margin is from where to remove (1 rows 2 columns)
  
  num_extinct <- NULL # Initialize the number of species going extinct
  # Loop through the extinction sequence
  for (e in extinction_sequence){
    # print(e)
    if (sum(m[,e])==0) {
      num_extinct <- c(num_extinct, 0)
      next
    }
    prev_extinctions <- sum(rowSums(m)==0)
    m <- single_extinct_unipartite(m,e)
    num_extinct <- c(num_extinct, record_2nd_extinctions_unipartite(m)+1-prev_extinctions)
  }
  # Produce the results
  results <- data.frame(num_removed=1:nrow(m), num_extinct=num_extinct)
  results$prop_removed <- results$num_removed/nrow(m)
  results %<>%
    mutate(prop_remain = 1-(cumsum(num_extinct)/nrow(m)))
  return(results)
}
extinct_unipartite(Farm_IT1_net,extinction_sequence)

play_game <- function(){
  # Read file with matrix
  A <- acast(Farm_IT1_edge_list_final, from~to, value.var = 'links', fill = 0)
  print(A) # Show the matrix
  df <- extinct_unipartite(A, extinction_sequence = order(colSums(A),decreasing = TRUE)) # Run extinction function
  
  # Calculate area under the curve
  y <- df[,2]
  y <- (sum(y) - cumsum(y))/sum(y)
  x <- df$prop_removed
  # x <- (object[, "no"]/max(object[, "no"]))
  ext.curve <- splinefun(x, y)
  ext.area <- integrate(ext.curve, 0, 1)
  R <- as.numeric(ext.area[[1]])
  # Plot
  p <- ggplot(df, aes(prop_removed, prop_remain))+
    geom_point(size=4, color='red')+
    geom_line(size=1, color='red')+
    labs(title=paste('Your score: ',round(100*R,1)), x='Proportion removed', y='Proportion remained') +
    theme_bw() # Make a plot
  print(p)
}

stability_output <- play_game()

stability_correlation_p <- stability_output$plot_env$p 
png(filename = 'local_output/figures/stability_IT1.png', width = 1300, height = 900, res = 300)
stability_correlation_p + 
  ggtitle("Stability of farm IT1 network") +
  theme_bw()+
  theme(panel.grid=element_blank(),
        axis.text = element_text(size=12, color='black'),
        axis.title = element_text(size=12, color='black'))
dev.off()
  
