
library(bipartite)
library(tidyverse)

# This function removes a single species from either the rows or the columns
single_extinct <- function(m, margin=NULL, x){
  if (margin==1){m[x,] <- 0}
  if (margin==2){m[,x] <- 0}
  return(m)
}

# This function counts the number of remaining species in either rows or columns
record_2nd_extinctions <- function(m, margin=NULL){
  # Notice that the margins are the opposite!!!
  # Because if removal is from rows we need to quantify columns and vice versa
  if (margin==2){x <- sum(rowSums(m)==0)}
  if (margin==1){x <- sum(colSums(m)==0)}
  return(x)
}

A <- readxl::read_excel('extinction_game.xlsx', col_names = paste('a',1:10,sep='_'), sheet = 1, range = 'A1:J10')
num_extinct <- 0
extinction_sequence <- order(colSums(A))
# Loop through the extinction sequence
for (e in extinction_sequence){
  # print(e)
  A <- single_extinct(A,2,e)
  num_extinct <- c(num_extinct, record_2nd_extinctions(A, 2))
}


# This is the main function
extinct <- function(m, margin, extinction_sequence){
  # m is the matrix, margin is from where to remove (1 rows 2 columns)
  
  num_extinct <- 0 # Initialize the number of species going extinct
  # Loop through the extinction sequence
  for (e in extinction_sequence){
    # print(e)
    m <- single_extinct(m,margin,e)
    num_extinct <- c(num_extinct, record_2nd_extinctions(m, margin))
  }
  # Produce the results
  if (margin==1){
    results <- data.frame(num_removed=0:nrow(m), num_extinct=num_extinct)
    results$prop_removed <- results$num_removed/nrow(m)
    results$prop_remain <- 1-results$num_extinct/ncol(m)
  }
  if (margin==2){
    results <- data.frame(num_removed=0:ncol(m), num_extinct=num_extinct)
    results$prop_removed <- results$num_removed/ncol(m)
    results$prop_remain <- 1-results$num_extinct/nrow(m)
  }
  return(results)
}

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
Farm_IT1_net <- acast(Farm_IT1_edge_list_final, from~to, value.var = 'links', fill = 0)
num_extinct <- 0
extinction_sequence <- order(colSums(Farm_IT1_net),decreasing = TRUE)
# Loop through the extinction sequence
for (e in extinction_sequence){
  # print(e)
  if (sum(Farm_IT1_net[,e])==0) {
    next
  }
  prev_extinctions <- sum(rowSums(Farm_IT1_net)==0)
  Farm_IT1_net <- single_extinct_unipartite(Farm_IT1_net,e)
  num_extinct <- c(num_extinct, record_2nd_extinctions_unipartite(Farm_IT1_net)-prev_extinctions)
  # if (e==9) {
  #   break
  # }
}

x <- sum(rowSums(A)==0) 

extinction_sequence[292]
