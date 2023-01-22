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

# If using google sheets
# drive_auth(email = 'shainova@gmail.com')
# drive_find(type = "spreadsheet") 
# (x <- drive_get('extinction_game'))
# (x <- drive_get(as_id('148NgXT3sdrNulcS81qkN12t3F_EuGEpx-OefMusXgiw')))

show_network <- function(){
  A <- readxl::read_excel('extinction_game.xlsx', col_names = paste('a',1:10,sep='_'), sheet = 1, range = 'A1:J10')
  # Convert to a matrix
  A <- data.matrix(A)
  A[is.na(A)] <- 0 # remove NA values
  plotweb(A)
}


play_game <- function(){
  # Read file with matrix
  # A <- read_sheet(ss = '148NgXT3sdrNulcS81qkN12t3F_EuGEpx-OefMusXgiw', col_names = paste('a',1:10,sep='_'), sheet = 1, range = 'A1:J10')
  A <- readxl::read_excel('extinction_game.xlsx', col_names = paste('a',1:10,sep='_'), sheet = 1, range = 'A1:J10')
  # Convert to a matrix
  A <- data.matrix(A)
  A[is.na(A)] <- 0 # remove NA values
  print(A) # Show the matrix
  df <- extinct(A, 2, extinction_sequence = order(colSums(A))) # Run extinction function
  
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
