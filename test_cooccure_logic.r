# here i will try to calculate P(lt) and P(gt) crom the cooccur index
# using: choose(n, r)

# functions:
calc_pj <- function(N, N1, N2, j) {
  n1_j <- choose(N1,j)
  nn1_n2j <- choose(N-N1,N2-j)
  n_n2 <- choose(N,N2)
  
  pj <- n1_j * nn1_n2j / n_n2
  
  return(pj)
}

calc_plt <- function(N, N1, N2, Qobs){
  sum_comb <- 0
  for (i in 0:(Qobs-1)) {
    #print(i)
    sum_comb <- sum_comb + calc_pj(N, N1, N2, i)
  }
  return(sum_comb)
}

calc_pgt <- function(N, N1, N2, Qobs){
  sum_comb <- 0
  for (i in (Qobs+1):N) {
    #print(i)
    sum_comb <- sum_comb + calc_pj(N, N1, N2, i)
  }
  return(sum_comb)
}

# calculate p-value
# test calculations
N=100 # cows 
N1=100
N2=20
Qobs=20

Pet <- calc_pj(N=N, N1=N1, N2=N2, j=Qobs)
Plt <- calc_plt(N=N, N1=N1, N2=N2, Qobs=Qobs) 
Pgt <- calc_pgt(N=N, N1=N1, N2=N2, Qobs=Qobs)

Plt + Pet # in package: Plt
Pgt + Pet # in package: Pgt
Pet
Plt+Pgt+Pet # == 1


# test calculating size effect
expct <- (N1/N)*(N2/N)*N
size_effect <- (Qobs-expct)/N





# test calculating size effect behavior ----
# case 1
N=17
N1=14
N2=14
Qobs<- 1:14

expct <- (N1/N)*(N2/N)*N
size_effect <- (Qobs-expct)/N
plot(Qobs, size_effect)

# case 2
N=4:15
N1=4
N2=4
Qobs=4
expct <- (N1/N)*(N2/N)*N
size_effect <- (Qobs-expct)/N

plot(N, size_effect)

# case 3
N=8
a=1:N # == N1 == N2 == Qobs

expct <- (a/N)*(a/N)*N
size_effect <- (a-expct)/N

plot(a, size_effect)


# test jaccard values?


