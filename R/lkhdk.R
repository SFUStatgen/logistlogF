lkhdk <- function(m,data_k,params_k,N) {
  # Input:
  # - m is the value of m
  # - data_k is the data consisting of the phenotype, SNV_k and confounders
  # - params_k is the output of MCEM()
  # - N is the number of Monte Carlo replicates
  # Output: Monte Carlo estimate of the profile likelihood
  
  betas <- log(rf(N,m,m))
  lvec <- rep(NA,N)
  for(j in 1:N) {
    s <- params_k[1]+data_k$X*betas[j]+as.matrix(data_k[,-c(1,2)])%*%(params_k[-1])
    lvec[j] <- prod(exp(data_k$Phenotype*s)/(1+exp(s)))
  }
  return(log(mean(lvec)))
}