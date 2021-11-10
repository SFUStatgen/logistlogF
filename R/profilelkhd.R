profilelkhd <- function(m,data,weight,Y_index,X_index,C_index,N) {
  # Input: 
  # - m isthe value of m
  # - data is the GWAS data
  # - weight is a list of weights for SNVs
  # - Y_index is the index of phenotype in data
  # - X_index is the index of SNVs in data
  # - C_index is the index of confounders in data
  # - N is the number of Monte Carlo replicates
  # Output: profile marginal likelihood of m
  
  ll <- 0
  K <- length(X_index)
  for(k in 1:K) {
    cat("k =",k,"\n")
    if (length(C_index) == 0) {data_k <- data[,c(Y_index,X_index[k])]}
    else {data_k <- data[,c(Y_index,X_index[k],C_index)]}
    colnames(data_k)[1:2]=c("Phenotype","X")
    params_k <- MCEM(m,data_k,C_index,N) 
    ll <- ll+weight[k]*lkhdk(m,data_k,params_k,N)
  }
  return(ll)
}