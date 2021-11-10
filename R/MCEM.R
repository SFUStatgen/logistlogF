MCEM <- function(m,data_k,C_index,N) {
  # Input:
  # - m is the value of m
  # - data_k is the data consisting of the phenotype, SNV_k and confounders
  # - C_index is the index of confounders in data
  # - N is the number of Monte Carlo replicates
  # Output: params_k
  
  model <- glm(Phenotype~.,data=data_k,family=binomial(link="logit"),maxit=100)
  initial_params <- model$coefficients[-2] # remove coefficient of X
  c <- length(C_index) # number of confounders
  params <- matrix(0,ncol=c+1)
  params <- rbind(params,initial_params)
  
  p <- 2
  threshold <- 1E-04
  
  Weight <- function(beta) {
    s <- params[p,1]+data_k$X*beta+as.matrix(data_k[,-c(1,2)])%*%(params[p,-1])
    prod(exp(data_k$Phenotype*s)/(1+exp(s)))
  }
  
  Y <- rep(data_k$Phenotype,times=N)
  Cov <- matrix(data=NA,ncol=c, nrow=N*dim(data_k)[1])
  if (c != 0) {
    for (i in 1:c) {
      Cov[,i] <- rep(data_k[,2+i])
    }
  }
  betas <- log(rf(N,m,m))
  
  O <- numeric() # offset
  for (j in 1:N) {
    O <- c(O,data_k$X*betas[j])
  }
  
  while(norm(as.matrix(params[p,]-params[p-1,]))>=threshold) {
    W_t <- numeric() # weight
    for (j in 1:N) {
      W_t[j] <- Weight(betas[j])
    }
    W <- rep(W_t,each=dim(data_k)[1])
    
    if (c != 0) {g <- glm(Y~offset(O)+Cov,weights=W,family=binomial(link="logit"),maxit=100)}
    else {g <- glm(Y~offset(O),weights=W,family=binomial(link="logit"),maxit=100)}
    #cat("EM iteration",p-1,":",g$coefficients,"\n")
    p <- p+1
    params <- rbind(params,g$coefficients)
  }
  return(params[p,])
}