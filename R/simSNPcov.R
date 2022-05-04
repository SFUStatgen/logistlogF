simSNPCov_popl <- function(n,Beta,MAF1,MAF2,beta_popl=1,scale=FALSE){
  # Input:
  # - n is the number of observations
  # - Beta is the effect sizes of SNVs
  # - MAF1 and MAF2 are minor allele frequencies in two populations
  # - beta_popl is the coefficient associated with population 
  # - scale indicates whether the covariate should be standardized
  
  ncase <- ncon <- n/2
  conX <- matrix(NA,ncol=K+1,nrow=ncon)
  caseX <- matrix(NA,ncol=K+1,nrow=ncase)
  
  # Simulate population status first (population is the K+1 column)
  f_0 <- f_1 <- 0.5
  conX[,K+1] <- rbinom(ncon,1,f_0)
  caseX[,K+1] <- sample(0:1,size=ncase,replace=TRUE,prob=c(f_0,f_1*exp(beta_popl)))
  
  # Simulate SNVs conditional on population status
  for (i in 1:K){ 
    maf0 <- MAF1[i]
    maf1 <- MAF2[i]
    beta <- Beta[i]
    
    ## control
    con_pop0 <- which(conX[,K+1]==0)
    con_pop1 <- which(conX[,K+1]==1)
    
    conX[con_pop0,i] <- rbinom(length(con_pop0),size=2,prob=maf0)
    conX[con_pop1,i] <- rbinom(length(con_pop1),size=2,prob=maf1)
    
    ## case
    case_pop0 <- which(caseX[,K+1]==0)
    case_pop1 <- which(caseX[,K+1]==1)
    
    pp0 <- c((1-maf0)^2,2*maf0*(1-maf0)*exp(beta),maf0^2*exp(2*beta))
    pp1 <- c((1-maf1)^2,2*maf1*(1-maf1)*exp(beta),maf1^2*exp(2*beta))
    caseX[case_pop0,i] <- sample(0:2,size=length(case_pop0),replace=TRUE,prob=pp0)
    caseX[case_pop1,i] <- sample(0:2,size=length(case_pop1),replace=TRUE,prob=pp1)
  }
  X <- rbind(caseX,conX)
  if(scale) X <- scale(X)
  colnames(X) <- c(paste0("X",1:K),"Population");rownames(X) <- NULL
  case <- c(rep(1,ncase),rep(0,ncon))
  return (data.frame(X,case))
}

