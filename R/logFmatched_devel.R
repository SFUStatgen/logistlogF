# log-F(m,m)-penalized conditional likelihood inference by data augmentation
augment.logFmatched = function(form,data,m) {
  # Input:
  # - form is an R formula
  # - dat is the data
  # - m is true value of m
  # Output: 
  # - augmented dataset
  
  # Step 1: Extract (i) the response and (ii) the design matrix 
  # from the input formula and data frame so that we can augment them. 
  mf = model.frame(form,data)
  D = model.response(mf)     # extract the response
  X = model.matrix(form,data) # extract the design matrix
  if(ncol(X)==1) { # intercept only model, no augmentation needed
    return(X)
  } else {
    X = model.matrix(form,data)[,-1,drop=FALSE] # we don't want the intercept
  }
  
  # Step 2 (augmentation): For an even degree of freedom m, add
  # m pseudo-matched sets of size 2 for each covariate:
  # In the first m/2 matched set, the case has a 1 at the covariate of interest
  # and 0 elsewhere, and the control has all covariates 0.
  # In the second m/2 matched set, the case has 0 at all covariates and the control
  # has a 1 at the covariate of interest and 0 elsewhere.
  ms = data$matchedset; curMS = max(ms)
  zeros = rep(0,ncol(X))
  pseudoD = rep(c(1,0),times=m)
  for (i in 1:ncol(X)) { # loop over covariates
    D = c(D,pseudoD)
    pseudoX = zeros
    pseudoX[i] = 1 # 1 at the covariate of interest
    augX1=c(); augX2=c()
    for (i in 1:(m/2)) {
      augX1 = rbind(augX1,pseudoX,zeros) # add m/2 pairs 
      augX2 = rbind(augX2,zeros,pseudoX) # add m/2 pairs
    }
    X = rbind(X,augX1,augX2)
    ms = c(ms,curMS+rep(seq(1,m),each=2))
    curMS = max(ms)
  }
  
  # Step 3: Set up data.frame with null rownames and correct colnames.
  rownames(X) = NULL
  aug_data = data.frame(D,X,ms)
  names(aug_data) = c(all.vars(form),"matchedset")
  return(aug_data)
}

# log-F(m,m)-penalized conditional likelihood inference by general optimization
logFmatched = function(data,m) {
  # Input:
  # - dat is the data
  # - m is true value of m
  # Output: 
  # - coef is the estimator coefficient
  # - ci is the 95% confidence interval of the estimator
  
  logFloglklh = function(beta) {
    lkhd = 0
    for (i in 1:max(data$matchedset)) {  # i is the number of matched set
      matchedset_data = data %>% filter(matchedset==i)
      numerator = sum(matchedset_data[1,3:ncol(data)]*beta)
      denominator = sum(exp(matchedset_data[,3:ncol(data)]*beta))
      lkhd = lkhd+numerator-log(denominator)
    }
    f_pen = m/2*beta-m*log(1+exp(beta))
    pen_lkhd = lkhd+f_pen
    return(pen_lkhd)
  }
  
  opt = optimize(logFloglklh,c(-5,5),maximum=T)
  coef = opt$maximum
  lkhdDrop = function(beta) {
    2*(opt$objective-logFloglklh(beta))-qchisq(1-0.05,1)
  }
  ci = cbind(uniroot(lkhdDrop,c(-15,coef))$root,uniroot(lkhdDrop,c(coef,15))$root) 
  # the CI is based on the profile penalized likelihood
  return(list(coef=coef,ci=ci))
}



# Testing: compare to Firth penalty
# DES = read.csv("DES.csv")
# form = formula(case~DES+matern.smoke)
# source("clogitf.R")
# # clogitf() needs the matched set variable to be named "matchedset"
# DES$matchedset = DES$matched.set 
# fit = clogitf(form,DES,pl=TRUE)
# coefficients(fit)
# cbind(log(fit$ci.lower),log(fit$ci.upper))
# 
# DESaug = augment.logFmatched(form,DES,m=2)
# fit = clogitf(form,DESaug,pl=TRUE,penalty=0)
# coefficients(fit)
# cbind(log(fit$ci.lower),log(fit$ci.upper))

