# log-F(m,m)-penalized likelihood inference (unconditional likelihood) by data augmentation



#' Logistic regression with log-F(m,m) penalty
#'
#' @param form R formula for the model
#' @param dat dataframe of data
#' @param m degree-of-freedome parameter
#' @param control control convergence
#'
#' @return The fitted penalized logistic regression
#' @export
#'
#' @examples
#' data(DES); DES$fmatched <- factor(DES$matched.set)
#' form = formula(case~fmatched+DES+matern.smoke)
#' fit = logF(form,DES,m=2)
#' coefficients(fit)
logistlogF = function(form,dat,m,control=glm.control()) {
  # form is an R formula, data is the data, m is the numerator and denominator
  # degrees of freedom for the log-F prior, control is algorithm control
  # arguments to be passed to glm().
  #---------------
  # Step 1: Extract (i) the response and (ii) the design matrix
  # from the input formula and data frame so that we can augment them.
  mf = model.frame(form,dat)
  D = model.response(mf)
  X = model.matrix(form,dat)
  xvars <- colnames(X)
  # Step 2 (augmentation): one pseudo-observation for each covariate,
  #  where the response is m/2 successes and m/2 failures (even if
  #  m is an odd number) and the covariates are all zeros except for
  #  a one indicating the index of the covariate.
  #  Following the recommendation of Greenland and Mansournia (2015; p. 3139)
  #  we do not penalize the intercept.
  n = rep(1,length(D))
  zeros = rep(0,ncol(X))
  pseudoD = m/2; pseudoN = m
  for(i in 2:ncol(X)) {
    D = c(D,pseudoD); n = c(n,pseudoN)
    pseudoX = zeros; pseudoX[i]=1; X = rbind(X,pseudoX)
  }
  # Step 3: Set up a response matrix with columns for number of successes
  # and number of failures.
  Y = cbind(D,n-D)
  # Step 4: Set up X's as a data.frame with null rownames and correct colnames.
  rownames(X) = NULL
  X = data.frame(X)
  # Seems that formulas ignore brackets, so  (Intercept) will look
  # for the variable Intercept. Remove brackets.
  xvars <- c("Intercept",xvars[-1])
  names(X) = xvars
  # Step 5: set up a formula and call glm()
  form = formula(paste("Y~ -1 + ",paste0(xvars,collapse="+")))
  out = glm(form,data=X,family=binomial(),control=control)
  return(out)
}

## Testing: compare to stratified Firth logistic regr
#library(logistf)
#DES = read.csv("DES.csv")
#DES$fmatched = factor(DES$matched.set)
#form = formula(case~fmatched+DES+matern.smoke)
#fit = logistf(form,data=DES)
#coefficients(fit)
#tem <- anova(fit,logistf(case~fmatched+matern.smoke,data=DES))
## Now logF
#fit = logF(form,DES,m=2)
#coefficients(fit)
