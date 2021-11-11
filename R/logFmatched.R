
#' log-F(m,m)-penalized conditional likelihood inference by data augmentation
#'
#' @param form
#' @param dat
#' @param m
#'
#' @return
#' @export
#'
#' @examples
logFmatched = function(form,dat,m=1) {
  # form is an R formula, dat is the data
  #---------------
  # Step 1: Extract (i) the response and (ii) the design matrix
  # from the input formula and data frame so that we can augment them.
  mf = model.frame(form,dat)
  D = model.response(mf)
  X = model.matrix(form,dat)
  if(ncol(X)==1) { #intercept only model, no augmentation needed
    return(X)
  } else {
    X = model.matrix(form,dat)[,-1,drop=FALSE] # we don't want the intercept
  }
  # Step 2 (augmentation): For each degree of freedom, add
  # two pseudo-matched sets of size 2 for each covariate.
  # In the first matched set, the case has a 1 at the covariate of interest
  # and 0 elsewhere, and the control has all covars 0.
  # In the second matched set, the case has 0 at all covars and the control
  # has a 1 at the covariate of interest and 0 elsewhere.
  ms = dat$matchedset; curMS = max(ms)
  zeros = rep(0,ncol(X))
  pseudoD = c(1,0,1,0)
  for(mm in 1:m) { # loop over degrees of freedom
    for(i in 1:ncol(X)) { # loop over covariates
      D = c(D,pseudoD)
      pseudoX = zeros; pseudoX[i]=1; X = rbind(X,pseudoX,zeros,zeros,pseudoX)
      ms = c(ms,rep(curMS+1,2),rep(curMS+2,2)); curMS = curMS+2
    }
  }
  # Step 3: Set up data.frame with null rownames and correct colnames.
  rownames(X) = NULL
  dat = data.frame(D,X,ms)
  names(dat) = c(all.vars(form),"matchedset")
  return(dat)
}

## Testing: compare to Firth penalty
#DES = read.csv("DES.csv")
#form = formula(case~DES+matern.smoke)
#source("clogitf.R")
## clogitf() needs the matched set variable to be named "matchedset"
#DES$matchedset = DES$matched.set
#fit = clogitf(form,DES,pl=TRUE)
#coefficients(fit)
#cbind(log(fit$ci.lower),log(fit$ci.upper))
#
#DESaug = augment.logFmatched(form,DES,m=2)
#fit = clogitf(form,DESaug,pl=TRUE,penalty=0)
#coefficients(fit)
#cbind(log(fit$ci.lower),log(fit$ci.upper))

