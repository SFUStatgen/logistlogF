#' Title
#'
#' @param mvals A vector of possible values of shrinkage parameter m to be considered
#' @param data A data.frame object containing response and covariates
#' @param weight The column name for the weights assigned to each covariate in data
#' @param Y_index The column index for the response
#' @param X_index The column index/indice for the covariate(s)
#' @param C_index The column index for the confounder in data
#' @param N The number of Monte Carlo replicates if method of "MCEM" is implemented
#' @param method The method to use, "MCEM"(default) or "LA"
#' @param ini_alpha The initial value for the intercept term
#'
#' @return the value of shrinkage parameter m with the largest approximate penalised likelihood
#' @export
#'
#' @examples
get_m <- function(mvals,data,weight=NULL,Y_index,X_index,C_index=NULL,N,method="MCEM",ini_alpha=NULL){
  require(parallel)
  # Input:
  # - mvals is a list of m values
  # - data is the GWAS data
  # - weight is a list of weights for SNVs
  # - Y_index is the index of phenotype in data
  # - X_index is the index of SNVs in data
  # - C_index is the index of confounders in data
  # - N is the number of Monte Carlo replicates
  # - method is either "MCEM" or "LA"
  # Output: an estimate of m
  if (method == "MCEM") {
    pp <- unlist(mclapply(X=mvals,FUN=profilelkhd,data=data,weight=weight,
                          Y_index=Y_index,X_index=X_index,C_index=C_index,N=N,
                          mc.preschedule=T,mc.cores=detectCores()))
    f <- function(x) {stats:::predict.smooth.spline(ss,x)$y}
    ss <- smooth.spline(pp)
    return(optimize(f,lower=first(mvals),upper=last(mvals),maximum=T)$maximum)}

  if (method == "LA"){
    logl_allm<-unlist(mclapply(X=mvals,FUN=LAapproxL,XM=as.matrix(data[,X_index]),
                               y=data[,Y_index],ini_alpha=ini_alpha,
                               confounding_factor=data[,C_index],weight_col=weight,mc.cores=detectCores()))
    return(mvals[which.max(logl_allm)])
  }
}
