#' Title
#'
#' @param mvals
#' @param data
#' @param weight
#' @param Y_index
#' @param X_index
#' @param C_index
#' @param N
#' @param method
#' @param ini_alpha
#'
#' @return
#' @export
#'
#' @examples
get_m <- function(mvals,data,weight,Y_index,X_index,C_index,N,method="MCEM",ini_alpha=NULL){
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
