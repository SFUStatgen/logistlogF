#' Title
#'
#' @param m
#' @param XM
#' @param y
#' @param ini_alpha
#' @param confounding_factor
#' @param weight_col
#'
#' @return
#' @export
#'
#' @examples
LAapproxL<-function(m,XM,y,ini_alpha,confounding_factor,weight_col){
    n<-nrow(XM)
    p<-ncol(XM)-1
    logl<-0
    if(is.null(ini_alpha)){ini_alpha=-4}
    ini.m<-m
    di=1
    for (di in 1:(p+1)){
      if(is.null(confounding_factor)){
      tracer1<-optimLA(ini_alpha,ini.m,XM[,di],y)}else{
      tracer1<-optimLA(ini_alpha,ini.m,XM[,di],y,0,confounding_factor)
      }
      if(is.null(weight_col)){
      logl<-logl+tracer1[nrow(tracer1),3]}else{
      logl<-logl+weight_col[di]*tracer1[nrow(tracer1),3]
      }
    }
  logl
}


