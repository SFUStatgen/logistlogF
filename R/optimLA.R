optimLA<-function(alpha,m,X,y,beta_int=NULL,confounding_factor=NULL){
  ftracer=0
  tracer<-matrix(0, nrow=1,ncol=3)
  i=1
  if (is.null(confounding_factor)){
    for (i in 1:50){
      dlogPenalisedL<-function(beta){
        sum(X*y-(X*exp(alpha+beta*as.numeric(X))/(1+exp(alpha+beta*as.numeric(X)))))-m/2+m*exp(-beta)/(1+exp(-beta))
      }
      beta_max<-uniroot(dlogPenalisedL, c(-20,20))$root
      logLP_betamax<-function(alpha0){
        temp1<-sum(X^2*exp(alpha0+beta_max*X)/(1+exp(alpha0+beta_max*X)))-sum(X^2*(exp(alpha0+beta_max*X)/(1+exp(alpha0+beta_max*X)))^2)
        temp2<-exp(-beta_max)/(1+exp(-beta_max))-(exp(-beta_max)/(1+exp(-beta_max)))^2
        c=temp1+temp2
        LP<-sum(y*(alpha0+beta_max*as.numeric(X))-log(1+exp(alpha0+beta_max*as.numeric(X))))-log(beta(m/2,m/2))-m/2*beta_max-m*log(1+exp(-beta_max))+0.5*log(c)
        LP
      }
      opt.result<-optimize(f=logLP_betamax,lower=-10,upper=10,
                           maximum=TRUE, tol=0.0001)
      if(abs(opt.result$objective-ftracer)>=0.0001*abs(ftracer)){
        alpha=opt.result$maximum
        tracer<-rbind(tracer,c(alpha,NA,opt.result$objective))
        ftracer<-opt.result$objective
      }else{
        break
      }
    }
  }else{
    for (i in 1:50){
      dlogPenalisedL<-function(beta){
        sum(X*y-(X*exp(alpha+confounding_factor*beta_int+beta*as.numeric(X))/(1+exp(alpha+beta_int*confounding_factor+beta*as.numeric(X)))))-m/2+m*exp(-beta)/(1+exp(-beta))
      }
      beta_max<-uniroot(dlogPenalisedL, c(-20,20))$root
      logLP_betamax<-function(alpha0beta0){
        alpha0<-alpha0beta0[1]
        beta0<-alpha0beta0[2]
        temp1<-sum(X^2*exp(alpha0+confounding_factor*beta0+beta_max*X)/(1+exp(alpha0+confounding_factor*beta0+beta_max*X)))-sum(X^2*(exp(alpha0+confounding_factor*beta0+beta_max*X)/(1+exp(alpha0+confounding_factor*beta0+beta_max*X)))^2)
        temp2<-exp(-beta_max)/(1+exp(-beta_max))-(exp(-beta_max)/(1+exp(-beta_max)))^2
        c=temp1+temp2
        LP<-sum(y*(alpha0+confounding_factor*beta0+beta_max*as.numeric(X))-log(1+exp(alpha0+confounding_factor*beta0+beta_max*as.numeric(X))))-log(beta(m/2,m/2))-m/2*beta_max-m*log(1+exp(-beta_max))+0.5*log(c)
        LP
      }
      opt.result<-optim(par=c(alpha,beta_int),
                        fn=logLP_betamax,
                        method="L-BFGS-B",
                        control=list(fnscale=-1))
      if(abs(opt.result$value-ftracer)>=0.0001*abs(ftracer)){
        alpha=opt.result$par[1]
        beta_int=opt.result$par[2]
        tracer<-rbind(tracer,c(alpha,beta_int,opt.result$value))
        ftracer<-opt.result$value
      }else{
        break
      }
    }
  }
  tracer
}




