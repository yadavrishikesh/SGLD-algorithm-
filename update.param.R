log.posterior<-function(Y, X2, X3, param, dist.mat, cens.ind, u.thr, cov){
  ntime <- nrow(Y)
  nsites <- ncol(Y)
  n <- ntime*nsites
  n.cov<-ncol(cov)
  alpha_hyper<-matrix(param[1:n.cov], ncol=1)
  alpha<-matrix(rep(exp(cov%*%alpha_hyper),ntime),nrow = ntime,ncol = nsites,byrow = TRUE)
  
  beta1 <- param[n.cov+1]
  beta2 <- param[n.cov+2]
  beta3<-param[n.cov+3]
  rho<-param[n.cov+4]
  
  log.prior.alphas <- sum(dnorm(param[1:n.cov],mean=0,sd=hyper.fixed[1],log=TRUE)) # log-prior for alpha0
  log.prior.beta1<- dbeta(beta1,shape1=hyper.fixed[2],shape2=hyper.fixed[3],log=TRUE) # log-prior for beta1.star
  log.prior.beta2 <- dbeta(beta2,shape1 = hyper.fixed[2],shape2=hyper.fixed[3],log=TRUE) # log-prior for beta1.ba
  log.prior.beta3 <- dgamma(beta3,shape = hyper.fixed[4],rate=hyper.fixed[5],log=TRUE) # log-prior for beta2
  # log.prior.beta1 <- log(PC.prior.beta1(beta=beta1,lambda=lambda1)) # log-prior for beta1
  # log.prior.beta2<-log(PC.prior.beta2(beta=beta2,lambda=lambda2)) # log-prior for beta2
  log.prior.rho <- dunif(rho,min=0,max=2*rho.upper.range,log=TRUE) # log-prior for rho
  
  X2_matrix<-matrix(rep(X2,times=nsites),nrow=ntime,ncol=nsites)
  log.density.X1 <- sum(ifelse(cens.ind, pweibull(u.thr, shape = 1/beta1, scale = (alpha*X2_matrix)/(X3*gamma(1+beta1)),log.p=TRUE),
                               dweibull(Y, shape = 1/beta1, scale =  (alpha*X2_matrix)/(X3*gamma(1+beta1)), log = TRUE))) # log-density for X1(s) (numerator) given X2(s) (denominator)
  
  log.density.X2 <- sum(dweibull(X2, shape = 1/beta2, scale = 1/gamma(1+beta2), log = TRUE)) # log-density for X1.bar (numerator) given X2(s) (denominator)
  
  x3.arg<-qnorm(matrix(pgamma(X3, shape = beta3, rate = (beta3-1)),nrow = ntime, ncol = nsites)) ## quantile at which we need to calculate the Gaussian density
  log.density.X3 <- sum(mvtnorm::dmvnorm(x3.arg, mean=rep(0,nsites), sigma=Sigma.function(rho=rho,dist.mat=dist.mat),log = TRUE))+sum(dgamma(X3, shape = beta3, rate = beta3-1, log = TRUE))-sum(dnorm(x3.arg, log=TRUE)) # log-density for X2(s) (denominator)
  
  log.post<-log.prior.alphas+log.prior.beta1+log.prior.beta2+log.prior.beta3+log.prior.rho+log.density.X1+log.density.X2+log.density.X3
  return(log.post)
}

#Tilde.LOG-POSTERIOR DENSITY
# Calculating log posterior density while input latent variables are on log scale, i.e., log.X2=log(X2)
tilde.log.posterior<-function(Y, log.X2, log.X3, log.tilde.param, dist.mat,cens.ind, u.thr, cov){
  X2<-exp(log.X2)
  X3 <- exp(log.X3)
  n.cov<-ncol(cov)
  param <- c(log.tilde.param[1:n.cov],
             (exp(log.tilde.param[n.cov+1]))/(exp(log.tilde.param[n.cov+1])+1),
             (exp(log.tilde.param[n.cov+2]))/(exp(log.tilde.param[n.cov+2])+1),
             1+exp(-log.tilde.param[n.cov+3]),
             2*rho.upper.range*exp(log.tilde.param[n.cov+4])/(1+exp(log.tilde.param[n.cov+4])))
             
  return(log.posterior(Y, X2, X3, param, dist.mat,cens.ind, u.thr, cov) + sum(log.X2)+sum(log.X3) + log.tilde.param[n.cov+1] +log.tilde.param[n.cov+2]-
           2*log(1+exp(log.tilde.param[n.cov+1]))-2*log(1+exp(log.tilde.param[n.cov+2])) - 
           log.tilde.param[n.cov+3]+log(2)+log(rho.upper.range)+log.tilde.param[n.cov+4]-2*log(1+exp(log.tilde.param[n.cov+4])))
}


# # #To check the tilde.function is correctly parametrized
# log.X2<-log(X2)
# log.X3<-log(X3)
# tilde.log.posterior(Y, log.X2=log.X2, log.X3=log.X3, log.tilde.param, dist.mat,cens.ind, u.thr, cov=cov)
# log.posterior(Y, X2=X2, X3=X3, param, dist.mat,cens.ind, u.thr, cov=cov)


###############################################################################################################
###### tilde.log.postererior directly not using the log-posteruiors and using the unbiasedness in grdaients ########################
##########################################################################################################
# batch<-ntime
# ntime.total<-ntime
# SG_index<-sample(1:ntime, size=batch, replace=FALSE)

tilde.log.posterior1<-function(Y, log.X2, log.X3, log.tilde.param, dist.mat,cens.ind, u.thr, cov){
 ntime <- nrow(Y)
  nsites <- ncol(Y)
  n <- ntime*nsites
  n.cov<-ncol(cov)
  
  X2<-exp(log.X2)
  X3 <- exp(log.X3)
  alpha_hyper<-matrix(log.tilde.param[1:n.cov], ncol=1)
  alpha<-matrix(rep(exp(cov%*%alpha_hyper),ntime),nrow = ntime,ncol = nsites,byrow = TRUE)
  
  beta1 <-(exp(log.tilde.param[n.cov+1]))/(exp(log.tilde.param[n.cov+1])+1)
  beta2 <- (exp(log.tilde.param[n.cov+2]))/(exp(log.tilde.param[n.cov+2])+1)
  beta3<-1+exp(-log.tilde.param[n.cov+3])
  2*rho.upper.range*exp(log.tilde.param[n.cov+4])/(1+exp(log.tilde.param[n.cov+4]))
  
  log.prior.alphas <- sum(dnorm(param[1:n.cov],mean=0,sd=hyper.fixed[1],log=TRUE)) # log-prior for alpha0
  log.prior.beta1<- dbeta(beta1,shape1=hyper.fixed[2],shape2=hyper.fixed[3],log=TRUE) # log-prior for beta1.star
  log.prior.beta2 <- dbeta(beta2,shape1 = hyper.fixed[2],shape2=hyper.fixed[3],log=TRUE) # log-prior for beta1.ba
  log.prior.beta3 <- dgamma(beta3,shape = hyper.fixed[4],rate=hyper.fixed[5],log=TRUE) # log-prior for beta2
  log.prior.rho <- dunif(rho,min=0,max=2*rho.upper.range,log=TRUE) # log-prior for rho
  
  X2_matrix<-matrix(rep(X2,times=nsites),nrow=ntime,ncol=nsites)
  log.density.X1 <- sum(ifelse(cens.ind, pweibull(u.thr, shape = 1/beta1, scale = (alpha*X2_matrix)/(X3*gamma(1+beta1)),log.p=TRUE),
                               dweibull(Y, shape = 1/beta1, scale =  (alpha*X2_matrix)/(X3*gamma(1+beta1)), log = TRUE))) # log-density for X1(s) (numerator) given X2(s) (denominator)
  
  log.density.X2 <- sum(dweibull(X2, shape = 1/beta2, scale = 1/gamma(1+beta2), log = TRUE)) # log-density for X1.bar (numerator) given X2(s) (denominator)
  
  x3.arg<-qnorm(matrix(pgamma(X3, shape = beta3, rate = (beta3-1)),nrow = ntime, ncol = nsites)) ## quantile at which we need to calculate the Gaussian density
  log.density.X3 <- sum(mvtnorm::dmvnorm(x3.arg, mean=rep(0,nsites), sigma=Sigma.function(rho=rho,dist.mat=dist.mat),log = TRUE))+sum(dgamma(X3, shape = beta3, rate = beta3-1, log = TRUE))-sum(dnorm(x3.arg, log=TRUE)) # log-density for X2(s) (denominator)
  
  log.post<-log.prior.alphas+log.prior.beta1+log.prior.beta2+log.prior.beta3+log.prior.rho+(ntime.total/batch)*log.density.X1+(ntime.total/batch)*log.density.X2+(ntime.total/batch)*log.density.X3+
    sum(log.X2)+sum(log.X3) + log.tilde.param[n.cov+1] +log.tilde.param[n.cov+2]-2*log(1+exp(log.tilde.param[n.cov+1]))-2*log(1+exp(log.tilde.param[n.cov+2])) - log.tilde.param[n.cov+3]+
    log(2)+log(rho.upper.range)+log.tilde.param[n.cov+4]-2*log(1+exp(log.tilde.param[n.cov+4])) ########### Jacobean
  return(log.post)
}


# # #To check the tilde.function is correctly parametrized
# log.X2<-log(X2)
# log.X3<-log(X3)
# tilde.log.posterior1(Y=Y[SG_index,], log.X2=log.X2[SG_index], log.X3=log.X3[SG_index,], log.tilde.param, dist.mat,cens.ind[SG_index,], u.thr[SG_index,], cov=cov)
# #log.posterior(Y, X2=X2, X3=X3, param, dist.mat,cens.ind, u.thr, cov=cov)



################ Gradients of the log-posterior wrt the hyperparameters
tilde.gradient.log.posterior.numer.param<-function(Y, log.X2, log.X3, log.tilde.param, dist.mat, delta, cens.ind, u.thr, cov){
  ntime<-nrow(Y)
  nsites<-ncol(Y)
  n<-ntime*nsites
  n.cov<-ncol(cov)
  gradient.log.post <- c()
  ind_num<-c(1:length(log.tilde.param))
  for(j in 1:length(ind_num)){
    i<-ind_num[j]
    log.tilde.param.plus <- log.tilde.param.minus <- log.tilde.param
    log.tilde.param.plus[i] <- log.tilde.param[i]+delta
    log.tilde.param.minus[i] <- log.tilde.param[i]-delta
    gradient.log.post[j]<-(tilde.log.posterior1(Y=Y, log.X2=log.X2, log.X3=log.X3, log.tilde.param=log.tilde.param.plus, dist.mat,cens.ind, u.thr, cov)-
                             tilde.log.posterior1(Y=Y, log.X2=log.X2, log.X3=log.X3, log.tilde.param=log.tilde.param.minus, dist.mat,cens.ind, u.thr, cov))/(2*delta)
    #gradient.log.post[i] <- (log.posterior(Y=Y,X2=X2,X3=X3,param=param.plus,dist.mat=dist.mat,cens.ind, u.thr, cov)-log.posterior(Y=Y,X2=X2,X3=X3,param=param.minus,dist.mat=dist.mat,cens.ind, u.thr, cov))/(2*delta)
  }
  return(gradient.log.post)
}

#####################################################################################################
########## Prposal distribution for scale: covariets coeffiecents using MALA ########################
#####################################################################################################
### log-posterior wrt scale 
tilde.log.post.scale<-function(Y, log.X2, log.X3, log.tilde.param.scale, log.tilde.param.beta1, log.tilde.param.beta2, cens.ind, u.thr, cov, batch){
  ntime <- nrow(Y)
  nsites <- ncol(Y)
  n <- ntime*nsites
  n.cov<-ncol(cov)
  X2<-exp(log.X2)
  X3<-exp(log.X3)
  
  alpha_hyper<-matrix(log.tilde.param.scale, ncol=1)
  alpha<-matrix(rep(exp(cov%*%alpha_hyper),ntime),nrow = ntime,ncol = nsites,byrow = TRUE)
  
  beta1 <- exp(log.tilde.param.beta1)/(1+exp(log.tilde.param.beta1))
  beta2 <- exp(log.tilde.param.beta2)/(1+exp(log.tilde.param.beta2))
  
  log.prior.alphas <- sum(dnorm(param[1:n.cov],mean=0,sd=hyper.fixed[1],log=TRUE)) # log-prior for alpha0
  X2_matrix<-matrix(rep(X2,times=nsites),nrow=ntime,ncol=nsites)
  log.density.X1 <- sum(ifelse(cens.ind, pweibull(u.thr, shape = 1/beta1, scale = (alpha*X2_matrix)/(X3*gamma(1+beta1)),log.p=TRUE),
                               dweibull(Y, shape = 1/beta1, scale =  (alpha*X2_matrix)/(X3*gamma(1+beta1)), log = TRUE))) # log-density for X1(s) (numerator) given X2(s) (denominator)
  log.post<- log.prior.alphas + (ntime.total/batch) * log.density.X1 
  return(log.post) ### ntime.total/batch is basically the adjustements in the SGD steps
}

##### Gradients wrt the scale alpha
gradient.log.posterior.theor.sacle<-function(Y, log.X2, log.X3, log.tilde.param.scale, log.tilde.param.beta1, log.tilde.param.beta2, cens.ind, u.thr, cov, batch){
  ntime <- nrow(Y)
  nsites <- ncol(Y)
  n <- ntime*nsites
  n.cov<-ncol(cov)
  X2<-exp(log.X2)
  X3<-exp(log.X3)
  alpha_hyper<-matrix(log.tilde.param.scale, ncol=1)
  alpha<-matrix(rep(exp(cov%*%alpha_hyper),ntime),nrow = ntime,ncol = nsites,byrow = TRUE)
  
  beta1 <- exp(log.tilde.param.beta1)/(1+exp(log.tilde.param.beta1))
  beta2 <- exp(log.tilde.param.beta2)/(1+exp(log.tilde.param.beta2))
  X2_matrix<-matrix(rep(X2,times=nsites),nrow=ntime,ncol=nsites)
  
  log.gradient.X1<-ifelse(cens.ind, -u.thr*(dweibull(u.thr, shape= 1/beta1, scale=(alpha*X2_matrix)/(X3*gamma(1+beta1))))/(pweibull(u.thr, shape=1/beta1, scale=(alpha*X2_matrix)/(X3*gamma(1+beta1)))),
                          ((((Y*X3*gamma(1+beta1))/(alpha*X2_matrix))^(1/beta1))/beta1)-(1/beta1))
  grad<-((ntime.total/batch) * ((matrix(colSums(log.gradient.X1), nrow = 1)) %*% cov)) - t(alpha_hyper/hyper.fixed[1]^2)
  grad.est<-grad
  # if(max(grad.est)>(1/sigma2.scale)) {
  #   grad.est<-grad.est/max(grad.est)
  # } else{
  #   grad.est<-grad.est
  # }
  return(grad.est)
}
#delta<-10^(-4)

### Check unbiasedness of gradients for scale parameters
# B<-10000
# ntime.total<-ntime
# grad.scale<-matrix(NA, nrow = B, ncol = 4)
# for (i in 1: B) {
#   batch<-10
#   SG_index<-sample(1:ntime, size=batch, replace=FALSE)
# grad.scale[i ,]<-gradient.log.posterior.theor.sacle(Y=Y[SG_index,], log.X2=log(X2)[SG_index], log.X3=log(X3)[SG_index,], log.tilde.param.scale=log.tilde.param[1:ncol(cov)], log.tilde.param.beta1=log.tilde.param[ncol(cov)+1],
#                                    log.tilde.param.beta2=log.tilde.param[ncol(cov)+2], cens.ind[SG_index,], u.thr[SG_index,], cov, batch=batch)
# }
# 
# apply(grad.scale, MARGIN = 2, FUN = mean)
# ### true gradients 
# batch<-ntime
# SG_index<-sample(1:ntime, size=batch, replace=FALSE)
# tilde.gradient.log.posterior.numer.param(Y[SG_index,], log.X2=log(X2)[SG_index], log.X3=log(X3)[SG_index,], log.tilde.param=log.tilde.param, dist.mat, delta, cens.ind[SG_index,], u.thr[SG_index,], cov)[1:ncol(cov)]
# 

### prop distribution of beta1
proposal.mala.scale<-function(Y, log.X2, log.X3, log.tilde.param.scale, log.tilde.param.beta1, log.tilde.param.beta2, cens.ind, u.thr, cov, sigma2.scale, batch){
  ntime <- nrow(Y)
  nsites <- ncol(Y)
  n <- ntime*nsites
  means.mala <- log.tilde.param.scale+(sigma2.scale/2)*gradient.log.posterior.theor.sacle(Y=Y, log.X2=log.X2, log.X3=log.X3, log.tilde.param.scale=log.tilde.param.scale, log.tilde.param.beta1=log.tilde.param.beta1, 
                                                                                          log.tilde.param.beta2=log.tilde.param.beta2, cens.ind=cens.ind, u.thr=u.thr, cov=cov, batch=batch)
  variances.mala <- rep(sigma2.scale, length(log.tilde.param.scale))
  proposals<-rnorm(n=length(log.tilde.param.scale), mean=as.numeric(means.mala), sd=sqrt(variances.mala))
  return(proposals)
}

## Proposal density for beta1
log.proposal.density.mala.scale<-function(Y, log.X2, log.X3, log.tilde.param.scale, log.tilde.param.scale.star, log.tilde.param.beta1, log.tilde.param.beta2, cens.ind, u.thr, cov, sigma2.scale, batch){
  ntime <- nrow(Y)
  nsites <- ncol(Y)
  n <- ntime*nsites
  means.mala <- log.tilde.param.scale+(sigma2.scale/2)*gradient.log.posterior.theor.sacle(Y=Y, log.X2=log.X2, log.X3=log.X3, log.tilde.param.scale=log.tilde.param.scale, log.tilde.param.beta1=log.tilde.param.beta1, 
                                                                                          log.tilde.param.beta2=log.tilde.param.beta2, cens.ind=cens.ind, u.thr=u.thr, cov=cov, batch=batch)  
  variances.mala <- rep(sigma2.scale, length(log.tilde.param.scale))
  densities.mala<-dnorm(x=log.tilde.param.scale.star, mean=as.numeric(means.mala), sd=sqrt(variances.mala),log=TRUE)
  return(sum(densities.mala))
}
###############################################################################################
########## Prposal distribution for beta1 using MALA ##########################################
###############################################################################################
## log-posteriors for beta1
tilde.log.post.beta1<-function(Y, log.X2, log.X3, log.tilde.param.scale, log.tilde.param.beta1, cens.ind, u.thr, cov, batch){
  ntime <- nrow(Y)
  nsites <- ncol(Y)
  n <- ntime*nsites
  n.cov<-ncol(cov)
  X2<-exp(log.X2)
  X3 <- exp(log.X3)
  alpha_hyper<-matrix(log.tilde.param.scale, ncol=1)
  alpha<-matrix(rep(exp(cov%*%alpha_hyper),ntime),nrow = ntime,ncol = nsites,byrow = TRUE)
  beta1 <-  exp(log.tilde.param.beta1)/(exp(log.tilde.param.beta1)+1)
  
  log.prior.beta1<- dbeta(beta1, shape1=hyper.fixed[2],shape2=hyper.fixed[3],log=TRUE) # log-prior for beta1.star
  
  X2_matrix<-matrix(rep(X2,times=nsites),nrow=ntime,ncol=nsites)
  log.density.X1 <- sum(ifelse(cens.ind, pweibull(u.thr, shape = 1/beta1, scale = (alpha*X2_matrix)/(X3*gamma(1+beta1)),log.p=TRUE),
                               dweibull(Y, shape = 1/beta1, scale =  (alpha*X2_matrix)/(X3*gamma(1+beta1)), log = TRUE))) # log-density for X1(s) (numerator) given X2(s) (denominator)
  
  log.post<- (ntime.total/batch) * log.density.X1 + log.prior.beta1 + log.tilde.param.beta1 - 2*log(1+exp(log.tilde.param.beta1))
  return(log.post)
}

###Gradient of beta1
tilde.gradient.log.posterior.numer.beta1<-function(Y, log.X2, log.X3, log.tilde.param.scale,  log.tilde.param.beta1, delta, cens.ind, u.thr, cov, batch){
  ntime<-nrow(Y)
  nsites<-ncol(Y)
  n<-ntime*nsites
  n.cov<-ncol(cov)
  log.tilde.param.beta1.plus <- log.tilde.param.beta1+delta
  log.tilde.param.beta1.minus<- log.tilde.param.beta1-delta
  gradient.log.post<-(tilde.log.post.beta1(Y=Y, log.X2=log.X2, log.X3=log.X3, log.tilde.param.scale=log.tilde.param.scale, log.tilde.param.beta1=log.tilde.param.beta1.plus, cens.ind=cens.ind, u.thr=u.thr, cov=cov, batch=batch)-
                        tilde.log.post.beta1(Y=Y, log.X2=log.X2, log.X3=log.X3, log.tilde.param.scale=log.tilde.param.scale, log.tilde.param.beta1=log.tilde.param.beta1.minus, cens.ind=cens.ind, u.thr=u.thr, cov=cov, batch=batch))/(2*delta)
  
  return(gradient.log.post)
}
# 
# ############## Comparing if the gradient is coded well and it is unbiased 
# delta<-0.0001
# B<-10000
# ntime.total<-ntime
# grad.beta1<-c()
# for (i in 1: B) {
#   batch<-10
#   SG_index<-sample(1:ntime, size=batch, replace=FALSE)
#   grad.beta1[i]<-tilde.gradient.log.posterior.numer.beta1(Y[SG_index,], log.X2[SG_index], log.X3[SG_index,], log.tilde.param.scale=log.tilde.param[1:ncol(cov)], 
#                                                           log.tilde.param.beta1=log.tilde.param[ncol(cov)+1], delta, cens.ind[SG_index,], u.thr[SG_index,], cov, batch)
# }
# mean(grad.beta1)
# 
# batch<-ntime
# SG_index<-sample(1:ntime, size=batch, replace=FALSE)
# tilde.gradient.log.posterior.numer.beta1(Y[SG_index,], log.X2[SG_index], log.X3[SG_index,], log.tilde.param.scale=log.tilde.param[1:ncol(cov)], 
#                                          log.tilde.param.beta1=log.tilde.param[ncol(cov)+1], delta, cens.ind[SG_index,], u.thr[SG_index,], cov, batch)
# tilde.gradient.log.posterior.numer.param(Y[SG_index,], log.X2=log(X2)[SG_index], log.X3=log(X3)[SG_index,], log.tilde.param=log.tilde.param, dist.mat, delta,
#                                          cens.ind[SG_index,], u.thr[SG_index,], cov)[ncol(cov)+1]

### prop distribution of beta1
proposal.mala.beta1<-function(Y, log.X2, log.X3, log.tilde.param.scale, log.tilde.param.beta1, delta, cens.ind, u.thr, cov, sigma2.beta1, batch){
  ntime <- nrow(Y)
  nsites <- ncol(Y)
  n <- ntime*nsites
  means.mala <- log.tilde.param.beta1+(sigma2.beta1/2)*tilde.gradient.log.posterior.numer.beta1(Y=Y, log.X2=log.X2, log.X3=log.X3, log.tilde.param.scale=log.tilde.param.scale,
                                                                                                log.tilde.param.beta1=log.tilde.param.beta1, delta=delta, cens.ind=cens.ind, u.thr=u.thr, cov=cov, batch=batch)
  variances.mala <- rep(sigma2.beta1, length(log.tilde.param.beta1))
  proposals<-rnorm(n=length(log.tilde.param.beta1), mean=as.numeric(means.mala), sd=sqrt(variances.mala))
  return(proposals)
}

## Proposal density for beta1
log.proposal.density.mala.beta1<-function(Y, log.X2, log.X3, log.tilde.param.scale, log.tilde.param.beta1, log.tilde.param.beta1.star, delta, cens.ind, u.thr, cov, sigma2.beta1, batch){
  ntime <- nrow(Y)
  nsites <- ncol(Y)
  n <- ntime*nsites
  means.mala <- log.tilde.param.beta1+(sigma2.beta1/2)*tilde.gradient.log.posterior.numer.beta1(Y=Y, log.X2=log.X2, log.X3=log.X3, log.tilde.param.scale=log.tilde.param.scale,
                                                                                                log.tilde.param.beta1=log.tilde.param.beta1, delta=delta, cens.ind=cens.ind, u.thr=u.thr, cov=cov, batch=batch)
  variances.mala <- rep(sigma2.beta1, length(log.tilde.param.beta1))
  densities.mala<-dnorm(x=log.tilde.param.beta1.star, mean=as.numeric(means.mala), sd=sqrt(variances.mala),log=TRUE)
  return(sum(densities.mala))
}

# gradient.log.posterior.theor.beta1<-function(Y, X2, X3, param, dist.mat, cens.ind, u.thr, cov){
#   ntime <- nrow(Y)
#   nsites <- ncol(Y)
#   n <- ntime*nsites
#   n.cov<-ncol(cov)
#   alpha_hyper<-matrix(param[1:n.cov], ncol=1)
#   alpha<-matrix(rep(exp(cov%*%alpha_hyper),ntime),nrow = ntime,ncol = nsites,byrow = TRUE)
# 
#   beta1 <- param[n.cov+1]
#   beta2 <- param[n.cov+2]
#   beta3<-param[n.cov+3]
#   rho<-param[n.cov+4]
#   X2_matrix<-matrix(rep(X2,nsites),nrow = ntime,ncol = nsites)
#   ll_u<- u.thr*X3/alpha*X2_matrix
#   ll_y<- Y*X3/alpha*X2_matrix
#   #
#   # log.density.X1 <- sum(ifelse(cens.ind, pweibull(u.thr, shape = 1/beta1, scale = (alpha*X2_matrix)/(X3*gamma(1+beta1)),log.p=TRUE),
#   #                              dweibull(Y, shape = 1/beta1, scale =  (alpha*X2_matrix)/(X3*gamma(1+beta1)), log = TRUE))) # log-density for X1(s) (numerator) given X2(s) (denominator)
#   #grad_lik<-sum((-1/beta2)+ (((gamma(1+beta2)*X2)^(1/beta2))-1)*(1/beta2^2)*(log(X2)+ lgamma(1+beta2)-(beta2*digamma(1+beta2))))
#   log.gradient.lik<- sum(ifelse(cens.ind, (((ll_u*gamma(1+beta1))^(1/beta1))*(beta1*digamma(1+beta1)-log(ll_u)-lgamma(1+beta1)))/((beta1^2)*(exp((ll_u*gamma(1+beta1))^(1/beta1))-1)),
#                            (-1/beta1)+ ((((gamma(1+beta1)*ll_y)^(1/beta1))-1)*(1/beta1^2)*(log(ll_y)+ lgamma(1+beta1)-(beta1*digamma(1+beta1))))))
#   log.gradient.prior<-((hyper.fixed[2]-1)/beta1)-((hyper.fixed[3]-1)/(1-beta1))
#   grad<-log.gradient.lik+log.gradient.prior
#   return(grad)
# }
# delta<-10^(-4)
# gradient.log.posterior.theor.beta1(Y, X2, X3, param, dist.mat, cens.ind, u.thr, cov)
##tilde.gradient.log.posterior.numer.param(Y, log.X2=log(X2), log.X3=log(X3), log.tilde.param=log.tilde.param, dist.mat, delta, cens.ind, u.thr, cov)[5]

###############################################################################################
########## Prposal distribution for beta2 using MALA ##########################################
###############################################################################################

##### log-posteriors for beta2
tilde.log.post.beta2<-function(Y, log.X2, log.tilde.param.beta2, batch){
  ntime <- nrow(Y)
  nsites <- ncol(Y)
  n <- ntime*nsites
  n.cov<-ncol(cov)
  X2<-exp(log.X2)
  beta2 <-  exp(log.tilde.param.beta2)/(exp(log.tilde.param.beta2)+1)
  log.prior.beta2<- dbeta(beta2, shape1=hyper.fixed[2],shape2=hyper.fixed[3],log=TRUE) # log-prior for beta1.star
  log.density.X2 <- sum(dweibull(X2, shape = 1/beta2, scale = 1/gamma(1+beta2), log = TRUE)) # log-density for X1.bar (numerator) given X2(s) (denominator)
  log.post<-(ntime.total/batch) * log.density.X2 + log.prior.beta2 + log.tilde.param.beta2 - 2*log(1+exp(log.tilde.param.beta2))
  return(log.post)
}

### Gradient
gradient.log.posterior.theor.beta2<-function(Y, log.X2, log.tilde.param.beta2, batch){
  ntime <- nrow(Y)
  nsites <- ncol(Y)
  n <- ntime*nsites
  n.cov<-ncol(cov)
  X2<-exp(log.X2)
  beta2 <- exp(log.tilde.param.beta2)/(1+exp(log.tilde.param.beta2))
  
  grad_lik<-(ntime.total/batch) * (sum((-1/beta2)+ (((gamma(1+beta2)*X2)^(1/beta2))-1)*(1/beta2^2)*(log(X2)+ lgamma(1+beta2)-(beta2*digamma(1+beta2)))))
  grad_prior<-((hyper.fixed[2]-1)/beta2)-((hyper.fixed[3]-1)/(1-beta2))
  grad_jacob<-(1-exp(log.tilde.param.beta2))/(1+exp(log.tilde.param.beta2))
  grad<- (grad_lik+grad_prior) * (exp(log.tilde.param.beta2)/(1+exp(log.tilde.param.beta2))^2) + grad_jacob
  return(grad)
}
# delta<-10^(-4)
# 
# delta<-0.0001
# B<-10000
# ntime.total<-ntime
# grad.beta2<-c()
# for (i in 1: B) {
#   batch<-10
#   SG_index<-sample(1:ntime, size=batch, replace=FALSE)
#   grad.beta2[i]<-gradient.log.posterior.theor.beta2(Y[SG_index,], log.X2=log(X2)[SG_index], log.tilde.param.beta2=log.tilde.param[6], batch=batch)
# }
# mean(grad.beta2)
# batch<-ntime
# SG_index<-sample(1:ntime, size=batch, replace=FALSE)
# tilde.gradient.log.posterior.numer.param(Y, log.X2=log(X2)[SG_index], log.X3=log(X3)[SG_index,], log.tilde.param=log.tilde.param, 
#                                          dist.mat, delta, cens.ind[SG_index,], u.thr[SG_index,], cov)[6]


### prop distribution of beta2
proposal.mala.beta2<-function(Y, log.X2, log.tilde.param.beta2, sigma2.beta2, batch){
  ntime <- nrow(Y)
  nsites <- ncol(Y)
  n <- ntime*nsites
  means.mala <- log.tilde.param.beta2+(sigma2.beta2/2)*gradient.log.posterior.theor.beta2(Y=Y, log.X2=log.X2, log.tilde.param.beta2=log.tilde.param.beta2, batch=batch)
  variances.mala <- rep(sigma2.beta2, length(log.tilde.param.beta2))
  proposals<-rnorm(n=length(log.tilde.param.beta2), mean=as.numeric(means.mala), sd=sqrt(variances.mala))
  return(proposals)
}

## Proposal density for beta2
log.proposal.density.mala.beta2<-function(Y, log.X2, log.tilde.param.beta2, log.tilde.param.beta2.star, sigma2.beta2, batch){
  ntime <- nrow(Y)
  nsites <- ncol(Y)
  n <- ntime*nsites
  means.mala <- log.tilde.param.beta2+(sigma2.beta2/2)*gradient.log.posterior.theor.beta2(Y=Y, log.X2=log.X2, log.tilde.param.beta2=log.tilde.param.beta2, batch=batch)
  variances.mala <- rep(sigma2.beta2, length(log.tilde.param.beta2))
  densities.mala<-dnorm(x=log.tilde.param.beta2.star, mean=as.numeric(means.mala), sd=sqrt(variances.mala),log=TRUE)
  return(sum(densities.mala))
}


##############################################################################################################
########## Proposal distribution for  beta3 and rho ############################
##############################################################################################################

tilde.log.post.beta3.rho<-function(Y, log.X3, log.tilde.param.beta3,  log.tilde.param.rho, dist.mat, batch){
  ntime <- nrow(Y)
  nsites <- ncol(Y)
  n <- ntime*nsites
  n.cov<-ncol(cov)
  beta3<-1+exp(-log.tilde.param.beta3)
  rho<-2*rho.upper.range*exp(log.tilde.param.rho)/(1+exp(log.tilde.param.rho))
  X3<-exp(log.X3)
  
  log.prior.beta3 <- dgamma(beta3,shape = hyper.fixed[4],rate=hyper.fixed[5],log=TRUE) # log-prior for beta2
  log.prior.rho <- dunif(rho,min = 0,max=2*rho.upper.range,log=TRUE) # log-prior for rho
  
  x3.arg<-qnorm(matrix(pgamma(X3, shape = beta3, rate = (beta3-1)),nrow = ntime, ncol = nsites)) ## quantile at which we need to calculate the Gaussian density
  log.density.X3 <- sum(mvtnorm::dmvnorm(x3.arg, mean=rep(0,nsites), sigma=Sigma.function(rho=rho,dist.mat=dist.mat),log = TRUE))+sum(dgamma(X3, shape = beta3, rate = beta3-1, log = TRUE))-sum(dnorm(x3.arg, log=TRUE)) # log-density for X2(s) (denominator)
  log.post<- ((ntime.total/batch) * log.density.X3) + log.prior.beta3 + log.prior.rho - log.tilde.param.beta3 +
              log(2)+log(rho.upper.range)+log.tilde.param.rho-2*log(1+exp(log.tilde.param.rho))
  return(log.post)
  
}

# SG_index<-c(179, 199,  7,  13, 126,  22, 173,  48,   1,  34, 104,  31,  11,  96, 170,  85,  82,  54, 10,  49)
# Y_1<-Y[SG_index,]
# log.X3_1<-log(X3)[SG_index,]
# log.tilde.param.beta3<-log(1/(5.164566-1))
# log.tilde.param.rho<-log(1.439711/(2*rho.upper.range-1.439711))
# tilde.log.post.beta3.rho(Y=Y_1, log.X3=log.X3_1, log.tilde.param.beta3=log.tilde.param.beta3,  log.tilde.param.rho=log.tilde.param.rho,
#                          dist.mat, batch)


### Gradient of the log-post wrt beta3.tilde and rho.tilde
gradient.log.posterior.numer.beta3.rho<-function(Y, log.X3, log.tilde.param.beta3,  log.tilde.param.rho, dist.mat, delta, batch){
  ntime<-nrow(Y)
  nsites<-ncol(Y)
  n<-ntime*nsites
  n.cov<-ncol(cov)
  gradient.log.post <- c()
  log.tilde.param<-c(log.tilde.param.beta3, log.tilde.param.rho)
  for(j in 1:length(log.tilde.param)){
    log.tilde.param.plus <- log.tilde.param.minus <- log.tilde.param
    log.tilde.param.plus[j] <- log.tilde.param[j]+delta
    log.tilde.param.minus[j] <- log.tilde.param[j]-delta
    gradient.log.post[j]<-(tilde.log.post.beta3.rho(Y=Y, log.X3=log.X3, log.tilde.param.beta3=log.tilde.param.plus[1], log.tilde.param.rho=log.tilde.param.plus[2], dist.mat, batch)-
                             tilde.log.post.beta3.rho(Y=Y, log.X3=log.X3, log.tilde.param.beta3=log.tilde.param.minus[1], log.tilde.param.rho=log.tilde.param.minus[2], dist.mat, batch))/(2*delta)
    #gradient.log.post[i] <- (log.posterior(Y=Y,X2=X2,X3=X3,param=param.plus,dist.mat=dist.mat,cens.ind, u.thr, cov)-log.posterior(Y=Y,X2=X2,X3=X3,param=param.minus,dist.mat=dist.mat,cens.ind, u.thr, cov))/(2*delta)
  }
  return(gradient.log.post)
}

# delta<-0.0001
# B<-10000
# ntime.total<-ntime
# grad.beta3.rho<-matrix(NA, ncol = 2, nrow = 5)
# SG_index<-1:ntime
# for (i in 1: 5) {
#   batch<-10
#   SG_index<-((i-1)*batch+1):(i*batch)
#   #SG_index<-sample(1:ntime, size=batch, replace=FALSE)
#   grad.beta3.rho[i,]<-gradient.log.posterior.numer.beta3.rho(Y[SG_index,], log.X3[SG_index,], log.tilde.param.beta3=log.tilde.param[ncol(cov)+3], 
#                                                              log.tilde.param.rho=log.tilde.param[ncol(cov)+4], dist.mat, delta, batch)
# print(i)
# }
# apply(grad.beta3.rho, MARGIN = 2, FUN = mean)
# ##### True gradients based on full log-posteriors
# tilde.gradient.log.posterior.numer.param(Y, log.X2=log(X2), log.X3=log(X3), log.tilde.param=log.tilde.param,
#                                          dist.mat, delta, cens.ind, u.thr, cov)[c(ncol(cov)+3, ncol(cov)+4)]

### prop distribution of beta3 and rho
proposal.mala.beta3.rho<-function(Y, log.X3, log.tilde.param.beta3,  log.tilde.param.rho, dist.mat, delta, sigma2.beta3.rho, batch){
  ntime <- nrow(Y)
  nsites <- ncol(Y)
  n <- ntime*nsites
  log.tilde.param.beta3.rho<-c(log.tilde.param.beta3, log.tilde.param.rho)
  means.mala <- log.tilde.param.beta3.rho + (sigma2.beta3.rho/2) * gradient.log.posterior.numer.beta3.rho(Y=Y, log.X3=log.X3, log.tilde.param.beta3=log.tilde.param.beta3,  
                                                                                                          log.tilde.param.rho=log.tilde.param.rho, dist.mat=dist.mat, delta=delta, batch=batch)
  variances.mala <- rep(sigma2.beta3.rho, length(log.tilde.param.beta3.rho))
  proposals<-rnorm(n=length(log.tilde.param.beta3.rho), mean=as.numeric(means.mala), sd=sqrt(variances.mala))
  return(proposals)
}

## Proposal density for beta3 and rho
log.proposal.density.mala.beta3.rho<-function(Y, log.X3, log.tilde.param.beta3, log.tilde.param.beta3.star, log.tilde.param.rho, log.tilde.param.rho.star, dist.mat, delta, sigma2.beta3.rho, batch){
  ntime <- nrow(Y)
  nsites <- ncol(Y)
  log.tilde.param.beta3.rho<-c(log.tilde.param.beta3, log.tilde.param.rho)
  log.tilde.param.beta3.rho.star<-c(log.tilde.param.beta3.star, log.tilde.param.rho.star)
  means.mala <- log.tilde.param.beta3.rho + (sigma2.beta3.rho/2) * gradient.log.posterior.numer.beta3.rho(Y=Y, log.X3=log.X3, log.tilde.param.beta3=log.tilde.param.beta3,  
                                                                                                          log.tilde.param.rho=log.tilde.param.rho, dist.mat=dist.mat, delta=delta, batch=batch)
  variances.mala <- rep(sigma2.beta3.rho, length(log.tilde.param.beta3.rho))
  densities.mala<-dnorm(x=log.tilde.param.beta3.rho.star, mean=as.numeric(means.mala), sd=sqrt(variances.mala),log=TRUE)
  return(sum(densities.mala))
}

#################################################################################
################################## Proposal of log.X2 ###########################
#################################################################################
## log-post wrt X2
tilde.log.post.X2<-function(Y, log.X2, log.X3, log.tilde.param.scale, log.tilde.param.beta1, log.tilde.param.beta2, cens.ind, u.thr, cov){
  ntime <- nrow(Y)
  nsites <- ncol(Y)
  n <- ntime*nsites
  n.cov<-ncol(cov)
  X2<-exp(log.X2)
  X3<-exp(log.X3)
  alpha_hyper<-matrix(log.tilde.param.scale, ncol=1)
  alpha<-matrix(rep(exp(cov%*%alpha_hyper),ntime),nrow = ntime,ncol = nsites,byrow = TRUE)
  
  beta1 <- exp(log.tilde.param.beta1)/(1+exp(log.tilde.param.beta1))
  beta2 <- exp(log.tilde.param.beta2)/(1+exp(log.tilde.param.beta2))
  
  X2_matrix<-matrix(rep(X2,times=nsites),nrow=ntime,ncol=nsites)
  log.density.X1 <- sum(ifelse(cens.ind, pweibull(u.thr, shape = 1/beta1, scale = (alpha*X2_matrix)/(X3*gamma(1+beta1)),log.p=TRUE),
                               dweibull(Y, shape = 1/beta1, scale =  (alpha*X2_matrix)/(X3*gamma(1+beta1)), log = TRUE))) # log-density for X1(s) (numerator) given X2(s) (denominator)
  
  log.density.X2 <- sum(dweibull(X2, shape = 1/beta2, scale = 1/gamma(1+beta2), log = TRUE)) # log-density for X1.bar (numerator) given X2(s) (denominator)
  
  log.post<- log.density.X1 + log.density.X2 + sum(log.X2)
  return(log.post)
}

######### Gradient wrt X2 theoratical
tilde.gradient.log.posterior.theor.X2<-function(Y, log.X2, log.X3, log.tilde.param.scale, log.tilde.param.beta1, log.tilde.param.beta2, cens.ind, u.thr, cov){ 
  ntime <- nrow(Y)
  nsites <- ncol(Y)
  n <- ntime*nsites
  n.cov<-ncol(cov)
  X2<-exp(log.X2)
  X3<-exp(log.X3)
  alpha_hyper<-matrix(log.tilde.param.scale, ncol=1)
  alpha<-matrix(rep(exp(cov%*%alpha_hyper),ntime),nrow = ntime,ncol = nsites,byrow = TRUE)
  
  beta1 <- exp(log.tilde.param.beta1)/(1+exp(log.tilde.param.beta1))
  beta2 <- exp(log.tilde.param.beta2)/(1+exp(log.tilde.param.beta2))
  
  X2_matrix<-matrix(rep(X2,nsites),nrow = ntime,ncol = nsites)
  
  log.gradient.data_label<- apply(ifelse(cens.ind, -((gamma(1+beta1)*X3*u.thr)/(alpha*X2_matrix^2))*(dweibull(u.thr*X3*gamma(1+beta1)/(alpha*X2_matrix), shape=1/beta1, scale=1)/(pweibull(u.thr*X3*gamma(1+beta1)/(alpha*X2_matrix), shape=1/beta1, scale=1))),
                                         -(1/(beta1*X2_matrix)) + (((gamma(1+beta1)*Y*X3)/alpha)^(1/beta1))*(1/beta1)*((1/X2_matrix)^((1/beta1)+1))), MARGIN = 1, FUN = sum)
  log.gradient.X2_label<-((1/beta2)-1)*(1/X2)-((gamma(1+beta2))^(1/beta2))*(1/beta2)*(X2^((1/beta2)-1))
  
  log.grad<- (log.gradient.data_label+log.gradient.X2_label)*X2 + 1
  
  return(log.grad)
}

### Propsal for X2
proposal.mala.X2<-function(Y, log.X2, log.X3, log.tilde.param.scale, log.tilde.param.beta1, log.tilde.param.beta2, cens.ind, u.thr, cov, sigma2.latentX2){
  ntime <- nrow(Y)
  nsites <- ncol(Y)
  n <- ntime*nsites
  means.mala <- log.X2+(sigma2.latentX2/2)*tilde.gradient.log.posterior.theor.X2(Y=Y, log.X2=log.X2, log.X3=log.X3, log.tilde.param.scale=log.tilde.param.scale, log.tilde.param.beta1=log.tilde.param.beta1,
                                                                                 log.tilde.param.beta2=log.tilde.param.beta2, cens.ind=cens.ind, u.thr=u.thr, cov=cov)
  variances.mala <- rep(sigma2.latentX2, ntime)
  proposals<-rnorm(ntime,mean=as.numeric(means.mala),sd=sqrt(variances.mala))
  return(proposals)
}

## Reverse proposal of X2
log.proposal.density.mala.X2<-function(Y, log.X2, log.X2.star, log.X3, log.tilde.param.scale, log.tilde.param.beta1, log.tilde.param.beta2, cens.ind, u.thr, cov, sigma2.latentX2){
  ntime <- nrow(Y)
  nsites <- ncol(Y)
  n <- ntime*nsites
  means.mala <- log.X2+(sigma2.latentX2/2)*tilde.gradient.log.posterior.theor.X2(Y=Y, log.X2=log.X2, log.X3=log.X3, log.tilde.param.scale=log.tilde.param.scale, log.tilde.param.beta1=log.tilde.param.beta1,
                                                                                 log.tilde.param.beta2=log.tilde.param.beta2, cens.ind=cens.ind, u.thr=u.thr, cov=cov)
  variances.mala <- rep(sigma2.latentX2, ntime)
  densities.mala<-dnorm(x=log.X2.star,mean=as.numeric(means.mala),sd=sqrt(variances.mala),log=TRUE)
  return(sum(densities.mala))
}

######### Gradient wrt X2 numerical
tilde.gradient.log.posterior.numer.X2<-function(Y, log.X2, log.X3, log.tilde.param, dist.mat, delta,cens.ind, u.thr, cov){
  gradient.log.post <- c()

  for(i in 1:ntime){
    log.X2.plus <- log.X2.minus <- log.X2
    log.X2.plus[i] <- log.X2[i]+delta
    log.X2.minus[i] <-log.X2[i]-delta
    gradient.log.post[i] <- (tilde.log.posterior(Y=Y,log.X2=log.X2.plus,log.X3=log.X3,log.tilde.param=log.tilde.param,dist.mat=dist.mat,cens.ind, u.thr, cov)-
                               tilde.log.posterior(Y=Y,log.X2=log.X2.minus,log.X3=log.X3,log.tilde.param=log.tilde.param,dist.mat=dist.mat,cens.ind, u.thr, cov))/(2*delta)
  }

  return(gradient.log.post)
}
#To check above function (should give the same result)
# delta <- 10^(-4)
# batch<-10
# SG_index<-sample(1:ntime, size = batch, replace = TRUE)
# tilde.gradient.log.posterior.numer.X2(Y[SG_index,],log.X2[SG_index],log.X3=log(X3)[SG_index,],
#                                       log.tilde.param,dist.mat=dist.mat,delta,cens.ind[SG_index,], 
#                                       u.thr[SG_index,], cov)
# 
# ### full grdients based on the 
# tilde.gradient.log.posterior.theor.X2(Y,log.X2,log.X3=log(X3),log.tilde.param.scale=log.tilde.param[1:ncol(cov)],
#                                       log.tilde.param.beta1=log.tilde.param[ncol(cov)+1], log.tilde.param.beta2=log.tilde.param[ncol(cov)+2], 
#                                       cens.ind, u.thr, cov)[SG_index]

#################################################################################
################################## Proposal of log.X3 ###########################
#################################################################################
## log-post wrt X3
tilde.log.post.X3<-function(Y, log.X2, log.X3, log.tilde.param.scale, log.tilde.param.beta1, log.tilde.param.beta3, log.tilde.param.rho, dist.mat, cens.ind, u.thr, cov){
  ntime <- nrow(Y)
  nsites <- ncol(Y)
  n <- ntime*nsites
  n.cov<-ncol(cov)
  X2<-exp(log.X2)
  X3<-exp(log.X3)
  alpha_hyper<-matrix(log.tilde.param.scale, ncol=1)
  alpha<-matrix(rep(exp(cov%*%alpha_hyper),ntime),nrow = ntime,ncol = nsites,byrow = TRUE)
  
  beta1 <- exp(log.tilde.param.beta1)/(1+exp(log.tilde.param.beta1))
  beta3<-1+exp(-(log.tilde.param.beta3))
  rho<-2*rho.upper.range*exp(log.tilde.param.rho)/(1+exp(log.tilde.param.rho))
  
  X2_matrix<-matrix(rep(X2,times=nsites),nrow=ntime,ncol=nsites)
  log.density.X1 <- sum(ifelse(cens.ind, pweibull(u.thr, shape = 1/beta1, scale = (alpha*X2_matrix)/(X3*gamma(1+beta1)),log.p=TRUE),
                               dweibull(Y, shape = 1/beta1, scale =  (alpha*X2_matrix)/(X3*gamma(1+beta1)), log = TRUE))) # log-density for X1(s) (numerator) given X2(s) (denominator)
  
  
  x3.arg<-qnorm(matrix(pgamma(X3, shape = beta3, rate = (beta3-1)),nrow = ntime, ncol = nsites)) ## quantile at which we need to calculate the Gaussian density
  log.density.X3 <- sum(mvtnorm::dmvnorm(x3.arg, mean=rep(0,nsites), sigma=Sigma.function(rho=rho,dist.mat=dist.mat),log = TRUE)) + sum(dgamma(X3, shape = beta3, rate = beta3-1, log = TRUE)) - sum(dnorm(x3.arg, log=TRUE)) # log-density for X2(s) (denominator)
  
  
  log.post<- log.density.X1 + log.density.X3 + sum(log.X3)
  return(log.post)
}

##### Gradient wrt X3 theoratical
tilde.gradient.log.posterior.theor.X3<-function(Y, log.X2, log.X3, log.tilde.param.scale, log.tilde.param.beta1, log.tilde.param.beta3, log.tilde.param.rho, dist.mat, cens.ind, u.thr, cov){ 
  ntime <- nrow(Y)
  nsites <- ncol(Y)
  n <- ntime*nsites
  n.cov<-ncol(cov)
  X2<-exp(log.X2)
  X3<-exp(log.X3)
  alpha_hyper<-matrix(log.tilde.param.scale[1:n.cov], ncol=1)
  alpha<-matrix(rep(exp(cov%*%alpha_hyper),ntime),nrow = ntime,ncol = nsites,byrow = TRUE)
  
  beta1 <- exp(log.tilde.param.beta1)/(1+exp(log.tilde.param.beta1))
  beta3<-1+exp(-(log.tilde.param.beta3))
  rho<-2*rho.upper.range*exp(log.tilde.param.rho)/(1+exp(log.tilde.param.rho))
  
  X2_matrix<-matrix(rep(X2,times=nsites),nrow=ntime,ncol=nsites)
  
  log.gradient.X1<-ifelse(cens.ind, ((gamma(1+beta1)*(u.thr))/(alpha*X2_matrix))*(dweibull(u.thr*X3*gamma(1+beta1)/(alpha*X2_matrix), shape=1/beta1, scale=1)/(pweibull(u.thr*X3*gamma(1+beta1)/(alpha*X2_matrix), shape=1/beta1, scale=1))),
                          (1/(beta1*X3))-(1/beta1)*(X3^((1/beta1)-1))*((Y*gamma(1+beta1))/(alpha*X2_matrix))^(1/beta1))
  quantile.norm<-qnorm(pgamma(X3, shape = beta3, rate = (beta3-1)))
  grad.quantile.norm<-dgamma(X3,shape=beta3, rate=(beta3-1))/dnorm(quantile.norm)
  latent.X3<-(((beta3-1)/X3)-(beta3-1))-(grad.quantile.norm)*t(solve(Sigma.function(rho,dist.mat))%*%t(quantile.norm))+quantile.norm*grad.quantile.norm
  
  grad<-(log.gradient.X1+latent.X3)*X3 + 1
  return(grad)
}

####### proposal distribution of X3
proposal.mala.X3<-function(Y, log.X2, log.X3, log.tilde.param.scale, log.tilde.param.beta1, log.tilde.param.beta3, log.tilde.param.rho, dist.mat, cens.ind, u.thr, cov, sigma2.latentX3){
  ntime <- nrow(Y)
  nsites <- ncol(Y)
  n <- ntime*nsites
  means.mala <- log.X3+(sigma2.latentX3/2)*tilde.gradient.log.posterior.theor.X3(Y=Y, log.X2=log.X2, log.X3=log.X3, log.tilde.param.scale=log.tilde.param.scale, log.tilde.param.beta1=log.tilde.param.beta1,
                                                                                 log.tilde.param.beta3=log.tilde.param.beta3, log.tilde.param.rho=log.tilde.param.rho, dist.mat=dist.mat, cens.ind=cens.ind, u.thr=u.thr, cov=cov)
  variances.mala <- rep(sigma2.latentX3,n)
  proposals <- matrix(rnorm(n,mean=as.numeric(means.mala), sd=as.numeric(sqrt(variances.mala))), nrow=ntime,ncol=nsites)
  return(proposals)
}

####### reversal proposal density of X3
log.proposal.density.mala.X3<-function(Y, log.X2, log.X3, log.X3.star, log.tilde.param.scale, log.tilde.param.beta1, log.tilde.param.beta3, log.tilde.param.rho, dist.mat, cens.ind, u.thr, cov, sigma2.latentX3){
  ntime <- nrow(Y)
  nsites <- ncol(Y)
  n <- ntime*nsites
  means.mala <- log.X3+(sigma2.latentX3/2)*tilde.gradient.log.posterior.theor.X3(Y=Y, log.X2=log.X2, log.X3=log.X3, log.tilde.param.scale=log.tilde.param.scale, log.tilde.param.beta1=log.tilde.param.beta1,
                                                                                 log.tilde.param.beta3=log.tilde.param.beta3, log.tilde.param.rho=log.tilde.param.rho, dist.mat=dist.mat, cens.ind=cens.ind, u.thr=u.thr, cov=cov)
  variances.mala <- rep(sigma2.latentX3,n)
  densities.mala <- matrix(dnorm(x=log.X3.star,mean=as.numeric(means.mala),sd=as.numeric(sqrt(variances.mala)),log=TRUE),nrow=ntime,ncol=nsites)
  return(sum(densities.mala))
}

# ####### numeriCAL GRADIENT OF X3
tilde.gradient.log.posterior.numer.X3<-function(Y, log.X2, log.X3,log.tilde.param, dist.mat, delta, cens.ind, u.thr, cov){
  ntime<-nrow(Y)
  nsites<-ncol(Y)
  n<-ntime*nsites
  gradient.log.post <- matrix(nrow=ntime,ncol=nsites)
  for(i in 1:ntime){
    for(j in 1:nsites){
      log.X3.plus <- log.X3.minus <- log.X3
      log.X3.plus[i,j] <- log.X3[i,j]+delta
      log.X3.minus[i,j] <- log.X3[i,j]-delta
      gradient.log.post[i,j] <- (tilde.log.posterior(Y=Y, log.X2=log.X2, log.X3=log.X3.plus, log.tilde.param = log.tilde.param, dist.mat=dist.mat, cens.ind, u.thr, cov)-
                                   tilde.log.posterior(Y=Y, log.X2=log.X2, log.X3=log.X3.minus, log.tilde.param = log.tilde.param, dist.mat=dist.mat, cens.ind, u.thr, cov))/(2*delta)
    }
  }

  return(gradient.log.post)
}

# delta <- 10^(-4)
# 
# ## gradients based on stochasticity
# batch<-10
# SG_index<-sample(1:ntime, size = batch, replace = FALSE)
# tilde.gradient.log.posterior.numer.X3(Y[SG_index,],log.X2=log.X2[SG_index],log.X3=log(X3)[SG_index,],log.tilde.param,dist.mat,delta, cens.ind[SG_index,], u.thr[SG_index,], cov)
# ### True gradients
# tilde.gradient.log.posterior.theor.X3(Y,log.X2=log.X2,log.X3=log(X3),log.tilde.param.scale=log.tilde.param[1:ncol(cov)], log.tilde.param.beta1=log.tilde.param[1+ncol(cov)], 
#                                       log.tilde.param.beta3=log.tilde.param[3+ncol(cov)],log.tilde.param.rho=log.tilde.param[4+ncol(cov)], dist.mat, cens.ind, u.thr, cov)[SG_index,]
# 
# range(grad_diff.X3)
