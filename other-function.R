###########################################
### USEFUL FUNCTIONS FOR MCMC ALGORITHM ###
###########################################

#CORRELATION MATRIX
#Sigma.function: function to calculate the correlation matrix of the latent process X2(s) at some fixed sites
# INPUTS:
#   rho: positive real number, range parameter of exponential correlation function exp(-h/rho), h>=0
#   dist.mat: distance matrix between locations (dimension: nsites x nsites)
# OUTPUT:
#   correlation matrix Sigma (dimension: nsites x nsites)
Sigma.function <- function(rho,dist.mat){ # correlation matrix for some fixed sites
  return(exp(-dist.mat/rho))
}

# PC priors for beta2
PC.prior.beta2 <- function(beta,lambda){
  if(beta>1){
    l.d<-sqrt(2)*lambda*exp(-sqrt(2)*lambda*((beta*(beta-1)))^(-1/2))*(beta-1/2)*(beta*(beta-1))^(-3/2)
  }else {l.d=10^(-50)}
  return(l.d)
}

# PC priors for beta1
KLD <- function(beta){
  return( (beta-1)*digamma(beta)-log(gamma(beta)) )
}
dKLD <- function(beta){
  return( (beta-1)*trigamma(beta) )
}
d.beta <- function(beta){
  return( sqrt(2*KLD(beta)) )
}
dd.beta <- function(beta){
  return( dKLD(beta)/sqrt(2*KLD(beta)) )
}
PC.prior.beta1 <- function(beta,lambda){
  res <- c()
  res[beta==0] <- 0
  res[beta==1] <- 0.5*lambda*mean(abs(c(dd.beta(1-10^(-6)),dd.beta(1+10^(-6)))))
  res[beta!=0 & beta!=1] <- 0.5*lambda*exp(-lambda*d.beta(beta[beta!=0 & beta!=1]))*abs(dd.beta(beta[beta!=0 & beta!=1])) ## The factor 0.5 is because we combine together the cases beta<1 and beta>1... (the PC prior is in fact a mixture between two priors defined over the intervals (0,1) and (1,Inf).)
  return( res )
}


# scale parameter in PC priors
lambda1=3  #for beta1
lambda2=3 #for beta2



