#rm(list=ls())
### loadrecquired library
library(mvtnorm)
library(fields)
#########################################################################
########### Main Examples of the product mixture model in Section 2.3 ####
#########################################################################

#Product mixture model: Y(s)= alpha(s) X1(s) * X2 * X3(s), where: gamma(x) is the gamma function
# X1(s): i.i.d Weibull(1/beta1,1/gamma(1+beta1)), with shape 1/beta1 and scale 1/gamma(1+beta1)
# X2 spatially constant values (varying with respect to time only) with marginal Weibull(1/beta2,1/gamma(1+beta2))
# X3(s): marginally Inv-Gamma(beta3,beta3-1), with scale beta3-1 >0 and shape beta3 > 1, dependence structure governed by Gaussian copula with exponential correlation function exp(-h/rho), h>=0, with range rho>0
# alpha(s): scale parameter of the model with exp(gamma0 + gamma1 * Z1(s) + gamma2 * Z2(s) + gamma3 * Z3(s)), where Zi(s), i=1,2,3 are the spatial covariates 

###################################
########### Simulation set-up ####
###################################
#DATA DIMENSION
nsites<-20#no of sites
ntime<-100# no of data at each site
n<-nsites*ntime #no of data points in space and time

#### PARAMETER COMBINATION
gamma0<-1   #intercept in alpha(s)
gamma1 <- 1  #first covariate in alpha(s)
gamma2 <- 1 #second covariate in alpha(s)
gamma3<-1 #third covariate in alpha(s)
beta1 <- 0.8 #Power for iid terms
beta2<-0.7 #Power for fully dependent terms
beta3 <- 5#gamma SHAPE for marginal distribution of X3, related to the tail index as xi=1/beta3 and scale beta3-1
rho <- 0.5 #RANGE parameter in exponential correlation function exp(-h/rho), h>=0, for the Gaussian copula in X2 (denominator)
## TILDE PARAMETRIZATION (LOG TRANSFROMATION TO MAKE THE SUPPORT ON REAL LINE) 
gamma0.tilde<-gamma0
gamma1.tilde<-gamma1
gamma2.tilde<-gamma2
gamma3.tilde<-gamma3
beta1.tilde<-log((beta1)/(1-beta1))
beta2.tilde<-log((beta2)/(1-beta2))
beta3.tilde<-log(1/(beta3-1))
rho.tilde<-log(rho/(2*rho.upper.range-rho))

#### SIMULATING THE SPATAIL LOCATIONS IN UNIFORM [0,1]^2 GRID 
set.seed(1)
loc<-cbind(runif(nsites),runif(nsites)) #locations generated uniformly in unit square [0,1]^2
Z1<-loc[,1] #first spatial covariate
Z2<-loc[,2] #second spatial covariate
Z3<-rnorm(nsites) #third spatial covariate
cov<-cbind(1,Z1,Z2,Z3) #degin matrix of dimesion 
#DISTANCE MATRIX BETWEEN OUR SITES
dist.mat <- as.matrix(dist(loc)) #distance matrix of all locations (dimension: nsites x nsites)
rho.upper.range<- 2*max(dist.mat) # range parameter the model

### STRING THE PARAMETERS IN A VECTORS
param <- c(gamma0,gamma1,gamma2,gamma3,beta1,beta2,beta3,rho) #Combines all hyperparameters in a single vector
log.tilde.param<-c(gamma0.tilde,gamma1.tilde,gamma2.tilde, gamma3.tilde, beta1.tilde,beta2.tilde,beta3.tilde, rho.tilde)

#PLOTTING SPATIAL LOCATIONS
par(mfrow=c(1,1))
plot(loc)
#######################
### DATA SIMULATION ###
#######################
#sim.data: function to simulated data
# INPUTS:
#   ntime: integer, number of independent time replicates
#   loc: matrix (dimension: nsites x 2), the i-th row gives the location of the i-th site
#   log.tilde.param: vector of tilde hyperparameters
#   cov: matrix(nsites x (p+1)), p:number of covariets
# OUTPUT: list with the following elements
#  Y: simulated data matrix (dimension: ntime x nsites), the i-th row gives the i-th temporal replicate at all the sites
#  X2: Simulated spatially constant and temporarily varying 
#  X3: Simulated matrix (dimension: ntime x nsites), the i-th row gives the i-th temporal replicate at all the sites
sim.data <- function(ntime, loc, log.tilde.param, cov){
  nsites <- nrow(loc) #number of sites
  n <- nsites*ntime #total number of observations in space and time
  n.cov<-ncol(cov)
  alpha_hyper<-matrix(log.tilde.param[1:n.cov], ncol=1)
  alpha<-matrix(rep(exp(cov%*%alpha_hyper),ntime),nrow = ntime,ncol = nsites,byrow = TRUE)
  
  beta1 <- exp(log.tilde.param[n.cov+1])/(1+exp(log.tilde.param[n.cov+1]))
  beta2 <- exp(log.tilde.param[n.cov+2])/(1+exp(log.tilde.param[n.cov+2]))
  beta3<-1+exp(-log.tilde.param[n.cov+3])
  rho<-2*rho.upper.range*exp(log.tilde.param[n.cov+4])/(1+exp(log.tilde.param[n.cov+4]))
  
  dist.mat <- as.matrix(dist(loc)) #distance matrix of all locations (dimension: nsites x nsites)
  Sigma <- exp(-dist.mat/rho) #correlation matrix for Gaussian copula
  
  X1<- ((matrix(rexp(n,rate=1),nrow=ntime,ncol=nsites)^beta1))/gamma(1+beta1)
  X2<-((matrix(rep(rexp(ntime,rate=1),times=nsites),nrow=ntime,ncol=nsites)^beta2))/gamma(1+beta2) 
  
  Gauss <-mvtnorm:: rmvnorm(n=ntime,mean=rep(0,nsites),Sigma) #simulation of ntime independent multivariate Gaussian vectors (dimension: ntime x nsites)
  X3<- (qgamma(pnorm(Gauss),shape=beta3,rate=1))/(beta3-1) #process X3(s) in denominator as independent vectors with gamma margins and Gaussian copula (dimension: ntime x nsites), this is further related to Inv-Gamma (as if X~Gamma, then 1/X ~ Inv-Gamma)
  Y <- (alpha*X1*X2)/X3 #data simulation as gamma ratio,
  return(list(Y=Y,X1=X1,X2=X2[,1],X3=X3))
}
set.seed(1)
data.simulated <- sim.data(ntime=ntime,loc=loc,log.tilde.param=log.tilde.param, cov = cov)

Y<-data.simulated$Y
X1<-data.simulated$X1
X2<-data.simulated$X2
X3<-data.simulated$X3
# 
# hist(Y, breaks = 100, prob=TRUE)
# lines(density(Y), col=2, lwd=2)
# 
mean(X1)
mean(X2)
mean(1/X3)

# ### Visulizing on the the spatail map for different parametrs combinatiation
# 
# #CHECK VALUES OF RANDOM FIELD SIMULATED ON A FINE GRID (20*20 grid)
# ngrid <- 30
# x.grid <- y.grid <- seq(0,1,length=ngrid)
# loc.grid <- as.matrix(expand.grid(x.grid,y.grid))
# data.simulated.grid <- sim.data(ntime=1,loc=loc.grid,log.tilde.param=log.tilde.param, cov=loc.grid)
# Y.grid <- data.simulated.grid$Y
# X1_iid_grid <- data.simulated.grid$X1
# X2_full_dep_grid<-data.simulated.grid$X2
# X3.grid <- data.simulated.grid$X3
# 
# #pdf(file="Project_presentation/Simbeta1tilde0.pdf",onefile=TRUE,width=11,height=5)
# par(mfrow=c(2,2),mgp=c(1.8,1,0),mar=c(4,3,1,1))
# image.plot(x.grid,y.grid,matrix(as.vector(log(t(Y.grid))),length(x.grid),length(y.grid)),main=expression(paste(Y(s), "on log scale",sep="")))
# image.plot(x.grid,y.grid,matrix(as.vector(log(X1_iid_grid)),length(x.grid),length(y.grid)),main=expression(paste(X[1][iid](s), "on log scale",sep="")))
# image.plot(x.grid,y.grid,matrix(as.vector(log(X2_full_dep_grid)),length(x.grid),length(y.grid)),main=expression(paste(X[2][full.dep](s), "on log scale",sep="")))
# image.plot(x.grid,y.grid,matrix(as.vector(log(t(X3.grid))),length(x.grid),length(y.grid)),main=expression(paste(X[3](s), "on log scale",sep="")))
# #dev.off()
# alpha<-matrix(rep(exp(param[1]+param[2]*Z1+param[3]*Z2),ntime),nrow = ntime,ncol = nsites,byrow = TRUE)
# ### QQplot to check if the data is simuplated coreectly marginally
# par(mfrow=c(1,2))
# qqplot(Y[,1],(alpha[,1]*(beta1.star+beta1.bar)/beta2)*qf(c(1:ntime)/(ntime+1),df1=2*(beta1.star+beta1.bar),df2=2*beta2),log="xy")
# abline(0,1,col="red")
# 
# #par(mfrow=c(1,1))
# true.q<-sort(Y[,1], decreasing = FALSE)
# fitted.q<-sort((alpha[,1]*(beta1.star+beta1.bar)/beta2)*qf(c(1:ntime)/(ntime+1), df1=2*(beta1.star+beta1.bar),df2=2*beta2), decreasing = FALSE)
# range.plot<-c(min(true.q,fitted.q),max(true.q,fitted.q))
# plot(true.q,fitted.q,xlim=range.plot,ylim=range.plot, xlab="Simulated quantile", ylab="Theroatical quantile (F dist.)")
# abline(0,1,col="red")


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

Y1<-Y
n.predict<-nsites/5
n.fit<-nsites-n.predict

prob<-0.50
fun.quan<-function(x){
  return(quantile(x[x>0],na.rm = TRUE,probs=prob))
}
quantile.matrix.fit<-matrix(rep(unname(apply(Y1[,-(1:n.predict)], 2, fun.quan)),ntime),nrow=ntime,ncol=n.fit,byrow = TRUE)

max.threhold<-10^2
## Creating the data full data Y (predict data + fit data) and threshold matrix u.thr (predict threshold + fit threshold)
## threhold for the stations where we fit the model 
u.thr.fit<-matrix(NA,nrow=ntime,ncol=n.fit)
## data for fitting
Y.fit<-matrix(NA,nrow=ntime,ncol=n.fit)
for(i in 1: ntime){
  for (j in 1: n.fit) {
    if(is.na(Y1[,-(1:n.predict)][i,j])==TRUE){
      u.thr.fit[i,j]<-max(Y1,na.rm=TRUE)+max.threhold
      Y.fit[i,j]<- NA#-10^6 
    } else{
      u.thr.fit[i,j]<-quantile.matrix.fit[i,j]
      Y.fit[i,j]<-Y1[,-(1:n.predict)][i,j]
    }
  }
}

## threshold for missing sites
u.thr.miss<-matrix(max(Y1,na.rm=TRUE)+max.threhold,nrow = ntime,ncol=n.predict)

## Threshold for all sites (predict+fit)
u.thr<-cbind(u.thr.miss,u.thr.fit)

## data for missing sites assumed to be +10^3 than the maximum oberved precipitation
data.to.predict<-matrix(NA,nrow = ntime,ncol=n.predict)

## data for all sites  (predict+fit)
Y<-cbind(data.to.predict,Y.fit)

#index for the cesnsored observations
ind1<-Y<u.thr
ind1[is.na(ind1)] <- TRUE

cens.ind<-ind1
rm(i)
rm(j)

############################################################################################
################ Initial-values ###############################################################
##############################################################################################
# INITIAL VALUES FOR FIRST CHAIN
gamma0.init1 <- 0.1 
gamma1.init1 <- 0
gamma2.init1 <- 0
gamma3.init1 <- 0
beta1.init1 <- 0.5
beta2.init1<-0.5
beta3.init1 <- 2.5
rho.init1<-0.1 

#tilde parametrization
gamma0.tilde.init1<-gamma0.init1
gamma1.tilde.init1<-gamma1.init1
gamma2.tilde.init1<-gamma2.init1
gamma3.tilde.init1<-gamma3.init1
beta1.tilde.init1<-log((beta1.init1)/(1-beta1.init1))
beta2.tilde.init1<-log((beta2.init1)/(1-beta2.init1))
beta3.tilde.init1<-log(1/(beta3.init1-1))
rho.tilde.init1<-log(rho.init1/(2*rho.upper.range-rho.init1))

tilde.param.init1<- c(gamma0.tilde.init1,gamma1.tilde.init1,gamma2.tilde.init1,gamma3.tilde.init1,
                      beta1.tilde.init1,
                      beta2.tilde.init1,
                      beta3.tilde.init1,
                      rho.tilde.init1) # "random" vector of all initial hyperparameters
set.seed(2)
#data.simulated.init1 <- sim.data(ntime=ntime,loc=loc,log.tilde.param=log.tilde.param, cov = cov)
data.simulated.init1 <- sim.data(ntime=ntime,loc=loc,log.tilde.param=tilde.param.init1, cov = cov)
# log.X2.init1<-log(data.simulated.init1$X2)
# log.X3.init1<-log(data.simulated.init1$X3)

log.X2.init1<-log(data.simulated.init1$X2)
log.X3.init1<-log(data.simulated.init1$X3)

init1<-c(tilde.param.init1, log.X2.init1, log.X3.init1)
#init1<-c(tilde.param.init1, log.X2.init1, log.X3.init1)


############################ Initial1 value for chain2

gamma0.init2<- 0.1 
gamma1.init2 <- 2
gamma2.init2 <- 2
gamma3.init2 <- 2
beta1.init2 <- 0.1 
beta2.init2<-0.1
beta3.init2 <- 10
rho.init2<-2 

# tilde parametrization
gamma0.tilde.init2<-gamma0.init2
gamma1.tilde.init2<-gamma1.init2
gamma2.tilde.init2<-gamma2.init2
gamma3.tilde.init2<-gamma3.init2
beta1.tilde.init2<-log((beta1.init2)/(1-beta1.init2))
beta2.tilde.init2<-log((beta2.init2)/(1-beta2.init2))
beta3.tilde.init2<-log(1/(beta3.init2-1))
rho.tilde.init2<-log(rho.init2/(2*rho.upper.range-rho.init2))

tilde.param.init2<- c(gamma0.tilde.init2, gamma1.tilde.init2, gamma2.tilde.init2, gamma3.tilde.init2,
                      beta1.tilde.init2,
                      beta2.tilde.init2,
                      beta3.tilde.init2,
                      rho.tilde.init2) # "random" vector of all initial hyperparameters
data.simulated.init2<- sim.data(ntime=ntime, loc=loc, log.tilde.param=tilde.param.init2, cov = cov)
log.X2.init2<-log(data.simulated.init2$X2)
log.X3.init2<-log(data.simulated.init2$X3)
init2<-c(tilde.param.init2, log.X2.init2, log.X3.init2)

 # save(Y, Y1, cov, cens.ind, log.tilde.param, X2, X3, init1, init2, ntime, nsites, u.thr, dist.mat=dist.mat, param, log.tilde.param, prob, rho.upper.range,
  #    file="/Users/yadavr/Dropbox (KAUST)/Paper2/Semi_conditional_model/Review_SpatialStats/Rcode/Simulation_GitHub/Data_nsites_20_ntime_100.Rdata")
