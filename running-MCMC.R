rm(list = ls())
setwd("/Users/yadavr/Dropbox (KAUST)/Paper2/Semi_conditional_model/Review_SpatialStats/Rcode/Simulation_GitHub/")

args <- commandArgs(TRUE)
for (arg in args) {
  eval(parse(text = arg))
  # print(arg)
}
rm(arg, args)
library(mvtnorm)
######### loading the data and source file
#load("Data.Rdata")
load("Data_nsites_20_ntime_100.Rdata")
#load("SimData_RealData_spaceDimesion_ntime_200.Rdata")
source("other-function.R")
source("update.param.R")
source("MCMC_main_function.R")
source("update_tuning_param_function.R")

############ deciding whethether to do simulation or data application
appltype<-"simulation"
if(appltype=="simulation"){
  true.param=c(log.tilde.param, log(X2)[1:4], log(X3)[1:4])
  param.names<-c(paste0("alpha","[",1:ncol(cov),"]"), paste0("beta","[",1:3,"]"), "rho", paste0("log(X2)","[",1:4,"]"), paste0("log(X3)","[",1:4,"]"))
} else if (appltype=="application"){
  true.param<-NULL
  param.names<-param.names<-c(paste0("alpha","[",1:n.cov,"]"), paste0("beta","[",1:n.cov,"]"), paste0("X2","[",1:n.cov,"]"), paste0("X3","[",1:n.cov,"]"))
}

### Types of proposals for scale

Alg.no<-"Alg1"

#### Algorithm 1
if(Alg.no=="Alg1"){
  Prop_type_1<-"MALA"
  Prop_type_2<-"MALA"
  Prop_type_3<-"MALA"
  Prop_type_4<-"MALA"
}

## Algorithm 2
if(Alg.no=="Alg2"){
Prop_type_1<-"RWM"
Prop_type_2<-"RWM"
Prop_type_3<-"RWM"
Prop_type_4<-"RWM"
} 



############## Fixed parameters
ntime.total<-ntime
delta<-10^(-5)
hyper.fixed<-c(10,2,2,1/3,1/100) ## fixed parameters in hyperparameter: sd of Gaussian in cov;  a and b in dbeta for beta1 and beta2; shape and rate in dgamma for beta3 and  rho
N.MCMC<- 1*10^5
burning1<-N.MCMC/4
burning2<-N.MCMC/2
adapt<-500
thin<-250  ## the thining 
sigma2.latentX2<-10^-6 # tuning parameter for mean/variance of the MALA proposals
sigma2.latentX3<-10^-6
sigma2.scale<-10^-6
sigma2.beta1<-10^-6
sigma2.beta2<-10^-4
sigma2.beta3.rho<-10^-4
theta.latentX2<-0.6
theta.latentX3<-0.6
theta.scale<-0.6
theta.beta1<-0.6
theta.beta2<-0.6
theta.beta3.rho<-0.6


######### Batch size 
batch<- 5
T.acc=floor(ntime.total/batch)

if(Prop_type_1=="MALA"){
  T.acc_1<-T.acc 
}
if(Prop_type_1=="RWM"){
  T.acc_1<-1
} 
if(Prop_type_2=="MALA"){
  T.acc_2<-T.acc
} 
if(Prop_type_2=="RWM"){
  T.acc_2<-1
} 
if (Prop_type_3=="MALA"){
  T.acc_3<-T.acc
} 
if (Prop_type_3=="RWM"){
  T.acc_3<-1
} 
if (Prop_type_4=="MALA"){
  T.acc_4<-T.acc
} 
if(Prop_type_4=="RWM"){
  T.acc_4<-1
}


epoch<-T.acc

## which chain to be used, for two different initial values
chains<-c("chain1","chain2")
chain<- 1
chain_s<- chains[chain]
if(chain_s=="chain1"){ 
  init<-init1
} else if (chain_s=="chain2") {
  init<-init2
}
############ The sequence after which we are updating the tuning parameter
n.block<-6
scale_adapt_seq<-seq(from=adapt, to=N.MCMC, by=n.block*adapt)
beta1_adapt_seq<-seq(from=2*adapt, to=N.MCMC, by=n.block*adapt)
beta2_adapt_seq<-seq(from=3*adapt, to=N.MCMC, by=n.block*adapt)
beta3.rho_adapt_seq<-seq(from=4*adapt, to=N.MCMC, by=n.block*adapt)
X2_adapt_seq<-seq(from=5*adapt, to=N.MCMC, by=n.block*adapt)
X3_adapt_seq<-seq(from=6*adapt, to=N.MCMC, by=n.block*adapt)

#### Running the MCMC algorithm
set.seed(1)
start_time1<-Sys.time()
MCMC.ouput<- MCMC_main_function(N.MCMC=N.MCMC, Y=Y, init=init, thin=thin, adapt=adapt, burning1 = burning1, burning2=burning2,
                         sigma2.latentX2=sigma2.latentX2, sigma2.latentX3=sigma2.latentX3, sigma2.scale=sigma2.scale, 
                         sigma2.beta1=sigma2.beta1, sigma2.beta2=sigma2.beta2, sigma2.beta3.rho=sigma2.beta3.rho,
                         theta.latentX2=theta.latentX2, theta.latentX3=theta.latentX3, theta.scale=theta.scale,
                         theta.beta1=theta.beta1, theta.beta2=theta.beta2, theta.beta3.rho=theta.beta3.rho,
                         cens.ind=cens.ind,  u.thr=u.thr, batch=batch, cov=cov, traceplots=TRUE, progress.show=TRUE, 
                         param.names=param.names, true.param=true.param, T.acc_1=T.acc_1, T.acc_2=T.acc_2, T.acc_3=T.acc_3, T.acc_4=T.acc_4,
                         Prop_type_1=Prop_type_1, Prop_type_2=Prop_type_2, Prop_type_3=Prop_type_3, Prop_type_4=Prop_type_4)
end_time1<-Sys.time()
run.time<-print(end_time1-start_time1)
 file.name<-paste0("MyData_Chain", chain, "_",batch,"_", Alg.no, ".RData")
 save.image(file.name)

