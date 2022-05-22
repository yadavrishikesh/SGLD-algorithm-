MCMC_main_function=function(N.MCMC, Y, init, thin, adapt, burning1, burning2, sigma2.latentX2, sigma2.latentX3, sigma2.scale, 
                            sigma2.beta1, sigma2.beta2, sigma2.beta3.rho, theta.latentX2, theta.latentX3, theta.scale,
                            theta.beta1, theta.beta2, theta.beta3.rho, cens.ind,  u.thr, batch, cov, traceplots, progress.show, 
                            param.names, true.param, T.acc_1, T.acc_2, T.acc_3, T.acc_4, Prop_type_1, Prop_type_2, Prop_type_3, Prop_type_4)
{
  ntime <- nrow(Y)
  nsites <- ncol(Y)
  n <- ntime*nsites
  n.cov<-ncol(cov)
  samples <- matrix(nrow=floor(N.MCMC/thin),ncol=n.cov+4+8,byrow=TRUE)
  sigma.matrix<-matrix(nrow=floor(N.MCMC/adapt),ncol=6,byrow=TRUE)
  sigma.matrix[1,]<-c(sigma2.scale, sigma2.beta1, sigma2.beta2, sigma2.beta3.rho, sigma2.latentX2, sigma2.latentX3)
  samples[1,]<-c(init[1:(n.cov+4)], init[(n.cov+5):(n.cov+8)],init[(ntime+n.cov+5):(ntime+n.cov+8)])
  
  cur.samples.scale <- init[1:n.cov]
  cur.samples.beta1 <- init[n.cov+1]
  cur.samples.beta2 <- init[n.cov+2]
  cur.samples.beta3 <- init[n.cov+3]
  cur.samples.rho<- init[n.cov+4]
  cur.samples.latentX2<-init[(n.cov+5):(ntime+n.cov+4)]
  cur.samples.latentX3<-init[(ntime+n.cov+5):(n+ntime+n.cov+4)]
  
  
  cur.samples.scale.zero<-cur.samples.scale
  cur.samples.beta1.zero<-cur.samples.beta1
  cur.samples.beta2.zero<-cur.samples.beta2
  cur.samples.beta3.zero<-cur.samples.beta3
  cur.samples.rho.zero<-cur.samples.rho
  
  j<-1
  l<-1
  m<-1
  k<-1
  
  for (i in 1:(N.MCMC-1)){
    if((i%%thin)-1==0){
      par(mfrow=c(4,4),oma=c(0,0,2,0),mar=c(4,5,1,1))
      if((i%%(adapt))-1==0 & i< (burning1+burning2+2)){ #to calculate the accpetance rate based on only current samples, burning+2 to calculate the accpertance rate after the burning samples
        rate.scale<-0;  rate.beta1<-0;  rate.beta2<-0;  rate.beta3.rho<-0;  rate.latentX2<-0; rate.latentX3<-0
      }
      ####### Show the progress ######
      if(progress.show==TRUE){
        if(i< (burning1+burning2+2)){
          print(paste("Iteration: ",i, "; Acc rate scale=",T.acc_1*rate.scale/((i%%(adapt))+1),"; Acc rate beta1=",T.acc_2*rate.beta1/((i%%(adapt))+1), "; Acc rate beta2=",T.acc_3*rate.beta2/((i%%(adapt))+1),
                      "; Acc rate beta3. rho=",T.acc_4*rate.beta3.rho/((i%%(adapt))+1),"; Acc rate latX2=",rate.latentX2/((i%%(adapt))+1),"; Acc rate latX3=",rate.latentX3/((i%%(adapt))+1),
                      "; sigma scale=", sigma2.scale, "; sigma beta1=", sigma2.beta1, "; sigma beta2=", sigma2.beta2,"; sigma beta3 rho=", sigma2.beta3.rho, 
                      "; sigma latX2=",sigma2.latentX2,"; sigma latX3=", sigma2.latentX3,sep=""))
        } else{
          print(paste("Iteration: ",i,"; Acceptance rate scale=",T.acc_1*rate.scale/(i-(burning1+burning2+2)),"; Acceptance rate  beta1=",T.acc_2*rate.beta1/(i-(burning1+burning2+2)),
                      "; Acceptance rate beta2=",T.acc_3*rate.beta2/(i-(burning1+burning2+2)),"; Acceptance rate beta3.rho=",T.acc_4*rate.beta3.rho/(i-(burning1+burning2+2)),
                      "; Acceptance rate latX2=",rate.latentX2/(i-(burning1+burning2+2)),"; Acceptance rate latX3=",rate.latentX3/(i-(burning1+burning2+2)),
                      "; sigma scale=", sigma2.scale, "; sigma beta1=", sigma2.beta1, "; sigma beta2=", sigma2.beta2,"; sigma beta3 rho=", sigma2.beta3.rho, 
                      "; sigma latX2=",sigma2.latentX2,"; sigma latX3=", sigma2.latentX3,sep=""))
        }
      }
      
      if(traceplots==TRUE){
        for (ll in 1: ncol(samples)) {
          plot(thin*c(0:(l-1))+1,samples[1:l,ll],type = "l",xlab="MCMC iteration",ylab=param.names[ll]) # Plot for alpha.tilde
          abline(h=true.param[ll],col="2")
        }
      }
      
      l<-l+1
    }
    
    ####### The SGD index ######
    #ind_sample<-sample(1:(ntime-batch+1), size=1)
    #SG_index<- (ind_sample) : ((ind_sample+ batch)-1)
    SG_index<- sample(1:ntime, size=batch, replace = FALSE)
    # ind.SGD<-i%%epoch
    # if(ind.SGD==0){ind.SGD=epoch}
    # SG_index<-((ind.SGD-1)*batch+1):(ind.SGD*batch)
    # cur.samples.latentX3<-c(log(X3))
    # cur.samples.latentX2<-c(log(X2))
    # cur.samples.scale<-log.tilde.param[1:n.cov]
    # cur.samples.beta1<-log.tilde.param[1+n.cov]
    # cur.samples.beta2<-log.tilde.param[2+n.cov]
    
    
    
    ###### Current parameters #######
    cur.log.scale<- cur.samples.scale 
    cur.log.beta1 <- cur.samples.beta1 
    cur.log.beta2 <- cur.samples.beta2 
    cur.log.beta3 <- cur.samples.beta3
    cur.log.rho <- cur.samples.rho
    cur.log.X3<-matrix(cur.samples.latentX3, nrow = ntime, ncol = nsites) # Current latent parameters X2(s) (on log scale)
    cur.log.X2<-cur.samples.latentX2
    
    cur.samples.scale.zero<-cur.samples.scale.zero ### the samples which is stored after the fixed number of steps
    cur.samples.beta1.zero<-cur.samples.beta1.zero ### the samples which is stored after the fixed number of steps
    cur.samples.beta2.zero<- cur.samples.beta2.zero ### the samples which is stored after the fixed number of steps
    cur.samples.beta3.zero<-cur.samples.beta3.zero ### the samples which is stored after the fixed number of steps
    cur.samples.rho.zero<-cur.samples.rho.zero ### the samples which is stored after the fixed number of steps
    
    ##########  Updating scale parameters #########
    if(Prop_type_1=="MALA"){
      prop.log.scale<-proposal.mala.scale(Y=Y[SG_index,], log.X2 = cur.log.X2[SG_index], log.X3 = cur.log.X3[SG_index,], log.tilde.param.scale = cur.log.scale, 
                                          log.tilde.param.beta1 = cur.log.beta1, log.tilde.param.beta2 = cur.log.beta2, cens.ind = cens.ind[SG_index,],
                                          u.thr = u.thr[SG_index,], cov=cov, sigma2.scale = sigma2.scale, batch=batch)
    } else if(Prop_type_1=="RWM"){
      prop.log.scale<-rnorm(n=length(cur.log.scale), mean=cur.log.scale, sd=sqrt(sigma2.scale))
    }
    cur.samples.scale<- c(prop.log.scale)
    
    if(i%%T.acc_1==0){ ## calculating the log-posteriors based on sample
      cur.log.posterior.scale<- tilde.log.post.scale(Y=Y, log.X2 = cur.log.X2, log.X3 = cur.log.X3, log.tilde.param.scale = cur.samples.scale.zero, 
                                                     log.tilde.param.beta1 = cur.log.beta1, log.tilde.param.beta2 = cur.log.beta2, cens.ind = cens.ind, 
                                                     u.thr = u.thr, cov=cov, batch=ntime.total)
      prop.log.posterior.scale <-  tilde.log.post.scale(Y=Y, log.X2 = cur.log.X2, log.X3 = cur.log.X3, log.tilde.param.scale = cur.samples.scale, 
                                                        log.tilde.param.beta1 = cur.log.beta1, log.tilde.param.beta2 = cur.log.beta2, cens.ind = cens.ind, 
                                                        u.thr = u.thr, cov=cov, batch=ntime.total)
      
      if(Prop_type_1=="MALA"){
        sum.cur.log.proposal.mala.scale.batch<-log.proposal.density.mala.scale(Y=Y, log.X2 = cur.log.X2, log.X3 = cur.log.X3, log.tilde.param.scale = prop.log.scale, 
                                                                               log.tilde.param.scale.star=cur.samples.scale.zero, log.tilde.param.beta1 = cur.log.beta1, log.tilde.param.beta2 = cur.log.beta2, 
                                                                               cens.ind = cens.ind, u.thr = u.thr, cov=cov, sigma2.scale=sigma2.scale, batch=ntime.total)
        sum.prop.log.proposal.mala.scale.batch<-log.proposal.density.mala.scale(Y=Y, log.X2 = cur.log.X2, log.X3 = cur.log.X3, log.tilde.param.scale = cur.samples.scale.zero, 
                                                                                log.tilde.param.scale.star=prop.log.scale, log.tilde.param.beta1 = cur.log.beta1, log.tilde.param.beta2 = cur.log.beta2, 
                                                                                cens.ind = cens.ind, u.thr = u.thr, cov=cov, sigma2.scale=sigma2.scale, batch=ntime.total)
        
        log.ratio.scale <- prop.log.posterior.scale + sum.cur.log.proposal.mala.scale.batch - cur.log.posterior.scale - sum.prop.log.proposal.mala.scale.batch
      } else if(Prop_type_1=="RWM"){
        log.ratio.scale<-prop.log.posterior.scale-cur.log.posterior.scale
      }
      
      if (log(runif(1)) < log.ratio.scale & is.na(log.ratio.scale)==FALSE){
        cur.samples.scale.zero<-cur.samples.scale
        rate.scale<-rate.scale+1
      } else{
        cur.samples.scale<- cur.samples.scale.zero
      }
    }
    
    sigma2.scale<-adpative_function(index_MCMC_iter = i, sigma2_adapt = sigma2.scale, target_accept = 0.57, rate_adapt = rate.scale, burning1 = burning1,
                                    burning2 = burning2, adapt_seq = scale_adapt_seq, adapt = adapt/T.acc_1, adpat_param = theta.scale, lower.acc = 0.5, upper.acc = 0.65)
    
    
    
    ############  Updating beta1 parameter ############
    if(Prop_type_2=="MALA"){
      prop.log.beta1<-proposal.mala.beta1(Y=Y[SG_index,], log.X2 = cur.log.X2[SG_index], log.X3 = cur.log.X3[SG_index,], log.tilde.param.scale  = cur.samples.scale.zero,
                                          log.tilde.param.beta1= cur.log.beta1, delta=delta, cens.ind=cens.ind[SG_index,], u.thr=u.thr[SG_index], cov=cov, sigma2.beta1=sigma2.beta1, batch=batch)
    } else if(Prop_type_2=="RWM"){
      prop.log.beta1<-rnorm(n=1, mean=cur.log.beta1, sd=sqrt(sigma2.beta1))
    }
    
    cur.samples.beta1<- prop.log.beta1
    
    if(i%%T.acc_2==0){ ## calculating the log-posteriors based on sample
      cur.log.posterior.beta1<- tilde.log.post.beta1(Y=Y, log.X2 = cur.log.X2, log.X3 = cur.log.X3, log.tilde.param.scale  = cur.samples.scale.zero,
                                                     log.tilde.param.beta1= cur.samples.beta1.zero, cens.ind=cens.ind, u.thr=u.thr, cov=cov, batch=ntime.total)
      prop.log.posterior.beta1 <- tilde.log.post.beta1(Y=Y, log.X2 = cur.log.X2, log.X3 = cur.log.X3, log.tilde.param.scale  = cur.samples.scale.zero,
                                                       log.tilde.param.beta1= prop.log.beta1, cens.ind=cens.ind, u.thr=u.thr, cov=cov, batch=ntime.total)
      
      #       cur.log.posterior.beta1<-tilde.log.posterior(Y=Y, log.X2=cur.log.X2, log.X3 = cur.log.X3, log.tilde.param=c(cur.samples.scale.zero, cur.samples.beta1.zero, cur.samples.beta2.zero,
      #                                                                                                                       cur.samples.beta3.zero, cur.samples.rho.zero), dist.mat=dist.mat, cens.ind=cens.ind, u.thr=u.thr, cov=cov )
      #       prop.log.posterior.beta1<-tilde.log.posterior(Y=Y, log.X2=cur.log.X2, log.X3 = cur.log.X3, log.tilde.param=c(cur.samples.scale.zero, prop.log.beta1, cur.samples.beta2.zero,
      #                                                                                                                    cur.samples.beta3.zero,  cur.samples.rho.zero), dist.mat=dist.mat, cens.ind=cens.ind, u.thr=u.thr, cov=cov )
      # # 
      if(Prop_type_2=="MALA"){
        sum.prop.log.proposal.mala.beta1.batch <-  log.proposal.density.mala.beta1(Y=Y, log.X2 = cur.log.X2, log.X3 = cur.log.X3, log.tilde.param.scale  = cur.samples.scale.zero,
                                                                                   log.tilde.param.beta1= cur.samples.beta1.zero, log.tilde.param.beta1.star=prop.log.beta1, delta=delta, cens.ind=cens.ind,
                                                                                   u.thr=u.thr, cov=cov, sigma2.beta1=sigma2.beta1, batch=ntime.total)
        sum.cur.log.proposal.mala.beta1.batch <-  log.proposal.density.mala.beta1(Y=Y, log.X2 = cur.log.X2, log.X3 = cur.log.X3, log.tilde.param.scale  = cur.samples.scale.zero,
                                                                                  log.tilde.param.beta1= prop.log.beta1, log.tilde.param.beta1.star=cur.samples.beta1.zero, delta=delta, cens.ind=cens.ind,
                                                                                  u.thr=u.thr, cov=cov, sigma2.beta1=sigma2.beta1, batch=ntime.total)
        log.ratio.beta1<- prop.log.posterior.beta1 + sum.cur.log.proposal.mala.beta1.batch - cur.log.posterior.beta1 - sum.prop.log.proposal.mala.beta1.batch
      } else if(Prop_type_2=="RWM"){
        log.ratio.beta1<- prop.log.posterior.beta1  - cur.log.posterior.beta1 
      }
      
      
      if (log(runif(1)) < log.ratio.beta1 & is.na(log.ratio.beta1)==FALSE){
        cur.samples.beta1.zero<-cur.samples.beta1
        rate.beta1<-rate.beta1+1
      } else{
        cur.samples.beta1<- cur.samples.beta1.zero
        
      }
    }
    sigma2.beta1<-adpative_function(index_MCMC_iter = i, sigma2_adapt = sigma2.beta1, target_accept = 0.57, rate_adapt = rate.beta1, burning1 = burning1,
                                    burning2 = burning2, adapt_seq = beta1_adapt_seq, adapt = adapt/T.acc_2, adpat_param = theta.beta1, lower.acc = 0.50, upper.acc = 0.65)
    
    ############  Updating beta2 parameter ############
    
    if(Prop_type_3=="MALA"){
      prop.log.beta2<-proposal.mala.beta2(Y=Y[SG_index,], log.X2 = cur.log.X2[SG_index],  log.tilde.param.beta2=cur.log.beta2, sigma2.beta2=sigma2.beta2, batch=batch)
    } else if(Prop_type_3=="RWM"){
      prop.log.beta2<-rnorm(n=1, mean = cur.log.beta2, sd=sqrt(sigma2.beta2))
    }
    
    cur.samples.beta2<- prop.log.beta2
    
    if(i%%T.acc_3==0){ ## calculating the log-posteriors based on sample
      cur.log.posterior.beta2<- tilde.log.post.beta2(Y=Y, log.X2 = cur.log.X2, log.tilde.param.beta2=cur.samples.beta2.zero, batch=ntime.total)
      prop.log.posterior.beta2 <- tilde.log.post.beta2(Y=Y, log.X2 = cur.log.X2, log.tilde.param.beta2=prop.log.beta2, batch=ntime.total)
      
      if(Prop_type_3=="MALA"){
        sum.prop.log.proposal.mala.beta2.batch <-  log.proposal.density.mala.beta2(Y=Y, log.X2 = cur.log.X2, log.tilde.param.beta2=cur.samples.beta2.zero, log.tilde.param.beta2.star=prop.log.beta2 ,sigma2.beta2=sigma2.beta2, batch=ntime.total)
        sum.cur.log.proposal.mala.beta2.batch <- log.proposal.density.mala.beta2(Y=Y, log.X2 = cur.log.X2, log.tilde.param.beta2=prop.log.beta2, log.tilde.param.beta2.star=cur.samples.beta2.zero ,sigma2.beta2=sigma2.beta2, batch=ntime.total)
        log.ratio.beta2<- prop.log.posterior.beta2 + sum.cur.log.proposal.mala.beta2.batch - cur.log.posterior.beta2 - sum.prop.log.proposal.mala.beta2.batch
      } else if(Prop_type_3=="RWM"){
        log.ratio.beta2<- prop.log.posterior.beta2 - cur.log.posterior.beta2 
        
      }
      
      if (log(runif(1)) < log.ratio.beta2 & is.na(log.ratio.beta2)==FALSE){
        
        cur.samples.beta2.zero<-cur.samples.beta2
        rate.beta2<-rate.beta2+1
      } else{
        cur.samples.beta2<- cur.samples.beta2.zero
        
      }
    }
    sigma2.beta2<-adpative_function(index_MCMC_iter = i, sigma2_adapt = sigma2.beta2, target_accept = 0.57, rate_adapt = rate.beta2, burning1 = burning1,
                                    burning2 = burning2, adapt_seq = beta2_adapt_seq, adapt = adapt/T.acc_3, adpat_param = theta.beta2, lower.acc = 0.50, upper.acc = 0.65)
    
    
    ############  Updating beta3 and rho parameter ############
    if(Prop_type_4=="MALA"){
      prop.log.beta3.rho<- proposal.mala.beta3.rho(Y=Y[SG_index,], log.X3 = cur.log.X3[SG_index,], log.tilde.param.beta3 = cur.log.beta3, log.tilde.param.rho = cur.log.rho,
                                                   dist.mat = dist.mat, delta = delta, sigma2.beta3.rho = sigma2.beta3.rho,  batch = batch)
    } else if(Prop_type_4=="RWM"){
      prop.log.beta3.rho<-rnorm(n=2, mean=c(cur.log.beta3, cur.log.rho),  sd= sqrt(sigma2.beta3.rho))
      
    }
    
    cur.samples.beta3<- prop.log.beta3.rho[1]
    cur.samples.rho<- prop.log.beta3.rho[2]
    
    if(i%%T.acc_4==0){ ## calculating the log-posteriors based on sample
      cur.log.posterior.beta3.rho<- tilde.log.post.beta3.rho(Y=Y, log.X3 = cur.log.X3, log.tilde.param.beta3 = cur.samples.beta3.zero,
                                                             log.tilde.param.rho = cur.samples.rho.zero,  dist.mat = dist.mat, batch = ntime.total)
      prop.log.posterior.beta3.rho <- tilde.log.post.beta3.rho(Y=Y, log.X3 = cur.log.X3, log.tilde.param.beta3 = prop.log.beta3.rho[1],
                                                               log.tilde.param.rho = prop.log.beta3.rho[2],  dist.mat = dist.mat, batch = ntime.total)
      
      # cur.log.posterior.beta3.rho<-tilde.log.posterior(Y=Y, log.X2=cur.log.X2, log.X3 = cur.log.X3, log.tilde.param=c(cur.samples.scale.zero, cur.samples.beta1.zero, cur.samples.beta2.zero,
      #                                                                                                                 cur.samples.beta3.zero, cur.samples.rho.zero), dist.mat=dist.mat, cens.ind=cens.ind, u.thr=u.thr, cov=cov )
      # prop.log.posterior.beta3.rho<-tilde.log.posterior(Y=Y, log.X2=cur.log.X2, log.X3 = cur.log.X3, log.tilde.param=c(cur.samples.scale.zero, cur.samples.beta1.zero, cur.samples.beta2.zero,
      #                                                                                                                  prop.log.beta3.rho[1],  prop.log.beta3.rho[2]), dist.mat=dist.mat, cens.ind=cens.ind, u.thr=u.thr, cov=cov )
      if(Prop_type_4=="MALA"){
        sum.prop.log.proposal.mala.beta3.rho.batch <-  log.proposal.density.mala.beta3.rho(Y=Y, log.X3 = cur.log.X3, log.tilde.param.beta3 = cur.samples.beta3.zero,
                                                                                           log.tilde.param.beta3.star=prop.log.beta3.rho[1], log.tilde.param.rho = cur.samples.rho.zero, log.tilde.param.rho.star=prop.log.beta3.rho[2],
                                                                                           dist.mat = dist.mat, delta = delta, sigma2.beta3.rho = sigma2.beta3.rho,  batch = ntime.total)
        sum.cur.log.proposal.mala.beta3.rho.batch <- log.proposal.density.mala.beta3.rho(Y=Y, log.X3 = cur.log.X3, log.tilde.param.beta3 = prop.log.beta3.rho[1],
                                                                                         log.tilde.param.beta3.star=cur.samples.beta3.zero, log.tilde.param.rho = prop.log.beta3.rho[2],log.tilde.param.rho.star=cur.samples.rho.zero,
                                                                                         dist.mat = dist.mat, delta = delta, sigma2.beta3.rho = sigma2.beta3.rho,  batch = ntime.total)
        
        log.ratio.beta3.rho<- prop.log.posterior.beta3.rho + sum.cur.log.proposal.mala.beta3.rho.batch - cur.log.posterior.beta3.rho - sum.prop.log.proposal.mala.beta3.rho.batch
      } else if(Prop_type_4=="RWM"){
        log.ratio.beta3.rho<- prop.log.posterior.beta3.rho - cur.log.posterior.beta3.rho 
        
      }
      
      
      if (log(runif(1)) < log.ratio.beta3.rho & is.na(log.ratio.beta3.rho)==FALSE){
        cur.samples.beta3.zero<-cur.samples.beta3
        cur.samples.rho.zero<-cur.samples.rho
        rate.beta3.rho<-rate.beta3.rho+1
      } else{
        cur.samples.beta3<- cur.samples.beta3.zero
        cur.samples.rho<-cur.samples.rho.zero
      }
    }
    sigma2.beta3.rho<-adpative_function(index_MCMC_iter = i, sigma2_adapt = sigma2.beta3.rho, target_accept = 0.57, rate_adapt = rate.beta3.rho, burning1 = burning1,
                                        burning2 = burning2, adapt_seq = beta3.rho_adapt_seq, adapt = adapt/T.acc_4, adpat_param = theta.beta3.rho, lower.acc = 0.50, upper.acc = 0.65)
    
    
    ############# Updating latent X2 #######################
    prop.log.X2 <- proposal.mala.X2(Y=Y[SG_index,], log.X2 = cur.log.X2[SG_index] ,log.X3=cur.log.X3[SG_index,],log.tilde.param.scale=cur.samples.scale.zero, log.tilde.param.beta1=cur.samples.beta1.zero,
                                    log.tilde.param.beta2=cur.samples.beta2.zero, cens.ind=cens.ind[SG_index,], u.thr=u.thr[SG_index,], cov = cov, sigma2.latentX2=sigma2.latentX2) # Proposed latent parameters X2(s) using the MALA
    
    cur.log.posterior.latentX2<- tilde.log.post.X2(Y=Y[SG_index,],log.X2=cur.log.X2[SG_index],log.X3=cur.log.X3[SG_index,], log.tilde.param.scale=cur.samples.scale.zero, log.tilde.param.beta1=cur.samples.beta1.zero,
                                                   log.tilde.param.beta2=cur.samples.beta2.zero, cens.ind=cens.ind[SG_index,], u.thr=u.thr[SG_index,], cov = cov)## log posterior at proposed values
    
    prop.log.posterior.latentX2 <- tilde.log.post.X2(Y=Y[SG_index,],log.X2=prop.log.X2,log.X3=cur.log.X3[SG_index,],log.tilde.param.scale=cur.samples.scale.zero, log.tilde.param.beta1=cur.samples.beta1.zero,
                                                     log.tilde.param.beta2=cur.samples.beta2.zero, cens.ind=cens.ind[SG_index,], u.thr=u.thr[SG_index,], cov = cov) ## log posterior at current values
    
    prop.log.proposal.mala.latentX2<- log.proposal.density.mala.X2(Y=Y[SG_index,],log.X2=cur.log.X2[SG_index],log.X2.star=prop.log.X2, log.X3=cur.log.X3[SG_index,],log.tilde.param.scale=cur.samples.scale.zero, log.tilde.param.beta1=cur.samples.beta1.zero,
                                                                   log.tilde.param.beta2=cur.samples.beta2.zero, cens.ind=cens.ind[SG_index,], u.thr=u.thr[SG_index,], cov = cov, sigma2.latentX2=sigma2.latentX2) #
    
    cur.log.proposal.mala.latentX2<- log.proposal.density.mala.X2(Y=Y[SG_index,],log.X2=prop.log.X2,log.X2.star=cur.log.X2[SG_index],log.X3=cur.log.X3[SG_index,],log.tilde.param.scale=cur.samples.scale.zero, log.tilde.param.beta1=cur.samples.beta1.zero,
                                                                  log.tilde.param.beta2=cur.samples.beta2.zero, cens.ind=cens.ind[SG_index,], u.thr=u.thr[SG_index,], cov = cov, sigma2.latentX2=sigma2.latentX2) #
    log.ratio.latentX2<- prop.log.posterior.latentX2 + cur.log.proposal.mala.latentX2 - cur.log.posterior.latentX2 - prop.log.proposal.mala.latentX2
    if (log(runif(1)) < log.ratio.latentX2 & is.na(log.ratio.latentX2)==FALSE){
      cur.log.X2[SG_index]<-prop.log.X2
      cur.samples.latentX2<- c(cur.log.X2)
      rate.latentX2<-rate.latentX2+1
    } else{
      cur.samples.latentX2<- cur.log.X2
    }
    
    sigma2.latentX2<-adpative_function(index_MCMC_iter = i, sigma2_adapt = sigma2.latentX2, target_accept = 0.57, rate_adapt = rate.latentX2, burning1 = burning1,
                                       burning2 = burning2, adapt_seq = X2_adapt_seq, adapt = adapt, adpat_param = theta.latentX2, lower.acc = 0.50, upper.acc = 0.65)
    
    ############# Updating latent X3 #######################
    ###################### Proposing X3 ##########################################
    prop.log.X3<- proposal.mala.X3(Y=Y[SG_index,], log.X2=cur.samples.latentX2[SG_index], log.X3=cur.log.X3[SG_index,], log.tilde.param.scale=cur.samples.scale.zero, log.tilde.param.beta1=cur.samples.beta1.zero, 
                                   log.tilde.param.beta3=cur.samples.beta3.zero, log.tilde.param.rho=cur.samples.rho.zero, dist.mat=dist.mat, cens.ind=cens.ind[SG_index,], u.thr=u.thr[SG_index,], cov = cov, sigma2.latentX3=sigma2.latentX3) # Proposed latent parameters X2(s) using the MALA
    
    prop.log.proposal.mala.latentX3<- log.proposal.density.mala.X3(Y=Y[SG_index,],log.X2=cur.samples.latentX2[SG_index],log.X3=cur.log.X3[SG_index,],log.X3.star=prop.log.X3, log.tilde.param.scale=cur.samples.scale.zero, log.tilde.param.beta1=cur.samples.beta1.zero, 
                                                                   log.tilde.param.beta3=cur.samples.beta3.zero, log.tilde.param.rho=cur.samples.rho.zero, dist.mat=dist.mat, cens.ind=cens.ind[SG_index,], u.thr=u.thr[SG_index,], cov = cov, sigma2.latentX3=sigma2.latentX3)
    
    cur.log.proposal.mala.latentX3 <- log.proposal.density.mala.X3(Y=Y[SG_index,],log.X2=cur.samples.latentX2[SG_index],log.X3=prop.log.X3,log.X3.star=cur.log.X3[SG_index,], log.tilde.param.scale=cur.samples.scale.zero, log.tilde.param.beta1=cur.samples.beta1.zero, 
                                                                   log.tilde.param.beta3=cur.samples.beta3.zero, log.tilde.param.rho=cur.samples.rho.zero, dist.mat=dist.mat, cens.ind=cens.ind[SG_index,], u.thr=u.thr[SG_index,], cov = cov, sigma2.latentX3=sigma2.latentX3)
    
    cur.log.posterior.latentX3<-  tilde.log.post.X3(Y=Y[SG_index,],log.X2=cur.samples.latentX2[SG_index],log.X3=cur.log.X3[SG_index,], log.tilde.param.scale=cur.samples.scale.zero, log.tilde.param.beta1=cur.samples.beta1.zero, 
                                                    log.tilde.param.beta3=cur.samples.beta3.zero, log.tilde.param.rho=cur.samples.rho.zero, dist.mat=dist.mat, cens.ind=cens.ind[SG_index,], u.thr=u.thr[SG_index,], cov = cov)
    
    prop.log.posterior.latentX3 <- tilde.log.post.X3(Y=Y[SG_index,],log.X2=cur.samples.latentX2[SG_index],log.X3=prop.log.X3, log.tilde.param.scale=cur.samples.scale.zero, log.tilde.param.beta1=cur.samples.beta1.zero, 
                                                     log.tilde.param.beta3=cur.samples.beta3.zero, log.tilde.param.rho=cur.samples.rho.zero, dist.mat=dist.mat, cens.ind=cens.ind[SG_index,], u.thr=u.thr[SG_index,], cov = cov)
    
    log.ratio.latentX3<- prop.log.posterior.latentX3 + cur.log.proposal.mala.latentX3 - cur.log.posterior.latentX3 - prop.log.proposal.mala.latentX3
    if (log(runif(1)) < log.ratio.latentX3 & is.na(log.ratio.latentX3)==FALSE){
      cur.log.X3[SG_index,]<-prop.log.X3
      cur.samples.latentX3<- c(cur.log.X3)
      rate.latentX3<-rate.latentX3+1
    } else{
      cur.samples.latentX3<- cur.log.X3
    }
    
    sigma2.latentX3<-adpative_function(index_MCMC_iter = i, sigma2_adapt = sigma2.latentX3, target_accept = 0.57, rate_adapt = rate.latentX3, burning1 = burning1,
                                       burning2 = burning2, adapt_seq = X3_adapt_seq, adapt = adapt, adpat_param = theta.latentX3, lower.acc = 0.50, upper.acc = 0.65)
    
    
    ### saving the samples for some all the hyperparameters and some selected latent parameters
    
    if((i%%thin)-1==0){
      samples[j,]<- c(cur.samples.scale.zero, cur.samples.beta1.zero, cur.samples.beta2.zero, cur.samples.beta3.zero, cur.samples.rho.zero, 
                      cur.samples.latentX2[1:4], cur.samples.latentX3[1:4])
      j=j+1
    }
    #print(i)
    
    ### saving the chains of tuning parameters 
    if((i%%adapt)-1==0){ # to save all exp the scale parameter of the MALA
      sigma.matrix[m,]<- c(sigma2.scale, sigma2.beta1, sigma2.beta2, sigma2.beta3.rho, sigma2.latentX2, sigma2.latentX3)
      m=m+1
    }
  }
  return(list("samples"=samples,"tuning_param"=sigma.matrix,"Acc.rate_hyper"=rate.scale/(N.MCMC-(burning1+burning2)),"Acc.rate_latentX2"=rate.latentX2/(N.MCMC-(burning1+burning2)),"Acc.rate.latentX3"=rate.latentX3/(N.MCMC-(burning1+burning2))))
}