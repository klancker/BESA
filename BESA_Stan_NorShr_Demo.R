#-------------------------------------------------------------------------------
## BESA IMPLEMENTATION FOR BARENTS SEA SHRIMPS
## 05.08.2022
## authors: Kira Lancker, Martin F. Quaas, Rudi Voss, Fabian Zimmermann
#-------------------------------------------------------------------------------
# clear workspace, and initial settings
 rm(list=ls(all=TRUE)) 
 graphics.off()
 gc()   
 set.seed(5) # use for comparing results between runs
 options(digits=3) # display numbers with three digits
 #options(mc.cores=parallel::detectCores())# To run in parallel on multiple cores
#-------------------------------------------------------------------------------
# install and load packages
# library(rstan)
# library(bayesplot)
# library(HDInterval)
# library("LearnBayes")
 for (package in c('rstan', 'bayesplot','HDInterval')) {
   if (!require(package, character.only=T, quietly=T)) {
     install.packages(package)
     library(package, character.only=T)
   }
 }

#-------------------------------------------------------------------------------
# read in codes needed
 source("CMSY_fun.R")
 mod_BESA <- stan_model('BESAstan.stan')
 rstan_options(auto_write = TRUE)
 # in case the line above does not compile, the problem might be related to a current compatibility issue
 # between the rstan package and R version 4.2, please visit:  
 # https://discourse.mc-stan.org/t/rstan-error-in-compilecode-f-code-language-language-verbose-verbose/28007
#-------------------------------------------------------------------------------
# define parameters for Stan runs:  
 nsamples<-20000   # number of samples to run per chain 
 nchains<-8      # number of chains to run
 pthin<-5   # thinning parameter
 pwarmup<-6000 # number of samples to include in warmup phase
 ncores<-8   # number of cores on computer to use for Stan runs
#-------------------------------------------------------------------------------
# read in data: 
 DATA<-read.csv("Data.csv")
 DATA_H<-DATA$Catch
 DATA_P<-DATA$Price_CPI
 Nt<-length(DATA_P)
 s<-1  # which year to start in with estimation? 
 Ht<-DATA_H[s:Nt]/1000  # optionally rescale or shorten time series
 pt<-log(DATA_P[s:Nt]/1000) # optionally rescale or shorten time series
 TT<-length(Ht)
#-------------------------------------------------------------------------------
# initialize: 
 result.depletion<-c()
 result.biomass<-c()
 hdi.depletionl<-c()
 hdi.depletionh<-c()
 hdi.biomassl<-c()
 hdi.biomassh<-c()
 BMEY=c()
 MEY=c()
#-------------------------------------------------------------------------------
# obtain priors for r, K from CMSY using the following settings:  
  res=2   # choose resilience class
# initial range of r based on resilience
  if(res == 4) {
    start_r <- c(0.6,1.5)} else if(res == 3) {
      start_r <- c(0.2,0.8)}    else if(res == 2) {
        start_r <- c(0.05,0.5)}  else { # i.e. res== 1
          start_r <- c(0.015,0.1)} 
  r.low<-start_r[1]
  r.hi<-start_r[2]
  
  output_list<-CMSY_fun_real(ct.raw=Ht,
                             r.low,r.hi,
                             nyr=TT,
                             int.yr.i=round(TT/2, digits = 0),
                             int.yr=round(TT/2, digits = 0),
                             end.yr=TT,
                             stb.low=0.5,stb.hi=0.95,
                             intb.low = 0.5,intb.hi = 0.95,
                             endb.low= 0.5,endb.hi=0.95)
 

# To connect the code parts, the CMSY code was slightly adapted 
 # from its original description in: 
 # Froese, R., Demirel, N., Coro, G., Kleisner, K. M., & Winker, H. (2016). 
 # Estimating fisheries reference points from catch and resilience. Fish and Fisheries. 
 # http://onlinelibrary.wiley.com/doi/10.1111/faf.12190/full. 
 # Original Code is available on http://oceanrep.geomar.de/33076/.
#-------------------------------------------------------------------------------
# define prior distribution parameters, fix parameters and multipliers for BESA: 
  mur=as.numeric(output_list[5])  # param. 1 for logn. prior for r from CMSY
  sigr=as.numeric(output_list[2]) # param. 2 for logn. prior for r from CMSY
  musigeps=0.5 # parameter 1 for normal uninformative prior sigma_epsilon
  sigsigeps=1  # parameter 2 for normal uninformative prior sigma_epsilon
  musigeta=0.5 # parameter 1 for normal uninformative prior sigma_eta
  sigsigeta=1  # parameter 2 for normal uninformative prior sigma_eta
  sigK=as.numeric(output_list[4]) # parameter 2 for logn. prior for K from CMSY
  sigc=1      # parameter 2 for lognormal prior for c
  F11=9.03#4.56    # parameter 1 for beta prior for starting depletion exp(x1)
  F12=2.74#4.56    # parameter 1 for beta prior for starting depletion exp(x1)
  m=as.numeric(output_list[3]) # multiplier for K: centers prior for K, from CMSY 
  mc=exp(pt[1]+(log(0.5)+log(m)))  # multiplier for c: based on priors x1 & K
  #quantile1=list(p=.025, x=0.5)     # 2.5% quantile should be 0.2
  #quantile2=list(p=.975, x=0.95)      # 97.5% quantile should be 0.8
  #beta.select(quantile1, quantile2)
#-------------------------------------------------------------------------------  
# run Bayesian estimation via Hamiltonian Monte Carlo MCMC (Stan)
  # input: 
  dataStan <- list(pt=pt,
                   Ht=Ht, 
                   TT=TT,
                   m=m,
                   mc=mc,
                   mur=mur,
                   sigr=sigr,
                   sigK=sigK,
                   sigc=sigc,
                   musigeps=musigeps,
                   sigsigeps=sigsigeps,
                   musigeta=musigeta,
                   sigsigeta=sigsigeta,
                   F11=F11,  
                   F12=F12)#, 
                   #sigd=1,
                   #muK=1)
 x <- sample(1:10000, 1) # draw random seed (can omit by deleting seed below)
 # run Stan
 BESAStan <- sampling(mod_BESA, 
                       data = dataStan, 
                       chains = nchains,
                       seed=x, 
                       iter=nsamples,
                       warmup=pwarmup,
                       thin=pthin, 
                       cores=ncores)
#-------------------------------------------------------------------------------
  # close down the Stan model and parallel implementation: 
  rm(mod_BESA)
  gc()
#-------------------------------------------------------------------------------  
# extract output: 
  fit_res<- extract(BESAStan, permuted = TRUE, inc_warmup = FALSE)
  print(BESAStan, max = 150) # for an overview 
  # extract the posterior distribution vectors for quantities of interest
  chat <- fit_res$truec
  rhat <- fit_res$truer
  Khat <- fit_res$truek
  sigepshat <- fit_res$sigepstrue
  sigetahat <- fit_res$sigetatrue
  MSYhat <- fit_res$MSY
  #lp_<- fit_res$lp_
  Xshat<- fit_res$Xs
  Z1hat<-fit_res$Z1
  Xt<-exp(fit_res$xt)
  Bt<-fit_res$Bt
  #diag_rhat<- fit_res$rhat
  #xtmeans=colMeans(fit_res$xt)
  
  # compute point estimates
  cmean<-exp(mean(log(chat))+0.5*sd(log(chat))^2)
  rmean<-exp(mean(log(rhat))+0.5*sd(log(rhat))^2)
  Kmean<-exp(mean(log(Khat))+0.5*sd(log(Khat))^2)
  sigepsmean<-mean(sigepshat)
  sigetamean<-mean(sigetahat)  
  MSYmean<-exp(mean(log(MSYhat))+0.5*sd(log(MSYhat))^2)
  Xsmean<-exp(mean(log(Xshat))+0.5*sd(log(Xshat))^2)
  Z1mean<-exp(mean(log(Z1hat))+0.5*sd(log(Z1hat))^2)
  for (t in 1:TT){
    result.depletion[t] <- exp(mean(log(Xt[,t])+0.5*sd(log(Xt[,t]))^2))
    result.biomass[t] <- exp(mean(log(Bt[,t])+0.5*sd(log(Bt[,t]))^2))
    HDIdepl=hdi(Xt[,t], credMass = 0.95)
    hdi.depletionl[t]=HDIdepl[1]
    hdi.depletionh[t]=HDIdepl[2]
    HDIbio=hdi(Bt[,t], credMass = 0.95)
    hdi.biomassl[t]=HDIbio[1]
    hdi.biomassh[t]=HDIbio[2]    
  }
#-------------------------------------------------------------------------------
# compute MEY
  pbar=mean(exp(pt))
  for (j in 1:length(rhat)){
    BMEY[j]=0.5*(Khat[j]+chat[j]/pbar)
    MEY[j]=rhat[j]*BMEY[j]*(1-BMEY[j]/Khat[j])
  }
  BMEYmean=mean(BMEY)
  MEYmean=mean(MEY)
#-------------------------------------------------------------------------------
# diagnostics, pairs plot, traceplots
 #check HMC specific diagnostics: Divergences, Tree depth, and Energy Bayesian Fraction of Missing Information.
 rstan::check_hmc_diagnostics(BESAStan)
 # pairs plots
 pairs(BESAStan, pars = c("r", "K", "c","Xs"), las = 1) # plots red for divergences, yellow for hitting the max treedepth
 pairs(BESAStan, pars = c("c", "sigeps","sigeta"), las = 1)
 # Traceplots (Bayesplot, includes warm up and should ensure whether warmup is long enough)
 # outside of warmup period (dark grey area), chains should remain in the same neighborhood, and different chains (colors) should 
 # overlap closely
 traceplot(BESAStan, pars="r", inc_warmup = TRUE)
 traceplot(BESAStan, pars="K", inc_warmup = TRUE)
 traceplot(BESAStan, pars="c", inc_warmup = TRUE)
 traceplot(BESAStan, pars="sigeps", inc_warmup = TRUE)
 traceplot(BESAStan, pars="sigeta", inc_warmup = TRUE)
 traceplot(BESAStan, pars="xt[1]", inc_warmup = TRUE)
 traceplot(BESAStan, pars="Z1", inc_warmup = TRUE)
 # number of divergent samples: 
 params <- as.data.frame(extract(BESAStan, permuted=FALSE))
 divergent=0
 for(w in 1:nchains){
   divergent <- divergent+sum(get_sampler_params(BESAStan, inc_warmup=FALSE)[[w]][,'divergent__'])
 }
 divergentsams<-divergent
 # Rhat
 fit_summary <- summary(BESAStan)
 vec_rhat<-as.numeric(fit_summary$summary[,10])
#-------------------------------------------------------------------------------
# plot posterior distributions
 dev.new()  # optionally open new window
 par(mfrow=c(4,2))  # figure with multiple tabs
 # r
 plot(density(rhat),main="",xlab="r",ylab="Density", xlim = c(0,2), ylim = c(0,7), col="blue")
 y<-rlnorm(100000,mean=log(as.numeric(output_list[5])),sd=as.numeric(output_list[2]))
 lines(density(y/4), ylim = c(0,7), col="red", lty = 2)
 # K
 plot(density(Khat),xlab="K",ylab="Density",main="", xlim = c(0,Kmean*2), col="blue")
 y<-rlnorm(100000,mean=log(1),sd=sigK)
 lines(density(y*m),main="", col="red", lty = 2)
 # c
 plot(density(chat),xlab="c",ylab="Density",main="", col="blue")
 x<-seq(from=0,to=200,length.out=300)
 lines(x,dlnorm(x,mean=log(mc),sd=sigc),main="", col="red", lty = 2)
 # Starting depletion x1 
 plot(density(Xshat),xlab="starting depletion, chi_1",ylab="Density",main="", col="blue")
 x<-seq(from=0,to=1,length.out=300)
 lines(x,dbeta(x,shape1=F11,shape2=F12),main="", col="red", lty = 2)
 # sigma epsilon
 plot(density(sigepshat),xlab="Observation std.err., sigma_epsilon",ylab="Density",main="", col="blue",xlim=c(0, 0.1))
 x<-seq(from=0,to=0.4,length.out=300)
 lines(x,dnorm(x,mean=musigeps,sd=sigsigeps),main="", col="red", lty = 2)
 # sigma eta
 plot(density(sigetahat),xlab="Process std.err., sigma_eta",ylab="Density",main="", col="blue")
 x<-seq(from=0,to=0.4,length.out=300)
 lines(x,dnorm(x,mean=musigeta,sd=sigsigeta),main="", col="red", lty = 2)
 # variance of starting depletion, Z1
 plot(density(Z1hat),xlab="Z_1",ylab="Density",main="", col="blue")
 x<-seq(from=0,to=0.4,length.out=300)
 lines(x,dnorm(x,mean=0.01,sd=1),main="", col="red", lty = 2)
 # legend
 plot(1, type = "n", axes=FALSE, xlab="", ylab="")
 legend(x="top", legend=c("prior", "posterior"),col=c("red", "blue"),inset=0, lty=c(2,1), cex=1.4, horiz=TRUE)

#-------------------------------------------------------------------------------
# plot time series outcomes:
 TIMES=DATA$Year[(length(DATA$Year)-Nt+1):length(DATA$Year)]

# depletion and biomass time series 
 par(mfrow=c(2,1),mai = c(1, 0.8, 0.1, 0.1))
 plot(TIMES,result.depletion[1:(45-s+1)],col="blue", cex.lab = 1.2, ty="o",ylim=c(0.3,1.27), lty = 3,xlab="time",ylab="depletion")
 polygon(c(rev(TIMES), TIMES), c(rev(hdi.depletionl), hdi.depletionh), col = 'grey', border = NA)
 lines(TIMES,result.depletion[1:(45-s+1)],col="blue", cex.lab = 1.2,lwd=2, ty="o",ylim=c(0.3,1.27), lty = 3,xlab="time",ylab="depletion")
 lines(TIMES,hdi.depletionl,lty = 'dashed')
 lines(TIMES,hdi.depletionh,lty = 'dashed')
 lines(TIMES, DATA$BBmsy[s:45]/2, col="red",lwd=2)
 legend("topright", legend=c("BESA", "NIPAG"),col=c("blue", "red"), pch = c(1, NA), lty=c(3,1), cex=0.8)
 plot(TIMES,result.biomass[1:Nt],col="blue", cex.lab = 1.2, ty="o",ylim=c(0, Kmean*2), lty = 3,xlab="time",ylab="biomass (tons)")
 polygon(c(rev(TIMES), TIMES), c(rev(hdi.biomassl), hdi.biomassh), col = 'grey', border = NA)
 lines(TIMES,result.biomass[1:Nt],col="blue", cex.lab = 1.2,lwd=2, ty="o",ylim=c(0, Kmean*2), lty = 3,xlab="time",ylab="depletion")
 lines(TIMES,hdi.biomassl,lty = 'dashed')
 lines(TIMES,hdi.biomassh,lty = 'dashed')
 lines(TIMES,  DATA$BBmsy[s:45]*DATA$K_Spict[s:45]/2, col="red",lwd=2)
 legend("topright", legend=c("BESA", "NIPAG"),col=c("blue", "red"), pch = c(1, NA), lty=c(3,1), cex=0.8)

#-------------------------------------------------------------------------------
# plot prices and harvest for paper
 par(mfrow=c(1,1))
 plot(TIMES, exp(pt[1:Nt])*1000,pch=3, ty="o", lty = 3, cex = 0.8,col="red", xlab = "year", ylim=c(10,35), ylab = "Price (NOK/ton, defl., red)",las = 1)
 par(new = TRUE) # second yaxis
 plot(TIMES, Ht,pch=3, cex = 0.8, ty="o", lty = 3, col="blue", xaxt = "n", yaxt = "n", ylab = "", xlab = "")
 axis(side = 4)
 mtext("catch (tons, blue)", side = 4)

#------------------------------------------------------------------------------- 
# price fit graph
 pred_Pt=cmean*(result.depletion[1:Nt]*Kmean)^(-d) 
 par(mfrow=c(1,1))
 plot(1:Nt, exp(pt[1:Nt])*1000,pch=3, cex = 0.8,col="red", ty="o", ylim=c(10,35), lty = 3,xlab = "year", ylab = "Price (NOK/ton, defl.)",xlim = c(0,Nt),las = 1)
 lines(1:Nt, pred_Pt[1:Nt]*1000, cex = 0.8, col="blue")
 legend("top", legend=c("observed","BESA predicted"),col=c("blue", "red"), lty=1, cex=0.8)
#------------------------------------------------------------------------------- 
# optionally export time series results:
 #df <- data.frame(result.depletion,hdi.depletionl,hdi.depletionh,result.biomass,hdi.biomassl,hdi.biomassh)
 #write.csv(df,"Results_Shrimps.csv", row.names = FALSE)
#-------------------------------------------------------------------------------
# highest density Bayesian 95\% credible intervals 
 hdi(fit_res$truec, credMass = 0.95)
 hdi(fit_res$truer, credMass = 0.95)
 hdi(fit_res$truek, credMass = 0.95)
 hdi(fit_res$sigepstrue, credMass = 0.95)
 hdi(fit_res$sigetatrue, credMass = 0.95)
 hdi(fit_res$MSY, credMass = 0.95)
 hdi(BMEY, credMass = 0.95)
 hdi(MEY, credMass = 0.95)
 
 
