# To connect the code parts, the CMSY code was slightly adapted 
# from its original description in: 
# Froese, R., Demirel, N., Coro, G., Kleisner, K. M., & Winker, H. (2016). 
# Estimating fisheries reference points from catch and resilience. Fish and Fisheries. 
# http://onlinelibrary.wiley.com/doi/10.1111/faf.12190/full. 
# Original Code is available on http://oceanrep.geomar.de/33076/.
##---------------------------------------------------------------------------------------------
## CMSY analysis with estimation of total biomass, including Bayesian Schaefer with recruitment 
## Written by Rainer Froese, Gianpaolo Coro and Henning Winker
##---------------------------------------------------------------------------------------------
CMSY_fun_real <- function(ct.raw,r.low,r.hi,
                     int.yr.i=25,
                     nyr=50,
                     start.yr=1,
                     end.yr=50,
                     stb.low=0.2,
                     stb.hi=0.8,
                     int.yr=25,
                     intb.low=0.2,
                     intb.hi=0.8, 
                     endb.low= 0.2,
                     endb.hi=0.6) {

#library(coda) 
library(EnvStats)

#-----------------------------------------
# General settings for the analysis
#-----------------------------------------
dataUncert   <- 0.1  # set observation error as uncertainty in catch - default is SD=0.1
sigmaR       <- 0.1 # overall process error for CMSY; SD=0.1 is the default
n            <- 10000 # initial number of r-k pairs
ni           <- 3 # iterations for r-k-startbiomass combinations, to test different variability patterns; no improvement seen above 3
mgraphs      <- F # set to TRUE to produce additional graphs for management
write.output <- T # set to true if table with results in output file is wanted
qyr          <- c(NA,NA) # c(1985,2000) #  set year range with stable catch/biomass ratio, min 5 years; else set NA
nab          <- 10 # default=10; minimum number of years with abundance data for full Schaefer model
u.CMSY.25p   <- NA # initialize variable; used to print 25th percentile in last year

#----------------------------------------------
# Monte Carlo filtering with Schaefer Function
#----------------------------------------------
SchaeferMC <- function(ri, ki, startbio, int.yr, intbio, endbio, sigR, pt, duncert, startbins, ni) {
   
  # create vector for initial biomasses
  startbt     <- seq(from =startbio[1], to=startbio[2], by = (startbio[2]-startbio[1])/startbins)
  nstartbt    <- length(startbt)
  npoints     <- length(ri)
  
  # create vectors for viable r, k and bt
  mdat        <- matrix(data=NA, nrow = npoints*nstartbt, ncol = 2+nyr+1) 
  
  
  #loop through r-k pairs
  for(i in 1 : npoints) {
    if (i%%1000==0)
    
    # create empty vector for annual biomasses
    bt <- vector()
    # set flag for viable point to FALSE
    vp <- F
    
    # loop through range of relative start biomasses
    for(j in startbt) {    
      # set initial biomass, including 0.1 process error to stay within bounds
      bt[1]=j*ki[i]*exp(rnorm(1,0, 0.1*sigR))  ## set biomass in first year
      
      # repeat test of r-k-startbt combination to allow for different random error
      for(re in 1:ni) {
        
        #loop through years in catch time series
        for(t in 1:nyr) {  # for all years in the time series
          xt=rnorm(1,0, sigR) # set new process error for every year  
          zlog.sd = sqrt(log(1+(duncert)^2))
          zt=rlnorm(1,meanlog = 0, sdlog = zlog.sd) # model the catch error as a log normal distribution.
          
          bt[t+1] <- bt[t]+ri[i]*bt[t]*(1-bt[t]/ki[i])*exp(xt)-ct[t]*zt
          
          # if biomass < 0.01 k, discard r-k-startbt combination
          if(bt[t+1] < 0.01*ki[i]) { break } # stop looping through years, go to next upper level
          
          # intermediate year check
          if ((t+1)==int.yr.i && (bt[t+1]>(intbio[2]*ki[i]) || bt[t+1]<(intbio[1]*ki[i]))) { break }  
          
        } # end of loop of years
        
        # if loop was broken or last biomass falls outside of expected ranges 
        # do not store results, go directly to next startbt
        if(t < nyr || bt[yr==end.yr] > (endbio[2]*ki[i]) || bt[yr==end.yr] < (endbio[1]*ki[i]) ) { next } else { 
          # store r, k, and bt, plot point, then go to next startbt
          mdat[((i-1)*nstartbt)+which(startbt==j),1]           <- ri[i]
          mdat[((i-1)*nstartbt)+which(startbt==j),2]           <- ki[i]
          mdat[((i-1)*nstartbt)+which(startbt==j),3:(2+nyr+1)] <- bt[1:(nyr+1)]/ki[i]
          vp <- T }
        
      } # end of repetition for random error
      
    } # end of j-loop of initial biomasses 

  } # end of i-loop of r-k pairs 
  
  # return results
  mdat <- na.omit(mdat) # ignore rows filled with NA from initial matrix   
  cat("\n")
  return(list(mdat))
} # end of SchaeferMC function

#---------------------------------
# Analyze stock(s)
#---------------------------------
stock=1

  # set global defaults for uncertainty
  duncert      <- dataUncert
  sigR         <- sigmaR

  # extract data on stock
  yr           <- c(1:50)
  bt <- NA


# change catch to 3 years moving average 
# for first year use reported catch and for second year 2 year average
ma                <- function(x){filter(x,rep(1/3,3),sides=1)}
  ct              <- ma(ct.raw)
  not.na          <- which(is.na(ct.raw)==F)
  ct[not.na[1]]   <- ct.raw[not.na[1]]
  ct[not.na[2]]   <- (ct.raw[not.na[1]]+ct.raw[not.na[2]])/2
  
  # initialize vectors for viable r, k, bt, and all in a matrix
  mdat.all    <- matrix(data=vector(),ncol=2+nyr+1)
  
  # initialize other vectors anew for each stock
  current.attempts <- NA
   
  #----------------------------------------------------
  # Determine initial ranges for parameters and biomass
  #----------------------------------------------------
  # initial range of r from input file

    start_r <- c(r.low,r.hi)

  
  # initial biomass range from input file
    startbio <- c(stb.low,stb.hi)
  
  # median of 5 highest catches = 3d highest catch
  m.max.ct     <- sort(ct.raw,decreasing=T)[3] 
  # index of years with lowest and highest catch
  min.yr     <- yr[which.min(ct.raw)]
  max.yr     <- yr[which.max(ct.raw)]
  
  # use year and biomass range for intermediate biomass from input file
  if(is.na(intb.low)==F & is.na(intb.hi)==F) {
    int.yr   <- int.yr
    intbio   <- c(intb.low,intb.hi)
  
  # if contrast in catch is low, use uninformative range 0-1 in mid-year
  } else if((m.max.ct-min(ct.raw))/m.max.ct < 0.5) {
    int.yr    <- as.integer(mean(c(start.yr, end.yr)))
    intbio    <- c(0,1) 
  
  # else if year of minimum catch is at least 3 years away from start.yr and end.yr of series, use min catch
  } else if((min.yr - start.yr) > 3 & (end.yr - min.yr) > 3 ) {
    # assume that biomass range in year before minimum catch was 0.01 - 0.4
    int.yr    <- min.yr-1
    intbio    <- c(0.01,0.4) 
  
  # else if year of max catch is at least 3 years away from start.yr and end.yr of series, use max catch  
  } else if((max.yr - start.yr) > 3 & (end.yr - max.yr) > 3 ) {
    # assume that biomass range in year before maximum catch was 0.2 - 0.6
    int.yr    <- max.yr-1
    intbio    <- c(0.2,0.6)  
  
  # else assume uninformative range 0-1 in mid-year 
  } else {
    int.yr     <- as.integer(mean(c(start.yr, end.yr)))
    intbio    <- c(0,1) }
  # end of intbio setting
  
  # final biomass range from input file
  if(is.na(endb.low)==F & is.na(endb.hi)==F) {
    endbio   <- c(endb.low,endb.hi)
  } else {
    # else use mean final catch/median max catch to estimate final biomass
    endbio  <- if(ct[nyr]/m.max.ct > 0.7) {c(0.5,0.9)} else if(ct[nyr]/m.max.ct < 0.3) {c(0.01,0.4)} else {c(0.2,0.6)}
  } # end of final biomass setting

  # initial prior range of k values, assuming min k will be larger than max catch / prior for r 
  if(mean(endbio) <= 0.5) {
     start_k <- c(m.max.ct/start_r[2],4*m.max.ct/start_r[1])} else {
      start_k <- c(2*m.max.ct/start_r[2],12*m.max.ct/start_r[1])} 
     cat("startbio=",startbio,", intbio=",int.yr,intbio,", endbio=",endbio,"\n")
    
  #------------------------------------------------------------------
  # Uniform sampling of the r-k space
  #------------------------------------------------------------------
  # get random set of r and k from log space distribution 
  ri1 = exp(runif(n, log(start_r[1]), log(start_r[2])))  
  ki1 = exp(runif(n, log(start_k[1]), log(start_k[2])))  

  
  #---------------------------------------------------------------------
  # 1 - Call CMSY-SchaeferMC function to preliminary explore the r-k space
  #---------------------------------------------------------------------
  cat("First Monte Carlo filtering of r-k space with ",n," points\n")
  MCA <-  SchaeferMC(ri=ri1, ki=ki1, startbio=startbio, int.yr=int.yr, intbio=intbio, endbio=endbio, sigR=sigR, 
                     pt=T, duncert=dataUncert, startbins=10, ni=ni)
   mdat.all <- rbind(mdat.all,MCA[[1]])
   rv.all   <- mdat.all[,1]
   kv.all   <- mdat.all[,2]
   btv.all  <- mdat.all[,3:(2+nyr+1)]
   # count viable trajectories and r-k pairs 
   n.viable.b   <- length(mdat.all[,1])
   n.viable.pt <- length(unique(mdat.all[,1]))
   cat("Found ",n.viable.b," viable trajectories for", n.viable.pt," r-k pairs\n")
 
  #----------------------------------------------------------------------- 
  # 2 - if the lower bound of k is too high, reduce it by half and rerun
  #-----------------------------------------------------------------------
  if(length(kv.all[kv.all < 1.1*start_k[1] & rv.all < mean(start_r)]) > 10) {
   cat("Reducing lower bound of k, resampling area with",n,"additional points\n")
   start_k <- c(0.5*start_k[1],start_k[2])
   ri1 = exp(runif(n, log(start_r[1]), log(start_r[2])))  
   ki1 = exp(runif(n, log(start_k[1]), log(start_k[2])))  
   MCA <-  SchaeferMC(ri=ri1, ki=ki1, startbio=startbio, int.yr=int.yr, intbio=intbio, endbio=endbio, sigR=sigR, 
                      pt=T, duncert=dataUncert, startbins=10, ni=ni)
   mdat.all <- rbind(mdat.all,MCA[[1]])
   rv.all   <- mdat.all[,1]
   kv.all   <- mdat.all[,2]
   btv.all  <- mdat.all[,3:(2+nyr+1)]
   n.viable.b   <- length(mdat.all[,1])
   n.viable.pt <- length(unique(mdat.all[,1]))
   cat("Found altogether",n.viable.b," viable trajectories for", n.viable.pt," r-k pairs\n")
 }

  #-------------------------------------------------------------------
  # 3 - if few points were found then resample and shrink the log k space
  #-------------------------------------------------------------------
  if (n.viable.b <= 1000){
    log.start_k.new  <- log(start_k) 
    max.attempts     <- 3
    current.attempts <- 1
    startbins        <- 10  
    while (n.viable.b <= 1000 && current.attempts <= max.attempts){
      if(n.viable.pt > 0) {
      log.start_k.new[1] <- mean(c(log(start_k[1]), min(log(kv.all))))
      log.start_k.new[2] <- mean(c(log.start_k.new[2], max(log(kv.all)))) }
      n.new <- n*current.attempts #add more points
      ri1 = exp(runif(n.new, log(start_r[1]), log(start_r[2])))  
      ki1 = exp(runif(n.new, log.start_k.new[1], log.start_k.new[2]))
      cat("Shrinking k space: repeating Monte Carlo in the interval [",exp(log.start_k.new[1]),",",exp(log.start_k.new[2]),"]\n")
      cat("Attempt ",current.attempts," of ",max.attempts," with ",n.new," additional points","\n")
      if(current.attempts==2 & n.viable.b < 50){
        duncert   <- 2*dataUncert
        sigR      <- 2*sigmaR
        startbins <- 20
        cat("Doubling startbins, catch and process error, and number of variability patterns \n")   
      }
      MCA <-  SchaeferMC(ri=ri1, ki=ki1, startbio=startbio, int.yr=int.yr, intbio=intbio, endbio=endbio, sigR=sigR, 
                         pt=T, duncert=duncert, startbins=startbins, ni=2*ni)
      mdat.all <- rbind(mdat.all,MCA[[1]])
      rv.all   <- mdat.all[,1]
      kv.all   <- mdat.all[,2]
      btv.all  <- mdat.all[,3:(2+nyr+1)]
      n.viable.b   <- length(mdat.all[,1])
      n.viable.pt <- length(unique(mdat.all[,1]))
      cat("Found altogether",n.viable.b," viable trajectories for", n.viable.pt," r-k pairs\n")
      current.attempts=current.attempts+1 #increment the number of attempts
    }
    if(n.viable.b < 5) {
      cat("Only",n.viable.pt,"viable r-k pairs found, check data and settings \n")
      next
    }  
  }
    
  #------------------------------------------------------------------
  # 4 - if tip of viable r-k pairs is 'thin', do extra sampling there
  #------------------------------------------------------------------
  if(length(rv.all[rv.all > 0.9*start_r[2]]) < 5) { 
    l.sample.r        <- quantile(rv.all,0.6)
    add.points        <- ifelse(is.na(current.attempts)==T,n,ifelse(current.attempts==2,2*n,ifelse(length(rv.all)>500,3*n,6*n)))
    cat("Final sampling in the tip area above r =",l.sample.r,"with",add.points,"additional points\n")
    log.start_k.new <- c(log(0.8*min(kv.all)),log(max(kv.all[rv.all > l.sample.r])))
  
    ri1 = exp(runif(add.points, log(l.sample.r), log(start_r[2])))  
    ki1 = exp(runif(add.points, log.start_k.new[1], log.start_k.new[2]))
    MCA <-  SchaeferMC(ri=ri1, ki=ki1, startbio=startbio, int.yr=int.yr, intbio=intbio, endbio=endbio, sigR=sigR, 
                     pt=T, duncert=duncert, startbins=10, ni=ni)
    mdat.all <- rbind(mdat.all,MCA[[1]])
    rv.all   <- mdat.all[,1]
    kv.all   <- mdat.all[,2]
    btv.all  <- mdat.all[,3:(2+nyr+1)]
    n.viable.b   <- length(mdat.all[,1])
    n.viable.pt <- length(unique(mdat.all[,1]))
   cat("Found altogether",n.viable.b," viable trajectories for", n.viable.pt," r-k pairs\n")
  }


#------------------------------------
# get results from CMSY
#------------------------------------
# get estimate of most probable r as 75th percentile of mid log.r-classes
# get unique combinations of r-k
unique.rk         <- unique(mdat.all[,1:2])
# get remaining viable log.r and log.k 
log.rs           <- log(unique.rk[,1])
log.ks           <- log(unique.rk[,2])
# get vectors with numbers of r and mid values in classes
# determine number of classes as a function of r-width
r.width         <- (max(unique.rk[,1])-start_r[1])/(start_r[2]-start_r[1])
classes         <- ifelse(r.width>0.8,100,ifelse(r.width>0.5,50,ifelse(r.width>0.3,25,12)))
hist.log.r      <- hist(x=log.rs, breaks=classes, plot=F)
log.r.counts    <- hist.log.r$counts
log.r.mids      <- hist.log.r$mids
# get most probable log.r as 75th percentile of mids with counts > 0
log.r.est       <- as.numeric(quantile(log.r.mids[which(log.r.counts > 0)],0.75))
lq.log.r        <- as.numeric(quantile(x=log.r.mids[which(log.r.counts > 0)], 0.625))
median.log.r    <- as.numeric(quantile(x=log.r.mids[which(log.r.counts > 0)], 0.50))
uq.log.r        <- as.numeric(quantile(x=log.r.mids[which(log.r.counts > 0)], 0.875))
lcl.log.r       <- as.numeric(quantile(x=log.r.mids[which(log.r.counts > 0)], 0.5125))
ucl.log.r       <- as.numeric(quantile(x=log.r.mids[which(log.r.counts > 0)], 0.9875))
r.est           <- exp(log.r.est)
lcl.r.est       <- exp(lcl.log.r)
ucl.r.est       <- exp(ucl.log.r)

# get r-k pairs above median of mids
rem            <- which(unique.rk[,1] > exp(median.log.r))
rem.log.r      <- log(unique.rk[,1][rem])
rem.log.k      <- log(unique.rk[,2][rem])
# do linear regression of log k ~ log r with slope fixed to -1 (from Schaefer)
reg         <- lm(rem.log.k ~ 1 + offset(-1*rem.log.r))
int.reg     <- as.numeric(reg[1])
sd.reg      <- sd(resid(reg))
# get estimate of log(k) from y where x = log.r.est
log.k.est   <- int.reg + (-1) * log.r.est
# get estimates of CL of log.k.est from y +/- SD where x = lcl.log.r or ucl.log.r 
lcl.log.k   <- int.reg + (-1) * ucl.log.r - sd.reg
ucl.log.k   <- int.reg + (-1) * lcl.log.r + sd.reg
k.est       <- exp(log.k.est)
lcl.k.est   <- exp(lcl.log.k)
ucl.k.est   <- exp(ucl.log.k)

#------------------------------------
# get results from CMSY
#------------------------------------
# get distribution for BESA:
allr      <- log.rs[which(log.rs>=lcl.log.r)]
allr      <- allr[which(allr<=ucl.log.r)]  # this is allr in logs
allk     <- (int.reg + (-1) * allr)       # this is allk in logs
distrCMSY<-elnorm(exp(allr), method = "mvue")
distkCMSY<-elnorm(exp(allk), method = "mvue")

mur=as.numeric(exp(distrCMSY$parameters[1]))
sigr=as.numeric(distrCMSY$parameters[2])
muk=as.numeric(exp(distkCMSY$parameters[1]))
sigk=as.numeric(distkCMSY$parameters[2])

allr1   <-log(exp(allr)*4)
distroneCMSY<-elnorm(exp(allr1), method = "mvue")
mur1=as.numeric(exp(distroneCMSY$parameters[1]))
sigr1=as.numeric(distroneCMSY$parameters[2])

output_list <- list(mur,sigr, muk, sigk,mur1)
# mur: CMSY estimate for r
# sigr: std. error associated with mur
# muk:  CMSY estimate for K --> m
# sigK: std. error associated with muk --> sigK
# mur1: transformed mur 

#detach(package:coda, unload = TRUE)
}

