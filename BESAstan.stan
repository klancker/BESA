/*----------------------- Data --------------------------*/
/* Data block: unpack the objects that will be input as data */
data {
int TT; // Length of state and observation time series
vector[TT] pt; // Log price observations
vector[TT] Ht; // Level Harvest observations
// lognormal and normal prior inputs (first mu, second sig) 
real mur;  
real sigr; 
real sigK; 
real musigeps;
real sigsigeps;
real musigeta;
real sigsigeta;
real sigc;
// inputs for beta distribution of starting depletion 
real F11;
real F12;
// fixed inputs 
real m;  // expected value for K, multiplier for better numerical performance
real mc; // expected value for c, multiplier for better numerical performance
}
/*----------------------- Parameters --------------------------*/
parameters {
real<lower=0>   r;    // parameter in process equation
real<lower=0>   K;     // parameter in process equation
real<lower=0> sigeps; // standard deviation of the observation equation error
real<lower=0> sigeta; // standard deviation of the process equation error
real<lower=0> c;     // parameter in observation equation
real<lower=0, upper=1> Xs;   // initial state of depletion
real<lower=0> Z1;		// initial state variance 
}

// transform some paramters for better numerical performance or to output (MSY)
transformed parameters {
real sigepstrue;
real sigetatrue;
real truer;
real truek; 
real truec;
real MSY;
sigepstrue=sigeps/4;
sigetatrue=sigeta/4;
truer=r/4;
truek=K*m;
truec=c*mc;
MSY=truer*truek/4; 
}
/*----------------------- Model --------------------------*/
/* Model block: defines and runs the model */
model {
// Initialization of variables
vector[TT+1] xt; // state variable time series
vector[TT+1] Zt; // variance of xt
vector[TT] xtt; // intermed state var time series
vector[TT] vt; // residual time series
vector[TT] Ft; // variance of vt
vector[TT] Ztt; // intermed state var variance 
vector[TT] Jtt; // Jacobian
real llik_obs[TT];
real llik;

// Prior distributions: 
target+=lognormal_lpdf(c|log(1),sigc);
target+=lognormal_lpdf(r|log(mur),sigr);
target+=lognormal_lpdf(K|log(1),sigK);
target+=normal_lpdf(sigeps|musigeps,sigsigeps);
target+=normal_lpdf(sigeta|musigeta,sigsigeta);
target+=beta_lpdf(Xs|F11,F12);
target+=lognormal_lpdf(Z1|log(0.01),1);

// initial: 
Zt[1]=Z1;    // initial variance
xt[1]=log(Xs);  // initial state
target+= -0.5*TT*log(2*3.14159);  // log likelihood constant part

// Extended Kalman filter part loop 
for (t in 1:TT){
vt[t]=pt[t]-log(truec)+log(K*m)+xt[t];   // residual
Ft[t]= Zt[t]+pow(sigepstrue,2);   // variance of residual 
// Each step, add this part to loglikelihood: 
llik_obs[t] = -0.5*(log(Ft[t])+pow(vt[t],2)/Ft[t]); // log likelihood

xtt[t]=xt[t]-(Zt[t]*vt[t])/Ft[t];   // corrected state estimate
Ztt[t]=pow(sigepstrue,2)*Zt[t]/Ft[t];    // variance of corrected state estimate
xt[t+1]=log((1+truer)*exp(xtt[t])-truer*exp(2*xtt[t])-Ht[t]/truek);  // a priori state estimate for next period
Jtt[t]=((1+truer)*exp(xtt[t])-2*truer*exp(2*xtt[t]))/((1+truer)*exp(xtt[t])-truer*exp(2*xtt[t])-Ht[t]/truek);  // Jacobian
Zt[t+1]=pow(Jtt[t],2)*Ztt[t]+pow(sigetatrue,2);  // a priori estimate for state variance for next period
}
llik = sum(llik_obs); // log likelihood is the sum over all the parts 
target+=llik;  // add the log likelihood to the target
}

// the next part only serves to make quantities available for R output afterwards
/*----------------------- Generated quantities (for output)--------------------------*/
generated quantities {
vector[TT+1] xt; 
vector[TT+1] Bt; // biomass time series
vector[TT] vt; 
vector[TT] Ft; 
vector[TT+1] Zt; 
vector[TT] xtt; 
vector[TT] Ztt; 
vector[TT] Jtt; 
real log_lik[TT];
real llik;

// initial: 
Zt[1]=Z1;
xt[1]=log(Xs);
Bt[1]=exp(xt[1])*truek;

// Kalman filter part loop 
for (t in 1:TT) {
vt[t]=pt[t] -log(truec)+log(truek)+xt[t];
Ft[t]= Zt[t]+pow(sigepstrue,2);
log_lik[t] = -0.5*(log(Ft[t])+pow(vt[t],2)/Ft[t]); 
xtt[t]=xt[t]-(Zt[t]*vt[t])/Ft[t]; 
Ztt[t]=pow(sigepstrue,2)*Zt[t]/Ft[t]; 
xt[t+1]=log((1+truer)*exp(xtt[t])-truer*exp(2*xtt[t])-Ht[t]/truek);
Jtt[t]=((1+truer)*exp(xtt[t])-2*truer*exp(2*xtt[t]))/((1+truer)*exp(xtt[t])-truer*exp(2*xtt[t])-Ht[t]/truek);
Zt[t+1]=pow(Jtt[t],2)*Ztt[t]+pow(sigetatrue,2);
Bt[t+1]=exp(xt[t+1])*truek;
}
llik = sum(log_lik); 

}
