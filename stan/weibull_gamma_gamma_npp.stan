//stan code to estimate the normalized power prior (NPP) posterior for the
//parameters of a weibull distribution given gamma priors 
//  beta: weibull shape parameter
//  eta:  weibull scale parameter (single prior, multiple posterior estimates)
//  npp_power: npp power parameter
//the approach follows that suggested in:
//  Ye, Keying, Zifei Han, Yuyan Duan, and Tianyu Bai. “Normalized Power Prior
//  Bayesian Analysis.” Journal of Statistical Planning and Inference 216
//  (January 1, 2022): 29–50. https://doi.org/10.1016/j.jspi.2021.05.005.

functions{
  //interpolate log(C) given a power and a vector of {power,log(C(power))} pairs
  real interpolateLogC(real nppPow, int nPowers, vector powerKnot, vector logCKnot){
    real logCest;
    for(id in 1:(nPowers-1)){
      if(nppPow >= powerKnot[id] && nppPow < powerKnot[id+1]){
        logCest = logCKnot[id]+ (nppPow-powerKnot[id])*(logCKnot[id+1]-logCKnot[id])/(powerKnot[id+1]-powerKnot[id]);
      }  // Interpolation function, given a sequence of logCKnot
    }
    return logCest;
  }
}
data {
  int<lower=0>  nData;
  int<lower=0>  nVariants;
  real<lower=0> mbsa[nData];
  int<lower=0>  is_cen[nData];
  int<lower=0>  variant[nData];
  int<lower=0,upper=1>  apply_pow[nData];
  real<lower=0> beta_shape;
  real<lower=0> beta_rate;
  real<lower=0> eta_shape;
  real<lower=0> eta_rate;
  int<lower=0>  nPowers;
  vector<lower=0>[nPowers] nppPowerKnot;
  vector[nPowers] nppLogCKnot;
}
parameters {
  real<lower=0> beta;
  real<lower=0> eta[nVariants];
  real<lower=0,upper=1> npp_power;
}
model {
  real logC; 

  // priors
  beta ~ gamma(beta_shape,beta_rate);
  for(j in 1:nVariants) {
    eta[j] ~ gamma(eta_shape, eta_rate);
  }
  // npp_power ~ beta(2,2);
  npp_power ~ uniform(0,1);
  
  // likelihood
  for(i in 1:nData) {
    // weibull scale parameter (stan uses R's parametrization)
    real scale  = eta[variant[i]];
    // effective power to apply (npp_power if apply_pow[i] is 1, 1 otherwise)
    //real pow = (apply_pow[i] == 1) ? npp_power : 1;
    real powr = 1 + apply_pow[i]*(npp_power-1);
    // if uncensored, compute (unnormalized) pdf at mbsa[i] directly
    if (is_cen[i] == 0)
      target += powr*weibull_lpdf( mbsa[i] | beta, scale );
    // if censored, compute probability that value is >= mbsa[i]
    if (is_cen[i] == 1)
      target += powr*weibull_lccdf( mbsa[i] | beta, scale );
  }
  
  logC =  interpolateLogC(npp_power, nPowers, nppPowerKnot, nppLogCKnot);
  target += -logC;
}
