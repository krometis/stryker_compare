//stan code to estimate the normalized power prior (NPP) normalization factor 
//associated with the parameters of a weibull distribution given gamma priors 
//  beta: weibull shape parameter
//  eta:  weibull scale parameter (single prior, multiple posterior estimates)
//  npp_power: npp power parameter
//the approach follows that suggested in:
//  Ye, Keying, Zifei Han, Yuyan Duan, and Tianyu Bai. “Normalized Power Prior
//  Bayesian Analysis.” Journal of Statistical Planning and Inference 216
//  (January 1, 2022): 29–50. https://doi.org/10.1016/j.jspi.2021.05.005.

data {
  int<lower=0>  nData;
  int<lower=0>  nVariants;
  real<lower=0> mbsa[nData];
  int<lower=0>  is_cen[nData];
  int<lower=0>  variant[nData];
  real<lower=0> beta_shape;
  real<lower=0> beta_rate;
  real<lower=0> eta_shape;
  real<lower=0> eta_rate;
  real<lower=0> npp_power;
}
parameters {
  real<lower=0> beta;
  real<lower=0> eta[nVariants];
}
transformed parameters {
  //compute log likelihood here so that it's saved at the end of the run
  //(we need it to compute the normalization factor for NPP)
  real llh; 

  //compute the log likelihood
  llh = 0;
  for(i in 1:nData) {
    // weibull scale parameter (stan uses R's parametrization)
    real scale  = eta[variant[i]];
    // if uncensored, compute (unnormalized) pdf at mbsa[i] directly
    if (is_cen[i] == 0)
      llh += weibull_lpdf( mbsa[i] | beta, scale );
      //mbsa[i] ~ weibull( beta, scale );
    // if censored, compute probability that value is >= mbsa[i]
    if (is_cen[i] == 1)
      llh += weibull_lccdf( mbsa[i] | beta, scale );
  }
}
model {
  // priors
  beta ~ gamma(beta_shape,beta_rate);
  for(j in 1:nVariants) {
    eta[j] ~ gamma(eta_shape, eta_rate);
  }

  // likelihood
  target += npp_power*llh;
}
