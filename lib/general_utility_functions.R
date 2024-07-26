
#function to compute the shape and rate parameters for a gamma distribution 
#with given mean mu and standard deviation sigma
gammaParams = function(mu,sigma) {
  a = (mu/sigma)^2; #scale
  b = mu/(sigma^2); #rate
  return( list(shape=a,rate=b) )
}

priorParameters = function(beta.mean.pr=1, beta.std.pr=2,
                           eta.mean.pr=1000, eta.std.pr=2000){
  #beta parameters
  beta.mean.pr = beta.mean.pr; #beta = 1 is sort of the default
  beta.std.pr  = beta.std.pr; #large-ish stdev relative to the value
  beta.pr.params = gammaParams(beta.mean.pr,beta.std.pr)
  beta.pr.shape  = beta.pr.params$shape
  beta.pr.rate   = beta.pr.params$rate
  cat(paste("beta prior parameters:\n  shape =",beta.pr.shape,"\n  rate = ",beta.pr.rate,"\n"))
  
  #eta parameters
  eta.mean.pr = eta.mean.pr; #mbsa = 1000 is the goal. try to be "neutral" about this
  eta.std.pr  = eta.std.pr; #large-ish stdev relative to the value
  eta.pr.params = gammaParams(eta.mean.pr,eta.std.pr)
  eta.pr.shape  = eta.pr.params$shape
  eta.pr.rate   = eta.pr.params$rate
  cat(paste("eta prior parameters:\n  shape =",eta.pr.shape,"\n  rate = ",eta.pr.rate,"\n"))
  
  return(list(betashape = beta.pr.shape, betarate = beta.pr.rate,
              etashape = eta.pr.shape, etarate = eta.pr.rate))
}

weibullMean = function(p_shape,p_scale) {
  return( p_scale*gamma(1+(1/p_shape)) )
}

summaryStats = function(df) {#,group_by=variable) {
  mcmc.stats = melt(as.data.frame(df),id.vars = NULL) %>% 
    group_by(variable) %>% 
    summarize(
      mean  = mean(value), 
      sd    = sd(value), 
      median= median(value),
      min   = min(value), 
      q025  = quantile(value,0.025, na.rm = TRUE), 
      q10   = quantile(value,0.1  , na.rm = TRUE),
      q25   = quantile(value,0.25 , na.rm = TRUE),
      q75   = quantile(value,0.75 , na.rm = TRUE), 
      q90   = quantile(value,0.9  , na.rm = TRUE), 
      q975  = quantile(value,0.975, na.rm = TRUE), 
      max   = max(value)
    )
}

mmbsaSummaryStats = function(p_shape,p_scale,cols) {
  mmbsa.stats = as.data.frame( weibullMean(as.vector(p_shape),as.matrix(p_scale)) )
  colnames(mmbsa.stats) = cols
  mmbsa.stats = summaryStats(mmbsa.stats)
  return(mmbsa.stats)
}


computeWeibullPosterior = function(
                                   data,
                                   nVariants,   #needed because of case where data doesn't include one
                                   priorParams, #list of prior parameters
                                   stanFile,    #stan file to use
                                   chains   = 4,
                                   samples  = 50000,
                                   refresh  = 25000
                                  ) {
  #beta ~ gamma
  if ("betashape" %in% names(priorParams)) {
    #cat("detected gamma prior for beta\n")
    stan.data <- list(
      nVariants  = nVariants,
      nData      = nrow(data),           #number of data points
      mbsa       = data$MBSA,            #miles to system abort
      is_cen     = data$rightCensor,     #censored (0 or 1)
      variant    = data$VehicleFactor,   #vehicle variant
      #prior parameters
      beta_shape = priorParams$betashape,      
      beta_rate  = priorParams$betarate,
      eta_shape  = priorParams$etashape,      
      eta_rate   = priorParams$etarate
    )
  #beta ~ normal
  } else {
    #cat("detected normal prior for beta\n")
    stan.data <- list(
      nVariants  = nVariants,
      nData      = nrow(data),           #number of data points
      mbsa       = data$MBSA,            #miles to system abort
      is_cen     = data$rightCensor,     #censored (0 or 1)
      variant    = data$VehicleFactor,   #vehicle variant
      #prior parameters
      beta_mean  = priorParams$betamean,      
      beta_std   = priorParams$betastd,
      eta_shape  = priorParams$etashape,      
      eta_rate   = priorParams$etarate
    )
  } 
  
  # Run Stan
  fit.stan <- stan(
                   stanFile,
                   data     = stan.data, 
                   chains   = chains,
                   iter     = samples,
                   refresh  = refresh
                  )

  return( extract(fit.stan) ) #return mcmc results
  #return( fit.stan ) #return mcmc results
}

computeWeibullPosteriorNpp = function(
                                   data,
                                   nVariants,   #needed because of case where data doesn't include one
                                   priorParams, #list of prior parameters
                                   normConst,   #NPP normalizing constants
                                   stanFile,    #stan file to use
                                   chains   = 4,
                                   samples  = 50000,
                                   refresh  = 25000
                                  ) {
  #assume beta ~ gamma
  stan.data <- list(
    nVariants  = nVariants,
    nData      = nrow(data),           #number of data points
    mbsa       = data$MBSA,            #miles to system abort
    is_cen     = data$rightCensor,     #censored (0 or 1)
    variant    = data$VehicleFactor,   #vehicle variant
    apply_pow  = df$applyPower,        #apply power if test phase is DT
    #prior parameters
    beta_shape = priorParams$betashape,      
    beta_rate  = priorParams$betarate,
    eta_shape  = priorParams$etashape,      
    eta_rate   = priorParams$etarate,
    #npp parameters
    nPowers      = length(normConst$npp.powers),
    nppPowerKnot = normConst$npp.powers,
    nppLogCKnot  = normConst$logC
  )
  
  # Run Stan
  fit.stan <- stan(
                   stanFile,
                   data     = stan.data, 
                   chains   = chains,
                   iter     = samples,
                   refresh  = refresh
                  )

  return( extract(fit.stan) ) #return mcmc results
  #return( fit.stan ) #return mcmc results
}

#trapezoidal rule
trapRule = function (x,y) {
  delX = x[2:length(x)] - x[1:length(x)-1]
  return( 0.5*sum(delX*(y[2:length(y)]+y[1:length(y)-1])) )
}
#cumulative trapezoidal rule
cumTrapRule = function (x,y) {
  intY = rep(0,length(x))
  for (n in 2:length(x)) {
    intY[n] = trapRule(x[1:n],y[1:n])
  }
  return (intY)
}
  
computeWeibullNormalizationNpp = function(
                                   data,
                                   nVariants,      #needed because of case where data doesn't include one
                                   priorParams,    #list of prior parameters
                                   stanFile,       #stan file to use
                                   nKnots    = 50, #number of knots
                                   knotPower = 2,  #power > 1 (larger pushes more near 0)
                                   chains    = 4,
                                   samples   = 50000,
                                   refresh   = -1,
                                   verbose   = 0
                                  ) {
  #assume beta ~ gamma
  stan.data <- list(
    nVariants  = nVariants,
    nData      = nrow(data),           #number of data points
    mbsa       = data$MBSA,            #miles to system abort
    is_cen     = data$rightCensor,     #censored (0 or 1)
    variant    = data$VehicleFactor,   #vehicle variant
    #prior parameters
    beta_shape = priorParams$betashape,      
    beta_rate  = priorParams$betarate,
    eta_shape  = priorParams$etashape,      
    eta_rate   = priorParams$etarate,
    npp_power  = 0
  )

  #get knots
  npp.powers = (seq(0,nKnots)/nKnots)^knotPower

  #loop over power parameters and compute expected value
  expLlh     = rep(0,length(npp.powers))
  for (i in 1:length(npp.powers)) {
    cat(paste("Starting on power =",npp.powers[i],"\n"))
    stan.data$npp_power = npp.powers[i]
    
    # Run STAN to sample from the prior
    fit.stan <- stan(
                     stanFile,
                     data     = stan.data, 
                     chains   = chains,
                     iter     = samples,
                     refresh  = refresh
                    )

    #get expected value of the log likelihood from the mcmc results
    expLlh[i]    = mean(extract(fit.stan)$llh)
  
    #sanity checks
    if (verbose > 0) {
      cat(paste("  Posterior statistics (delta =",npp.powers[i],"):\n"))
      cat(paste("    beta: mean =",mean(extract(fit.stan)$beta),", stdev =",sd(extract(fit.stan)$beta),", median =",median(extract(fit.stan)$beta),", min =",min(extract(fit.stan)$beta),", max =",max(extract(fit.stan)$beta),"\n"))
      mcmc.eta = extract(fit.stan)$n_j
      for (j in 1:ncol(mcmc.eta)) {
        cat(paste("    eta:  mean =",mean(mcmc.eta[,j]),", stdev =",sd(mcmc.eta[,j]),", median =",median(mcmc.eta[,j]),", min =",min(mcmc.eta[,j]),", max =",max(mcmc.eta[,j]),"\n"))
      }
      cat(paste("    llh:  mean =",mean(extract(fit.stan)$llh),", stdev =",sd(extract(fit.stan)$llh),", median =",median(extract(fit.stan)$llh),", min =",min(extract(fit.stan)$llh),", max =",max(extract(fit.stan)$llh),"\n"))
      cat(paste("  expLlh = ",expLlh[i],"\n"))
    }
  }
  
  #integrate to get approximation to normalization factor ( C(\delta) from Ye 2022 )
  logC = cumTrapRule(npp.powers,expLlh)
  
  #return the knots and constants
  return( list(npp.powers = npp.powers, logC = logC) )
} 
