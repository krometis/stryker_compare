#Functions to run Stryker


#Load the dataset
strykerData <- function(file) {
  #Insert file as string leading to stryker csv
  df = read.csv(file)
  
  df$VehicleID = 1
  df$AbortID = 1
  for(i in seq(2,nrow(df))) {
    # if TestPhase or VehicleVariant are the same then we might have the same vehicle
    if ( (df[i,"TestPhase"] == df[i-1,"TestPhase"]) &
         (df[i,"VehicleVariant"] == df[i-1,"VehicleVariant"]) ) {
      # if MilesBeforeSystemAbort = TestLifeMiles then we've started a new vehicle
      if (df[i,"MilesBeforeSystemAbort"] != df[i,"TestLifeMiles"]) {
        df[i,"VehicleID"] = df[i-1,"VehicleID"]
        df[i,"AbortID"] = df[i-1,"AbortID"]+1
      } else {
        df[i,"VehicleID"] = df[i-1,"VehicleID"]+1
        df[i,"AbortID"] = 1
      }
    }
    # if TestPhase or VehicleVariant change then we've started a new vehicle
    # in this case restart numbering at 1
    else {
      df[i,"VehicleID"] = 1
      df[i,"AbortID"] = 1
    }
  }
  
  #make a variable for right censored
  df$rightCensor = 1-df$ExactFailure
  
  #make a variable for MBSA to deal with missing (miles=0) data
  #copy from Pct50
  #have to do this down here because doing it earlier would mess up the vehicle id computation
  mbsaCol="Pct50"
  #MEV "PctXX" columns are NA for whatever reason, so initialize column and then copy over 
  #df$MBSA = df$MilesBeforeSystemAbort 
  #df$MBSA[df$VehicleVariant != "MEV"] = df[df$VehicleVariant != "MEV",mbsaCol]
  df$MBSA = df[,mbsaCol]
  df$MBSA[df$VehicleVariant == "MEV"] = df[df$VehicleVariant == "MEV","MilesBeforeSystemAbort"]
  
  #TODO: might also want to re-compute "TestLifeMiles" to match
  
  ## Additional columns that we didn't bother with in the hierarchical Bayesian model case ##
  
  # Add two factors for use in JAGS (add before split to ensure consistency since MEV isn't in DT)
  df$TestPhaseFactor = as.numeric(factor(df$TestPhase))      #test phase
  df$VehicleFactor   = as.numeric(factor(df$VehicleVariant)) #vehicle variant
  
  
  return(df)
}


plotMbsaHist = function(df,filename=NULL,variant=NULL,break.width=250) {
  #filter by variant if specified
  if (!is.null(variant)) {
    dfPlot = df[df$VehicleVariant==variant,]
    title.suffix = paste("(",variant,")",sep="")
  } else {
    dfPlot = df
    title.suffix = ""
  }
  #set filename if specified
  if (!is.null(filename)) png(filename = filename)
  #set histogram breaks
  brks=seq(0,max(dfPlot$MBSA)+break.width,break.width)
  xlab="Miles Before System Abort"
  
  plot1 = custom_hist(data = subset(dfPlot, ExactFailure == 1 & TestPhase == "DT")$MBSA,
                    xlab = xlab, ylab = 'Frequency', title = paste("DT Exact Failures",title.suffix), 
                    breaks = brks)
  
  plot2 = custom_hist(data = subset(dfPlot, ExactFailure== 1 & TestPhase == "OT")$MBSA,
                    xlab = xlab, ylab = 'Frequency', title = paste("OT Exact Failures",title.suffix),breaks=brks)
  
  plot3 = custom_hist(data = subset(dfPlot, ExactFailure== 0 & TestPhase == "DT")$MBSA,
                    xlab = xlab, ylab = 'Frequency', title = paste("DT Right-Censored Points",title.suffix),breaks=brks)
  
  plot4 = custom_hist(data = subset(dfPlot, ExactFailure== 0 & TestPhase == "OT")$MBSA,
                    xlab = xlab, ylab = 'Frequency', title = paste("OT-Right Censored Points",title.suffix),breaks=brks)
  
  grid.arrange(plot1, plot2, plot3, plot4, ncol = 2)
  
  
  #save to filename if specified
  if (!is.null(filename)) dev.off()
}

#sample MBSA (not mean!) from Weibull distributions given samples of beta and eta
sampleMbsa = function(data,etaIndex,nDraws,thinFactor=1) {
  nSamples = nrow(data$beta)
  thinIdx = seq(1,nSamples,thinFactor)
  return(as.vector(apply(as.data.frame(data)[thinIdx,],1,function(x) rweibull(nDraws,x[1],x[etaIndex+1]))))
}

#impliedMmbsa(): Given gamma priors on beta and eta, output draws from the MMBSA
#associated with Weibull(beta,eta)
impliedMmbsa = function(n = 1000, prior.params) {
  y_beta=rgamma(n, shape = prior.params$betashape, rate = prior.params$betarate)
  y_eta =rgamma(n, shape = prior.params$etashape , rate = prior.params$etarate )
  y_mmbsa = weibullMean(y_beta,y_eta)
  return(y_mmbsa)
}
impliedMmbsaStats = function(n = 1000, prior.params, variantNames) {
  #draw samples
  mmbsa.prior = impliedMmbsa(n,prior.params)

  #make a column for each variant
  df.prior = data.frame(mmbsa.prior)
  for (i in 2:length(variantNames)) { df.prior=cbind(df.prior,mmbsa.prior) }
  colnames(df.prior)=variantNames
  
  #compute stats and return
  return( summaryStats(df.prior) )
}

