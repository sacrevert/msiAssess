## JAGS code from https://github.com/BiologicalRecordsCentre/BRCindicators for reference only

smoothing <- function() {'
    ######################### Smoothing #######################

    beta[1] ~ dnorm(0, 0.000001)
    beta[2] ~ dnorm(0, 0.000001)
    taub ~ dgamma(0.000001, 0.000001)
    for(k in 1:num.knots){b[k]~dnorm(0,taub)}
  
    for (t in 1:(nyears - 1)){
      logLambda[t] <- m[t]
      m[t] <- mfe[t] + mre[t]
      mfe[t] <- beta[1] * X[t,1] + beta[2] * X[t,2]
      for (k in 1:num.knots){
        temp[t,k] <- b[k] * Z[t,k]
      }
      mre[t]<-sum(temp[t,1:num.knots])
    }  
  '
}

################################################################################

bma_model_Smooth <- function(incl.2deriv = FALSE, 
                             seFromData = FALSE,
                             Y1perfect = TRUE){
  
  # BUGS code assembled by metacoding
  # Indicator defined by Growth rates, with Ruppert smoother (deterministic version)
  # Defined by equation 7 in Steve Freeman's document of 2/11/17
  # this version is suitable for datasets where species have zero SE in year 1
  # and some species area allowed to join late
  # Therefore the indicator must be plotted with zero error in year 1 
  # uncertainty in M[1] doesn't measure the same thing as uncertainty in other years 
  # tau.spi is on the growth rates, not the index
  # Mprime is now estimated without uncertainty due to interspecific variation
  # standard errors can be read from a file or fixed
  
  ########### functions
  
  priors <- function(seFromData = FALSE){
    process_errors <- '
  # process errors
  tau.spi <- pow(sigma.spi,-2)
  sigma.spi ~ dunif(0.0001,30)'
    
    if(seFromData){
      obsErrors <- '
  # observation errors: one value per site-species
  for (s in 1:nsp){
   for (t in 1:nyears){
    sigma.obs[s,t] ~ dunif(0.0001, max_se) # for the missing values
  }}
    '
    } else {
      obsErrors <- '
  # observation error is constant
  theta ~ dunif(0.0001,10)'
    }
    
    return(paste(c(
      "    ###################  Define priors  ###########################",
      process_errors, 
      obsErrors),
      collapse = "\n"
    ))
  }
  
  likelihood <- function(seFromData = FALSE, Y1perfect = TRUE) {
    # option: is the first year assumed to be known without error, or not?
    part1 <- "
    M[1] <- 0 
  
    for (t in 2:nyears){
      M[t] <- M[t-1] + logLambda[t-1]
    }"
    
    part2 <-"
    for (s in 1:nsp){
      for (t in 1:(nyears-1)){
        spgrowth[s,t] ~ dnorm(logLambda[t], tau.spi)
    }}"
    
    if(Y1perfect){ # y1perfect = T
      part3 <- "
    for (s in 1:nsp){
      for (t in 1:(FY[s]-1)){ # e.g. FY[s] = 10 and nyears = 20, loop then 1:9
        spindex[s,t] <- spindex[s,t+1] - spgrowth[s,t] # t = 9 means spindex for year 10 minus backwards 1 growth step from estimate[,10]
    }

    spindex[s,FY[s]] <- estimate[s,FY[s]] # we assume the first year is known without error

    for (t in (FY[s]+1):(nyears)){ # from year 11 in example
      spindex[s,t] <- estimate[s,FY[s]] + sum(spgrowth[s,FY[s]:(t-1)]) # index in 10 = year 10 anchor plus sum of growths from 10 to 9
      estimate[s,t] ~ dnorm(spindex[s,t], tau.obs[s,t]) # estimate yr 11 drawn from spindex for yr 11 etc.
      tau.obs[s,t] <- pow(theta, -2)
    }}"
    } else { # y1perfect = F
      part3 <- "
    for (s in 1:nsp){
      for (t in 1:FY[s]){ # years 1 to 10 in example
        spindex[s,t] <- spindex[s,t+1] - spgrowth[s,t] # t = 10 means spindex for yr 11 - growth step for year 10
      }

      for (t in (FY[s]+1):(nyears)){ # years 11 to nyears (20)
        spindex[s,t] <- estimate[s,FY[s]] + sum(spgrowth[s,FY[s]:(t-1)]) # year 10 estimate + sum of growth from 10 to 9
        estimate[s,t] ~ dnorm(spindex[s,t], tau.obs[s,t])
        tau.obs[s,t] <- pow(theta, -2)
      }}"
    }
    if(seFromData) part3 <- gsub(part3, pattern ="theta", replacement ="sigma.obs[s,t]")
    
    return(paste(c(
      "    ###################  Define likelihood  ###########################",
      part1, 
      part2,
      part3),
      collapse = "\n"
    ))
  }
  
  geomeanExpected <- function(){'
  ####################  geomean of expected values ######################
  
  for (t in 1:nyears){
    Mprime[t] <- sum(spindex[,t])/nsp
  }'
  }
  
  ##########
  
  derivatives <- ifelse(!incl.2deriv, "", {'
    #########################  second derivatives #######################
  
  I <- M
  t2dash[1] <- 0
  t2dash[2] <- (I[2+1] - 2*I[2] + I[2-1])/1
  t2dash[nyears-1] <- (I[nyears] - 2*I[nyears-1] + I[nyears-2])/1
  t2dash[3] <- (-I[5]+16*I[4]-30*I[3]+16*I[2]-I[1] )/12
  t2dash[nyears-2] <- (-I[nyears]+16*I[nyears-1]-30*I[nyears-2]+16*I[nyears-3]-I[nyears-4])/12
  for (t in 4:(nyears-3)){
    t2dash[t]<-(2*I[t+3]-27*I[t+2]+270*I[t+1]-490*I[t]+270*I[t-1]-27*I[t-2]+2*I[t-3])/180
  }
  
    #####################################################################
  
  '})
  
  model <- paste(c("model {",
                   priors(seFromData), 
                   smoothing(), 
                   likelihood(seFromData, Y1perfect),
                   geomeanExpected(),
                   derivatives,
                   "}"), collapse = "\n")
  
  return(model)
}


foo <- bma_model_Smooth(incl.2deriv = FALSE, seFromData = FALSE, Y1perfect = TRUE)
