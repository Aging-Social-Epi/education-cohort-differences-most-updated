
mediation_analysis <- function(outcome,mediator, dataframe, sim.df, mc_iteration, reference_cohort , strata = FALSE){
  
  ## check to make sure the correct dataset is being used
  # if (nrow(dataframe) != 163760) {return(paste0("CHECK IF DATASET IS CORRECT! Current dataset has ", nrow(dataframe)," rows"))}
  
  if (strata == FALSE){
    ## sample from dataframe
    # bootstrap_sample  <- cfdecomp::cluster.resample(dataframe, cluster.name = "hhidpn",
    #                      size = length(unique(dataframe$hhidpn)))
    
    bootstrap_sample <-resample_data(dataframe, N = length(unique(dataframe$hhidpn)), ID_labels = "hhidpn")
    
    
    ## fit outcome model with APC specifications (Carstensen approach)
    bootstrap_sample$age <- as.numeric(bootstrap_sample$age)
    bootstrap_sample$period <- as.numeric(bootstrap_sample$period)
    #bootstrap_sample$cohort <- as.numeric(bootstrap_sample$cohort)
    outcome.model     <- glm(paste(outcome, "~ Ns(age, df=6) + Ns(period, df=5,detrend=TRUE) + cohort +", mediator, "+ gender + raceth+pedu_cat",sep=" "),
                             family ="binomial",
                             data = bootstrap_sample)
    
    ## copy simulated df
    sim.df.nc <- sim.df.cf <- data.frame(1:nrow(sim.df))
    
    
    ## fit mediator model
    if (class(dataframe[,mediator, drop = TRUE]) == "factor" & length(unique(dataframe[,mediator, drop = TRUE])) == 2){
      mediator.model <- glm(paste(mediator, "~ Ns(age, df=6) + Ns(period, df=5,detrend=TRUE) + cohort + gender  + raceth+pedu_cat",sep=""), 
                            family ="binomial",
                            data = bootstrap_sample)
      
      sim.df$age <- as.numeric(sim.df$age)
      sim.df$period <- as.numeric(sim.df$period)
      #sim.df$cohort <- as.numeric(sim.df$cohort)
      
      
      probs <- predict(mediator.model, type = "response", newdata = sim.df)
      
      ## start of Monte Carlo loop
      for (m in 1:mc_iteration){
        ## natural course
        ## predict mediator values
        sim.df[,mediator] <- as.numeric(rbinom(n = nrow(sim.df), size = 1,
                                              prob = probs))
        levels(sim.df[,mediator]) <- levels(dataframe[,mediator,drop = TRUE])
        sim.df.nc[,m]   <- predict(outcome.model, type="response", newdata = sim.df) 
        
        ## counterfactual
        ## our counterfactual is the education distribution of the 1945 cohort
        ## so we predict just with that cohort
        sim.df[,mediator] <- sim.df[,mediator][sim.df$cohort==reference_cohort]
        sim.df.cf[,m]  <- predict(outcome.model, type="response", newdata = sim.df) 
      }
    }   else if (class(dataframe[,mediator, drop = TRUE]) == "factor" & length(unique(dataframe[,mediator, drop = TRUE])) > 2){
      mediator.model <- multinom(paste(mediator, " ~ Ns(age, df=6) + Ns(period, df=5,detrend=TRUE) + cohort + gender + raceth+pedu_cat", sep=""),
                                 data = bootstrap_sample)
      sim.df$age <- as.numeric(sim.df$age)
      sim.df$period <- as.numeric(sim.df$period)
      #sim.df$cohort <- as.numeric(sim.df$cohort)
      simulated_mediator <- predict(mediator.model, type="probs",newdata = sim.df)
      
      for (m in 1:mc_iteration){
        #natural course
        #predict mediator values
        sim.df[,mediator] <- rMultinom(m = 1, probs = simulated_mediator)
        sim.df.nc[,m]   <- predict(outcome.model, type="response", newdata = sim.df) 
        
        #counterfactual
        sim.df[,mediator] <- sim.df[,mediator][sim.df$cohort==reference_cohort]
        sim.df.cf[,m]  <- predict(outcome.model, type="response", newdata = sim.df) 
      }
      
    } else if (class(dataframe[,mediator]) == "numeric"){
      mediator.model <- lm(bootstrap_sample[,mediator] ~ Ns(age, df=6) + 
                             Ns(period, df=5,detrend=TRUE) + 
                             cohort + gender + raceth+pedu_cat, 
                           data = bootstrap_sample)
      simulated_mediator  <- predict(mediator.model, type = "response", newdata = sim.df)
      
      for (m in 1:mc_iteration){
        #natural course
        #predict mediator values
        sim.df[,mediator]     <- rnorm(n=nrow(sim.df), mean = mean(simulated_mediator), sd=sd(mediator.model$residuals))
        sim.df.nc[,m]   <- predict(outcome.model, type="response", newdata = sim.df) 
        
        #counterfactual
        sim.df[,mediator] <- sim.df[,mediator][sim.df$cohort==reference_cohort]
        sim.df.cf[,m]  <- predict(outcome.model, type="response", newdata = sim.df) 
      }
    } else {
      return("Mediator has the wrong class")
    }
    
    ## aggregate MC results for natural course and counterfactual
    df.BS.nc <- rowMeans(sim.df.nc)
    df.BS.cf <- rowMeans(sim.df.cf)
    return(c(df.BS.nc, df.BS.cf))
  }
  
  if (strata =="raceth"){
    ## sample from dataframe
    # bootstrap_sample  <- cfdecomp::cluster.resample(dataframe, cluster.name = "hhidpn",
    #                      size = length(unique(dataframe$hhidpn)))
    
    bootstrap_sample <-resample_data(dataframe, N = length(unique(dataframe$hhidpn)), ID_labels = "hhidpn")
    
    
    ## fit outcome model with APC specifications (Carstensen approach)
    bootstrap_sample$age <- as.numeric(bootstrap_sample$age)
    bootstrap_sample$period <- as.numeric(bootstrap_sample$period)
    #bootstrap_sample$cohort <- as.numeric(bootstrap_sample$cohort)
    outcome.model     <- glm(paste(outcome, "~ Ns(age, df=6) + Ns(period, df=5,detrend=TRUE) + cohort +", mediator, "+ gender+pedu_cat",sep=" "),
                             family ="binomial",
                             data = bootstrap_sample)
    
    ## copy simulated df
    sim.df.nc <- sim.df.cf <- data.frame(1:nrow(sim.df))
    
    
    ## fit mediator model
    if (class(dataframe[,mediator, drop = TRUE]) == "factor" & length(unique(dataframe[,mediator, drop = TRUE])) == 2){
      mediator.model <- glm(paste(mediator, "~ Ns(age, df=6) + Ns(period, df=5,detrend=TRUE) + cohort + gender+pedu_cat",sep=""), 
                            family ="binomial",
                            data = bootstrap_sample)
      
      sim.df$age <- as.numeric(sim.df$age)
      sim.df$period <- as.numeric(sim.df$period)
      #sim.df$cohort <- as.numeric(sim.df$cohort)
      
      
      probs <- predict(mediator.model, type = "response", newdata = sim.df)
      
      ## start of Monte Carlo loop
      for (m in 1:mc_iteration){
        ## natural course
        ## predict mediator values
        sim.df[,mediator] <- as.numeric(rbinom(n = nrow(sim.df), size = 1,
                                               prob = probs))
        levels(sim.df[,mediator]) <- levels(dataframe[,mediator,drop = TRUE])
        sim.df.nc[,m]   <- predict(outcome.model, type="response", newdata = sim.df) 
        
        ## counterfactual
        ## our counterfactual is the health behavior factor distribution of the 1945 cohort
        ## so we predict just with that cohort
        sim.df[,mediator] <- sim.df[,mediator][sim.df$cohort==reference_cohort]
        sim.df.cf[,m]  <- predict(outcome.model, type="response", newdata = sim.df) 
      }
    }   else if (class(dataframe[,mediator, drop = TRUE]) == "factor" & length(unique(dataframe[,mediator, drop = TRUE])) > 2){
      mediator.model <- multinom(paste(mediator, " ~ Ns(age, df=6) + Ns(period, df=5,detrend=TRUE) + cohort + gender+pedu_cat", sep=""),
                                 data = bootstrap_sample)
      sim.df$age <- as.numeric(sim.df$age)
      sim.df$period <- as.numeric(sim.df$period)
      #sim.df$cohort <- as.numeric(sim.df$cohort)
      simulated_mediator <- predict(mediator.model, type="probs",newdata = sim.df)
      
      for (m in 1:mc_iteration){
        #natural course
        #predict mediator values
        sim.df[,mediator] <- rMultinom(m = 1, probs = simulated_mediator)
        sim.df.nc[,m]   <- predict(outcome.model, type="response", newdata = sim.df) 
        
        #counterfactual
        sim.df[,mediator] <- sim.df[,mediator][sim.df$cohort==reference_cohort]
        sim.df.cf[,m]  <- predict(outcome.model, type="response", newdata = sim.df) 
      }
      
    } else if (class(dataframe[,mediator]) == "numeric"){
      mediator.model <- lm(bootstrap_sample[,mediator] ~ Ns(age, df=6) + 
                             Ns(period, df=5,detrend=TRUE) + 
                             cohort + gender+pedu_cat, 
                           data = bootstrap_sample)
      simulated_mediator  <- predict(mediator.model, type = "response", newdata = sim.df)
      
      for (m in 1:mc_iteration){
        #natural course
        #predict mediator values
        sim.df[,mediator]     <- rnorm(n=nrow(sim.df), mean = mean(simulated_mediator), sd=sd(mediator.model$residuals))
        sim.df.nc[,m]   <- predict(outcome.model, type="response", newdata = sim.df) 
        
        #counterfactual
        sim.df[,mediator] <- sim.df[,mediator][sim.df$cohort==reference_cohort]
        sim.df.cf[,m]  <- predict(outcome.model, type="response", newdata = sim.df) 
      }
    } else {
      return("Mediator has the wrong class")
    }
    
    ## aggregate MC results for natural course and counterfactual
    df.BS.nc <- rowMeans(sim.df.nc)
    df.BS.cf <- rowMeans(sim.df.cf)
    return(c(df.BS.nc, df.BS.cf))
  }
  
  if (strata =="gender"){
    ## sample from dataframe
    # bootstrap_sample  <- cfdecomp::cluster.resample(dataframe, cluster.name = "hhidpn",
    #                      size = length(unique(dataframe$hhidpn)))
    
    bootstrap_sample <-resample_data(dataframe, N = length(unique(dataframe$hhidpn)), ID_labels = "hhidpn")
    
    
    ## fit outcome model with APC specifications (Carstensen approach)
    bootstrap_sample$age <- as.numeric(bootstrap_sample$age)
    bootstrap_sample$period <- as.numeric(bootstrap_sample$period)
    #bootstrap_sample$cohort <- as.numeric(bootstrap_sample$cohort)
    outcome.model     <- glm(paste(outcome, "~ Ns(age, df=6) + Ns(period, df=5,detrend=TRUE) + cohort +", mediator, "+ raceth+pedu_cat",sep=" "),
                             family ="binomial",
                             data = bootstrap_sample)
    
    ## copy simulated df
    sim.df.nc <- sim.df.cf <- data.frame(1:nrow(sim.df))
    
    
    ## fit mediator model
    if (class(dataframe[,mediator, drop = TRUE]) == "factor" & length(unique(dataframe[,mediator, drop = TRUE])) == 2){
      mediator.model <- glm(paste(mediator, "~ Ns(age, df=6) + Ns(period, df=5,detrend=TRUE) + cohort + raceth+pedu_cat",sep=""), 
                            family ="binomial",
                            data = bootstrap_sample)
      
      sim.df$age <- as.numeric(sim.df$age)
      sim.df$period <- as.numeric(sim.df$period)
      #sim.df$cohort <- as.numeric(sim.df$cohort)
      
      
      probs <- predict(mediator.model, type = "response", newdata = sim.df)
      
      ## start of Monte Carlo loop
      for (m in 1:mc_iteration){
        ## natural course
        ## predict mediator values
        sim.df[,mediator] <- as.numeric(rbinom(n = nrow(sim.df), size = 1,
                                               prob = probs))
        levels(sim.df[,mediator]) <- levels(dataframe[,mediator,drop = TRUE])
        sim.df.nc[,m]   <- predict(outcome.model, type="response", newdata = sim.df) 
        
        ## counterfactual
        ## our counterfactual is the health behavior factor distribution of the 1945 cohort
        ## so we predict just with that cohort
        sim.df[,mediator] <- sim.df[,mediator][sim.df$cohort==reference_cohort]
        sim.df.cf[,m]  <- predict(outcome.model, type="response", newdata = sim.df) 
      }
    }   else if (class(dataframe[,mediator, drop = TRUE]) == "factor" & length(unique(dataframe[,mediator, drop = TRUE])) > 2){
      mediator.model <- multinom(paste(mediator, " ~ Ns(age, df=6) + Ns(period, df=5,detrend=TRUE) + cohort + raceth+pedu_cat", sep=""),
                                 data = bootstrap_sample)
      sim.df$age <- as.numeric(sim.df$age)
      sim.df$period <- as.numeric(sim.df$period)
      #sim.df$cohort <- as.numeric(sim.df$cohort)
      simulated_mediator <- predict(mediator.model, type="probs",newdata = sim.df)
      
      for (m in 1:mc_iteration){
        #natural course
        #predict mediator values
        sim.df[,mediator] <- rMultinom(m = 1, probs = simulated_mediator)
        sim.df.nc[,m]   <- predict(outcome.model, type="response", newdata = sim.df) 
        
        #counterfactual
        sim.df[,mediator] <- sim.df[,mediator][sim.df$cohort==reference_cohort]
        sim.df.cf[,m]  <- predict(outcome.model, type="response", newdata = sim.df) 
      }
      
    } else if (class(dataframe[,mediator]) == "numeric"){
      mediator.model <- lm(bootstrap_sample[,mediator] ~ Ns(age, df=6) + 
                             Ns(period, df=5,detrend=TRUE) + 
                             cohort + raceth+pedu_cat, 
                           data = bootstrap_sample)
      simulated_mediator  <- predict(mediator.model, type = "response", newdata = sim.df)
      
      for (m in 1:mc_iteration){
        #natural course
        #predict mediator values
        sim.df[,mediator]     <- rnorm(n=nrow(sim.df), mean = mean(simulated_mediator), sd=sd(mediator.model$residuals))
        sim.df.nc[,m]   <- predict(outcome.model, type="response", newdata = sim.df) 
        
        #counterfactual
        sim.df[,mediator] <- sim.df[,mediator][sim.df$cohort==reference_cohort]
        sim.df.cf[,m]  <- predict(outcome.model, type="response", newdata = sim.df) 
      }
    } else {
      return("Mediator has the wrong class")
    }
    
    ## aggregate MC results for natural course and counterfactual
    df.BS.nc <- rowMeans(sim.df.nc)
    df.BS.cf <- rowMeans(sim.df.cf)
    return(c(df.BS.nc, df.BS.cf))
  }


}

# Byte-compile for speedup (10-30% on tight MC loop, bit-for-bit identical)
library(compiler)
mediation_analysis <- cmpfun(mediation_analysis)
