mediation_analysis <- function(outcome, mediator, dataframe, sim.df, mc_iteration, reference_cohort, strata = FALSE){
  
  if (strata == FALSE){
    
    bootstrap_sample <- resample_data(dataframe, N = length(unique(dataframe$hhidpn)), ID_labels = "hhidpn")
    
    # Convert to numeric for splines (keep cohort as factor)
    bootstrap_sample$age <- as.numeric(bootstrap_sample$age)
    bootstrap_sample$period <- as.numeric(bootstrap_sample$period)
    
    # Fit outcome model
    outcome.model <- glm(paste(outcome, "~ Ns(age, df=6) + Ns(period, df=5, detrend=TRUE) + cohort +", 
                               mediator, "+ gender + raceth", sep=" "),
                         family = "binomial",
                         data = bootstrap_sample)
    
    # Copy simulated df
    sim.df.nc <- sim.df.cf <- data.frame(1:nrow(sim.df))
    
    # Fit mediator model (multinomial only)
    mediator.model <- multinom(paste(mediator, "~ Ns(age, df=6) + Ns(period, df=5, detrend=TRUE) + cohort + gender + raceth", 
                                     sep=""),
                               data = bootstrap_sample)
    
    sim.df$age <- as.numeric(sim.df$age)
    sim.df$period <- as.numeric(sim.df$period)
    
    simulated_mediator <- predict(mediator.model, type = "probs", newdata = sim.df)
    
    for (m in 1:mc_iteration){
      # Natural course
      sim.df[,mediator] <- rMultinom(m = 1, probs = simulated_mediator)
      sim.df.nc[,m] <- predict(outcome.model, type = "response", newdata = sim.df) 
      
      # Counterfactual
      sim.df[,mediator] <- sim.df[,mediator][sim.df$cohort == reference_cohort]
      sim.df.cf[,m] <- predict(outcome.model, type = "response", newdata = sim.df) 
    }
    
    # Aggregate MC results for natural course and counterfactual
    df.BS.nc <- rowMeans(sim.df.nc)
    df.BS.cf <- rowMeans(sim.df.cf)
    return(c(df.BS.nc, df.BS.cf))
  }
  
  if (strata == "raceth"){
    
    bootstrap_sample <- resample_data(dataframe, N = length(unique(dataframe$hhidpn)), ID_labels = "hhidpn")
    
    # Convert to numeric for splines (keep cohort as factor)
    bootstrap_sample$age <- as.numeric(bootstrap_sample$age)
    bootstrap_sample$period <- as.numeric(bootstrap_sample$period)
    
    # Fit outcome model
    outcome.model <- glm(paste(outcome, "~ Ns(age, df=6) + Ns(period, df=5, detrend=TRUE) + cohort +", 
                               mediator, "+ gender", sep=" "),
                         family = "binomial",
                         data = bootstrap_sample)
    
    # Copy simulated df
    sim.df.nc <- sim.df.cf <- data.frame(1:nrow(sim.df))
    
    # Fit mediator model (multinomial only)
    mediator.model <- multinom(paste(mediator, "~ Ns(age, df=6) + Ns(period, df=5, detrend=TRUE) + cohort + gender", 
                                     sep=""),
                               data = bootstrap_sample)
    
    sim.df$age <- as.numeric(sim.df$age)
    sim.df$period <- as.numeric(sim.df$period)
    
    simulated_mediator <- predict(mediator.model, type = "probs", newdata = sim.df)
    
    for (m in 1:mc_iteration){
      # Natural course
      sim.df[,mediator] <- rMultinom(m = 1, probs = simulated_mediator)
      sim.df.nc[,m] <- predict(outcome.model, type = "response", newdata = sim.df) 
      
      # Counterfactual
      sim.df[,mediator] <- sim.df[,mediator][sim.df$cohort == reference_cohort]
      sim.df.cf[,m] <- predict(outcome.model, type = "response", newdata = sim.df) 
    }
    
    # Aggregate MC results for natural course and counterfactual
    df.BS.nc <- rowMeans(sim.df.nc)
    df.BS.cf <- rowMeans(sim.df.cf)
    return(c(df.BS.nc, df.BS.cf))
  }
  
  if (strata == "gender"){
    
    bootstrap_sample <- resample_data(dataframe, N = length(unique(dataframe$hhidpn)), ID_labels = "hhidpn")
    
    # Convert to numeric for splines (keep cohort as factor)
    bootstrap_sample$age <- as.numeric(bootstrap_sample$age)
    bootstrap_sample$period <- as.numeric(bootstrap_sample$period)
    
    # Fit outcome model
    outcome.model <- glm(paste(outcome, "~ Ns(age, df=6) + Ns(period, df=5, detrend=TRUE) + cohort +", 
                               mediator, "+ raceth", sep=" "),
                         family = "binomial",
                         data = bootstrap_sample)
    
    # Copy simulated df
    sim.df.nc <- sim.df.cf <- data.frame(1:nrow(sim.df))
    
    # Fit mediator model (multinomial only)
    mediator.model <- multinom(paste(mediator, "~ Ns(age, df=6) + Ns(period, df=5, detrend=TRUE) + cohort + raceth", 
                                     sep=""),
                               data = bootstrap_sample)
    
    sim.df$age <- as.numeric(sim.df$age)
    sim.df$period <- as.numeric(sim.df$period)
    
    simulated_mediator <- predict(mediator.model, type = "probs", newdata = sim.df)
    
    for (m in 1:mc_iteration){
      # Natural course
      sim.df[,mediator] <- rMultinom(m = 1, probs = simulated_mediator)
      sim.df.nc[,m] <- predict(outcome.model, type = "response", newdata = sim.df) 
      
      # Counterfactual
      sim.df[,mediator] <- sim.df[,mediator][sim.df$cohort == reference_cohort]
      sim.df.cf[,m] <- predict(outcome.model, type = "response", newdata = sim.df) 
    }
    
    # Aggregate MC results for natural course and counterfactual
    df.BS.nc <- rowMeans(sim.df.nc)
    df.BS.cf <- rowMeans(sim.df.cf)
    return(c(df.BS.nc, df.BS.cf))
  }
}

library(compiler)
mediation_analysis <- cmpfun(mediation_analysis)