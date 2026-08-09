# =============================================================================
# Sensitivity (NH White) -- additionally adjusted for parental education.
# Runs LW + GP. Joint distribution samples (gender, educat, pedu_cat) tuples.
# =============================================================================

if (dir.exists("~/R/library")) .libPaths("~/R/library")

options(error = function() {
  traceback(2)
  dump.frames("error_dump", to.file = TRUE)
  q("no", 1, FALSE)
})

library(foreach); library(parallel); library(doParallel); library(foreign)
library(dplyr); library(tidyselect); library(lattice); library(Epi); library(ggplot2)
library(nnet); library(MASS); library(Hmisc); library(cfdecomp); library(fabricatr)
library(survey); library(questionr)
options(scipen = 999)

source('/restricted/project/gly-ukb/Ruijia/function_education_sen.R')


data_date <- "2026_06_08"
out_date  <- "2026_06_08"
dataset_dir <- paste0("/restricted/project/gly-ukb/Ruijia/data/hrs_all_clean_", data_date, ".rds")

APC_df <- readRDS(dataset_dir) %>%
  filter(birthcohort2 %in% c("AHEAD <1924", "CODA 1924-30", "HRS 1931-41", "WB 1942-47")) %>%
  filter(raceth == "NH White/Caucasian")

APC_df$cohort <- droplevels(factor(APC_df$birthcohort2))

APC_df <- APC_df %>% dplyr::select(hhidpn, langa_di, power_dem, iwstat,
                                   age, period, cohort, gender, educat, raceth,
                                   wtresp, pedu_cat) %>%
  filter(!is.na(age), !is.na(period), !is.na(cohort), !is.na(gender),
         !is.na(educat), !is.na(raceth), !is.na(pedu_cat))


APC_df_full <- APC_df


# =============================================================================
# (1) Langa-Weir analysis (primary)
# =============================================================================

APC_df <- APC_df_full %>% filter(!is.na(langa_di))

n_obs <- 100
n_cohort <- length(unique(APC_df$cohort))
n_period <- length(unique(APC_df$period))
n_age <- length(unique(APC_df$age))

sim <- expand.grid(age=sort(unique(APC_df$age)),
                   period=sort(unique(APC_df$period)),
                   cohort=sort(unique(APC_df$cohort)))
sim <- sim[rep(seq_len(nrow(sim)), each = n_obs), ]

n_cohort<-length(unique(APC_df$cohort))
n_age<-length(unique(APC_df$age))
n_period<-length(unique(APC_df$period))

groups <- n_cohort*n_age*n_period
sim_length <- n_obs*n_cohort*n_age*n_period
sim_length_cohort <- n_obs*n_age*n_period

gender <- educat <- pedu_cat <- rep(NA_character_, sim_length)
cohort_levels <- levels(APC_df$cohort)

set.seed(123)
for (i in 1:n_cohort){
  a <-(1+((i-1)*sim_length_cohort))
  b <-(sim_length_cohort+((i-1)*sim_length_cohort))

  cohort_data <- APC_df[APC_df$cohort == cohort_levels[i],
                        c("gender", "educat", "pedu_cat")]
  idx <- sample(nrow(cohort_data), size = sim_length_cohort, replace = TRUE)

  gender[a:b]   <- as.character(cohort_data$gender)[idx]
  educat[a:b]   <- as.character(cohort_data$educat)[idx]
  pedu_cat[a:b] <- as.character(cohort_data$pedu_cat)[idx]
}

sim1 <- as.data.frame(cbind(gender, educat, pedu_cat), stringsAsFactors = FALSE)
sim2 <- cbind(sim, sim1)
sim2$gender   <- factor(sim2$gender,   levels = levels(APC_df$gender))
sim2$educat   <- factor(sim2$educat,   levels = levels(APC_df$educat))
sim2$pedu_cat <- factor(sim2$pedu_cat, levels = levels(APC_df$pedu_cat))

edu.sim <- sim2


cl <- makeCluster(20)
registerDoParallel(cl)

bs_it = 1000
mc_it = 1000

set.seed(123)
start.time <- Sys.time()

start.time
edu <- foreach(bs = 1:bs_it, .combine='cbind', .inorder=FALSE, .packages=c("cfdecomp","nnet","Epi","Hmisc","fabricatr")) %dopar% {
  mediation_analysis(outcome="langa_di", mediator="educat", dataframe=APC_df, sim.df=edu.sim, mc_iteration=mc_it, reference_cohort = "WB 1942-47",
                     strata="raceth")
}
end.time <- Sys.time()
end.time

saveRDS(edu, paste0("/restricted/project/gly-ukb/Ruijia/Output/sen/parental/output_all_white_LW_", out_date, ".RDS"))

edu.nc <- edu[1:(nrow(sim)),]
edu.cf <- edu[((nrow(sim)+1)):(2*nrow(sim)),]


df_mean <- cbind(sim[1:n_obs,], edu.cf[1:n_obs,])
edu.cf2<- df_mean %>% group_by(age, period, cohort) %>%
  summarise_at(.vars = names(.)[4:(bs_it+3)], .funs = c(mean="mean"))

for(i in 2:groups){
  a <- (1+n_obs*(i-1))
  b <- (n_obs+n_obs*(i-1))
  df_mean <- cbind(sim[a:b,], edu.cf[a:b,])
  df_sum<- df_mean %>% group_by(age, period, cohort) %>%
    summarise_at(.vars = names(.)[4:(bs_it+3)], .funs = c(mean="mean"))
  edu.cf2 <-rbind(edu.cf2,df_sum)
}

df_mean <- cbind(sim[1:n_obs,], edu.nc[1:n_obs,])
edu.nc2<- df_mean %>% group_by(age, period, cohort) %>%
  summarise_at(.vars = names(.)[4:(bs_it+3)], .funs = c(mean="mean"))

for(i in 2:groups){
  a <- (1+n_obs*(i-1))
  b <- (n_obs+n_obs*(i-1))
  df_mean <- cbind(sim[a:b,], edu.nc[a:b,])
  df_sum<- df_mean %>% group_by(age, period, cohort) %>%
    summarise_at(.vars = names(.)[4:(bs_it+3)], .funs = c(mean="mean"))
  edu.nc2 <-rbind(edu.nc2,df_sum)
}

saveRDS(edu.nc2, paste0("/restricted/project/gly-ukb/Ruijia/Output/sen/parental/output_nc_white_LW_", out_date, ".RDS"))
saveRDS(edu.cf2, paste0("/restricted/project/gly-ukb/Ruijia/Output/sen/parental/output_cf_white_LW_", out_date, ".RDS"))


# =============================================================================
# (2) Gianattasio-Power analysis (co-primary)
# =============================================================================

APC_df <- APC_df_full %>% filter(!is.na(power_dem))

n_cohort <- length(unique(APC_df$cohort))
n_period <- length(unique(APC_df$period))
n_age <- length(unique(APC_df$age))

sim <- expand.grid(age=sort(unique(APC_df$age)),
                   period=sort(unique(APC_df$period)),
                   cohort=sort(unique(APC_df$cohort)))
sim <- sim[rep(seq_len(nrow(sim)), each = n_obs), ]

n_cohort<-length(unique(APC_df$cohort))
n_age<-length(unique(APC_df$age))
n_period<-length(unique(APC_df$period))

groups <- n_cohort*n_age*n_period
sim_length <- n_obs*n_cohort*n_age*n_period
sim_length_cohort <- n_obs*n_age*n_period

gender <- educat <- pedu_cat <- rep(NA_character_, sim_length)
cohort_levels <- levels(APC_df$cohort)

set.seed(123)
for (i in 1:n_cohort){
  a <-(1+((i-1)*sim_length_cohort))
  b <-(sim_length_cohort+((i-1)*sim_length_cohort))

  cohort_data <- APC_df[APC_df$cohort == cohort_levels[i],
                        c("gender", "educat", "pedu_cat")]
  idx <- sample(nrow(cohort_data), size = sim_length_cohort, replace = TRUE)

  gender[a:b]   <- as.character(cohort_data$gender)[idx]
  educat[a:b]   <- as.character(cohort_data$educat)[idx]
  pedu_cat[a:b] <- as.character(cohort_data$pedu_cat)[idx]
}

sim1 <- as.data.frame(cbind(gender, educat, pedu_cat), stringsAsFactors = FALSE)
sim2 <- cbind(sim, sim1)
sim2$gender   <- factor(sim2$gender,   levels = levels(APC_df$gender))
sim2$educat   <- factor(sim2$educat,   levels = levels(APC_df$educat))
sim2$pedu_cat <- factor(sim2$pedu_cat, levels = levels(APC_df$pedu_cat))

edu.sim <- sim2


set.seed(123)
start.time <- Sys.time()

start.time
edu <- foreach(bs = 1:bs_it, .combine='cbind', .inorder=FALSE, .packages=c("cfdecomp","nnet","Epi","Hmisc","fabricatr")) %dopar% {
  mediation_analysis(outcome="power_dem", mediator="educat", dataframe=APC_df, sim.df=edu.sim, mc_iteration=mc_it, reference_cohort = "WB 1942-47",
                     strata="raceth")
}
end.time <- Sys.time()
end.time

saveRDS(edu, paste0("/restricted/project/gly-ukb/Ruijia/Output/sen/parental/output_all_white_GP_", out_date, ".RDS"))

edu.nc <- edu[1:(nrow(sim)),]
edu.cf <- edu[((nrow(sim)+1)):(2*nrow(sim)),]


df_mean <- cbind(sim[1:n_obs,], edu.cf[1:n_obs,])
edu.cf2<- df_mean %>% group_by(age, period, cohort) %>%
  summarise_at(.vars = names(.)[4:(bs_it+3)], .funs = c(mean="mean"))

for(i in 2:groups){
  a <- (1+n_obs*(i-1))
  b <- (n_obs+n_obs*(i-1))
  df_mean <- cbind(sim[a:b,], edu.cf[a:b,])
  df_sum<- df_mean %>% group_by(age, period, cohort) %>%
    summarise_at(.vars = names(.)[4:(bs_it+3)], .funs = c(mean="mean"))
  edu.cf2 <-rbind(edu.cf2,df_sum)
}

df_mean <- cbind(sim[1:n_obs,], edu.nc[1:n_obs,])
edu.nc2<- df_mean %>% group_by(age, period, cohort) %>%
  summarise_at(.vars = names(.)[4:(bs_it+3)], .funs = c(mean="mean"))

for(i in 2:groups){
  a <- (1+n_obs*(i-1))
  b <- (n_obs+n_obs*(i-1))
  df_mean <- cbind(sim[a:b,], edu.nc[a:b,])
  df_sum<- df_mean %>% group_by(age, period, cohort) %>%
    summarise_at(.vars = names(.)[4:(bs_it+3)], .funs = c(mean="mean"))
  edu.nc2 <-rbind(edu.nc2,df_sum)
}

saveRDS(edu.nc2, paste0("/restricted/project/gly-ukb/Ruijia/Output/sen/parental/output_nc_white_GP_", out_date, ".RDS"))
saveRDS(edu.cf2, paste0("/restricted/project/gly-ukb/Ruijia/Output/sen/parental/output_cf_white_GP_", out_date, ".RDS"))


