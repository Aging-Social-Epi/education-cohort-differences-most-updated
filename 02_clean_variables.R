
library(table1)
library(ggplot2)
library(ggsci)
library(dplyr)
library(tidyverse)


########read data (output from 01_merge_file.R)
hrs_all <- readRDS("data/hrs_all_2026_06_08.rds")

####rename raw RAND HRS parental education variables to align with downstream
####pedu construction (rafeduc -> fedu, rameduc -> medu)
hrs_all <- hrs_all %>%
  rename(fedu = rafeduc,
         medu = rameduc)

######reformat cohort, gender, race, education variables
hrs_all_v2 <- hrs_all %>%
  mutate(
    gender = factor(
      ragender,
      levels = c(1, 2),
      labels = c("Male", "Female")
    ),
    race = factor(
      raracem,
      levels = c(1, 2, 3),
      labels = c("White/Caucasian", "Black/African American", "Other")
    ),
    raceth = case_when(raracem == 1 & rahispan == 0 ~ 1,
                       raracem == 2 & rahispan == 0 ~ 2,
                       raracem == 3 & rahispan == 0 ~ 3,
                       rahispan == 1 ~ 4),
    raceth = factor(
      raceth,
      levels = c(1, 2, 3, 4),
      labels = c("NH White/Caucasian",
                 "NH Black/African American",
                 "NH Other",
                 "Hispanic")
    ),
    educat = factor(
      raeduc,
      levels = c(1, 2, 3, 4, 5),
      labels = c("Lt High-school",
                 "GED",
                 "High-school graduate",
                 "Some college",
                 "College and above")
    )
  ) %>% ungroup() %>%
  complete(hhidpn, period) %>%
  group_by(hhidpn) %>%
  arrange(hhidpn, period) %>%
  tidyr::fill(birthcohort, .direction = "updown") %>%
  ungroup() %>%
  mutate(birthcohort2 = cut(
    rabyear,
    breaks = c(1890, 1924, 1931, 1942, 1948, 1954, 1960, 1966, 1972, Inf),
    labels = c("AHEAD <1924", "CODA 1924-30", "HRS 1931-41", "WB 1942-47",
               "EBB 1948-53", "MBB 1954-59", "LBB 1960-65", "EGENX 1966-71", "Other"),
    right = FALSE
  ))


###### Langa-Weir dementia indicator (langa_di)
hrs_all_v2$langa_di <- NA
hrs_all_v2$langa_di[hrs_all_v2$langa_cogfunction %in% c(3)] <- 1
hrs_all_v2$langa_di[hrs_all_v2$langa_cogfunction %in% c(1, 2)] <- 0


###### Parental education years (pedu) -- max of father's (fedu) and mother's (medu)
hrs_all_v2 <- hrs_all_v2 %>%
  mutate(pedu = case_when(
    !is.na(fedu) & !is.na(medu) ~ pmax(fedu, medu),
    !is.na(fedu) &  is.na(medu) ~ fedu,
     is.na(fedu) & !is.na(medu) ~ medu,
    TRUE ~ NA_real_
  )) %>%
  group_by(hhidpn) %>%
  mutate(pedu = ifelse(all(is.na(pedu)), NA, max(pedu, na.rm = TRUE))) %>%
  ungroup()


###### Parental education category (pedu_cat) -- 5 levels based on years
hrs_all_v2$pedu_cat <- cut(
  as.numeric(hrs_all_v2$pedu),
  breaks = c(-Inf, 11, 12, 13, 15, 17),
  labels = c("11 or fewer years",
             "12 years",
             "13 years",
             "14-15 years",
             "16 or more years"),
  right = TRUE
)


###### Set factor and reference category for educat 
hrs_all_v2$educat <- factor(hrs_all_v2$educat)
hrs_all_v2$educat <- relevel(hrs_all_v2$educat, ref = "Lt High-school")


###### Save
saveRDS(hrs_all_v2, "data/hrs_all_clean_2026_06_08.rds")
