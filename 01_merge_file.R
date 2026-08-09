
options(scipen = 999)

library(tidyr)
library(haven)
library(dplyr)
library(stringr)
library(tidyselect)


#######set directory
setwd("C:/Users/rchen3/BOSTON UNIVERSITY Dropbox/Ruijia Chen/causal decomposition project/")

#######merge the data together
###read the rand sas wide file
hrs_wd <- read_sas("raw data/RAND HRS Longitudinal file/randhrs1992_2022v1.sas7bdat")
###convert the names to lower case
names(hrs_wd) <- tolower(names(hrs_wd))
####do not select spouse variables
hrs_wd <- hrs_wd %>%
  dplyr::select(vars_select(names(hrs_wd), !starts_with('s')), -hhidpn) %>%
  rename("hhidpn" = "rahhidpn")


## create a new dataset with variables of interest
base <- hrs_wd[, c("hhidpn", "ragender", "racohbyr", "rahispan", "raracem",
                   "rabmonth", "rabyear", "rabdate", "rabplace", "radmonth", "radyear", "raddate",
                   "rafeduc", "rameduc", "raeduc", "raedyrs", "raedegrm", "rarelig")]


##########read the mediator variables
hrs_wd_med <- hrs_wd %>% dplyr::select(hhidpn, starts_with("r") &
  ends_with(c("proxy", "iwstat", "wtresp", "age_y", "agey_b",
              "mrct", "shlt", "deprex", "depres", "cesd", "hibpe",
              "diabe", "cancre", "lunge", "lungf", "hearte",
              "stroke", "psyche", "psychf", "arthre", "sleepe",
              "drink", "drinkr", "drinkd", "drinkn", "smokev", "smoken", "vigact", "vgactf", "vgactx",
              "vigactn", "vigctp", "mdactx", "conde", "bmi", "adlf",
              "nsscre", "slfmem", "imrc20", "dlrc", "dlr10", "dlrc20",
              "ser7", "bwc20", "mo", "height", "wtresp")))


long_med <- hrs_wd_med %>%
  pivot_longer(
    cols = starts_with("r"),
    names_to = c("time", ".value"),
    names_pattern = "^r(\\d+)(.*)$",
    values_to = "Value"
  )


hrs_inw <- hrs_wd %>% dplyr::select(hhidpn, starts_with("inw"))

long_inw <- hrs_inw %>%
  pivot_longer(
    cols = starts_with("inw"),
    names_to = "time",
    names_pattern = "inw(\\d+)",
    values_to = "inw"
  ) %>% dplyr::select(hhidpn, time, inw)

long_med_df <- full_join(long_med, long_inw, by = c("hhidpn", "time"))


#####merge hrs_wd_med and base
hrs_long_df <- base %>% full_join(long_med_df, by = "hhidpn")

hrs_long_df <- hrs_long_df %>%
  mutate(
    time = as.numeric(time),
    year = 1990 + (time * 2)
  )


########read the dementia probability file
dem_dir <- "raw data/Dementia probability"

# List all subdirectories in the main directory
subdirs <- list.dirs(path = dem_dir, full.names = TRUE, recursive = FALSE)

###read all the sas files (no save action; metadata loop)
for (dir in subdirs) {
  sas_files <- list.files(path = dir, pattern = "\\.sas7bdat$", full.names = TRUE)
  for (sas_file in sas_files) {
    data <- read_sas(sas_file)
    file_base_name <- tools::file_path_sans_ext(basename(sas_file))
    rds_file <- file.path(dem_dir, paste0(file_base_name, ".rds"))
    # saveRDS(data, rds_file)
  }
}


#############################read the langa algorithm#############################
langa <- readRDS("raw data/Dementia probability/cogfinalimp_9520wide.rds") %>% rename_all(tolower)
langa$hhidpn <- paste(langa$hhid, langa$pn, sep = "")
langa$hhidpn <- as.character(langa$hhidpn)

langa_df <- langa %>%
  pivot_longer(
    cols = -c(hhid, pn, hhidpn),
    names_to = c(".value", "year"),
    names_pattern = "(.*?)(\\d{4})$"
  )

langa_df <- langa_df %>%
  rename_with(~ paste0("langa_", .), -c(hhid, pn, hhidpn, year))

langa_df$year <- as.character(langa_df$year)
langa_df <- langa_df %>% filter(!is.na(langa_cogfunction))
langa_df$year[langa_df$year == 1995] <- 1996


################################read melinda power's algorithm####################################################
power <- read_sas("raw data/Dementia probability/hrsdementia_2024_1217.sas7bdat") %>%
  rename_all(tolower) %>%
  rename(year = hrs_year)
power$hhidpn <- paste(power$hhid, power$pn, sep = "")
power$hhidpn <- as.character(power$hhidpn)
power <- power %>%
  rename_with(~ paste0("power_", .), -c(hhid, pn, hhidpn, year)) %>%
  rename(power_dem = power_expert_dem)


########################################merge the data ########################################

# merge dementia file
dementia_df <- power %>%
  full_join(langa_df, by = c("hhidpn", "year")) %>%
  rename(hhid = hhid.y, pn = pn.y)

####recode year
hrs_long_df$year <- as.character(hrs_long_df$year)


##############################merge dementia_df with the long data frame
all_df <- hrs_long_df %>%
  full_join(dementia_df, by = c("hhidpn", "year")) %>%
  rename(period = year, age = agey_b) %>%
  mutate(birthcohort = as.numeric(period) - as.numeric(age)) %>%
  dplyr::select(-hhid, -pn) %>%
  dplyr::select(hhidpn, period, age, birthcohort, everything())

all_df <- all_df %>% filter(!is.na(period))


# Filter to study period 2000-2020 and drop deaths / lost-to-follow-up
all_df1 <- all_df %>% filter(period >= 2000)
all_df1 <- all_df1 %>% filter(period < 2022) %>%
  group_by(hhidpn) %>%
  mutate(year_death = max(radyear),
         death = case_when(year_death <= period & is.na(age) ~ 1,
                           TRUE ~ 0))

all_df2 <- all_df1 %>% filter(death != 1)
all_df3 <- all_df2 %>% filter(!is.na(age))


###################### Combined inclusion: person-years with either Langa-Weir or
###################### Gianattasio-Power dementia indicator available
all_df4 <- all_df3 %>% filter(!is.na(langa_cogfunction) | !is.na(power_dem))

# Age limit 70-99 
all_df5 <- all_df4 %>% filter(age >= 70)
all_df5 <- all_df5 %>% filter(age < 100)


######save the combined dataset
saveRDS(all_df5, "data/hrs_all_2026_06_08.rds")
