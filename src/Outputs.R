# The whole script takes just over 3 hours

# Housekeeping
# Load libraries
library(dplyr)
library(readxl)
library(ggplot2)
library(tidyr)
library(stringr)
options(dplyr.summarise.inform = FALSE)

# Path to files
path_output <- "/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/PhD project/Project 2/Code base/RSV_MV_mAB_catchup/Outputs/" # Where model outputs are stored
path_data <- "/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/PhD project/Project 2/Code base/Data/" # Where the data is stored
path_code <- "/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/PhD project/Project 2/Code base/RSV_MV_mAB_catchup/src/"
path_fitting <- "/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/PhD project/Project 2/Code base/RSV_MV_mAB_catchup/Fitting outputs/"

# --------- Read in files ------------------------------------------------------
# Read in birth numbers for 2023
suppressMessages(births_de <- read_excel(paste0(path_data, "Births.xlsx")))
colnames(births_de)[2] <- "Month"
colnames(births_de)[5] <- "Value"
births_de <- births_de %>%
  select(Month, Value) %>%
  filter(grepl(".0", Value)) %>%
  mutate(Value = as.numeric(Value))

births_de <- births_de %>%
  mutate(season = case_when(
    Month == "January" | Month == "February" | Month == "December" ~ "winter",
    Month == "March" | Month == "April" | Month == "May" ~ "spring",
    Month == "June" | Month == "July" | Month == "August" ~ "summer",
    Month == "September" | Month == "October" | Month == "November" ~ "autumn"
  ))

# Get % births by season
births_de <- births_de %>%
  mutate(N = sum(Value)) %>%
  group_by(season, N) %>%
  reframe(
    total = sum(Value),
    prop = total / N
  ) %>%
  unique()

births_de <- births_de %>%
  rbind(
    births_de %>%
      summarise(
        season = "all",
        total = sum(total),
        prop = sum(prop),
        N = N[1]
      )
  )

# Read in hospitalisation posterior distribution from Mahmud et al
hosp_posterior <- readRDS(paste0(path_data, "inpatients_posterior_estimates.RDS"))

# Read-in ICU posterior distribution from Mahmud et al
icu_posterior <- readRDS(paste0(path_data, "icu_posterior_estimates.RDS"))


# ----------- Seroprevalence data -------------------------------------------------------------
# Read in the data from Andeweg et al and put it in the right format
data <- read.csv("https://raw.githubusercontent.com/Stijn-A/RSV_serology/refs/heads/master/data/infection_status.csv")
# Group age into intervals
# bi-monthly for 0-2 years and 6-monthly for 2-5 years
data$age_grp <- cut(data$age_days,
  breaks = c(
    seq(0, 730, by = 30.25 * 2),
    seq(909, 2000, by = 30.25 * 6)
  ),
  include.lowest = T, right = F
)

# Divide by season of birth
spring <- c(3, 4, 5)
summer <- c(6, 7, 8)
autumn <- c(9, 10, 11)
winter <- c(1, 2, 12)

data <- data %>%
  mutate(
    Birth_mo = format(as.Date(Birth_doy, origin = "2020-12-31"), "%m") %>% as.numeric(),
    season_birth = case_when(
      Birth_mo %in% spring ~ "spring",
      Birth_mo %in% summer ~ "summer",
      Birth_mo %in% autumn ~ "autumn",
      Birth_mo %in% winter ~ "winter"
    )
  )

# Define the midpoint of each age group
get_midpoint <- function(cut_label) {
  round(mean(as.numeric(unlist(strsplit(gsub("\\(|\\)|\\[|\\]", "", as.character(cut_label)), ",")))))
}

data$xMidpoint <- sapply(data$age_grp, get_midpoint)

# Get number of seroconversions by age and season of birth
incidence_data_season <- data %>%
  select(age_grp, age_days, infection, season_birth, xMidpoint) %>%
  mutate(N_tot = n()) %>%
  group_by(age_grp, xMidpoint, season_birth) %>%
  summarise(
    age_mid = round(median(age_days)),
    N = n(),
    n_infection = sum(infection),
    prop_seroconv = n_infection / N
  ) %>%
  ungroup() %>%
  distinct() %>%
  group_by(season_birth) %>%
  mutate(
    cum_pop = cumsum(N),
    incidence = n_infection / cum_pop
  )

# Calculate seroprevalence and 95% CI using binomial confidence intervals
incidence_data_season[, c("seroprev_mean", "seroprev_low95", "seroprev_up95")] <- binom::binom.confint(incidence_data_season$n_infection,
  incidence_data_season$N,
  method = "exact"
)[, c("mean", "lower", "upper")]

# Calculate incidence of seroconversion and binomial confidence intervals
incidence_data_season[, c("incidence_mean", "incidence_low95", "incidence_up95")] <- binom::binom.confint(incidence_data_season$n_infection,
  incidence_data_season$cum_pop,
  method = "exact"
)[, c("mean", "lower", "upper")]

# Pivot to wide format for fitting
incidence_data_season_wide <- incidence_data_season %>%
  select(!c(age_mid, age_grp, seroprev_mean, incidence_mean, cum_pop)) %>%
  pivot_wider(
    names_from = season_birth,
    values_from = c(
      N, n_infection, prop_seroconv, seroprev_low95, seroprev_up95,
      incidence, incidence_low95, incidence_up95
    ),
    values_fill = 0
  ) %>%
  rename(time = "xMidpoint")

# --------- Define our functions ----------------------------------------------
# Newly seroconverted at a given age
new_seroconv_age <- function(index) {
  # Extract a random iteration
  converted_all_rand <- as.numeric(t(converted_all[index, 0:366]))
  converted_sp_rand <- as.numeric(t(converted_sp[index, 0:366]))
  converted_sm_rand <- as.numeric(t(converted_sm[index, 0:366]))
  converted_au_rand <- as.numeric(t(converted_au[index, 0:366]))
  converted_wt_rand <- as.numeric(t(converted_wt[index, 0:366]))

  converted <- data.frame(
    age_midpoint = 0:365,
    sero_all = converted_all_rand,
    sero_sp = converted_sp_rand,
    sero_sm = converted_sm_rand,
    sero_au = converted_au_rand,
    sero_wt = converted_wt_rand
  )

  # Get incidence of seroconversion
  incidence_df <- data.frame(
    age_midpoint = converted$age_midpoint,
    incidence = c(0, diff(converted$sero_sp)),
    season_birth = "spring"
  ) %>%
    rbind(
      data.frame(
        age_midpoint = converted$age_midpoint,
        incidence = c(0, diff(converted$sero_sm)),
        season_birth = "summer"
      )
    ) %>%
    rbind(
      data.frame(
        age_midpoint = converted$age_midpoint,
        incidence = c(0, diff(converted$sero_au)),
        season_birth = "autumn"
      )
    ) %>%
    rbind(
      data.frame(
        age_midpoint = converted$age_midpoint,
        incidence = c(0, diff(converted$sero_wt)),
        season_birth = "winter"
      )
    ) %>%
    rbind(
      data.frame(
        age_midpoint = converted$age_midpoint,
        incidence = c(0, diff(converted$sero_all)),
        season_birth = "all"
      )
    )
}

# Number of infections
infections <- function(conversion, births) {
  # Convert ages to the same units
  conversion_formated <- conversion %>%
    mutate(
      age_months = floor(age_midpoint / 30), # turn age into months
      incidence = ifelse(incidence < 0, 0, incidence)
    ) %>% # incidence can't be negative
    arrange(age_midpoint) # arrange by age instead of by season of birth

  # Combine seroconversion data with birth data
  infections <- conversion_formated %>%
    merge(births %>% select(total, season),
      by.x = "season_birth", by.y = "season"
    )

  # Get number of infections
  infections <- infections %>%
    mutate(n_infections = incidence * total)

  return(infections)
}

# Aggregate infections into weekly infections to match disease progression data
infections_weekly <- function(infections) {
  infections_week <- infections %>%
    filter(season_birth != "all") %>%
    select(age_midpoint, age_months, n_infections, season_birth) %>%
    # We need to get age in weeks and the corresponding infections
    mutate(age_weeks = floor(age_midpoint / 7))

  return(infections_week)
}

# Save weeks to day mapping
weeks_to_day <- function(inf_week) {
  weeks_day <- inf_week %>%
    select(age_weeks, age_midpoint, age_months) %>%
    distinct() %>%
    arrange(age_midpoint)

  return(weeks_day)
}

# Get seasonal distribution of infections in each age group
inf_dist_seasonal <- function(infections_week) {
  # Total infections by age
  dist_infections_season <- infections_week %>%
    select(!age_midpoint) %>%
    group_by(age_weeks, season_birth) %>%
    mutate(n_infections = sum(n_infections)) %>%
    ungroup() %>%
    distinct() %>%
    # Get total infections by age in weeks
    group_by(age_weeks) %>%
    mutate(total_infections = sum(n_infections)) %>%
    ungroup() %>%
    # Get proportion of infections by season in a given week
    group_by(age_weeks) %>%
    mutate(
      prop_infections = ifelse(n_infections > 0, n_infections / total_infections, 0)
    ) %>%
    arrange(age_weeks) %>%
    ungroup() %>%
    select(age_weeks, age_months, season_birth, prop_infections)

  return(dist_infections_season)
}

# Get distribution of hospitalisations by age in weeks for <5 yo
hosp_prop <- function(hosp_posterior, i) {
  # sample from Mahmud et al distribution
  k <- hosp_posterior[i, "mu_k"]
  c <- hosp_posterior[i, "mu_c"]
  lambda <- hosp_posterior[i, "mu_lambda"]

  age_weeks <- seq(0, 260, by = 1)
  prop_hosp <- actuar::dburr(x = age_weeks, shape1 = k, shape2 = c, scale = lambda)
  prop_hosp <- data.frame(
    age_weeks = age_weeks,
    prop_hosp = prop_hosp
  )
  return(prop_hosp)
}

# Get number of hosp <1 in Germany by age
hosp_age_de <- function(prop_hosp) {
  # 21,000 in 2023 according to InEKDatenBrowser for <5 year, we need to divide
  # them by age group
  hosp_age_de <- prop_hosp %>%
    mutate(n_hosp = prop_hosp * 21000) %>%
    # Restrict to < 1 yo
    filter(age_weeks <= 52) %>%
    return(hosp_age_de)
}

# Get number of hospitalisations by age group and by season
hosp_age_season <- function(hosp_age_de, dist_infections_season, weeks_day) {
  hosp_age_season <- hosp_age_de %>%
    merge(dist_infections_season, by = "age_weeks") %>%
    group_by(age_weeks, season_birth) %>%
    mutate(n_hosp_season = n_hosp * prop_infections) %>%
    ungroup() %>%
    select(age_weeks, age_months, season_birth, n_hosp_season)

  # Add total by age
  hosp_age_season <- hosp_age_season %>%
    rbind(
      hosp_age_season %>%
        group_by(age_weeks, age_months) %>%
        summarise(
          season_birth = "all",
          n_hosp_season = sum(n_hosp_season)
        ) %>%
        ungroup()
    )

  # Put this back as daily numbers to match VE and IE distributions
  hosp_age_season <- hosp_age_season %>%
    merge(weeks_day, by = c("age_weeks", "age_months")) %>%
    select(age_midpoint, age_months, season_birth, n_hosp_season) %>%
    mutate(n_hosp_season = n_hosp_season / 7)

  return(hosp_age_season)
}

# Get distribution of ICU admissions in children <5 yo by age in weeks
icu_dist <- function(icu_posterior, i) {
  # sample from Mahmud et al distribution
  k <- icu_posterior[i, "mu_k"]
  c <- icu_posterior[i, "mu_c"]
  lambda <- icu_posterior[i, "mu_lambda"]

  age_weeks <- seq(0, 260, by = 1)
  icu_prop <- actuar::dburr(x = age_weeks, shape1 = k, shape2 = c, scale = lambda)
  icu_prop <- data.frame(
    age_weeks = age_weeks,
    icu_prop = icu_prop
  )
  return(icu_prop)
}

# Get number of ICU admissions <1 in Germany by age
icu_age_de <- function(prop_icu) {
  # 1,400 in 2023 according to InEKDatenBrowser for <5 year, we need to get
  # them by age group
  icu_age_de <- prop_icu %>%
    mutate(n_icu = icu_prop * 1400) %>%
    # Restrict them to <1
    filter(age_weeks <= 52)

  return(icu_age_de)
}

# Get number of ICU admissions by age group and by season
icu_age_season <- function(icu_age_de, dist_infections_season, weeks_day) {
  icu_age_season <- icu_age_de %>%
    merge(dist_infections_season, by = "age_weeks") %>%
    group_by(age_weeks, season_birth) %>%
    # The seasonal proportion of hospitalisations in each age group will be the same
    # as the proportion of infections
    mutate(n_icu_season = n_icu * prop_infections) %>%
    ungroup() %>%
    select(age_weeks, age_months, season_birth, n_icu_season)

  # Add total by age
  icu_age_season <- icu_age_season %>%
    rbind(
      icu_age_season %>%
        group_by(age_weeks, age_months) %>%
        summarise(
          season_birth = "all",
          n_icu_season = sum(n_icu_season)
        ) %>%
        ungroup()
    )

  # Put this back into daily numbers
  icu_age_season <- icu_age_season %>%
    merge(weeks_day, by = c("age_weeks", "age_months")) %>%
    select(age_midpoint, age_months, season_birth, n_icu_season) %>%
    mutate(n_icu_season = n_icu_season / 7)

  return(icu_age_season)
}

# ----- Interventions
# Effect of maternal vaccination on hospitalisations and ICU admissions
outcome_vacc <- function(cases, icu, VE_distribution) {
  # Add the estimated VE to the number of cases
  outcome_prevented_vacc <- cases %>%
    filter(!is.na(age_months)) %>%
    select(season_birth, age_months, age_midpoint, n_hosp_season) %>%
    merge(
      VE_distribution %>% filter(group == "severe") %>% select(t, VE_t),
      by.x = "age_midpoint", by.y = "t"
    ) %>%
    merge(
      icu %>% select(season_birth, age_months, age_midpoint, n_icu_season),
      by = c("season_birth", "age_months", "age_midpoint")
    )

  # Get number of prevented cases
  outcome_prevented_vacc <- outcome_prevented_vacc %>%
    group_by(season_birth, age_months, age_midpoint, VE_t) %>%
    summarise(
      n_vacc = n_hosp_season * (1 - VE_t), # number of cases despite vaccination
      n_averted = n_hosp_season * VE_t,
      n_icu = n_icu_season * (1 - VE_t), # number of ICU cases despite vaccination
      n_icu_averted = n_icu_season * VE_t
    ) %>%
    # Add age_bracket column for later
    mutate(age_bracket = as.character(age_months))

  # Get totals by season of birth
  outcome_prevented_vacc <- outcome_prevented_vacc %>%
    filter(season_birth != "all") %>%
    rbind(
      outcome_prevented_vacc %>%
        filter(season_birth != "all") %>%
        group_by(season_birth) %>%
        summarise(
          age_bracket = "all",
          age_months = NA,
          age_midpoint = NA,
          n_vacc = sum(n_vacc),
          n_averted = sum(n_averted),
          n_icu = sum(n_icu),
          n_icu_averted = sum(n_icu_averted)
        )
    )
  return(outcome_prevented_vacc)
}

# Effect of mAB on hospitalisations and ICU admissions
outcome_mAB_season <- function(cases, icu, VE_distribution, IE_distribution, 
                                   admin_age, season){
  # Only keep the IE from the time of administration to 12 months of age
  days_to_12_m <- round((12 - admin_age)*30.41)
  IE_subset <- IE_distribution %>%
    filter(group == "hospitalisation" &
             t <= days_to_12_m) %>%
    mutate(age = round(admin_age*30.41 + t),
           season_birth = season) %>%
    select(age, IE_t, season_birth)
  
  # Create a dataframe with IE = 0 for all previous ages (no effect as not given)
  IE_younger <- data.frame(age = seq(0, round(admin_age*30.41), by = 1),
                           IE_t = 0,
                           season_birth = season)
  
  # Combine the two dataframes
  IE_subset <- bind_rows(IE_younger, IE_subset)
  
  # Merge outcomes and VE and IE estimates
  outcomes <- cases %>%
    filter(!is.na(age_months)) %>%
    ungroup() %>%
    select(
      season_birth, age_months, age_midpoint, n_hosp_season) %>%
    # Add ICU admissions
    merge(
      icu %>% select(season_birth, age_months, age_midpoint, n_icu_season),
      by = c("season_birth", "age_months", "age_midpoint")
    ) 
  
  # Add VE and IE estimates
  # Immunised cohort first
  outcome_prevented <- outcomes %>%
    # Add VE estimates
    merge(
      VE_distribution %>% filter(group == "severe") %>% select(t, VE_t),
      by.x = "age_midpoint", by.y = "t"
    ) %>%
    # Add IE estimates
    merge(
      IE_subset,
      by.x = c("age_midpoint", "season_birth"), by.y = c("age", "season_birth"), keep_all.x = T
    ) 
  
    # Add IE = 0 for other seasons
  outcome_prevented <- outcome_prevented %>%
    rbind(
      outcomes %>%
        filter(!is.na(age_months) & season_birth != season) %>%
        ungroup() %>%
        select(season_birth, age_months, age_midpoint, n_hosp_season, n_icu_season) %>%
        merge(
          VE_distribution %>% filter(group == "severe") %>% select(t, VE_t),
          by.x = "age_midpoint", by.y = "t"
        ) %>%
        mutate(IE_t = 0)
      ) %>%
    arrange(age_midpoint)
  
  # Change VE to 0 when we have IE (effects don't add up)
  outcome_prevented <- outcome_prevented %>%
    mutate(VE_t = ifelse(IE_t > 0, 0, VE_t))
  
  # Add age_bracket for future DFs
  outcome_prevented <- outcome_prevented %>%
    mutate(age_bracket = as.character(age_months))
  
  # Apply estimates
  outcome_prevented_mAB <- outcome_prevented %>%
    group_by(season_birth, age_bracket, age_months, age_midpoint, IE_t) %>%
    mutate(
      n_nirs = n_hosp_season * (1 - VE_t) * (1 - IE_t),
      n_icu_nirs = n_icu_season * (1 - VE_t) * (1 - IE_t),
      # Cases averted compared to no intervention
      n_averted_nirs = n_hosp_season - (n_hosp_season * (1 - VE_t) * (1 - IE_t)),
      n_icu_averted_nirs = n_icu_season - (n_icu_season * (1 - VE_t) * (1 - IE_t))
    ) %>%
    arrange(season_birth, age_midpoint) %>%
    ungroup()
  
  
  return(outcome_prevented_mAB)
}

# Summarising the outputs of both interventions
outcome_intervention <- function(cases, icu, outcome_prevented_vacc, outcome_prevented_mAB) {
  # Get the numbers
  total_outcome_intervention <- cases %>%
    mutate(age_bracket = as.character(age_months)) %>%
    select(season_birth, age_bracket, n_hosp_season, age_midpoint) %>%
    rename(n_hospitalisations = n_hosp_season) %>%
    mutate(
      prevented_hospitalisations = 0,
      prevented_icu = 0
    ) %>%
    merge(
      icu %>%
        mutate(age_bracket = as.character(age_months)) %>%
        select(season_birth, age_bracket, n_icu_season, age_midpoint) %>%
        rename(n_icu = n_icu_season),
      by = c("season_birth", "age_bracket", "age_midpoint")
    ) %>%
    select(!age_midpoint) %>%
    mutate(intervention = "No immunisation") %>%
    rbind(
      outcome_prevented_vacc %>% filter(age_bracket == "all") %>%
        ungroup() %>%
        select(season_birth, age_bracket, n_vacc, n_averted, n_icu, n_icu_averted) %>%
        group_by(season_birth) %>%
        summarise(
          age_bracket = "all",
          n_hospitalisations = sum(n_vacc),
          prevented_hospitalisations = sum(n_averted),
          n_icu = sum(n_icu),
          prevented_icu = sum(n_icu_averted),
          intervention = "MV"
        )
    ) %>%
    rbind(
      outcome_prevented_vacc %>%
        filter(age_bracket != "all") %>%
        ungroup() %>%
        select(season_birth, age_bracket, n_vacc, n_averted, n_icu, n_icu_averted) %>%
        group_by(season_birth, age_bracket) %>%
        summarise(
          n_hospitalisations = sum(n_vacc),
          prevented_hospitalisations = sum(n_averted),
          n_icu = sum(n_icu),
          prevented_icu = sum(n_icu_averted),
          intervention = "MV"
        )
    ) %>%
    rbind(
      outcome_prevented_vacc %>% filter(age_bracket == "all") %>%
        ungroup() %>%
        select(age_bracket, n_vacc, n_averted, n_icu, n_icu_averted) %>%
        summarise(
          season_birth = "all",
          age_bracket = "all",
          n_hospitalisations = sum(n_vacc),
          prevented_hospitalisations = sum(n_averted),
          n_icu = sum(n_icu),
          prevented_icu = sum(n_icu_averted),
          intervention = "MV"
        )
    ) %>%
    rbind(
      outcome_prevented_mAB %>%
        select(season_birth, age_bracket, n_nirs, n_icu_nirs, n_averted_nirs,
               n_icu_averted_nirs, intervention) %>%
        group_by(season_birth, age_bracket, intervention) %>%
        summarise(
          n_hospitalisations = sum(n_nirs),
          prevented_hospitalisations = sum(n_averted_nirs),
          n_icu = sum(n_icu_nirs),
          prevented_icu = sum(n_icu_averted_nirs)
        )
    ) %>%
    rbind(
      outcome_prevented_mAB %>%
        select(age_bracket, n_nirs, n_icu_nirs, n_averted_nirs,
               n_icu_averted_nirs, intervention) %>%
        group_by(age_bracket, intervention) %>%
        summarise(
          season_birth = "all",
          n_hospitalisations = sum(n_nirs),
          prevented_hospitalisations = sum(n_averted_nirs),
          n_icu = sum(n_icu_nirs),
          prevented_icu = sum(n_icu_averted_nirs)
        )
    ) %>%
    rbind(
      outcome_prevented_mAB %>%
        select(season_birth, n_nirs, n_icu_nirs, n_averted_nirs,
               n_icu_averted_nirs, intervention) %>%
        group_by(season_birth, intervention) %>%
        summarise(
          age_bracket = "all",
          n_hospitalisations = sum(n_nirs),
          prevented_hospitalisations = sum(n_averted_nirs),
          n_icu = sum(n_icu_nirs),
          prevented_icu = sum(n_icu_averted_nirs)
        )
    ) %>%
    rbind(
      outcome_prevented_mAB %>%
        filter(intervention == "+ mAB for all") %>%
        select(n_nirs, n_icu_nirs, n_averted_nirs,
               n_icu_averted_nirs, intervention) %>%
        group_by(intervention) %>%
        summarise(
          age_bracket = "all",
          season_birth = "all",
          n_hospitalisations = sum(n_nirs),
          prevented_hospitalisations = sum(n_averted_nirs),
          n_icu = sum(n_icu_nirs),
          prevented_icu = sum(n_icu_averted_nirs)
        )
    ) %>%
    mutate(intervention = factor(intervention,
                                 levels = c(
                                   "No immunisation", "MV",
                                   "+ mAB for spring births",
                                   "+ mAB for summer births",
                                   "+ mAB for autumn births",
                                   "+ mAB for winter births",
                                   "+ mAB for all"
                                 )
    ))
  
  return(total_outcome_intervention)
}
# -----------------------------------------------------------------------------

# --------- Run the models once -----------------------------------------------
# Get VE
# source(paste0(path_code, "VE estimates.R")) # takes about 10 min
# or read in VE estimates if saved
model_ve_runs_t <- readRDS(paste0(path_fitting, "VE estimates.rds"))
# Get IE of mAB
# source(paste0(path_code, "mAB IE estimates.R")) # takes about 5 min
# or read in if saved
model_mab_runs_t <- readRDS(paste0(path_fitting, "mAB IE estimates.rds"))
# Get seroconversion by age
# source(paste0(path_code, "Main.R")) # takes about 3 hours
# or read in if saved
converted_summary <- readRDS(paste0(path_fitting, "seroconversion simulation outputs.rds"))
converted_all <- do.call(rbind.data.frame, converted_summary[1])
converted_sp <- do.call(rbind.data.frame, converted_summary[2])
converted_sm <- do.call(rbind.data.frame, converted_summary[3])
converted_au <- do.call(rbind.data.frame, converted_summary[4])
converted_wt <- do.call(rbind.data.frame, converted_summary[5])

# Replace negative values with 0
converted_all[converted_all < 0] <- 0
converted_sp[converted_sp < 0] <- 0
converted_sm[converted_sm < 0] <- 0
converted_au[converted_au < 0] <- 0
converted_wt[converted_wt < 0] <- 0


# -----------------------------------------------------------------------------

# --------- Combine the estimates ---------------------------------------------
# Define where we're going to store the outputs
total_outcome_intervention <- list()
infection_list <- list()

for (i in 1:100) {
  # ------ Calculate number of RSV cases by age in Germany ------------
  # Proportion of newly seroconverted children at a given age
  # Take one random iteration of the model
  set.seed(i)
  rand <- base::sample(1:10000, 1)
  # And extract the values for it
  conversion <- new_seroconv_age(rand)

  # Number of infections
  n_infections <- infections(conversion, births_de)
  # Number of infections by week
  n_infections_week <- infections_weekly(n_infections)
  # Save days/weeks/months mapping
  days_weeks_months <- weeks_to_day(n_infections_week)
  # Distribution of infections by season
  dist_inf_season <- inf_dist_seasonal(n_infections_week)

  # Number of hospitalisations
  index <- base::sample(1:500, 1) # take one random iteration of the posterior
  hosp_prob <- hosp_prop(hosp_posterior, index)
  hosp_cases_de <- hosp_age_de(hosp_prob)
  hosp_age_season_de <- hosp_age_season(hosp_cases_de, dist_inf_season, days_weeks_months)

  # Number of ICU admissions
  index <- base::sample(1:500, 1) # take one random iteration of the posterior
  icu_prob <- icu_dist(icu_posterior, index)
  icu_admissions_de <- icu_age_de(icu_prob)
  icu_age_season_de <- icu_age_season(icu_admissions_de, dist_inf_season, days_weeks_months)

  # ------ Maternal vaccination -----------------------------------------------
  # Take one random iteration of the model
  set.seed(42 * i)
  draw <- base::sample(1:10000, 1)
  # And extract the values for it
  VE_distribution <- model_ve_runs_t %>% filter(iter == draw)

  # Calculate number of outcomes despite the vaccine
  outcome_prevented_vacc <- outcome_vacc(hosp_age_season_de, icu_age_season_de, VE_distribution)

  # ------ Monoclonal antibody -----------------------------------------------
  # Take one random iteration of the model
  draw <- base::sample(1:10000, 1)
  # And extract the values for it
  IE_distribution <- model_mab_runs_t %>% filter(iter == draw)
  # Calculate number of outcomes despite the vaccine and the mAB
  outcome_prevented_mAB_autumn <- outcome_mAB_season(hosp_age_season_de, 
                                                             icu_age_season_de, 
                                                             VE_distribution, 
                                                             IE_distribution, 
                                                             admin_age = 1, 
                                                             season = "autumn") %>%
    mutate(intervention = "+ mAB for autumn births")
  
  outcome_prevented_mAB_summer <- outcome_mAB_season(hosp_age_season_de, 
                                                             icu_age_season_de, 
                                                             VE_distribution, 
                                                             IE_distribution, 
                                                             admin_age = 4, 
                                                             season = "summer") %>%
    mutate(intervention = "+ mAB for summer births")
  
  outcome_prevented_mAB_spring <- outcome_mAB_season(hosp_age_season_de, 
                                                             icu_age_season_de, 
                                                             VE_distribution, 
                                                             IE_distribution, 
                                                             admin_age = 8, 
                                                             season = "spring") %>%
    mutate(intervention = "+ mAB for spring births")
  
  outcome_prevented_mAB_winter <- outcome_mAB_season(hosp_age_season_de, 
                                                             icu_age_season_de, 
                                                             VE_distribution, 
                                                             IE_distribution, 
                                                             admin_age = 11, 
                                                             season = "winter") %>%
    mutate(intervention = "+ mAB for winter births")
  
  outcome_prevented_mAB_all <- outcome_prevented_mAB_autumn %>%
    filter(season_birth == "autumn") %>%
    rbind(
      outcome_prevented_mAB_spring %>% filter(season_birth == "spring"),
      outcome_prevented_mAB_summer %>% filter(season_birth == "summer"),
      outcome_prevented_mAB_winter %>% filter(season_birth == "winter")
    ) %>%
    mutate(intervention = "+ mAB for all")
  
  # Bind the mAB outcomes for all seasons together
  outcome_prevented_mAB <- outcome_prevented_mAB_autumn %>%
    rbind(
      outcome_prevented_mAB_spring,
      outcome_prevented_mAB_summer,
      outcome_prevented_mAB_winter,
      outcome_prevented_mAB_all
    )

  # Extract number of infections
  infection_list[[i]] <- n_infections %>%
    select(season_birth, age_midpoint, n_infections)
  infection_list[[i]]$iter <- rep(i, each = nrow(infection_list[[i]])) # save index

  # Get total outcomes by intervention
  total_outcome_intervention[[i]] <- outcome_intervention(hosp_age_season_de,
                                                          icu_age_season_de,
                                                          outcome_prevented_vacc,
                                                          outcome_prevented_mAB_alt)
  total_outcome_intervention[[i]]$iter <- rep(i, each = nrow(total_outcome_intervention[[i]])) # save index
}

total_outcome_intervention_df <- do.call("rbind", total_outcome_intervention)
infection_df <- do.call("rbind", infection_list)

total_outcome_intervention_df %>% readr::write_rds(paste0(path_output, "total outcome intervention.rds"))
infection_df %>% readr::write_rds(paste0(path_output, "infections.rds"))