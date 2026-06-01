# ------------- Compartmental model --------------------------------------------
## ------------------ Percentage is each compartment -----------------------------
traj <- readRDS("/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/LSTHM project/Extension/RSV_infection_NL/RDS files/trajectories monty season.rds")
converted_summary <- readRDS("/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/LSTHM project/Extension/RSV_infection_NL/RDS files/seroconversion simulation outputs.rds")

# Extract each compartment
M1_sp <- do.call(rbind.data.frame, traj[1])
M2_sp <- do.call(rbind.data.frame, traj[5])
M_sp <- M1_sp + M2_sp
S_sp <- do.call(rbind.data.frame, traj[9])
R_sp <- do.call(rbind.data.frame, converted_summary[2])

M1_sm <- do.call(rbind.data.frame, traj[2])
M2_sm <- do.call(rbind.data.frame, traj[6])
M_sm <- M1_sm + M2_sm
S_sm <- do.call(rbind.data.frame, traj[10])
R_sm <- do.call(rbind.data.frame, converted_summary[3])

M1_au <- do.call(rbind.data.frame, traj[3])
M2_au <- do.call(rbind.data.frame, traj[7])
M_au <- M1_au + M2_au
S_au <- do.call(rbind.data.frame, traj[11])
R_au <- do.call(rbind.data.frame, converted_summary[4])

M1_wt <- do.call(rbind.data.frame, traj[4])
M2_wt <- do.call(rbind.data.frame, traj[8])
M_wt <- M1_wt + M2_wt
S_wt <- do.call(rbind.data.frame, traj[12])
R_wt <- do.call(rbind.data.frame, converted_summary[5])

M1_all <- 0.26 * traj$M1_sp +
  0.29 * traj$M1_sm +
  0.24 * traj$M1_au +
  0.20 * traj$M1_wt

M2_all <- 0.26 * traj$M2_sp +
  0.29 * traj$M2_sm +
  0.24 * traj$M2_au +
  0.20 * traj$M2_wt

M_all <- M1_all + M2_all

S_all <- 0.26 * traj$S_sp +
  0.29 * traj$S_sm +
  0.24 * traj$S_au +
  0.20 * traj$S_wt

R_all <- do.call(rbind.data.frame, converted_summary[1])

# Transpose to have time as rows and convert to data frame, focusing on first 366 days
M1_sp <- as.data.frame(t(M1_sp[, 0:366]))
M2_sp <- as.data.frame(t(M2_sp[, 0:366]))
M_sp <- as.data.frame(t(M_sp[, 0:366]))
S_sp <- as.data.frame(t(S_sp[, 0:366]))
R_sp <- as.data.frame(t(R_sp[, 0:366]))

M1_sm <- as.data.frame(t(M1_sm[, 0:366]))
M2_sm <- as.data.frame(t(M2_sm[, 0:366]))
M_sm <- as.data.frame(t(M_sm[, 0:366]))
S_sm <- as.data.frame(t(S_sm[, 0:366]))
R_sm <- as.data.frame(t(R_sm[, 0:366]))

M1_au <- as.data.frame(t(M1_au[, 0:366]))
M2_au <- as.data.frame(t(M2_au[, 0:366]))
M_au <- as.data.frame(t(M_au[, 0:366]))
S_au <- as.data.frame(t(S_au[, 0:366]))
R_au <- as.data.frame(t(R_au[, 0:366]))

M1_wt <- as.data.frame(t(M1_wt[, 0:366]))
M2_wt <- as.data.frame(t(M2_wt[, 0:366]))
M_wt <- as.data.frame(t(M_wt[, 0:366]))
S_wt <- as.data.frame(t(S_wt[, 0:366]))
R_wt <- as.data.frame(t(R_wt[, 0:366]))

M1_all <- as.data.frame(t(M1_all[, 0:366]))
M2_all <- as.data.frame(t(M2_all[, 0:366]))
M_all <- as.data.frame(t(M_all[, 0:366]))
S_all <- as.data.frame(t(S_all[, 0:366]))
R_all <- as.data.frame(t(R_all[, 0:366]))

# Add time
time <- seq(0, 365, by = 1)
M1_sp <- cbind(time, M1_sp)
M2_sp <- cbind(time, M2_sp)
M_sp <- cbind(time, M_sp)
S_sp <- cbind(time, S_sp)
R_sp <- cbind(time, R_sp)
M1_sm <- cbind(time, M1_sm)
M2_sm <- cbind(time, M2_sm)
M_sm <- cbind(time, M_sm)
S_sm <- cbind(time, S_sm)
R_sm <- cbind(time, R_sm)
M1_au <- cbind(time, M1_au)
M2_au <- cbind(time, M2_au)
M_au <- cbind(time, M_au)
S_au <- cbind(time, S_au)
R_au <- cbind(time, R_au)
M1_wt <- cbind(time, M1_wt)
M2_wt <- cbind(time, M2_wt)
M_wt <- cbind(time, M_wt)
S_wt <- cbind(time, S_wt)
R_wt <- cbind(time, R_wt)
M1_all <- cbind(time, M1_all)
M2_all <- cbind(time, M2_all)
M_all <- cbind(time, M_all)
S_all <- cbind(time, S_all)
R_all <- cbind(time, R_all)

# Get 95% CIs
M1_sp_CI <- M1_sp %>%
  as.data.frame() %>%
  pivot_longer(-time, names_to = "sim", values_to = "value") %>%
  group_by(time) %>%
  summarise(
    lower = quantile(value, 0.025),
    median = quantile(value, 0.5),
    upper = quantile(value, 0.975)
  )

M2_sp_CI <- M2_sp %>%
  as.data.frame() %>%
  pivot_longer(-time, names_to = "sim", values_to = "value") %>%
  group_by(time) %>%
  summarise(
    lower = quantile(value, 0.025),
    median = quantile(value, 0.5),
    upper = quantile(value, 0.975)
  )

M_sp_CI <- M_sp %>%
  as.data.frame() %>%
  pivot_longer(-time, names_to = "sim", values_to = "value") %>%
  group_by(time) %>%
  summarise(
    lower = quantile(value, 0.025),
    median = quantile(value, 0.5),
    upper = quantile(value, 0.975)
  )

S_sp_CI <- S_sp %>%
  as.data.frame() %>%
  pivot_longer(-time, names_to = "sim", values_to = "value") %>%
  group_by(time) %>%
  summarise(
    lower = quantile(value, 0.025),
    median = quantile(value, 0.5),
    upper = quantile(value, 0.975)
  )

R_sp_CI <- R_sp %>%
  as.data.frame() %>%
  pivot_longer(-time, names_to = "sim", values_to = "value") %>%
  group_by(time) %>%
  summarise(
    lower = quantile(value, 0.025),
    median = quantile(value, 0.5),
    upper = quantile(value, 0.975)
  )

M1_sm_CI <- M1_sm %>%
  as.data.frame() %>%
  pivot_longer(-time, names_to = "sim", values_to = "value") %>%
  group_by(time) %>%
  summarise(
    lower = quantile(value, 0.025),
    median = quantile(value, 0.5),
    upper = quantile(value, 0.975)
  )

M2_sm_CI <- M2_sm %>%
  as.data.frame() %>%
  pivot_longer(-time, names_to = "sim", values_to = "value") %>%
  group_by(time) %>%
  summarise(
    lower = quantile(value, 0.025),
    median = quantile(value, 0.5),
    upper = quantile(value, 0.975)
  )

M_sm_CI <- M_sm %>%
  as.data.frame() %>%
  pivot_longer(-time, names_to = "sim", values_to = "value") %>%
  group_by(time) %>%
  summarise(
    lower = quantile(value, 0.025),
    median = quantile(value, 0.5),
    upper = quantile(value, 0.975)
  )

S_sm_CI <- S_sm %>%
  as.data.frame() %>%
  pivot_longer(-time, names_to = "sim", values_to = "value") %>%
  group_by(time) %>%
  summarise(
    lower = quantile(value, 0.025),
    median = quantile(value, 0.5),
    upper = quantile(value, 0.975)
  )

R_sm_CI <- R_sm %>%
  as.data.frame() %>%
  pivot_longer(-time, names_to = "sim", values_to = "value") %>%
  group_by(time) %>%
  summarise(
    lower = quantile(value, 0.025),
    median = quantile(value, 0.5),
    upper = quantile(value, 0.975)
  )

M1_au_CI <- M1_au %>%
  as.data.frame() %>%
  pivot_longer(-time, names_to = "sim", values_to = "value") %>%
  group_by(time) %>%
  summarise(
    lower = quantile(value, 0.025),
    median = quantile(value, 0.5),
    upper = quantile(value, 0.975)
  )

M2_au_CI <- M2_au %>%
  as.data.frame() %>%
  pivot_longer(-time, names_to = "sim", values_to = "value") %>%
  group_by(time) %>%
  summarise(
    lower = quantile(value, 0.025),
    median = quantile(value, 0.5),
    upper = quantile(value, 0.975)
  )

M_au_CI <- M_au %>%
  as.data.frame() %>%
  pivot_longer(-time, names_to = "sim", values_to = "value") %>%
  group_by(time) %>%
  summarise(
    lower = quantile(value, 0.025),
    median = quantile(value, 0.5),
    upper = quantile(value, 0.975)
  )

S_au_CI <- S_au %>%
  as.data.frame() %>%
  pivot_longer(-time, names_to = "sim", values_to = "value") %>%
  group_by(time) %>%
  summarise(
    lower = quantile(value, 0.025),
    median = quantile(value, 0.5),
    upper = quantile(value, 0.975)
  )

R_au_CI <- R_au %>%
  as.data.frame() %>%
  pivot_longer(-time, names_to = "sim", values_to = "value") %>%
  group_by(time) %>%
  summarise(
    lower = quantile(value, 0.025),
    median = quantile(value, 0.5),
    upper = quantile(value, 0.975)
  )

M1_wt_CI <- M1_wt %>%
  as.data.frame() %>%
  pivot_longer(-time, names_to = "sim", values_to = "value") %>%
  group_by(time) %>%
  summarise(
    lower = quantile(value, 0.025),
    median = quantile(value, 0.5),
    upper = quantile(value, 0.975)
  )

M2_wt_CI <- M2_wt %>%
  as.data.frame() %>%
  pivot_longer(-time, names_to = "sim", values_to = "value") %>%
  group_by(time) %>%
  summarise(
    lower = quantile(value, 0.025),
    median = quantile(value, 0.5),
    upper = quantile(value, 0.975)
  )

M_wt_CI <- M_wt %>%
  as.data.frame() %>%
  pivot_longer(-time, names_to = "sim", values_to = "value") %>%
  group_by(time) %>%
  summarise(
    lower = quantile(value, 0.025),
    median = quantile(value, 0.5),
    upper = quantile(value, 0.975)
  )

S_wt_CI <- S_wt %>%
  as.data.frame() %>%
  pivot_longer(-time, names_to = "sim", values_to = "value") %>%
  group_by(time) %>%
  summarise(
    lower = quantile(value, 0.025),
    median = quantile(value, 0.5),
    upper = quantile(value, 0.975)
  )

R_wt_CI <- R_wt %>%
  as.data.frame() %>%
  pivot_longer(-time, names_to = "sim", values_to = "value") %>%
  group_by(time) %>%
  summarise(
    lower = quantile(value, 0.025),
    median = quantile(value, 0.5),
    upper = quantile(value, 0.975)
  )

M1_all_CI <- M1_all %>%
  as.data.frame() %>%
  pivot_longer(-time, names_to = "sim", values_to = "value") %>%
  group_by(time) %>%
  summarise(
    lower = quantile(value, 0.025),
    median = quantile(value, 0.5),
    upper = quantile(value, 0.975)
  )

M2_all_CI <- M2_all %>%
  as.data.frame() %>%
  pivot_longer(-time, names_to = "sim", values_to = "value") %>%
  group_by(time) %>%
  summarise(
    lower = quantile(value, 0.025),
    median = quantile(value, 0.5),
    upper = quantile(value, 0.975)
  )

M_all_CI <- M_all %>%
  as.data.frame() %>%
  pivot_longer(-time, names_to = "sim", values_to = "value") %>%
  group_by(time) %>%
  summarise(
    lower = quantile(value, 0.025),
    median = quantile(value, 0.5),
    upper = quantile(value, 0.975)
  )

S_all_CI <- S_all %>%
  as.data.frame() %>%
  pivot_longer(-time, names_to = "sim", values_to = "value") %>%
  group_by(time) %>%
  summarise(
    lower = quantile(value, 0.025),
    median = quantile(value, 0.5),
    upper = quantile(value, 0.975)
  )

R_all_CI <- R_all %>%
  as.data.frame() %>%
  pivot_longer(-time, names_to = "sim", values_to = "value") %>%
  group_by(time) %>%
  summarise(
    lower = quantile(value, 0.025),
    median = quantile(value, 0.5),
    upper = quantile(value, 0.975)
  )

# Save outputs
write.csv(M1_sp_CI, paste0(path_output, "M1_sp_CI.csv"), row.names = FALSE)
write.csv(M2_sp_CI, paste0(path_output, "M2_sp_CI.csv"), row.names = FALSE)
write.csv(M_sp_CI, paste0(path_output, "M_sp_CI.csv"), row.names = FALSE)
write.csv(S_sp_CI, paste0(path_output, "S_sp_CI.csv"), row.names = FALSE)
write.csv(R_sp_CI, paste0(path_output, "R_sp_CI.csv"), row.names = FALSE)
write.csv(M1_sm_CI, paste0(path_output, "M1_sm_CI.csv"), row.names = FALSE)
write.csv(M2_sm_CI, paste0(path_output, "M2_sm_CI.csv"), row.names = FALSE)
write.csv(M_sm_CI, paste0(path_output, "M_sm_CI.csv"), row.names = FALSE)
write.csv(S_sm_CI, paste0(path_output, "S_sm_CI.csv"), row.names = FALSE)
write.csv(R_sm_CI, paste0(path_output, "R_sm_CI.csv"), row.names = FALSE)
write.csv(M1_au_CI, paste0(path_output, "M1_au_CI.csv"), row.names = FALSE)
write.csv(M2_au_CI, paste0(path_output, "M2_au_CI.csv"), row.names = FALSE)
write.csv(M_au_CI, paste0(path_output, "M_au_CI.csv"), row.names = FALSE)
write.csv(S_au_CI, paste0(path_output, "S_au_CI.csv"), row.names = FALSE)
write.csv(R_au_CI, paste0(path_output, "R_au_CI.csv"), row.names = FALSE)
write.csv(M1_wt_CI, paste0(path_output, "M1_wt_CI.csv"), row.names = FALSE)
write.csv(M2_wt_CI, paste0(path_output, "M2_wt_CI.csv"), row.names = FALSE)
write.csv(M_wt_CI, paste0(path_output, "M_wt_CI.csv"), row.names = FALSE)
write.csv(S_wt_CI, paste0(path_output, "S_wt_CI.csv"), row.names = FALSE)
write.csv(R_wt_CI, paste0(path_output, "R_wt_CI.csv"), row.names = FALSE)
write.csv(M1_all_CI, paste0(path_output, "M1_all_CI.csv"), row.names = FALSE)
write.csv(M2_all_CI, paste0(path_output, "M2_all_CI.csv"), row.names = FALSE)
write.csv(M_all_CI, paste0(path_output, "M_all_CI.csv"), row.names = FALSE)
write.csv(S_all_CI, paste0(path_output, "S_all_CI.csv"), row.names = FALSE)
write.csv(R_all_CI, paste0(path_output, "R_all_CI.csv"), row.names = FALSE)


# --------------------- Disease model -----------------------------------------
total_outcome_intervention_df <- readRDS(paste0(path_output, "total outcome intervention.rds"))
infection_df <- readRDS(paste0(path_output, "infections.rds"))

# Get incidence per 100,000
incidence_hosp <- total_outcome_intervention_df %>%
  # filter(age_bracket == "all") %>%
  merge(births_de %>% select(season, total), by.x = "season_birth", by.y = "season") %>%
  mutate(
    incidence_per_100000 = n_hospitalisations / total * 100000,
    icu_per_100000 = n_icu / total * 100000,
    prevented_hospitalisations_per_100000 = prevented_hospitalisations / total * 100000,
    prevented_icu_per_100000 = prevented_icu / total * 100000
  )

# Get reduction in incidence compared to no immunisation
incidence_reduction_no_imm <- incidence_hosp %>%
  filter(intervention != "No immunisation" & age_bracket == "all") %>%
  merge(
    incidence_hosp %>%
      filter(intervention == "No immunisation" & age_bracket == "all") %>%
      select(season_birth, age_bracket, n_hospitalisations, n_icu, iter) %>%
      rename(
        n_hospitalisations_no_int = n_hospitalisations,
        n_icu_no_int = n_icu
      ),
    by = c("season_birth", "age_bracket", "iter")
  ) %>%
  mutate(
    incidence_reduction_perc = (n_hospitalisations_no_int - n_hospitalisations) / n_hospitalisations_no_int * 100,
    icu_reduction_perc = (n_icu_no_int - n_icu) / n_icu_no_int * 100
  )

# Get reduction in incidence for mAB compared to MV
incidence_reduction_mab_mv <- incidence_hosp %>%
  filter(intervention != "No immunisation" & intervention != "MV" & age_bracket == "all") %>%
  merge(
    incidence_hosp %>%
      filter(intervention == "MV" & age_bracket == "all") %>%
      select(season_birth, age_bracket, incidence_per_100000, icu_per_100000, iter) %>%
      rename(
        incidence_per_100000_mv = incidence_per_100000,
        icu_per_100000_mv = icu_per_100000
      ),
    by = c("season_birth", "age_bracket", "iter")
  ) %>%
  mutate(
    incidence_reduction_perc = (incidence_per_100000_mv - incidence_per_100000) / incidence_per_100000_mv * 100,
    icu_reduction_perc = (icu_per_100000_mv - icu_per_100000) / icu_per_100000_mv * 100
  )

# Compute 95% CI
total_outcome_intervention_int <- total_outcome_intervention_df %>%
  filter(age_bracket == "all") %>%
  group_by(intervention, season_birth, age_bracket) %>%
  summarise(
    hosp_low_95 = quantile(n_hospitalisations, 0.025, na.rm = T),
    hosp_median = quantile(n_hospitalisations, 0.5, na.rm = T),
    hosp_up_95 = quantile(n_hospitalisations, 0.975, na.rm = T),
    prev_low_95 = quantile(prevented_hospitalisations, 0.025, na.rm = T),
    prev_median = quantile(prevented_hospitalisations, 0.5, na.rm = T),
    prev_up_95 = quantile(prevented_hospitalisations, 0.975, na.rm = T),
    icu_low_95 = quantile(n_icu, 0.025, na.rm = T),
    icu_median = quantile(n_icu, 0.5, na.rm = T),
    icu_up_95 = quantile(n_icu, 0.975, na.rm = T),
    icu_prev_low_95 = quantile(prevented_icu, 0.025, na.rm = T),
    icu_prev_median = quantile(prevented_icu, 0.5, na.rm = T),
    icu_prev_up_95 = quantile(prevented_icu, 0.975, na.rm = T)
  ) %>%
  ungroup() %>%
  mutate(season_birth = str_to_sentence(season_birth)) %>%
  mutate(season_birth = factor(season_birth,
                               levels = c(
                                 "Winter", "Spring", "Summer", "Autumn",
                                 "All"
                               )
  ))

incidence_outcome_int <- incidence_hosp %>%
  group_by(intervention, season_birth, age_bracket) %>%
  summarise(
    incidence_up_95 = quantile(incidence_per_100000, 0.025, na.rm = T),
    incidence_median = quantile(incidence_per_100000, 0.5, na.rm = T),
    incidence_low_95 = quantile(incidence_per_100000, 0.975, na.rm = T),
    icu_inc_up_95 = quantile(icu_per_100000, 0.025, na.rm = T),
    icu_inc_median = quantile(icu_per_100000, 0.5, na.rm = T),
    icu_inc_low_95 = quantile(icu_per_100000, 0.975, na.rm = T),
    prev_inc_up_95 = quantile(prevented_hospitalisations_per_100000, 0.025, na.rm = T),
    prev_inc_median = quantile(prevented_hospitalisations_per_100000, 0.5, na.rm = T),
    prev_inc_low_95 = quantile(prevented_hospitalisations_per_100000, 0.975, na.rm = T),
    icu_prev_inc_up_95 = quantile(prevented_icu_per_100000, 0.025, na.rm = T),
    icu_prev_inc_median = quantile(prevented_icu_per_100000, 0.5, na.rm = T),
    icu_prev_inc_low_95 = quantile(prevented_icu_per_100000, 0.975, na.rm = T)
  ) %>%
  ungroup() %>%
  mutate(season_birth = str_to_sentence(season_birth)) %>%
  mutate(season_birth = factor(season_birth,
                               levels = c(
                                 "Winter", "Spring", "Summer", "Autumn",
                                 "All"
                               )
  ))

incidence_reduction_mab_mv_int <- incidence_reduction_mab_mv %>%
  group_by(intervention, season_birth, age_bracket) %>%
  summarise(
    incidence_reduction_low_95 = quantile(incidence_reduction_perc, 0.025, na.rm = T),
    incidence_reduction_median = quantile(incidence_reduction_perc, 0.5, na.rm = T),
    incidence_reduction_up_95 = quantile(incidence_reduction_perc, 0.975, na.rm = T),
    icu_reduction_low_95 = quantile(icu_reduction_perc, 0.025, na.rm = T),
    icu_reduction_median = quantile(icu_reduction_perc, 0.5, na.rm = T),
    icu_reduction_up_95 = quantile(icu_reduction_perc, 0.975, na.rm = T)
  ) %>%
  ungroup() %>%
  mutate(season_birth = str_to_sentence(season_birth)) %>%
  mutate(season_birth = factor(season_birth,
                               levels = c(
                                 "Winter", "Spring", "Summer", "Autumn",
                                 "All"
                               )
  ))


incidence_reduction_no_imm_int <- incidence_reduction_no_imm %>%
  group_by(intervention, season_birth, age_bracket) %>%
  summarise(
    incidence_reduction_low_95 = quantile(incidence_reduction_perc, 0.025, na.rm = T),
    incidence_reduction_median = quantile(incidence_reduction_perc, 0.5, na.rm = T),
    incidence_reduction_up_95 = quantile(incidence_reduction_perc, 0.975, na.rm = T),
    icu_inc_reduction_low_95 = quantile(icu_reduction_perc, 0.025, na.rm = T),
    icu_inc_reduction_median = quantile(icu_reduction_perc, 0.5, na.rm = T),
    icu_inc_reduction_up_95 = quantile(icu_reduction_perc, 0.975, na.rm = T)
  ) %>%
  ungroup() %>%
  mutate(season_birth = str_to_sentence(season_birth)) %>%
  mutate(season_birth = factor(season_birth,
                               levels = c(
                                 "Winter", "Spring", "Summer", "Autumn",
                                 "All"
                               )
  ))

# Get absolute reduction in incidence for mAB compared to MV
abs_inc_reduction <- incidence_hosp %>%
  filter(intervention == "MV" & age_bracket == "all" & season_birth == "spring") %>%
  select(season_birth, age_bracket, incidence_per_100000, icu_per_100000, iter) %>%
  merge(
    incidence_hosp %>%
      filter(intervention == "+ mAB for spring births") %>%
      select(season_birth, age_bracket, incidence_per_100000, icu_per_100000, iter) %>%
      rename(
        incidence_mAB = incidence_per_100000,
        icu_inc_mAB = icu_per_100000
      ),
    by = c("season_birth", "age_bracket", "iter")
  ) %>%
  mutate(
    abs_incidence_reduction = incidence_per_100000 - incidence_mAB,
    abs_icu_incidence_reduction = icu_per_100000 - icu_inc_mAB,
    season_birth = "spring"
  ) %>%
  rbind(
    incidence_hosp %>%
      filter(intervention == "MV" & age_bracket == "all" & season_birth == "summer") %>%
      select(season_birth, age_bracket, incidence_per_100000, icu_per_100000, iter) %>%
      merge(
        incidence_hosp %>%
          filter(intervention == "+ mAB for summer births") %>%
          select(season_birth, age_bracket, incidence_per_100000, icu_per_100000, iter) %>%
          rename(
            incidence_mAB = incidence_per_100000,
            icu_inc_mAB = icu_per_100000
          ),
        by = c("season_birth", "age_bracket", "iter")
      ) %>%
      mutate(
        abs_incidence_reduction = incidence_per_100000 - incidence_mAB,
        abs_icu_incidence_reduction = icu_per_100000 - icu_inc_mAB,
        season_birth = "summer"
      )
  ) %>%
  rbind(
    incidence_hosp %>%
      filter(intervention == "MV" & age_bracket == "all" & season_birth == "autumn") %>%
      select(season_birth, age_bracket, incidence_per_100000, icu_per_100000, iter) %>%
      merge(
        incidence_hosp %>%
          filter(intervention == "+ mAB for autumn births") %>%
          select(season_birth, age_bracket, incidence_per_100000, icu_per_100000, iter) %>%
          rename(
            incidence_mAB = incidence_per_100000,
            icu_inc_mAB = icu_per_100000
          ),
        by = c("season_birth", "age_bracket", "iter")
      ) %>%
      mutate(
        abs_incidence_reduction = incidence_per_100000 - incidence_mAB,
        abs_icu_incidence_reduction = icu_per_100000 - icu_inc_mAB,
        season_birth = "autumn"
      )
  ) %>%
  rbind(
    incidence_hosp %>%
      filter(intervention == "MV" & age_bracket == "all" & season_birth == "winter") %>%
      select(season_birth, age_bracket, incidence_per_100000, icu_per_100000, iter) %>%
      merge(
        incidence_hosp %>%
          filter(intervention == "+ mAB for winter births") %>%
          select(season_birth, age_bracket, incidence_per_100000, icu_per_100000, iter) %>%
          rename(
            incidence_mAB = incidence_per_100000,
            icu_inc_mAB = icu_per_100000
          ),
        by = c("season_birth", "age_bracket", "iter")
      ) %>%
      mutate(
        abs_incidence_reduction = incidence_per_100000 - incidence_mAB,
        abs_icu_incidence_reduction = icu_per_100000 - icu_inc_mAB,
        season_birth = "winter"
      )
  ) %>%
  rbind(
    incidence_hosp %>%
      filter(intervention == "MV" & age_bracket == "all" & season_birth == "all") %>%
      select(season_birth, age_bracket, incidence_per_100000, icu_per_100000, iter) %>%
      merge(
        incidence_hosp %>%
          filter(intervention == "+ mAB for all") %>%
          select(season_birth, age_bracket, incidence_per_100000, icu_per_100000, iter) %>%
          rename(
            incidence_mAB = incidence_per_100000,
            icu_inc_mAB = icu_per_100000
          ),
        by = c("season_birth", "age_bracket", "iter")
      ) %>%
      mutate(
        abs_incidence_reduction = incidence_per_100000 - incidence_mAB,
        abs_icu_incidence_reduction = icu_per_100000 - icu_inc_mAB,
        season_birth = "all"
      )
  )

abs_inc_reduction_ci <- abs_inc_reduction %>%
  group_by(season_birth, age_bracket) %>%
  summarise(
    abs_incidence_reduction_low_95 = quantile(abs_incidence_reduction, 0.025, na.rm = T),
    abs_incidence_reduction_median = quantile(abs_incidence_reduction, 0.5, na.rm = T),
    abs_incidence_reduction_up_95 = quantile(abs_incidence_reduction, 0.975, na.rm = T),
    abs_icu_incidence_reduction_low_95 = quantile(abs_icu_incidence_reduction, 0.025, na.rm = T),
    abs_icu_incidence_reduction_median = quantile(abs_icu_incidence_reduction, 0.5, na.rm = T),
    abs_icu_incidence_reduction_up_95 = quantile(abs_icu_incidence_reduction, 0.975, na.rm = T)
  ) %>%
  ungroup() %>%
  mutate(season_birth = str_to_sentence(season_birth)) %>%
  mutate(season_birth = factor(season_birth,
                               levels = c(
                                 "Winter", "Spring", "Summer", "Autumn",
                                 "All"
                               )
  ))


# We need to recalculate the CI for No immunisation
total_outcome_intervention_plt <- total_outcome_intervention_df %>%
  filter(!(intervention == "No immunisation" & season_birth == "all")) %>%
  rbind(
    total_outcome_intervention_df %>%
      filter(intervention == "No immunisation" & season_birth != "all" & age_bracket != "all") %>%
      group_by(intervention, season_birth, iter) %>%
      summarise(
        age_bracket = "all",
        n_hospitalisations = sum(n_hospitalisations),
        prevented_hospitalisations = 0,
        n_icu = sum(n_icu),
        prevented_icu = 0
      ) %>%
      ungroup()
  )

total_outcome_intervention_plt <- total_outcome_intervention_plt %>%
  select(
    intervention, season_birth, age_bracket, n_hospitalisations, n_icu,
    prevented_hospitalisations, prevented_icu,
    iter
  ) %>%
  rbind(
    total_outcome_intervention_plt %>%
      filter(intervention == "No immunisation" & age_bracket != "all") %>%
      group_by(intervention, iter) %>%
      summarise(
        season_birth = "all",
        age_bracket = "all",
        n_hospitalisations = sum(n_hospitalisations),
        prevented_hospitalisations = sum(prevented_hospitalisations),
        n_icu = sum(n_icu),
        prevented_icu = sum(prevented_icu)
      ) %>%
      ungroup()
  )

# Compute 95% CI
total_outcome_intervention_plt_int <- total_outcome_intervention_plt %>%
  filter(age_bracket == "all") %>%
  group_by(intervention, season_birth, age_bracket) %>%
  summarise(
    hosp_low_95 = quantile(n_hospitalisations, 0.025, na.rm = T),
    hosp_median = quantile(n_hospitalisations, 0.5, na.rm = T),
    hosp_up_95 = quantile(n_hospitalisations, 0.975, na.rm = T),
    prev_low_95 = quantile(prevented_hospitalisations, 0.025, na.rm = T),
    prev_median = quantile(prevented_hospitalisations, 0.5, na.rm = T),
    prev_up_95 = quantile(prevented_hospitalisations, 0.975, na.rm = T),
    icu_low_95 = quantile(n_icu, 0.025, na.rm = T),
    icu_median = quantile(n_icu, 0.5, na.rm = T),
    icu_up_95 = quantile(n_icu, 0.975, na.rm = T),
    prev_icu_low_95 = quantile(prevented_icu, 0.025, na.rm = T),
    prev_icu_median = quantile(prevented_icu, 0.5, na.rm = T),
    prev_icu_up_95 = quantile(prevented_icu, 0.975, na.rm = T)
  ) %>%
  ungroup() %>%
  mutate(season_birth = str_to_sentence(season_birth)) %>%
  mutate(season_birth = factor(season_birth,
                               levels = c(
                                 "Winter", "Spring", "Summer", "Autumn",
                                 "All"
                               )
  ))

# Total number of hospitalisations and ICU admissions by age
total_hosp_no_int_age <- total_outcome_intervention_df %>%
  filter(age_bracket != "all" & intervention == "No immunisation") %>%
  group_by(age_bracket, iter, intervention) %>%
  summarise(
    season_birth = "all",
    n_hospitalisations = sum(n_hospitalisations)
  ) %>%
  rbind(total_outcome_intervention_df %>%
          filter(age_bracket != "all" & intervention == "No immunisation") %>%
          group_by(age_bracket, iter, season_birth, intervention) %>%
          summarise(n_hospitalisations = sum(n_hospitalisations, na.rm = T))) %>%
  ungroup() %>%
  group_by(age_bracket, intervention, season_birth) %>%
  summarise(
    hosp_low_95 = quantile(n_hospitalisations, 0.025, na.rm = T),
    hosp_median = quantile(n_hospitalisations, 0.5, na.rm = T),
    hosp_up_95 = quantile(n_hospitalisations, 0.975, na.rm = T)
  ) %>%
  ungroup() %>%
  mutate(
    season_birth = str_to_sentence(season_birth),
    season_birth = factor(season_birth,
                          levels = c(
                            "Winter", "Spring", "Summer", "Autumn",
                            "All"
                          )
    ),
    age_bracket = factor(age_bracket,
                         levels = c(
                           "0", "1", "2", "3", "4", "5", "6", "7", "8",
                           "9", "10", "11", "12"
                         )
    )
  )

total_icu_no_int_age <- total_outcome_intervention_df %>%
  filter(age_bracket != "all" & intervention == "No immunisation") %>%
  group_by(age_bracket, iter, intervention) %>%
  summarise(
    season_birth = "all",
    n_icu = sum(n_icu)
  ) %>%
  rbind(total_outcome_intervention_df %>%
          filter(age_bracket != "all" & intervention == "No immunisation") %>%
          group_by(age_bracket, iter, season_birth, intervention) %>%
          summarise(n_icu = sum(n_icu, na.rm = T))) %>%
  ungroup() %>%
  group_by(age_bracket, intervention, season_birth) %>%
  summarise(
    icu_low_95 = quantile(n_icu, 0.025, na.rm = T),
    icu_median = quantile(n_icu, 0.5, na.rm = T),
    icu_up_95 = quantile(n_icu, 0.975, na.rm = T)
  ) %>%
  ungroup() %>%
  mutate(
    season_birth = str_to_sentence(season_birth),
    season_birth = factor(season_birth,
                          levels = c(
                            "Winter", "Spring", "Summer", "Autumn",
                            "All"
                          )
    ),
    age_bracket = factor(age_bracket,
                         levels = c(
                           "0", "1", "2", "3", "4", "5", "6", "7", "8",
                           "9", "10", "11", "12"
                         )
    )
  )

# Look at total burden by season of birth
total_hosp_1_y <- total_outcome_intervention_df %>%
  filter(age_bracket == "all") %>%
  group_by(season_birth, intervention) %>%
  summarise(
    hosp_low_95 = quantile(n_hospitalisations, 0.025),
    hosp_median = quantile(n_hospitalisations, 0.5),
    hosp_up_95 = quantile(n_hospitalisations, 0.975)
  ) %>%
  ungroup() %>%
  mutate(season_birth = str_to_sentence(season_birth))

total_icu_1_y <- total_outcome_intervention_df %>%
  filter(age_bracket == "all") %>%
  group_by(season_birth, intervention) %>%
  summarise(
    icu_low_95 = quantile(n_icu, 0.025),
    icu_median = quantile(n_icu, 0.5),
    icu_up_95 = quantile(n_icu, 0.975)
  ) %>%
  ungroup() %>%
  mutate(season_birth = str_to_sentence(season_birth))

# Save files
write.csv(total_outcome_intervention_int, paste0(path_output, "Final outputs hosp and icu.csv"), row.names = F)
write.csv(incidence_outcome_int, paste0(path_output, "Final outputs hosp incidence.csv"), row.names = F)
write.csv(incidence_reduction_no_imm_int, paste0(path_output, "Final outputs incidence vs no.csv"), row.names = F)
write.csv(incidence_reduction_mab_mv_int, paste0(path_output, "Final outputs incidence vs MV.csv"), row.names = F)
write.csv(abs_inc_reduction_ci, paste0(path_output, "Final outputs absolute incidence reduction.csv"), row.names = F)


# --------------------- Reduction in hospitalisations and ICU visits -----------
outcome_reduction_perc_MV <- total_outcome_intervention_df %>%
  filter(intervention == "No immunisation") %>%
  group_by(season_birth, age_bracket, iter) %>%
  summarise(
    n_hosp_none = sum(n_hospitalisations),
    n_icu_non = sum(n_icu)
  ) %>%
  merge(
    total_outcome_intervention_df %>%
      filter(intervention == "MV") %>%
      group_by(season_birth, age_bracket, iter) %>%
      summarise(
        n_hosp_MV = n_hospitalisations,
        hosp_prevented_MV = prevented_hospitalisations,
        n_icu_MV = n_icu
      ),
    by = c("season_birth", "age_bracket", "iter")
  ) %>%
  mutate(
    perc_reduction_hosp_MV = (n_hosp_none - n_hosp_MV) / n_hosp_none * 100,
    perc_reduction_icu_MV = (n_icu_non - n_icu_MV) / n_icu_non * 100
  )

outcome_reduction_perc_spring <- total_outcome_intervention_df %>%
  filter(intervention == "No immunisation") %>%
  group_by(season_birth, age_bracket, iter) %>%
  summarise(
    n_hosp_none = sum(n_hospitalisations),
    n_icu_non = sum(n_icu)
  ) %>%
  merge(
    total_outcome_intervention_df %>%
      filter(season_birth == "spring" & intervention == "+ mAB for spring births") %>%
      group_by(season_birth, age_bracket, iter) %>%
      summarise(
        n_hosp_spring = n_hospitalisations,
        hosp_prevented_spring = prevented_hospitalisations,
        n_icu_spring = n_icu
      ),
    by = c("season_birth", "age_bracket", "iter")
  ) %>%
  mutate(
    perc_reduction_hosp_spring = (n_hosp_none - n_hosp_spring) / n_hosp_none * 100,
    perc_reduction_icu_spring = (n_icu_non - n_icu_spring) / n_icu_non * 100
  )

outcome_reduction_perc_summer <- total_outcome_intervention_df %>%
  filter(intervention == "No immunisation") %>%
  group_by(season_birth, age_bracket, iter) %>%
  summarise(
    n_hosp_none = sum(n_hospitalisations),
    n_icu_non = sum(n_icu)
  ) %>%
  merge(
    total_outcome_intervention_df %>%
      filter(season_birth == "summer" & intervention == "+ mAB for summer births") %>%
      group_by(season_birth, age_bracket, iter) %>%
      summarise(
        n_hosp_summer = n_hospitalisations,
        hosp_prevented_summer = prevented_hospitalisations,
        n_icu_summer = n_icu
      ),
    by = c("season_birth", "age_bracket", "iter")
  ) %>%
  mutate(
    perc_reduction_hosp_summer = (n_hosp_none - n_hosp_summer) / n_hosp_none * 100,
    perc_reduction_icu_summer = (n_icu_non - n_icu_summer) / n_icu_non * 100
  )

outcome_reduction_perc_autumn <- total_outcome_intervention_df %>%
  filter(intervention == "No immunisation") %>%
  group_by(season_birth, age_bracket, iter) %>%
  summarise(
    n_hosp_none = sum(n_hospitalisations),
    n_icu_non = sum(n_icu)
  ) %>%
  merge(
    total_outcome_intervention_df %>%
      filter(season_birth == "autumn" & intervention == "+ mAB for autumn births") %>%
      group_by(season_birth, age_bracket, iter) %>%
      summarise(
        n_hosp_autumn = n_hospitalisations,
        hosp_prevented_autumn = prevented_hospitalisations,
        n_icu_autumn = n_icu
      ),
    by = c("season_birth", "age_bracket", "iter")
  ) %>%
  mutate(
    perc_reduction_hosp_autumn = (n_hosp_none - n_hosp_autumn) / n_hosp_none * 100,
    perc_reduction_icu_autumn = (n_icu_non - n_icu_autumn) / n_icu_non * 100
  )

outcome_reduction_perc_winter <- total_outcome_intervention_df %>%
  filter(intervention == "No immunisation") %>%
  group_by(season_birth, age_bracket, iter) %>%
  summarise(
    n_hosp_none = sum(n_hospitalisations),
    n_icu_non = sum(n_icu)
  ) %>%
  merge(
    total_outcome_intervention_df %>%
      filter(season_birth == "winter" & intervention == "+ mAB for winter births") %>%
      group_by(season_birth, age_bracket, iter) %>%
      summarise(
        n_hosp_winter = n_hospitalisations,
        hosp_prevented_winter = prevented_hospitalisations,
        n_icu_winter = n_icu
      ),
    by = c("season_birth", "age_bracket", "iter")
  ) %>%
  mutate(
    perc_reduction_hosp_winter = (n_hosp_none - n_hosp_winter) / n_hosp_none * 100,
    perc_reduction_icu_winter = (n_icu_non - n_icu_winter) / n_icu_non * 100
  )

outcome_reduction_perc <- outcome_reduction_perc_MV %>%
  select(
    season_birth, age_bracket, iter,
    perc_reduction_hosp_MV, perc_reduction_icu_MV
  ) %>%
  merge(
    outcome_reduction_perc_spring %>%
      select(
        season_birth, age_bracket, iter,
        perc_reduction_hosp_spring, perc_reduction_icu_spring
      ),
    by = c("season_birth", "age_bracket", "iter")
  ) %>%
  merge(
    outcome_reduction_perc_summer %>%
      select(
        season_birth, age_bracket, iter,
        perc_reduction_hosp_summer, perc_reduction_icu_summer
      ),
    by = c("season_birth", "age_bracket", "iter")
  ) %>%
  merge(
    outcome_reduction_perc_autumn %>%
      select(
        season_birth, age_bracket, iter,
        perc_reduction_hosp_autumn, perc_reduction_icu_autumn
      ),
    by = c("season_birth", "age_bracket", "iter")
  ) %>%
  merge(
    outcome_reduction_perc_winter %>%
      select(
        season_birth, age_bracket, iter,
        perc_reduction_hosp_winter, perc_reduction_icu_winter
      ),
    by = c("season_birth", "age_bracket", "iter")
  )

# Get 95% CIs for reduction in hospitalisations and ICU visits
outcome_reduction_perc_CI <- outcome_reduction_perc %>%
  group_by(season_birth, age_bracket) %>%
  summarise(
    hosp_reduction_MV_lower = quantile(perc_reduction_hosp_MV, 0.025, na.rm = TRUE),
    hosp_reduction_MV_median = quantile(perc_reduction_hosp_MV, 0.5, na.rm = TRUE),
    hosp_reduction_MV_upper = quantile(perc_reduction_hosp_MV, 0.975, na.rm = TRUE),
    hosp_reduction_spring_lower = quantile(perc_reduction_hosp_spring, 0.025, na.rm = TRUE),
    hosp_reduction_spring_median = quantile(perc_reduction_hosp_spring, 0.5, na.rm = TRUE),
    hosp_reduction_spring_upper = quantile(perc_reduction_hosp_spring, 0.975, na.rm = TRUE),
    hosp_reduction_summer_lower = quantile(perc_reduction_hosp_summer, 0.025, na.rm = TRUE),
    hosp_reduction_summer_median = quantile(perc_reduction_hosp_summer, 0.5, na.rm = TRUE),
    hosp_reduction_summer_upper = quantile(perc_reduction_hosp_summer, 0.975, na.rm = TRUE),
    hosp_reduction_autumn_lower = quantile(perc_reduction_hosp_autumn, 0.025, na.rm = TRUE),
    hosp_reduction_autumn_median = quantile(perc_reduction_hosp_autumn, 0.5, na.rm = TRUE),
    hosp_reduction_autumn_upper = quantile(perc_reduction_hosp_autumn, 0.975, na.rm = TRUE),
    hosp_reduction_winter_lower = quantile(perc_reduction_hosp_winter, 0.025, na.rm = TRUE),
    hosp_reduction_winter_median = quantile(perc_reduction_hosp_winter, 0.5, na.rm = TRUE),
    hosp_reduction_winter_upper = quantile(perc_reduction_hosp_winter, 0.975, na.rm = TRUE),
    icu_reduction_MV_lower = quantile(perc_reduction_icu_MV, 0.025, na.rm = TRUE),
    icu_reduction_MV_median = quantile(perc_reduction_icu_MV, 0.5, na.rm = TRUE),
    icu_reduction_MV_upper = quantile(perc_reduction_icu_MV, 0.975, na.rm = TRUE),
    icu_reduction_spring_lower = quantile(perc_reduction_icu_spring, 0.025, na.rm = TRUE),
    icu_reduction_spring_median = quantile(perc_reduction_icu_spring, 0.5, na.rm = TRUE),
    icu_reduction_spring_upper = quantile(perc_reduction_icu_spring, 0.975, na.rm = TRUE),
    icu_reduction_summer_lower = quantile(perc_reduction_icu_summer, 0.025, na.rm = TRUE),
    icu_reduction_summer_median = quantile(perc_reduction_icu_summer, 0.5, na.rm = TRUE),
    icu_reduction_summer_upper = quantile(perc_reduction_icu_summer, 0.975, na.rm = TRUE),
    icu_reduction_autumn_lower = quantile(perc_reduction_icu_autumn, 0.025, na.rm = TRUE),
    icu_reduction_autumn_median = quantile(perc_reduction_icu_autumn, 0.5, na.rm = TRUE),
    icu_reduction_autumn_upper = quantile(perc_reduction_icu_autumn, 0.975, na.rm = TRUE),
    icu_reduction_winter_lower = quantile(perc_reduction_icu_winter, 0.025, na.rm = TRUE),
    icu_reduction_winter_median = quantile(perc_reduction_icu_winter, 0.5, na.rm = TRUE),
    icu_reduction_winter_upper = quantile(perc_reduction_icu_winter, 0.975, na.rm = TRUE)
  )

# Reduction by season of birth
outcome_reduction_season <- total_outcome_intervention_df %>%
  filter(intervention == "MV" & season_birth == "all") %>%
  mutate(
    prev_total = prevented_hospitalisations,
    prev_total_icu = prevented_icu
  ) %>%
  select(age_bracket, iter, prev_total, prev_total_icu) %>%
  merge(
    total_outcome_intervention_df %>%
      filter(intervention == "MV" & season_birth == "spring") %>%
      mutate(
        prev_spring = prevented_hospitalisations,
        prev_spring_icu = prevented_icu
      ) %>%
      select(age_bracket, iter, prev_spring, prev_spring_icu),
    by = c("age_bracket", "iter")
  ) %>%
  merge(
    total_outcome_intervention_df %>%
      filter(intervention == "MV" & season_birth == "summer") %>%
      mutate(
        prev_summer = prevented_hospitalisations,
        prev_summer_icu = prevented_icu
      ) %>%
      select(age_bracket, iter, prev_summer, prev_summer_icu),
    by = c("age_bracket", "iter")
  ) %>%
  merge(
    total_outcome_intervention_df %>%
      filter(intervention == "MV" & season_birth == "autumn") %>%
      mutate(
        prev_autumn = prevented_hospitalisations,
        prev_autumn_icu = prevented_icu
      ) %>%
      select(age_bracket, iter, prev_autumn, prev_autumn_icu),
    by = c("age_bracket", "iter")
  ) %>%
  merge(
    total_outcome_intervention_df %>%
      filter(intervention == "MV" & season_birth == "winter") %>%
      mutate(
        prev_winter = prevented_hospitalisations,
        prev_winter_icu = prevented_icu
      ) %>%
      select(age_bracket, iter, prev_winter, prev_winter_icu),
    by = c("age_bracket", "iter")
  ) %>%
  mutate(
    perc_reduction_spring = prev_spring / prev_total * 100,
    perc_reduction_summer = prev_summer / prev_total * 100,
    perc_reduction_autumn = prev_autumn / prev_total * 100,
    perc_reduction_winter = prev_winter / prev_total * 100,
    perc_reduction_spring_icu = prev_spring_icu / prev_total_icu * 100,
    perc_reduction_summer_icu = prev_summer_icu / prev_total_icu * 100,
    perc_reduction_autumn_icu = prev_autumn_icu / prev_total_icu * 100,
    perc_reduction_winter_icu = prev_winter_icu / prev_total_icu * 100
  ) %>%
  group_by(age_bracket) %>%
  summarise(
    hosp_reduction_spring_lower = quantile(perc_reduction_spring, 0.025, na.rm = TRUE),
    hosp_reduction_spring_median = quantile(perc_reduction_spring, 0.5, na.rm = TRUE),
    hosp_reduction_spring_upper = quantile(perc_reduction_spring, 0.975, na.rm = TRUE),
    hosp_reduction_summer_lower = quantile(perc_reduction_summer, 0.025, na.rm = TRUE),
    hosp_reduction_summer_median = quantile(perc_reduction_summer, 0.5, na.rm = TRUE),
    hosp_reduction_summer_upper = quantile(perc_reduction_summer, 0.975, na.rm = TRUE),
    hosp_reduction_autumn_lower = quantile(perc_reduction_autumn, 0.025, na.rm = TRUE),
    hosp_reduction_autumn_median = quantile(perc_reduction_autumn, 0.5, na.rm = TRUE),
    hosp_reduction_autumn_upper = quantile(perc_reduction_autumn, 0.975, na.rm = TRUE),
    hosp_reduction_winter_lower = quantile(perc_reduction_winter, 0.025, na.rm = TRUE),
    hosp_reduction_winter_median = quantile(perc_reduction_winter, 0.5, na.rm = TRUE),
    hosp_reduction_winter_upper = quantile(perc_reduction_winter, 0.975, na.rm = TRUE),
    icu_reduction_spring_lower = quantile(perc_reduction_spring_icu, 0.025, na.rm = TRUE),
    icu_reduction_spring_median = quantile(perc_reduction_spring_icu, 0.5, na.rm = TRUE),
    icu_reduction_spring_upper = quantile(perc_reduction_spring_icu, 0.975, na.rm = TRUE),
    icu_reduction_summer_lower = quantile(perc_reduction_summer_icu, 0.025, na.rm = TRUE),
    icu_reduction_summer_median = quantile(perc_reduction_summer_icu, 0.5, na.rm = TRUE),
    icu_reduction_summer_upper = quantile(perc_reduction_summer_icu, 0.975, na.rm = TRUE),
    icu_reduction_autumn_lower = quantile(perc_reduction_autumn_icu, 0.025, na.rm = TRUE),
    icu_reduction_autumn_median = quantile(perc_reduction_autumn_icu, 0.5, na.rm = TRUE),
    icu_reduction_autumn_upper = quantile(perc_reduction_autumn_icu, 0.975, na.rm = TRUE),
    icu_reduction_winter_lower = quantile(perc_reduction_winter_icu, 0.025, na.rm = TRUE),
    icu_reduction_winter_median = quantile(perc_reduction_winter_icu, 0.5, na.rm = TRUE),
    icu_reduction_winter_upper = quantile(perc_reduction_winter_icu, 0.975, na.rm = TRUE)
  )

# Save outputs
write.csv(outcome_reduction_perc_CI, paste0(path_output, "outcome_reduction_perc_CI.csv"), row.names = FALSE)
write.csv(outcome_reduction_season, paste0(path_output, "outcome_reduction_season.csv"), row.names = FALSE)


# --------------- Incidence vs winter -----------------------
incidence_outcome_wint <- incidence_hosp %>%
  filter(season_birth == "winter" & age_bracket == "all") %>%
  select(
    intervention, age_bracket, incidence_per_100000, icu_per_100000,
    iter
  ) %>%
  rename(
    hosp_winter = incidence_per_100000,
    icu_winter = icu_per_100000
  ) %>%
  merge(
    incidence_hosp %>%
      filter(season_birth == "spring" & age_bracket == "all") %>%
      select(
        intervention, age_bracket, incidence_per_100000, icu_per_100000,
        iter
      ) %>%
      rename(
        hosp_spring = incidence_per_100000,
        icu_spring = icu_per_100000
      ),
    by = c("intervention", "age_bracket", "iter")
  ) %>%
  merge(
    incidence_hosp %>%
      filter(season_birth == "summer" & age_bracket == "all") %>%
      select(
        intervention, age_bracket, incidence_per_100000, icu_per_100000,
        iter
      ) %>%
      rename(
        hosp_summer = incidence_per_100000,
        icu_summer = icu_per_100000
      ),
    by = c("intervention", "age_bracket", "iter")
  ) %>%
  merge(
    incidence_hosp %>%
      filter(season_birth == "autumn" & age_bracket == "all") %>%
      select(
        intervention, age_bracket, incidence_per_100000, icu_per_100000,
        iter
      ) %>%
      rename(
        hosp_autumn = incidence_per_100000,
        icu_autumn = icu_per_100000
      ),
    by = c("intervention", "age_bracket", "iter")
  ) %>%
  mutate(
    perc_increase_hosp_spring = (hosp_spring - hosp_winter) / hosp_winter * 100,
    perc_increase_hosp_summer = (hosp_summer - hosp_winter) / hosp_winter * 100,
    perc_increase_hosp_autumn = (hosp_autumn - hosp_winter) / hosp_winter * 100,
    perc_increase_icu_spring = (icu_spring - icu_winter) / icu_winter * 100,
    perc_increase_icu_summer = (icu_summer - icu_winter) / icu_winter * 100,
    perc_increase_icu_autumn = (icu_autumn - icu_winter) / icu_winter * 100
  )

# Get 95% CIs
incidence_outcome_wint_CI <- incidence_outcome_wint %>%
  group_by(intervention, age_bracket) %>%
  summarise(
    hosp_increase_spring_lower = quantile(perc_increase_hosp_spring, 0.025, na.rm = TRUE),
    hosp_increase_spring_median = quantile(perc_increase_hosp_spring, 0.5, na.rm = TRUE),
    hosp_increase_spring_upper = quantile(perc_increase_hosp_spring, 0.975, na.rm = TRUE),
    hosp_increase_summer_lower = quantile(perc_increase_hosp_summer, 0.025, na.rm = TRUE),
    hosp_increase_summer_median = quantile(perc_increase_hosp_summer, 0.5, na.rm = TRUE),
    hosp_increase_summer_upper = quantile(perc_increase_hosp_summer, 0.975, na.rm = TRUE),
    hosp_increase_autumn_lower = quantile(perc_increase_hosp_autumn, 0.025, na.rm = TRUE),
    hosp_increase_autumn_median = quantile(perc_increase_hosp_autumn, 0.5, na.rm = TRUE),
    hosp_increase_autumn_upper = quantile(perc_increase_hosp_autumn, 0.975, na.rm = TRUE),
    icu_increase_spring_lower = quantile(perc_increase_icu_spring, 0.025, na.rm = TRUE),
    icu_increase_spring_median = quantile(perc_increase_icu_spring, 0.5, na.rm = TRUE),
    icu_increase_spring_upper = quantile(perc_increase_icu_spring, 0.975, na.rm = TRUE),
    icu_increase_summer_lower = quantile(perc_increase_icu_summer, 0.025, na.rm = TRUE),
    icu_increase_summer_median = quantile(perc_increase_icu_summer, 0.5, na.rm = TRUE),
    icu_increase_summer_upper = quantile(perc_increase_icu_summer, 0.975, na.rm = TRUE),
    icu_increase_autumn_lower = quantile(perc_increase_icu_autumn, 0.025, na.rm = TRUE),
    icu_increase_autumn_median = quantile(perc_increase_icu_autumn, 0.5, na.rm = TRUE),
    icu_increase_autumn_upper = quantile(perc_increase_icu_autumn, 0.975, na.rm = TRUE)
  )

# -------------------- Average age at infection --------------------------------
# Add total infections by season of birth
infection_age <- infection_df %>%
  mutate(age_bracket = trunc(age_midpoint / 30)) %>%
  group_by(season_birth, age_bracket, iter) %>%
  summarise(n_infections = sum(n_infections)) %>%
  select(season_birth, age_bracket, iter, n_infections)

infection_age_tot <- infection_age %>%
  mutate(age_bracket = as.character(age_bracket)) %>%
  rbind(
    infection_age %>%
      group_by(season_birth, iter) %>%
      summarise(
        age_bracket = "all",
        n_infections = sum(n_infections)
      )
  )

# Calculate average age at infection
avg_age_infection <- infection_age %>%
  ungroup() %>%
  group_by(iter, season_birth) %>%
  summarise(av_age = sum(age_bracket * n_infections) / sum(n_infections)) %>%
  ungroup() %>%
  group_by(season_birth) %>%
  summarise(
    av_age_low_95 = quantile(av_age, 0.025),
    av_age_median = quantile(av_age, 0.5),
    av_age_up_95 = quantile(av_age, 0.975)
  )

# Save outputs
write.csv(avg_age_infection, paste0(path_output, "avg_age_infection.csv"), row.names = FALSE)

# -------------------- Infection incidence rate --------------------------------
incidence_infection <- infection_age_tot %>%
  merge(
    births_de %>% select(season, total),
    by.x = "season_birth", by.y = "season"
  ) %>%
  mutate(incidence_per_100000 = n_infections / total * 100000)

# Get 95% CIs
incidence_infection_CI <- incidence_infection %>%
  group_by(season_birth, age_bracket) %>%
  summarise(
    incidence_lower = quantile(incidence_per_100000, 0.025),
    incidence_median = quantile(incidence_per_100000, 0.5),
    incidence_upper = quantile(incidence_per_100000, 0.975)
  )

# Save outputs
write.csv(incidence_infection_CI, paste0(path_output, "incidence_infection.csv"), row.names = FALSE)

# ---------------------------- Average age with mAB ----------------------------
avg_age_mAB <- total_outcome_age_mAB_df %>%
  ungroup() %>%
  group_by(iter, season_birth, intervention) %>%
  summarise(
    av_age_hosp = sum(age_months * hosp) / sum(hosp),
    av_age_icu = sum(age_months * icu) / sum(icu)
  ) %>%
  ungroup() %>%
  group_by(season_birth, intervention) %>%
  summarise(
    av_age_hosp_low_95 = quantile(av_age_hosp, 0.025),
    av_age_hosp_median = quantile(av_age_hosp, 0.5),
    av_age_hosp_up_95 = quantile(av_age_hosp, 0.975),
    av_age_icu_low_95 = quantile(av_age_icu, 0.025),
    av_age_icu_median = quantile(av_age_icu, 0.5),
    av_age_icu_up_95 = quantile(av_age_icu, 0.975)
  )

# Save outputs
write.csv(avg_age_mAB, paste0(path_output, "avg_age_mAB.csv"), row.names = FALSE)

# Get average age at hospitalisation
av_age_hosp <- total_outcome_intervention_df %>%
  ungroup() %>%
  filter(age_bracket != "all" & age_bracket != "12-14") %>%
  mutate(age_bracket = case_when(
    age_bracket == "<1" ~ 0,
    T ~ as.numeric(age_bracket)
  )) %>%
  rbind(
    total_outcome_intervention_df %>%
      ungroup() %>%
      filter(age_bracket != "all" & age_bracket != "12-14" & intervention == "MV") %>%
      mutate(age_bracket = case_when(
        age_bracket == "<1" ~ 0,
        T ~ as.numeric(age_bracket)
      )) %>%
      group_by(iter, intervention, age_bracket) %>%
      summarise(
        season_birth = "all",
        n_hospitalisations = sum(n_hospitalisations),
        prevented_hospitalisations = sum(prevented_hospitalisations),
        prevented_icu = sum(prevented_icu),
        n_icu = sum(n_icu)
      ) %>%
      ungroup()
  ) %>%
  group_by(iter, intervention, season_birth) %>%
  summarise(av_age = sum(age_bracket * n_hospitalisations) / sum(n_hospitalisations)) %>%
  ungroup() %>%
  group_by(intervention, season_birth) %>%
  summarise(
    av_age_low_95 = quantile(av_age, 0.025),
    av_age_median = quantile(av_age, 0.5),
    av_age_up_95 = quantile(av_age, 0.975)
  ) %>%
  ungroup()

av_age_hosp %>% write.csv(paste0(path_output, "Avg age at hospitalisation.csv"), row.names = F)

# Get average age at ICU admission
av_age_icu <- total_outcome_intervention_df %>%
  ungroup() %>%
  filter(age_bracket != "all" & age_bracket != "12-14") %>%
  mutate(age_bracket = case_when(
    age_bracket == "<1" ~ 0,
    T ~ as.numeric(age_bracket)
  )) %>%
  group_by(iter, intervention, season_birth) %>%
  summarise(av_age = sum(age_bracket * n_icu) / sum(n_icu)) %>%
  ungroup() %>%
  group_by(intervention, season_birth) %>%
  summarise(
    av_age_low_95 = quantile(av_age, 0.025),
    av_age_median = quantile(av_age, 0.5),
    av_age_up_95 = quantile(av_age, 0.975)
  ) %>%
  ungroup() %>%
  rbind(
    total_outcome_intervention_df %>%
      ungroup() %>%
      filter(age_bracket != "all" & age_bracket != "12-14" &
        intervention == "MV") %>%
      mutate(age_bracket = case_when(
        age_bracket == "<1" ~ 0,
        T ~ as.numeric(age_bracket)
      )) %>%
      group_by(iter, age_bracket, intervention) %>%
      summarise(
        season_birth = "all",
        n_hospitalisations = sum(n_hospitalisations),
        prevented_hospitalisations = sum(prevented_hospitalisations),
        prevented_icu = sum(prevented_icu),
        n_icu = sum(n_icu)
      ) %>%
      ungroup() %>%
      group_by(iter, intervention, season_birth) %>%
      summarise(av_age = sum(age_bracket * n_icu) / sum(n_icu)) %>%
      ungroup() %>%
      group_by(intervention, season_birth) %>%
      summarise(
        av_age_low_95 = quantile(av_age, 0.025),
        av_age_median = quantile(av_age, 0.5),
        av_age_up_95 = quantile(av_age, 0.975)
      ) %>%
      ungroup()
  ) %>%
  arrange(season_birth)

av_age_icu %>% write.csv(paste0(path_output, "Avg age at ICU admission.csv"), row.names = F)

# ------------ NNVs --------------------
nnv_hosp <- total_outcome_intervention_df %>%
  filter(age_bracket == "all") %>%
  select(intervention, season_birth, prevented_hospitalisations, iter) %>%
  merge(births_de %>% select(season, total),
    by.x = "season_birth", by.y = "season"
  ) %>%
  mutate(
    NNV = case_when(
      intervention == "No immunisation" ~ NA,
      T ~ total / prevented_hospitalisations
    )
  ) %>%
  group_by(intervention, season_birth) %>%
  summarise(
    NNV_median = median(NNV, na.rm = T),
    NNV_low95 = quantile(NNV, 0.025, na.rm = T),
    NNV_up95 = quantile(NNV, 0.975, na.rm = T)
  ) %>%
  mutate(
    NNV_median = ifelse(is.infinite(NNV_median), NA, NNV_median),
    NNV_low95 = ifelse(is.infinite(NNV_low95), NA, NNV_low95),
    NNV_up95 = ifelse(is.infinite(NNV_up95), NA, NNV_up95)
  )

nnv_icu <- total_outcome_intervention_df %>%
  filter(age_bracket == "all") %>%
  select(intervention, season_birth, prevented_icu, iter) %>%
  merge(births_de %>% select(season, total),
    by.x = "season_birth", by.y = "season"
  ) %>%
  mutate(
    NNV = case_when(
      intervention == "No immunisation" ~ NA,
      T ~ total / prevented_icu
    )
  ) %>%
  group_by(intervention, season_birth) %>%
  summarise(
    NNV_median = median(NNV, na.rm = T),
    NNV_low95 = quantile(NNV, 0.025, na.rm = T),
    NNV_up95 = quantile(NNV, 0.975, na.rm = T)
  ) %>%
  mutate(
    NNV_median = ifelse(is.infinite(NNV_median), NA, NNV_median),
    NNV_low95 = ifelse(is.infinite(NNV_low95), NA, NNV_low95),
    NNV_up95 = ifelse(is.infinite(NNV_up95), NA, NNV_up95)
  )

# Save outputs
write.csv(nnv_hosp, file = paste0(path_output, "NNV_hosp.csv"), row.names = F)
write.csv(nnv_icu, file = paste0(path_output, "NNV_icu.csv"), row.names = F)

# ------------- Incidence by age group ----------------------------------------
# We need to divide births into 12 age groups
births_monthly <- births_de %>%
  select(season, total) %>%
  mutate(total = total / 12) %>%
  # repeat each season 12 times to get monthly numbers
  slice(rep(1:n(), each = 12)) %>%
  group_by(season) %>%
  mutate(age_months = row_number()) %>%
  ungroup()

incidence_outcome_all_inter <- total_outcome_intervention_df %>%
  filter(intervention == "No immunisation" | intervention == "MV") %>%
  rbind(total_outcome_intervention_df %>%
    filter(intervention == "+ mAB for autumn births" &
      season_birth == "autumn") %>%
    mutate(intervention = "MV + mAB")) %>%
  rbind(total_outcome_intervention_df %>%
    filter(intervention == "+ mAB for winter births" &
      season_birth == "winter") %>%
    mutate(intervention = "MV + mAB")) %>%
  rbind(total_outcome_intervention_df %>%
    filter(intervention == "+ mAB for summer births" &
      season_birth == "summer") %>%
    mutate(intervention = "MV + mAB")) %>%
  rbind(total_outcome_intervention_df %>%
    filter(intervention == "+ mAB for spring births" &
      season_birth == "spring") %>%
    mutate(intervention = "MV + mAB")) %>%
  filter(age_bracket != "all") %>%
  group_by(iter, intervention, season_birth, age_bracket) %>%
  summarise(
    n_hospitalisations = sum(n_hospitalisations),
    n_icu = sum(n_icu)
  )

incidence_hosp_monthly <- incidence_outcome_all_inter %>%
  merge(births_monthly,
    by.x = c("season_birth", "age_bracket"),
    by.y = c("season", "age_months")
  ) %>%
  mutate(
    incidence_per_100000 = n_hospitalisations / total * 100000,
    icu_per_100000 = n_icu / total * 100000,
    age_bracket = as.numeric(age_bracket)
  )

incidence_outcome_monthly_int <- incidence_hosp_monthly %>%
  group_by(intervention, season_birth, age_bracket) %>%
  summarise(
    incidence_up_95 = quantile(incidence_per_100000, 0.025, na.rm = T),
    incidence_median = quantile(incidence_per_100000, 0.5, na.rm = T),
    incidence_low_95 = quantile(incidence_per_100000, 0.975, na.rm = T),
    icu_inc_up_95 = quantile(icu_per_100000, 0.025, na.rm = T),
    icu_inc_median = quantile(icu_per_100000, 0.5, na.rm = T),
    icu_inc_low_95 = quantile(icu_per_100000, 0.975, na.rm = T)
  ) %>%
  ungroup() %>%
  mutate(season_birth = str_to_sentence(season_birth)) %>%
  mutate(season_birth = factor(season_birth,
    levels = c(
      "Winter", "Spring", "Summer",
      "Autumn", "All"
    )
  ))

# Cumulative incidence by age group
incidence_outcome_monthly_cum <- incidence_hosp_monthly %>%
  arrange(intervention, season_birth, age_bracket, iter) %>%
  group_by(intervention, season_birth, iter) %>%
  mutate(
    cum_incidence_per_100000 = cumsum(incidence_per_100000),
    cum_icu_per_100000 = cumsum(icu_per_100000)
  ) %>%
  ungroup()

incidence_outcome_monthly_cum_int <- incidence_outcome_monthly_cum %>%
  group_by(intervention, season_birth, age_bracket) %>%
  summarise(
    cum_incidence_up_95 = quantile(cum_incidence_per_100000, 0.025, na.rm = T),
    cum_incidence_median = quantile(cum_incidence_per_100000, 0.5, na.rm = T),
    cum_incidence_low_95 = quantile(cum_incidence_per_100000, 0.975, na.rm = T),
    cum_icu_inc_up_95 = quantile(cum_icu_per_100000, 0.025, na.rm = T),
    cum_icu_inc_median = quantile(cum_icu_per_100000, 0.5, na.rm = T),
    cum_icu_inc_low_95 = quantile(cum_icu_per_100000, 0.975, na.rm = T)
  ) %>%
  ungroup() %>%
  mutate(season_birth = str_to_sentence(season_birth)) %>%
  mutate(season_birth = factor(season_birth, levels = c("Winter", "Spring", "Summer", "Autumn")))

# ---------- Proportion of disease progression ---------------------------------
# Bind the two DFs
infection_to_disease <- total_outcome_intervention_df %>%
  group_by(season_birth, age_bracket, intervention, iter) %>%
  summarise(
    n_hospitalisations = sum(n_hospitalisations),
    n_icu = sum(n_icu)
  ) %>%
  ungroup() %>%
  merge(
    infection_df %>%
      mutate(age_bracket = trunc(age_midpoint / 30)) %>%
      group_by(season_birth, age_bracket, iter) %>%
      summarise(n_infections = sum(n_infections)) %>%
      ungroup(),
    by.x = c("age_bracket", "season_birth", "iter"),
    by.y = c("age_bracket", "season_birth", "iter")
  )

# Add a row for the first 2 months of life
infection_to_disease <- infection_to_disease %>%
  rbind(
    infection_to_disease %>%
      filter(age_bracket == 0 | age_bracket == 1) %>%
      group_by(season_birth, intervention, iter) %>%
      summarise(
        age_bracket = "0-1",
        n_hospitalisations = sum(n_hospitalisations),
        n_icu = sum(n_icu),
        n_infections = sum(n_infections)
      ) %>%
      ungroup()
  )

# Add totals
infection_to_disease <- infection_to_disease %>%
  rbind(
    infection_to_disease %>%
      filter(age_bracket != "0-1") %>%
      group_by(season_birth, intervention, iter) %>%
      summarise(
        n_hospitalisations = sum(n_hospitalisations),
        n_icu = sum(n_icu),
        n_infections = sum(n_infections)
      ) %>%
      ungroup() %>%
      mutate(age_bracket = "all")
  )

infection_to_disease_prop <- infection_to_disease %>%
  group_by(season_birth, age_bracket, intervention, iter) %>%
  summarise(
    prop_hosp = n_hospitalisations / n_infections,
    prop_icu = n_icu / n_hospitalisations
  ) %>%
  ungroup() %>%
  group_by(season_birth, age_bracket, intervention) %>%
  summarise(
    prop_hosp_up_95 = quantile(prop_hosp, 0.025, na.rm = T),
    prop_hosp_median = quantile(prop_hosp, 0.5, na.rm = T),
    prop_hosp_low_95 = quantile(prop_hosp, 0.975, na.rm = T),
    prop_icu_up_95 = quantile(prop_icu, 0.025, na.rm = T),
    prop_icu_median = quantile(prop_icu, 0.5, na.rm = T),
    prop_icu_low_95 = quantile(prop_icu, 0.975, na.rm = T)
  ) %>%
  ungroup() %>%
  mutate(season_birth = str_to_sentence(season_birth)) %>%
  mutate(season_birth = factor(season_birth, levels = c("Winter", "Spring", "Summer", "Autumn", "All")))

infection_to_disease_prop %>% write.csv(paste0(path_output, "Infection to disease progression.csv"), row.names = F)

# Pivot wider for combined plotting
infection_to_disease_prop_wide <- infection_to_disease_prop %>%
  filter(age_bracket != "0-1") %>%
select(!c(prop_icu_up_95, prop_icu_median, prop_icu_low_95)) %>%
  rename(
    prop_up_95 = prop_hosp_up_95,
    prop_median = prop_hosp_median,
    prop_low_95 = prop_hosp_low_95
  ) %>%
  mutate(progression = "Infections to hospitalisations") %>%
  rbind(
    infection_to_disease_prop %>%
      select(!c(prop_hosp_up_95, prop_hosp_median, prop_hosp_low_95)) %>%
      rename(
        prop_up_95 = prop_icu_up_95,
        prop_median = prop_icu_median,
        prop_low_95 = prop_icu_low_95
      ) %>%
      mutate(progression = "Hospitalisations to ICU admissions")
  ) %>%
  mutate(progression = factor(progression, levels = c(
    "Infections to hospitalisations",
    "Hospitalisations to ICU admissions"
  )))

# ------------ Reduction in hosp and ICU visits with mAB compared to MV -------
reduction_mAB_MV <- total_outcome_intervention_df %>%
  filter(intervention == "MV" |
    (season_birth == "all" & intervention == "+ mAB for all") |
    (season_birth == "spring" & intervention == "+ mAB for spring births") |
    (season_birth == "summer" & intervention == "+ mAB for summer births") |
    (season_birth == "autumn" & intervention == "+ mAB for autumn births") |
    (season_birth == "winter" & intervention == "+ mAB for winter births")) %>%
  group_by(iter, season_birth, age_bracket, intervention) %>%
  summarise(
    n_hospitalisations = (n_hospitalisations),
    n_icu = (n_icu)
  ) %>%
  mutate(intervention = case_when(
    intervention != "MV" ~ "MV + mAB",
    T ~ intervention
  )) %>%
  ungroup() %>%
  pivot_wider(
    names_from = intervention,
    values_from = c(n_hospitalisations, n_icu)
  ) %>%
  group_by(iter, season_birth, age_bracket) %>%
  summarise(
    prevented_hospitalisations_vs_MV = (n_hospitalisations_MV - `n_hospitalisations_MV + mAB`),
    prevented_icu_vs_MV = (n_icu_MV - `n_icu_MV + mAB`)
  ) %>%
  ungroup()

# Reduction as incidence
reduction_mAB_MV_inc <- total_outcome_intervention_df %>%
  filter(intervention == "MV" |
    (season_birth == "all" & intervention == "+ mAB for all") |
    (season_birth == "spring" & intervention == "+ mAB for spring births") |
    (season_birth == "summer" & intervention == "+ mAB for summer births") |
    (season_birth == "autumn" & intervention == "+ mAB for autumn births") |
    (season_birth == "winter" & intervention == "+ mAB for winter births")) %>%
  group_by(iter, season_birth, age_bracket, intervention) %>%
  summarise(
    n_hospitalisations = (n_hospitalisations),
    n_icu = (n_icu)
  ) %>%
  mutate(intervention = case_when(
    intervention != "MV" ~ "MV + mAB",
    T ~ intervention
  )) %>%
  filter(age_bracket == "all") %>%
  ungroup() %>%
  merge(births_de,
    by.x = c("season_birth"),
    by.y = c("season")
  ) %>%
  mutate(
    incidence_per_100000 = n_hospitalisations / total * 100000,
    icu_per_100000 = n_icu / total * 100000
  ) %>%
  select(!c(N, total, prop)) %>%
  pivot_wider(
    names_from = intervention,
    values_from = c(
      incidence_per_100000, icu_per_100000,
      n_hospitalisations, n_icu
    )
  ) %>%
  group_by(iter, season_birth) %>%
  summarise(
    prevented_hospitalisations_vs_MV = (incidence_per_100000_MV - `incidence_per_100000_MV + mAB`),
    prevented_icu_vs_MV = (icu_per_100000_MV - `icu_per_100000_MV + mAB`),
    incidence_per_100000_MV = incidence_per_100000_MV,
    icu_per_100000_MV = icu_per_100000_MV,
    incidence_per_100000_MV_mAB = `incidence_per_100000_MV + mAB`,
    icu_per_100000_MV_mAB = `icu_per_100000_MV + mAB`
  ) %>%
  ungroup()

reduction_mAB_MV_inc_CI <- reduction_mAB_MV_inc %>%
  group_by(season_birth) %>%
  summarise(
    hosp_reduction_hosp_lower = quantile(prevented_hospitalisations_vs_MV, 0.025, na.rm = TRUE),
    hosp_reduction_hosp_median = quantile(prevented_hospitalisations_vs_MV, 0.5, na.rm = TRUE),
    hosp_reduction_hosp_upper = quantile(prevented_hospitalisations_vs_MV, 0.975, na.rm = TRUE),
    hosp_reduction_icu_lower = quantile(prevented_icu_vs_MV, 0.025, na.rm = TRUE),
    hosp_reduction_icu_median = quantile(prevented_icu_vs_MV, 0.5, na.rm = TRUE),
    hosp_reduction_icu_upper = quantile(prevented_icu_vs_MV, 0.975, na.rm = TRUE),
    incidence_per_100000_MV_mAB_lower = quantile(incidence_per_100000_MV_mAB, 0.025, na.rm = TRUE),
    incidence_per_100000_MV_mAB_median = quantile(incidence_per_100000_MV_mAB, 0.5, na.rm = TRUE),
    incidence_per_100000_MV_mAB_upper = quantile(incidence_per_100000_MV_mAB, 0.975, na.rm = TRUE),
    icu_per_100000_MV_mAB_lower = quantile(icu_per_100000_MV_mAB, 0.025, na.rm = TRUE),
    icu_per_100000_MV_mAB_median = quantile(icu_per_100000_MV_mAB, 0.5, na.rm = TRUE),
    icu_per_100000_MV_mAB_upper = quantile(icu_per_100000_MV_mAB, 0.975, na.rm = TRUE)
  )

# Percentage prevented compared to MV
prev_mAB_MV_perc_spring <- reduction_mAB_MV %>%
  filter(season_birth == "spring") %>%
  merge(
    total_outcome_intervention_df %>%
      filter(intervention == "MV" & season_birth == "spring") %>%
      mutate(
        prev_MV = prevented_hospitalisations,
        prev_icu_MV = prevented_icu
      ) %>%
      select(age_bracket, iter, season_birth, prev_MV, prev_icu_MV),
    by = c("age_bracket", "iter", "season_birth")
  )


prev_mAB_MV_perc_summer <- reduction_mAB_MV %>%
  filter(season_birth == "summer") %>%
  merge(
    total_outcome_intervention_df %>%
      filter(intervention == "MV" & season_birth == "summer") %>%
      mutate(
        prev_MV = prevented_hospitalisations,
        prev_icu_MV = prevented_icu
      ) %>%
      select(age_bracket, iter, season_birth, prev_MV, prev_icu_MV),
    by = c("age_bracket", "iter", "season_birth")
  )


prev_mAB_MV_perc_autumn <- reduction_mAB_MV %>%
  filter(season_birth == "autumn") %>%
  merge(
    total_outcome_intervention_df %>%
      filter(intervention == "MV" & season_birth == "autumn") %>%
      mutate(
        prev_MV = prevented_hospitalisations,
        prev_icu_MV = prevented_icu
      ) %>%
      select(age_bracket, iter, season_birth, prev_MV, prev_icu_MV),
    by = c("age_bracket", "iter", "season_birth")
  )

prev_mAB_MV_perc_winter <- reduction_mAB_MV %>%
  filter(season_birth == "winter") %>%
  merge(
    total_outcome_intervention_df %>%
      filter(intervention == "MV" & season_birth == "winter") %>%
      mutate(
        prev_MV = prevented_hospitalisations,
        prev_icu_MV = prevented_icu
      ) %>%
      select(age_bracket, iter, season_birth, prev_MV, prev_icu_MV),
    by = c("age_bracket", "iter", "season_birth")
  )

prev_mAB_MV_perc_all <- reduction_mAB_MV %>%
  filter(season_birth == "all") %>%
  merge(
    total_outcome_intervention_df %>%
      filter(intervention == "MV" & season_birth == "all") %>%
      mutate(
        prev_MV = prevented_hospitalisations,
        prev_icu_MV = prevented_icu
      ) %>%
      select(age_bracket, iter, season_birth, prev_MV, prev_icu_MV),
    by = c("age_bracket", "iter", "season_birth")
  )

prev_mAB_MV_perc <- prev_mAB_MV_perc_spring %>%
  rbind(prev_mAB_MV_perc_summer) %>%
  rbind(prev_mAB_MV_perc_autumn) %>%
  rbind(prev_mAB_MV_perc_winter) %>%
  rbind(prev_mAB_MV_perc_all)

prev_mAB_MV_perc <- prev_mAB_MV_perc %>%
  mutate(
    perc_reduction_hosp = prevented_hospitalisations_vs_MV / prev_MV * 100,
    perc_reduction_icu = prevented_icu_vs_MV / prev_icu_MV * 100
  )

prev_mAB_MV_perc_CI <- prev_mAB_MV_perc %>%
  group_by(age_bracket, season_birth) %>%
  summarise(
    hosp_reduction_hosp_lower = quantile(perc_reduction_hosp, 0.025, na.rm = TRUE),
    hosp_reduction_hosp_median = quantile(perc_reduction_hosp, 0.5, na.rm = TRUE),
    hosp_reduction_hosp_upper = quantile(perc_reduction_hosp, 0.975, na.rm = TRUE),
    hosp_reduction_icu_lower = quantile(perc_reduction_icu, 0.025, na.rm = TRUE),
    hosp_reduction_icu_median = quantile(perc_reduction_icu, 0.5, na.rm = TRUE),
    hosp_reduction_icu_upper = quantile(perc_reduction_icu, 0.975, na.rm = TRUE)
  )

prev_mAB_MV_perc_CI %>% write.csv(paste0(path_output, "Prevented mAB vs MV percentage.csv"), row.names = F)
reduction_mAB_MV_inc_CI %>% write.csv(paste0(path_output, "Prevented mAB vs MV incidence.csv"), row.names = F)

# Percentage reduction compared to MV
reduction_mAB_MV_perc_spring <- reduction_mAB_MV %>%
  filter(season_birth == "spring") %>%
  merge(
    total_outcome_intervention_df %>%
      filter(intervention == "MV" & season_birth == "spring") %>%
      mutate(
        hosp_MV = n_hospitalisations,
        n_icu_MV = n_icu
      ) %>%
      select(age_bracket, iter, season_birth, hosp_MV, n_icu_MV),
    by = c("age_bracket", "iter", "season_birth")
  )

reduction_mAB_MV_perc_summer <- reduction_mAB_MV %>%
  filter(season_birth == "summer") %>%
  merge(
    total_outcome_intervention_df %>%
      filter(intervention == "MV" & season_birth == "summer") %>%
      mutate(
        hosp_MV = n_hospitalisations,
        n_icu_MV = n_icu
      ) %>%
      select(age_bracket, iter, season_birth, hosp_MV, n_icu_MV),
    by = c("age_bracket", "iter", "season_birth")
  )

reduction_mAB_MV_perc_autumn <- reduction_mAB_MV %>%
  filter(season_birth == "autumn") %>%
  merge(
    total_outcome_intervention_df %>%
      filter(intervention == "MV" & season_birth == "autumn") %>%
      mutate(
        hosp_MV = n_hospitalisations,
        n_icu_MV = n_icu
      ) %>%
      select(age_bracket, iter, season_birth, hosp_MV, n_icu_MV),
    by = c("age_bracket", "iter", "season_birth")
  )

reduction_mAB_MV_perc_winter <- reduction_mAB_MV %>%
  filter(season_birth == "winter") %>%
  merge(
    total_outcome_intervention_df %>%
      filter(intervention == "MV" & season_birth == "winter") %>%
      mutate(
        hosp_MV = n_hospitalisations,
        n_icu_MV = n_icu
      ) %>%
      select(age_bracket, iter, season_birth, hosp_MV, n_icu_MV),
    by = c("age_bracket", "iter", "season_birth")
  )

reduction_mAB_MV_perc_all <- reduction_mAB_MV %>%
  filter(season_birth == "all") %>%
  merge(
    total_outcome_intervention_df %>%
      filter(intervention == "MV" & season_birth == "all") %>%
      mutate(
        hosp_MV = n_hospitalisations,
        n_icu_MV = n_icu
      ) %>%
      select(age_bracket, iter, season_birth, hosp_MV, n_icu_MV),
    by = c("age_bracket", "iter", "season_birth")
  )

reduction_mAB_MV_perc_2 <- reduction_mAB_MV_perc_spring %>%
  rbind(reduction_mAB_MV_perc_summer) %>%
  rbind(reduction_mAB_MV_perc_autumn) %>%
  rbind(reduction_mAB_MV_perc_winter) %>%
  rbind(reduction_mAB_MV_perc_all)

reduction_mAB_MV_perc_2 <- reduction_mAB_MV_perc_2 %>%
  mutate(
    perc_reduction_hosp = prevented_hospitalisations_vs_MV / hosp_MV * 100,
    perc_reduction_icu = prevented_icu_vs_MV / n_icu_MV * 100
  )

reduction_mAB_MV_perc_2_CI <- reduction_mAB_MV_perc_2 %>%
  group_by(age_bracket, season_birth) %>%
  summarise(
    hosp_reduction_hosp_lower = quantile(perc_reduction_hosp, 0.025, na.rm = TRUE),
    hosp_reduction_hosp_median = quantile(perc_reduction_hosp, 0.5, na.rm = TRUE),
    hosp_reduction_hosp_upper = quantile(perc_reduction_hosp, 0.975, na.rm = TRUE),
    hosp_reduction_icu_lower = quantile(perc_reduction_icu, 0.025, na.rm = TRUE),
    hosp_reduction_icu_median = quantile(perc_reduction_icu, 0.5, na.rm = TRUE),
    hosp_reduction_icu_upper = quantile(perc_reduction_icu, 0.975, na.rm = TRUE)
  )

reduction_mAB_MV_perc_2_CI %>% write.csv(paste0(path_output, "Prevented mAB vs MV percentage by age.csv"), row.names = F)

# NNV for mAB compared to MV
nnv_mAB_MV <- reduction_mAB_MV %>%
  filter(age_bracket == "all") %>%
  merge(births_de, by.x = "season_birth", by.y = "season") %>%
  mutate(
    NNV_hosp = case_when(
      prevented_hospitalisations_vs_MV == 0 ~ NA,
      T ~ total / prevented_hospitalisations_vs_MV
    ),
    NNV_icu = case_when(
      prevented_icu_vs_MV == 0 ~ NA,
      T ~ total / prevented_icu_vs_MV
    )
  ) %>%
  group_by(season_birth) %>%
  summarise(
    NNV_hosp_median = median(NNV_hosp, na.rm = T),
    NNV_hosp_low95 = quantile(NNV_hosp, 0.025, na.rm = T),
    NNV_hosp_up95 = quantile(NNV_hosp, 0.975, na.rm = T),
    NNV_icu_median = median(NNV_icu, na.rm = T),
    NNV_icu_low95 = quantile(NNV_icu, 0.025, na.rm = T),
    NNV_icu_up95 = quantile(NNV_icu, 0.975, na.rm = T)
  ) %>%
  mutate(
    NNV_hosp_median = ifelse(is.infinite(NNV_hosp_median), NA, NNV_hosp_median),
    NNV_hosp_low95 = ifelse(is.infinite(NNV_hosp_low95), NA, NNV_hosp_low95),
    NNV_hosp_up95 = ifelse(is.infinite(NNV_hosp_up95), NA, NNV_hosp_up95),
    NNV_icu_median = ifelse(is.infinite(NNV_icu_median), NA, NNV_icu_median),
    NNV_icu_low95 = ifelse(is.infinite(NNV_icu_low95), NA, NNV_icu_low95),
    NNV_icu_up95 = ifelse(is.infinite(NNV_icu_up95), NA, NNV_icu_up95)
  )
nnv_mAB_MV %>% write.csv(paste0(path_output, " NNV mAB vs MV.csv"), row.names = F)

# Get FOI by season compared to summer
samples_tuned <- readRDS("/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/PhD project/Project 2/Code base/RSV_MV_mAB_catchup/Fitting output/samples season.rds")
samples_thinned <- monty_samples_thin(samples_tuned, burnin = 3000)
draws_thinned <- as_draws_df(samples_thinned) %>%
  rename_variables(
    "Spring component of the FOI" = "spring_comp",
    "Summer component of FOI" = "summer_comp",
    "Autumn component of the FOI" = "autumn_comp",
    "Winter component of the FOI" = "winter_comp",
    "ω" = "mu",
    "π" = "prop"
  )

params_relative <- draws_thinned %>%
  select(
    `Spring component of the FOI`,
    `Summer component of FOI`,
    `Autumn component of the FOI`,
    `Winter component of the FOI`,
    `.chain`, `.iteration`, `.draw`
  ) %>%
  group_by(`.chain`, `.iteration`, `.draw`) %>%
  mutate(
    foi_spring = `Spring component of the FOI` + `Summer component of FOI`,
    foi_summer = `Summer component of FOI`,
    foi_autumn = `Autumn component of the FOI` + `Summer component of FOI`,
    foi_winter = `Winter component of the FOI` + `Summer component of FOI`
  ) %>%
  ungroup() %>%
  group_by(`.chain`, `.iteration`, `.draw`) %>%
  mutate(
    rel_foi_spring = (foi_spring - foi_summer) / foi_summer,
    rel_foi_summer = (foi_summer - foi_summer) / foi_summer,
    rel_foi_autumn = (foi_autumn - foi_summer) / foi_summer,
    rel_foi_winter = (foi_winter - foi_summer) / foi_summer
  ) %>%
  ungroup()

params_relative_summary <- params_relative %>%
  summarise(
    foi_spring_low95 = quantile(foi_spring, 0.025, na.rm = T),
    foi_spring_median = quantile(foi_spring, 0.5, na.rm = T),
    foi_spring_up95 = quantile(foi_spring, 0.975, na.rm = T),
    foi_summer_low95 = quantile(foi_summer, 0.025, na.rm = T),
    foi_summer_median = quantile(foi_summer, 0.5, na.rm = T),
    foi_summer_up95 = quantile(foi_summer, 0.975, na.rm = T),
    foi_autumn_low95 = quantile(foi_autumn, 0.025, na.rm = T),
    foi_autumn_median = quantile(foi_autumn, 0.5, na.rm = T),
    foi_autumn_up95 = quantile(foi_autumn, 0.975, na.rm = T),
    foi_winter_low95 = quantile(foi_winter, 0.025, na.rm = T),
    foi_winter_median = quantile(foi_winter, 0.5, na.rm = T),
    foi_winter_up95 = quantile(foi_winter, 0.975, na.rm = T),
    rel_foi_spring_low95 = quantile(rel_foi_spring, 0.025, na.rm = T),
    rel_foi_spring_median = quantile(rel_foi_spring, 0.5, na.rm = T),
    rel_foi_spring_up95 = quantile(rel_foi_spring, 0.975, na.rm = T),
    rel_foi_summer_low95 = quantile(rel_foi_summer, 0.025, na.rm = T),
    rel_foi_summer_median = quantile(rel_foi_summer, 0.5, na.rm = T),
    rel_foi_summer_up95 = quantile(rel_foi_summer, 0.975, na.rm = T),
    rel_foi_autumn_low95 = quantile(rel_foi_autumn, 0.025, na.rm = T),
    rel_foi_autumn_median = quantile(rel_foi_autumn, 0.5, na.rm = T),
    rel_foi_autumn_up95 = quantile(rel_foi_autumn, 0.975, na.rm = T),
    rel_foi_winter_low95 = quantile(rel_foi_winter, 0.025, na.rm = T),
    rel_foi_winter_median = quantile(rel_foi_winter, 0.5, na.rm = T),
    rel_foi_winter_up95 = quantile(rel_foi_winter, 0.975, na.rm = T)
  )

# Pivot longer
params_relative_summary_1 <- params_relative_summary %>%
  select(
    foi_spring_low95, foi_spring_median, foi_spring_up95,
    foi_summer_low95, foi_summer_median, foi_summer_up95,
    foi_autumn_low95, foi_autumn_median, foi_autumn_up95,
    foi_winter_low95, foi_winter_median, foi_winter_up95
  ) %>%
  pivot_longer(
    cols = everything(),
    names_to = c("season", ".value"),
    names_pattern = "foi_(.*)_(.*)"
  )


params_relative_summary_2 <- params_relative_summary %>%
  select(
    rel_foi_spring_low95, rel_foi_spring_median, rel_foi_spring_up95,
    rel_foi_summer_low95, rel_foi_summer_median, rel_foi_summer_up95,
    rel_foi_autumn_low95, rel_foi_autumn_median, rel_foi_autumn_up95,
    rel_foi_winter_low95, rel_foi_winter_median, rel_foi_winter_up95
  ) %>%
  pivot_longer(
    cols = everything(),
    names_to = c("season", ".value"),
    names_pattern = "rel_foi_(.*)_(.*)"
  ) %>%
  rename(
    rel_low95 = low95,
    rel_median = median,
    rel_up95 = up95
  )

# Join both
params_relative_summary_wide <- left_join(params_relative_summary_1,
  params_relative_summary_2,
  by = "season"
)

# Save file
params_relative_summary_wide %>% write.csv(paste0(path_output, "Relative FOI compared to summer.csv"))

# ----------- Proportion of total outcomes by season of birth ----------
proportion_outcomes_season_birth <- total_outcome_intervention_df %>%
  filter(age_bracket == "all" & season_birth != "all") %>%
  group_by(iter, intervention, season_birth) %>%
  summarise(
    n_hospitalisations = sum(n_hospitalisations),
    n_icu = sum(n_icu)
  ) %>%
  ungroup() %>%
  group_by(iter, intervention) %>%
  mutate(
    total_hosp = sum(n_hospitalisations),
    total_icu = sum(n_icu)
  ) %>%
  ungroup() %>%
  mutate(
    proportion_hosp = n_hospitalisations / total_hosp,
    proportion_icu = n_icu / total_icu
  ) %>%
  group_by(intervention, season_birth) %>%
  summarise(
    prop_hosp_low_95 = quantile(proportion_hosp, 0.025, na.rm = T),
    prop_hosp_median = quantile(proportion_hosp, 0.5, na.rm = T),
    prop_hosp_up_95 = quantile(proportion_hosp, 0.975, na.rm = T),
    prop_icu_low_95 = quantile(proportion_icu, 0.025, na.rm = T),
    prop_icu_median = quantile(proportion_icu, 0.5, na.rm = T),
    prop_icu_up_95 = quantile(proportion_icu, 0.975, na.rm = T)
  ) %>%
  ungroup() %>%
  mutate(season_birth = str_to_sentence(season_birth)) %>%
  mutate(season_birth = factor(season_birth,
    levels = c(
      "Winter", "Spring", "Summer",
      "Autumn", "Winter + autumn"
    )
  ))

proportion_outcomes_season_birth %>% write.csv(paste0(path_output, "Proportion outcomes by season of birth.csv"), row.names = F)

# -------------- In and out-of-season cohorts for no intervention and MV -------
# Infection incidence 
infection_age <- infection_df %>%
  mutate(age_bracket = trunc(age_midpoint / 30)) %>%
  group_by(season_birth, age_bracket, iter) %>%
  summarise(n_infections = sum(n_infections)) %>%
  select(season_birth, age_bracket, iter, n_infections)

infection_age_tot <- infection_age %>%
  mutate(age_bracket = as.character(age_bracket)) %>%
  rbind(
    infection_age %>%
      group_by(season_birth, iter) %>%
      summarise(
        age_bracket = "all",
        n_infections = sum(n_infections)
      )
  )

incidence_infection_in_out <- infection_age_tot %>%
  filter(age_bracket == "all") %>%
  merge(
    births_de %>% select(season, total),
    by.x = "season_birth", by.y = "season"
  ) %>%
  mutate(RSV_season = case_when (season_birth == "winter" | season_birth == "autumn" ~ "In season",
                                 season_birth == "summer" | season_birth == "spring" ~ "Out of season",
                                 season_birth == "all" ~ "All")) %>%
  group_by(RSV_season, age_bracket, iter) %>%
  summarise(n_infections = sum(n_infections),
            total = sum(total)) %>%
  ungroup() %>%
  mutate(incidence_per_100000 = n_infections / total * 100000)

# Get 95% CIs
incidence_infection_in_out_CI <- incidence_infection_in_out %>%
  group_by(RSV_season, age_bracket) %>%
  summarise(
    incidence_lower = quantile(incidence_per_100000, 0.025),
    incidence_median = quantile(incidence_per_100000, 0.5),
    incidence_upper = quantile(incidence_per_100000, 0.975)
  ) %>%
  ungroup()

incidence_infection_in_out_CI %>% write.csv(paste0(path_output, "Incidence infection in and out of season.csv"), row.names = F)

# Get outcome incidence per 100,000
incidence_hosp_in_out <- total_outcome_intervention_df %>%
  filter((intervention == "No immunisation" | intervention == "MV") &
           age_bracket == "all") %>%
  merge(births_de %>% select(season, total), by.x = "season_birth", by.y = "season") %>%
  mutate(RSV_season = case_when (season_birth == "winter" | season_birth == "autumn" ~ "In season",
                                 season_birth == "spring" | season_birth == "summer" ~ "Out of season",
                                 T ~ NA)) %>%
  group_by(RSV_season, age_bracket, intervention, iter) %>%
  summarise(
    n_hospitalisations = sum(n_hospitalisations),
    prevented_hospitalisations = sum(prevented_hospitalisations),
    prevented_icu = sum(prevented_icu),
    n_icu = sum(n_icu),
    total = sum(total)
  ) %>%
  ungroup() %>%
  mutate(
    incidence_per_100000 = n_hospitalisations / total * 100000,
    icu_per_100000 = n_icu / total * 100000,
    prevented_hospitalisations_per_100000 = prevented_hospitalisations / total * 100000,
    prevented_icu_per_100000 = prevented_icu / total * 100000
  )

# Get reduction in incidence compared to no immunisation
incidence_reduction_no_imm_in_out <- incidence_hosp_in_out %>%
  filter(intervention != "No immunisation" & age_bracket == "all") %>%
  merge(
    incidence_hosp_in_out %>%
      filter(intervention == "No immunisation" & age_bracket == "all") %>%
      select(RSV_season, age_bracket, n_hospitalisations, n_icu, iter) %>%
      rename(
        n_hospitalisations_no_int = n_hospitalisations,
        n_icu_no_int = n_icu
      ),
    by = c("RSV_season", "age_bracket", "iter")
  ) %>%
  mutate(
    incidence_reduction_perc = (n_hospitalisations_no_int - n_hospitalisations) / n_hospitalisations_no_int * 100,
    icu_reduction_perc = (n_icu_no_int - n_icu) / n_icu_no_int * 100
  )


# Compute 95% CI
incidence_outcome_in_out_int <- incidence_hosp_in_out %>%
  group_by(intervention, RSV_season, age_bracket) %>%
  summarise(
    incidence_up_95 = quantile(incidence_per_100000, 0.025, na.rm = T),
    incidence_median = quantile(incidence_per_100000, 0.5, na.rm = T),
    incidence_low_95 = quantile(incidence_per_100000, 0.975, na.rm = T),
    icu_inc_up_95 = quantile(icu_per_100000, 0.025, na.rm = T),
    icu_inc_median = quantile(icu_per_100000, 0.5, na.rm = T),
    icu_inc_low_95 = quantile(icu_per_100000, 0.975, na.rm = T),
    prev_inc_up_95 = quantile(prevented_hospitalisations_per_100000, 0.025, na.rm = T),
    prev_inc_median = quantile(prevented_hospitalisations_per_100000, 0.5, na.rm = T),
    prev_inc_low_95 = quantile(prevented_hospitalisations_per_100000, 0.975, na.rm = T),
    icu_prev_inc_up_95 = quantile(prevented_icu_per_100000, 0.025, na.rm = T),
    icu_prev_inc_median = quantile(prevented_icu_per_100000, 0.5, na.rm = T),
    icu_prev_inc_low_95 = quantile(prevented_icu_per_100000, 0.975, na.rm = T)
  ) %>%
  ungroup()

incidence_outcome_in_out_int %>% write.csv(paste0(path_output, "Incidence hosp icu in and out.csv"), row.names = F)

incidence_reduction_no_imm_in_out_int <- incidence_reduction_no_imm_in_out %>%
  group_by(intervention, RSV_season, age_bracket) %>%
  summarise(
    incidence_reduction_low_95 = quantile(incidence_reduction_perc, 0.025, na.rm = T),
    incidence_reduction_median = quantile(incidence_reduction_perc, 0.5, na.rm = T),
    incidence_reduction_up_95 = quantile(incidence_reduction_perc, 0.975, na.rm = T),
    icu_inc_reduction_low_95 = quantile(icu_reduction_perc, 0.025, na.rm = T),
    icu_inc_reduction_median = quantile(icu_reduction_perc, 0.5, na.rm = T),
    icu_inc_reduction_up_95 = quantile(icu_reduction_perc, 0.975, na.rm = T)
  ) %>%
  ungroup()
incidence_outcome_in_out_int %>% write.csv(paste0(path_output, "Incidence reduction no in and out.csv"), row.names = F)