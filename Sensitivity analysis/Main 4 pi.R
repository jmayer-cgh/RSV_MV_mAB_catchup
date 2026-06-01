# Fitting the odin model to the data

# Housekeeping
rm(list=ls())

# Load libraries
library(dust2)
library(odin2)
library(monty)
library(dplyr)
library(tidyr)
library(ggplot2)
library(posterior)
library(bayesplot)
library(magrittr)

# Path for outputs
path <- "/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/PhD project/Project 2/Code base/RSV_MV_mAB_catchup/Fitting outputs/Sensitivity analysis/4 pi/"

# ----------- Model ------------------------------------------------------------
# Load model
msr <- odin("/Users/juliamayer/Library/CloudStorage/OneDrive-Charité-UniversitätsmedizinBerlin/PhD project/Project 2/Code base/RSV_MV_mAB_catchup/Sensitivity analysis/Odin2 model 4 pi.R")

# ----------- Data -------------------------------------------------------------
# Read in and the data and put it in the right format
data <- read.csv("https://raw.githubusercontent.com/Stijn-A/RSV_serology/refs/heads/master/data/infection_status.csv")

# Group age into intervals 
# bi-monthly for 0-2 years and 6-monthly for 2-5 years
data$age_grp <- cut(data$age_days,
                    breaks = c(seq(0,730, by = 30.25*2),
                               seq(909,2000, by = 30.25*6)), 
                    include.lowest = T, right = F)

# Divide by season of birth
spring <- c(3, 4, 5)
summer <- c(6, 7, 8)
autumn <- c (9, 10, 11)
winter <- c(1, 2, 12)

data <- data %>%
  mutate(
    Birth_mo = format(as.Date(Birth_doy, origin = "2020-12-31"), "%m") %>% as.numeric(),
    season_birth = case_when (Birth_mo %in% spring ~ "spring",
                              Birth_mo %in% summer ~ "summer",
                              Birth_mo %in% autumn ~ "autumn",
                              Birth_mo %in% winter ~ "winter"))

get_midpoint <- function(cut_label) {
  round(mean(as.numeric(unlist(strsplit(gsub("\\(|\\)|\\[|\\]", "", as.character(cut_label)), ",")))))
}

data$xMidpoint <- sapply(data$age_grp, get_midpoint)

# Different groupings
# Get number of cases by age
incidence_data_season[,c("seroprev_mean","seroprev_low95","seroprev_up95")] <- binom::binom.confint(incidence_data_season$n_infection, 
                                                                                                    incidence_data_season$N, 
                                                                                                    method="exact")[,c("mean","lower","upper")]

incidence_data_season_wide <- incidence_data_season %>% 
  select (!c(age_mid, age_grp, seroprev_mean)) %>%
  pivot_wider(
    names_from = season_birth,
    values_from = c(N , n_infection, prop_seroconv, seroprev_low95, seroprev_up95, cum_infection),
    values_fill = 0
  ) %>%
  rename(time = "xMidpoint")

# ----------- Fitting ----------------------------------------------------------
# Build a filter and test it
filter <- dust_filter_create(msr, time_start = 0, data = incidence_data_season_wide, n_particles = 1000)

# Now running an MCMC
# List of inputs
packer <- monty_packer(c("spring_comp", "summer_comp", "autumn_comp", "winter_comp", 
                         "mu", 
                         "prop_sp", "prop_sm", "prop_au", "prop_wt"))

# Likelihood
likelihood <- dust_likelihood_monty(filter, packer, save_trajectories = T)
dust_likelihood_run(filter, list(spring_comp = 0.003478,
                                 summer_comp = 0.00060,
                                 autumn_comp = 0.00264,
                                 winter_comp = 0.00692, 
                                 mu = 0.09,
                                 prop_sp = 0.3,
                                 prop_sm = 0.01,
                                 prop_au = 0.3,
                                 prop_wt = 0.3),
                    save_trajectories = T)

# Priors
prior <- monty_dsl({
  spring_comp ~ Uniform(0, 0.1)
  summer_comp ~ Uniform(0, 0.1)
  autumn_comp ~ Uniform(0, 0.1)
  winter_comp ~ Uniform(0, 0.1)
  mu ~ Normal(1/181, 0.0011546) # centered around 200 days with 95% CI [137 - 365]
  prop_sp ~ Uniform(0, 1)
  prop_sm ~ Uniform(0, 1)
  prop_au ~ Uniform(0, 1)
  prop_wt ~ Uniform(0, 1)
})

# Posterior
posterior <- likelihood + prior

# Define a sampler (adaptive MCMC)
vcv <- diag(9) * 0.0004
sampler <- monty_sampler_adaptive(vcv)

# Parallelise the process
runner <- monty_runner_callr(2, progress = T)

# Sample
properties <- monty_model_properties(is_stochastic = F)
model <- monty_model(posterior, properties = properties)
samples <- monty_sample(model, sampler, 80000, 
                        initial = c(1e-05, 0.02002, 3e-05, 4e-05, 0.004, 
                                    0.5, 0.5, 0.5, 0.5),
                        runner = runner,
                        n_chains = 2)

# Tune the sampler
draws <- as_draws_df(samples)
vcv_tuned <- cov(draws[1:9])

# Sample
runner <- monty_runner_callr(6, progress = T)
sampler_tuned <- monty_sampler_adaptive(vcv_tuned)
samples_tuned <- monty_sample(model, sampler_tuned, 80000, 
                              initial = c(1e-05, 0.02002, 3e-05, 4e-05, 0.004, 
                                          0.5, 0.5, 0.5, 0.5),
                              runner = runner,
                              n_chains = 6)

# ----------- Diagnostics -----------------------------------------------------
# Check mixing
matplot(samples_tuned$density, type = "l", lty = 1,
        xlab = "Sample", ylab = "Posterior probability density")

# Thin and check mixing
samples_thinned <- monty_samples_thin(samples_tuned, burnin = 3000)
png(filename = paste0(path, "Trace log posterior.png"),
    width = 2200, height = 1500, res = 300)
matplot(samples_thinned$density, type = "l", lty = 1,
        xlab = "Iteration", ylab = "Posterior probability density",
        xaxt = "n",
        cex.lab = 1.4,
        cex.axis = 1.2,
        lwd = 2,
        col = 1:ncol(samples_thinned$density),  # optional colors
        panel.first = {
          rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4],
               col = "grey95")  # grey background
          abline(h = axTicks(2), col = "white")   # horizontal grid lines
          abline(v = axTicks(1), col = "white")   # vertical grid lines
        })

ticks <- axTicks(1)
axis(1,
     at = ticks,
     labels = format(ticks, big.mark = ","),
     cex.axis = 1.2)
dev.off ()

# Check trace
draws_thinned <- as_draws_df(samples_thinned) %>%
  rename_variables(
    "Spring component of the FOI" = "spring_comp",
    "Summer component of FOI" = "summer_comp",
    "Autumn component of the FOI" = "autumn_comp",
    "Winter component of the FOI" = "winter_comp",
    "μ" = "mu",
    "π" = "prop")

png(filename = paste0(path, "Trace parameters.png"),
    width = 2200, height = 3000, res = 300)
mcmc_trace(draws_thinned, facet_args = list(ncol = 1, strip.position = "top")) +
  theme(strip.text = element_text(size = 14),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) +
  scale_x_continuous(labels = scales::comma) +
  labs(y = "Posterior probability density\n", x = "\nIteration")
dev.off ()

# Summary
params_est <- summarise_draws(draws_thinned)
params_est %>% write.csv(paste0(path, "Parameter estimates.csv"))

draws_array_tuned <- as_draws_array(draws_thinned)
# Posterior uncertainty intervals
# Plot them separately because the values are very different
png(filename = paste0(path, "Uncertainty FOI season.png"),
    width = 2200, height = 1500, res = 300)
# Wrapped labels for display
param_labels_wrapped <- stringr::str_wrap(
  c("Spring component of the FOI",
    "Summer component of FOI",
    "Autumn component of the FOI",
    "Winter component of the FOI"),
  width = 16)

mcmc_intervals(draws_array_tuned,
               pars = rev(c("Spring component of the FOI", "Summer component of FOI",
                            "Autumn component of the FOI",  "Winter component of the FOI"))) +
  scale_y_discrete(labels = rev(param_labels_wrapped)) +  
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    axis.text.y = element_text(margin = margin(r = 6))) + 
  labs(y = "Parameter\n", x = "\nPosterior uncertainty interval")
dev.off ()

png(filename = paste0(path, "Uncertainty mu season.png"),
    width = 2200, height = 1500, res = 300)
mcmc_intervals(draws_array_tuned, pars = "μ") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) +
  labs(y = "Parameter\n", x = "\nPosterior uncertainty interval")
dev.off ()

png(filename = paste0(path, "Uncertainty prop season.png"),
    width = 2200, height = 1500, res = 300)
mcmc_intervals(draws_array_tuned, pars = "π") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)) +
  labs(y = "Parameter\n", x = "\nPosterior uncertainty interval")
dev.off ()

# Univariate marginal posterior distributions
png(filename = paste0(path, "Marginal posterior distributions season.png"),
    width = 2200, height = 1500, res = 300)
mcmc_hist(draws_array_tuned,
          facet_args = list(strip.position = "left",
                            labeller = labeller(
                              .variable = label_wrap_gen(width = 12)))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.title = element_text(size = 14),
        strip.text.x = element_text(size = 14,
                                    margin = margin(t = 4, b = 10)),
        panel.spacing  = unit(1.5, "lines")) +
  labs(y = "Posterior density\n") 
dev.off ()

# Trajectories
trajectories <- dust_unpack_state(filter,
                                  samples_thinned$observations$trajectories)
sero_conv <- array(trajectories$R_all, c(19, 25000))

# Plot the fit
matplot(incidence_data$time, sero_conv, type = "l", col = "#00000044", lty = 1,
        xlab = "Age (days)", ylab = "Proportion of seroconverted children")
points(x = incidence_data$time, y = incidence_data$prop_seroconv, pch = 19, col = "red")
arrows(x0 = incidence_data$time, y0 = incidence_data$seroprev_low95, 
       x1 = incidence_data$time, y1 = incidence_data$seroprev_up95, 
       angle = 90, code = 3, length = 0.1, col = "red")

dev.copy(jpeg,filename = paste0(path, "Model fit season.png"))
dev.off ()

# By season
# Spring
sero_conv_spring <- array(trajectories$R_sp, c(19, 1000))
png(filename = paste0(path, "Model fit season spring.png"),
    width = 2200, height = 1500, res = 300)
matplot(incidence_data_season_wide$time, sero_conv_spring, type = "l", col = "#00000044", lty = 1,
        xlab = "Age (days)", ylab = "Proportion of seroconverted children", main = "Spring cohort")
points(x = incidence_data_season_wide$time, y = incidence_data_season_wide$prop_seroconv_spring, pch = 19, col = "red")
arrows(x0 = incidence_data_season_wide$time, y0 = incidence_data_season_wide$seroprev_low95_spring, 
       x1 = incidence_data_season_wide$time, y1 = incidence_data_season_wide$seroprev_up95_spring, 
       angle = 90, code = 3, length = 0.1, col = "red")
legend("bottomright",  # or "topright", "topleft", etc.
       legend = c("Modelled seroconversion", "Observed seroconversion (95% CI)"),
       col = c("#00000044", "red"),
       lty = c(1, NA),
       pch = c(NA, 19),
       bty = "n")
dev.off ()

# Summer
sero_conv_summer <- array(trajectories$R_sm, c(19, 1000))
png(filename = paste0(path, "Model fit season summer.png"),
    width = 2200, height = 1500, res = 300)
matplot(incidence_data_season_wide$time, sero_conv_summer, type = "l", col = "#00000044", lty = 1,
        xlab = "Age (days)", ylab = "Proportion of seroconverted children", main = "Summer cohort")
points(x = incidence_data_season_wide$time, y = incidence_data_season_wide$prop_seroconv_summer, pch = 19, col = "red")
arrows(x0 = incidence_data_season_wide$time, y0 = incidence_data_season_wide$seroprev_low95_summer, 
       x1 = incidence_data_season_wide$time, y1 = incidence_data_season_wide$seroprev_up95_summer, 
       angle = 90, code = 3, length = 0.1, col = "red")
legend("bottomright",  # or "topright", "topleft", etc.
       legend = c("Modelled seroconversion", "Observed seroconversion (95% CI)"),
       col = c("#00000044", "red"),
       lty = c(1, NA),
       pch = c(NA, 19),
       bty = "n")
dev.off ()

# Autumn
sero_conv_autumn <- array(trajectories$R_au, c(19, 1000))
png(filename = paste0(path, "Model fit season autumn.png"),
    width = 2200, height = 1500, res = 300)
matplot(incidence_data_season_wide$time, sero_conv_autumn, type = "l", col = "#00000044", lty = 1,
        xlab = "Age (days)", ylab = "Proportion of seroconverted children", main = "Autumn cohort")
points(x = incidence_data_season_wide$time, y = incidence_data_season_wide$prop_seroconv_autumn, pch = 19, col = "red")
arrows(x0 = incidence_data_season_wide$time, y0 = incidence_data_season_wide$seroprev_low95_autumn, 
       x1 = incidence_data_season_wide$time, y1 = incidence_data_season_wide$seroprev_up95_autumn, 
       angle = 90, code = 3, length = 0.1, col = "red")
legend("bottomright",  # or "topright", "topleft", etc.
       legend = c("Modelled seroconversion", "Observed seroconversion (95% CI)"),
       col = c("#00000044", "red"),
       lty = c(1, NA),
       pch = c(NA, 19),
       bty = "n")
dev.off ()

# Winter
sero_conv_winter <- array(trajectories$R_wt, c(19, 1000))
png(filename = paste0(path, "Model fit season winter.png"),
    width = 2200, height = 1500, res = 300)
matplot(incidence_data_season_wide$time, sero_conv_winter, type = "l", col = "#00000044", lty = 1,
        xlab = "Age (days)", ylab = "Proportion of seroconverted children", main = "Winter cohort")
points(x = incidence_data_season_wide$time, y = incidence_data_season_wide$prop_seroconv_winter, pch = 19, col = "red")
arrows(x0 = incidence_data_season_wide$time, y0 = incidence_data_season_wide$seroprev_low95_winter, 
       x1 = incidence_data_season_wide$time, y1 = incidence_data_season_wide$seroprev_up95_winter, 
       angle = 90, code = 3, length = 0.1, col = "red")
legend("bottomright",  # or "topright", "topleft", etc.
       legend = c("Modelled seroconversion", "Observed seroconversion (95% CI)"),
       col = c("#00000044", "red"),
       lty = c(1, NA),
       pch = c(NA, 19),
       bty = "n")
dev.off ()