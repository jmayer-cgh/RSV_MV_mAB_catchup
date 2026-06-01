library(patchwork)
library(ggtext)

# Define palette
palette_season <- c(
  "Autumn" = "#E8743B", "Winter" = "#3B82F6",
  "Spring" = "#10B981", "Summer" = "#FFC300",
  "All" = "grey50"
)

palette_season_w_space <- c(
  "Autumn " = "#E8743B", "Winter " = "#3B82F6",
  "Spring " = "#10B981", "Summer " = "#FFC300",
  "All " = "#C0392B"
)

# --------------------- Figure 2 -----------------------------------------------
# Plot seroconversions with CI
converted_ci <- list()

for (index in seq(1:1000)) {
  converted_all_rand <- as.numeric(t(converted_all[index, ]))
  converted_sp_rand <- as.numeric(t(converted_sp[index, ]))
  converted_sm_rand <- as.numeric(t(converted_sm[index, ]))
  converted_au_rand <- as.numeric(t(converted_au[index, ]))
  converted_wt_rand <- as.numeric(t(converted_wt[index, ]))

  converted <- data.frame(
    age_midpoint = 0:1825,
    All = converted_all_rand,
    Spring = converted_sp_rand,
    Summer = converted_sm_rand,
    Autumn = converted_au_rand,
    Winter = converted_wt_rand
  )

  converted_ci[[index]] <- converted
  converted_ci[[index]]$iter <- index # save index in case we want to check one iteration
}

converted_ci_df <- do.call("rbind", converted_ci)

converted_ci_long <- converted_ci_df %>%
  pivot_longer(
    cols = !c(age_midpoint, iter),
    names_to = "season_birth",
    values_to = "seroconverted"
  ) %>%
  mutate(season_birth = paste0(season_birth, " ")) %>%
  mutate(season_birth = factor(season_birth, levels = c(
    "Autumn ", "Winter ",
    "Spring ", "Summer ",
    "All "
  )))

converted_ci <- converted_ci_long %>%
  group_by(season_birth, age_midpoint) %>%
  summarise(
    sero_low_95 = quantile(seroconverted, 0.025, na.rm = T),
    sero_median = quantile(seroconverted, 0.5, na.rm = T),
    sero_up_95 = quantile(seroconverted, 0.975, na.rm = T)
  ) %>%
  ungroup()

# Format seroprevalence data for plotting
incidence_data_season_plt <- incidence_data_season %>%
  filter(season_birth != "all" & xMidpoint <= 1826) %>%
  mutate(season_birth = paste0(str_to_sentence(season_birth), " ")) %>%
  mutate(season_birth = factor(season_birth,
    levels = c("Autumn ", "Winter ", "Spring ", "Summer ")
  ))

# Panel A
converted_plt <- converted_ci %>%
  filter(season_birth != "All ") %>%
  ggplot() +
  geom_line(aes(x = age_midpoint, y = sero_median, col = season_birth),
    linewidth = 1.4
  ) +
  geom_ribbon(
    aes(
      x = age_midpoint, ymin = sero_low_95, ymax = sero_up_95,
      fill = season_birth
    ),
    alpha = 0.2
  ) +
  geom_point(
    data = incidence_data_season_plt,
    aes(x = xMidpoint, y = prop_seroconv),
    col = "#8B1A1A", size = 4
  ) +
  geom_errorbar(
    data = incidence_data_season_plt,
    aes(x = xMidpoint, ymin = seroprev_low95, ymax = seroprev_up95),
    width = 45, col = "#8B1A1A", linewidth = 0.5
  ) +
  facet_wrap(~season_birth) +
  labs(
    x = "Age (months)",
    y = "Proportion seroconverted\n",
    colour = "Season of birth",
    fill = "Season of birth",
    tag = "A"
  ) +
  theme_light() +
  theme(
    axis.ticks.y = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 40),
    axis.title.x = element_text(size = 45),
    axis.text.y = element_text(size = 45),
    axis.title.y = element_text(size = 45),
    legend.text = element_text(size = 30),
    panel.grid.minor.y = element_blank(),
    strip.text.x = element_text(size = 45, color = "black"),
    strip.background = element_rect(fill = "white"),
    plot.tag = element_text(size = 50, face = "bold"),
    plot.tag.position = c(0.08, 0.99)
  ) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_x_continuous(
    breaks = seq(0, 1826, 365),
    labels = c("0", "12", "24", "36", "48", "60")
  ) +
  scale_colour_manual(values = palette_season_w_space) +
  scale_fill_manual(values = palette_season_w_space)

converted_plt

# Panel B
converted_plt_B <- converted_ci %>%
  filter(season_birth != "All " & age_midpoint <= 365 * 2) %>%
  ggplot() +
  geom_line(aes(x = age_midpoint, y = sero_median, col = season_birth),
    linewidth = 1.4
  ) +
  geom_ribbon(
    aes(
      x = age_midpoint, ymin = sero_low_95, ymax = sero_up_95,
      fill = season_birth
    ),
    alpha = 0.2
  ) +
  labs(
    x = "Age (months)",
    y = "\n\nProportion seroconverted\n",
    colour = "Season of birth",
    fill = "Season of birth",
    tag = "B"
  ) +
  theme_light() +
  theme(
    axis.ticks.y = element_blank(),
    legend.position = "bottom",
    axis.text.x = element_text(size = 45, hjust = c(1, 0.5, 0.6), vjust = c(1.1, 1, 1)),
    axis.title.x = element_text(size = 45),
    axis.text.y = element_text(size = 45),
    axis.title.y = element_text(size = 45),
    legend.title = element_text(size = 50, margin = margin(b = 20)),
    legend.text = element_text(size = 45),
    legend.key.height = unit(50, "pt"),
    legend.spacing.y = unit(10, "pt"),
    panel.grid.minor.y = element_blank(),
    strip.text.x = element_text(size = 33, color = "black"),
    strip.background = element_rect(fill = "white"),
    plot.tag = element_text(size = 50, face = "bold"),
    plot.tag.position = c(0.08, 0.99)
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25),
    labels = c("", "25%", "50%", "75%", "100%"),
    expand = c(0, 0)
  ) +
  scale_x_continuous(
    breaks = seq(0, 730, 365),
    labels = c("0 ", "12", "24 "),
    expand = c(0, 0)
  ) +
  guides(
    colour = guide_legend(title.position = "top", title.hjust = 0.5),
    fill = guide_legend(title.position = "top", title.hjust = 0.5)
  ) +
  scale_color_manual(values = palette_season_w_space) +
  scale_fill_manual(values = palette_season_w_space)

converted_plt_B

# --------------------- Combined Figure 2 (A + B) ------------------------------
# Combine both panels into one plot
combined_plot <- converted_plt + converted_plt_B + guide_area() +
  plot_layout(
    ncol = 2, guides = "collect", design = "AB\nCC",
    heights = c(10, 1)
  )

combined_plot

combined_plot %>% ggsave(
  filename = paste0(path_output, "Figure_2_combined.png"),
  width = 36, height = 20, units = "in",
  device = "png"
)


# --------------------- Figure 4 -----------------------------------------------
palette_intervention <- c(
  "No immunisation" = "#E91E8C",
  "MV" = "#00BCD4",
  "MV + mAB" = "#9B59B6"
)

# Plot incidence by age group
cum_incidence_plt <- incidence_outcome_monthly_cum_int %>%
  filter(season_birth != "All") %>%
  mutate(season_birth = factor(season_birth, levels = c("Autumn", "Winter", "Spring", "Summer"))) %>%
  ggplot(aes(
    x = age_bracket, y = cum_incidence_median,
    ymin = cum_incidence_low_95, ymax = cum_incidence_up_95,
    colour = intervention, fill = intervention
  )) +
  geom_line(linewidth = 1.5) +
  geom_ribbon(alpha = 0.2, color = NA) +
  labs(
    y = "**Hospitalisations** <br>per 100,000<br>",
    colour = "Intervention",
    fill = "Intervention",
    tag = "A"
  ) +
  facet_wrap(~season_birth, nrow = 1) +
  theme_light() +
  theme(
    axis.ticks.y = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 45),
    axis.text.y = element_text(size = 45),
    axis.title.x = element_blank(),
    axis.title.y = element_markdown(size = 45),
    legend.title = element_text(size = 45),
    legend.text = element_text(size = 40),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    strip.text.x = element_text(size = 45, color = "black"),
    strip.background = element_rect(fill = "white"),
    title = element_text(size = 45),
    plot.tag = element_text(size = 50, face = "bold"),
    plot.tag.position = c(0.08, 0.99)
  ) +
  scale_y_continuous(labels = scales::comma) +
  scale_x_continuous(breaks = seq(0, 12, by = 2)) +
  scale_color_manual(values = palette_intervention) +
  scale_fill_manual(values = palette_intervention)

cum_incidence_plt

cum_incidence_icu_plt <- incidence_outcome_monthly_cum_int %>%
  filter(season_birth != "All") %>%
  mutate(season_birth = factor(season_birth, levels = c("Autumn", "Winter", "Spring", "Summer"))) %>%
  ggplot(aes(
    x = age_bracket, y = cum_icu_inc_median,
    ymin = cum_icu_inc_low_95, ymax = cum_icu_inc_up_95,
    colour = intervention, fill = intervention
  )) +
  geom_line(linewidth = 1.5) +
  geom_ribbon(alpha = 0.2, color = NA) +
  labs(
    x = "\nAge (months)",
    y = "**ICU admissions** <br>per 100,000<br>",
    colour = "Intervention",
    fill = "Intervention",
    tag = "B"
  ) +
  facet_wrap(~season_birth, nrow = 1) +
  theme_light() +
  theme(
    axis.ticks.y = element_blank(),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 45),
    axis.title.x = element_text(size = 45),
    axis.text.y = element_text(size = 45),
    axis.title.y = element_markdown(size = 45),
    legend.title = element_text(size = 45, face = "bold"),
    legend.text = element_text(size = 45),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    strip.text.x = element_text(size = 45, color = "black"),
    strip.background = element_rect(fill = "white"),
    title = element_text(size = 45),
    plot.tag = element_text(size = 50, face = "bold"),
    plot.tag.position = c(0.08, 0.99)
  ) +
  guides(
    colour = guide_legend(title.position = "top", title.hjust = 0.5),
    fill = guide_legend(title.position = "top", title.hjust = 0.5)
  ) +
  scale_y_continuous(labels = scales::comma) +
  scale_x_continuous(breaks = seq(0, 12, by = 2)) +
  scale_color_manual(values = palette_intervention) +
  scale_fill_manual(values = palette_intervention)

cum_incidence_icu_plt

combined_cum_inc_plot <- cum_incidence_plt + cum_incidence_icu_plt +
  plot_layout(nrow = 2)

combined_cum_inc_plot

combined_cum_inc_plot %>% ggsave(
  filename = paste0(path_output, "Figure_4_combined.png"),
  width = 24, height = 16, units = "in",
  device = "png"
)

# Add stacked bar chart to this
plt <- total_outcome_intervention_plt_int %>%
  filter(season_birth != "All" &
    (intervention == "No immunisation" |
      intervention == "MV" |
      intervention == "+ mAB for all")) %>%
  ggplot() +
  geom_col(aes(x = intervention, y = hosp_median, fill = season_birth), position = "stack", alpha = 0.5) +
  geom_errorbar(
    data = subset(
      total_outcome_intervention_plt_int,
      season_birth == "All" &
        (intervention == "No immunisation" |
          intervention == "MV" |
          intervention == "+ mAB for all")
    ),
    aes(x = intervention, ymin = hosp_low_95, ymax = hosp_up_95),
    width = 0.2, position = "dodge"
  ) +
  labs(
    x = "\nIntervention",
    y = "<br>**Hospitalisations** <br>Total<br>",
    fill = "Season of birth",
    tag = "C"
  ) +
  theme_light() +
  theme(
    axis.ticks.y = element_blank(),
    legend.justification = "top",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 45),
    axis.text.y = element_text(size = 45),
    axis.title.x = element_blank(),
    axis.title.y = element_markdown(size = 45),
    legend.title = element_markdown(size = 45, face = "bold"),
    legend.text = element_text(size = 45),
    plot.tag = element_text(size = 50, face = "bold"),
    plot.tag.position = c(0.08, 0.99)
  ) +
  scale_fill_manual(values = palette_season) +
  scale_y_continuous(
    labels = scales::comma,
    breaks = scales::pretty_breaks(3)
  )

plt

plt_icu <- total_outcome_intervention_plt_int %>%
  filter(season_birth != "All" &
    (intervention == "No immunisation" |
      intervention == "MV" |
      intervention == "+ mAB for all")) %>%
  ggplot() +
  geom_col(aes(x = intervention, y = icu_median, fill = season_birth), position = "stack", alpha = 0.5) +
  geom_errorbar(
    data = subset(
      total_outcome_intervention_plt_int,
      season_birth == "All" &
        (intervention == "No immunisation" |
          intervention == "MV" |
          intervention == "+ mAB for all")
    ),
    aes(x = intervention, ymin = icu_low_95, ymax = icu_up_95),
    width = 0.2, position = "dodge"
  ) +
  labs(
    x = "\nIntervention",
    y = "<br>**ICU admissions** <br>Total<br>",
    fill = "Season of birth",
    tag = "D"
  ) +
  theme_light() +
  theme(
    axis.ticks.y = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 45),
    axis.text.y = element_text(size = 45),
    axis.title.x = element_text(size = 45),
    axis.title.y = element_markdown(size = 45),
    legend.title = element_markdown(size = 45, face = "bold"),
    legend.text = element_text(size = 45),
    plot.tag = element_text(size = 50, face = "bold"),
    plot.tag.position = c(0.08, 0.99)
  ) +
  scale_fill_manual(values = palette_season) +
  scale_y_continuous(
    labels = scales::comma,
    breaks = scales::pretty_breaks(3)
  )

plt_icu

combined_inc_plot <- plt + plt_icu +
  plot_layout(nrow = 2)

combined_inc_plot

figure_4_combined <- (cum_incidence_plt + plt) / (cum_incidence_icu_plt + plt_icu) +
  plot_layout(widths = c(30, 1))

figure_4_combined

left_col <- wrap_plots(cum_incidence_plt, cum_incidence_icu_plt, ncol = 1)
right_col <- wrap_plots(plt, plt_icu, ncol = 1)

figure_4_combined <- wrap_plots(left_col, right_col, ncol = 2, widths = c(4, 1))

figure_4_combined %>% ggsave(
  filename = paste0(path_output, "Figure 4 combined.png"),
  width = 46, height = 28, units = "in",
  device = "png"
)

# --------------------- Figure 3 -----------------------------------------------
incidence_3_outcomes <- incidence_infection_CI %>%
  filter(season_birth != "all") %>%
  mutate(
    outcome = "Infection",
    season_birth = str_to_sentence(season_birth)
  ) %>%
  rename(
    incidence_up_95 = incidence_upper,
    incidence_low_95 = incidence_lower
  ) %>%
  rbind(
    incidence_outcome_int %>%
      filter(intervention == "No immunisation" & season_birth != "All") %>%
      select(season_birth, age_bracket, incidence_up_95, incidence_median, incidence_low_95) %>%
      mutate(outcome = "Hospitalisation")
  ) %>%
  rbind(
    incidence_outcome_int %>%
      filter(intervention == "No immunisation" & season_birth != "All") %>%
      select(season_birth, age_bracket, icu_inc_up_95, icu_inc_median, icu_inc_low_95) %>%
      rename(
        incidence_up_95 = icu_inc_up_95,
        incidence_median = icu_inc_median,
        incidence_low_95 = icu_inc_low_95
      ) %>%
      mutate(outcome = "ICU admission")
  ) %>%
  filter(age_bracket != "all") %>%
  mutate(
    age_bracket = as.numeric(age_bracket),
    season_birth = factor(season_birth, levels = c("Spring", "Summer", "Autumn", "Winter", "All")),
    outcome = factor(outcome, levels = c("Infection", "Hospitalisation", "ICU admission"))
  )

# Plot on a panel
incidence_plt <- incidence_3_outcomes %>%
  ggplot(aes(x = age_bracket, y = incidence_median)) +
  geom_line(aes(col = season_birth), linewidth = 1.5) +
  geom_ribbon(aes(ymin = incidence_low_95, ymax = incidence_up_95, fill = season_birth), alpha = 0.15) +
  facet_wrap(~outcome, scales = "free") +
  labs(
    x = "\nAge (months)\n\n", y = "Incidence per 100,000 children\n",
    colour = "Season of birth",
    fill = "Season of birth",
    tag = "A"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 45),
    axis.text.y = element_text(size = 45),
    axis.title.x = element_text(size = 45),
    axis.title.y = element_text(size = 45),
    legend.title = element_text(size = 45, face = "bold"),
    legend.text = element_text(size = 45),
    strip.text = element_text(size = 40, face = "bold", color = "black"),
    plot.tag = element_text(size = 50, face = "bold"),
    plot.tag.position = c(0.08, 0.99)
  ) +
  scale_x_continuous(breaks = seq(0, 12, by = 2)) +
  scale_y_continuous(labels = scales::comma) +
  scale_color_manual(values = palette_season) +
  scale_fill_manual(values = palette_season) +
  ggh4x::facetted_pos_scales(y = list(
    Infection = scale_y_continuous(limits = c(0, 35000), labels = scales::comma),
    Hospitalisation = scale_y_continuous(limits = c(0, 35), labels = scales::comma),
    `ICU admission` = scale_y_continuous(limits = c(0, 7), labels = scales::comma)
  ))

incidence_plt

#  Disease progression
infection_to_disease_plt <- infection_to_disease_prop %>%
  filter(age_bracket != "all" & season_birth == "All" &
    (intervention == "No immunisation")) %>%
  mutate(
    age_bracket = as.numeric(age_bracket),
    season_birth = factor(season_birth,
      levels = c("Spring", "Summer", "Autumn", "Winter", "All")
    )
  ) %>%
  ggplot(aes(x = age_bracket)) +
  geom_line(aes(y = prop_hosp_median, colour = season_birth), linewidth = 1.5) +
  geom_ribbon(aes(ymin = prop_hosp_low_95, ymax = prop_hosp_up_95, fill = season_birth),
    alpha = 0.2
  ) +
  labs(
    x = "\nAge (months)",
    y = "Proportion",
    title = "Infections leading to hospitalisation\n",
    tag = "B"
  ) +
  theme_light() +
  theme(
    axis.ticks.y = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 45),
    axis.title.x = element_text(size = 45),
    axis.text.y = element_text(size = 45),
    axis.title.y = element_blank(),
    plot.title = element_text(size = 45, hjust = 0.5, face = "bold"),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    strip.text.x = element_text(size = 45, color = "black"),
    plot.title.position = "plot",
    strip.background = element_rect(fill = "white"),
    plot.tag = element_text(size = 50, face = "bold"),
    plot.tag.position = c(0.08, 0.99)
  ) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_x_continuous(breaks = seq(0, 12, by = 2)) +
  scale_color_manual(values = palette_season) +
  scale_fill_manual(values = palette_season)

infection_to_disease_plt


hosp_to_icu_plt <- infection_to_disease_prop %>%
  filter(age_bracket != "all" & season_birth == "All" &
    (intervention == "No immunisation")) %>%
  mutate(
    age_bracket = as.numeric(age_bracket),
    season_birth = factor(season_birth,
      levels = c("Spring", "Summer", "Autumn", "Winter", "All")
    )
  ) %>%
  ggplot(aes(x = age_bracket)) +
  geom_line(aes(y = prop_icu_median, colour = season_birth), linewidth = 1.5) +
  geom_ribbon(aes(ymin = prop_icu_low_95, ymax = prop_icu_up_95, fill = season_birth),
    alpha = 0.2
  ) +
  labs(
    x = "\nAge (months)",
    y = "\n",
    title = "Hospitalisations leading to ICU admission\n",
    color = "Season of birth",
    fill = "Season of birth"
  ) +
  theme_light() +
  theme(
    axis.ticks.y = element_blank(),
    legend.position = "right",
    legend.title = element_markdown(size = 45, face = "bold" ),
    legend.text = element_text(size = 45),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 45),
    axis.title.x = element_text(size = 45),
    axis.text.y = element_text(size = 45),
    axis.title.y = element_text(size = 45),
    plot.title = element_text(size = 45, hjust = 0.5, face = "bold"),
    plot.title.position = "plot",
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    strip.text.x = element_text(size = 45, color = "black"),
    strip.background = element_rect(fill = "white")
  ) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_x_continuous(breaks = seq(0, 12, by = 2)) +
  scale_color_manual(values = palette_season) +
  scale_fill_manual(values = palette_season)

hosp_to_icu_plt

infection_to_disease_prop_wide %>%
  filter(age_bracket != "all" & season_birth == "All" &
    (intervention == "No immunisation" | intervention == "MV" | intervention == "+ mAB for all")) %>%
  mutate(
    age_bracket = as.numeric(age_bracket),
    intervention = case_when(
      intervention == "+ mAB for all" ~ "MV + mAB",
      T ~ intervention
    ),
    intervention = factor(intervention, levels = c("No immunisation", "MV", "MV + mAB"))
  ) %>%
  ggplot(aes(x = age_bracket)) +
  geom_line(aes(y = prop_median, colour = intervention)) +
  geom_ribbon(aes(ymin = prop_low_95, ymax = prop_up_95, fill = intervention),
    alpha = 0.2
  ) +
  labs(
    x = "\nAge (months)",
    color = "Intervention",
    fill = "Intervention"
  ) +
  facet_wrap(~progression) +
  theme_light() +
  theme(
    axis.ticks.y = element_blank(),
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
    axis.title.x = element_text(size = 25),
    axis.text.y = element_text(size = 20),
    axis.title.y = element_blank(),
    title = element_text(size = 25),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    strip.text.x = element_text(size = 25, color = "black"),
    strip.background = element_rect(fill = "white")
  ) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_x_continuous(breaks = seq(0, 12, by = 1)) +
  scale_color_manual(values = palette_intervention) +
  scale_fill_manual(values = palette_intervention)

# Combine for figure 3
figure_3_combined <- incidence_plt / (infection_to_disease_plt + hosp_to_icu_plt) +
  plot_layout(heights = c(2, 1))

figure_3_combined %>% ggsave(
  filename = paste0(path_output, "Figure_3_combined.png"),
  width = 38, height = 32, units = "in",
  device = "png"
)

# Try using a bar chart instead
disease_prog_plt <- infection_to_disease_prop_wide %>%
  filter(age_bracket == "all" & season_birth != "All" &
    (intervention == "No immunisation" | intervention == "MV" | intervention == "+ mAB for all")) %>%
  mutate(
    intervention = case_when(
      intervention == "+ mAB for all" ~ "MV + mAB",
      T ~ intervention
    ),
    intervention = factor(intervention, levels = c("No immunisation", "MV", "MV + mAB"))
  ) %>%
  ggplot(aes(x = intervention, y = prop_median, fill = season_birth)) +
  geom_col(position = "stack") +
  labs(
    x = "\nAge (months)",
    y = "Disease progression\n",
    fill = "Season of birth"
  ) +
  facet_wrap(~progression, scale = "free") +
  theme_light() +
  theme(
    axis.ticks.y = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 25),
    legend.text = element_text(size = 20),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 20),
    axis.title.x = element_text(size = 25),
    axis.text.y = element_text(size = 20),
    axis.title.y = element_text(size = 25),
    strip.text.x = element_text(size = 25, color = "black"),
    strip.background = element_rect(fill = "white"),
    title = element_text(size = 25)
  ) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = palette_season)

disease_prog_plt

# --------------------- Plot NNVs ----------------------------------
nnv_data <- nnv_hosp %>% 
  filter((intervention == "MV" & season_birth == "all")) %>%
  mutate(outcome = "Hospitalisation") %>%
  rbind(
    nnv_icu %>%
      filter((intervention == "MV" & season_birth == "all")) %>%
      mutate(outcome = "ICU admission"))

# Format NNV for mAB compared to MV
nnv_mAB <- nnv_mAB_MV %>%
  select(season_birth, NNV_hosp_median, NNV_hosp_low95, NNV_hosp_up95) %>%
  rename(
    NNV_low95 = NNV_hosp_low95,
    NNV_up95 = NNV_hosp_up95,
    NNV_median = NNV_hosp_median
  ) %>%
  mutate(outcome = "Hospitalisation") %>%
  rbind(
    nnv_mAB_MV %>%
      select(season_birth, NNV_icu_median, NNV_icu_low95, NNV_icu_up95) %>%
      rename(
        NNV_low95 = NNV_icu_low95,
        NNV_up95 = NNV_icu_up95,
        NNV_median = NNV_icu_median
      ) %>%
      mutate(outcome = "ICU admission")
  ) %>%
  # Turn negative values into 0
  mutate(
    NNV_low95 = ifelse(NNV_low95 < 0, 0, NNV_low95),
    NNV_up95 = ifelse(NNV_up95 < 0, 0, NNV_up95),
    NNV_median = ifelse(NNV_median < 0, 0, NNV_median)
  )

# Bind both and order
nnv_data <- nnv_data %>%
  rbind(nnv_mAB %>% mutate(intervention = "MV + mAB")) %>%
  mutate(season_birth = str_to_sentence(season_birth)) %>%
  mutate(season_birth = factor(season_birth, 
                               levels = c("Spring", "Summer", "Autumn", 
                                          "Winter", "All"))) %>%
  mutate(intervention = factor(intervention, 
                               levels = c("MV", "MV + mAB")))

# Plot
nnv_plot <- nnv_data %>%
  ggplot(aes(x = intervention, y = NNV_median, fill = season_birth)) +
  geom_col(position = "dodge", alpha = 0.5) +
  geom_errorbar(aes(ymin = NNV_low95, ymax = NNV_up95), width = 0.2, 
                position = position_dodge(0.9)) +
  facet_wrap(~outcome, scales = "free") +
  labs(
    x = "\nIntervention",
    y = "NNV\n",
    fill = "Season of birth"
  ) +
  theme_light() +
  theme(
    axis.ticks.y = element_blank(),
    legend.position = "right",
    legend.title = element_markdown(size = 30, face = "bold"),
    legend.text = element_text(size = 30),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 35),
    axis.title.x = element_text(size = 35),
    axis.text.y = element_text(size = 35),
    axis.title.y = element_text(size = 35),
    strip.text.x = element_text(size = 35, color = "black"),
    strip.background = element_rect(fill = "white")
  ) +
  scale_y_continuous(labels = scales::comma) +
  scale_fill_manual(values = palette_season) +
  scale_y_log10() 


nnv_plot
nnv_plot %>% ggsave(
  filename = paste0(path_output, "NNV plot.png"),
  width = 18, height = 16, units = "in",
  device = "png"
)
