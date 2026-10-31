# Driveline PULSE Sensor Statistics & Data Analysis
# Six Sigma Black Belt Project

# Team Big Red Crew: Alexis Staveski, Ally Schwartz, Collin Shiang, Sophi Ramirez, Vedant Chavan
# SYSEN 5300, Cornell University
# Current Revision Date: 10/28/26

# This code will ... 
# - importing data and doing any post-processing + formatting, noting what the ideal location distance and angle is from the medial bony bone on the elbow for each thrower. also note the order in which treatments were applied (relevant for block randomization)
# - doing inital data visualization. plot a distribution of torque and distribution of arm speed as a histogram.
# - make control charts to visualize stability of the process. each thrower is a subgroup. make a averages chart, deviation chart, and ranges chart just to see. set up a framework for calculating the cp/cpk/pp/ppk values as applicable and calculating the confidence intervals. make a note that if confidence interval is uncertain, bootstrapping is needed. set up a framework for comparing them to the upper and lower control limits that we decide on.
# - perform the factorial design analysis. set up a framwork for determing direct effects of longitudinal deviation, rotational deviation, and combined deviation, each controlling for the effects of fatigue from the throw order. calculate confidence intervals on the dbar values where necessary.
# ^ rough descriptions - EDIT LATER

# ------------------------------------------------------------------------------


# Import Packages
# We load the tidyverse and statistical analysis packages used throughout SYSEN 5300.
library(dplyr)    # For data manipulation and piping (%>%)
library(readr)    # For reading CSV files
library(tidyr)    # For reshaping data (pivot_wider)
library(tibble)   # For tidy data frames
library(ggplot2)  # For visualization conforming to the course theme
library(broom)    # For tidying regression and test outputs (tidy, glance)
library(moments)  # For skewness and kurtosis
library(ggpubr)   # For arranging multi-panel plots (ggarrange)

# Source Helper Functions - see file for descriptions.
source("helper_functions.R")


# ------------------------------------------------------------------------------

# Import Experimental Data - data expected in CSV format
raw_data = read_csv("throwData_10.31.csv", col_types = cols())
# throw_data %>% head()

# thrower_ID is who is throwing
# throw_ID is the throw number (1-50 for most trials)
# sensor_placement indicates what treatment is being applied
# ideal_long_dist is the ideal distance from the medial epicondyle for that thrower (usually 4.5cm)
# arm_speed is in revolutions per minute (RPM)
# valgus_torque is in N*m
# longitudinal_pos is the distance from the medial epicondyle bone (from bone towards the hand)
# rotational_pos is the rotational angle from the ideal position (positive = CW looking from the shoulder to the hand)

# ideal placement is at ideal_long_dist, 0 deg rotation.

# Experimental Design Notes:
#   In this preliminary single-thrower test, treatments were conducted in sequential sets of 10 throws.
#   For the multi-thrower study, each Thrower will act as a BLOCK, with the 5 treatment sets 
#   randomized in order (5! = 120 possible permutations) to prevent fatigue confounding.


# ------------------------------------------------------------------------------
# Process Data
# ------------------------------------------------------------------------------

pulse_data = raw_data %>%
  # Select and order core measurement columns
  select(
    thrower_ID,
    throw_ID,
    arm_speed,
    valgus_torque,
    longitudinal_pos,
    rotational_pos,
    ideal_long_dist,
    sensor_placement
  ) %>%
  # Filter out bad data if it exists
  filter(!is.na(throw_ID), !is.na(valgus_torque), !is.na(arm_speed)) %>%
  mutate(
    # Identify thrower as a factor to simplify the label + to be used in blocking
    thrower = factor(paste("Thrower", thrower_ID)),

    # Calculate anatomical deviations from ideal baseline (ideal_long_dist and 0 deg rotation)
    longitudinal_dev = longitudinal_pos - ideal_long_dist,
    rotational_dev   = rotational_pos,

    # Convert sensor_placement to a factor ordered logically:
    # ideal -> long_displaced -> rot_90 -> rot_180 -> long_and_rot90
    sensor_placement = factor(
      sensor_placement,
      levels = c("ideal", "long_displaced", "rot_90", "rot_180", "long_and_rot90")
    ),

    # Combined condition key for grouping, matching sensor_placement order
    condition_group = factor(
      paste0(longitudinal_dev, "cm_", rotational_dev, "deg"),
      levels = c("0cm_0deg", "4.5cm_0deg", "0cm_90deg", "0cm_180deg", "4.5cm_90deg")
    )
  ) %>%

  # select certain columns to keep for actual data analysis
  select(
    thrower,
    sensor_placement,
    arm_speed,
    valgus_torque,
    longitudinal_dev,
    rotational_dev,
    condition_group
  )

# Data Ingestion Summary
cat("--- Data Ingestion Summary ---\n")
cat("Total valid throws imported:", nrow(pulse_data), "\n")
cat("Throwers included:", paste(levels(pulse_data$thrower), collapse = ", "), "\n")
cat("Sensor placement conditions:", paste(levels(pulse_data$sensor_placement), collapse = ", "), "\n\n")

# Treatment breakdown per thrower
# print(
#   pulse_data %>%
#     count(thrower, sensor_placement, condition_group, name = "throws_count")
# )

# look at organized data
glimpse(pulse_data)


# ------------------------------------------------------------------------------
# Visualize raw data with distributions
# ------------------------------------------------------------------------------

# compute mean, SD, skewness, kurtosis with helper function
torque_desc = pulse_data$valgus_torque %>% describe()
speed_desc  = pulse_data$arm_speed %>% describe()

# create plots:

# apply same ggplot theme across all plots
set_theme()

# histogram of peak valgus torque (all throws) for each thrower
plot_torque_hist = pulse_data %>%
  ggplot() +
  geom_histogram(
    mapping = aes(x = valgus_torque, fill = thrower),
    bins = 15,
    color = "white",
    alpha = 0.7,
    show.legend = FALSE
  ) +
  facet_wrap(~ thrower, scales = "free_y") +
  labs(
    title = "Raw Distribution of Valgus Torque (All Throws)",
    x = "Peak Valgus Torque (N*m)",
    y = "Count"
  )

# histogram of peak valgus torque (ideal placement) for each thrower
plot_torque_hist_ideal = pulse_data %>%
  filter(sensor_placement == "ideal") %>%
  ggplot() +
  geom_histogram(
    mapping = aes(x = valgus_torque, fill = thrower),
    bins = 15,
    color = "white",
    alpha = 0.7,
    show.legend = FALSE
  ) +
  facet_wrap(~ thrower, scales = "free_y") +
  labs(
    title = "Raw Distribution of Valgus Torque (Ideal Placement Throws)",
    x = "Peak Valgus Torque (N*m)",
    y = "Count"
  )

# histogram of arm speed for each thrower
plot_speed_hist = pulse_data %>%
  ggplot() +
  geom_histogram(
    mapping = aes(x = arm_speed, fill = thrower),
    bins = 15,
    color = "white",
    alpha = 0.7,
    show.legend = FALSE
  ) +
  facet_wrap(~ thrower, scales = "free_y") +
  labs(
    title = "Raw Distribution of Arm Speed (All Throws)",
    x = "Arm Speed (rpm)",
    y = "Count"
  )

# display plots
print(plot_torque_hist)
print(plot_torque_hist_ideal)
print(plot_speed_hist)

# ggsave("torque_histogram.png", plot_torque_hist)
# ggsave("torque_histogram_ideal.png", plot_torque_hist_ideal)
# ggsave("speed_histogram.png", plot_speed_hist)


# ------------------------------------------------------------------------------
# Table of average torques for each sensor placement
# ------------------------------------------------------------------------------

# compute average torque at ideal placement for each thrower
ideal_means = pulse_data %>%
  filter(sensor_placement == "ideal") %>%
  group_by(thrower) %>%
  summarize(mean_torque_ideal = mean(valgus_torque)) %>%
  ungroup()

# create summary table of average torque at each sensor placement
avg_torque_table = pulse_data %>%
  group_by(thrower, sensor_placement) %>%
  summarize(mean_torque = mean(valgus_torque)) %>%
  ungroup() %>%
  pivot_wider(names_from = sensor_placement, values_from = mean_torque)
# add bottom row providing the average of each column
avg_torque_summary = avg_torque_table %>%
  bind_rows(summarize(., thrower = "Average", across(where(is.numeric), ~ mean(.x, na.rm = TRUE)))) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))

# display summary table
print(avg_torque_summary)


# Table of percent deviation from ideal valgus torque

# calculate percent deviation for each throw relative to each thrower's ideal mean (add to main data frame)
pulse_data = pulse_data %>%
  left_join(ideal_means, by = "thrower") %>%  # add the ideal means as a column
  mutate(pct_deviation = (valgus_torque - mean_torque_ideal) / mean_torque_ideal * 100)   # calculate % deviation

# create summary table of average percent error at each sensor placement
pct_dev_table = pulse_data %>%
  filter(sensor_placement != "ideal") %>%
  group_by(thrower, sensor_placement) %>%
  summarize(mean_pct_error = mean(pct_deviation)) %>%
  ungroup() %>%
  pivot_wider(names_from = sensor_placement, values_from = mean_pct_error)

avg_row = pct_dev_table %>%  # provide the average of each column at the bottom
  summarize(
    thrower = "Average", across(where(is.numeric), ~ mean(.x, na.rm = TRUE))
  )

# append whole table
pct_dev_summary = bind_rows(pct_dev_table, avg_row) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2))) # round data

# display summary table
print(pct_dev_summary)

# boxplot with overlaid jitter of percent deviations across sensor placements
# plot_pct_dev_box_jitter = pulse_data %>%
#   ggplot() +
#   geom_hline(yintercept = 0, color = "darkgrey", linetype = "dashed", linewidth = 0.8) +
#   geom_boxplot(
#     mapping = aes(x = sensor_placement, y = pct_deviation, fill = thrower),
#     position = position_dodge(width = 0.8),
#     outlier.shape = NA,
#     alpha = 0.35,
#     width = 0.7
#   ) +
#   geom_point(
#     mapping = aes(x = sensor_placement, y = pct_deviation, color = thrower),
#     position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
#     size = 2,
#     alpha = 0.8
#   ) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   labs(
#     title = "Percent Deviation from Ideal Torque by Sensor Placement",
#     x = "Sensor Placement",
#     y = "Valgus Torque Percent Deviation (%)",
#     color = "Thrower",
#     fill  = "Thrower"
#   )
# # display plot
# print(plot_pct_dev_box_jitter)



# ------------------------------------------------------------------------------
# Factorial Design Analysis: Direct Effects & Interactions
# ------------------------------------------------------------------------------

# note: longitudinal displacement and rotational displacement at 180 was not tested for (would have completed the 2x3 design). so, 
# a 2x2 factorial experiment is performed and the rot_180 placement is viewed as simply an extra direct effect to look at.

# create factorial factor columns for 2x2 design (ideal, long_displaced, rot_90, long_and_rot90)
factorial_data_2x2 = pulse_data %>%
  filter(sensor_placement %in% c("ideal", "long_displaced", "rot_90", "long_and_rot90")) %>%
  mutate(  # make the factor levels (2 long., 2 rotational)
    long_factor = factor(
      if_else(longitudinal_dev == 0, "Ideal (0cm)", "Displaced (4.5cm)"),
      levels = c("Ideal (0cm)", "Displaced (4.5cm)")
    ),
    rot_factor = factor(
      if_else(rotational_dev == 0, "0deg", "90deg"),
      levels = c("0deg", "90deg")
    )
  )

# create separate rotational subset for comparing ideal vs 180 deg rotation
rot180_data = pulse_data %>%
  filter(sensor_placement %in% c("ideal", "rot_180")) %>%
  mutate(
    rot180_factor = factor(
      if_else(sensor_placement == "ideal", "0deg", "180deg"),
      levels = c("0deg", "180deg")
    )
  )

thrower_list = levels(pulse_data$thrower) # make a list of throwers for applying analysis to each

# Table 1: Direct effect of longitudinal displacement on valgus torque (raw N*m and percent deviation)
factorial_effects_long = lapply(c(thrower_list, "Overall"), function(th) {
  # handle the summary row case (average of all dbars) vs. actual data
  sub = if (th == "Overall") {
    factorial_data_2x2 %>% filter(rot_factor == "0deg")
  } else {
    factorial_data_2x2 %>% filter(thrower == th, rot_factor == "0deg")
  }
  mean_ideal = if (th == "Overall") {
    mean(pulse_data %>% filter(sensor_placement == "ideal") %>% pull(valgus_torque))
  } else {
    mean(pulse_data %>% filter(thrower == th, sensor_placement == "ideal") %>% pull(valgus_torque))
  }

  # calculate parameters using helper functions + CIs
  d   = dbar_oneway(valgus_torque ~ long_factor, data = sub)
  se  = se_factorial(valgus_torque ~ long_factor, data = sub)
  low = d - 1.96 * se
  upp = d + 1.96 * se

  # generate percent errors which are more useful
  d_pct   = (d / mean_ideal) * 100
  se_pct  = (se / mean_ideal) * 100
  low_pct = d_pct - 1.96 * se_pct
  upp_pct = d_pct + 1.96 * se_pct

  # package data
  tibble(
    thrower        = th,
    "dbar (N*m)"   = round(d, 2),
    "se (N*m)"     = round(se, 2),
    "95% CI (N*m)" = paste0("[", round(low, 2), ", ", round(upp, 2), "]"),
    "dbar (%)"     = round(d_pct, 2),
    "se (%)"       = round(se_pct, 2),
    "95% CI (%)"   = paste0("[", round(low_pct, 2), "%, ", round(upp_pct, 2), "%]")
  )
}) %>% bind_rows()

# Table 2: Direct effects of pure rotation (rot_90 and rot_180) in raw N*m and percent deviation
factorial_effects_rot = lapply(c(thrower_list, "Overall"), function(th) {
  # handle the summary row case (average of all dbars) vs. actual data
  mean_ideal = if (th == "Overall") {
    mean(pulse_data %>% filter(sensor_placement == "ideal") %>% pull(valgus_torque))
  } else {
    mean(pulse_data %>% filter(thrower == th, sensor_placement == "ideal") %>% pull(valgus_torque))
  }

  # 90 deg rotation vs ideal
  sub_90  = if (th == "Overall") {
    factorial_data_2x2 %>% filter(long_factor == "Ideal (0cm)")
  } else {
    factorial_data_2x2 %>% filter(thrower == th, long_factor == "Ideal (0cm)")
  }
  # calculate params + CIs for torque and percent deviation
  d_90       = dbar_oneway(valgus_torque ~ rot_factor, data = sub_90)
  se_90      = se_factorial(valgus_torque ~ rot_factor, data = sub_90)
  low_90     = d_90 - 1.96 * se_90
  upp_90     = d_90 + 1.96 * se_90
  d_90_pct   = (d_90 / mean_ideal) * 100
  se_90_pct  = (se_90 / mean_ideal) * 100
  low_90_pct = d_90_pct - 1.96 * se_90_pct
  upp_90_pct = d_90_pct + 1.96 * se_90_pct

  # 180 deg rotation vs ideal - separate case
  sub_180 = if (th == "Overall") {
    rot180_data
  } else {
    rot180_data %>% filter(thrower == th)
  }
  d_180       = dbar_oneway(valgus_torque ~ rot180_factor, data = sub_180)
  se_180      = se_factorial(valgus_torque ~ rot180_factor, data = sub_180)
  low_180     = d_180 - 1.96 * se_180
  upp_180     = d_180 + 1.96 * se_180
  d_180_pct   = (d_180 / mean_ideal) * 100
  se_180_pct  = (se_180 / mean_ideal) * 100
  low_180_pct = d_180_pct - 1.96 * se_180_pct
  upp_180_pct = d_180_pct + 1.96 * se_180_pct

  tibble( # package data
    thrower               = th,
    "dbar rot90 (N*m)"    = round(d_90, 2),
    "se rot90 (N*m)"      = round(se_90, 2),
    "95% CI rot90 (N*m)"  = paste0("[", round(low_90, 2), ", ", round(upp_90, 2), "]"),
    "dbar rot90 (%)"      = round(d_90_pct, 2),
    "se rot90 (%)"        = round(se_90_pct, 2),
    "95% CI rot90 (%)"    = paste0("[", round(low_90_pct, 2), "%, ", round(upp_90_pct, 2), "%]"),
    "dbar rot180 (N*m)"   = round(d_180, 2),
    "se rot180 (N*m)"     = round(se_180, 2),
    "95% CI rot180 (N*m)" = paste0("[", round(low_180, 2), ", ", round(upp_180, 2), "]"),
    "dbar rot180 (%)"     = round(d_180_pct, 2),
    "se rot180 (%)"       = round(se_180_pct, 2),
    "95% CI rot180 (%)"   = paste0("[", round(low_180_pct, 2), "%, ", round(upp_180_pct, 2), "%]")
  )
}) %>% bind_rows()

# Table 3: Interaction effects using lm() linear regression models (raw N*m and percent)
factorial_effects_interaction = lapply(c(thrower_list, "Overall (Blocked)"), function(th) {
  # handle the summary row case (average of all dbars) vs. actual data
  mean_ideal = if (th == "Overall (Blocked)") {
    mean(pulse_data %>% filter(sensor_placement == "ideal") %>% pull(valgus_torque))
  } else {
    mean(pulse_data %>% filter(thrower == th, sensor_placement == "ideal") %>% pull(valgus_torque))
  }
  m = if (th == "Overall (Blocked)") {
    lm(valgus_torque ~ long_factor * rot_factor + thrower, data = factorial_data_2x2)
  } else {
    lm(valgus_torque ~ long_factor * rot_factor, data = factorial_data_2x2 %>% filter(thrower == th))
  }
  td = tidy(m, conf.int = TRUE) %>%
    filter(term == "long_factorDisplaced (4.5cm):rot_factor90deg")

  # WAITING ON REMAINING COMBINATION TREATMENT DATA (E.G. THROWER 4): POPULATE NA IF MISSING CELL
  if (nrow(td) == 0 || is.na(td$estimate[1])) {
    return(tibble(
      thrower             = th,
      "interaction (N*m)" = NA_real_,
      "se (N*m)"          = NA_real_,
      "95% CI (N*m)"      = NA_character_,
      "interaction (%)"   = NA_real_,
      "se (%)"            = NA_real_,
      "95% CI (%)"        = NA_character_,
      p_value             = NA_real_
    ))
  }

  est_pct = (td$estimate / mean_ideal) * 100
  se_pct  = (td$std.error / mean_ideal) * 100
  low_pct = (td$conf.low / mean_ideal) * 100
  upp_pct = (td$conf.high / mean_ideal) * 100

  tibble(  # package data
    thrower             = th,
    "interaction (N*m)" = round(td$estimate, 2),
    "se (N*m)"          = round(td$std.error, 2),
    "95% CI (N*m)"      = paste0("[", round(td$conf.low, 2), ", ", round(td$conf.high, 2), "]"),
    "interaction (%)"   = round(est_pct, 2),
    "se (%)"            = round(se_pct, 2),
    "95% CI (%)"        = paste0("[", round(low_pct, 2), "%, ", round(upp_pct, 2), "%]"),
    p_value             = round(td$p.value, 4)
  )
}) %>% bind_rows()

# display summary tables
print(factorial_effects_long, width = Inf)
print(factorial_effects_rot, width = Inf)
print(factorial_effects_interaction, width = Inf)



# plot dbars with 95% confidence intervals for each treatment across throwers
plot_effects_data = bind_rows(
  lapply(thrower_list, function(th) {
    # longitudinal
    sub_l = factorial_data_2x2 %>% filter(thrower == th, rot_factor == "0deg")
    res_l = calc_direct_effects_table(
      tibble(factor = "Longitudinal (+4.5 cm)", dbar = dbar_oneway(valgus_torque ~ long_factor, data = sub_l)),
      se_pooled = se_factorial(valgus_torque ~ long_factor, data = sub_l)
    )
    # rot 90
    sub_r90 = factorial_data_2x2 %>% filter(thrower == th, long_factor == "Ideal (0cm)")
    res_r90 = calc_direct_effects_table(
      tibble(factor = "Rotational (90 deg)", dbar = dbar_oneway(valgus_torque ~ rot_factor, data = sub_r90)),
      se_pooled = se_factorial(valgus_torque ~ rot_factor, data = sub_r90)
    )
    # rot 180
    sub_r180 = rot180_data %>% filter(thrower == th)
    res_r180 = calc_direct_effects_table(
      tibble(factor = "Rotational (180 deg)", dbar = dbar_oneway(valgus_torque ~ rot180_factor, data = sub_r180)),
      se_pooled = se_factorial(valgus_torque ~ rot180_factor, data = sub_r180)
    )
    # combined (long + rot 90 vs ideal)
    # WAITING ON REMAINING DATA FOR COMBINATION TREATMENT IF MISSING
    sub_comb = pulse_data %>% filter(thrower == th, sensor_placement %in% c("ideal", "long_and_rot90")) %>%
      mutate(comb_f = factor(sensor_placement, levels = c("ideal", "long_and_rot90")))
    res_comb = if (n_distinct(sub_comb$comb_f) == 2 && sum(sub_comb$comb_f == "long_and_rot90") > 0) {
      calc_direct_effects_table(
        tibble(factor = "Combined (Long + Rot 90)", dbar = dbar_oneway(valgus_torque ~ comb_f, data = sub_comb)),
        se_pooled = se_factorial(valgus_torque ~ comb_f, data = sub_comb)
      )
    } else {
      tibble(factor = "Combined (Long + Rot 90)", dbar = NA_real_, se = NA_real_, lower_95 = NA_real_, upper_95 = NA_real_, label = NA_character_)
    }

    bind_rows(res_l, res_r90, res_r180, res_comb) %>%
      mutate(thrower = th)
  }) %>% bind_rows()
) %>%
  mutate(
    treatment_name = factor(
      factor,
      levels = c("Longitudinal (+4.5 cm)", "Rotational (90 deg)", "Rotational (180 deg)", "Combined (Long + Rot 90)")
    )
  )

plot_dbars_ci = plot_effects_data %>%
  ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgrey", linewidth = 0.8) +
  geom_errorbar(
    mapping = aes(x = treatment_name, ymin = lower_95, ymax = upper_95, color = thrower),
    position = position_dodge(width = 0.6),
    width = 0.2,
    linewidth = 0.8
  ) +
  geom_point(
    mapping = aes(x = treatment_name, y = dbar, color = thrower),
    position = position_dodge(width = 0.6),
    size = 3
  ) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  labs(
    title = "Direct Treatment Effects (dbar) on Valgus Torque with 95% CI",
    subtitle = "Differences relative to ideal placement (dashed line = 0 baseline)",
    x = "Treatment Condition",
    y = "dbar: Valgus Torque Difference (N*m)",
    color = "Thrower"
  )

# plot 2-way factorial interaction between longitudinal and rotational placement
interaction_plot_data = factorial_data_2x2 %>%
  group_by(thrower, long_factor, rot_factor) %>%
  summarize(
    mean_torque = mean(valgus_torque),
    se_torque   = sd(valgus_torque) / sqrt(n()),
    .groups     = "drop"
  )

plot_interaction = interaction_plot_data %>%
  ggplot() +
  geom_line(
    mapping = aes(x = long_factor, y = mean_torque, group = rot_factor, color = rot_factor),
    linewidth = 1
  ) +
  geom_point(
    mapping = aes(x = long_factor, y = mean_torque, color = rot_factor),
    size = 3
  ) +
  geom_errorbar(
    mapping = aes(x = long_factor, ymin = mean_torque - se_torque, ymax = mean_torque + se_torque, color = rot_factor),
    width = 0.15,
    linewidth = 0.7
  ) +
  facet_wrap(~ thrower) +
  labs(
    title = "Factorial Interaction Plot: Longitudinal vs. Rotational Placement",
    subtitle = "Error bars indicate +/- 1 SE | Non-parallel lines indicate interaction",
    x = "Longitudinal Placement",
    y = "Mean Peak Valgus Torque (N*m)",
    color = "Rotational\nPlacement"
  )

# display plots
print(plot_dbars_ci)
print(plot_interaction)



# ------------------------------------------------------------------------------
# Regression Analysis: Direct Effects & Interactions
# ------------------------------------------------------------------------------

