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
library(tibble)   # For tidy data frames
library(ggplot2)  # For visualization conforming to the course theme
library(broom)    # For tidying regression and test outputs (tidy, glance)
library(moments)  # For skewness and kurtosis
library(ggpubr)   # For arranging multi-panel plots (ggarrange)

# Source Helper Functions - see file for descriptions.
source("helper_functions.R")

# Apply course ggplot theme across all visualizations
set_theme()


# ------------------------------------------------------------------------------

# Import Experimental Data - data expected in CSV format
raw_data = read_csv("throwData_10.29.csv", col_types = cols())
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

    # add sensor_placement labels directly if any entries are missing in raw CSV    FIX THIS LATER TO BE JUST IN THE CSV
    sensor_placement = case_when(
      !is.na(sensor_placement) & sensor_placement != "" ~ sensor_placement,
      longitudinal_dev == 0 & rotational_pos == 0       ~ "ideal",
      rotational_pos == 90                              ~ "rot_90",
      rotational_pos == 180                             ~ "rot_180",
      longitudinal_dev == 4.5                           ~ "long_9cm",
      longitudinal_dev == 9                             ~ "long_13.5cm",
      longitudinal_dev == -9                            ~ "long_-4.5cm",
      longitudinal_dev == -13.5                         ~ "long_-9cm",
      TRUE                                              ~ paste0("pos_", longitudinal_pos)
    ),

    # Convert sensor_placement to a factor with "ideal" as the reference (first) level
    sensor_placement = factor(
      sensor_placement,
      levels = c("ideal", setdiff(unique(sensor_placement), "ideal"))
    ),

    # Combined condition key for grouping
    condition_group = paste0(longitudinal_dev, "cm_", rotational_dev, "deg")
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

# compute mean, SD, skewness, kurtosis with helper function
torque_desc = pulse_data$valgus_torque %>% describe()
speed_desc  = pulse_data$arm_speed %>% describe()

# create plots:

# 1. Histogram for each person's torque
plot_torque_hist_ind = pulse_data %>%
  ggplot() +
  geom_histogram(mapping = aes(x = valgus_torque, fill = thrower), bins = 15, color = "white", alpha = 0.7, show.legend = FALSE) +
  facet_wrap(~ thrower, scales = "free_y") +
  labs(
    title = "Distribution of Peak Valgus Torque by Thrower",
    x = "Peak Valgus Torque (N*m)",
    y = "Count"
  )

# 2. Scatter plot for each person's torque across sensor placements
plot_torque_scatter_ind = pulse_data %>%
  ggplot() +
  geom_jitter(mapping = aes(x = sensor_placement, y = valgus_torque, color = thrower), width = 0.2, height = 0, size = 2, alpha = 0.7, show.legend = FALSE) +
  facet_wrap(~ thrower, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Peak Valgus Torque by Sensor Placement for Each Thrower",
    x = "Sensor Placement",
    y = "Peak Valgus Torque (N*m)"
  )

# 3. Histogram for each person's arm speed
plot_speed_hist_ind = pulse_data %>%
  ggplot() +
  geom_histogram(mapping = aes(x = arm_speed, fill = thrower), bins = 15, color = "white", alpha = 0.7, show.legend = FALSE) +
  facet_wrap(~ thrower, scales = "free_y") +
  labs(
    title = "Distribution of Arm Speed by Thrower",
    x = "Arm Speed (rpm)",
    y = "Count"
  )

# 4. Scatter plot for all throwers' torque across sensor placements
plot_torque_scatter_all = pulse_data %>%
  ggplot() +
  geom_jitter(mapping = aes(x = sensor_placement, y = valgus_torque, color = thrower), width = 0.2, height = 0, size = 2.5, alpha = 0.8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Peak Valgus Torque across Sensor Placements (All Throwers)",
    x = "Sensor Placement",
    y = "Peak Valgus Torque (N*m)",
    color = "Thrower"
  )

# 5. Scatter plot for all throwers' arm speed across sensor placements
plot_speed_scatter_all = pulse_data %>%
  ggplot() +
  geom_jitter(mapping = aes(x = sensor_placement, y = arm_speed, color = thrower), width = 0.2, height = 0, size = 2.5, alpha = 0.8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Arm Speed across Sensor Placements (All Throwers)",
    x = "Sensor Placement",
    y = "Arm Speed (rpm)",
    color = "Thrower"
  )

# Display plots
print(plot_torque_hist_ind)
print(plot_torque_scatter_ind)
print(plot_speed_hist_ind)
print(plot_torque_scatter_all)
print(plot_speed_scatter_all)


# 4. Statistical Process Control (SPC) & Control Charts
# ------------------------------------------------------------------------------
# Course SPC concepts:
# - In the final multi-thrower experiment, each THROWER is a natural subgroup.
# - For this preliminary test set (single thrower), we treat each 10-throw condition SET
#   as a subgroup to examine stability across placements.
# - sigma_short (sigma_s) measures within-subgroup common cause variation: sqrt(mean(sd^2)).
# - sigma_total (sigma_t) measures overall variation across the entire dataset.

stat_subgroups = pulse_data %>%
  group_by(set_order) %>%
  summarize(
    condition = first(condition_group),
    xbar      = mean(valgus_torque),
    r         = max(valgus_torque) - min(valgus_torque),
    sd        = sd(valgus_torque),
    nw        = n(),
    df        = nw - 1,
    .groups   = "drop"
  ) %>%
  mutate(
    # Short-term variation (pooled within-group variation)
    sigma_s = sqrt(sum(df * sd^2) / sum(df)),
    # Standard error for 6-sigma control limits
    se      = sigma_s / sqrt(nw),
    # Grand mean (xbbar)
    xbbar   = mean(xbar),
    # Control limits for Xbar chart (+/- 3 standard errors)
    ucl_x   = xbbar + 3 * se,
    lcl_x   = xbbar - 3 * se,
    # S chart limits (using course B3/B4 constants approximation via chi/normal simulation or formula)
    sbar    = mean(sd),
    ucl_s   = sbar + 3 * (sigma_s / sqrt(2 * first(df))),
    lcl_s   = pmax(0, sbar - 3 * (sigma_s / sqrt(2 * first(df)))),
    # R chart limits
    rbar    = mean(r),
    ucl_r   = rbar + 3 * sd(r),
    lcl_r   = pmax(0, rbar - 3 * sd(r))
  )

sigma_s_val = unique(stat_subgroups$sigma_s)
sigma_t_val = sd(pulse_data$valgus_torque)
xbbar_val   = unique(stat_subgroups$xbbar)

# 4A. Average (Xbar) Chart
labels_xbar = stat_subgroups %>%
  reframe(
    set_order = max(set_order),
    name      = c("UCL (+3s)", "Grand Mean", "LCL (-3s)"),
    val       = c(first(ucl_x), first(xbbar), first(lcl_x)),
    label     = paste(name, " = ", round(val, 2), sep = "")
  )

chart_xbar = stat_subgroups %>%
  ggplot(mapping = aes(x = set_order, y = xbar)) +
  geom_hline(aes(yintercept = xbbar), color = "lightgrey", size = 2) +
  geom_ribbon(aes(ymin = lcl_x, ymax = ucl_x), fill = "steelblue", alpha = 0.15) +
  geom_line(color = "steelblue", size = 1) +
  geom_point(color = "steelblue", size = 4) +
  geom_label(data = labels_xbar, aes(x = set_order, y = val, label = label), hjust = 1) +
  labs(
    title = "Xbar Control Chart: Average Valgus Torque",
    subtitle = "Subgroups = 10-Throw Sets | Limits = Grand Mean +/- 3 SE",
    x = "Subgroup (Set Order)",
    y = "Subgroup Mean Torque (N*m)"
  )
print(chart_xbar)

# 4B. Standard Deviation (S) Chart
chart_s = stat_subgroups %>%
  ggplot(mapping = aes(x = set_order, y = sd)) +
  geom_hline(aes(yintercept = sbar), color = "lightgrey", size = 2) +
  geom_ribbon(aes(ymin = lcl_s, ymax = ucl_s), fill = "forestgreen", alpha = 0.15) +
  geom_line(color = "forestgreen", size = 1) +
  geom_point(color = "forestgreen", size = 4) +
  labs(
    title = "S Control Chart: Within-Subgroup Standard Deviation",
    subtitle = "Tracking process variation across placement conditions",
    x = "Subgroup (Set Order)",
    y = "Subgroup SD (N*m)"
  )
print(chart_s)

# 4C. Range (R) Chart
chart_r = stat_subgroups %>%
  ggplot(mapping = aes(x = set_order, y = r)) +
  geom_hline(aes(yintercept = rbar), color = "lightgrey", size = 2) +
  geom_ribbon(aes(ymin = lcl_r, ymax = ucl_r), fill = "darkorange", alpha = 0.15) +
  geom_line(color = "darkorange", size = 1) +
  geom_point(color = "darkorange", size = 4) +
  labs(
    title = "R Control Chart: Within-Subgroup Range",
    subtitle = "Tracking max - min spread within each 10-throw set",
    x = "Subgroup (Set Order)",
    y = "Subgroup Range (N*m)"
  )
print(chart_r)


# 5. Process Capability & Performance Framework (Cp, Cpk, Pp, Ppk)
# ------------------------------------------------------------------------------
# Distinction from class:
# - Control Limits come from the PROCESS itself (empirical variation).
# - Specification Limits come from the CUSTOMER / BIOMECHANICS REQUIREMENTS.
#
# Placeholder Spec Limits (To be finalized by project team):
# E.g., based on literature and baseline ideal placement (mean ~33 Nm),
# an acceptable tolerance band might be +/- 6 Nm from baseline (27 Nm to 39 Nm).
spec_lower = 27.0
spec_upper = 39.0

# Calculate point estimates for the four indices (cp, pp, cpk, ppk from helper_functions.R)
est_cp  = cp(sigma_s = sigma_s_val, upper = spec_upper, lower = spec_lower)
est_pp  = pp(sigma_t = sigma_t_val, upper = spec_upper, lower = spec_lower)
est_cpk = cpk(mu = xbbar_val, sigma_s = sigma_s_val, lower = spec_lower, upper = spec_upper)
est_ppk = ppk(mu = xbbar_val, sigma_t = sigma_t_val, lower = spec_lower, upper = spec_upper)

# Theoretical 95% Confidence Intervals for Capability Indices
# (Computed via calc_capability_summary() from helper_functions.R)
capability_summary = calc_capability_summary(
  mu            = xbbar_val,
  sigma_s       = sigma_s_val,
  sigma_t       = sigma_t_val,
  lower         = spec_lower,
  upper         = spec_upper,
  n_subgroups   = nrow(stat_subgroups),
  subgroup_size = first(stat_subgroups$nw),
  n_total       = nrow(pulse_data)
)
print(capability_summary)

# NOTE ON UNCERTAINTY & BOOTSTRAPPING:
# If theoretical confidence intervals are wide, straddle 1.0, or if the underlying
# normality assumption is violated (e.g., severe skewness in rotated placements),
# BOOTSTRAPPING with replacement is required to obtain empirical percentiles.
#
# Available from helper_functions.R:
# get_index(x, y, index, upper, lower, bootstrap_reps = 1000, ci_level = 0.95, by_subgroup = TRUE)
#
# Example bootstrap evaluation for Cpk:
# boot_cpk = get_index(
#   x = pulse_data$set_order,
#   y = pulse_data$valgus_torque,
#   index = "cpk",
#   upper = spec_upper,
#   lower = spec_lower,
#   bootstrap_reps = 1000
# )
# print(boot_cpk)


# 6. Factorial Design Analysis: Direct Effects & Interactions
# ------------------------------------------------------------------------------
# 6A. Group-Means Summary of Factorial Conditions
# Calculate mean (xbar), sd (s), and count (n) per condition group
groups = pulse_data %>%
  group_by(longitudinal_dev_f, rotational_pos_f) %>%
  summarize(
    xbar = mean(valgus_torque),
    s    = sd(valgus_torque),
    n    = n(),
    .groups = "drop"
  )
print(groups)

# 6B. Calculate Differences of Grand Means (d-double-bar) for Direct Effects
# Reference baseline: Longitudinal Deviation = 0 cm, Rotational Position = 0 deg
baseline_xbar = groups %>%
  filter(longitudinal_dev_f == "0" & rotational_pos_f == "0") %>%
  pull(xbar)

# Direct effect of longitudinal deviations (isolated against ideal baseline 0 cm)
# Direct effect of rotational deviations (isolated against ideal baseline 0 deg)
dbbar_direct = tibble(
  factor = c(
    "Longitudinal: +4.5 cm vs 0 cm",
    "Longitudinal: +9.0 cm vs 0 cm",
    "Longitudinal: -9.0 cm vs 0 cm",
    "Longitudinal: -13.5 cm vs 0 cm",
    "Rotational: 90 deg vs 0 deg",
    "Rotational: 180 deg vs 0 deg"
  ),
  dbar = c(
    mean(groups$xbar[groups$longitudinal_dev_f == "4.5"]) - baseline_xbar,
    mean(groups$xbar[groups$longitudinal_dev_f == "9"])   - baseline_xbar,
    mean(groups$xbar[groups$longitudinal_dev_f == "-9"])  - baseline_xbar,
    mean(groups$xbar[groups$longitudinal_dev_f == "-13.5"]) - baseline_xbar,
    mean(groups$xbar[groups$rotational_pos_f == "90"])   - baseline_xbar,
    mean(groups$xbar[groups$rotational_pos_f == "180"])  - baseline_xbar
  )
)

# 6C. Pooled Standard Error Across Factorial Treatments
# Following Week 9 LC 01/02 formula:
error_pooled = groups %>%
  mutate(df = n - 1) %>%
  summarize(
    var_pooled = sum(s^2 * df) / sum(df),
    s_pooled   = sqrt(var_pooled),
    se         = sqrt(s_pooled^2 * sum(1 / n))
  )

# Direct Effects Table with 95% Confidence Intervals
# (Computed using calc_direct_effects_table() from helper_functions.R)
direct_effects = calc_direct_effects_table(
  dbar_df   = dbbar_direct,
  se_pooled = error_pooled$se
)
print(direct_effects)

# Note: For standard 2^k factorial designs, functions from helper_functions.R can be used:
# - se_factorial(formula, data) to estimate pooled SE across factorial groups
# - dbar_oneway(), dbar_twoway(), dbar_threeway() for main and interaction effects

# 6D. Linear Model Controlling for Fatigue (Throw Sequence / Block)
# In SYSEN 5300, lm() models both direct effects and two-way interaction (* operator),
# while controlling for throw_ID (order) to account for joint fatigue over time.
m_factorial = lm(
  valgus_torque ~ throw_ID + longitudinal_dev_f * rotational_pos_f,
  data = pulse_data
)

# Tidy summary of model coefficients
m_tidy = m_factorial %>% tidy(conf.int = TRUE, conf.level = 0.95)
print(m_tidy)

# Model goodness of fit (R-squared, sigma / residual SE, F-statistic)
m_fit = m_factorial %>% glance()
print(m_fit)

# 6E. Interaction Plot Framework
# Visualize predicted torque across conditions with 95% confidence intervals
# Parallel lines indicate NO interaction (effects are additive).
# Non-parallel / crossing lines indicate an INTERACTION between displacement and rotation.
pred_grid = pulse_data %>%
  distinct(longitudinal_dev_f, rotational_pos_f) %>%
  mutate(throw_ID = mean(pulse_data$throw_ID)) # hold throw order constant at average

pred_fit = predict(
  m_factorial,
  newdata = pred_grid,
  interval = "confidence",
  level = 0.95
) %>% as_tibble()

plot_interaction = bind_cols(pred_grid, pred_fit) %>%
  ggplot(mapping = aes(
    x = longitudinal_dev_f,
    y = fit,
    group = rotational_pos_f,
    color = rotational_pos_f
  )) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.2) +
  labs(
    title = "Factorial Interaction Plot: Predicted Valgus Torque",
    subtitle = "Controlled for throw sequence order | Error bars = 95% CI",
    x = "Longitudinal Deviation (cm from Ideal 4.5 cm)",
    y = "Predicted Peak Torque (N*m)",
    color = "Rotational\nPosition (deg)"
  )
print(plot_interaction)
