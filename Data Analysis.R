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


# ------------------------------------------------------------------------------

# Import Experimental Data - data expected in CSV format
raw_data = read_csv("throwData_10.30.csv", col_types = cols())
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

# compute mean, SD, skewness, kurtosis with helper function
torque_desc = pulse_data$valgus_torque %>% describe()
speed_desc  = pulse_data$arm_speed %>% describe()

# create plots:

# apply same ggplot theme across all plots
set_theme()

# Histogram of peak valgus torque for each thrower
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
    title = "Distribution of Peak Valgus Torque by Thrower",
    x = "Peak Valgus Torque (N*m)",
    y = "Count"
  )

# Histogram of arm speed for each thrower
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
    title = "Distribution of Arm Speed by Thrower",
    x = "Arm Speed (rpm)",
    y = "Count"
  )

# Boxplot with overlaid jitter of valgus torque across sensor placements by thrower
plot_torque_box_jitter = pulse_data %>%
  ggplot() +
  geom_boxplot(
    mapping = aes(x = sensor_placement, y = valgus_torque, fill = thrower),
    position = position_dodge(width = 0.8),
    outlier.shape = NA,
    alpha = 0.35,
    width = 0.7
  ) +
  geom_point(
    mapping = aes(x = sensor_placement, y = valgus_torque, color = thrower),
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
    size = 2,
    alpha = 0.8
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Peak Valgus Torque by Sensor Placement across Throwers",
    x = "Sensor Placement",
    y = "Peak Valgus Torque (N*m)",
    color = "Thrower",
    fill  = "Thrower"
  )

# Boxplot with overlaid jitter of arm speed across sensor placements by thrower
plot_speed_box_jitter = pulse_data %>%
  ggplot() +
  geom_boxplot(
    mapping = aes(x = sensor_placement, y = arm_speed, fill = thrower),
    position = position_dodge(width = 0.8),
    outlier.shape = NA,
    alpha = 0.35,
    width = 0.7
  ) +
  geom_point(
    mapping = aes(x = sensor_placement, y = arm_speed, color = thrower),
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
    size = 2,
    alpha = 0.8
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Arm Speed by Sensor Placement across Throwers",
    x = "Sensor Placement",
    y = "Arm Speed (rpm)",
    color = "Thrower",
    fill  = "Thrower"
  )

# Display plots
print(plot_torque_hist)
print(plot_speed_hist)
print(plot_torque_box_jitter)
print(plot_speed_box_jitter)



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
