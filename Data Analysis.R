# Driveline PULSE Sensor Statistics & Data Analysis
# Six Sigma Black Belt Project

# Team Big Red Crew: Alexis Staveski, Ally Schwartz, Collin Shiang, Sophi Ramirez, Vedant Chavan
# SYSEN 5300, Cornell University
# Current Revision Date: 10/25/26

# This code will ... 
# - importing data and doing any post-processing + formatting, noting what the ideal location distance and angle is from the medial bony bone on the elbow for each thrower. also note the order in which treatments were applied (relevant for block randomization)
# - doing inital data visualization. plot a distribution of torque and distribution of arm speed as a histogram.
# - make control charts to visualize stability of the process. each thrower is a subgroup. make a averages chart, deviation chart, and ranges chart just to see. set up a framework for calculating the cp/cpk/pp/ppk values as applicable and calculating the confidence intervals. make a note that if confidence interval is uncertain, bootstrapping is needed. set up a framework for comparing them to the upper and lower control limits that we decide on.
# - perform the factorial design analysis. set up a framwork for determing direct effects of longitudinal deviation, rotational deviation, and combined deviation, each controlling for the effects of fatigue from the throw order. calculate confidence intervals on the dbar values where necessary.
# ^ rough descriptions - EDIT LATER


# 1. Setup & Package Loading
# ------------------------------------------------------------------------------
# We load the tidyverse and statistical analysis packages used throughout SYSEN 5300.
library(dplyr)    # For data manipulation and piping (%>%)
library(readr)    # For reading CSV files
library(tibble)   # For tidy data frames
library(ggplot2)  # For visualization conforming to the course theme
library(broom)    # For tidying regression and test outputs (tidy, glance)
library(moments)  # For skewness and kurtosis
library(ggpubr)   # For arranging multi-panel plots (ggarrange)

# Course SPC functions can be sourced directly if available:
if (file.exists("sysen/functions/functions_process_control.R")) {
  source("sysen/functions/functions_process_control.R")
}

# Universal course plot theme (theme_classic with centered titles & clean borders)
theme_set(
  theme_classic(base_size = 14) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      plot.caption = element_text(hjust = 0.5),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      panel.border = element_rect(fill = NA, color = "grey"),
      plot.margin = margin(r = 10, t = 10, b = 10, l = 10)
    )
)

# In-house course descriptive summary function
describe = function(x) {
  tibble(x) %>%
    summarize(
      mean     = mean(x, na.rm = TRUE),
      sd       = sd(x, na.rm = TRUE),
      skew     = moments::skewness(x, na.rm = TRUE),
      kurtosis = moments::kurtosis(x, na.rm = TRUE)
    ) %>%
    mutate(
      caption = paste(
        "Mean: ", round(mean, 2), " | ",
        "SD: ", round(sd, 2), " | ",
        "Skewness: ", round(skew, 2), " | ",
        "Kurtosis: ", round(kurtosis, 2),
        sep = ""
      )
    )
}

# In-house capability and performance index functions (from Week 6)
cp  = function(sigma_s, upper, lower) { abs(upper - lower) / (6 * sigma_s) }
pp  = function(sigma_t, upper, lower) { abs(upper - lower) / (6 * sigma_t) }
cpk = function(mu, sigma_s, lower = NULL, upper = NULL) {
  a = if (!is.null(lower)) abs(mu - lower) / (3 * sigma_s) else NULL
  b = if (!is.null(upper)) abs(upper - mu) / (3 * sigma_s) else NULL
  if (!is.null(lower) & !is.null(upper)) min(a, b) else if (is.null(lower)) b else a
}
ppk = function(mu, sigma_t, lower = NULL, upper = NULL) {
  a = if (!is.null(lower)) abs(mu - lower) / (3 * sigma_t) else NULL
  b = if (!is.null(upper)) abs(upper - mu) / (3 * sigma_t) else NULL
  if (!is.null(lower) & !is.null(upper)) min(a, b) else if (is.null(lower)) b else a
}


# 2. Data Ingestion, Cleaning & Preprocessing
# ------------------------------------------------------------------------------
# Load testData2.csv (path adapts to project or workspace working directory)
data_path = if (file.exists("testData2.csv")) "testData2.csv" else "Project/testData2.csv"
raw_pulse = read_csv(data_path, col_types = cols())

# Select relevant columns and clean column names
pulse_data = raw_pulse %>%
  select(
    throw_ID,
    arm_speed,
    valgus_torque,
    longitudinal_pos,
    rotational_pos
  ) %>%
  filter(!is.na(throw_ID)) %>%
  # Data Hygiene Correction:
  # In throws 69 and 70, arm_speed and valgus_torque were transposed during data entry
  mutate(
    temp_speed = arm_speed,
    arm_speed = if_else(throw_ID %in% c(69, 70), valgus_torque, arm_speed),
    valgus_torque = if_else(throw_ID %in% c(69, 70), temp_speed, valgus_torque)
  ) %>%
  select(-temp_speed)

# Anatomical Reference Notes:
# - Ideal sensor placement: 4.5 cm distal to medial epicondyle (medial bone), 0 deg rotation.
# - Longitudinal Deviation = longitudinal_pos - 4.5 cm (0 cm = Ideal).
# - Rotational Deviation   = rotational_pos (0 deg = Ideal).
#
# Experimental Design & Block Randomization Notes:
# - In this preliminary single-thrower test, treatments were conducted in sequential sets of 10 throws.
# - For the multi-thrower study, each Thrower will act as a BLOCK, with the 5 treatment sets 
#   randomized in order (5! = 120 possible permutations) to prevent fatigue confounding.
pulse_data = pulse_data %>%
  mutate(
    # Identify the thrower (currently single thrower for test data)
    thrower = factor("Thrower 1"),

    # Track the chronological treatment set (1 to 7)
    set_order = as.integer(ceiling(throw_ID / 10)),
    
    # Calculate longitudinal deviation from ideal (4.5 cm)
    longitudinal_dev = longitudinal_pos - 4.5,
    
    # Convert deviations to factors with the IDEAL / baseline condition as the FIRST level.
    # Following course conventions, the reference group must be first so differences (dbar)
    # and linear model contrasts directly measure changes against ideal baseline.
    longitudinal_dev_f = factor(
      as.character(longitudinal_dev),
      levels = c("0", "4.5", "9", "-9", "-13.5")
    ),
    rotational_pos_f = factor(
      as.character(rotational_pos),
      levels = c("0", "90", "180")
    ),
    
    # Combined condition label for subgrouping
    condition_group = paste(longitudinal_dev, "cm_", rotational_pos, "deg", sep = "")
  )

# Data Completeness Check / Missing Data Notes:
# - Missing multi-thrower variation: currently n = 1 thrower (70 throws).
# - Missing full factorial interaction combinations: this dataset tests deviations along one 
#   axis at a time (e.g., no joint condition where longitudinal = 9 cm AND rotational = 90 deg).
# - Post-throw displacement was removed from scope as displacement was observed to be minimal.
glimpse(pulse_data)


# 3. Initial Data Visualization & Distributional Analysis
# ------------------------------------------------------------------------------
# Compute distribution statistics for peak torque and arm speed using course describe()
torque_desc = pulse_data$valgus_torque %>% describe()
speed_desc  = pulse_data$arm_speed %>% describe()

# 3A. Distribution of Peak Valgus Torque (Histogram + Density Curve)
plot_torque_dist = pulse_data %>%
  ggplot(mapping = aes(x = valgus_torque)) +
  geom_histogram(aes(y = after_stat(density)), bins = 15, fill = "steelblue", color = "white", alpha = 0.7) +
  geom_density(color = "darkblue", size = 1) +
  labs(
    title = "Distribution of Peak Valgus Torque",
    subtitle = "Observed across all sensor placement trials",
    x = "Peak Valgus Torque (N*m)",
    y = "Density",
    caption = torque_desc$caption
  )

# 3B. Distribution of Arm Speed (Histogram + Density Curve)
plot_speed_dist = pulse_data %>%
  ggplot(mapping = aes(x = arm_speed)) +
  geom_histogram(aes(y = after_stat(density)), bins = 15, fill = "darkseagreen", color = "white", alpha = 0.7) +
  geom_density(color = "darkgreen", size = 1) +
  labs(
    title = "Distribution of Arm Speed",
    subtitle = "Observed across all sensor placement trials",
    x = "Arm Speed (rpm)",
    y = "Density",
    caption = speed_desc$caption
  )

# Display distribution plots
print(plot_torque_dist)
print(plot_speed_dist)

# 3C. Fatigue / Time Trend Check: Torque vs. Throw Order
# A scatter plot with a linear best fit to assess whether arm fatigue causes a systematic drift.
plot_fatigue_trend = pulse_data %>%
  ggplot(mapping = aes(x = throw_ID, y = valgus_torque)) +
  geom_point(color = "grey30", size = 2.5, alpha = 0.8) +
  geom_smooth(method = "lm", color = "firebrick", se = TRUE, fill = "lightpink", alpha = 0.3) +
  labs(
    title = "Fatigue Trend: Valgus Torque vs. Throw Sequence",
    subtitle = "Testing for downward drift across sequential throws",
    x = "Throw Order (Sequential Throw ID)",
    y = "Peak Valgus Torque (N*m)"
  )
print(plot_fatigue_trend)

# 3D. Process Overview Diagram (Jittered Boxplot + Marginal Flipped Histogram)
# Following the course ggprocess() design from Recitation 5 & SPC Learning Check
box_panel = pulse_data %>%
  ggplot(mapping = aes(x = factor(set_order), y = valgus_torque)) +
  geom_hline(aes(yintercept = mean(valgus_torque)), color = "lightgrey", size = 2) +
  geom_jitter(height = 0, width = 0.2, color = "steelblue", alpha = 0.6) +
  geom_boxplot(fill = NA, color = "black", outlier.shape = NA) +
  labs(
    title = "Process Overview: Torque Across 10-Throw Sets",
    x = "Treatment Set (10 Throws Each)",
    y = "Valgus Torque (N*m)",
    caption = torque_desc$caption
  )

hist_panel = pulse_data %>%
  ggplot(mapping = aes(x = valgus_torque)) +
  geom_histogram(bins = 15, color = "white", fill = "grey70") +
  theme_void() +
  coord_flip()

process_overview = ggarrange(box_panel, hist_panel, widths = c(4, 1), align = "h")
print(process_overview)


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

# Calculate point estimates for the four indices
est_cp  = cp(sigma_s = sigma_s_val, upper = spec_upper, lower = spec_lower)
est_pp  = pp(sigma_t = sigma_t_val, upper = spec_upper, lower = spec_lower)
est_cpk = cpk(mu = xbbar_val, sigma_s = sigma_s_val, lower = spec_lower, upper = spec_upper)
est_ppk = ppk(mu = xbbar_val, sigma_t = sigma_t_val, lower = spec_lower, upper = spec_upper)

# Theoretical 95% Confidence Intervals for Capability Indices
# Degrees of freedom:
# v_short = k * (nw - 1)
# v_total = total_n - 1
k_subgroups = nrow(stat_subgroups)
nw_size     = first(stat_subgroups$nw)
v_short     = k_subgroups * (nw_size - 1)
v_total     = nrow(pulse_data) - 1

se_cp  = est_cp * sqrt(1 / (2 * v_short))
se_pp  = est_pp * sqrt(1 / (2 * v_total))
se_cpk = est_cpk * sqrt(1 / (2 * v_short) + 1 / (9 * nrow(pulse_data) * (est_cpk^2)))
se_ppk = est_ppk * sqrt(1 / (2 * v_total) + 1 / (9 * nrow(pulse_data) * (est_ppk^2)))

capability_summary = tibble(
  index    = c("Cp", "Cpk", "Pp", "Ppk"),
  sigma    = c("sigma_short", "sigma_short", "sigma_total", "sigma_total"),
  estimate = c(est_cp, est_cpk, est_pp, est_ppk),
  se       = c(se_cp, se_cpk, se_pp, se_ppk),
  lower_95 = estimate - qnorm(0.975) * se,
  upper_95 = estimate + qnorm(0.975) * se
)
print(capability_summary)

# NOTE ON UNCERTAINTY & BOOTSTRAPPING:
# If the theoretical confidence intervals are wide, straddle 1.0, or if the underlying
# normality assumption is violated (e.g., severe skewness in rotated placements),
# BOOTSTRAPPING with replacement (sample_n(size = n(), replace = TRUE)) is required
# to obtain empirical percentiles (2.5% and 97.5%) for the capability indices.
#
# Skeleton for Bootstrapping Capability Indices:
# set.seed(5300)
# boot_indices = tibble(rep = 1:1000) %>%
#   group_by(rep) %>%
#   reframe(pulse_data) %>%
#   group_by(rep) %>%
#   sample_n(size = n(), replace = TRUE) %>%
#   summarize(
#     boot_mu = mean(valgus_torque),
#     boot_sigma_t = sd(valgus_torque),
#     boot_ppk = ppk(boot_mu, boot_sigma_t, spec_lower, spec_upper)
#   )


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
direct_effects = dbbar_direct %>%
  mutate(
    se       = error_pooled$se,
    lower_95 = dbar - qnorm(0.975) * se,
    upper_95 = dbar + qnorm(0.975) * se,
    label    = paste(round(dbar, 2), " [", round(lower_95, 2), ", ", round(upper_95, 2), "]", sep = "")
  )
print(direct_effects)

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
