# Driveline PULSE Sensor Statistics & Data Analysis
# Six Sigma Black Belt Project

# Team Big Red Crew: Alexis Staveski, Ally Schwartz, Collin Shiang, Sophi Ramirez, Vedant Chavan
# SYSEN 5300, Cornell University
# Current Revision Date: 10/27/26

# This code will ... 
  # Import data; perform post-processing & formatting, noting what the ideal location distance and angle is from the medial bony bone on the elbow for each thrower.
  # Note the order in which treatments were applied for block randomization
  # Perform inital data visualization. plot a distribution of torque and distribution of arm speed as a histogram.
  # Create control charts to visualize stability of the process: each thrower is a subgroup. 
      # Make a averages chart, SD chart, and ranges chart. Set up a framework for calculating the performance indicies of cp/cpk/pp/ppk values as applicable and calculating the confidence intervals. 
      # Make a note that if confidence interval is uncertain, bootstrapping is needed.
  # Set up a framework for comparing them to the upper and lower control limits.
  # Perform the factorial design analysis. 
    # Set up a framwork for determing direct effects of longitudinal deviation, rotational deviation, and combined deviation, each controlling for the effects of fatigue from the throw order. 
  # Calculate confidence intervals on the dbar values where necessary.
  # ^ rough descriptions - EDIT LATER"""


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
# Import collected data (CSVs)
# Data information:
# throw_ID is the throw number (1-50 for actual data)
# arm_speed is in revolutions per minute (RPM)
# valgus_torque is in N*m
# longitudinal_pos is the distance from the medial epicondyle bone (from bone towards the hand)
# rotational_pos is the rotational angle from the ideal position (positive = CW looking from the shoulder to the hand)

raw_pulse = read_csv("testData2.csv", col_types = cols())


# Select relevant columns and clean column names
pulse_data = raw_pulse %>%
  select(
    throw_ID,
    thrower_ID,
    arm_speed,
    valgus_torque,
    longitudinal_pos,
    rotational_pos
  ) %>%
  # The raw CSV already contains the corrected speed and torque columns.
  filter(!is.na(throw_ID))

# Former data-entry correction, disabled because the raw CSV is now corrected.
# This swapped arm_speed and valgus_torque for throws 69 and 70.
# pulse_data = pulse_data %>%
#   mutate(
#     temp_speed = arm_speed,
#     arm_speed = if_else(throw_ID %in% c(69, 70), valgus_torque, arm_speed),
#     valgus_torque = if_else(throw_ID %in% c(69, 70), temp_speed, valgus_torque)
#   ) %>%
#   select(-temp_speed)

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
    # Identify each thrower from the recorded ID
    thrower = factor(paste("Thrower", thrower_ID)),

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
# - Some planned throws have missing measurements or placement values.
#   Keep them in pulse_data, but exclude them from analyses requiring those values.
# - Missing full factorial interaction combinations: this dataset tests deviations along one 
#   axis at a time (e.g., no joint condition where longitudinal = 9 cm AND rotational = 90 deg).
# - Post-throw displacement was removed from scope as displacement was observed to be minimal.
glimpse(pulse_data)

# Report missing values rather than silently treating planned throws as observations.
missing_data_summary = pulse_data %>%
  group_by(thrower_ID) %>%
  summarize(across(c(arm_speed, valgus_torque, longitudinal_pos, rotational_pos),
                   ~ sum(is.na(.x))), .groups = "drop")
print(missing_data_summary)
# Other analyses require complete records for every planned row of a thrower.
analysis_data = pulse_data %>%
  group_by(thrower_ID) %>%
  filter(all(complete.cases(pick(arm_speed, valgus_torque,
                                longitudinal_pos, rotational_pos)))) %>%
  ungroup() %>%
  droplevels()


# 3. Initial Data Visualization & Distributional Analysis
# ------------------------------------------------------------------------------
# Compute distribution statistics for peak torque and arm speed using course describe()
torque_desc = analysis_data$valgus_torque %>% describe()
speed_desc  = analysis_data$arm_speed %>% describe()

# 3A. Distribution of Peak Valgus Torque (Histogram + Density Curve)
plot_torque_dist = analysis_data %>%
  ggplot(mapping = aes(x = valgus_torque)) +
  geom_histogram(aes(y = after_stat(density)), bins = 15, fill = "steelblue", color = "white", alpha = 0.7) +
  geom_density(color = "darkblue", size = 1) +
  labs(
    title = "Distribution of Peak Valgus Torque",
    subtitle = "Throwers with complete records only",
    x = "Peak Valgus Torque (N*m)",
    y = "Density",
    caption = torque_desc$caption
  )

# 3B. Distribution of Arm Speed (Histogram + Density Curve)
plot_speed_dist = analysis_data %>%
  ggplot(mapping = aes(x = arm_speed)) +
  geom_histogram(aes(y = after_stat(density)), bins = 15, fill = "darkseagreen", color = "white", alpha = 0.7) +
  geom_density(color = "darkgreen", size = 1) +
  labs(
    title = "Distribution of Arm Speed",
    subtitle = "Throwers with complete records only",
    x = "Arm Speed (rpm)",
    y = "Density",
    caption = speed_desc$caption
  )

# Display distribution plots
print(plot_torque_dist)
print(plot_speed_dist)

# 3C. Fatigue / Time Trend Check: Torque vs. Throw Order
# A scatter plot with a linear best fit to assess whether arm fatigue causes a systematic drift.
plot_fatigue_trend = analysis_data %>%
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
box_panel = analysis_data %>%
  ggplot(mapping = aes(x = factor(set_order), y = valgus_torque)) +
  geom_hline(aes(yintercept = mean(valgus_torque, na.rm = TRUE)), color = "lightgrey", size = 2) +
  geom_jitter(height = 0, width = 0.2, color = "steelblue", alpha = 0.6) +
  geom_boxplot(fill = NA, color = "black", outlier.shape = NA) +
  labs(
    title = "Process Overview: Torque Across 10-Throw Sets",
    x = "Treatment Set (10 Throws Each)",
    y = "Valgus Torque (N*m)",
    caption = torque_desc$caption
  )

hist_panel = analysis_data %>%
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

# 4D. Average Chart: Nine Thrower x Shared 10-Throw Block Subgroups
# Keep a full block only when each of the three throwers has all ten torques.
# All filtering here is local to Sections 4D and 4E.
required_throwers = c(1, 2, 3)
thrower_block_candidates = pulse_data %>%
  filter(thrower_ID %in% required_throwers, !is.na(throw_ID)) %>%
  mutate(block = as.integer(ceiling(throw_ID / 10)))

complete_thrower_blocks = thrower_block_candidates %>%
  group_by(thrower_ID, block) %>%
  summarize(complete = n() == 10 && n_distinct(throw_ID) == 10 &&
              all(is.finite(valgus_torque)), .groups = "drop")
shared_blocks = complete_thrower_blocks %>%
  group_by(block) %>%
  summarize(complete = n_distinct(thrower_ID) == 3 && all(complete),
            .groups = "drop") %>%
  filter(complete)

thrower_block_data = thrower_block_candidates %>%
  semi_join(shared_blocks, by = "block")

if (nrow(shared_blocks) == 0) {
  message("Skipping 4D and 4E: no complete 10-throw blocks shared by all three throwers.")
} else {
  stat_thrower_block = thrower_block_data %>%
    group_by(block, thrower_ID) %>%
    summarize(xbar = mean(valgus_torque), s = sd(valgus_torque),
              nw = n(), .groups = "drop") %>%
    group_by(block) %>%
    mutate(
      # Compare the three thrower subgroups within the same block.
      sigma_s = sqrt(sum((nw - 1) * s^2) / sum(nw - 1)),
      xbbar = mean(xbar),
      ucl_x = xbbar + 3 * sigma_s / sqrt(nw),
      lcl_x = xbbar - 3 * sigma_s / sqrt(nw),
      # Use the same SD standard-error approximation as Section 4B.
      sbar = mean(s),
      ucl_s = sbar + 3 * (sigma_s / sqrt(2 * first(nw - 1))),
      lcl_s = pmax(0, sbar - 3 * (sigma_s / sqrt(2 * first(nw - 1)))),
      thrower_label = paste("Thrower", thrower_ID),
      block_label = paste0("Throws ", (block - 1) * 10 + 1, "–", block * 10)
    ) %>%
    ungroup() %>%
    arrange(block, thrower_ID) %>%
    mutate(subgroup_index = row_number(),
           subgroup_label = paste0(thrower_label, "\n", block_label))

  # Each row compares throwers within one shared throw block.
  thrower_block_limits = stat_thrower_block %>%
    group_by(block) %>%
    summarize(left = min(thrower_ID) - 0.4,
              right = max(thrower_ID) + 0.4,
              across(c(xbbar, ucl_x, lcl_x, sbar, ucl_s, lcl_s), first),
              .groups = "drop") %>%
    mutate(block_label = paste0("Throws ", (block - 1) * 10 + 1, "–", block * 10))

  xbar_thrower = ggplot(stat_thrower_block, aes(x = thrower_ID, y = xbar)) +
    geom_segment(data = thrower_block_limits,
                 aes(x = left, xend = right, y = xbbar, yend = xbbar),
                 inherit.aes = FALSE, color = "grey50") +
    geom_segment(data = thrower_block_limits,
                 aes(x = left, xend = right, y = ucl_x, yend = ucl_x),
                 inherit.aes = FALSE, color = "firebrick", linetype = "dashed") +
    geom_segment(data = thrower_block_limits,
                 aes(x = left, xend = right, y = lcl_x, yend = lcl_x),
                 inherit.aes = FALSE, color = "firebrick", linetype = "dashed") +
    geom_line(aes(group = block), color = "steelblue", linewidth = 0.8) +
    geom_point(color = "steelblue", size = 3) +
    geom_text(data = thrower_block_limits,
              aes(x = right, y = ucl_x, label = paste0("UCL = ", round(ucl_x, 2))),
              inherit.aes = FALSE, hjust = 1, vjust = -0.4, color = "firebrick", size = 3) +
    geom_text(data = thrower_block_limits,
              aes(x = right, y = lcl_x, label = paste0("LCL = ", round(lcl_x, 2))),
              inherit.aes = FALSE, hjust = 1, vjust = 1.3, color = "firebrick", size = 3) +
    scale_y_continuous(expand = expansion(mult = c(0.18, 0.18))) +
    geom_text(aes(label = paste0("Mean = ", round(xbar, 2))), vjust = -1, size = 3.5) +
    facet_grid(block_label ~ ., scales = "free_y") +
    scale_x_continuous(breaks = required_throwers,
                       labels = paste("Thrower", required_throwers)) +
    labs(
      title = "Average Valgus Torque by Thrower and 10-Throw Block",
      subtitle = "Complete blocks shared by all three throwers | 10 throws per point",
      x = "Thrower subgroup", y = "Mean Valgus Torque (N*m)",
      caption = "Y-axis scales and control limits differ by row.\nGrey line: block grand mean across throwers. Dashed lines: +/- 3 SE using pooled within-thrower SD."
    )
  print(xbar_thrower)
}

# 4E. SD Chart: Same Nine Thrower x Shared 10-Throw Block Subgroups
if (nrow(shared_blocks) > 0) {
  s_thrower = ggplot(stat_thrower_block, aes(x = thrower_ID, y = s)) +
    geom_segment(data = thrower_block_limits,
                 aes(x = left, xend = right, y = sbar, yend = sbar),
                 inherit.aes = FALSE, color = "grey50") +
    geom_segment(data = thrower_block_limits,
                 aes(x = left, xend = right, y = ucl_s, yend = ucl_s),
                 inherit.aes = FALSE, color = "firebrick", linetype = "dashed") +
    geom_segment(data = thrower_block_limits,
                 aes(x = left, xend = right, y = lcl_s, yend = lcl_s),
                 inherit.aes = FALSE, color = "firebrick", linetype = "dashed") +
    geom_line(aes(group = block), color = "forestgreen", linewidth = 0.8) +
    geom_point(color = "forestgreen", size = 3) +
    geom_text(data = thrower_block_limits,
              aes(x = right, y = ucl_s, label = paste0("UCL = ", round(ucl_s, 2))),
              inherit.aes = FALSE, hjust = 1, vjust = -0.4, color = "firebrick", size = 3) +
    geom_text(data = thrower_block_limits,
              aes(x = right, y = lcl_s, label = paste0("LCL = ", round(lcl_s, 2))),
              inherit.aes = FALSE, hjust = 1, vjust = 1.3, color = "firebrick", size = 3) +
    scale_y_continuous(expand = expansion(mult = c(0.18, 0.18))) +
    geom_text(aes(label = paste0("SD = ", round(s, 2))), vjust = -1, size = 3.5) +
    facet_grid(block_label ~ ., scales = "free_y") +
    scale_x_continuous(breaks = required_throwers,
                       labels = paste("Thrower", required_throwers)) +
    labs(
      title = "Valgus Torque SD by Thrower and 10-Throw Block",
      subtitle = "Complete blocks shared by all three throwers | 10 throws per point",
      x = "Thrower subgroup", y = "Within-Subgroup SD (N*m)",
      caption = "Y-axis scales and control limits differ by row.\nGrey line: block average SD across throwers. Dashed lines: average SD +/- 3 approximate SE, using the same formula as Section 4B."
    )
  print(s_thrower)

}


# 4F. All Nine Subgroups Together: Averages and SD
# Calculate one set of limits across all nine subgroups, using Section 4 formulas.
if (nrow(shared_blocks) > 0) {
  limits_all_subgroups = stat_thrower_block %>%
    summarize(
      left = min(subgroup_index) - 0.4,
      right = max(subgroup_index) + 0.4,
      sigma_s = sqrt(sum((nw - 1) * s^2) / sum(nw - 1)),
      xbbar = mean(xbar),
      ucl_x = xbbar + 3 * sigma_s / sqrt(first(nw)),
      lcl_x = xbbar - 3 * sigma_s / sqrt(first(nw)),
      sbar = mean(s),
      ucl_s = sbar + 3 * sigma_s / sqrt(2 * (first(nw) - 1)),
      lcl_s = pmax(0, sbar - 3 * sigma_s / sqrt(2 * (first(nw) - 1)))
    )

  # Build both views with the same subgroup order and labels.
  charts_all_subgroups = list()
  for (metric in c("xbar", "s")) {
    is_average = metric == "xbar"
    center = if (is_average) "xbbar" else "sbar"
    upper = if (is_average) "ucl_x" else "ucl_s"
    lower = if (is_average) "lcl_x" else "lcl_s"
    chart = ggplot(stat_thrower_block,
                   aes(x = subgroup_index, y = .data[[metric]])) +
      geom_segment(data = limits_all_subgroups,
                   aes(x = left, xend = right, y = .data[[center]], yend = .data[[center]]),
                   inherit.aes = FALSE, color = "grey50") +
      geom_segment(data = limits_all_subgroups,
                   aes(x = left, xend = right, y = .data[[upper]], yend = .data[[upper]]),
                   inherit.aes = FALSE, color = "firebrick", linetype = "dashed") +
      geom_segment(data = limits_all_subgroups,
                   aes(x = left, xend = right, y = .data[[lower]], yend = .data[[lower]]),
                   inherit.aes = FALSE, color = "firebrick", linetype = "dashed") +
      geom_line(aes(group = 1), color = "steelblue", linewidth = 0.8) +
      geom_point(aes(color = block_label), size = 3) +
      geom_text(aes(label = round(.data[[metric]], 2)), vjust = -1, size = 3) +
      scale_x_continuous(breaks = stat_thrower_block$subgroup_index,
                         labels = stat_thrower_block$subgroup_label) +
      scale_y_continuous(expand = expansion(mult = c(0.15, 0.15))) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) +
      labs(
        title = if (is_average) "Average Valgus Torque: All Nine Subgroups" else
          "Valgus Torque SD: All Nine Subgroups",
        subtitle = "10 throws per subgroup | Only complete blocks shared by all three throwers",
        x = "Thrower and throw block", color = "Throw block",
        y = if (is_average) "Mean Valgus Torque (N*m)" else "Within-Subgroup SD (N*m)",
        caption = "Grey: overall center. Dashed red: common UCL/LCL calculated across all nine subgroups.\nPoints are ordered by block, then thrower; connections do not represent a time series."
      )
    charts_all_subgroups[[metric]] = chart
  }
  xbar_all_subgroups = charts_all_subgroups$xbar
  s_all_subgroups = charts_all_subgroups$s
  # Stack the averages and SD charts with matching subgroup order.
  combined_all_subgroups = ggarrange(
    xbar_all_subgroups, s_all_subgroups,
    ncol = 1, nrow = 2, align = "v",
    common.legend = TRUE, legend = "bottom"
  )
  print(combined_all_subgroups)
}


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
v_short     = sum(stat_subgroups$df)
v_total     = nrow(analysis_data) - 1

se_cp  = est_cp * sqrt(1 / (2 * v_short))
se_pp  = est_pp * sqrt(1 / (2 * v_total))
se_cpk = est_cpk * sqrt(1 / (2 * v_short) + 1 / (9 * nrow(analysis_data) * (est_cpk^2)))
se_ppk = est_ppk * sqrt(1 / (2 * v_total) + 1 / (9 * nrow(analysis_data) * (est_ppk^2)))

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
# Placement summaries use only throwers with complete records.
model_data = analysis_data
groups = model_data %>%
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
    mean(groups$xbar[groups$longitudinal_dev_f == "4.5" & groups$rotational_pos_f == "0"]) - baseline_xbar,
    mean(groups$xbar[groups$longitudinal_dev_f == "9" & groups$rotational_pos_f == "0"])   - baseline_xbar,
    mean(groups$xbar[groups$longitudinal_dev_f == "-9" & groups$rotational_pos_f == "0"])  - baseline_xbar,
    mean(groups$xbar[groups$longitudinal_dev_f == "-13.5" & groups$rotational_pos_f == "0"]) - baseline_xbar,
    mean(groups$xbar[groups$rotational_pos_f == "90" & groups$longitudinal_dev_f == "0"])   - baseline_xbar,
    mean(groups$xbar[groups$rotational_pos_f == "180" & groups$longitudinal_dev_f == "0"])  - baseline_xbar
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
# Include thrower as a block when multiple throwers have complete data.
model_formula = valgus_torque ~ throw_ID + longitudinal_dev_f * rotational_pos_f
if (nlevels(model_data$thrower) > 1) {
  model_formula = update(model_formula, . ~ . + thrower)
}
interaction_formula = update(model_formula, . ~ . + longitudinal_dev_f:rotational_pos_f)
interaction_matrix = model.matrix(interaction_formula, data = model_data)
interaction_estimable = qr(interaction_matrix)$rank == ncol(interaction_matrix)
if (!interaction_estimable) {
  message("Skipping interaction model and predictions: placement combinations are missing.")
} else {
m_factorial = lm(interaction_formula, data = model_data, na.action = na.fail)

# Tidy summary of model coefficients
m_tidy = m_factorial %>% tidy(conf.int = TRUE, conf.level = 0.95)
print(m_tidy)

# Model goodness of fit (R-squared, sigma / residual SE, F-statistic)
m_fit = m_factorial %>% glance()
print(m_fit)

# 6E. Interaction Plot Framework
# Visualize predicted torque across conditions with 95% confidence intervals
# Run only when the interaction model is estimable from the observed design.
pred_grid = model_data %>%
  distinct(longitudinal_dev_f, rotational_pos_f) %>%
  mutate(thrower = factor(levels(model_data$thrower)[1], levels = levels(model_data$thrower)),
         throw_ID = mean(model_data$throw_ID)) # hold throw order constant at average

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
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.2) +
  labs(
    title = "Placement Effects: Predicted Valgus Torque",
    subtitle = paste("At mean throw order for", levels(model_data$thrower)[1],
                     "| Error bars = 95% CI"),
    x = "Longitudinal Deviation (cm from Ideal 4.5 cm)",
    y = "Predicted Peak Torque (N*m)",
    color = "Rotational\nPosition (deg)"
  )
print(plot_interaction)
}
