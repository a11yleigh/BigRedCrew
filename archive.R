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

