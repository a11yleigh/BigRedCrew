# Driveline PULSE Sensor Statistics & Data Analysis - Helper Functions
# Six Sigma Black Belt Project

# Team Big Red Crew: Alexis Staveski, Ally Schwartz, Collin Shiang, Sophi Ramirez, Vedant Chavan
# SYSEN 5300, Cornell University

# Citation:
# Many functions in this script originate from and build upon the Cornell SYSEN 5300
# course repository (https://github.com/timothyfraser/sysen), (c) 2022 Timothy Fraser, PhD

# ------------------------------------------------------------------------------
# Available functions to use:
# set_theme(): Sets universal course ggplot theme with centered titles and clean borders.
# describe(x): Computes sample mean, SD, skewness, kurtosis, and returns a formatted caption string.
# ggprocess(x, y, xlab, ylab): Creates process overview plot with a jittered boxplot and flipped marginal histogram.
# get_stat_s(x, y): Calculates within-subgroup statistics (xbar, range, SD, nw, df, sigma_s, se, control limits).
# get_stat_t(x, y): Calculates overall process statistics (xbbar, rbar, sbar, sigma_s, sigma_t, total sample size).
# get_labels(data): Extracts plot coordinates and text labels for chart centerlines and +/- 3-sigma control limits.
# dn(n, reps): Simulates normal ranges to compute range control constants (d2, d3, D3, D4).
# bn(n, reps): Simulates normal standard deviations to compute S-chart control constants (b2, b3, C4, A3, B3, B4).
# limits_avg(x, y): Computes upper and lower control limits for an Averages (Xbar) chart using the A3 constant.
# limits_s(x, y): Computes upper and lower control limits for a Standard Deviation (S) chart using B3 and B4.
# limits_r(x, y): Computes upper and lower control limits for a Range (R) chart using D3 and D4.
# limits_mr(x, y): Computes upper and lower control limits for a Moving Range (MR) chart (subgroup size n = 1).
# ggxbar(x, y, xlab, ylab): Generates an Average (Xbar) control chart with 3-sigma limits and labels.
# ggs(x, y, xlab, ylab): Generates a Standard Deviation (S) control chart using B3/B4 control limits.
# ggr(x, y, xlab, ylab): Generates a Range (R) control chart using D3/D4 control limits.
# ggmr(x, y, xlab, ylab): Generates a Moving Range (MR) control chart for individual observations.
# ggp(t, x, n, xlab, ylab): Generates a Fraction Defective (p) control chart.
# ggnp(t, x, n, xlab, ylab): Generates a Number of Defectives (np) control chart.
# ggu(t, x, xlab, ylab): Generates a Defects per Unit (u) control chart.
# cp(sigma_s, upper, lower): Computes potential capability index Cp using short-term within-subgroup variation.
# pp(sigma_t, upper, lower): Computes potential performance index Pp using total overall variation.
# cpk(mu, sigma_s, lower, upper): Computes actual capability index Cpk accounting for process centering.
# ppk(mu, sigma_t, lower, upper): Computes actual performance index Ppk accounting for process centering.
# get_index(x, y, index, upper, lower, reps, ci_level, by_subgroup): Resamples data to bootstrap empirical 95% CIs for capability indices.
# calc_capability_summary(mu, sigma_s, sigma_t, lower, upper, ...): Computes Cp, Cpk, Pp, Ppk with analytical SEs and 95% CIs into a tidy table.
# format_ci_label(estimate, lower, upper, digits): Formats point estimate and bounds into "estimate [lower, upper]".
# se_factorial(formula, data): Computes pooled standard error across treatment combinations for a factorial model.
# dbar_oneway(formula, data): Computes one-way main treatment effect (dbar) between high and low factor levels.
# dbar_twoway(formula, data): Computes two-way interaction effect (dbar) between two factors.
# dbar_threeway(formula, data): Computes three-way interaction effect (dbar) across three factors.
# calc_direct_effects_table(dbar_df, se_pooled, conf_level): Formats a tidy direct effects table with pooled SE and 95% CIs.
# ------------------------------------------------------------------------------


library(dplyr)
library(readr)
library(tibble)
library(ggplot2)
library(ggpubr)
library(moments)


# ==============================================================================
# Universal Plot Theme
# ==============================================================================

# Sets the universal course theme for all ggplots:
# theme_classic with centered title/subtitle/caption, clean borders, bottom legend.
set_theme = function() {
  theme_set(
    theme_bw(base_size = 14) +
      theme(
        legend.position = "right",
        plot.title      = element_text(hjust = 0.5),
        plot.subtitle   = element_text(hjust = 0.5),
        plot.caption    = element_text(hjust = 0.5),
        panel.border    = element_rect(fill = NA, color = "grey"),
        plot.margin     = margin(r = 0)
      )
  )
}


# ==============================================================================
# General Descriptive Statistics
# ==============================================================================

# Computes sample mean, standard deviation, skewness, and kurtosis for a numeric vector x.
# Returns a single-row tibble with a pre-formatted caption string.
describe = function(x) {
  tibble(x) %>%
    summarize(
      mean     = mean(x, na.rm = TRUE),
      sd       = sd(x, na.rm = TRUE),
      skew     = moments::skewness(x, na.rm = TRUE),
      kurtosis = moments::kurtosis(x, na.rm = TRUE)
    ) %>%
    mutate(
      caption = paste0(
        "Mean: ", round(mean, 2), " | ",
        "SD: ", round(sd, 2), " | ",
        "Skewness: ", round(skew, 2), " | ",
        "Kurtosis: ", round(kurtosis, 2)
      )
    )
}

# Creates a process overview plot combining a jittered boxplot over subgroups
# with a marginal flipped histogram on the right (via ggarrange).
ggprocess = function(x, y, xlab = "Subgroup", ylab = "Metric") {
  data = tibble(x = x, y = y)
  stat = data %>% summarize(mu = mean(y))
  tab  = describe(data$y)

  g1 = ggplot() +
    geom_hline(data = stat, mapping = aes(yintercept = mu), color = "lightgrey", linewidth = 3) +
    geom_jitter(data = data, mapping = aes(x = x, y = y), height = 0, width = 0.25) +
    geom_boxplot(data = data, mapping = aes(x = x, y = y, group = x)) +
    labs(x = xlab, y = ylab, subtitle = "Process Overview", caption = tab$caption)

  g2 = ggplot() +
    geom_histogram(data = data, mapping = aes(x = y), bins = 15, color = "white", fill = "grey") +
    theme_void() +
    coord_flip()

  output = ggpubr::ggarrange(plotlist = list(g1, g2), nrow = 1, ncol = 2, widths = c(4, 1), align = "h")
  return(output)
}


# ==============================================================================
# Statistical Process Control (SPC) Helper & Limit Functions
# ==============================================================================

# Computes within-subgroup statistics (xbar, range, SD, nw, df) along with
# grand mean (xbbar), pooled short-term variation (sigma_s), total variation (sigma_t),
# standard error (se), and standard 3-sigma control limits.
get_stat_s = function(x, y) {
  data    = tibble(x = x, y = y)
  sigma_t = sd(y)

  output = data %>%
    group_by(x) %>%
    summarize(
      xbar = mean(y),
      r    = max(y) - min(y),
      s    = sd(y),
      nw   = n(),
      df   = nw - 1,
      .groups = "drop"
    ) %>%
    mutate(
      xbbar   = mean(xbar),
      sigma_s = sqrt(mean(s^2)),
      sigma_t = sigma_t,
      se      = sigma_s / sqrt(nw),
      upper   = mean(xbar) + 3 * se,
      lower   = mean(xbar) - 3 * se
    )

  return(output)
}

# Summarizes overall (total) process statistics across all subgroups:
# grand mean (xbbar), average range (rbar), average SD (sbar), sigma_s, sigma_t, and total n.
get_stat_t = function(x, y) {
  data    = tibble(x = x, y = y)
  sigma_t = sd(y)

  stat_s = data %>%
    group_by(x) %>%
    summarize(
      xbar = mean(y),
      r    = max(y) - min(y),
      s    = sd(y),
      nw   = n(),
      df   = nw - 1,
      .groups = "drop"
    )

  output = stat_s %>%
    summarize(
      xbbar   = mean(xbar),
      rbar    = mean(r),
      sbar    = mean(s),
      sigma_s = sqrt(sum(s^2) / sum(nw)),
      sigma_t = sigma_t,
      n       = sum(nw)
    )

  return(output)
}

# Extracts label coordinates and formatted text for chart centerlines and +/- 3s limits.
get_labels = function(data) {
  data %>%
    reframe(
      x     = c(max(x), max(x), max(x)),
      type  = c("xbbar", "upper", "lower"),
      name  = c("mean", "+3 s", "-3 s"),
      value = c(mean(xbar), max(upper), min(lower))
    ) %>%
    mutate(
      value = round(value, 2),
      text  = paste0(name, " = ", value)
    )
}

# Simulates normal ranges to compute range control constants: d2, d3, D3, D4.
dn = function(n = 12, reps = 1e4) {
  tibble(rep = 1:reps) %>%
    group_by(rep) %>%
    summarize(r = rnorm(n = n, mean = 0, sd = 1) %>% range() %>% diff() %>% abs(), .groups = "drop") %>%
    summarize(
      d2 = mean(r),
      d3 = sd(r),
      D3 = 1 - 3 * (d3 / d2),
      D4 = 1 + 3 * (d3 / d2),
      D3 = if_else(D3 < 0, true = 0, false = D3)
    )
}

# Simulates normal standard deviations to compute S-chart control constants: b2, b3, C4, A3, B3, B4.
bn = function(n, reps = 1e4) {
  tibble(rep = 1:reps) %>%
    group_by(rep) %>%
    summarize(s = rnorm(n, mean = 0, sd = 1) %>% sd(), .groups = "drop") %>%
    summarize(
      b2 = mean(s),
      b3 = sd(s),
      C4 = b2,
      A3 = 3 / (b2 * sqrt(n)),
      B3 = 1 - 3 * b3 / b2,
      B4 = 1 + 3 * b3 / b2,
      B3 = if_else(B3 < 0, true = 0, false = B3)
    )
}

# Computes control limits for an Averages (Xbar) chart using the A3 control constant.
limits_avg = function(x, y) {
  data = tibble(x = x, y = y)

  stat_s = data %>%
    group_by(x) %>%
    summarize(xbar = mean(y), s = sd(y), nw = n(), df = nw - 1, .groups = "drop")

  constants = stat_s %>%
    dplyr::select(nw) %>%
    distinct() %>%
    group_by(nw) %>%
    summarize(A3 = bn(n = nw, reps = 1e4)$A3, .groups = "drop")

  stat_s = stat_s %>%
    left_join(by = "nw", y = constants) %>%
    mutate(
      sbar  = sqrt(sum(df * s^2) / sum(df)),
      xbbar = mean(xbar),
      lower = xbbar - A3 * sbar,
      upper = xbbar + A3 * sbar
    )

  return(stat_s)
}

# Computes control limits for a Standard Deviation (S) chart using B3 and B4 constants.
limits_s = function(x, y) {
  data = tibble(x = x, y = y)

  stat_s = data %>%
    group_by(x) %>%
    summarize(s = sd(y), nw = n(), df = nw - 1, .groups = "drop")

  constants = stat_s %>%
    dplyr::select(nw) %>%
    distinct() %>%
    group_by(nw) %>%
    summarize(bn(n = nw, reps = 1e4), .groups = "drop")

  stat_s = stat_s %>%
    left_join(by = "nw", y = constants) %>%
    mutate(
      sbar  = sqrt(sum(df * s^2) / sum(df)),
      lower = B3 * sbar,
      upper = B4 * sbar
    )

  return(stat_s)
}

# Computes control limits for a Range (R) chart using D3 and D4 constants.
limits_r = function(x, y) {
  data = tibble(x = x, y = y)

  stat_s = data %>%
    group_by(x) %>%
    summarize(r = max(y) - min(y), nw = n(), df = nw - 1, .groups = "drop")

  constants = stat_s %>%
    dplyr::select(nw) %>%
    distinct() %>%
    group_by(nw) %>%
    summarize(dn(n = nw, reps = 1e4), .groups = "drop")

  stat_s = stat_s %>%
    left_join(by = "nw", y = constants) %>%
    mutate(
      rbar  = mean(r),
      lower = D3 * rbar,
      upper = D4 * rbar
    )

  return(stat_s)
}

# Computes control limits for a Moving Range (MR) chart when subgroup size n = 1.
limits_mr = function(x, y) {
  data = tibble(x = x, y = y)

  data2 = data %>%
    reframe(
      x  = x[-1],
      mr = y %>% diff() %>% abs()
    )

  d2 = rnorm(n = 10000, mean = 0, sd = 1) %>% diff() %>% abs() %>% mean()

  stat = data2 %>%
    mutate(
      mrbar   = mean(mr),
      d2      = d2,
      sigma_s = mrbar / d2,
      n       = 1,
      se      = sigma_s / sqrt(n),
      upper   = mrbar + 3 * se,
      lower   = 0
    )

  return(stat)
}


# ==============================================================================
# Shewhart Control Charts (ggplot)
# ==============================================================================

# Plots an Average (Xbar) Control Chart with ribbon, line, points, and limit labels.
ggxbar = function(x, y, xlab = "Time (Subgroups)", ylab = "Average") {
  data   = tibble(x = x, y = y)
  stat_s = get_stat_s(x = data$x, y = data$y)

  labels = stat_s %>%
    reframe(
      x     = c(max(x), max(x), max(x)),
      type  = c("xbbar", "upper", "lower"),
      name  = c("xbbar", "+3 s", "-3 s"),
      value = c(mean(xbbar), max(upper), min(lower))
    ) %>%
    mutate(value = round(value, 2), text = paste0(name, " = ", value))

  stat_t = get_stat_t(x = x, y = y)

  gg = ggplot() +
    geom_hline(data = stat_t, mapping = aes(yintercept = xbbar), color = "lightgrey") +
    geom_ribbon(data = stat_s, mapping = aes(x = x, ymin = lower, ymax = upper), fill = "steelblue", alpha = 0.2) +
    geom_line(data = stat_s, mapping = aes(x = x, y = xbar), size = 1) +
    geom_point(data = stat_s, mapping = aes(x = x, y = xbar), size = 5) +
    geom_label(data = labels, mapping = aes(x = x, y = value, label = text), hjust = 1) +
    labs(x = xlab, y = ylab, subtitle = "Average Chart")

  return(gg)
}

# Plots a Standard Deviation (S) Control Chart with limits derived from B3/B4 constants.
ggs = function(x, y, xlab = "Time (Subgroups)", ylab = "Standard Deviation") {
  data   = tibble(x = x, y = y)
  stat_s = limits_s(x = data$x, y = data$y)
  stat_t = get_stat_t(x = data$x, y = data$y)

  labels = stat_s %>%
    reframe(
      x     = c(max(x), max(x), max(x)),
      type  = c("sbar", "upper", "lower"),
      name  = c("sbar", "+3 s", "-3 s"),
      value = c(mean(sbar), max(upper), min(lower))
    ) %>%
    mutate(value = round(value, 2), text = paste0(name, " = ", value))

  gg = ggplot() +
    geom_hline(data = stat_t, mapping = aes(yintercept = sbar), color = "lightgrey") +
    geom_ribbon(data = stat_s, mapping = aes(x = x, ymin = lower, ymax = upper), fill = "steelblue", alpha = 0.2) +
    geom_line(data = stat_s, mapping = aes(x = x, y = s), linewidth = 1) +
    geom_point(data = stat_s, mapping = aes(x = x, y = s), size = 5) +
    geom_label(data = labels, mapping = aes(x = x, y = value, label = text), hjust = 1) +
    labs(x = xlab, y = ylab, subtitle = "Standard Deviation Chart")

  return(gg)
}

# Plots a Range (R) Control Chart with limits derived from D3/D4 constants.
ggr = function(x, y, xlab = "Time (Subgroups)", ylab = "Range") {
  data   = tibble(x = x, y = y)
  stat_s = limits_r(x = data$x, y = data$y)
  stat_t = get_stat_t(x = data$x, y = data$y)

  labels = stat_s %>%
    reframe(
      x     = c(max(x), max(x), max(x)),
      type  = c("rbar", "upper", "lower"),
      name  = c("rbar", "+3 s", "-3 s"),
      value = c(mean(rbar), max(upper), min(lower))
    ) %>%
    mutate(value = round(value, 2), text = paste0(name, " = ", value))

  gg = ggplot() +
    geom_hline(data = stat_t, mapping = aes(yintercept = rbar), color = "lightgrey") +
    geom_ribbon(data = stat_s, mapping = aes(x = x, ymin = lower, ymax = upper), fill = "steelblue", alpha = 0.2) +
    geom_line(data = stat_s, mapping = aes(x = x, y = r), linewidth = 1) +
    geom_point(data = stat_s, mapping = aes(x = x, y = r), size = 5) +
    geom_label(data = labels, mapping = aes(x = x, y = value, label = text), hjust = 1) +
    labs(x = xlab, y = ylab, subtitle = "Range Chart")

  return(gg)
}

# Plots a Moving Range (MR) Control Chart for individual observations (subgroup size = 1).
ggmr = function(x, y, xlab = "Time (Subgroups)", ylab = "Moving Range") {
  data = tibble(x = x, y = y)

  data2 = data %>%
    reframe(x = x[-1], mr = y %>% diff() %>% abs())

  stat_s = limits_mr(x = data$x, y = data$y)

  labels = stat_s %>%
    reframe(
      x     = c(max(x), max(x), max(x)),
      type  = c("mrbar", "upper", "lower"),
      name  = c("mrbar", "+3 s", "-3 s"),
      value = c(mean(mrbar), max(upper), min(lower))
    ) %>%
    mutate(value = round(value, 2), text = paste0(name, " = ", value))

  stat_t = stat_s %>% summarize(mrbar = mean(mr))

  gg = ggplot() +
    geom_hline(data = stat_t, mapping = aes(yintercept = mrbar), color = "lightgrey") +
    geom_ribbon(data = stat_s, mapping = aes(x = x, ymin = lower, ymax = upper), fill = "steelblue", alpha = 0.2) +
    geom_line(data = data2, mapping = aes(x = x, y = mr), linewidth = 1) +
    geom_point(data = stat_s, mapping = aes(x = x, y = mr), size = 5) +
    geom_label(data = labels, mapping = aes(x = x, y = value, label = text), hjust = 1) +
    labs(x = xlab, y = ylab, subtitle = "Moving Range Chart")

  return(gg)
}


# ==============================================================================
# Attribute Control Charts (p, np, u)
# ==============================================================================

# Plots a Fraction Defective (p) Control Chart.
ggp = function(t, x, n, xlab = "Time (Subgroup)", ylab = "Fraction Defective") {
  data = tibble(t = t, x = x, n = n)

  stat_s = data %>%
    group_by(t) %>%
    mutate(p = x / n, mu = n * p, sigma = sqrt(n * p * (1 - p))) %>%
    ungroup() %>%
    mutate(
      xsum  = sum(x),
      nsum  = sum(n),
      pbar  = xsum / nsum,
      se    = sqrt(pbar * (1 - pbar) / n),
      lower = pmax(0, pbar - 3 * se),
      upper = pbar + 3 * se
    )

  gg = ggplot() +
    geom_ribbon(data = stat_s, mapping = aes(x = t, ymin = lower, ymax = upper), fill = "steelblue", alpha = 0.2) +
    geom_hline(data = stat_s, mapping = aes(yintercept = pbar), linewidth = 1.5, color = "darkgrey") +
    geom_line(data = stat_s, mapping = aes(x = t, y = p)) +
    geom_point(data = stat_s, mapping = aes(x = t, y = p)) +
    labs(x = xlab, y = ylab, subtitle = "Fraction Defective (p) Chart")

  return(gg)
}

# Plots a Number of Defectives (np) Control Chart.
ggnp = function(t, x, n, xlab = "Time (Subgroups)", ylab = "Number of Defectives (np)") {
  data = tibble(t = t, x = x, n = n)

  stat_s = data %>%
    group_by(t) %>%
    mutate(p = x / n, np = n * p) %>%
    ungroup() %>%
    mutate(
      xsum  = sum(x),
      nsum  = sum(n),
      npbar = sum(n * p) / n(),
      pbar  = sum(n * p) / sum(n),
      se    = sqrt(npbar * (1 - pbar)),
      lower = pmax(0, npbar - 3 * se),
      upper = npbar + 3 * se
    )

  labels = stat_s %>%
    reframe(
      t     = c(max(t), max(t), max(t)),
      type  = c("npbar", "upper", "lower"),
      name  = c("npbar", "+3 s", "-3 s"),
      value = c(mean(npbar), max(upper), min(lower))
    ) %>%
    mutate(value = round(value, 2), text = paste0(name, " = ", value))

  gg = ggplot() +
    geom_ribbon(data = stat_s, mapping = aes(x = t, ymin = lower, ymax = upper), fill = "steelblue", alpha = 0.2) +
    geom_hline(data = stat_s, mapping = aes(yintercept = npbar), linewidth = 1.5, color = "darkgrey") +
    geom_line(data = stat_s, mapping = aes(x = t, y = np)) +
    geom_point(data = stat_s, mapping = aes(x = t, y = np)) +
    geom_label(data = labels, mapping = aes(x = t, y = value, label = text), hjust = 1) +
    labs(x = xlab, y = ylab, subtitle = "Mean Defective (np) Chart")

  return(gg)
}

# Plots a Defects per Unit (u) Control Chart.
ggu = function(t, x, xlab = "Time (Subgroups)", ylab = "Number of Defects (u)") {
  data = tibble(t = t, x = x)

  stat_s = data %>%
    group_by(t) %>%
    mutate(u = sum(x), nw = n()) %>%
    ungroup() %>%
    mutate(
      ubar  = sum(u) / sum(nw),
      se    = sqrt(ubar / nw),
      lower = pmax(0, ubar - 3 * se),
      upper = ubar + 3 * se
    )

  labels = stat_s %>%
    reframe(
      t     = c(max(t), max(t), max(t)),
      type  = c("ubar", "upper", "lower"),
      name  = c("ubar", "+3 s", "-3 s"),
      value = c(mean(ubar), max(upper), min(lower))
    ) %>%
    mutate(value = round(value, 2), text = paste0(name, " = ", value))

  gg = ggplot() +
    geom_ribbon(data = stat_s, mapping = aes(x = t, ymin = lower, ymax = upper), fill = "steelblue", alpha = 0.2) +
    geom_hline(data = stat_s, mapping = aes(yintercept = ubar), linewidth = 1.5, color = "darkgrey") +
    geom_line(data = stat_s, mapping = aes(x = t, y = u)) +
    geom_point(data = stat_s, mapping = aes(x = t, y = u)) +
    geom_label(data = labels, mapping = aes(x = t, y = value, label = text), hjust = 1) +
    labs(x = xlab, y = ylab, subtitle = "Number of Defects (u) Chart")

  return(gg)
}


# ==============================================================================
# Capability & Performance Indices
# ==============================================================================

# Potential Process Capability Index (Cp): tolerance width / (6 * sigma_s)
# Measures potential capability of a centered, stable process using short-term variation.
cp = function(sigma_s, upper, lower) {
  abs(upper - lower) / (6 * sigma_s)
}

# Potential Process Performance Index (Pp): tolerance width / (6 * sigma_t)
# Measures potential performance of a centered process using total overall variation.
pp = function(sigma_t, upper, lower) {
  abs(upper - lower) / (6 * sigma_t)
}

# Actual Process Capability Index (Cpk): accounts for mean centering using short-term variation.
cpk = function(mu, sigma_s, lower = NULL, upper = NULL) {
  a = if (!is.null(lower)) abs(mu - lower) / (3 * sigma_s) else NULL
  b = if (!is.null(upper)) abs(upper - mu) / (3 * sigma_s) else NULL
  if (!is.null(lower) & !is.null(upper)) min(a, b) else if (is.null(lower)) b else a
}

# Actual Process Performance Index (Ppk): accounts for mean centering using total variation.
ppk = function(mu, sigma_t, lower = NULL, upper = NULL) {
  a = if (!is.null(lower)) abs(mu - lower) / (3 * sigma_t) else NULL
  b = if (!is.null(upper)) abs(upper - mu) / (3 * sigma_t) else NULL
  if (!is.null(lower) & !is.null(upper)) min(a, b) else if (is.null(lower)) b else a
}

# Resamples data to calculate empirical bootstrap confidence intervals for cp, pp, cpk, or ppk.
# Resampling can occur by subgroup (by_subgroup = TRUE) or individually (by_subgroup = FALSE).
get_index = function(x, y, index = "cp", upper, lower,
                      bootstrap_reps = 1000, ci_level = 0.95,
                      by_subgroup = TRUE) {
  if (!index %in% c("cp", "pp", "cpk", "ppk")) {
    stop("index must be one of: cp, pp, cpk, ppk")
  }
  if (index %in% c("cp", "pp")) {
    if (missing(upper) || missing(lower)) {
      stop("upper and lower specification limits are required for cp and pp")
    }
  } else {
    if (missing(upper) && missing(lower)) {
      stop("at least one of upper or lower specification limit is required for cpk and ppk")
    }
  }
  if (bootstrap_reps < 500) {
    warning("bootstrap_reps < 500 may result in unreliable confidence intervals")
  }

  data = tibble(x = x, y = y)

  # Observed point estimate
  if (index == "cp") {
    stat_s   = get_stat_s(x = x, y = y)
    estimate = cp(sigma_s = stat_s$sigma_s[1], upper = upper, lower = lower)
  } else if (index == "pp") {
    stat_t   = get_stat_t(x = x, y = y)
    estimate = pp(sigma_t = stat_t$sigma_t[1], upper = upper, lower = lower)
  } else if (index == "cpk") {
    stat_s   = get_stat_s(x = x, y = y)
    stat_t   = get_stat_t(x = x, y = y)
    estimate = cpk(mu = stat_t$xbbar[1], sigma_s = stat_s$sigma_s[1],
                   upper = if (missing(upper)) NULL else upper,
                   lower = if (missing(lower)) NULL else lower)
  } else if (index == "ppk") {
    stat_t   = get_stat_t(x = x, y = y)
    estimate = ppk(mu = stat_t$xbbar[1], sigma_t = stat_t$sigma_t[1],
                   upper = if (missing(upper)) NULL else upper,
                   lower = if (missing(lower)) NULL else lower)
  }

  boot_values = numeric(bootstrap_reps)

  for (i in 1:bootstrap_reps) {
    if (by_subgroup) {
      subgroups         = unique(data$x)
      sampled_subgroups = sample(subgroups, size = length(subgroups), replace = TRUE)
      boot_data         = tibble()
      for (sg in sampled_subgroups) {
        boot_data = bind_rows(boot_data, data %>% filter(x == sg))
      }
    } else {
      boot_data = data %>% slice_sample(n = nrow(data), replace = TRUE)
    }

    if (index == "cp") {
      boot_stat_s    = get_stat_s(x = boot_data$x, y = boot_data$y)
      boot_values[i] = cp(sigma_s = boot_stat_s$sigma_s[1], upper = upper, lower = lower)
    } else if (index == "pp") {
      boot_stat_t    = get_stat_t(x = boot_data$x, y = boot_data$y)
      boot_values[i] = pp(sigma_t = boot_stat_t$sigma_t[1], upper = upper, lower = lower)
    } else if (index == "cpk") {
      boot_stat_s    = get_stat_s(x = boot_data$x, y = boot_data$y)
      boot_stat_t    = get_stat_t(x = boot_data$x, y = boot_data$y)
      boot_values[i] = cpk(mu = boot_stat_t$xbbar[1], sigma_s = boot_stat_s$sigma_s[1],
                           upper = if (missing(upper)) NULL else upper,
                           lower = if (missing(lower)) NULL else lower)
    } else if (index == "ppk") {
      boot_stat_t    = get_stat_t(x = boot_data$x, y = boot_data$y)
      boot_values[i] = ppk(mu = boot_stat_t$xbbar[1], sigma_t = boot_stat_t$sigma_t[1],
                           upper = if (missing(upper)) NULL else upper,
                           lower = if (missing(lower)) NULL else lower)
    }
  }

  se        = sd(boot_values)
  alpha     = 1 - ci_level
  ci_bounds = quantile(boot_values, probs = c(alpha / 2, 1 - alpha / 2))

  output = tibble(
    term     = index,
    estimate = estimate,
    se       = se,
    lower    = ci_bounds[1],
    upper    = ci_bounds[2]
  )

  return(output)
}


# ==============================================================================
# Capability Summary & Confidence Interval Helpers
# ==============================================================================

# Computes Cp, Cpk, Pp, Ppk point estimates along with theoretical standard errors
# and 95% confidence intervals based on normal distribution theory.
calc_capability_summary = function(mu, sigma_s, sigma_t, lower, upper,
                                    n_subgroups, subgroup_size, n_total,
                                    conf_level = 0.95) {
  est_cp  = cp(sigma_s = sigma_s, upper = upper, lower = lower)
  est_pp  = pp(sigma_t = sigma_t, upper = upper, lower = lower)
  est_cpk = cpk(mu = mu, sigma_s = sigma_s, lower = lower, upper = upper)
  est_ppk = ppk(mu = mu, sigma_t = sigma_t, lower = lower, upper = upper)

  v_short = n_subgroups * (subgroup_size - 1)
  v_total = n_total - 1

  se_cp  = est_cp * sqrt(1 / (2 * v_short))
  se_pp  = est_pp * sqrt(1 / (2 * v_total))
  se_cpk = est_cpk * sqrt(1 / (2 * v_short) + 1 / (9 * n_total * (est_cpk^2)))
  se_ppk = est_ppk * sqrt(1 / (2 * v_total) + 1 / (9 * n_total * (est_ppk^2)))

  z_crit    = qnorm(1 - (1 - conf_level) / 2)
  estimates = c(est_cp, est_cpk, est_pp, est_ppk)
  ses       = c(se_cp, se_cpk, se_pp, se_ppk)
  lowers    = estimates - z_crit * ses
  uppers    = estimates + z_crit * ses

  tibble(
    index      = c("Cp", "Cpk", "Pp", "Ppk"),
    sigma_type = c("sigma_short", "sigma_short", "sigma_total", "sigma_total"),
    estimate   = round(estimates, 3),
    se         = round(ses, 3),
    lower_ci   = round(lowers, 3),
    upper_ci   = round(uppers, 3),
    label      = paste0(round(estimates, 2), " [", round(lowers, 2), ", ", round(uppers, 2), "]")
  )
}

# Formats numerical point estimates and confidence intervals as: "estimate [lower, upper]".
format_ci_label = function(estimate, lower, upper, digits = 2) {
  paste0(
    round(estimate, digits), " [",
    round(lower, digits), ", ",
    round(upper, digits), "]"
  )
}


# ==============================================================================
# Factorial Design Analysis Functions
# ==============================================================================

# Calculates the pooled standard error for a factorial model formula (e.g. y ~ a + b + c).
se_factorial = function(formula, data) {
  frame = model.frame(formula, data) %>% rename(y = 1)
  xvars = names(frame)[-1]

  output = frame %>%
    group_by(across(any_of(xvars))) %>%
    summarize(s = sd(y), n = n(), .groups = "drop") %>%
    summarize(se = sqrt(sum(s^2 / n))) %>%
    with(se)

  return(output)
}

# Calculates one-way main treatment effect (dbar) between High (level 2) and Low (level 1) groups.
dbar_oneway = function(formula, data) {
  frame = model.frame(formula, data) %>%
    dplyr::select(y = 1, a = 2) %>%
    mutate(a = factor(a)) %>%
    mutate(a = as.integer(a) - 1)

  frame %>%
    reframe(dbar = mean(y[a == 1] - y[a == 0])) %>%
    with(dbar)
}

# Calculates two-way interaction effect (dbar) between two factors.
# Compares combinations where factors are aligned (HH or LL) vs opposite (HL or LH).
dbar_twoway = function(formula, data) {
  frame = model.frame(formula, data) %>%
    dplyr::select(y = 1, a = 2, b = 3) %>%
    mutate(a = factor(a), b = factor(b)) %>%
    mutate(a = as.integer(a), b = as.integer(b))

  output = frame %>%
    reframe(
      x1 = y[(a == 2 & b == 2) | (a == 1 & b == 1)],
      x0 = y[(a == 1 & b == 2) | (a == 2 & b == 1)]
    ) %>%
    reframe(dbar = mean(x1) - mean(x0)) %>%
    with(dbar)

  return(output)
}

# Calculates three-way interaction effect (dbar) across three factors.
dbar_threeway = function(formula, data) {
  frame = model.frame(formula, data) %>%
    dplyr::select(y = 1, a = 2, b = 3, c = 4) %>%
    mutate(a = factor(a), b = factor(b), c = factor(c)) %>%
    mutate(a = as.integer(a) - 1, b = as.integer(b) - 1, c = as.integer(c) - 1)

  differences = frame %>%
    reframe(
      d1a = y[a == 1 & b == 0 & c == 1] - y[a == 0 & b == 0 & c == 1],
      d0a = y[a == 1 & b == 0 & c == 0] - y[a == 0 & b == 0 & c == 0],
      d1b = y[a == 1 & b == 1 & c == 1] - y[a == 0 & b == 1 & c == 1],
      d0b = y[a == 1 & b == 1 & c == 0] - y[a == 0 & b == 1 & c == 0]
    )

  output = differences %>%
    reframe(
      dbar_a = (mean(d1a) - mean(d0a)) / 2,
      dbar_b = (mean(d1b) - mean(d0b)) / 2,
      dbar   = (dbar_b - dbar_a) / 2
    ) %>%
    dplyr::select(dbar) %>%
    with(dbar)

  return(output)
}

# Constructs a formatted direct effects table from dbar differences and pooled SE.
calc_direct_effects_table = function(dbar_df, se_pooled, conf_level = 0.95) {
  z_crit = qnorm(1 - (1 - conf_level) / 2)
  dbar_df %>%
    mutate(
      se       = se_pooled,
      lower_95 = dbar - z_crit * se,
      upper_95 = dbar + z_crit * se,
      label    = format_ci_label(dbar, lower_95, upper_95, digits = 2)
    )
}
