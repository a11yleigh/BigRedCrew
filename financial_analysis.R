# financial_analysis.R -- Big Red Crew
# Matches the Financial Impacts section: placement correction factors with bootstrap 95% CIs,
# applied to an MLB-level throw. Raw Nm from our low-effort throws are never used directly.

library(dplyr); library(readr)
set.seed(5300)

# ---- 0. DATA: 3 pitchers x 5 placements x 10 throws = 150 throws ----
data = read_csv("throwData_10_30.csv", show_col_types = FALSE) %>%
  transmute(pitcher = factor(thrower_ID), placement = sensor_placement, arm_speed, torque = valgus_torque)
places = c("ideal", "rot_90", "rot_180", "long_displaced", "long_and_rot90")

# ---- 1. STATISTIC: reading-to-true ratio for each placement ----
# Step 1: from CORRECTLY placed throws only, model torque ~ pitcher + arm speed.
#         This predicts what each throw SHOULD have read.
# Step 2: ratio = actual reading / predicted true torque, averaged per placement.
#         ratio 0.78 means the sensor reads 22% low. Correction factor = 1 / ratio.
f = function(d){
  m = lm(torque ~ pitcher + arm_speed, data = filter(d, placement == "ideal"))
  d %>% mutate(ratio = torque / predict(m, newdata = d)) %>%
    group_by(placement) %>% summarize(ratio = mean(ratio), .groups = "drop")
}

# ---- 2. BOOTSTRAP: 1,000 resamples ----
# Within each pitcher x placement group (10 throws), redraw 10 throws WITH replacement,
# then rerun step 1 (including refitting the model).
boot = bind_rows(lapply(1:1000, function(i)
  data %>% group_by(pitcher, placement) %>% slice_sample(prop = 1, replace = TRUE) %>%
    ungroup() %>% f() %>% mutate(rep = i)))

# ---- 3. CORRECTION FACTORS with 95% CIs (middle 95% of the 1,000 results) ----
factors = boot %>% mutate(cf = 1 / ratio) %>% group_by(placement) %>%
  summarize(pct_error = median(ratio - 1),          # how far off the reading is
            factor    = median(cf),                 # multiply the reading by this
            lower     = quantile(cf, 0.025),
            upper     = quantile(cf, 0.975),
            .groups = "drop") %>%
  mutate(needs_correction = !(lower <= 1 & upper >= 1)) %>%   # CI excludes 1 = real bias
  arrange(match(placement, places))
print(factors, width = Inf)

# ---- 4. MLB EXAMPLE: what a true 120 Nm throw would read at each placement ----
# 120 Nm = top of the published 50-120 Nm range (Scientific Reports, 2023); >100 Nm hazardous (Slowik et al., 2019)
true_t = 120; hazard = 100
boot %>% mutate(reads = true_t * ratio) %>% group_by(placement) %>%
  summarize(reads_as = median(reads), lower = quantile(reads, 0.025), upper = quantile(reads, 0.975),
            looks_safe = median(reads) < hazard, .groups = "drop") %>%
  arrange(match(placement, places)) %>% print()

# ---- 5. INJURY COST CONTEXT (published values, scenario only) ----
salary = 5.34e6              # avg MLB salary 2026 (Blum, 2026)
surg_per_team = 276 / 10 / 30  # UCL surgeries per team-season, 2014-2024
tibble(at_risk_per_team_season = surg_per_team * salary,
       save_if_5pct_prevented  = 0.05 * surg_per_team * salary,
       save_if_10pct_prevented = 0.10 * surg_per_team * salary)
