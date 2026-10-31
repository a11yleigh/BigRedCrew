# pulse_financial_analysis.R
# Financial impact of improved PULSE sensor placement -- Big Red Crew
# Structure follows Fraser's ZZ_financial_analysis.R:
#   statistic f() -> bootstrap -> cost -> CI -> compare scenarios
# Data: throwData_10_30.csv -- 3 pitchers x 5 placements x 10 throws = 150 throws

library(dplyr)
library(readr)
set.seed(5300)

# ---------------------------------------------------------------
# 0. DATA
data = read_csv("throwData_10.30.csv", show_col_types = FALSE) %>%
  transmute(pitcher = factor(thrower_ID),
            condition = sensor_placement,
            arm_speed, torque = valgus_torque)

# ---------------------------------------------------------------
# 1. ASSUMPTIONS (cite or justify each one in the report)
k           = 1.96          # failure band: +/- 1.96 residual SDs of ideal-placement throws
sensor_cost = 245           # PULSE Throw price, USD (Driveline Baseball, n.d.)
at_wage     = 60250 / 2080  # athletic trainer median wage, USD/hr (BLS, 2025)
sessions    = 100           # tracked sessions per pitcher-season         [ASSUMPTION]
min_session = 5             # AT minutes per session to don/sync/review   [ASSUMPTION]
fix_min     = 1             # extra donning minutes per session for a fix [ASSUMPTION]
fix_cost    = 10            # one-time cost of jig / alignment mark       [ASSUMPTION]
staff       = 15            # pitchers on staff                           [ASSUMPTION]

spend    = sensor_cost + sessions * (min_session / 60) * at_wage   # $ monitoring spend / pitcher-season
cost_fix = fix_cost   + sessions * (fix_min / 60) * at_wage        # $ added by a fix / pitcher-season

# ---------------------------------------------------------------
# 2. SCENARIOS: how often each placement occurs in practice
# BEFORE: uniform over the 5 tested placements (no field data on real donning error yet)
conds = c("ideal", "rot_90", "rot_180", "long_displaced", "long_and_rot90")
scen = list(
  before       = c(ideal=1, rot_90=1, rot_180=1, long_displaced=1, long_and_rot90=1),
  fix_rotation = c(ideal=1, rot_90=0, rot_180=0, long_displaced=1, long_and_rot90=0), # keyed sleeve / mark
  fix_length   = c(ideal=1, rot_90=1, rot_180=1, long_displaced=0, long_and_rot90=0), # epicondyle jig
  fix_both     = c(ideal=1, rot_90=0, rot_180=0, long_displaced=0, long_and_rot90=0)
)
scen = lapply(scen, function(w) w / sum(w))

# ---------------------------------------------------------------
# 3. STATISTIC
# Torque rises with arm speed, and arm speed differed across conditions, so the
# "true" reading for each throw is predicted from IDEAL-placement throws only:
#   torque ~ pitcher + arm_speed   (pitcher-specific baseline, shared arm-speed slope)
# A throw FAILS if its reading falls outside +/- k * residual SD of that ideal model.
# By construction ~5% of ideal throws fail (the measurement noise floor).
f = function(d, k){
  m     = lm(torque ~ pitcher + arm_speed, data = filter(d, condition == "ideal"))
  sigma = summary(m)$sigma
  pc = d %>%
    mutate(pred = predict(m, newdata = d),
           fail = abs(torque - pred) > k * sigma) %>%
    group_by(condition) %>%
    summarize(p = mean(fail), .groups = "drop")
  p = setNames(pc$p, pc$condition)[conds]
  bind_cols(as_tibble(as.list(p)), as_tibble(lapply(scen, function(w) sum(w * p))))
}
point = data %>% f(k = k)
print(point, width = Inf)

# ---------------------------------------------------------------
# 4. BOOTSTRAP: resample throws within each pitcher x condition cell, refit everything
B = 1000
boot = bind_rows(lapply(1:B, function(i){
  data %>% group_by(pitcher, condition) %>%
    slice_sample(prop = 1, replace = TRUE) %>% ungroup() %>%
    f(k = k) %>% mutate(rep = i)
}))

# ---------------------------------------------------------------
# 5. COSTS & 95% CIs
qi = function(x) tibble(estimate = median(x), se = sd(x),
                        lower = quantile(x, 0.025), upper = quantile(x, 0.975))
cond_rates = bind_rows(lapply(conds, function(c) qi(boot[[c]]) %>% mutate(condition = c))) %>%
  select(condition, everything())
print(cond_rates)

results = bind_rows(lapply(names(scen), function(s){
  p    = boot[[s]]
  cost = p * spend + if (s == "before") 0 else cost_fix
  save = boot$before * spend - cost
  bind_rows(
    qi(p)            %>% mutate(metric = "P(mis-measured throw)"),
    qi(cost)         %>% mutate(metric = "Cost per pitcher-season ($)"),
    qi(save)         %>% mutate(metric = "Net savings per pitcher-season ($)"),
    qi(save * staff) %>% mutate(metric = "Net savings per staff-season ($)")
  ) %>% mutate(scenario = s)
})) %>% select(scenario, metric, everything())
print(results, n = Inf, width = Inf)

# Share of bootstrap replicates in which each fix beats doing nothing
boot %>% summarize(across(c(fix_rotation, fix_length, fix_both),
                          ~ mean((before - .x) * spend - cost_fix > 0)))

# ---------------------------------------------------------------
# 6. SENSITIVITY to the failure band
bind_rows(lapply(c(1.645, 1.96, 2.576), function(kk)
  data %>% f(k = kk) %>% select(before:fix_both) %>% mutate(k = kk)))

# ---------------------------------------------------------------
# 7. Mean bias by placement (arm-speed adjusted), for the Results section
summary(lm(torque ~ pitcher + arm_speed + relevel(factor(condition), "ideal"), data = data))$coef

# ---------------------------------------------------------------
# 8. BREAK-EVEN vs UCL surgery cost (NOT a claimed saving)
tibble(source = c("Meldau et al. 2020 (2004-14)", "2014-24 study: $1.26B / 276"),
       cor = c(1.9e6, 1.26e9 / 276)) %>%
  mutate(staff_seasons_per_surgery_to_break_even = cor / (cost_fix * staff))
