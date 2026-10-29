# REGRESSION ANALYSIS - DRIVELINE PULSE
#
# Research question:
# How does rotational sensor placement affect measured torque?
#
# We use raw torque in Nm, without normalization.
# Each thrower acts as a blocking variable.

# Load packages
library(dplyr)   # Data manipulation
library(readr)   # Read CSV files
library(broom)   # Regression results and confidence intervals
library(texreg)  # Regression comparison tables
library(ggplot2) # Optional visualizations

# Import the team's original data
raw_pulse = read_csv(
  "testData2.csv",
  na = c("", "NA")
)

# Check the dataset
glimpse(raw_pulse)

# Select the rotational experiments shared by all 3 throwers
regression_data = raw_pulse %>%

  filter(
    thrower_ID %in% c(1, 2, 3),

    # Throw IDs for the three rotational treatments
    throw_ID %in% c(1:10, 51:70),

    # Only the tested rotation angles
    rotational_pos %in% c(0, 90, 180),

    # Only observations with actual measurements
    !is.na(valgus_torque)
  ) %>%

  mutate(
    # Convert rotation to a category.
    # 0 degrees is the reference / ideal position.
    rotation = factor(
      rotational_pos,
      levels = c(0, 90, 180)
    ),

    # Identify each thrower.
    # Ally is the reference category.
    thrower = factor(
      thrower_ID,
      levels = c(1, 2, 3),
      labels = c("Ally", "Sophia", "Alexis")
    )
  )

# Check how many throws we have per person and rotation
regression_data %>%
  count(thrower, rotation)

# Confirm total observations
nrow(regression_data)

# Just to make sure it has 10 throws
# Stop the analysis if any group has fewer than 10 throws
group_check = regression_data %>%
  count(thrower, rotation)

stopifnot(
  nrow(regression_data) == 90,
  nrow(group_check) == 9,
  all(group_check$n == 10)
)

# Summarize torque for each thrower at each rotation
torque_summary = regression_data %>%

  group_by(thrower, rotation) %>%

  summarize(
    n = n(),

    mean_torque = mean(valgus_torque),

    sd_torque = sd(valgus_torque),

    mean_speed = mean(arm_speed),

    # Standard error of the torque mean
    se_torque = sd_torque / sqrt(n),

    # 95% confidence interval for the mean
    ci_lower = mean_torque -
      qt(0.975, df = n - 1) * se_torque,

    ci_upper = mean_torque +
      qt(0.975, df = n - 1) * se_torque,

    .groups = "drop"
  )

print(torque_summary, n = Inf)

# MODEL 1: Control for individual differences between throwers
m1 = regression_data %>%
  lm(formula = valgus_torque ~ rotation + thrower)

# View model results
summary(m1)

# Get coefficients, p-values and 95% confidence intervals
m1 %>%
  tidy(conf.int = TRUE)

# Get R-squared, adjusted R-squared and F-statistic
m1 %>%
  glance()

# MODEL 2:
# Does rotation affect each thrower differently?
m2 = regression_data %>%
  lm(formula = valgus_torque ~ rotation * thrower)

# Full model summary
summary(m2)

# Coefficients and confidence intervals
m2 %>%
  tidy(conf.int = TRUE)

# Model-level statistics
m2 %>%
  glance()

# Does adding interaction terms improve the model fit?
anova(m1, m2)


# --------------------------------------------------------
# SECONDARY ANALYSIS: ROTATION VS IDEAL POSITION
#
# Compare each thrower's torque at 90 and 180 degrees
# against their own ideal position at 0 degrees.
#
# Use Welch's unpaired t-test because the throws are not
# matched pairs and the groups may have different variances.
#
# These comparisons supplement our regression models.

# Function to compare one thrower at one rotation
compare_rotation = function(person, angle) {

  # Select the measurements for this thrower
  person_data = regression_data %>%
    filter(thrower == person)

  # Get the 10 torque measurements at ideal placement (0 degrees)
  ideal = person_data %>%
    filter(rotation == "0") %>%
    pull(valgus_torque)

  # Get the 10 torque measurements at the selected rotation
  rotated = person_data %>%
    filter(rotation == angle) %>%
    pull(valgus_torque)

  # Compare average torque between the rotated and ideal positions
  # Welch's test allows the two groups to have different variances
  result = t.test(
    rotated,
    ideal,
    paired = FALSE,
    var.equal = FALSE,
    conf.level = 0.95
  )

  # Save the mean difference, confidence interval, and p-value
  tibble(
    thrower = person,
    rotation = angle,

    # Positive means rotated torque is higher than ideal
    mean_difference = mean(rotated) - mean(ideal),

    # 95% confidence interval for the mean difference
    ci_lower = result$conf.int[1],
    ci_upper = result$conf.int[2],

    # P-value from Welch's test
    p_value = result$p.value
  )
}


# Compare 90 and 180 degrees against ideal for each thrower
rotation_effects = bind_rows(
  compare_rotation("Ally", "90"),
  compare_rotation("Ally", "180"),
  compare_rotation("Sophia", "90"),
  compare_rotation("Sophia", "180"),
  compare_rotation("Alexis", "90"),
  compare_rotation("Alexis", "180")
)


# Adjust p-values because we performed six comparisons
# Holm's method reduces the risk of false-positive findings
rotation_effects = rotation_effects %>%
  mutate(
    p_adjusted = p.adjust(p_value, method = "holm")
  )


# Display all six rotation comparisons
# The 95% confidence intervals are individual Welch intervals;
# they are not adjusted for multiple comparisons
print(rotation_effects, n = Inf)



# REGRESSION DIAGNOSTICS
#
# Check whether the assumptions of our regression model
# are reasonable before interpreting its results.
#
# Focus on unequal variances and unusual observations.


# Create a separate group for each thrower and rotation
# This gives us nine experimental groups
diagnostic_data = regression_data %>%
  mutate(
    treatment_group = interaction(thrower, rotation)
  )


# Test whether the nine groups have similar variances
# A small p-value suggests that variances differ across groups
fligner.test(
  valgus_torque ~ treatment_group,
  data = diagnostic_data
)
