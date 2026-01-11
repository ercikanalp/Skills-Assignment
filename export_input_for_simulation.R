source("StatisticsPart.R")

# Translate Results into Simulation Inputs 
# Slot length recommendations 
# You can switch between q90/q95:
SLOT_POLICY <- "q95"   # choose "q90" or "q95"

# Point-estimate quantiles from the observed data (minutes)
type1_q90 <- as.numeric(quantile(dur1_min, 0.90))
type1_q95 <- as.numeric(quantile(dur1_min, 0.95))
type2_q90 <- as.numeric(quantile(dur2_min, 0.90))
type2_q95 <- as.numeric(quantile(dur2_min, 0.95))

slot_type1_min <- if (SLOT_POLICY == "q90") type1_q90 else type1_q95
slot_type2_min <- if (SLOT_POLICY == "q90") type2_q90 else type2_q95

# Round to practical slot lengths (edit rounding if hospital uses different time blocks)
round_to <- function(x, base = 5) base * ceiling(x / base)

slot_recommendations <- tibble(
  PatientType = c("Type 1", "Type 2"),
  policy = SLOT_POLICY,
  slot_minutes_raw = c(slot_type1_min, slot_type2_min),
  slot_minutes_recommended = c(round_to(slot_type1_min, 5), round_to(slot_type2_min, 5))
)

cat("\n--- Phase 6.1 Slot recommendations (minutes) ---\n")
print(slot_recommendations)

# q90/q95 side-by-side table
slot_quantiles_table <- tibble(
  PatientType = c("Type 1", "Type 2"),
  q90_min = c(type1_q90, type2_q90),
  q95_min = c(type1_q95, type2_q95),
  q90_slot_rounded = c(round_to(type1_q90, 5), round_to(type2_q90, 5)),
  q95_slot_rounded = c(round_to(type1_q95, 5), round_to(type2_q95, 5))
)

cat("\n--- Phase 6.1 Quantile-based slot options (minutes) ---\n")
print(slot_quantiles_table)

# Reusable stochastic input generators
WORK_START_H <- 8
WORK_END_H   <- 17
WORK_MINUTES <- (WORK_END_H - WORK_START_H) * 60  # 540

# Type 1 arrivals: Poisson process over working day
# Uses exponential gaps with rate = lambda_per_day / WORK_MINUTES
generate_arrivals_type1 <- function(lambda_per_day,
                                    work_minutes = WORK_MINUTES) {
  stopifnot(lambda_per_day > 0, work_minutes > 0)
  
  rate_per_min <- lambda_per_day / work_minutes
  
  t <- 0
  arrivals <- numeric(0)
  while (TRUE) {
    t <- t + rexp(1, rate = rate_per_min)
    if (t > work_minutes) break
    arrivals <- c(arrivals, t)
  }
  arrivals  # minutes since 08:00
}
# Type 1 durations: Normal (truncate at small positive)

generate_duration_type1 <- function(n,
                                    mu_min,
                                    sigma_min) {
  stopifnot(n >= 0, sigma_min > 0)
  
  x <- rnorm(n, mean = mu_min, sd = sigma_min)
  x[x <= 0] <- 0.1
  x
}

# Type 2 arrivals: resample daily count + uniform times
# 1) sample a daily count from empirical daily counts
# 2) place arrivals uniformly across the working day
generate_arrivals_type2 <- function(empirical_daily_counts,
                                    work_minutes = WORK_MINUTES) {
  stopifnot(length(empirical_daily_counts) > 0, work_minutes > 0)
  
  n <- sample(empirical_daily_counts, size = 1, replace = TRUE)
  if (n == 0) return(numeric(0))
  sort(runif(n, min = 0, max = work_minutes))
}

#Type 2 durations: resample from empirical durations
generate_duration_type2 <- function(n,
                                    empirical_durations_min) {
  stopifnot(n >= 0, length(empirical_durations_min) > 0)
  
  sample(empirical_durations_min, size = n, replace = TRUE)
}

# Build empirical inputs once for convience 
# empirical inputs for Type 2
type2_daily_emp <- type2_daily 
dur2_emp_min <- dur2_min 

# estimated parametric inputs for Type 1
lambda1_hat <- lambda_hat # from 4.1 ("StatisticsPart.R)
mu1_hat <- mu_hat # from 4.2
sigma1_hat <- sigma_hat # from 4.2

# Sanity check: simulate one day
simulate_one_day_inputs <- function() {
  a1 <- generate_arrivals_type1(lambda1_hat)
  d1 <- generate_duration_type1(length(a1), mu1_hat, sigma1_hat)
  
  a2 <- generate_arrivals_type2(type2_daily_emp)
  d2 <- generate_duration_type2(length(a2), dur2_emp_min)
  
  tibble(
    PatientType = c(rep("Type 1", length(a1)), rep("Type 2", length(a2))),
    arrival_min_since_open = c(a1, a2),
    duration_min = c(d1, d2)
  ) %>% arrange(arrival_min_since_open)
}

cat("\n--- Phase 6.2 Sanity check: one simulated day of inputs ---\n")
print(simulate_one_day_inputs())
#View(simulate_one_day_inputs())

######################################################
# EXPORT SIMULATION INPUTS 
######################################################
# Safety checks
stopifnot(exists("lambda_hat"), exists("mu_hat"), exists("sigma_hat"))
stopifnot(exists("type2_daily"), exists("dur2_min"))

# 1. Type 1 parametric inputs
type1_params <- tibble(
  parameter = c("lambda_per_day", "mu_duration_min", "sigma_duration_min"),
  value = c(
    unname(lambda_hat),
    unname(mu_hat),
    unname(sigma_hat)
  )
)

write_csv(type1_params, "type1_params.csv")

# 2. Type 2 empirical daily arrivals
type2_daily_df <- tibble(
  daily_arrivals = as.integer(type2_daily)
)

write_csv(type2_daily_df, "type2_daily_arrivals.csv")

# 3. Type 2 empirical durations (minutes)
type2_durations_df <- tibble(
  duration_min = as.numeric(dur2_min)
)

write_csv(type2_durations_df, "type2_durations_min.csv")

# 4. Metadata
metadata <- tibble(
  key = c(
    "work_start_hour",
    "work_end_hour",
    "work_minutes",
    "slot_policy",
    "slot_minutes_type1",
    "slot_minutes_type2"
  ),
  value = c(
    8,
    17,
    540,
    SLOT_POLICY,
    round_to(slot_type1_min, 5),
    round_to(slot_type2_min, 5)
  )
)

write_csv(metadata, "simulation_metadata.csv")

cat("\nCSV files written\n")
