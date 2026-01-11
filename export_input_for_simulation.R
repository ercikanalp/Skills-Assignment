source("StatisticsPart.R")

# EXPORT SIMULATION INPUTS 
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
