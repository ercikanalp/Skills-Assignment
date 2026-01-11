# MRI Scheduling Case Study
# By: Alp Erçikan, Qiyuan Li, Laura Maas, Mart van der Vleuten
# Computational Research Skills (E&OR) (2025-2026-002-EBS4043)

library(tidyverse)
library(lubridate)

setwd("") # fill in

############################################################################
# 1. Setup, load and clean data 
############################################################################
path <- "" # fill in
raw <- readr::read_csv(path, show_col_types = FALSE)
stopifnot(all(c("Date","Time","Duration","PatientType") %in% names(raw)))

dat <- raw %>%
  mutate(
    Date = as.Date(Date),
    Time = as.numeric(Time),
    Duration = as.numeric(Duration),
    PatientType = as.factor(PatientType)
  ) %>%
  filter(
    !is.na(Date), !is.na(Time), !is.na(Duration), !is.na(PatientType),
    Duration > 0
  )

dat <- dat %>%
  mutate(
    seconds_since_midnight = round(Time * 3600),
    Timestamp = as.POSIXct(Date) + seconds(seconds_since_midnight),
    wday = wday(Date, label = TRUE, week_start = 1)
  ) %>%
  arrange(Timestamp)

# working-hour filter
dat <- dat %>%
  mutate(time_in_hours = hour(Timestamp) + minute(Timestamp)/60 + second(Timestamp)/3600) %>%
  filter(time_in_hours >= 8, time_in_hours <= 17)

data_type1 <- dat %>% filter(PatientType == "Type 1")
data_type2 <- dat %>% filter(PatientType == "Type 2")

cat("Rows total:", nrow(dat), "\n")
cat("Type 1:", nrow(data_type1), " | Type 2:", nrow(data_type2), "\n")

#Working-minutes clock
to_working_minutes <- function(df, day_start = 8, day_end = 17) {
  df <- df %>%
    arrange(Date, Time)
  
  date_tbl <- df %>%
    distinct(Date) %>%
    arrange(Date) %>%
    mutate(
      is_workday = !(wday(Date, week_start = 1) %in% c(6, 7)),
      workday_incr = as.integer(is_workday),
      workday_index = cumsum(workday_incr) - workday_incr
    )
  
  df %>%
    left_join(date_tbl, by = "Date") %>%
    mutate(
      minutes_in_day = pmin(pmax((Time - day_start) * 60, 0), (day_end - day_start) * 60),
      working_minutes = workday_index * ((day_end - day_start) * 60) + minutes_in_day
    )
}
############################################################################
# 2. EDA 
############################################################################
daily_counts <- dat %>%
  count(PatientType, Date, name = "n_calls") %>%
  arrange(PatientType, Date)

arrival_summary <- daily_counts %>%
  group_by(PatientType) %>%
  summarise(
    days_observed = n(),
    mean_calls_per_day = mean(n_calls),
    var_calls_per_day  = var(n_calls),
    sd_calls_per_day   = sd(n_calls),
    q50 = quantile(n_calls, 0.50),
    q90 = quantile(n_calls, 0.90),
    q95 = quantile(n_calls, 0.95),
    .groups = "drop"
  )
print(arrival_summary)

p_daily_counts <- ggplot(daily_counts, aes(x = Date, y = n_calls)) +
  geom_col() +
  facet_wrap(~ PatientType, scales = "free_y") +
  labs(title = "Daily Number of MRI Scan Requests Per Patient Type", x = NULL, y = "Calls per day")
print(p_daily_counts)

p_hist_counts <- ggplot(daily_counts, aes(x = n_calls)) +
  geom_histogram(bins = 20) +
  facet_wrap(~ PatientType, scales = "free_y") +
  labs(title = "Histogram of daily calls", x = "Calls per day", y = "Frequency")
print(p_hist_counts)

interarrival_by_type <- dat %>%
  group_split(PatientType) %>%
  purrr::map_dfr(~{
    d <- to_working_minutes(.x) %>% arrange(Date, Time)
    d %>%
      mutate(interarrival_min = c(NA, diff(working_minutes))) %>%
      dplyr::select(Date, Time, Timestamp, PatientType, interarrival_min)
  })


interarrival_summary <- interarrival_by_type %>%
  filter(!is.na(interarrival_min), interarrival_min >= 0) %>%
  group_by(PatientType) %>%
  summarise(
    n_gaps = n(),
    mean_gap_min = mean(interarrival_min),
    median_gap_min = median(interarrival_min),
    q90_gap_min = quantile(interarrival_min, 0.90),
    q95_gap_min = quantile(interarrival_min, 0.95),
    .groups = "drop"
  )
print(interarrival_summary)

# Histogram of Arrival times 
p_hist_gaps <- ggplot(
  interarrival_by_type %>% filter(!is.na(interarrival_min), interarrival_min >= 0),
  aes(x = interarrival_min)
) +
  geom_histogram(bins = 40) +
  facet_wrap(~ PatientType, scales = "free_y") +
  labs(title = "Inter-arrival times (minutes, excluding non-working hours)",
       x = "Minutes between consecutive calls", y = "Frequency")
print(p_hist_gaps) 

#Empirical Cumulative Distribution Function
p_ecdf_gaps <- ggplot(
  interarrival_by_type %>% filter(!is.na(interarrival_min), interarrival_min >= 0),
  aes(x = interarrival_min)
) +
  stat_ecdf(geom = "step") +
  facet_wrap(~ PatientType) +
  labs(title = "Empirical Distribution of interarrival times per patient type",
       x = "Minutes between consecutive calls", y = "F(t)")
print(p_ecdf_gaps)

duration_summary <- dat %>%
  group_by(PatientType) %>%
  summarise(
    n = n(),
    mean_h = mean(Duration),
    sd_h   = sd(Duration),
    q50_h  = quantile(Duration, 0.50),
    q75_h  = quantile(Duration, 0.75),
    q90_h  = quantile(Duration, 0.90),
    q95_h  = quantile(Duration, 0.95),
    max_h  = max(Duration),
    .groups = "drop"
  ) %>%
  mutate(across(ends_with("_h"), ~ .x * 60, .names = "{.col}_min"))
print(duration_summary)

p_hist_dur <- ggplot(dat, aes(x = Duration * 60)) +
  geom_histogram(bins = 40) +
  facet_wrap(~ PatientType, scales = "free_y") +
  labs(title = "Scan duration histogram", x = "Duration (minutes)", y = "Frequency")
print(p_hist_dur)

p_dens_dur <- ggplot(dat, aes(x = Duration * 60)) +
  geom_density() +
  facet_wrap(~ PatientType, scales = "free_y") +
  labs(title = "Scan duration density", x = "Duration (minutes)", y = "Density")
print(p_dens_dur)

p_box_dur <- ggplot(dat, aes(x = PatientType, y = Duration * 60)) +
  geom_boxplot() +
  labs(title = "Scan duration boxplot", x = NULL, y = "Duration (minutes)")
print(p_box_dur)

############################################################################
# 3. Hospital-relevant quantities 
############################################################################
hospital_spec <- list(
  duration_thresholds_min = c(30, 45, 60),
  duration_quantiles = c(0.50, 0.90, 0.95),
  arrivals_quantiles = c(0.50, 0.90, 0.95)
)

compute_hospital_metrics <- function(df, spec = hospital_spec) {
  daily <- df %>% count(Date, name = "arrivals")
  mean_arrivals <- mean(daily$arrivals)
  q_arrivals <- quantile(daily$arrivals, probs = spec$arrivals_quantiles, names = TRUE)
  
  dur_min <- df$Duration * 60
  mean_dur <- mean(dur_min)
  q_dur <- quantile(dur_min, probs = spec$duration_quantiles, names = TRUE)
  
  tail_probs <- map_dbl(spec$duration_thresholds_min, ~ mean(dur_min > .x)) %>%
    setNames(paste0("P(Duration>", spec$duration_thresholds_min, "min)"))
  
  list(
    arrivals = list(mean_per_day = mean_arrivals, quantiles = q_arrivals),
    durations = list(mean_min = mean_dur, quantiles_min = q_dur, tail_probs = tail_probs)
  )
}

metrics_type1 <- compute_hospital_metrics(data_type1, hospital_spec)
metrics_type2 <- compute_hospital_metrics(data_type2, hospital_spec)

metrics_table <- bind_rows(
  tibble(
    PatientType = "Type 1",
    mean_arrivals = metrics_type1$arrivals$mean_per_day,
    q50_arrivals = as.numeric(metrics_type1$arrivals$quantiles["50%"]),
    q90_arrivals = as.numeric(metrics_type1$arrivals$quantiles["90%"]),
    q95_arrivals = as.numeric(metrics_type1$arrivals$quantiles["95%"]),
    mean_duration_min = metrics_type1$durations$mean_min,
    q50_duration_min = as.numeric(metrics_type1$durations$quantiles_min["50%"]),
    q90_duration_min = as.numeric(metrics_type1$durations$quantiles_min["90%"]),
    q95_duration_min = as.numeric(metrics_type1$durations$quantiles_min["95%"]),
    p_gt_30 = as.numeric(metrics_type1$durations$tail_probs["P(Duration>30min)"]),
    p_gt_45 = as.numeric(metrics_type1$durations$tail_probs["P(Duration>45min)"]),
    p_gt_60 = as.numeric(metrics_type1$durations$tail_probs["P(Duration>60min)"])
  ),
  tibble(
    PatientType = "Type 2",
    mean_arrivals = metrics_type2$arrivals$mean_per_day,
    q50_arrivals = as.numeric(metrics_type2$arrivals$quantiles["50%"]),
    q90_arrivals = as.numeric(metrics_type2$arrivals$quantiles["90%"]),
    q95_arrivals = as.numeric(metrics_type2$arrivals$quantiles["95%"]),
    mean_duration_min = metrics_type2$durations$mean_min,
    q50_duration_min = as.numeric(metrics_type2$durations$quantiles_min["50%"]),
    q90_duration_min = as.numeric(metrics_type2$durations$quantiles_min["90%"]),
    q95_duration_min = as.numeric(metrics_type2$durations$quantiles_min["95%"]),
    p_gt_30 = as.numeric(metrics_type2$durations$tail_probs["P(Duration>30min)"]),
    p_gt_45 = as.numeric(metrics_type2$durations$tail_probs["P(Duration>45min)"]),
    p_gt_60 = as.numeric(metrics_type2$durations$tail_probs["P(Duration>60min)"])
  )
)
print(metrics_table)

############################################################################
# 4 & 5: Bootstrap Inference 
############################################################################
percentile_ci <- function(x, level = 0.95) {
  alpha <- (1 - level) / 2
  quantile(x, probs = c(alpha, 1 - alpha), na.rm = TRUE, names = FALSE)
}

daily_arrivals_vec <- function(df) {
  df %>% count(Date, name = "n_calls") %>% pull(n_calls)
}

mean_interarrival_working_min <- function(df) {
  d <- to_working_minutes(df) %>% arrange(Date, Time)
  gaps <- diff(d$working_minutes)
  mean(gaps[gaps >= 0], na.rm = TRUE)
}

set.seed(123)
B <- 5000

############################################################################
# 4.1 Type 1 arrivals (Poisson parametric bootstrap)
############################################################################
type1_daily <- daily_arrivals_vec(data_type1)
n_days1 <- length(type1_daily)
lambda_hat <- mean(type1_daily)

type1_arrivals_point <- tibble(
  PatientType = "Type 1",
  lambda_hat = lambda_hat,
  mean_arrivals_per_day = lambda_hat,
  mean_interarrival_min_model = 540 / lambda_hat,
  mean_interarrival_min_empirical = mean_interarrival_working_min(data_type1)
)
print(type1_arrivals_point)

boot_lambda <- replicate(B, {
  sim_daily <- rpois(n_days1, lambda_hat)
  mean(sim_daily)
})
boot_mean_gap <- 540 / boot_lambda

type1_arrivals_ci <- tibble(
  PatientType = "Type 1",
  metric = c("Arrivals/day (lambda)", "Mean inter-arrival (min, model-based)"),
  estimate = c(lambda_hat, 540 / lambda_hat),
  ci_low = c(percentile_ci(boot_lambda)[1], percentile_ci(boot_mean_gap)[1]),
  ci_high = c(percentile_ci(boot_lambda)[2], percentile_ci(boot_mean_gap)[2])
)
print(type1_arrivals_ci)

############################################################################
# 4.2: Type 1 durations (Normal parametric bootstrap)
############################################################################
dur1_min <- data_type1$Duration * 60
n1 <- length(dur1_min)
mu_hat <- mean(dur1_min)
sigma_hat <- sd(dur1_min)

boot_dur1 <- replicate(B, {
  sim <- rnorm(n1, mean= mu_hat, sd = sigma_hat)
  c(
    mean = mean(sim),
    q50 = as.numeric(quantile(sim, 0.50)),
    q90 = as.numeric(quantile(sim, 0.90)),
    q95 = as.numeric(quantile(sim, 0.95)),
    p30 = mean(sim > 30),
    p45 = mean(sim > 45),
    p60 = mean(sim > 60)
  )
})
boot_dur1 <- as.data.frame(t(boot_dur1))

type1_duration_ci <- tibble(
  PatientType = "Type 1",
  metric = c("Mean duration (min)", "q50 duration (min)", "q90 duration (min)", "q95 duration (min)",
             "P(duration>30)", "P(duration>45)", "P(duration>60)"),
  estimate = c(mu_hat,
              as.numeric(quantile(dur1_min, 0.50)),
              as.numeric(quantile(dur1_min, 0.90)),
              as.numeric(quantile(dur1_min, 0.95)),
              mean(dur1_min > 30),
              mean(dur1_min > 45),
              mean(dur1_min > 60)),
  ci_low = c(percentile_ci(boot_dur1$mean)[1],
             percentile_ci(boot_dur1$q50)[1],
             percentile_ci(boot_dur1$q90)[1],
             percentile_ci(boot_dur1$q95)[1],
             percentile_ci(boot_dur1$p30)[1],
             percentile_ci(boot_dur1$p45)[1],
             percentile_ci(boot_dur1$p60)[1]),
  ci_high = c(percentile_ci(boot_dur1$mean)[2],
              percentile_ci(boot_dur1$q50)[2],
              percentile_ci(boot_dur1$q90)[2],
              percentile_ci(boot_dur1$q95)[2],
              percentile_ci(boot_dur1$p30)[2],
              percentile_ci(boot_dur1$p45)[2],
              percentile_ci(boot_dur1$p60)[2])
)
print(type1_duration_ci)

############################################################################
# 5.1: Type 2 nonparametric bootstrap (days for arrivals, scans for duration)
############################################################################
type2_daily <- daily_arrivals_vec(data_type2)
n_days2 <- length(type2_daily)

boot_arr2 <- replicate(B, {
  resamp <- sample(type2_daily, size = n_days2, replace = TRUE)
  c(
    mean = mean(resamp),
    q50 = as.numeric(quantile(resamp, 0.50)),
    q90 = as.numeric(quantile(resamp, 0.90)),
    q95 =as.numeric(quantile(resamp, 0.95))
  )
})
boot_arr2 <- as.data.frame(t(boot_arr2))

type2_arrivals_ci <- tibble(
  PatientType = "Type 2",
  metric = c("Mean arrivals/day", "q50 arrivals/day", "q90 arrivals/day", "q95 arrivals/day"),
  estimate = c(mean(type2_daily),
               as.numeric(quantile(type2_daily, 0.50)),
               as.numeric(quantile(type2_daily, 0.90)),
               as.numeric(quantile(type2_daily, 0.95))),
  ci_low = c(percentile_ci(boot_arr2$mean)[1],
             percentile_ci(boot_arr2$q50)[1],
             percentile_ci(boot_arr2$q90)[1],
             percentile_ci(boot_arr2$q95)[1]),
  ci_high = c(percentile_ci(boot_arr2$mean)[2],
              percentile_ci(boot_arr2$q50)[2],
              percentile_ci(boot_arr2$q90)[2],
              percentile_ci(boot_arr2$q95)[2])
)
print(type2_arrivals_ci)

dur2_min <- data_type2$Duration * 60
n2 <- length(dur2_min)

boot_dur2 <- replicate(B, {
  resamp <- sample(dur2_min, size = n2, replace = TRUE)
  c(
    mean = mean(resamp),
    q50 = as.numeric(quantile(resamp, 0.50)),
    q90 = as.numeric(quantile(resamp, 0.90)),
    q95 = as.numeric(quantile(resamp, 0.95)),
    p30 = mean(resamp > 30),
    p45 = mean(resamp > 45),
    p60= mean(resamp > 60)
  )
})
boot_dur2 <- as.data.frame(t(boot_dur2))

type2_duration_ci <- tibble(
  PatientType = "Type 2",
  metric = c("Mean duration (min)", "q50 duration (min)", "q90 duration (min)", "q95 duration (min)",
             "P(duration>30)", "P(duration>45)", "P(duration>60)"),
  estimate = c(mean(dur2_min),
               as.numeric(quantile(dur2_min, 0.50)),
               as.numeric(quantile(dur2_min, 0.90)),
               as.numeric(quantile(dur2_min, 0.95)),
               mean(dur2_min > 30),
               mean(dur2_min > 45),
               mean(dur2_min > 60)),
  ci_low = c(percentile_ci(boot_dur2$mean)[1],
             percentile_ci(boot_dur2$q50)[1],
             percentile_ci(boot_dur2$q90)[1],
             percentile_ci(boot_dur2$q95)[1],
             percentile_ci(boot_dur2$p30)[1],
             percentile_ci(boot_dur2$p45)[1],
             percentile_ci(boot_dur2$p60)[1]),
  ci_high = c(percentile_ci(boot_dur2$mean)[2],
              percentile_ci(boot_dur2$q50)[2],
              percentile_ci(boot_dur2$q90)[2],
              percentile_ci(boot_dur2$q95)[2],
              percentile_ci(boot_dur2$p30)[2],
              percentile_ci(boot_dur2$p45)[2],
              percentile_ci(boot_dur2$p60)[2])
)
print(type2_duration_ci)

############################################################################
# 5.2: parametric comparison for Type 2 durations
############################################################################
do_parametric_compare <- TRUE
if (do_parametric_compare) {
  
  fit_lnorm <- MASS::fitdistr(dur2_min, densfun = "lognormal")
  meanlog_hat <- fit_lnorm$estimate["meanlog"]
  sdlog_hat <- fit_lnorm$estimate["sdlog"]
  
  fit_gamma <- MASS::fitdistr(dur2_min, densfun = "gamma")
  shape_hat <- fit_gamma$estimate["shape"]
  rate_hat <- fit_gamma$estimate["rate"]
  
  qq_plot <- function(x, qfun, title) {
    n <- length(x)
    p <- ppoints(n)
    theo <- qfun(p)
    dfqq <- tibble(theoretical = sort(theo), empirical = sort(x))
    ggplot(dfqq, aes(theoretical, empirical)) +
      geom_point() +
      geom_abline(slope = 1, intercept = 0) +
      labs(title = title, x = "Theoretical quantiles", y = "Empirical quantiles")
  }
  
  print(qq_plot(dur2_min, \(p) qlnorm(p, meanlog = meanlog_hat, sdlog = sdlog_hat),
                "Type 2 duration Q-Q plot: Lognormal fit"))
  print(qq_plot(dur2_min, \(p) qgamma(p, shape = shape_hat, rate = rate_hat),
                "Type 2 duration Q-Q plot: Gamma fit"))
  
  ks_distance <- function(x, cdf_fun) {
    xs <- sort(x)
    Fn <- (1:length(xs)) / length(xs)
    max(abs(Fn - cdf_fun(xs)))
  }
  
  ks_lnorm <- ks_distance(dur2_min, \(z) plnorm(z, meanlog = meanlog_hat, sdlog = sdlog_hat))
  ks_gamma <- ks_distance(dur2_min, \(z) pgamma(z, shape = shape_hat, rate = rate_hat))
  
  type2_parametric_fit_summary <- tibble(
    distribution = c("Lognormal", "Gamma"),
    ks_distance = c(ks_lnorm, ks_gamma)
  ) %>% arrange(ks_distance)
  
  print(type2_parametric_fit_summary)
}

############################################################################
# 6: Monte Carlo Robustness Study (Type 2 durations)
############################################################################
set.seed(2026)

metrics_names <- c("mean", "q90", "q95", "p_gt_45", "p_gt_60")

plugin_estimator <- function(x) {
  c(
    mean = mean(x),
    q90 = as.numeric(quantile(x, 0.90, type = 7)),
    q95 = as.numeric(quantile(x, 0.95, type = 7)),
    p_gt_45 = mean(x > 45),
    p_gt_60 = mean(x > 60)
  )
}

bootstrap_ci_np <- function(x, B = 800, level = 0.95) {
  alpha <- (1 - level) / 2
  n <- length(x)
  
  reps <- replicate(B, {
    xs <- sample(x, size = n, replace = TRUE)
    plugin_estimator(xs)
  })
  reps <- t(reps)
  
  est <- plugin_estimator(x)
  ci_low <- apply(reps, 2, quantile, probs = alpha, na.rm = TRUE)
  ci_high<- apply(reps, 2, quantile, probs = 1 - alpha, na.rm = TRUE)
  
  list(est = est, ci_low = ci_low, ci_high = ci_high)
}

# Calibrate truth models from empirical dur2_min
fit_gamma_truth <- MASS::fitdistr(dur2_min, densfun = "gamma")
gamma_shape <- unname(fit_gamma_truth$estimate["shape"])
gamma_rate  <- unname(fit_gamma_truth$estimate["rate"])

fit_lnorm_truth <- MASS::fitdistr(dur2_min, densfun = "lognormal")
ln_meanlog <- unname(fit_lnorm_truth$estimate["meanlog"])
ln_sdlog   <- unname(fit_lnorm_truth$estimate["sdlog"])

# Mixture (simple: 2 normals, truncated >0)
mix_w  <- 0.85
mix_mu1 <- as.numeric(quantile(dur2_min, 0.50))
mix_sd1 <- sd(dur2_min[dur2_min <= as.numeric(quantile(dur2_min, 0.80))])
mix_mu2 <- as.numeric(quantile(dur2_min, 0.95))
mix_sd2 <- sd(dur2_min[dur2_min >= as.numeric(quantile(dur2_min, 0.80))])

rmixture_pos <- function(n) {
  z <- rbinom(n, 1, mix_w)
  x <- numeric(n)
  n1 <- sum(z == 1)
  n0 <- n - n1
  if (n1 > 0) x[z == 1] <- rnorm(n1, mix_mu1, mix_sd1)
  if (n0 > 0) x[z == 0] <- rnorm(n0, mix_mu2, mix_sd2)
  x[x <= 0] <- 0.1
  x
}

truth_models <- list(
  Gamma = list(
    r = function(n) rgamma(n, shape = gamma_shape, rate = gamma_rate),
    true = function() plugin_estimator(rgamma(200000, shape = gamma_shape, rate = gamma_rate))
  ),
  Lognormal = list(
    r = function(n) rlnorm(n, meanlog = ln_meanlog, sdlog = ln_sdlog),
    true = function() plugin_estimator(rlnorm(200000, meanlog = ln_meanlog, sdlog = ln_sdlog))
  ),
  Mixture = list(
    r = function(n) rmixture_pos(n),
    true = function() plugin_estimator(rmixture_pos(200000))
  )
)

true_values <- lapply(truth_models, \(m) m$true())

M <- 300
Bboot <- 400
level <- 0.95

run_one_truth <- function(model_name) {
  gen <- truth_models[[model_name]]$r
  theta_true <- true_values[[model_name]]
  
  cover_mat <- matrix(FALSE, nrow = M, ncol = length(metrics_names),
                      dimnames = list(NULL, metrics_names))
  bias_mat  <- matrix(NA_real_, nrow = M, ncol = length(metrics_names),
                      dimnames = list(NULL, metrics_names))
  
  for (m in 1:M) {
    x <- gen(n2)
    out <- bootstrap_ci_np(x, B = Bboot, level = level)
    
    est <- out$est[metrics_names]
    lo  <- out$ci_low[metrics_names]
    hi  <- out$ci_high[metrics_names]
    tru <- theta_true[metrics_names]
    
    cover_mat[m, ] <- (tru >= lo) & (tru <= hi)
    bias_mat[m, ] <- est - tru
  }
  
  tibble(
    truth = model_name,
    metric = metrics_names,
    coverage = colMeans(cover_mat),
    mean_bias = colMeans(bias_mat),
    rmse = sqrt(colMeans(bias_mat^2))
  )
}

results <- bind_rows(lapply(names(truth_models), run_one_truth))

results_table <- results %>%
  mutate(
    coverage = round(coverage, 3),
    mean_bias = round(mean_bias, 3),
    rmse = round(rmse, 3)
  ) %>%
  arrange(truth, metric)

print(results_table)

coverage_wide <- results_table %>%
  dplyr::select(truth, metric, coverage) %>%
  tidyr::pivot_wider(names_from = metric, values_from = coverage)

print(coverage_wide)
