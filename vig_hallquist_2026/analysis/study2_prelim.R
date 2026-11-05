library(dplyr)
library(ggplot2)
library(ggh4x)

source_helpers <- function() {
  analysis_dir <- here::here("vig_hallquist_2026", "analysis")
  source(file.path(analysis_dir, "data_import.R"))
  source(here::here("vig_hallquist_2026", "vh_study2.R"))
}

source_helpers()

# condition_summary <- get_condition_summary("vig_hallquist_slurm", filter = study == "study2")
# manifest <- get_manifest("vig_hallquist_slurm", study_arg = "all")
study2_results <- read_study_cache("study2", this_run = "vig_hallquist_slurm") %>%
  mutate(
      status10_failure = !is.na(status_code) & status_code == 10L,
      converged = !status10_failure & !is.na(estimate),
      bias = estimate - truth,
      sq_error = (estimate - truth)^2,
      covered = ci_low <= truth & ci_high >= truth
  ) %>%
  add_study2_legacy_eligibility() %>%
  add_study2_analysis_eligibility() %>%
  mutate(analysis_ready = converged & analysis_eligible)

unique(study2_results$method)
names(study2_results)

study2_results %>%
  filter(analysis_ready) %>%
  filter(
    #marginal_rho == 0.5,
    standardized_beta_target == 0.4,
    method %in% c("fuller_closed_form", "lai_2spa", "closed_form_dual")
  ) %>%
  ggplot(aes(x = target_reliability, y = sq_error, fill = method)) +
  stat_summary(fun = function(x) sqrt(mean(x, na.rm = TRUE)), geom = "col", position = "dodge") +
  facet_nested("Slope-Intercept Correlation" + marginal_rho ~ "Number of Clusters" + as.factor(num_clus)) +
  labs(x = "Target Reliability", y = "RMSE", fill = "Method") +
  labs(title = "Study 2: RMSE", subtitle = "Beta = 0.4, Converged Only") +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom")

study2_results %>%
  filter(analysis_ready) %>%
  filter(
    #marginal_rho == 0.5,
    standardized_beta_target == 0.6,
    method %in% c("fuller_closed_form", "lai_2spa")
  ) %>%
  ggplot(aes(x = target_reliability, y = bias, fill = method)) +
  stat_summary(fun = safe_mean, geom = "col", position = "dodge") +
  facet_nested("Slope-Intercept Correlation" + marginal_rho ~ "Number of Clusters" + as.factor(num_clus)) +
  labs(x = "Target Reliability", y = "Bias", fill = "Method") +
  labs(title = "Study 2: Bias", subtitle = "Beta = 0.6, Converged Only") +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom")


study2_results %>%
  filter(
    method %in% c("fuller_closed_form", "lai_2spa", "msem")
  ) %>%
  ggplot(
    aes(
      x = target_reliability,
      y = analysis_ready,
      color = method,
      group = method
    )
  ) +
  stat_summary(fun = safe_mean, geom = "line", linewidth = 1) +
  stat_summary(fun = safe_mean, geom = "point", size = 2) +
  scale_y_continuous(
    labels = scales::label_percent(),
    limits = c(0, 1)
  ) +
  facet_nested(
    "Slope-Intercept Correlation" + marginal_rho ~ "Number of Clusters" + as.factor(num_clus)
  ) +
  labs(
    x = "Target Reliability",
    y = "Convergence rate",
    color = "Method"
  ) +
  theme_bw()
