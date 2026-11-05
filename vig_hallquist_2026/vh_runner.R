#' Generic runner and output aggregation for Vig-Hallquist (2026).
#'
#' Study-specific modules define the data generators and estimators; the runner
#' is responsible for selecting design conditions, executing replications,
#' writing per-condition artifacts, and rebuilding aggregate summaries from
#' completed condition files.

vh_pipeline_version <- function() {
  "openmx_one_row_stage1_covariance_v6_20260817"
}

#' Return the deterministic seed for one simulation replication.
#'
#' Ordinary conditions use their stable `condition_id`. Matched bridge and
#' sensitivity conditions instead use `simulation_seed_group`, which is shared
#' by the arms being compared. This implements common random numbers: the ICC
#' bridge has identical samples at its m = 10 anchor, and variance-ratio arms
#' reuse the same underlying random draws after their covariance-scale changes.
vh_replication_seed <- function(condition, rep_id,
                                base_seed = 20260612L,
                                condition_stride = 100000L) {
  seed_group <- if (
    "simulation_seed_group" %in% names(condition) &&
      length(condition$simulation_seed_group) > 0L &&
      is.finite(condition$simulation_seed_group[[1]])
  ) {
    as.numeric(condition$simulation_seed_group[[1]])
  } else {
    as.numeric(condition$condition_id[[1]])
  }
  seed <- as.numeric(base_seed) + seed_group * as.numeric(condition_stride) +
    as.numeric(rep_id)
  if (!is.finite(seed) || seed < 0 || seed > .Machine$integer.max) {
    stop("Derived replication seed is outside R's supported integer range.")
  }
  as.integer(seed)
}

#' Parse an optional command-line integer argument.
#'
#' @details
#' Shell scripts and batch schedulers often pass optional values as empty
#' strings or literal missing-value tokens. This helper normalizes those
#' sentinel values to `NA_integer_` so downstream validation can treat omitted
#' and explicitly missing inputs the same way.
#'
#' @param x A scalar-like value, usually a character vector element from
#'   `commandArgs()`.
#'
#' @return An integer scalar, or `NA_integer_` when `x` is absent or one of the
#'   recognized missing-value tokens.
parse_optional_integer_arg <- function(x) {
  # Accept the common string encodings used by shell wrappers for omitted args.
  if (is.null(x) || length(x) == 0L || is.na(x) || x %in% c("", "NA", "NaN", "NULL", "null")) {
    return(NA_integer_)
  }
  as.integer(x)
}

#' Parse a command-line logical flag.
#'
#' @param x A scalar-like value from `commandArgs()`.
#' @param default Logical value returned for omitted or missing-value tokens.
#'
#' @return A logical scalar.
parse_logical_arg <- function(x, default = FALSE) {
  if (is.null(x) || length(x) == 0L || is.na(x) || x %in% c("", "NA", "NaN", "NULL", "null")) {
    return(default)
  }
  if (is.logical(x)) {
    return(isTRUE(x))
  }
  x_chr <- tolower(trimws(as.character(x[[1]])))
  if (x_chr %in% c("1", "true", "t", "yes", "y")) {
    return(TRUE)
  }
  if (x_chr %in% c("0", "false", "f", "no", "n")) {
    return(FALSE)
  }
  stop("Could not parse logical argument: ", x)
}

#' Restrict a simulation design to one chunk.
#'
#' @details
#' Chunking lets large condition grids be split across independent jobs without
#' changing the underlying design identifiers. The input design is sliced by
#' row position after `select_design()` has applied any study and
#' `max_conditions` filters. The returned data frame carries a `chunk_meta`
#' attribute used for file naming and progress reporting.
#'
#' If both chunk arguments are `NA`, the design is returned unchanged and the
#' metadata describes the full selected design. If one chunk argument is
#' supplied, the other must be supplied too.
#'
#' @param design A data frame or tibble returned by `select_design()`. Must
#'   contain `condition_id`.
#' @param chunk_index One-based chunk number to run, or `NA_integer_` for the
#'   full selected design.
#' @param chunk_size Number of design rows per chunk, or `NA_integer_` for the
#'   full selected design.
#'
#' @return A sliced design with a `chunk_meta` attribute containing
#'   `chunk_index`, `chunk_size`, `condition_start`, `condition_end`, and
#'   `n_conditions`.
slice_design_chunk <- function(design, chunk_index = NA_integer_, chunk_size = NA_integer_) {
  if (is.na(chunk_index) && is.na(chunk_size)) {
    attr(design, "chunk_meta") <- list(
      chunk_index = NA_integer_,
      chunk_size = NA_integer_,
      condition_start = min(design$condition_id),
      condition_end = max(design$condition_id),
      n_conditions = nrow(design)
    )
    return(design)
  }

  if (is.na(chunk_index) || is.na(chunk_size)) {
    stop("`chunk_index` and `chunk_size` must be supplied together.")
  }
  if (chunk_index < 1L || chunk_size < 1L) {
    stop("`chunk_index` and `chunk_size` must be positive integers.")
  }

  # Chunks are contiguous 1-based windows over the selected design, not over
  # the unfiltered full factorial grid.
  start_idx <- ((chunk_index - 1L) * chunk_size) + 1L
  end_idx <- min(nrow(design), chunk_index * chunk_size)
  if (start_idx > nrow(design)) {
    out <- design[0, , drop = FALSE]
    attr(out, "chunk_meta") <- list(
      chunk_index = chunk_index,
      chunk_size = chunk_size,
      condition_start = NA_integer_,
      condition_end = NA_integer_,
      n_conditions = 0L
    )
    return(out)
  }

  out <- design %>%
    dplyr::slice(start_idx:end_idx)

  attr(out, "chunk_meta") <- list(
    chunk_index = chunk_index,
    chunk_size = chunk_size,
    condition_start = min(out$condition_id),
    condition_end = max(out$condition_id),
    n_conditions = nrow(out)
  )

  out
}

#' Create the stable label used for chunk-specific aggregate files.
#'
#' @param chunk_meta The `chunk_meta` attribute produced by
#'   `slice_design_chunk()`.
#'
#' @return A character scalar. Full-design runs return `"full_selection"`;
#'   chunked runs return a label containing the zero-padded chunk number and
#'   first/last condition identifiers.
make_chunk_label <- function(chunk_meta) {
  if (is.null(chunk_meta) || is.na(chunk_meta$chunk_index)) {
    "full_selection"
  } else {
    sprintf(
      "chunk_%03d_conditions_%04d_%04d",
      chunk_meta$chunk_index,
      chunk_meta$condition_start,
      chunk_meta$condition_end
    )
  }
}

#' Build per-condition output paths.
#'
#' @details
#' Condition-level outputs live under `out_dir/conditions`. The replication
#' results are compressed because they can be large; summaries are small CSV
#' files intended for quick inspection.
#'
#' @param out_dir Root output directory for the simulation run.
#' @param condition_id Integer condition identifier.
#'
#' @return A named list with paths for `replications`, `summary`,
#'   `issue_summary`, and `stage1_summary`.
get_condition_file_paths <- function(out_dir, condition_id) {
  condition_dir <- file.path(out_dir, "conditions")
  dir.create(condition_dir, recursive = TRUE, showWarnings = FALSE)

  list(
    replications = file.path(condition_dir, sprintf("condition_%04d_replication_results.csv.gz", condition_id)),
    summary = file.path(condition_dir, sprintf("condition_%04d_summary.csv", condition_id)),
    issue_summary = file.path(condition_dir, sprintf("condition_%04d_issue_summary.csv", condition_id)),
    stage1_summary = file.path(condition_dir, sprintf("condition_%04d_stage1_problem_summary.csv", condition_id))
  )
}

# NA-stable helpers for the inferential Monte Carlo summaries below.
vh_safe_sd <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2L) NA_real_ else stats::sd(x)
}

vh_safe_quantile <- function(x, probability) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) {
    return(NA_real_)
  }
  unname(stats::quantile(x, probs = probability, names = FALSE, na.rm = TRUE))
}

vh_safe_correlation <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 2L || stats::sd(x[keep]) <= sqrt(.Machine$double.eps) ||
      stats::sd(y[keep]) <= sqrt(.Machine$double.eps)) {
    return(NA_real_)
  }
  stats::cor(x[keep], y[keep])
}

vh_binomial_rate <- function(event, eligible) {
  keep <- !is.na(eligible) & eligible
  if (!any(keep) || anyNA(event[keep])) {
    return(NA_real_)
  }
  mean(event[keep])
}

vh_binomial_mc_se <- function(event, eligible) {
  keep <- !is.na(eligible) & eligible
  n <- sum(keep)
  rate <- vh_binomial_rate(event, eligible)
  if (n == 0L || !is.finite(rate)) {
    return(NA_real_)
  }
  sqrt(rate * (1 - rate) / n)
}

vh_binomial_wilson_bound <- function(event, eligible, bound = c("low", "high"),
                                     confidence = 0.95) {
  bound <- match.arg(bound)
  keep <- !is.na(eligible) & eligible
  n <- sum(keep)
  rate <- vh_binomial_rate(event, eligible)
  if (n == 0L || !is.finite(rate)) {
    return(NA_real_)
  }
  z <- stats::qnorm(1 - (1 - confidence) / 2)
  denominator <- 1 + z^2 / n
  center <- (rate + z^2 / (2 * n)) / denominator
  half_width <- z * sqrt(rate * (1 - rate) / n + z^2 / (4 * n^2)) /
    denominator
  if (identical(bound, "low")) {
    max(0, center - half_width)
  } else {
    min(1, center + half_width)
  }
}

vh_compact_reason_counts <- function(reason, eligible) {
  excluded <- !is.na(eligible) & !eligible
  reason <- as.character(reason[excluded])
  reason[is.na(reason) | !nzchar(reason)] <- "unspecified"
  if (length(reason) == 0L) {
    return(NA_character_)
  }
  counts <- sort(table(reason), decreasing = TRUE)
  paste0(names(counts), "=", as.integer(counts), collapse = "; ")
}

#' Summarize replication-level estimator results.
#'
#' @details
#' The summary is computed at the condition-study-method level while preserving
#' the design parameters needed to compare rows across the study grids.
#' Point-estimate summaries use `point_eligible`; interval summaries use the
#' stricter `interval_eligible`. Conditional coverage/rejection rates and joint
#' success probabilities are both reported so a method cannot appear well
#' calibrated solely because difficult replications failed to yield intervals.
#' Bias, RMSE, coverage, Type I error, and power use the standardized estimand
#' stored in `truth`.
#'
#' @param results Replication-level result tibble produced by
#'   `run_condition_replications()` or loaded from per-condition files. Expected
#'   columns include `condition_id`, `study`, `method`, `estimate`, `truth`,
#'   confidence limits, status diagnostics, explicit point/interval eligibility,
#'   and the design descriptors. Older inputs without the explicit fields are
#'   classified conservatively from their status, estimates, SEs, and limits.
#'
#' @return A tibble with one row per condition-study-method combination and
#'   Monte Carlo summary columns for eligibility, bias and its MC interval,
#'   empirical-versus-model SE calibration on the common interval-eligible
#'   subset, interval-width tails, conditional and joint coverage,
#'   rejection/Type-I error/power, Wilson intervals for binomial rates, and
#'   estimator diagnostic averages.
summarize_results_df <- function(results) {
  required_point_columns <- c("estimate", "truth", "status_code")
  missing_point_columns <- setdiff(required_point_columns, names(results))
  if (length(missing_point_columns) > 0L) {
    stop(
      "Cannot summarize results without: ",
      paste(missing_point_columns, collapse = ", ")
    )
  }
  # Point-only diagnostic methods and older fixtures may omit interval fields.
  # Materialize them as missing so those rows remain available for point
  # summaries and are transparently excluded from interval summaries.
  for (column in c("se", "ci_low", "ci_high")) {
    if (!(column %in% names(results))) {
      results[[column]] <- rep(NA_real_, nrow(results))
    }
  }

  # The explicit fields are written by current Studies 1--4. Reclassifying here
  # also makes aggregate rebuilds from older condition files conservative and
  # gives them the same denominator contract.
  results <- add_vh_analysis_eligibility(results)

  design_cols <- intersect(
    c(
      "condition_id", "study", "study_version", "calibration_version",
      "simulation_module", "sensitivity_block", "sensitivity_block_label",
      "variance_ratio_arm", "variance_ratio_pair_id",
      "variance_ratio_pair_label", "sensitivity_reference_ratio",
      "sensitivity_is_low_reliability_stress",
      "method", "method_role", "fuller_variant",
      "fuller_preliminary_moment", "fuller_variance_bread",
      "fuller_predictor_outcome_covariance_source",
      "num_clus", "mean_clus_size",
      "target_reliability", "achieved_reliability", "marginal_rho",
      "calibration_arm", "calibration_metric", "calibration_target",
      "bridge_pair_id", "bridge_pair_label", "simulation_seed_group",
      "calibration_target_value", "achieved_calibration_value",
      "covariance_shape_fixed", "posterior_reliability_anchor", "icc_anchor",
      "target_calibration_reliability", "achieved_calibration_reliability",
      "achieved_partial_reliability", "slope_variance_marginal",
      "residualized_slope_variance", "slope_intercept_variance_ratio",
      "G_condition_number", "intercept_icc", "marginal_slope_icc",
      "conditional_slope_icc",
      "standardized_beta_target", "structural_target", "structural_r2",
      "focal_unique_r2", "mean_clus_size_y", "mean_clus_size_q",
      "target_reliability_y", "target_reliability_q",
      "achieved_reliability_y", "achieved_reliability_q",
      "achieved_partial_reliability_y", "achieved_partial_reliability_q",
      "intercept_icc_y", "intercept_icc_q",
      "conditional_slope_icc_y", "conditional_slope_icc_q",
      "balance_mode", "r_structure", "r_rho", "sigma", "sigma_y", "sigma_q",
      "information_profile", "is_falsification_control", "information_matched",
      "profile_min_clus_size", "profile_max_clus_size",
      "reliability_sd", "reliability_iqr", "reliability_min",
      "reliability_max", "reliability_small", "reliability_large",
      "population_lambda22_mean", "population_lambda22_min",
      "population_lambda22_max",
      "population_lambda_matrix_frobenius_rms_dispersion",
      "population_theta22_mean", "population_theta22_min",
      "population_theta22_max",
      "population_theta_matrix_frobenius_rms_dispersion",
      "population_ols_var22_mean", "population_ols_var22_min",
      "population_ols_var22_max",
      "population_ols_cov_matrix_frobenius_rms_dispersion"
    ),
    names(results)
  )
  replication_diagnostic_cols <- intersect(
    c(
      "realized_reliability_mean", "realized_reliability_sd",
      "realized_reliability_iqr", "realized_reliability_min",
      "realized_reliability_max",
      "lambda22_mean", "lambda22_sd", "lambda22_min", "lambda22_max",
      "lambda_matrix_frobenius_rms_dispersion",
      "theta22_mean", "theta22_sd", "theta22_min", "theta22_max",
      "theta_matrix_frobenius_rms_dispersion",
      "ols_var22_mean", "ols_var22_sd", "ols_var22_min", "ols_var22_max",
      "ols_cov_matrix_frobenius_rms_dispersion",
      "blup_slope_bias_small", "blup_slope_rmse_small",
      "blup_slope_bias_large", "blup_slope_rmse_large",
      "corrected_slope_bias_small", "corrected_slope_rmse_small",
      "corrected_slope_bias_large", "corrected_slope_rmse_large",
      "average_measurement_slope_bias_small",
      "average_measurement_slope_rmse_small",
      "average_measurement_slope_bias_large",
      "average_measurement_slope_rmse_large",
      "fuller_measurement_weight_used",
      "fuller_predictor_outcome_covariance_max_abs",
      "fuller_alpha_step1_used", "fuller_alpha_step3_used",
      "fuller_alpha_scaling_used",
      "fuller_correction1", "fuller_correction_c",
      "fuller_correction_scaling",
      "fuller_sx1_star_condition", "fuller_sx1_star_min_eigen",
      "fuller_sx1_star_relative_min_eigen",
      "fuller_sx_star_condition", "fuller_sx_star_min_eigen",
      "fuller_sx_star_relative_min_eigen",
      "fuller_scaling_condition", "fuller_scaling_min_eigen",
      "fuller_scaling_relative_min_eigen",
      "fuller_auto_guard_pass", "fuller_auto_full_weight_guard_pass",
      "fuller_auto_search_evaluations"
    ),
    names(results)
  )
  # Automatically retain the expanded replication-level audit fields. Keeping
  # this prefix-based prevents a newly added matrix entry from silently being
  # written to condition files but omitted from aggregate summaries.
  expanded_diagnostic_cols <- grep(
    paste0(
      "^stage1(_[yq])?_(lambda|theta|posterior_variance|",
      "corrected_score_error_covariance|fitted_|true_blup|blup_slope|",
      "corrected_score_slope|latent_covariance)|",
      "^mx_(raw_focal|latent_|predictor_latent_|outcome_latent_)|",
      "^mplus_"
    ),
    names(results),
    value = TRUE
  )
  expanded_diagnostic_cols <- expanded_diagnostic_cols[
    vapply(results[expanded_diagnostic_cols], function(x) {
      is.numeric(x) || is.integer(x) || is.logical(x)
    }, logical(1L))
  ]
  replication_diagnostic_cols <- unique(c(
    replication_diagnostic_cols,
    expanded_diagnostic_cols
  ))
  results %>%
    dplyr::mutate(
      # Status 10 is an OpenMx optimizer failure class; keep it visible rather
      # than folding it into generic missing-estimate behavior.
      status10_failure = !is.na(status_code) & status_code == 10L,
      converged = point_eligible,
      signed_error = dplyr::if_else(
        point_eligible,
        estimate - truth,
        NA_real_
      ),
      abs_error = abs(signed_error),
      sq_error = signed_error^2,
      eligible_estimate = dplyr::if_else(point_eligible, estimate, NA_real_),
      interval_eligible_estimate = dplyr::if_else(
        interval_eligible,
        estimate,
        NA_real_
      ),
      eligible_se = dplyr::if_else(interval_eligible, se, NA_real_),
      interval_width = dplyr::if_else(
        interval_eligible,
        ci_high - ci_low,
        NA_real_
      ),
      standardized_error = dplyr::if_else(
        interval_eligible,
        (estimate - truth) / se,
        NA_real_
      ),
      covered = ci_low <= truth & ci_high >= truth,
      reject_zero = ci_low > 0 | ci_high < 0,
      reject_positive = ci_low > 0,
      reject_negative = ci_high < 0,
      success_and_cover_event = interval_eligible & dplyr::coalesce(covered, FALSE),
      success_and_reject_event = interval_eligible & dplyr::coalesce(reject_zero, FALSE)
    ) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(design_cols))) %>%
    dplyr::summarise(
      truth = dplyr::first(truth),
      null_condition = is.finite(dplyr::first(truth)) &&
        abs(dplyr::first(truth)) <= sqrt(.Machine$double.eps),
      n_rep = dplyr::n(),
      n_point_eligible = sum(point_eligible, na.rm = TRUE),
      n_interval_eligible = sum(interval_eligible, na.rm = TRUE),
      point_eligibility = safe_mean(point_eligible),
      interval_eligibility = safe_mean(interval_eligible),
      point_eligibility_mc_se = vh_binomial_mc_se(point_eligible, rep(TRUE, dplyr::n())),
      interval_eligibility_mc_se = vh_binomial_mc_se(interval_eligible, rep(TRUE, dplyr::n())),
      point_exclusion_reasons = vh_compact_reason_counts(
        point_exclusion_reason,
        point_eligible
      ),
      interval_exclusion_reasons = vh_compact_reason_counts(
        interval_exclusion_reason,
        interval_eligible
      ),
      # Backward-compatible aliases: "success" and "analysis eligibility"
      # now unambiguously mean availability of a defensible point estimate.
      convergence = point_eligibility,
      n_success = n_point_eligible,
      n_analysis_eligible = n_point_eligible,
      analysis_eligibility = point_eligibility,
      mean_estimate = safe_mean(eligible_estimate),
      empirical_sd = vh_safe_sd(eligible_estimate),
      mc_se_mean = if (n_point_eligible > 1L) {
        empirical_sd / sqrt(n_point_eligible)
      } else {
        NA_real_
      },
      mc_se_bias = mc_se_mean,
      bias = safe_mean(signed_error),
      bias_mc_low = if (is.finite(bias) && is.finite(mc_se_bias)) {
        bias - stats::qnorm(0.975) * mc_se_bias
      } else {
        NA_real_
      },
      bias_mc_high = if (is.finite(bias) && is.finite(mc_se_bias)) {
        bias + stats::qnorm(0.975) * mc_se_bias
      } else {
        NA_real_
      },
      rmse = if (all(is.na(sq_error))) NA_real_ else sqrt(mean(sq_error, na.rm = TRUE)),
      median_absolute_error = vh_safe_quantile(abs_error, 0.50),
      p95_absolute_error = vh_safe_quantile(abs_error, 0.95),
      p99_absolute_error = vh_safe_quantile(abs_error, 0.99),
      empirical_sd_interval_subset = vh_safe_sd(interval_eligible_estimate),
      mean_se = safe_mean(eligible_se),
      median_se = vh_safe_quantile(eligible_se, 0.50),
      p95_se = vh_safe_quantile(eligible_se, 0.95),
      p99_se = vh_safe_quantile(eligible_se, 0.99),
      max_se = if (all(is.na(eligible_se))) NA_real_ else max(eligible_se, na.rm = TRUE),
      mean_se_to_empirical_sd = if (is.finite(empirical_sd_interval_subset) &&
          empirical_sd_interval_subset > sqrt(.Machine$double.eps)) {
        mean_se / empirical_sd_interval_subset
      } else {
        NA_real_
      },
      mean_se_to_interval_empirical_sd = mean_se_to_empirical_sd,
      mean_interval_width = safe_mean(interval_width),
      median_interval_width = vh_safe_quantile(interval_width, 0.50),
      p95_interval_width = vh_safe_quantile(interval_width, 0.95),
      p99_interval_width = vh_safe_quantile(interval_width, 0.99),
      max_interval_width = if (all(is.na(interval_width))) {
        NA_real_
      } else {
        max(interval_width, na.rm = TRUE)
      },
      mean_standardized_error = safe_mean(standardized_error),
      sd_standardized_error = vh_safe_sd(standardized_error),
      q025_standardized_error = vh_safe_quantile(standardized_error, 0.025),
      q975_standardized_error = vh_safe_quantile(standardized_error, 0.975),
      signed_error_se_correlation = vh_safe_correlation(signed_error, eligible_se),
      absolute_error_se_correlation = vh_safe_correlation(abs_error, eligible_se),
      coverage = vh_binomial_rate(covered, interval_eligible),
      coverage_mc_se = vh_binomial_mc_se(covered, interval_eligible),
      coverage_mc_low = vh_binomial_wilson_bound(covered, interval_eligible, "low"),
      coverage_mc_high = vh_binomial_wilson_bound(covered, interval_eligible, "high"),
      success_and_cover = safe_mean(success_and_cover_event),
      conditional_rejection_rate = vh_binomial_rate(reject_zero, interval_eligible),
      rejection_rate_mc_se = vh_binomial_mc_se(reject_zero, interval_eligible),
      rejection_rate_mc_low = vh_binomial_wilson_bound(
        reject_zero,
        interval_eligible,
        "low"
      ),
      rejection_rate_mc_high = vh_binomial_wilson_bound(
        reject_zero,
        interval_eligible,
        "high"
      ),
      conditional_positive_rejection_rate = vh_binomial_rate(
        reject_positive,
        interval_eligible
      ),
      conditional_negative_rejection_rate = vh_binomial_rate(
        reject_negative,
        interval_eligible
      ),
      success_and_reject = safe_mean(success_and_reject_event),
      type1_error = if (null_condition) conditional_rejection_rate else NA_real_,
      type1_error_positive = if (null_condition) {
        conditional_positive_rejection_rate
      } else {
        NA_real_
      },
      type1_error_negative = if (null_condition) {
        conditional_negative_rejection_rate
      } else {
        NA_real_
      },
      type1_error_mc_se = if (null_condition) rejection_rate_mc_se else NA_real_,
      type1_error_mc_low = if (null_condition) rejection_rate_mc_low else NA_real_,
      type1_error_mc_high = if (null_condition) rejection_rate_mc_high else NA_real_,
      operational_type1_error = if (null_condition) success_and_reject else NA_real_,
      power = if (!null_condition) conditional_rejection_rate else NA_real_,
      power_mc_se = if (!null_condition) rejection_rate_mc_se else NA_real_,
      power_mc_low = if (!null_condition) rejection_rate_mc_low else NA_real_,
      power_mc_high = if (!null_condition) rejection_rate_mc_high else NA_real_,
      operational_power = if (!null_condition) success_and_reject else NA_real_,
      n_status10_fail = sum(status10_failure, na.rm = TRUE),
      prop_status10_fail = safe_mean(status10_failure),
      dplyr::across(
        dplyr::all_of(replication_diagnostic_cols),
        safe_mean,
        .names = "mean_{.col}"
      ),
      .groups = "drop"
    )
}

#' Summarize results by first-stage singularity/problem status.
#'
#' @details
#' Some simulation conditions can be dominated by first-stage mixed-model
#' singularity or design-matrix instability. When the replication results carry
#' those diagnostics, this helper stratifies the standard performance summaries
#' by `stage1_singular_problem` and adds aggregate first-stage diagnostic
#' measures. If the diagnostics are not present, it returns an empty tibble so
#' callers can always write a stable output file.
#'
#' @param results Replication-level result tibble.
#'
#' @return A tibble grouped by condition-study-method and
#'   `stage1_singular_problem`, or an empty tibble when first-stage problem
#'   diagnostics are unavailable.
summarize_stage1_problem_df <- function(results) {
  if (!("stage1_singular_problem" %in% names(results))) {
    return(tibble::tibble())
  }

  results <- add_vh_analysis_eligibility(results)

  design_cols <- intersect(
    c(
      "condition_id", "study", "study_version", "calibration_version",
      "simulation_module", "sensitivity_block", "sensitivity_block_label",
      "variance_ratio_arm", "variance_ratio_pair_id",
      "variance_ratio_pair_label", "sensitivity_reference_ratio",
      "sensitivity_is_low_reliability_stress",
      "method", "method_role", "fuller_variant",
      "fuller_preliminary_moment", "fuller_variance_bread",
      "fuller_predictor_outcome_covariance_source",
      "stage1_singular_problem",
      "num_clus", "mean_clus_size", "target_reliability",
      "achieved_reliability", "marginal_rho", "standardized_beta_target",
      "calibration_arm", "calibration_metric", "calibration_target",
      "bridge_pair_id", "bridge_pair_label", "simulation_seed_group",
      "calibration_target_value", "achieved_calibration_value",
      "covariance_shape_fixed", "posterior_reliability_anchor", "icc_anchor",
      "target_calibration_reliability", "achieved_calibration_reliability",
      "achieved_partial_reliability", "slope_variance_marginal",
      "residualized_slope_variance", "slope_intercept_variance_ratio",
      "G_condition_number", "intercept_icc", "marginal_slope_icc",
      "conditional_slope_icc",
      "structural_target", "structural_r2", "focal_unique_r2",
      "mean_clus_size_y", "mean_clus_size_q",
      "target_reliability_y", "target_reliability_q",
      "achieved_reliability_y", "achieved_reliability_q",
      "balance_mode", "r_structure", "r_rho", "sigma", "sigma_y", "sigma_q",
      "information_profile", "is_falsification_control", "information_matched",
      "profile_min_clus_size", "profile_max_clus_size",
      "reliability_sd", "reliability_min", "reliability_max"
    ),
    names(results)
  )

  results <- results %>%
    dplyr::mutate(
      status10_failure = !is.na(status_code) & status_code == 10L,
      converged = point_eligible,
      eligible_estimate = dplyr::if_else(point_eligible, estimate, NA_real_),
      bias = dplyr::if_else(point_eligible, estimate - truth, NA_real_),
      sq_error = bias^2,
      covered = ci_low <= truth & ci_high >= truth,
      success_and_cover_event = interval_eligible & dplyr::coalesce(covered, FALSE)
    ) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(design_cols))) %>%
    dplyr::summarise(
      truth = dplyr::first(truth),
      n_rep = dplyr::n(),
      n_point_eligible = sum(point_eligible, na.rm = TRUE),
      n_interval_eligible = sum(interval_eligible, na.rm = TRUE),
      point_eligibility = safe_mean(point_eligible),
      interval_eligibility = safe_mean(interval_eligible),
      convergence = point_eligibility,
      n_success = n_point_eligible,
      mean_estimate = safe_mean(eligible_estimate),
      bias = safe_mean(bias),
      coverage = vh_binomial_rate(covered, interval_eligible),
      success_and_cover = safe_mean(success_and_cover_event),
      rmse = if (all(is.na(sq_error))) NA_real_ else sqrt(mean(sq_error, na.rm = TRUE)),
      prop_status10_fail = safe_mean(status10_failure),
      prop_stage1_lmer_singular = safe_mean(stage1_lmer_singular),
      mean_stage1_re_corr = safe_mean(stage1_re_corr),
      mean_stage1_eb_corr = safe_mean(stage1_eb_corr),
      mean_stage1_design_kappa = safe_mean(stage1_design_kappa),
      sample_stage1_problem_detail = compact_message(stage1_problem_detail[!is.na(stage1_problem_detail)]),
      # Random-effects summaries: first parameterization
      mean_stage1_intercept_variance = if ("stage1_intercept_variance" %in% names(results)) {
        safe_mean(stage1_intercept_variance)
      } else {
        NULL
      },
      mean_stage1_slope_variance = if ("stage1_slope_variance" %in% names(results)) {
        safe_mean(stage1_slope_variance)
      } else {
        NULL
      },
      mean_stage1_intercept_slope_covariance = if ("stage1_intercept_slope_covariance" %in% names(results)) {
        safe_mean(stage1_intercept_slope_covariance)
      } else {
        NULL
      },
      mean_stage1_intercept_slope_correlation = if ("stage1_intercept_slope_correlation" %in% names(results)) {
        safe_mean(stage1_intercept_slope_correlation)
      } else {
        NULL
      },

      # Random-effects summaries: y/q parameterization
      mean_stage1_y_intercept_variance = if ("stage1_y_intercept_variance" %in% names(results)) {
        safe_mean(stage1_y_intercept_variance)
      } else {
        NULL
      },
      mean_stage1_y_slope_variance = if ("stage1_y_slope_variance" %in% names(results)) {
        safe_mean(stage1_y_slope_variance)
      } else {
        NULL
      },
      mean_stage1_y_intercept_slope_covariance = if ("stage1_y_intercept_slope_covariance" %in% names(results)) {
        safe_mean(stage1_y_intercept_slope_covariance)
      } else {
        NULL
      },
      mean_stage1_y_intercept_slope_correlation = if ("stage1_y_intercept_slope_correlation" %in% names(results)) {
        safe_mean(stage1_y_intercept_slope_correlation)
      } else {
        NULL
      },
      mean_stage1_q_intercept_variance = if ("stage1_q_intercept_variance" %in% names(results)) {
        safe_mean(stage1_q_intercept_variance)
      } else {
        NULL
      },
      mean_stage1_q_slope_variance = if ("stage1_q_slope_variance" %in% names(results)) {
        safe_mean(stage1_q_slope_variance)
      } else {
        NULL
      },
      mean_stage1_q_intercept_slope_covariance = if ("stage1_q_intercept_slope_covariance" %in% names(results)) {
        safe_mean(stage1_q_intercept_slope_covariance)
      } else {
        NULL
      },
      mean_stage1_q_intercept_slope_correlation = if ("stage1_q_intercept_slope_correlation" %in% names(results)) {
        safe_mean(stage1_q_intercept_slope_correlation)
      } else {
        NULL
      },
      .groups = "drop"
    )
  results
}

#' Summarize OpenMx and estimator issue classes.
#'
#' @details
#' Estimator wrappers classify OpenMx or standard-error problems into compact
#' `mx_issue_class` labels. This helper keeps only non-`ok` issue classes and
#' reports their frequency within each condition-study-method group, along with
#' a compact example detail string for inspection.
#'
#' @param results Replication-level result tibble.
#'
#' @return A tibble with issue-class counts and proportions, or an empty tibble
#'   when `mx_issue_class` is absent.
summarize_issue_df <- function(results) {
  if (!("mx_issue_class" %in% names(results))) {
    return(tibble::tibble())
  }

  design_cols <- intersect(
    c(
      "condition_id", "study", "study_version", "calibration_version",
      "simulation_module", "sensitivity_block", "sensitivity_block_label",
      "variance_ratio_arm", "variance_ratio_pair_id",
      "variance_ratio_pair_label", "sensitivity_reference_ratio",
      "sensitivity_is_low_reliability_stress",
      "method", "method_role", "fuller_variant",
      "fuller_preliminary_moment", "fuller_variance_bread",
      "fuller_predictor_outcome_covariance_source",
      "mx_issue_class", "num_clus",
      "mean_clus_size", "target_reliability", "achieved_reliability",
      "marginal_rho", "standardized_beta_target", "structural_target",
      "calibration_arm", "calibration_metric", "calibration_target",
      "calibration_target_value", "achieved_calibration_value",
      "covariance_shape_fixed", "posterior_reliability_anchor", "icc_anchor",
      "target_calibration_reliability", "achieved_calibration_reliability",
      "achieved_partial_reliability", "slope_variance_marginal",
      "residualized_slope_variance", "slope_intercept_variance_ratio",
      "G_condition_number", "intercept_icc", "marginal_slope_icc",
      "conditional_slope_icc",
      "structural_r2", "focal_unique_r2", "balance_mode", "r_structure",
      "r_rho", "sigma", "mean_clus_size_y", "mean_clus_size_q",
      "target_reliability_y", "target_reliability_q",
      "achieved_reliability_y", "achieved_reliability_q", "sigma_y", "sigma_q",
      "information_profile", "is_falsification_control", "information_matched",
      "profile_min_clus_size", "profile_max_clus_size",
      "reliability_sd", "reliability_min", "reliability_max"
    ),
    names(results)
  )

  results %>%
    dplyr::group_by(condition_id, study, method) %>%
    dplyr::mutate(n_rep_total = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::filter(
      !is.na(mx_issue_class),
      nzchar(mx_issue_class),
      mx_issue_class != "ok"
    ) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(design_cols))) %>%
    dplyr::summarise(
      n_rep = dplyr::n(),
      n_rep_total = dplyr::first(n_rep_total),
      prop_rep = dplyr::n() / dplyr::first(n_rep_total),
      n_distinct_issue_details = dplyr::n_distinct(mx_issue_detail[!is.na(mx_issue_detail)]),
      sample_issue_detail = compact_message(mx_issue_detail[!is.na(mx_issue_detail)]),
      .groups = "drop"
    )
}

#' Append one row to a progress CSV.
#'
#' @details
#' The progress file is intentionally append-only so interrupted jobs retain a
#' record of completed and skipped conditions. The first write includes headers;
#' later writes append rows without column names.
#'
#' @param progress_path Path to the run-level progress CSV.
#' @param row_df One-row data frame or tibble with progress metadata.
#'
#' @return Invisibly returns the result of the underlying write call.
write_progress_row <- function(progress_path, row_df) {
  if (!file.exists(progress_path)) {
    utils::write.csv(row_df, file = progress_path, row.names = FALSE)
  } else {
    utils::write.table(
      row_df,
      file = progress_path,
      sep = ",",
      row.names = FALSE,
      col.names = FALSE,
      append = TRUE
    )
  }
}

#' Collect a complete, one-row provenance record for a simulation invocation.
#'
#' The run manifest deliberately separates invocation-level provenance from the
#' condition manifest and replication rows. This makes it possible to identify
#' the exact software state that produced a chunk without copying the same
#' package metadata into every estimator row.
collect_vh_run_provenance <- function(
    design, n_sim, study_arg, chunk_meta, n_cores,
    repo_root_path = if (exists("repo_root", inherits = TRUE)) {
      get("repo_root", inherits = TRUE)
    } else {
      getwd()
    },
    started_at = Sys.time()) {
  git_output <- function(args) {
    tryCatch(
      system2(
        "git",
        c("-C", shQuote(normalizePath(repo_root_path)), args),
        stdout = TRUE,
        stderr = FALSE
      ),
      error = function(e) character()
    )
  }
  git_commit <- git_output(c("rev-parse", "HEAD"))
  git_status <- git_output(c("status", "--porcelain", "--untracked-files=normal"))

  package_names <- c("OpenMx", "MplusAutomation", "lme4", "nlme", "geigen")
  package_versions <- vapply(package_names, function(package_name) {
    if (!requireNamespace(package_name, quietly = TRUE)) {
      return(NA_character_)
    }
    as.character(utils::packageVersion(package_name))
  }, character(1L))

  mplus_executable <- unname(Sys.which("mplus"))
  mplus_version <- Sys.getenv("MPLUS_VERSION", unset = NA_character_)
  if (!nzchar(mplus_executable)) {
    mplus_executable <- NA_character_
  }
  if (is.na(mplus_version) && !is.na(mplus_executable)) {
    executable_strings <- tryCatch(
      system2(
        "strings",
        shQuote(mplus_executable),
        stdout = TRUE,
        stderr = FALSE
      ),
      error = function(e) character()
    )
    version_lines <- grep(
      "^Mplus VERSION ", executable_strings,
      value = TRUE, ignore.case = TRUE
    )
    if (length(version_lines) > 0L) {
      mplus_version <- trimws(version_lines[[1L]])
    }
  }

  out <- tibble::tibble(
    metadata_schema_version = "vh_run_provenance_v1",
    pipeline_version = vh_pipeline_version(),
    git_commit = if (length(git_commit) > 0L) git_commit[[1L]] else NA_character_,
    git_dirty = length(git_status) > 0L,
    git_dirty_entry_count = as.integer(length(git_status)),
    requested_n_sim = as.integer(n_sim),
    study_selector = as.character(study_arg[[1L]]),
    chunk_index = as.integer(chunk_meta$chunk_index),
    chunk_size = as.integer(chunk_meta$chunk_size),
    condition_id_start = as.integer(chunk_meta$condition_start),
    condition_id_end = as.integer(chunk_meta$condition_end),
    selected_condition_count = as.integer(nrow(design)),
    requested_cores = as.integer(n_cores),
    r_version = R.version.string,
    mplus_executable = mplus_executable,
    mplus_version = mplus_version,
    slurm_job_id = Sys.getenv("SLURM_JOB_ID", unset = NA_character_),
    slurm_array_job_id = Sys.getenv("SLURM_ARRAY_JOB_ID", unset = NA_character_),
    slurm_array_task_id = Sys.getenv("SLURM_ARRAY_TASK_ID", unset = NA_character_),
    run_started_at = format(started_at, "%Y-%m-%dT%H:%M:%S%z"),
    hostname = Sys.info()[["nodename"]]
  )
  for (package_name in names(package_versions)) {
    out[[paste0("package_", tolower(package_name), "_version")]] <-
      package_versions[[package_name]]
  }
  out
}

#' Atomically write a small CSV artifact.
#'
#' Atomic replacement prevents concurrent array tasks from exposing a partial
#' deterministic crosswalk or provenance file to downstream readers.
write_vh_csv_atomic <- function(data, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary_path <- tempfile(
    pattern = paste0(".", basename(path), "-"),
    tmpdir = dirname(path),
    fileext = ".tmp"
  )
  on.exit(unlink(temporary_path, force = TRUE), add = TRUE)
  utils::write.csv(data, temporary_path, row.names = FALSE)
  if (!file.rename(temporary_path, path)) {
    stop("Could not atomically write CSV artifact: ", path)
  }
  invisible(path)
}

#' Save the versioned deterministic ICC/reliability crosswalk.
save_vh_icc_crosswalk <- function(out_dir) {
  crosswalk <- make_icc_posterior_reliability_crosswalk() %>%
    dplyr::mutate(
      crosswalk_version = "icc_posterior_reliability_crosswalk_v1",
      pipeline_version = vh_pipeline_version(),
      .before = 1L
    )
  path <- file.path(
    out_dir,
    "icc_posterior_reliability_crosswalk_v1.csv"
  )
  write_vh_csv_atomic(crosswalk, path)
  invisible(path)
}

#' Read and type-normalize a replication result file.
#'
#' @details
#' `data.table::fread()` is fast, but resumed aggregate builds can encounter
#' plain CSV files from older runs as well as current `.csv.gz` files. This
#' helper reads either form and coerces known numeric, integer, and logical
#' columns back to their expected types before row-binding across conditions.
#'
#' @param path Path to a per-condition replication result CSV or CSV.GZ file.
#'
#' @return A tibble containing the replication-level rows from `path`.
read_replication_results_file <- function(path) {
  out <- tibble::as_tibble(data.table::fread(path))
  fuller_logical_cols <- c(
    "fuller_auto_tempered", "fuller_auto_guard_pass",
    "fuller_auto_full_weight_guard_pass", "fuller_auto_search_nonmonotone"
  )
  fuller_integer_cols <- "fuller_auto_search_evaluations"
  fuller_character_cols <- c(
    "fuller_auto_guard_reason", "fuller_auto_full_weight_guard_reason",
    "fuller_variant", "fuller_preliminary_moment", "fuller_variance_bread",
    "fuller_predictor_outcome_covariance_source"
  )
  fuller_numeric_cols <- setdiff(
    grep("^fuller_", names(out), value = TRUE),
    c(fuller_logical_cols, fuller_integer_cols, fuller_character_cols)
  )
  expanded_logical_cols <- grep(
    "^(stage1(_[yq])?|mx(_[a-z]+)?|mplus(_[a-z]+)?)_.*boundary$",
    names(out),
    value = TRUE
  )
  expanded_integer_cols <- intersect(
    c(
      "replication_seed", "bridge_pair_id", "variance_ratio_pair_id",
      "simulation_seed_group",
      "mplus_warning_count", "mplus_target_parameter_count",
      "mplus_boundary_variance_count", "mplus_nonpositive_variance_count"
    ),
    names(out)
  )
  expanded_character_cols <- intersect(
    c(
      "bridge_pair_label", "simulation_module", "sensitivity_block",
      "sensitivity_block_label", "variance_ratio_arm",
      "variance_ratio_pair_label"
    ),
    names(out)
  )
  expanded_numeric_cols <- setdiff(
    grep(
      paste0(
        "^stage1(_[yq])?_(lambda|theta|posterior_variance|",
        "corrected_score_error_covariance|fitted_|true_blup|blup_slope|",
        "corrected_score_slope|latent_covariance)|",
        "^mx_(raw_focal|latent_|predictor_latent_|outcome_latent_)|",
        "^mplus_"
      ),
      names(out),
      value = TRUE
    ),
    c(
      expanded_logical_cols, expanded_integer_cols,
      "mplus_critical_warning", "mplus_critical_warning_detail"
    )
  )

  # Keep these casts centralized so aggregate rebuilds do not depend on the
  # exact type inference chosen by fread for any single condition file.
  numeric_cols <- unique(c(intersect(
    # TODO: update these columns
    c(
      "estimate", "se", "ci_low", "ci_high", "truth",
      "mx_condition_number", "mplus_fitted_latent_slope_sd",
      "stage1_re_corr", "stage1_eb_corr", "stage1_design_kappa",
      "stage1_y_re_corr", "stage1_y_eb_corr", "stage1_y_design_kappa",
      "stage1_q_re_corr", "stage1_q_eb_corr", "stage1_q_design_kappa",
      "target_reliability", "achieved_reliability", "marginal_rho",
      "calibration_target_value", "achieved_calibration_value",
      "posterior_reliability_anchor", "icc_anchor",
      "target_calibration_reliability", "achieved_calibration_reliability",
      "achieved_partial_reliability", "slope_variance_marginal",
      "residualized_slope_variance", "slope_intercept_variance_ratio",
      "sensitivity_reference_ratio",
      "G_condition_number", "intercept_icc", "marginal_slope_icc",
      "conditional_slope_icc", "calibration_tau0", "calibration_tau0_sq",
      "calibration_tau1_sq", "sigma",
      "standardized_beta_target", "structural_r2", "focal_unique_r2",
      "tau1", "beta1z", "beta2z", "outcome_residual_variance",
      "target_reliability_y", "target_reliability_q",
      "achieved_reliability_y", "achieved_reliability_q",
      "achieved_partial_reliability_y", "achieved_partial_reliability_q",
      "intercept_icc_y", "intercept_icc_q",
      "marginal_slope_icc_y", "marginal_slope_icc_q",
      "conditional_slope_icc_y", "conditional_slope_icc_q",
      "tau1_y", "tau1_q", "theta0", "theta1",
      "standardized_theta0", "slope_variance_marginal_y",
      "slope_variance_marginal_q", "slope_variance_residual_q",
      "tau1_residual_q", "rho_residual_q", "sigma_y", "sigma_q",
      "structural_residual_G_min_eigen", "structural_joint_G_min_eigen",
      "sigma2", "var_u1", "sigma_z", "fuller_lambda1", "fuller_lambda2",
      "fuller_sigma2", "fuller_weight_min", "fuller_weight_max",
      "fuller_correction_c", "profile_small_weight", "profile_large_weight",
      "reliability_sd", "reliability_iqr", "reliability_min", "reliability_max",
      "reliability_small", "reliability_large",
      "partial_reliability_mean", "partial_reliability_min",
      "partial_reliability_max", "partial_reliability_small",
      "partial_reliability_large",
      "population_lambda22_mean", "population_lambda22_min",
      "population_lambda22_max",
      "population_lambda_matrix_frobenius_rms_dispersion",
      "population_theta22_mean", "population_theta22_min",
      "population_theta22_max",
      "population_theta_matrix_frobenius_rms_dispersion",
      "population_ols_var22_mean", "population_ols_var22_min",
      "population_ols_var22_max",
      "population_ols_cov_matrix_frobenius_rms_dispersion",
      "realized_reliability_mean", "realized_reliability_sd",
      "realized_reliability_iqr", "realized_reliability_min",
      "realized_reliability_max",
      "mean_realized_trials", "min_realized_trials", "max_realized_trials",
      "prop_ids_leq_3_trials",
      "lambda22_mean", "lambda22_sd", "lambda22_min", "lambda22_max",
      "lambda_matrix_frobenius_rms_dispersion",
      "theta22_mean", "theta22_sd", "theta22_min", "theta22_max",
      "theta_matrix_frobenius_rms_dispersion",
      "ols_var22_mean", "ols_var22_sd", "ols_var22_min", "ols_var22_max",
      "ols_cov_matrix_frobenius_rms_dispersion",
      "blup_slope_bias_small", "blup_slope_rmse_small",
      "blup_slope_bias_large", "blup_slope_rmse_large",
      "corrected_slope_bias_small", "corrected_slope_rmse_small",
      "corrected_slope_bias_large", "corrected_slope_rmse_large",
      "average_measurement_slope_bias_small",
      "average_measurement_slope_rmse_small",
      "average_measurement_slope_bias_large",
      "average_measurement_slope_rmse_large",
      "stage1_intercept_variance",
      "stage1_slope_variance",
      "stage1_y_intercept_variance",
      "stage1_y_slope_variance",
      "stage1_q_intercept_variance",
      "stage1_q_slope_variance",
      "stage1_intercept_slope_covariance",
      "stage1_intercept_slope_correlation",
      "stage1_y_intercept_slope_covariance",
      "stage1_y_intercept_slope_correlation",
      "stage1_q_intercept_slope_covariance",
      "stage1_q_intercept_slope_correlation"
    ),
    names(out)
  ), fuller_numeric_cols, expanded_numeric_cols))
  integer_cols <- intersect(
    c(
      "condition_id", "rep", "num_clus", "mean_clus_size",
      "mean_clus_size_y", "mean_clus_size_q", "status_code",
      "profile_min_clus_size", "profile_max_clus_size",
      "profile_small_clus_size", "profile_large_clus_size",
      "score_small_cluster_size", "score_large_cluster_size",
      "mplus_warning_count", "mplus_target_parameter_count",
      fuller_integer_cols, expanded_integer_cols
    ),
    names(out)
  )
  logical_cols <- intersect(
    c(
      "stage1_singular_problem", "stage1_lmer_singular",
      "stage1_y_singular_problem", "stage1_y_lmer_singular",
      "stage1_q_singular_problem", "stage1_q_lmer_singular",
      "is_falsification_control", "information_matched",
      "sensitivity_is_low_reliability_stress",
      "covariance_shape_fixed",
      "mx_info_definite", "mplus_critical_warning",
      "point_eligible", "interval_eligible", "analysis_eligible",
      fuller_logical_cols, expanded_logical_cols
    ),
    names(out)
  )
  character_cols <- intersect(
    c(
      "point_exclusion_reason", "interval_exclusion_reason",
      "analysis_exclusion_reason", "mx_issue_class", "mx_issue_detail",
      "mx_status_msg", "mplus_critical_warning_detail",
      "stage1_problem_detail", fuller_character_cols,
      expanded_character_cols
    ),
    names(out)
  )

  if (length(numeric_cols) > 0L) {
    out <- out %>% dplyr::mutate(dplyr::across(dplyr::all_of(numeric_cols), ~ suppressWarnings(as.numeric(.x))))
  }
  if (length(integer_cols) > 0L) {
    out <- out %>% dplyr::mutate(dplyr::across(dplyr::all_of(integer_cols), ~ suppressWarnings(as.integer(.x))))
  }
  if (length(logical_cols) > 0L) {
    out <- out %>% dplyr::mutate(dplyr::across(dplyr::all_of(logical_cols), as.logical))
  }
  if (length(character_cols) > 0L) {
    out <- out %>% dplyr::mutate(dplyr::across(dplyr::all_of(character_cols), as.character))
  }

  out
}

#' Load all completed condition-level replication results for a design.
#'
#' @details
#' This is used after the condition loop to rebuild aggregate outputs from the
#' per-condition artifacts on disk. Rebuilding from disk makes `resume_existing`
#' runs deterministic: skipped conditions and newly completed conditions are
#' combined through the same code path.
#'
#' @param design Design rows that should be represented in the aggregate output.
#'   Must contain `condition_id`.
#' @param out_dir Root output directory containing the `conditions` subdirectory.
#'
#' @return A tibble formed by row-binding all discovered condition replication
#'   files, or an empty tibble if none are present.
load_completed_condition_results <- function(design, out_dir) {
  files <- purrr::map_chr(design$condition_id, function(condition_id) {
    gz_path <- get_condition_file_paths(out_dir, condition_id)$replications
    csv_path <- sub("\\.gz$", "", gz_path)
    if (file.exists(gz_path)) {
      gz_path
    } else if (file.exists(csv_path)) {
      csv_path
    } else {
      NA_character_
    }
  })
  existing <- stats::na.omit(files)
  if (length(existing) == 0L) {
    return(tibble::tibble())
  }

  purrr::map_dfr(existing, read_replication_results_file)
}

load_completed_condition_csv <- function(design, out_dir, artifact) {
  if (!(artifact %in% c("summary", "issue_summary", "stage1_summary"))) {
    stop("Unsupported condition artifact: ", artifact)
  }
  files <- purrr::map_chr(design$condition_id, function(condition_id) {
    path <- get_condition_file_paths(out_dir, condition_id)[[artifact]]
    if (file.exists(path)) path else NA_character_
  })
  existing <- stats::na.omit(files)
  if (length(existing) == 0L) {
    return(tibble::tibble())
  }

  purrr::map_dfr(existing, function(path) {
    tryCatch(
      tibble::as_tibble(utils::read.csv(path, check.names = FALSE)),
      error = function(e) tibble::tibble()
    )
  })
}

make_condition_replication_index <- function(design, out_dir) {
  purrr::map_dfr(design$condition_id, function(condition_id) {
    path <- get_condition_file_paths(out_dir, condition_id)$replications
    tibble::tibble(
      condition_id = condition_id,
      replication_file = path,
      file_exists = file.exists(path),
      file_size_bytes = if (file.exists(path)) file.info(path)$size else NA_real_
    )
  })
}

estimate_replication_result_rows <- function(design, n_sim) {
  method_counts <- vapply(
    seq_len(nrow(design)),
    function(i) length(study_methods_for_condition(design[i, , drop = FALSE])),
    integer(1)
  )
  as.double(n_sim) * sum(method_counts)
}

condition_output_is_complete <- function(
    path,
    expected_methods,
    expected_n_sim,
    expected_pipeline_version = vh_pipeline_version()) {
  if (!file.exists(path)) {
    return(FALSE)
  }
  out <- tryCatch(read_replication_results_file(path), error = function(e) NULL)
  required <- c("method", "rep", "pipeline_version")
  if (is.null(out) || !all(required %in% names(out))) {
    return(FALSE)
  }
  if (!all(out$pipeline_version == expected_pipeline_version)) {
    return(FALSE)
  }
  # A row-complete file can still represent an interrupted Mplus executable or
  # filesystem failure. Such conditions must be rerun instead of being accepted
  # by resume solely because every method/replication row was written.
  if ("mx_issue_class" %in% names(out) &&
      any(out$mx_issue_class == "mplus_null_fit", na.rm = TRUE)) {
    return(FALSE)
  }
  if (!setequal(unique(out$method), expected_methods) ||
      !setequal(unique(out$rep), seq_len(expected_n_sim))) {
    return(FALSE)
  }
  method_counts <- table(factor(out$method, levels = expected_methods))
  all(method_counts == expected_n_sim)
}

#' Dispatch one replication to the correct study module.
#'
#' @param condition One-row condition tibble with a `study` column equal to
#'  a legacy, amended-v2, Study 5, or ICC-bridge study key. Covariance-shape
#'  sensitivity rows retain the amended-v2 key of their source DGM.
#'
#' @return A replication-level result tibble from the corresponding
#'   `run_study*_rep()` function.
run_study_rep <- function(condition) {
  study_key <- as.character(condition$study[[1]])
  switch(
    study_key,
    study0 = run_study1_rep(condition),
    study1 = run_study1_rep(condition),
    study2 = run_study2_rep(condition),
    study3 = run_study3_rep(condition),
    study4 = run_study4_rep(condition),
    study5 = run_study5_rep(condition),
    study1v2 = run_study1_rep(condition),
    study2v2 = run_study2_rep(condition),
    study3v2 = run_study3_rep(condition),
    study4v2 = run_study4_rep(condition),
    iccbridge = run_study2_rep(condition),
    stop("Unsupported VH study key: ", study_key)
  )
}

#' List objects that must be exported to parallel workers.
#'
#' @details
#' `foreach` workers run in separate R processes and therefore need explicit
#' access to locally defined functions that are not package exports. This list
#' includes study simulators, estimator helpers, diagnostics, and OpenMx
#' utilities used anywhere below `run_study_rep()`.
#'
#' @return A character vector of object names passed to `foreach(..., .export)`.
vig_hallquist_parallel_exports <- function() {
  candidates <- c(
    "run_study_rep", "run_study1_rep", "run_study2_rep", "run_study3_rep", "run_study4_rep", "run_study5_rep",
    "vh_pipeline_version", "vh_replication_seed",
    "add_vh_analysis_eligibility",
    "add_study1_analysis_eligibility",
    "add_study2_method_roles", "add_study2_analysis_eligibility",
    "add_study3_analysis_eligibility",
    "add_study4_method_roles", "add_study4_analysis_eligibility",
    "add_stage1_estimates", "summarize_stage1_measurement_diagnostics",
    "rescale_fuller_to_population_sd",
    "add_zero_fuller_predictor_outcome_covariance",
    "prepare_fuller_average_measurement", "fit_fuller_average_measurement",
    "simulate_study1", "simulate_study2", "simulate_study3", "simulate_study4", "simulate_study5",
    "simulate_data_blup_as_outcome", "simulate_data_blup_as_predictor",
    "simulate_data_dual_blup", "simulate_data_study4",
    "make_failed_result", "add_study_result_context",
    "combine_dual_stage1_diagnostics",
    "fit_stage1", "condition_uses_non_iid_R", "condition_to_r_spec",
    "condition_to_nlme_correlation", "draw_level1_residuals",
    "balance_mode_to_sim_arg", "make_reliability_time_design",
    "study4_profile_spec", "study4_time_design", "make_study4_cluster_sizes",
    "study4_weighted_quantile", "study4_matrix_rms_dispersion",
    "study4_measurement_matrix_summary", "study4_measurement_diagnostics",
    "draw_random_effects",
    "matched_study_methods", "disparate_study_methods", "study_methods_for_condition", "tempered_eiv_methods",
    "condition_includes_tempered_eiv", "lai_condition_to_r_spec", "lai_condition_uses_non_iid_R",
    "lai_condition_to_nlme_correlation", "add_lai_trial_index", "draw_lai_level1_residuals",
    "fit_lai_stage1", "extract_centered_slope_eb", "fit_tempered_eiv_dual_set", "lai_truth", "make_covu",
    "make_study2_cluster_sizes", "make_study3_cluster_sizes", "fixed_params",
    "study1_methods", "study2_methods", "study3_methods", "study4_methods", "study5_methods",
    "safe_lmer", "safe_lme", "empty_stage1_diagnostics",
    "assess_stage1_random_effect_covariance", "get_stage1_diagnostics",
    "normalize_r_spec", "make_R_matrix", "draw_residuals_from_R",
    "make_random_effect_covariance", "posterior_random_effect_covariance",
    "get_closed_form_corrected_scores", "get_stage1_eb_components",
    "extract_stage1_components", "extract_stage1_components.merMod", "extract_stage1_components.lme",
    "extract_stage1_components.default", "normalize_R_list", "as_plain_vcov_matrix", "stage1_fixef",
    "format_stage1_eb_row", "default_re_design", "make_eb_output_row", "select_lai_measurement_columns",
    "openmx_bivariate_loading_columns",
    "compute_eb_measurement_inputs", "compute_bivariate_eb_inputs",
    "compute_univariate_eb_inputs", "compute_lai_2spa_inputs",
    "assess_dual_ols_design", "fit_observed_single", "fit_observed_dual",
    "finalize_ols_se_variants", "fit_eiv_dual", "fit_ridge_dual", "fit_fuller",
    "fuller_dual_result_columns", "fuller_matrix_diagnostics",
    "fuller_relative_min_eigen", "fuller_row_max_finite",
    "fuller_guard_penalty", "score_fuller_auto_candidates",
    "fuller_reference_dual_se", "fit_fuller_dual_core", "fit_fuller_dual",
    "fit_fuller_dual_variants",
    "fit_fuller_dual_stepdown", "fit_fuller_dual_alpha_stepdown",
    "fit_lai_2spa", "fit_lai_2spa_observed_outcome", "fit_lai_2spa_disparate",
    "fit_lai_2spa_dual_process",
    "run_mx_safe", "extract_mx_stats", "extract_mx_se_details",
    "mx_numeric_scalar", "mx_diagnostics_tibble",
    "mx_latent_covariance_diagnostics",
    "classify_mx_issue", "compact_message", "project_to_pd",
    "run_mplus_modeler_writable_tmp", "fit_mplus_blup_predictor",
    "extract_mplus_stats", "mplus_diagnostics_template", "mplus_message_lines",
    "populate_mplus_latent_diagnostics",
    "mplus_warning_diagnostics"
  )
  candidates[vapply(candidates, exists, logical(1), mode = "function", inherits = TRUE) |
    candidates == "fixed_params"]
}

#' Run all replications for one design condition.
#'
#' @details
#' Each replication receives a deterministic seed derived from a fixed base
#' seed, a condition seed group, and the replication index. The seed group is
#' normally the condition identifier; paired ICC-bridge and variance-ratio
#' sensitivity arms instead share a `simulation_seed_group` so their Monte
#' Carlo draws are paired.
#'
#' Parallel execution uses `foreach` and expects a registered backend, which is
#' set up by `run_simulation()` when `n_cores > 1`.
#'
#' @param condition One-row condition tibble.
#' @param n_sim Number of replications to run for the condition.
#' @param n_cores Number of cores requested. Values greater than one select the
#'   `foreach` path.
#'
#' @return A tibble with one row per estimator/method result per replication,
#'   augmented with condition descriptors and the replication id.
run_condition_replications <- function(condition, n_sim, n_cores = 1L) {
  rep_ids <- seq_len(n_sim)

  run_single_rep <- function(rep_id) {
    replication_seed <- vh_replication_seed(condition, rep_id)
    set.seed(replication_seed)
    rep_out <- run_study_rep(condition)
    dplyr::bind_cols(
      rep_out,
      condition[rep(1L, nrow(rep_out)), , drop = FALSE] %>% dplyr::select(-study),
      tibble::tibble(
        rep = rep_id,
        replication_seed = replication_seed,
        pipeline_version = vh_pipeline_version()
      )
    )
  }

  if (n_cores > 1L) {
    foreach::foreach(
      rep_id = rep_ids,
      .combine = dplyr::bind_rows,
      # TODO: check these exports
      .packages = c("data.table", "lme4", "MASS", "dplyr", "tidyr", "purrr", "tibble", "OpenMx",
                    "glmnet", "sandwich", "geigen", "MplusAutomation", "glue"),
      .export = vig_hallquist_parallel_exports()
    ) %dopar% {
      run_single_rep(rep_id)
    }
  } else {
    purrr::map_dfr(rep_ids, run_single_rep)
  }
}

#' Run a Vig-Hallquist simulation batch.
#'
#' @details
#' This is the top-level entry point used by the replication script. It
#' selects the requested study design, optionally slices the design into a
#' chunk, runs each condition, writes condition-level artifacts immediately, and
#' finally rebuilds aggregate replication and summary files from the completed
#' condition outputs.
#'
#' The per-condition write pattern is deliberate: long simulation jobs can be
#' interrupted and later resumed with `resume_existing = TRUE`. A condition is
#' skipped only when both its replication file and main summary file already
#' exist. After the loop, aggregate outputs are rebuilt from all condition files
#' found for the selected design/chunk.
#'
#' Output files include:
#' - A manifest CSV with the selected design rows.
#' - A run-metadata CSV with code, software, scheduler, and host provenance.
#' - For ICC-bridge runs, a standalone versioned deterministic crosswalk CSV.
#' - An append-only progress CSV.
#' - Per-condition replication, summary, issue-summary, and first-stage
#'   problem-summary files under `conditions/`.
#' - Chunk-specific aggregate files. Full-selection runs also write legacy
#'   aggregate filenames without the chunk label.
#'
#' @param n_sim Positive integer number of replications per condition.
#' @param study_arg Study selector passed to `select_design()`, typically
#'   `"all"` for legacy Studies 1--5, `"allv2"` for amended Studies
#'   1--4, an individual legacy/v2 key, `"iccbridge"`, or the compact
#'   `"ratiosensitivity"` module.
#' @param out_dir Root output directory. Defaults to
#'   `file.path(vig_hallquist_dir, "outputs", "vig_hallquist")`.
#' @param n_cores Positive integer number of worker cores. Values greater than
#'   one register a `doParallel` backend and use the parallel replication path.
#' @param max_conditions Optional cap passed to `select_design()` after study
#'   selection.
#' @param chunk_index Optional one-based chunk index for distributed runs. Must
#'   be supplied with `chunk_size`.
#' @param chunk_size Optional number of selected design rows per chunk. Must be
#'   supplied with `chunk_index`.
#' @param resume_existing Logical. If `TRUE`, skip conditions whose
#'   per-condition replication and summary files already exist.
#' @param max_aggregate_replication_rows Maximum expected rows that may be
#'   loaded to materialize a combined replication table. Larger runs aggregate
#'   condition summaries and write a replication-file index instead.
#'
#' @return Invisibly returns a list with aggregate `results`, `summary`,
#'   `issue_summary`, and `stage1_summary` tibbles.
run_simulation <- function(n_sim = 100L,
                           study_arg = "all",
                           out_dir = file.path(vig_hallquist_dir, "outputs", "vig_hallquist"),
                           n_cores = 1L,
                           max_conditions = NA_integer_,
                           chunk_index = NA_integer_,
                           chunk_size = NA_integer_,
                           resume_existing = TRUE,
                           max_aggregate_replication_rows = 2e6) {
  run_started_at <- Sys.time()
  if (is.na(n_sim) || n_sim < 1L) {
    stop("`n_sim` must be a positive integer.")
  }
  if (is.na(n_cores) || n_cores < 1L) {
    stop("`n_cores` must be a positive integer.")
  }
  if (is.na(max_aggregate_replication_rows) ||
      max_aggregate_replication_rows <= 0) {
    stop("`max_aggregate_replication_rows` must be positive or `Inf`.")
  }

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  design <- select_design(study_arg = study_arg, max_conditions = max_conditions)
  design <- slice_design_chunk(design, chunk_index = chunk_index, chunk_size = chunk_size)
  if (nrow(design) == 0L) {
    message(
      "No conditions in requested chunk ", chunk_index,
      " for study selection `", study_arg, "`; exiting successfully."
    )
    return(invisible(list(
      results = tibble::tibble(),
      summary = tibble::tibble(),
      issue_summary = tibble::tibble(),
      stage1_summary = tibble::tibble()
    )))
  }
  chunk_meta <- attr(design, "chunk_meta")
  chunk_label <- make_chunk_label(chunk_meta)
  selection_label <- gsub(
    "[^a-z0-9]+",
    "_",
    tolower(as.character(study_arg[[1]]))
  )
  file_prefix <- sprintf(
    "vig_hallquist_%s_%s",
    selection_label,
    chunk_label
  )

  # Chunk labels keep concurrently run jobs from overwriting each other's
  # aggregate files while preserving condition-level paths shared by resume.
  manifest_path <- file.path(out_dir, sprintf("%s_manifest.csv", file_prefix))
  run_metadata_path <- file.path(
    out_dir,
    sprintf("%s_run_metadata.csv", file_prefix)
  )
  progress_path <- file.path(out_dir, sprintf("%s_progress.csv", file_prefix))
  aggregate_replications_path <- file.path(out_dir, sprintf("%s_replication_results.csv.gz", file_prefix))
  aggregate_summary_path <- file.path(out_dir, sprintf("%s_summary.csv", file_prefix))
  aggregate_issue_summary_path <- file.path(out_dir, sprintf("%s_issue_summary.csv", file_prefix))
  aggregate_stage1_summary_path <- file.path(out_dir, sprintf("%s_stage1_problem_summary.csv", file_prefix))
  replication_index_path <- file.path(
    out_dir,
    sprintf("%s_condition_replication_index.csv", file_prefix)
  )

  utils::write.csv(design, file = manifest_path, row.names = FALSE)
  run_metadata <- collect_vh_run_provenance(
    design = design,
    n_sim = n_sim,
    study_arg = study_arg,
    chunk_meta = chunk_meta,
    n_cores = n_cores,
    started_at = run_started_at
  )
  write_vh_csv_atomic(run_metadata, run_metadata_path)
  if (any(design$study == "iccbridge")) {
    save_vh_icc_crosswalk(out_dir)
  }

  if (n_cores > 1L) {
    doParallel::registerDoParallel(cores = n_cores)
    on.exit(doParallel::stopImplicitCluster(), add = TRUE)
  }

  for (i in seq_len(nrow(design))) {
    condition <- design[i, , drop = FALSE]
    condition_id <- condition$condition_id[[1]]
    paths <- get_condition_file_paths(out_dir, condition_id)
    started_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

    # A condition is resumable only when both the detailed rows and primary
    # summary exist; partial artifacts are overwritten by rerunning the condition.
    if (isTRUE(resume_existing) && file.exists(paths$replications) && file.exists(paths$summary)) {
      expected_methods <- study_methods_for_condition(condition)
      if (condition_output_is_complete(
          paths$replications,
          expected_methods = expected_methods,
          expected_n_sim = n_sim)) {
        write_progress_row(
          progress_path,
          tibble::tibble(
            condition_id = condition_id,
            chunk_label = chunk_label,
            status = "skipped_existing",
            started_at = started_at,
            finished_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
            elapsed_seconds = 0,
            n_methods = NA_integer_,
            n_success = NA_integer_,
            replications_file = paths$replications,
            summary_file = paths$summary
          )
        )
        next
      }
    }

    timing <- system.time({
      condition_results <- run_condition_replications(condition, n_sim = n_sim, n_cores = n_cores)
    })
    condition_summary <- summarize_results_df(condition_results)
    condition_issue_summary <- summarize_issue_df(condition_results)
    condition_stage1_summary <- summarize_stage1_problem_df(condition_results)

    data.table::fwrite(condition_results, file = paths$replications)
    utils::write.csv(condition_summary, file = paths$summary, row.names = FALSE)
    utils::write.csv(condition_issue_summary, file = paths$issue_summary, row.names = FALSE)
    utils::write.csv(condition_stage1_summary, file = paths$stage1_summary, row.names = FALSE)

    write_progress_row(
      progress_path,
      tibble::tibble(
        condition_id = condition_id,
        chunk_label = chunk_label,
        status = "completed",
        started_at = started_at,
        finished_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
        elapsed_seconds = unname(timing[["elapsed"]]),
        n_methods = nrow(condition_results),
        n_success = sum(!is.na(condition_results$estimate)),
        replications_file = paths$replications,
        summary_file = paths$summary
      )
    )
  }

  expected_result_rows <- estimate_replication_result_rows(design, n_sim)
  materialize_replications <- is.infinite(max_aggregate_replication_rows) ||
    expected_result_rows <= max_aggregate_replication_rows

  if (materialize_replications) {
    results <- load_completed_condition_results(design, out_dir)
    summary_df <- if (nrow(results) > 0L) summarize_results_df(results) else tibble::tibble()
    issue_summary_df <- if (nrow(results) > 0L) summarize_issue_df(results) else tibble::tibble()
    stage1_summary_df <- if (nrow(results) > 0L) summarize_stage1_problem_df(results) else tibble::tibble()
    data.table::fwrite(results, file = aggregate_replications_path)
  } else {
    results <- tibble::tibble()
    summary_df <- load_completed_condition_csv(design, out_dir, "summary")
    issue_summary_df <- load_completed_condition_csv(design, out_dir, "issue_summary")
    stage1_summary_df <- load_completed_condition_csv(design, out_dir, "stage1_summary")
    replication_index <- make_condition_replication_index(design, out_dir)
    utils::write.csv(replication_index, replication_index_path, row.names = FALSE)
    if (file.exists(aggregate_replications_path)) {
      unlink(aggregate_replications_path)
    }
    message(
      "Skipped in-memory replication aggregation for approximately ",
      format(expected_result_rows, big.mark = ",", scientific = FALSE),
      " rows. Condition-level replication files remain authoritative."
    )
  }

  utils::write.csv(summary_df, file = aggregate_summary_path, row.names = FALSE)
  utils::write.csv(issue_summary_df, file = aggregate_issue_summary_path, row.names = FALSE)
  utils::write.csv(stage1_summary_df, file = aggregate_stage1_summary_path, row.names = FALSE)
  # Preserve the original aggregate filenames for downstream scripts that do
  # not know about chunk-specific output naming.
  if (identical(chunk_label, "full_selection") &&
      identical(tolower(as.character(study_arg[[1]])), "all")) {
    if (materialize_replications) {
      data.table::fwrite(results, file = file.path(out_dir, "vig_hallquist_replication_results.csv.gz"))
    } else {
      legacy_replication_path <- file.path(
        out_dir,
        "vig_hallquist_replication_results.csv.gz"
      )
      if (file.exists(legacy_replication_path)) {
        unlink(legacy_replication_path)
      }
      utils::write.csv(
        replication_index,
        file = file.path(out_dir, "vig_hallquist_condition_replication_index.csv"),
        row.names = FALSE
      )
    }
    utils::write.csv(summary_df, file = file.path(out_dir, "vig_hallquist_summary.csv"), row.names = FALSE)
    utils::write.csv(issue_summary_df, file = file.path(out_dir, "vig_hallquist_issue_summary.csv"), row.names = FALSE)
    utils::write.csv(stage1_summary_df, file = file.path(out_dir, "vig_hallquist_stage1_problem_summary.csv"), row.names = FALSE)
  }

  message("Saved outputs to: ", normalizePath(out_dir))
  message("Run provenance: ", run_metadata_path)
  if (materialize_replications) {
    message("Aggregate replication results: ", aggregate_replications_path)
  } else {
    message("Condition replication index: ", replication_index_path)
  }
  message("Aggregate summary results: ", aggregate_summary_path)
  message("Aggregate issue summary results: ", aggregate_issue_summary_path)
  message("Aggregate stage1 problem summary results: ", aggregate_stage1_summary_path)

  invisible(list(results = results, summary = summary_df, issue_summary = issue_summary_df, stage1_summary = stage1_summary_df))
}
