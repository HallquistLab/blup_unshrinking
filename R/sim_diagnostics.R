#' Diagnostics for first-stage mixed-model fits used in stage-2 score analyses.

#' Return an empty stage-1 diagnostics template.
#'
#' @details
#' When upstream model fitting fails, this helper returns a row filled with
#' missing values. This keeps downstream calls to `dplyr::bind_rows()` stable
#' by ensuring the expected columns are always present.
#'
#' @return A `tibble` with a single row of NA values for diagnostic columns.
empty_stage1_diagnostics <- function() {
  tibble::tibble(
    stage1_singular_problem = NA,
    stage1_problem_detail = NA_character_,
    stage1_lmer_singular = NA,
    stage1_re_corr = NA_real_,
    stage1_eb_corr = NA_real_,
    stage1_design_kappa = NA_real_
  )
}


#' Apply the common VH point- and interval-eligibility contract.
#'
#' @details
#' Simulation performance has two distinct denominators. A replication can
#' provide a usable point estimate while failing to provide a defensible
#' standard error or confidence interval. This helper applies the same contract
#' to VH Studies 1--4:
#'
#' - `point_eligible` requires a successful status, a finite estimate, and any
#'   method-specific design/convergence checks supplied by the estimator.
#' - `interval_eligible` additionally requires a finite positive standard error
#'   and finite, ordered confidence limits.
#'
#' Existing point- or interval-eligibility fields produced by an estimator are
#' consumed as inputs, including legacy `analysis_eligible` and
#' `analysis_exclusion_reason` fields from dual-OLS design checks. Thus, calling
#' the classifier again can add exclusions but cannot reinstate a result that an
#' estimator explicitly rejected. On output, the legacy `analysis_*` names are
#' retained as aliases for point eligibility so older analysis scripts keep a
#' stable point-estimate meaning. New code should use the explicit point and
#' interval names.
#'
#' Lai 2S-PA rows are screened for classified OpenMx problems, non-definite
#' information matrices, and extreme reported condition numbers. MSEM rows are
#' screened for critical Mplus warnings and ambiguous target-parameter matches.
#' First-stage singularity and Fuller matrix diagnostics remain reportable
#' diagnostics rather than automatic exclusions unless they cause estimation to
#' return a nonzero status or nonfinite result.
#'
#' @param results Replication-level estimator rows.
#' @param mx_condition_limit Finite upper bound for an available OpenMx
#'   condition-number diagnostic.
#'
#' @return `results` augmented with explicit point/interval eligibility fields
#'   and backward-compatible `analysis_*` aliases.
add_vh_analysis_eligibility <- function(results, mx_condition_limit = 1e12) {
  if (!is.finite(mx_condition_limit) || mx_condition_limit <= 0) {
    stop("`mx_condition_limit` must be finite and positive.")
  }

  get_column <- function(name, default) {
    if (name %in% names(results)) results[[name]] else rep(default, nrow(results))
  }
  set_reason <- function(reason, when, value) {
    index <- is.na(reason) & !is.na(when) & when
    if (length(value) == 1L) {
      reason[index] <- value
    } else {
      reason[index] <- value[index]
    }
    reason
  }

  method <- as.character(get_column("method", NA_character_))
  status_code <- suppressWarnings(as.integer(get_column("status_code", NA_integer_)))
  estimate <- suppressWarnings(as.numeric(get_column("estimate", NA_real_)))
  se <- suppressWarnings(as.numeric(get_column("se", NA_real_)))
  ci_low <- suppressWarnings(as.numeric(get_column("ci_low", NA_real_)))
  ci_high <- suppressWarnings(as.numeric(get_column("ci_high", NA_real_)))

  # Preserve method-specific eligibility emitted upstream (currently the
  # complete-case dual-OLS rank/VIF check) before replacing the legacy alias.
  upstream_point_eligible <- if ("point_eligible" %in% names(results)) {
    as.logical(results$point_eligible)
  } else {
    as.logical(get_column("analysis_eligible", NA))
  }
  upstream_point_reason <- if ("point_exclusion_reason" %in% names(results)) {
    as.character(results$point_exclusion_reason)
  } else {
    as.character(get_column("analysis_exclusion_reason", NA_character_))
  }
  upstream_interval_eligible <- as.logical(
    get_column("interval_eligible", NA)
  )
  upstream_interval_reason <- as.character(
    get_column("interval_exclusion_reason", NA_character_)
  )

  mx_issue_class <- as.character(get_column("mx_issue_class", NA_character_))
  mx_info_definite <- as.logical(get_column("mx_info_definite", NA))
  mx_condition_number <- suppressWarnings(as.numeric(
    get_column("mx_condition_number", NA_real_)
  ))
  mplus_critical_warning <- as.logical(
    get_column("mplus_critical_warning", NA)
  )
  mplus_target_parameter_count <- suppressWarnings(as.integer(
    get_column("mplus_target_parameter_count", NA_integer_)
  ))

  point_reason <- rep(NA_character_, nrow(results))
  point_reason <- set_reason(
    point_reason,
    is.na(status_code) | is.na(estimate),
    "estimation_unavailable"
  )
  point_reason <- set_reason(
    point_reason,
    !is.na(status_code) & status_code != 0L,
    "estimation_status_nonzero"
  )
  point_reason <- set_reason(
    point_reason,
    !is.na(estimate) & !is.finite(estimate),
    "nonfinite_estimate"
  )
  point_reason <- set_reason(
    point_reason,
    !is.na(upstream_point_eligible) & !upstream_point_eligible,
    ifelse(
      is.na(upstream_point_reason) | upstream_point_reason == "",
      "stage2_design_ineligible",
      upstream_point_reason
    )
  )

  openmx_method <- method %in% c("lai_2spa", "lai_2spaa")
  point_reason <- set_reason(
    point_reason,
    openmx_method & !is.na(mx_issue_class) & mx_issue_class != "ok",
    "openmx_issue"
  )
  point_reason <- set_reason(
    point_reason,
    openmx_method & !is.na(mx_info_definite) & !mx_info_definite,
    "openmx_information_not_definite"
  )
  point_reason <- set_reason(
    point_reason,
    openmx_method & is.finite(mx_condition_number) &
      mx_condition_number > mx_condition_limit,
    "openmx_condition_number_excessive"
  )

  mplus_method <- method %in% c("msem", "sem")
  point_reason <- set_reason(
    point_reason,
    mplus_method & !is.na(mplus_critical_warning) & mplus_critical_warning,
    "mplus_critical_warning"
  )
  point_reason <- set_reason(
    point_reason,
    mplus_method & !is.na(mplus_target_parameter_count) &
      mplus_target_parameter_count != 1L,
    "mplus_target_parameter_not_unique"
  )

  interval_reason <- point_reason
  interval_reason <- set_reason(
    interval_reason,
    !is.na(upstream_interval_eligible) & !upstream_interval_eligible,
    ifelse(
      is.na(upstream_interval_reason) | upstream_interval_reason == "",
      "estimator_interval_ineligible",
      upstream_interval_reason
    )
  )
  interval_reason <- set_reason(
    interval_reason,
    !is.finite(se) | se <= 0,
    "invalid_standard_error"
  )
  interval_reason <- set_reason(
    interval_reason,
    !is.finite(ci_low) | !is.finite(ci_high) | ci_low > ci_high,
    "invalid_confidence_interval"
  )

  dplyr::mutate(
    results,
    point_eligible = is.na(point_reason),
    point_exclusion_reason = point_reason,
    interval_eligible = is.na(interval_reason),
    interval_exclusion_reason = interval_reason,
    # Backward-compatible aliases. Their meaning is now explicitly the
    # point-estimate denominator in every VH study.
    analysis_eligible = point_eligible,
    analysis_exclusion_reason = point_exclusion_reason
  )
}


#' Assess a fitted two-dimensional random-effect covariance matrix.
#'
#' A failed or boundary mixed-model fit can contain `NA`, non-finite, or
#' non-positive marginal variance estimates. Calling `stats::cov2cor()` on
#' that matrix emits a warning and does not produce an interpretable
#' correlation. Treat those matrices as singular, retain `NA` for the
#' correlation, and let the replication-level diagnostic record the event.
#'
#' @param vcov_re Candidate intercept/slope covariance matrix.
#' @param model_singular Optional singularity flag supplied by the mixed-model
#'   implementation (for example, `lme4::isSingular()`).
#' @param tolerance Eigenvalue tolerance used to identify boundary fits.
#'
#' @return A list containing `correlation`, `invalid`, and `singular`.
assess_stage1_random_effect_covariance <- function(
    vcov_re, model_singular = FALSE, tolerance = 1e-5) {
  matrix_valid <- is.matrix(vcov_re) &&
    nrow(vcov_re) >= 2L && ncol(vcov_re) >= 2L
  covariance <- if (matrix_valid) {
    as.matrix(vcov_re[1:2, 1:2, drop = FALSE])
  } else {
    matrix(NA_real_, nrow = 2L, ncol = 2L)
  }
  covariance <- (covariance + t(covariance)) / 2
  finite_entries <- all(is.finite(covariance))
  positive_variances <- finite_entries && all(diag(covariance) > 0)
  invalid <- !matrix_valid || !finite_entries || !positive_variances

  correlation <- if (!invalid) {
    covariance[1L, 2L] /
      sqrt(covariance[1L, 1L] * covariance[2L, 2L])
  } else {
    NA_real_
  }
  min_eigenvalue <- if (!invalid) {
    tryCatch(
      min(eigen(covariance, symmetric = TRUE, only.values = TRUE)$values),
      error = function(e) NA_real_
    )
  } else {
    NA_real_
  }
  singular <- invalid || isTRUE(model_singular) ||
    !is.finite(min_eigenvalue) || min_eigenvalue < tolerance

  list(
    correlation = correlation,
    invalid = invalid,
    singular = singular
  )
}


#' Extract diagnostics from a first-stage mixed-model fit.
#'
#' @details
#' This function extracts random-effect correlations and checks for boundary
#' or near-singular estimates that could cause numerical issues in stage-2
#' score estimators. For `lme4` fits, it uses `lme4::isSingular`; for
#' `nlme::lme` fits, it checks whether the fitted random-effect covariance has
#' a near-zero eigenvalue. It also flags nearly perfect fitted random-effect or
#' empirical Bayes correlations.
#'
#' @param fit_obj Fitted first-stage model from `lme4` or `nlme`. Can be `NULL`
#' if the fit failed.
#' @param stage2_df A data frame containing the cluster-level EB scores or
#' predictors used in the stage-2 model.
#' @param predictor_u0 The column name in `stage2_df` for the intercept-like score.
#' @param predictor_u1 The column name in `stage2_df` for the slope-like score.
#'
#' @return A `tibble` with 1 row containing summary diagnostics for the fit,
#' including singularity flags and correlation bounds.
get_stage1_diagnostics <- function(fit_obj, stage2_df, predictor_u0 = "u0_eb", predictor_u1 = "u1_eb") {
  if (is.null(fit_obj)) {
    return(empty_stage1_diagnostics())
  }

  vcov_re <- tryCatch({
    if (inherits(fit_obj, "lme")) {
      as_plain_vcov_matrix(nlme::getVarCov(fit_obj, type = "random.effects"))
    } else {
      as.matrix(lme4::VarCorr(fit_obj)[[1]])
    }
  }, error = function(e) matrix(NA_real_, nrow = 2L, ncol = 2L))
  model_singular <- if (inherits(fit_obj, "lme")) {
    FALSE
  } else {
    tryCatch(lme4::isSingular(fit_obj, tol = 1e-5), error = function(e) TRUE)
  }
  covariance_assessment <- assess_stage1_random_effect_covariance(
    vcov_re,
    model_singular = model_singular,
    tolerance = 1e-5
  )
  re_corr <- covariance_assessment$correlation
  lmer_singular <- covariance_assessment$singular

  # Filter to complete cases to ensure the empirical Bayes correlation and
  # design kappa calculations match the exact data available for stage-2 analysis.
  dat <- stage2_df[, c(predictor_u0, predictor_u1), drop = FALSE]
  dat <- dat[stats::complete.cases(dat), , drop = FALSE]

  eb_corr <- if (nrow(dat) >= 2L &&
    is.finite(stats::sd(dat[[predictor_u0]])) &&
    is.finite(stats::sd(dat[[predictor_u1]])) &&
    stats::sd(dat[[predictor_u0]]) > sqrt(.Machine$double.eps) &&
    stats::sd(dat[[predictor_u1]]) > sqrt(.Machine$double.eps)) {
    stats::cor(dat[[predictor_u0]], dat[[predictor_u1]])
  } else {
    NA_real_
  }

  # Compute the condition number (kappa) of the design matrix. A kappa > 1e6
  # suggests extreme multicollinearity that can destabilize stage-2 estimators.
  design_kappa <- if (nrow(dat) >= 3L) {
    tryCatch({
      k <- kappa(cbind(1, as.matrix(dat)))
      if (k > 1/.Machine$double.eps) { Inf } else { k }
    }, error = function(e) NA_real_
    )
  } else {
    NA_real_
  }

  # Check for near-perfect correlations which often indicate structural
  # non-identifiability or zero variance components in the simulation.
  boundary_re_corr <- is.finite(re_corr) && abs(re_corr) > 0.999
  boundary_eb_corr <- is.finite(eb_corr) && abs(eb_corr) > 0.999
  singular_problem <- isTRUE(lmer_singular) || boundary_re_corr || boundary_eb_corr

  detail_parts <- c(
    if (isTRUE(covariance_assessment$invalid)) {
      "stage1_random_effect_cov_invalid"
    } else if (isTRUE(lmer_singular)) {
      "stage1_random_effect_cov_singular"
    } else {
      NA_character_
    },
    if (boundary_re_corr) sprintf("random_effect_corr=%0.6f", re_corr) else NA_character_,
    if (boundary_eb_corr) sprintf("eb_corr=%0.6f", eb_corr) else NA_character_,
    if (is.finite(design_kappa) && design_kappa > 1e6) sprintf("design_kappa=%0.3e", design_kappa) else NA_character_
  )

  tibble::tibble(
    stage1_singular_problem = singular_problem,
    stage1_problem_detail = if (isTRUE(singular_problem)) compact_message(detail_parts) else "ok",
    stage1_lmer_singular = lmer_singular,
    stage1_re_corr = re_corr,
    stage1_eb_corr = eb_corr,
    stage1_design_kappa = design_kappa
  )
}
