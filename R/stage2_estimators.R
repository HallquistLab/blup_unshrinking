#' Stage-2 estimators for observed EB/corrected-score regressions.
#'
#' The functions in this file fit second-stage regressions where the predictor
#' of substantive interest is a cluster-level random slope estimate. To keep
#' simulation summaries directly comparable across estimator variants, the
#' reported slope is scaled by the empirical standard deviation of the slope
#' predictor. A one-unit estimate therefore corresponds to a one-SD change in
#' the stage-2 slope predictor.

#' Construct the empty OLS result template.
#'
#' @details
#' OLS-based estimators return both model-based (`naive`) and HC3 robust
#' standard-error variants. This helper centralizes the NA-filled return shape
#' used when the model cannot be estimated, keeping downstream row-binding
#' stable.
#'
#' @return
#' A tibble with one row for each `se_type` (`naive`, `hc3`) and columns
#' `estimate`, `se`, `ci_low`, `ci_high`, and `status_code`, all filled with
#' missing values.
empty_lm_variants <- function(analysis_diagnostics = NULL) {
  if (is.null(analysis_diagnostics)) {
    analysis_diagnostics <- tibble::tibble(
      analysis_eligible = FALSE,
      analysis_exclusion_reason = "estimation_unavailable",
      stage2_predictor_correlation = NA_real_,
      stage2_vif = NA_real_,
      stage2_design_rank = NA_integer_,
      stage2_design_rcond = NA_real_
    )
  }
  dplyr::bind_rows(
    naive = dplyr::bind_cols(tibble::tibble(
      estimate = NA_real_,
      se = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      status_code = NA_integer_
    ), analysis_diagnostics),
    hc3 = dplyr::bind_cols(tibble::tibble(
      estimate = NA_real_,
      se = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      status_code = NA_integer_
    ), analysis_diagnostics),
    .id = "se_type"
  )
}

#' Assess whether a dual-predictor OLS design is suitable for performance summaries.
#'
#' @details
#' A finite OLS coefficient can still be uninformative when its two predictors
#' are nearly collinear. This helper preserves such estimates for auditability,
#' but marks them ineligible for primary bias, variance, RMSE, and coverage
#' summaries. The diagnostic is calculated on the exact complete-case Stage-2
#' data used by the OLS fit and is invariant to the predictors' units.
#'
#' A VIF of 100 is the default near-collinearity threshold. With two
#' predictors, this is equivalent to an absolute predictor correlation of about
#' 0.995. Exact rank deficiency is always ineligible regardless of the VIF.
#'
#' @param dat Complete-case Stage-2 data containing the two predictors.
#' @param predictor_u0 Character name of the intercept-like predictor.
#' @param predictor_u1 Character name of the slope-like predictor.
#' @param vif_threshold Positive VIF threshold used to flag near collinearity.
#'
#' @return A one-row tibble with analysis-eligibility and Stage-2 design
#' diagnostics.
assess_dual_ols_design <- function(
    dat,
    predictor_u0,
    predictor_u1,
    vif_threshold = 100) {
  unavailable <- function(reason) {
    tibble::tibble(
      analysis_eligible = FALSE,
      analysis_exclusion_reason = reason,
      stage2_predictor_correlation = NA_real_,
      stage2_vif = NA_real_,
      stage2_design_rank = NA_integer_,
      stage2_design_rcond = NA_real_
    )
  }

  if (!is.finite(vif_threshold) || vif_threshold <= 1) {
    stop("`vif_threshold` must be finite and greater than one.")
  }
  if (nrow(dat) < 4L) {
    return(unavailable("insufficient_complete_cases"))
  }

  x <- as.matrix(dat[, c(predictor_u0, predictor_u1), drop = FALSE])
  if (any(!is.finite(x))) {
    return(unavailable("nonfinite_stage2_predictors"))
  }
  x_sd <- apply(x, 2L, stats::sd)
  if (any(!is.finite(x_sd)) || any(x_sd <= sqrt(.Machine$double.eps))) {
    return(unavailable("degenerate_stage2_predictor"))
  }

  # Standardization makes the rank, reciprocal condition number, and VIF
  # comparable across the very different intercept and slope scales.
  x_standardized <- scale(x, center = TRUE, scale = TRUE)
  x_design <- cbind(1, x_standardized)
  rank <- qr(x_design)$rank
  correlation <- stats::cor(x_standardized[, 1L], x_standardized[, 2L])
  vif <- if (is.finite(correlation) && abs(correlation) < 1) {
    1 / (1 - correlation^2)
  } else {
    Inf
  }
  design_rcond <- tryCatch(
    rcond(crossprod(x_design)),
    error = function(e) NA_real_
  )

  exclusion_reason <- if (rank < ncol(x_design)) {
    "stage2_rank_deficient"
  } else if (!is.finite(correlation) || !is.finite(vif) || !is.finite(design_rcond)) {
    "stage2_nonfinite_design_diagnostic"
  } else if (vif >= vif_threshold) {
    "stage2_near_collinear"
  } else {
    NA_character_
  }

  tibble::tibble(
    analysis_eligible = is.na(exclusion_reason),
    analysis_exclusion_reason = exclusion_reason,
    stage2_predictor_correlation = correlation,
    stage2_vif = vif,
    stage2_design_rank = as.integer(rank),
    stage2_design_rcond = design_rcond
  )
}

#' Fit a single-predictor observed-score stage-2 OLS regression.
#'
#' @details
#' Fits `outcome ~ predictor` after complete-case filtering. The coefficient
#' for `predictor` and its standard errors are multiplied by the observed
#' standard deviation of `predictor`, so the reported estimand is the expected
#' change in `outcome` per one observed SD of the predictor.
#'
#' Both conventional model-based and HC3 robust standard errors are returned.
#' Failures are represented by the `empty_lm_variants()` template rather than
#' by thrown errors, because these functions are used inside simulation loops.
#'
#' @param stage2_df Data frame containing the outcome and stage-2 predictor.
#' @param outcome Character scalar naming the outcome column.
#' @param predictor Character scalar naming the single predictor column.
#' @param reporting_scale Optional positive scale applied to the raw
#'   coefficient and standard error. Defaults to the observed predictor SD.
#'
#' @return
#' A two-row tibble keyed by `se_type`, with scaled `estimate`, `se`,
#' Wald-normal confidence limits, and `status_code` (`0L` for successfully
#' estimated rows).
fit_observed_single <- function(stage2_df, outcome, predictor, reporting_scale = NULL) {
  empty_out <- empty_lm_variants()

  # Restrict to the variables used by this estimator before applying the
  # complete-case rule; missingness in unrelated simulation columns is ignored.
  dat <- stage2_df[, c(outcome, predictor), drop = FALSE]
  dat <- dat[stats::complete.cases(dat), , drop = FALSE]

  # Need at least three rows for an intercept, slope, and residual degree of
  # freedom, and a non-degenerate predictor to make the scaled slope meaningful.
  if (nrow(dat) < 3L || !is.finite(stats::sd(dat[[predictor]])) ||
    stats::sd(dat[[predictor]]) <= sqrt(.Machine$double.eps)) {
    return(empty_out)
  }

  fit <- tryCatch(stats::lm(stats::as.formula(paste(outcome, "~", predictor)), data = dat), error = function(e) NULL)
  if (is.null(fit)) {
    return(empty_out)
  }

  coef_tab <- summary(fit)$coefficients
  if (!(predictor %in% rownames(coef_tab))) {
    return(empty_out)
  }

  # By default report the effect per observed predictor SD. Simulation callers
  # can supply a fixed population scale so all proxy methods target the same
  # standardized latent effect.
  scale_u1 <- if (is.null(reporting_scale)) {
    stats::sd(dat[[predictor]])
  } else {
    as.numeric(reporting_scale[[1]])
  }
  if (!is.finite(scale_u1) || scale_u1 <= 0) {
    return(empty_out)
  }
  est <- unname(coef_tab[predictor, "Estimate"]) * scale_u1
  se_naive <- unname(coef_tab[predictor, "Std. Error"]) * scale_u1

  # HC3 is the leverage-adjusted sandwich variance used for the robust OLS
  # variant. If sandwich fails, only the HC3 row is marked unavailable.
  vcov_hc3 <- tryCatch(sandwich::vcovHC(fit, type = "HC3"), error = function(e) NULL)
  se_hc3 <- if (!is.null(vcov_hc3) && predictor %in% rownames(vcov_hc3)) {
    sqrt(unname(vcov_hc3[predictor, predictor])) * scale_u1
  } else {
    NA_real_
  }

  dplyr::bind_rows(
    naive = tibble::tibble(
      estimate = est,
      se = se_naive,
      ci_low = est - stats::qnorm(0.975) * se_naive,
      ci_high = est + stats::qnorm(0.975) * se_naive,
      status_code = 0L
    ),
    hc3 = tibble::tibble(
      estimate = est,
      se = se_hc3,
      ci_low = est - stats::qnorm(0.975) * se_hc3,
      ci_high = est + stats::qnorm(0.975) * se_hc3,
      status_code = ifelse(is.finite(se_hc3), 0L, NA_integer_)
    ),
    .id = "se_type"
  )
}

#' Fit a dual-predictor observed-score stage-2 OLS regression.
#'
#' @details
#' Fits `outcome ~ predictor_u0 + predictor_u1` after complete-case filtering
#' and reports the coefficient for `predictor_u1`. The intercept-like random
#' effect proxy (`predictor_u0`) is included as an adjustment variable, while
#' the slope-like random effect proxy (`predictor_u1`) is the target estimand.
#'
#' The returned estimate and standard errors are scaled by the observed
#' standard deviation of `predictor_u1`. This mirrors `fit_observed_single()`
#' and supports direct comparisons across simulation conditions.
#'
#' @param stage2_df Data frame containing the outcome and both predictors.
#' @param outcome Character scalar naming the outcome column.
#' @param predictor_u0 Character scalar naming the intercept-like predictor.
#' @param predictor_u1 Character scalar naming the slope-like predictor whose
#' coefficient is reported.
#' @param reporting_scale Optional positive scale applied to the raw slope
#'   coefficient and standard error. Defaults to the observed slope-proxy SD.
#'
#' @return
#' A two-row tibble keyed by `se_type`, with scaled `estimate`, `se`,
#' Wald-normal confidence limits, and `status_code` (`0L` for successfully
#' estimated rows).
fit_observed_dual <- function(
    stage2_df,
    outcome,
    predictor_u0,
    predictor_u1,
    reporting_scale = NULL) {
  # Complete-case filtering is estimator-specific: rows only need to be
  # observed on the outcome and the two predictors used here.
  dat <- stage2_df[, c(outcome, predictor_u0, predictor_u1), drop = FALSE]
  dat <- dat[stats::complete.cases(dat), , drop = FALSE]
  analysis_diagnostics <- assess_dual_ols_design(
    dat = dat,
    predictor_u0 = predictor_u0,
    predictor_u1 = predictor_u1
  )
  empty_out <- empty_lm_variants(analysis_diagnostics)

  # Four rows are the minimum for a two-predictor model with residual degrees
  # of freedom. Degenerate slope predictors cannot be converted to one-SD units.
  if (nrow(dat) < 4L || !is.finite(stats::sd(dat[[predictor_u1]])) ||
    stats::sd(dat[[predictor_u1]]) <= sqrt(.Machine$double.eps)) {
    return(empty_out)
  }

  fit <- tryCatch(
    stats::lm(stats::as.formula(paste(outcome, "~", predictor_u0, "+", predictor_u1)), data = dat),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(empty_out)
  }

  coef_tab <- summary(fit)$coefficients
  if (!(predictor_u1 %in% rownames(coef_tab))) {
    return(empty_out)
  }

  # Only the slope-like predictor is reported; the u0 proxy is an adjustment
  # covariate used to reduce confounding between random intercepts and slopes.
  scale_u1 <- if (is.null(reporting_scale)) {
    stats::sd(dat[[predictor_u1]])
  } else {
    as.numeric(reporting_scale[[1]])
  }
  if (!is.finite(scale_u1) || scale_u1 <= 0) {
    return(empty_out)
  }
  est <- unname(coef_tab[predictor_u1, "Estimate"]) * scale_u1
  se_naive <- unname(coef_tab[predictor_u1, "Std. Error"]) * scale_u1

  # Compute the HC3 standard error for the target coefficient on the same
  # one-SD scale as the conventional OLS standard error.
  vcov_hc3 <- tryCatch(sandwich::vcovHC(fit, type = "HC3"), error = function(e) NULL)
  se_hc3 <- if (!is.null(vcov_hc3) && predictor_u1 %in% rownames(vcov_hc3)) {
    sqrt(unname(vcov_hc3[predictor_u1, predictor_u1])) * scale_u1
  } else {
    NA_real_
  }

  dplyr::bind_rows(
    naive = dplyr::bind_cols(tibble::tibble(
      estimate = est,
      se = se_naive,
      ci_low = est - stats::qnorm(0.975) * se_naive,
      ci_high = est + stats::qnorm(0.975) * se_naive,
      status_code = 0L
    ), analysis_diagnostics),
    hc3 = dplyr::bind_cols(tibble::tibble(
      estimate = est,
      se = se_hc3,
      ci_low = est - stats::qnorm(0.975) * se_hc3,
      ci_high = est + stats::qnorm(0.975) * se_hc3,
      status_code = ifelse(is.finite(se_hc3), 0L, NA_integer_)
    ), analysis_diagnostics),
    .id = "se_type"
  )
}

#' Convert OLS standard-error variants to method-labelled estimator rows.
#'
#' @details
#' `fit_observed_single()` and `fit_observed_dual()` return rows keyed by
#' `se_type` so the model fit and the robust variance calculation can stay
#' together. Simulation summaries instead group by `method`, so this helper
#' maps the `naive` row to `base_method` and the `hc3` row to
#' `paste0(base_method, "_hc3")`.
#'
#' @param fit_tbl Tibble returned by an OLS stage-2 fitting helper.
#' @param base_method Character scalar naming the non-HC3 method.
#'
#' @return
#' A tibble with `method`, `estimate`, `se`, `ci_low`, `ci_high`, and
#' `status_code`.
finalize_ols_se_variants <- function(fit_tbl, base_method) {
  fit_tbl <- dplyr::mutate(
    fit_tbl,
    method = ifelse(se_type == "hc3", paste0(base_method, "_hc3"), base_method)
  )
  dplyr::select(
    fit_tbl,
    method, estimate, se, ci_low, ci_high, status_code,
    dplyr::any_of(c(
      "analysis_eligible", "analysis_exclusion_reason",
      "stage2_predictor_correlation", "stage2_vif",
      "stage2_design_rank", "stage2_design_rcond"
    ))
  )
}

#' Fit a ridge-regularized dual-predictor stage-2 regression.
#'
#' @details
#' Fits a Gaussian ridge regression with `glmnet::cv.glmnet()` using
#' `predictor_u0` and `predictor_u1`, then reports the cross-validated
#' coefficient for `predictor_u1` scaled by the observed SD of `predictor_u1`.
#' Ridge is used as a point-estimation fallback for collinear stage-2 predictors;
#' no analytic standard error or confidence interval is computed here.
#'
#' The number of folds is capped at 10 and floored at 3, with the sample size
#' used when fewer than 10 complete cases are available.
#'
#' @param stage2_df Data frame containing the outcome and both predictors.
#' @param outcome Character scalar naming the outcome column.
#' @param predictor_u0 Character scalar naming the intercept-like predictor.
#' @param predictor_u1 Character scalar naming the slope-like predictor whose
#' coefficient is reported.
#' @param lambda_rule Character scalar naming the `cv.glmnet` lambda selector,
#' typically `"lambda.1se"` or `"lambda.min"`.
#'
#' @return
#' A one-row tibble with scaled `estimate`, missing `se` and confidence limits,
#' `status_code`, and the selected `ridge_lambda` plus `ridge_lambda_rule`.
fit_ridge_dual <- function(stage2_df, outcome, predictor_u0, predictor_u1, lambda_rule = "lambda.1se") {
  # Ridge uses the same complete-case and non-degeneracy checks as the dual OLS
  # estimator so comparisons are driven by estimator behavior, not row sets.
  dat <- stage2_df[, c(outcome, predictor_u0, predictor_u1), drop = FALSE]
  dat <- dat[stats::complete.cases(dat), , drop = FALSE]

  if (nrow(dat) < 4L || !is.finite(stats::sd(dat[[predictor_u1]])) ||
    stats::sd(dat[[predictor_u1]]) <= sqrt(.Machine$double.eps)) {
    return(tibble::tibble(
      estimate = NA_real_, se = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
      status_code = NA_integer_, ridge_lambda = NA_real_, ridge_lambda_rule = lambda_rule
    ))
  }

  x_mat <- as.matrix(dat[, c(predictor_u0, predictor_u1), drop = FALSE])
  y_vec <- dat[[outcome]]

  # Cross-validation needs at least three folds; with small simulation samples,
  # use leave-one-ish folds up to the usual 10-fold cap.
  nfolds <- max(3L, min(10L, nrow(dat)))
  cv_fit <- tryCatch(
    glmnet::cv.glmnet(
      x = x_mat,
      y = y_vec,
      family = "gaussian",
      alpha = 0,
      nfolds = nfolds,
      standardize = TRUE,
      intercept = TRUE,
      type.measure = "mse"
    ),
    error = function(e) NULL
  )

  if (is.null(cv_fit)) {
    return(tibble::tibble(
      estimate = NA_real_, se = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
      status_code = NA_integer_, ridge_lambda = NA_real_, ridge_lambda_rule = lambda_rule
    ))
  }

  lambda_value <- cv_fit[[lambda_rule]]

  # Extracting by the lambda rule rather than by the numeric lambda preserves
  # glmnet's handling of named selectors such as "lambda.1se".
  coef_mat <- tryCatch(as.matrix(stats::coef(cv_fit, s = lambda_rule)), error = function(e) NULL)
  if (is.null(coef_mat) || !(predictor_u1 %in% rownames(coef_mat))) {
    return(tibble::tibble(
      estimate = NA_real_, se = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
      status_code = NA_integer_, ridge_lambda = lambda_value, ridge_lambda_rule = lambda_rule
    ))
  }

  scale_u1 <- stats::sd(dat[[predictor_u1]])
  tibble::tibble(
    estimate = unname(coef_mat[predictor_u1, 1]) * scale_u1,
    se = NA_real_,
    ci_low = NA_real_,
    ci_high = NA_real_,
    status_code = 0L,
    ridge_lambda = lambda_value,
    ridge_lambda_rule = lambda_rule
  )
}

#' Columns returned by EIV estimators for simulation output.
#'
#' @return Character vector of EIV result columns to preserve in replication
#'   outputs.
eiv_result_columns <- function() {
  c(
    "method", "estimate", "se", "ci_low", "ci_high", "status_code",
    "se_type",
    "mx_issue_class", "mx_issue_detail",
    "eiv_measurement_weight_requested", "eiv_measurement_weight_used",
    "eiv_regularized", "eiv_latent_cov_min_eigen",
    "eiv_latent_cov_condition_number"
  )
}

#' Select the stable EIV result schema after assigning a method label.
#'
#' @param result One-row EIV result tibble with a `method` column.
#'
#' @return `result` restricted to the EIV simulation-output columns available
#'   in the input.
select_eiv_result_columns <- function(result) {
  dplyr::select(result, dplyr::any_of(eiv_result_columns()))
}

#' Convert EIV standard-error variants to method-labelled estimator rows.
#'
#' @details
#' `fit_eiv_dual()` returns rows keyed by `se_type`. This helper mirrors
#' `finalize_ols_se_variants()`: the conventional model-based row keeps
#' `base_method`, while HC0 and HC3 rows append `"_hc0"` and `"_hc3"`.
#'
#' @param fit_tbl Tibble returned by `fit_eiv_dual()`.
#' @param base_method Character scalar naming the model-based EIV method.
#'
#' @return
#' A tibble with method labels and the stable EIV result columns.
finalize_eiv_se_variants <- function(fit_tbl, base_method) {
  fit_tbl <- dplyr::mutate(
    fit_tbl,
    method = dplyr::case_when(
      .data$se_type == "naive" ~ base_method,
      TRUE ~ paste0(base_method, "_", .data$se_type)
    )
  )
  select_eiv_result_columns(fit_tbl)
}

#' Fit a dual-predictor errors-in-variables corrected-score estimator.
#'
#' @details
#' This estimator treats the observed stage-2 predictors as noisy measurements
#' of latent random-effect scores. For each complete case it subtracts the
#' supplied measurement-error covariance matrix from `x_i x_i'` before solving
#' the corrected normal equations. The reported target is the `predictor_u1`
#' coefficient multiplied by an estimated latent SD for `predictor_u1`.
#'
#' The estimating equation uses predictors ordered as
#' `(intercept, predictor_u0, predictor_u1)`. The measurement-error variances
#' and covariance (`meas11`, `meas12`, `meas22`) fill the predictor block only;
#' the intercept is treated as measured without error. `measurement_weight`
#' scales the supplied covariance terms. Values below 1 define a tempered EIV
#' sensitivity path: they deliberately under-subtract sampling-error covariance
#' to regularize the deconvolution, so they should not be interpreted as the
#' classical full EIV correction unless there is a separate calibration reason.
#'
#' A sandwich variance is computed from the empirical estimating-function
#' residuals for the HC0 and HC3 rows. The `naive` row uses a conventional
#' homoskedastic model-based variance for the corrected normal equations. If
#' `stabilize_a_mat` is `TRUE`, the full corrected cross-product matrix is
#' projected to positive definite before solving. If `ridge_predictor_block` is
#' `TRUE`, only the two-predictor block is ridged up to `ridge_min_eigen`.
#'
#' @param stage2_df Data frame containing outcome, predictors, and measurement
#' error covariance columns.
#' @param outcome Character scalar naming the outcome column.
#' @param predictor_u0 Character scalar naming the intercept-like observed
#' predictor.
#' @param predictor_u1 Character scalar naming the slope-like observed
#' predictor whose coefficient is reported.
#' @param meas11 Character scalar naming the measurement-error variance column
#' for `predictor_u0`.
#' @param meas12 Character scalar naming the measurement-error covariance column
#' between `predictor_u0` and `predictor_u1`.
#' @param meas22 Character scalar naming the measurement-error variance column
#' for `predictor_u1`.
#' @param outcome_meas_var Optional character scalar naming an outcome
#' measurement-error variance column. The current estimator does not use this
#' variance in the estimating equations, but includes it in complete-case
#' filtering when supplied so corrected-outcome workflows use a consistent row
#' set.
#' @param stabilize_a_mat Logical; project the corrected normal-equation matrix
#' to positive definite before solving.
#' @param min_eigen Numeric lower bound used by positive-definite projections.
#' @param ridge_predictor_block Logical; add ridge only to the corrected
#' two-predictor cross-product block.
#' @param ridge_min_eigen Numeric minimum eigenvalue targeted by the predictor
#' block ridge.
#' @param measurement_weight Numeric multiplier applied to the supplied
#' predictor measurement-error covariance terms. The default `1` is full EIV;
#' values between `0` and `1` are tempered/regularized EIV variants.
#' @param regularize Logical; if `TRUE`, adaptively reduces
#' `measurement_weight` until the corrected latent predictor covariance is
#' positive definite. If `FALSE`, non-admissible corrected covariance returns
#' `NA` with `mx_issue_class = "corrected_predictor_cov_not_pd"`.
#' @param regularize_tol Numeric tolerance for the binary search used when
#' `regularize = TRUE`.
#' @param se_types Character vector of EIV standard-error variants to return.
#' Options are `"naive"`, `"hc0"`, and `"hc3"`. The HC0 variant corresponds to
#' the empirical sandwich used by earlier one-row versions of this function.
#'
#' @return
#' A tibble with one row per requested `se_type`, containing scaled `estimate`,
#' `se`, Wald-normal confidence limits, `status_code`, EIV diagnostics, and
#' issue-class fields.
#' Status `0L` indicates success, `1L` indicates that the corrected normal
#' equations could not be solved, `2L` indicates a non-finite scaled estimate or
#' standard error, and `3L` indicates that the corrected latent predictor
#' covariance was not positive definite.
fit_eiv_dual <- function(stage2_df,
                         outcome,
                         predictor_u0,
                         predictor_u1,
                         meas11,
                         meas12,
                         meas22,
                         outcome_meas_var = NULL,
                         stabilize_a_mat = FALSE,
                         min_eigen = 1e-6,
                         ridge_predictor_block = FALSE,
                         ridge_min_eigen = 1e-4,
                         measurement_weight = 1,
                         regularize = FALSE,
                         regularize_tol = 1e-6,
                         se_types = c("naive", "hc0", "hc3")) {
  se_types <- unique(match.arg(se_types, choices = c("naive", "hc0", "hc3"), several.ok = TRUE))

  latent_cov_diagnostics <- function(weight) {
    sum_s <- weight * matrix(c(sum(s11), sum(s12), sum(s12), sum(s22)), nrow = 2L, byrow = TRUE)
    sigma_x <- (crossprod(w_centered) - sum_s) / max(1, nrow(dat) - 1L)
    sigma_x <- (sigma_x + t(sigma_x)) / 2
    eig_values <- tryCatch(eigen(sigma_x, symmetric = TRUE, only.values = TRUE)$values, error = function(e) rep(NA_real_, 2L))
    min_eig <- min(eig_values, na.rm = TRUE)
    max_eig <- max(eig_values, na.rm = TRUE)
    condition_number <- if (is.finite(min_eig) && is.finite(max_eig) && min_eig > 0) {
      max_eig / min_eig
    } else {
      Inf
    }
    list(
      sigma_x = sigma_x,
      eig_values = eig_values,
      min_eig = min_eig,
      condition_number = condition_number,
      admissible = all(is.finite(eig_values)) && min_eig > min_eigen
    )
  }

  choose_measurement_weight <- function(initial_weight) {
    full_diag <- latent_cov_diagnostics(initial_weight)
    if (isTRUE(full_diag$admissible)) {
      return(list(weight = initial_weight, regularized = FALSE, diag = full_diag))
    }

    if (!isTRUE(regularize)) {
      return(list(weight = initial_weight, regularized = FALSE, diag = full_diag))
    }

    zero_diag <- latent_cov_diagnostics(0)
    if (!isTRUE(zero_diag$admissible)) {
      return(list(weight = 0, regularized = TRUE, diag = zero_diag))
    }

    lo <- 0
    hi <- initial_weight
    best_weight <- lo
    best_diag <- zero_diag

    for (iter in seq_len(60L)) {
      mid <- (lo + hi) / 2
      mid_diag <- latent_cov_diagnostics(mid)
      if (isTRUE(mid_diag$admissible)) {
        best_weight <- mid
        best_diag <- mid_diag
        lo <- mid
      } else {
        hi <- mid
      }
      if (abs(hi - lo) <= regularize_tol * max(1, abs(initial_weight))) {
        break
      }
    }

    list(weight = best_weight, regularized = TRUE, diag = best_diag)
  }

  out_fail <- tibble::tibble(
    se_type = se_types,
    estimate = NA_real_,
    se = NA_real_,
    ci_low = NA_real_,
    ci_high = NA_real_,
    status_code = NA_integer_,
    mx_issue_class = NA_character_,
    mx_issue_detail = NA_character_,
    eiv_measurement_weight_requested = measurement_weight,
    eiv_measurement_weight_used = NA_real_,
    eiv_regularized = FALSE,
    eiv_latent_cov_min_eigen = NA_real_,
    eiv_latent_cov_condition_number = NA_real_
  )

  # Include the optional outcome measurement-variance column only to enforce a
  # consistent complete-case sample when upstream simulations create it.
  cols_needed <- c(outcome, predictor_u0, predictor_u1, meas11, meas12, meas22)
  if (!is.null(outcome_meas_var)) {
    cols_needed <- c(cols_needed, outcome_meas_var)
  }

  dat <- stage2_df[, cols_needed, drop = FALSE]
  dat <- dat[stats::complete.cases(dat), , drop = FALSE]

  # EIV needs enough rows to estimate a three-parameter equation and its
  # sandwich variance with some stability. The target predictor must also vary.
  if (nrow(dat) < 8L || !is.finite(stats::sd(dat[[predictor_u1]])) ||
    stats::sd(dat[[predictor_u1]]) <= sqrt(.Machine$double.eps)) {
    return(out_fail)
  }

  w_mat <- as.matrix(dat[, c(predictor_u0, predictor_u1), drop = FALSE])
  y_vec <- dat[[outcome]]

  # Variance components should not reduce observed second moments because of
  # negative numerical estimates; clamp variances at zero but retain covariance.
  s11 <- pmax(dat[[meas11]], 0)
  s12 <- dat[[meas12]]
  s22 <- pmax(dat[[meas22]], 0)

  # Build the per-row measurement-error covariance array for
  # (intercept, predictor_u0, predictor_u1). The intercept row/column remains 0.
  x_mat <- cbind(1, w_mat)
  w_centered <- scale(w_mat, center = TRUE, scale = FALSE)

  weight_choice <- choose_measurement_weight(measurement_weight)
  measurement_weight_used <- weight_choice$weight
  latent_cov_diag <- weight_choice$diag

  if (!isTRUE(latent_cov_diag$admissible)) {
    return(dplyr::mutate(
      out_fail,
      status_code = 3L,
      mx_issue_class = "corrected_predictor_cov_not_pd",
      mx_issue_detail = sprintf(
        "min_eigen=%0.6e; condition_number=%0.6e; requested_weight=%0.6f; used_weight=%0.6f; regularize=%s",
        latent_cov_diag$min_eig,
        latent_cov_diag$condition_number,
        measurement_weight,
        measurement_weight_used,
        isTRUE(regularize)
      ),
      eiv_measurement_weight_used = measurement_weight_used,
      eiv_regularized = isTRUE(weight_choice$regularized),
      eiv_latent_cov_min_eigen = latent_cov_diag$min_eig,
      eiv_latent_cov_condition_number = latent_cov_diag$condition_number
    ))
  }

  # Corrected normal equations:
  #   sum_i (x_i x_i' - S_i) beta = sum_i x_i y_i
  # where S_i is the supplied measurement-error covariance for the predictors.
  a_mat <- crossprod(x_mat)
  a_mat[2, 2] <- a_mat[2, 2] - measurement_weight_used * sum(s11)
  a_mat[2, 3] <- a_mat[2, 3] - measurement_weight_used * sum(s12)
  a_mat[3, 2] <- a_mat[3, 2] - measurement_weight_used * sum(s12)
  a_mat[3, 3] <- a_mat[3, 3] - measurement_weight_used * sum(s22)

  b_vec <- as.vector(crossprod(x_mat, y_vec))

  a_mat_use <- a_mat
  if (isTRUE(stabilize_a_mat)) {
    a_mat_use <- project_to_pd(a_mat_use, min_eigen = min_eigen)
  }
  if (isTRUE(ridge_predictor_block)) {
    # Ridge only the observed predictor block, leaving the intercept coupling
    # intact while preventing near-singular corrected predictor covariance.
    pred_block <- (a_mat_use[2:3, 2:3, drop = FALSE] + t(a_mat_use[2:3, 2:3, drop = FALSE])) / 2
    eig_min <- min(eigen(pred_block, symmetric = TRUE, only.values = TRUE)$values)
    ridge_lambda <- max(0, ridge_min_eigen - eig_min)
    a_mat_use[2:3, 2:3] <- pred_block + diag(ridge_lambda, nrow = 2L)
  }

  beta_hat <- tryCatch(solve(a_mat_use, b_vec), error = function(e) NULL)
  if (is.null(beta_hat) || any(!is.finite(beta_hat))) {
    return(dplyr::mutate(
      out_fail,
      status_code = 1L,
      mx_issue_class = "corrected_normal_equations_solve_failed",
      eiv_measurement_weight_used = measurement_weight_used,
      eiv_regularized = isTRUE(weight_choice$regularized),
      eiv_latent_cov_min_eigen = latent_cov_diag$min_eig,
      eiv_latent_cov_condition_number = latent_cov_diag$condition_number
    ))
  }

  # If stabilization changed the estimating matrix, include the average matrix
  # perturbation in each influence contribution so the sandwich variance matches
  # the equation actually solved.
  a_adjustment <- (a_mat_use - a_mat) / nrow(dat)
  
  term1 <- x_mat * y_vec
  term2 <- x_mat * as.vector(x_mat %*% beta_hat)
  
  term3 <- matrix(0, nrow = nrow(dat), ncol = 3L)
  term3[, 2] <- measurement_weight_used * (s11 * beta_hat[2] + s12 * beta_hat[3])
  term3[, 3] <- measurement_weight_used * (s12 * beta_hat[2] + s22 * beta_hat[3])
  
  adj_beta <- as.vector(a_adjustment %*% beta_hat)
  term4 <- matrix(adj_beta, nrow = nrow(dat), ncol = 3L, byrow = TRUE)
  
  psi_mat <- term1 - term2 + term3 - term4

  # Empirical sandwich variance for beta_hat. The generalized inverse fallback
  # preserves a finite variance estimate in borderline stabilized cases.
  bread_inv <- tryCatch(solve(a_mat_use), error = function(e) MASS::ginv(a_mat_use))

  # Naive/model-based EIV variance treats the corrected normal-equation matrix
  # as fixed and uses a homoskedastic residual variance for Var(X'y).
  resid_vec <- y_vec - as.vector(x_mat %*% beta_hat)
  sigma2_hat <- sum(resid_vec^2) / max(1L, nrow(dat) - ncol(x_mat))
  vcov_naive <- sigma2_hat * bread_inv %*% crossprod(x_mat) %*% bread_inv

  meat_hc0 <- crossprod(psi_mat)
  vcov_hc0 <- bread_inv %*% meat_hc0 %*% bread_inv

  # HC3-style finite-sample correction for the EIV estimating functions. The
  # leverage is computed from the corrected bread and capped for numerical
  # stability, matching the defensive treatment used in the stacked sandwich.
  h_ii <- rowSums((x_mat %*% bread_inv) * x_mat)
  h_ii <- pmin(pmax(h_ii, 0), 0.999)
  psi_mat_hc3 <- psi_mat / pmax(1 - h_ii, sqrt(.Machine$double.eps))
  meat_hc3 <- crossprod(psi_mat_hc3)
  vcov_hc3 <- bread_inv %*% meat_hc3 %*% bread_inv

  # Estimate the latent predictor covariance by subtracting the aggregate
  # measurement-error covariance from the centered observed cross-product.
  sigma_x_hat <- latent_cov_diag$sigma_x
  scale_u1 <- sqrt(sigma_x_hat[2, 2])

  # Convert the raw latent-slope coefficient and its sandwich SE to the same
  # one-SD target scale used by the observed-score estimators.
  est <- unname(beta_hat[[3]]) * scale_u1
  vcov_by_type <- list(naive = vcov_naive, hc0 = vcov_hc0, hc3 = vcov_hc3)
  issue_detail <- if (isTRUE(weight_choice$regularized)) {
    sprintf(
      "regularized_weight=%0.6f; requested_weight=%0.6f; min_eigen=%0.6e; condition_number=%0.6e",
      measurement_weight_used,
      measurement_weight,
      latent_cov_diag$min_eig,
      latent_cov_diag$condition_number
    )
  } else {
    "ok"
  }

  purrr::map_dfr(se_types, function(se_type) {
    vcov_beta <- vcov_by_type[[se_type]]
    var_beta1 <- if (!is.null(vcov_beta) && nrow(vcov_beta) >= 3L) {
      unname(vcov_beta[3, 3])
    } else {
      NA_real_
    }
    se_beta1 <- if (is.finite(var_beta1) && var_beta1 >= 0) sqrt(var_beta1) else NA_real_
    se <- se_beta1 * scale_u1
    ok <- is.finite(est) && is.finite(se)

    tibble::tibble(
      se_type = se_type,
      estimate = if (ok) est else NA_real_,
      se = if (ok) se else NA_real_,
      ci_low = if (ok) est - stats::qnorm(0.975) * se else NA_real_,
      ci_high = if (ok) est + stats::qnorm(0.975) * se else NA_real_,
      status_code = if (ok) 0L else 2L,
      mx_issue_class = if (ok) "ok" else "nonfinite_eiv_estimate_or_se",
      mx_issue_detail = if (ok) issue_detail else NA_character_,
      eiv_measurement_weight_requested = measurement_weight,
      eiv_measurement_weight_used = measurement_weight_used,
      eiv_regularized = isTRUE(weight_choice$regularized),
      eiv_latent_cov_min_eigen = latent_cov_diag$min_eig,
      eiv_latent_cov_condition_number = latent_cov_diag$condition_number
    )
  })
}

fuller_dual_result_columns <- function() {
  c(
    "estimate", "se", "ci_low", "ci_high", "status_code",
    "fuller_raw_estimate", "fuller_raw_se",
    "mx_issue_class", "mx_issue_detail",
    "fuller_lambda1", "fuller_lambda2", "fuller_lambda_scaling",
    "fuller_sigma2", "fuller_weight_min", "fuller_weight_max",
    "fuller_correction1", "fuller_correction_c",
    "fuller_correction_scaling",
    "fuller_measurement_weight_requested",
    "fuller_measurement_weight_used",
    "fuller_preliminary_moment", "fuller_variance_bread",
    "fuller_predictor_outcome_covariance_source",
    "fuller_predictor_outcome_covariance_max_abs",
    "fuller_alpha_step1_requested", "fuller_alpha_step1_used",
    "fuller_alpha_step3_requested", "fuller_alpha_step3_used",
    "fuller_alpha_scaling_requested", "fuller_alpha_scaling_used",
    "fuller_auto_tempered",
    "fuller_sx1_star_condition", "fuller_sx1_star_min_eigen",
    "fuller_sx1_observed_max_eigen",
    "fuller_sx1_star_relative_min_eigen",
    "fuller_sx_star_condition", "fuller_sx_star_min_eigen",
    "fuller_sx_observed_max_eigen",
    "fuller_sx_star_relative_min_eigen",
    "fuller_scaling_condition", "fuller_scaling_min_eigen",
    "fuller_scaling_observed_max_eigen",
    "fuller_scaling_relative_min_eigen",
    "fuller_reference_se", "fuller_se_ratio",
    "fuller_auto_guard_pass", "fuller_auto_guard_reason",
    "fuller_auto_guard_score",
    "fuller_auto_full_weight_guard_pass",
    "fuller_auto_full_weight_guard_reason",
    "fuller_auto_full_weight_se_ratio",
    "fuller_auto_search_evaluations", "fuller_auto_search_nonmonotone"
  )
}

rescale_fuller_to_population_sd <- function(result, reporting_scale) {
  reporting_scale <- as.numeric(reporting_scale[[1]])
  if (!is.finite(reporting_scale) || reporting_scale <= 0) {
    stop("`reporting_scale` must be finite and positive.")
  }
  if (!all(c("fuller_raw_estimate", "fuller_raw_se") %in% names(result))) {
    stop("Fuller result is missing raw estimate or standard-error columns.")
  }

  dplyr::mutate(
    result,
    estimate = fuller_raw_estimate * reporting_scale,
    se = fuller_raw_se * reporting_scale,
    ci_low = estimate - stats::qnorm(0.975) * se,
    ci_high = estimate + stats::qnorm(0.975) * se
  )
}

fuller_matrix_diagnostics <- function(mat) {
  mat <- (mat + t(mat)) / 2
  vals <- tryCatch(eigen(mat, symmetric = TRUE, only.values = TRUE)$values, error = function(e) NA_real_)
  vals <- vals[is.finite(vals)]
  if (length(vals) == 0L) {
    return(list(min_eigen = NA_real_, max_eigen = NA_real_, condition_number = Inf))
  }
  min_eig <- min(vals)
  max_eig <- max(vals)
  condition_number <- if (is.finite(min_eig) && is.finite(max_eig) && min_eig > 0) {
    max_eig / min_eig
  } else {
    Inf
  }
  list(min_eigen = min_eig, max_eigen = max_eig, condition_number = condition_number)
}

fuller_relative_min_eigen <- function(corrected_diag, observed_diag) {
  if (is.finite(corrected_diag$min_eigen) &&
    is.finite(observed_diag$max_eigen) &&
    observed_diag$max_eigen > sqrt(.Machine$double.eps)) {
    corrected_diag$min_eigen / observed_diag$max_eigen
  } else {
    NA_real_
  }
}

fuller_row_max_finite <- function(...) {
  mat <- cbind(...)
  apply(mat, 1L, function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0L) Inf else max(x)
  })
}

fuller_guard_penalty <- function(value, floor) {
  ifelse(
    is.finite(value) & is.finite(floor) & floor > 0 & value > 0,
    pmax(0, log(floor / value)),
    Inf
  )
}

#' Fit the shared Fuller errors-in-variables (EIV) core.
#'
#' @details
#' Implements the three-step procedure from Section \ref{sec-fuller-method} of
#' `documentation/re_regression_vig_5-12-26.tex` (after Fuller, 1987) for the
#' case where both the stage-2 outcome and the stage-2 predictors are measured
#' with error.
#'
#' The stage-2 regression uses predictors ordered as
#' `(intercept, predictor_u0, predictor_u1)` when `predictor_u0` is supplied.
#' If `predictor_u0` is `NULL`, the regression uses `(intercept, predictor_u1)`.
#' The intercept is treated as measured without error. The supplied predictor
#' measurement-error variances and covariance (`meas11`, `meas12`, `meas22`)
#' populate the predictor block. Outcome measurement error is supplied via
#' `outcome_meas_var`. Predictor--outcome measurement-error covariances may be
#' supplied explicitly via `predictor_outcome_meas_cov_u0` and
#' `predictor_outcome_meas_cov_u1`. Omitted covariance columns are represented
#' by exact zeros. This is the appropriate VH Studies 1--4 default because the
#' Stage-1 predictor and outcome errors arise from separate processes (or the
#' Stage-2 outcome is observed without error).
#'
#' There are two independent algebra switches; neither changes the assumed
#' measurement-error model. `preliminary_moment` controls the preliminary
#' coefficient used to estimate the structural residual variance and case
#' weights. `variance_bread` controls the matrix inverted in the final Wald
#' variance. The final point estimator always uses Fuller's lambda/alpha-modified
#' weighted moment equation. The VH primary estimator uses `"modified"` for
#' both switches. `fit_fuller_dual_variants()` exposes a nested comparison that
#' first replaces the preliminary moment with Fuller (1987, Eq. 3.1.19) and
#' then replaces the variance bread with Eq. 3.1.28. This ordering identifies
#' whether preliminary-weight regularization or bread regularization controls
#' the extreme-SE tail.
#'
#' The returned estimate and standard error are scaled by an estimated latent
#' SD for `predictor_u1`, computed by subtracting the aggregate predictor
#' measurement-error covariance from the centered observed cross-product.
#'
#' @param stage2_df Data frame containing outcome, predictors, and
#' measurement-error covariance columns. This function assumes one row per
#' cluster/unit contributing to the second-stage regression.
#' @param outcome Character scalar naming the outcome column.
#' @param predictor_u0 Optional character scalar naming the intercept-like
#' observed predictor. If `NULL`, the model uses only `predictor_u1`.
#' @param predictor_u1 Character scalar naming the slope-like observed
#' predictor whose coefficient is reported.
#' @param meas11 Optional character scalar naming the predictor measurement-error
#' variance column for `predictor_u0`. Required when `predictor_u0` is supplied.
#' @param meas12 Optional character scalar naming the predictor measurement-error
#' covariance column between `predictor_u0` and `predictor_u1`. Required when
#' `predictor_u0` is supplied.
#' @param meas22 Character scalar naming the predictor measurement-error
#' variance column for `predictor_u1`.
#' @param outcome_meas_var Optional character scalar naming an outcome
#' measurement-error variance column. If `NULL`, outcome measurement error is
#' treated as zero.
#' @param predictor_outcome_meas_cov_u0 Optional character scalar naming the
#' measurement-error covariance between `predictor_u0` and the outcome. It must
#' be omitted when `predictor_u0` is `NULL`. If omitted, this covariance is zero.
#' @param predictor_outcome_meas_cov_u1 Optional character scalar naming the
#' measurement-error covariance between `predictor_u1` and the outcome. If
#' omitted, this covariance is zero.
#' @param preliminary_moment Either `"modified"` (the VH default), which uses
#' the lambda/alpha-modified Step-1 moments, or `"fuller"`, which uses the
#' fully measurement-error-subtracted preliminary moments in Fuller (1987,
#' Eq. 3.1.19). Both choices use the predictor--outcome covariance in the
#' corrected cross-moment; outcome error variance is never substituted for it.
#' @param variance_bread Either `"modified"` (the VH default), which uses the
#' lambda/alpha-modified final predictor moment as the sandwich bread, or
#' `"fuller"`, which uses the fully measurement-error-subtracted weighted
#' predictor moment in Fuller (1987, Eq. 3.1.28). The two switches expose the
#' stabilizations identified in the Fuller algebra audit as method variants.
#' @param alpha_step1 Optional numeric Step 1 correction factor. Defaults to
#' `p + 1` after data-dependent `p` is known.
#' @param alpha_step3 Optional numeric Step 3 correction factor. Defaults to
#' `p + 1` after data-dependent `p` is known.
#'
#' @return
#' A one-row tibble containing the scaled `estimate`, `se`, Wald-normal
#' confidence limits, `status_code`, and Fuller step diagnostics. Status `0L`
#' indicates success, `1L` indicates a linear-system solve failure in the
#' Fuller correction steps, `2L` indicates non-finite estimates/standard errors
#' or a negative/non-finite variance estimate, and `3L` indicates an
#' admissibility failure such as invalid `measurement_weight`, insufficient
#' complete cases, invalid alpha values, or an invalid scaling denominator
fit_fuller_dual_core <- function(stage2_df,
                                 outcome,
                                 predictor_u0 = NULL,
                                 predictor_u1,
                                 meas11 = NULL,
                                 meas12 = NULL,
                                 meas22,
                                 outcome_meas_var = NULL,
                                 predictor_outcome_meas_cov_u0 = NULL,
                                 predictor_outcome_meas_cov_u1 = NULL,
                                 measurement_weight = 1,
                                 alpha_step1 = NULL,
                                 alpha_step3 = NULL,
                                 skip_internal_scaling = TRUE,
                                 alpha_scaling = NULL,
                                 preliminary_moment = c("modified", "fuller"),
                                 variance_bread = c("modified", "fuller"),
                                 auto_tempered = FALSE) {
  if (!requireNamespace("geigen", quietly = TRUE)) {
    stop(
      "The `geigen` package is required for Fuller EIV estimation (generalized eigenvalues). ",
      "Install it with `install.packages(\"geigen\")`."
    )
  }

  preliminary_moment <- match.arg(preliminary_moment)
  variance_bread <- match.arg(variance_bread)

  # Record whether the caller supplied predictor--outcome error covariance.
  # Missing columns mean exact zero covariance, never a request to reuse
  # `outcome_meas_var`. VH Studies 1--4 pass named zero columns so that this
  # assumption remains visible in saved replication results.
  covariance_source <- if (
    is.null(predictor_outcome_meas_cov_u0) &&
      is.null(predictor_outcome_meas_cov_u1)
  ) {
    "zero_default"
  } else {
    "supplied"
  }

  out_fail <- tibble::tibble(
    estimate = NA_real_,
    se = NA_real_,
    ci_low = NA_real_,
    ci_high = NA_real_,
    status_code = NA_integer_,
    fuller_raw_estimate = NA_real_,
    fuller_raw_se = NA_real_,
    mx_issue_class = NA_character_,
    mx_issue_detail = NA_character_,
    fuller_lambda1 = NA_real_,
    fuller_lambda2 = NA_real_,
    fuller_lambda_scaling = NA_real_,
    fuller_sigma2 = NA_real_,
    fuller_weight_min = NA_real_,
    fuller_weight_max = NA_real_,
    fuller_correction_c = NA_real_,
    fuller_correction1 = NA_real_,
    fuller_correction_scaling = NA_real_,
    fuller_measurement_weight_requested = measurement_weight,
    fuller_measurement_weight_used = measurement_weight,
    fuller_preliminary_moment = preliminary_moment,
    fuller_variance_bread = variance_bread,
    fuller_predictor_outcome_covariance_source = covariance_source,
    fuller_predictor_outcome_covariance_max_abs = NA_real_,
    fuller_alpha_step1_requested = NA_real_,
    fuller_alpha_step1_used = NA_real_,
    fuller_alpha_step3_requested = NA_real_,
    fuller_alpha_step3_used = NA_real_,
    fuller_alpha_scaling_requested = NA_real_,
    fuller_alpha_scaling_used = NA_real_,
    fuller_auto_tempered = auto_tempered,
    fuller_sx1_star_condition = NA_real_,
    fuller_sx1_star_min_eigen = NA_real_,
    fuller_sx1_star_relative_min_eigen = NA_real_,
    fuller_sx1_observed_max_eigen = NA_real_,
    fuller_sx_star_condition = NA_real_,
    fuller_sx_star_min_eigen = NA_real_,
    fuller_sx_star_relative_min_eigen = NA_real_,
    fuller_sx_observed_max_eigen = NA_real_,
    fuller_scaling_condition = NA_real_,
    fuller_scaling_min_eigen = NA_real_,
    fuller_scaling_observed_max_eigen = NA_real_,
    fuller_scaling_relative_min_eigen = NA_real_,
    fuller_reference_se = NA_real_,
    fuller_se_ratio = NA_real_,
    fuller_auto_guard_pass = NA,
    fuller_auto_guard_reason = NA_character_,
    fuller_auto_guard_score = NA_real_,
    fuller_auto_full_weight_guard_pass = NA,
    fuller_auto_full_weight_guard_reason = NA_character_,
    fuller_auto_full_weight_se_ratio = NA_real_,
    fuller_auto_search_evaluations = NA_integer_,
    fuller_auto_search_nonmonotone = NA
  )

  fuller_fail <- function(status_code, issue_class, issue_detail = NA_character_) {
    dplyr::mutate(
      out_fail,
      status_code = status_code,
      mx_issue_class = issue_class,
      mx_issue_detail = issue_detail
    )
  }

  has_u0 <- !is.null(predictor_u0)
  if (has_u0 && (is.null(meas11) || is.null(meas12))) {
    stop("`meas11` and `meas12` must be supplied when `predictor_u0` is provided.")
  }
  if (!has_u0 && (!is.null(meas11) || !is.null(meas12))) {
    stop("`meas11` and `meas12` should be omitted when `predictor_u0` is NULL.")
  }
  if (!has_u0 && !is.null(predictor_outcome_meas_cov_u0)) {
    stop(
      "`predictor_outcome_meas_cov_u0` should be omitted when `predictor_u0` is NULL."
    )
  }

  cols_needed <- c(outcome, predictor_u1, meas22)
  if (has_u0) {
    cols_needed <- c(outcome, predictor_u0, predictor_u1, meas11, meas12, meas22)
  }
  if (!is.finite(measurement_weight) || measurement_weight < 0 || measurement_weight > 1) {
    return(fuller_fail(
      3L,
      "fuller_invalid_measurement_weight",
      sprintf("measurement_weight=%0.6f", measurement_weight)
    ))
  }

  if (!is.null(outcome_meas_var)) {
    cols_needed <- c(cols_needed, outcome_meas_var)
  }
  if (!is.null(predictor_outcome_meas_cov_u0)) {
    cols_needed <- c(cols_needed, predictor_outcome_meas_cov_u0)
  }
  if (!is.null(predictor_outcome_meas_cov_u1)) {
    cols_needed <- c(cols_needed, predictor_outcome_meas_cov_u1)
  }
  cols_needed <- unique(cols_needed)

  dat <- stage2_df[, cols_needed, drop = FALSE]
  dat <- dat[stats::complete.cases(dat), , drop = FALSE]

  # P counts the intercept in the second-stage regression.
  p <- if (has_u0) 3L else 2L
  m <- nrow(dat)
  # EIV needs enough rows to estimate a p-parameter equation and its sandwich
  # variance with some stability. The target predictor must also vary.
  if (m < 8L || m <= p || !is.finite(stats::sd(dat[[predictor_u1]])) ||
    stats::sd(dat[[predictor_u1]]) <= sqrt(.Machine$double.eps)) {
    return(fuller_fail(
      3L,
      "fuller_insufficient_complete_cases",
      sprintf("m=%d; p=%d", m, p)
    ))
  }

  alpha_step1 <- if (is.null(alpha_step1)) p + 1 else alpha_step1
  alpha_step3 <- if (is.null(alpha_step3)) p + 1 else alpha_step3
  
  out_fail <- dplyr::mutate(
    out_fail,
    fuller_alpha_step1_requested = alpha_step1,
    fuller_alpha_step3_requested = alpha_step3
  )
  if (!is.finite(alpha_step1) || !is.finite(alpha_step3) ||
    alpha_step1 <= 0 || alpha_step3 <= 0) {
    return(fuller_fail(
      3L,
      "fuller_invalid_alpha",
      sprintf(
        "alpha_step1=%s; alpha_step3=%s",
        alpha_step1,
        alpha_step3
      )
    ))
  }
  
  if (!skip_internal_scaling) {
    alpha_scaling <- if (is.null(alpha_scaling)) p + 1 else alpha_scaling
    out_fail <- dplyr::mutate(
      out_fail,
      fuller_alpha_scaling_requested = alpha_scaling
    )
    if (!is.finite(alpha_scaling) || alpha_scaling <= 0) {
      return(fuller_fail(
        3L,
        "fuller_invalid_alpha",
        sprintf(
          "alpha_scaling=%s",
          alpha_scaling
        )
      ))
    }
    alpha_scaling <- max(alpha_scaling, p + 1)
  } else {
    alpha_scaling <- NA
  }
  
  alpha_step1 <- max(alpha_step1, p + 1)
  alpha_step3 <- max(alpha_step3, p + 1)
  out_fail <- dplyr::mutate(
    out_fail,
    fuller_alpha_step1_used = alpha_step1,
    fuller_alpha_step3_used = alpha_step3,
    fuller_alpha_scaling_used = alpha_scaling
  )

  # Extract observed scores. `has_u0` changes only the dimension and ordering of
  # the predictor vector; it is unrelated to the modified-versus-Fuller algebra
  # choices controlled by `preliminary_moment` and `variance_bread`.
  y_vec <- dat[[outcome]]
  u1_vec <- dat[[predictor_u1]]
  u0_vec <- if (has_u0) dat[[predictor_u0]] else NULL
  x_mat <- if (has_u0) cbind(1, u0_vec, u1_vec) else cbind(1, u1_vec)
  pred_block_idx <- if (has_u0) 2:3 else 2:2

  # Build the predictor error block Omega_xj. The dual-predictor path uses
  # (intercept, u0, u1); the single-predictor path uses (intercept, u1).
  # The design-matrix intercept is measured without error in both paths. Clamp
  # variances at zero, but retain the signed u0--u1 covariance.
  if (has_u0) {
    s11 <- measurement_weight * pmax(dat[[meas11]], 0)
    s12 <- measurement_weight * dat[[meas12]]
    s22 <- measurement_weight * pmax(dat[[meas22]], 0)
  } else {
    s11 <- rep(0, m)
    s12 <- rep(0, m)
    s22 <- measurement_weight * pmax(dat[[meas22]], 0)
  }

  # Build the remaining pieces of the full measurement-error covariance:
  # outcome variance omega_yj and predictor--outcome vector omega_xyj. These
  # are distinct: positive outcome error variance does not imply nonzero
  # predictor--outcome covariance.
  # `measurement_weight` tempers every predictor-involving correction together;
  # it is one in the primary and algebra-comparison estimators.
  omega_y <- if (!is.null(outcome_meas_var)) pmax(dat[[outcome_meas_var]], 0) else rep(0, m)
  omega_xy_u0 <- if (!is.null(predictor_outcome_meas_cov_u0)) {
    measurement_weight * dat[[predictor_outcome_meas_cov_u0]]
  } else {
    rep(0, m)
  }
  omega_xy_u1 <- if (!is.null(predictor_outcome_meas_cov_u1)) {
    measurement_weight * dat[[predictor_outcome_meas_cov_u1]]
  } else {
    rep(0, m)
  }
  # Stack each cross-covariance vector in x_mat order. The leading zero belongs
  # to the error-free design-matrix intercept.
  omega_xy_mat <- if (has_u0) {
    cbind(0, omega_xy_u0, omega_xy_u1)
  } else {
    cbind(0, omega_xy_u1)
  }
  omega_xy_sum <- colSums(omega_xy_mat)
  omega_y_sum <- sum(omega_y)
  out_fail <- dplyr::mutate(
    out_fail,
    fuller_predictor_outcome_covariance_max_abs = max(abs(omega_xy_mat))
  )

  # Sum Omega_xj in x_mat order. These branches are dimensional bookkeeping for
  # dual versus single predictors, not competing Fuller estimators.
  if (has_u0) {
    omega_x_sum <- matrix(
      c(
        0, 0, 0,
        0, sum(s11), sum(s12),
        0, sum(s12), sum(s22)
      ),
      nrow = p,
      byrow = TRUE
    )
  } else {
    omega_x_sum <- matrix(c(0, 0, 0, sum(s22)), nrow = p, byrow = TRUE)
  }

  # Helper: smallest finite generalized eigenvalue of (A, B) with B
  # positive-semidefinite but possibly singular.
  smallest_det_root <- function(a_mat, b_mat) {
    # We want the smallest root of det(A - lambda B) = 0.
    # For symmetric A and PSD B, these roots are the generalized eigenvalues
    # solving A v = lambda B v.
    #
    # Directly computing det(A - lambda B) is numerically unstable, so we use
    # a generalized eigen-solver. The `geigen` package wraps LAPACK routines
    # for the generalized eigenproblem.
    #
    # Note: `geigen(..., symmetric=TRUE)` requires B to be positive definite.
    # In Fuller, B is often only positive semidefinite (and can be singular),
    # so we explicitly set `symmetric = FALSE`.
    ge_out <- tryCatch(
      geigen::geigen(a_mat, b_mat, symmetric = FALSE, only.values = TRUE),
      error = function(e) NULL
    )
    if (is.null(ge_out) || is.null(ge_out$values) || length(ge_out$values) == 0L) {
      return(NA_real_)
    }

    vals <- ge_out$values
    # In non-PD/singular cases LAPACK can return complex-valued eigenvalues.
    # Treat tiny imaginary parts as numerical noise, otherwise drop.
    if (is.complex(vals)) {
      imag_tol <- sqrt(.Machine$double.eps) * max(1, max(abs(vals), na.rm = TRUE))
      vals <- ifelse(abs(Im(vals)) <= imag_tol, Re(vals), NA_real_)
    }

    vals <- vals[is.finite(vals)]
    if (length(vals) == 0L) {
      return(NA_real_)
    }

    # Returning the smallest matches the "smallest determinant root" used for
    # lambda_1 and lambda_2 in the Fuller procedure.
    min(vals)
  }

  # Step 1: lambda_1.
  b_mat <- cbind(y_vec, x_mat)
  bb_sum <- crossprod(b_mat)
  omega_sum <- matrix(0, nrow = p + 1L, ncol = p + 1L)
  omega_sum[1, 1] <- omega_y_sum
  omega_sum[1, 2:(p + 1L)] <- omega_xy_sum
  omega_sum[2:(p + 1L), 1] <- omega_xy_sum
  omega_sum[2:(p + 1L), 2:(p + 1L)] <- omega_x_sum

  lambda1_hat <- smallest_det_root(bb_sum, omega_sum)

  # Step 1 obtains the lambda/alpha correction used by the stabilized
  # preliminary path (and retained in the diagnostics for either path).
  alpha_step1_scaled <- alpha_step1 / m
  c_correction1 <- if (!is.na(lambda1_hat) && is.finite(lambda1_hat) && lambda1_hat <= 1 + 1 / m) {
    lambda1_hat - 1 / m - alpha_step1_scaled
  } else {
    1 - alpha_step1_scaled
  }

  # PRELIMINARY-MOMENT SWITCH
  # * "modified" (VH primary) uses S1* and subtracts c1 times the full
  #   predictor/cross-error moments. This regularizes gamma0_hat, which then
  #   enters the residual-variance estimate, case weights, and variance meat.
  # * "fuller" uses the unmodified moment equation in Fuller (1987),
  #   Eq. 3.1.19, with coefficient 1 on those same error moments.
  # Both paths subtract omega_xy; neither substitutes outcome variance omega_y.
  if (identical(preliminary_moment, "modified")) {
    a0_mat <- crossprod(x_mat) - c_correction1 * omega_x_sum
    b0_vec <- as.vector(
      crossprod(x_mat, y_vec) - c_correction1 * omega_xy_sum
    )
  } else {
    a0_mat <- crossprod(x_mat) - omega_x_sum
    b0_vec <- as.vector(crossprod(x_mat, y_vec) - omega_xy_sum)
  }

  sx1_observed <- crossprod(x_mat) / m
  sx1_observed_diag <- fuller_matrix_diagnostics(sx1_observed[pred_block_idx, pred_block_idx, drop = FALSE])
  sx1_star <- (a0_mat / m)[pred_block_idx, pred_block_idx, drop = FALSE]
  sx1_star_diag <- fuller_matrix_diagnostics(sx1_star)
  sx1_star_relative_min_eigen <- fuller_relative_min_eigen(sx1_star_diag, sx1_observed_diag)

  out_fail <- dplyr::mutate(
    out_fail,
    fuller_lambda1 = lambda1_hat,
    fuller_correction1 = c_correction1,
    fuller_sx1_star_relative_min_eigen = sx1_star_relative_min_eigen,
    fuller_sx1_observed_max_eigen = sx1_observed_diag$max_eigen,
    fuller_sx1_star_condition = sx1_star_diag$condition_number,
    fuller_sx1_star_min_eigen = sx1_star_diag$min_eigen
  )
  gamma0_hat <- tryCatch(as.vector(solve(a0_mat, b0_vec)), error = function(e) NULL)
  if (is.null(gamma0_hat) || any(!is.finite(gamma0_hat))) {
    return(fuller_fail(1L, "fuller_gamma0_solve_failed"))
  }
  gamma0_u0 <- if (has_u0) gamma0_hat[2] else 0
  gamma0_u1 <- gamma0_hat[p]

  # Regression error variance estimate: SSE/(M-P) minus average measurement
  # error in the composite residual (u_y - gamma0' u_x).
  resid0 <- y_vec - as.vector(x_mat %*% gamma0_hat)
  sigma2_ols <- sum(resid0^2) / max(1, m - p)
  predictor_outcome_cross <- gamma0_u0 * omega_xy_u0 +
    gamma0_u1 * omega_xy_u1
  quad_x <- (gamma0_u0^2) * s11 +
    2 * gamma0_u0 * gamma0_u1 * s12 +
    (gamma0_u1^2) * s22
  composite_measurement_variance <- omega_y -
    2 * predictor_outcome_cross + quad_x
  sigma2_corr <- mean(composite_measurement_variance)

  sigma2_hat <- if (!is.na(lambda1_hat) && is.finite(lambda1_hat) && lambda1_hat < 1) {
    0
  } else {
    sigma2_ols - sigma2_corr
  }
  if (!is.finite(sigma2_hat)) {
    return(fuller_fail(2L, "fuller_sigma2_nonfinite"))
  }
  sigma2_hat <- max(0, sigma2_hat)

  # Step 3: weights, lambda_2, corrected S* matrix, and final gamma.
  w_j <- sigma2_hat + composite_measurement_variance
  out_fail <- dplyr::mutate(
    out_fail,
    fuller_lambda1 = lambda1_hat,
    fuller_sigma2 = sigma2_hat,
    fuller_weight_min = suppressWarnings(min(w_j, na.rm = TRUE)),
    fuller_weight_max = suppressWarnings(max(w_j, na.rm = TRUE))
  )

  if (any(!is.finite(w_j)) || any(w_j <= sqrt(.Machine$double.eps))) {
    return(fuller_fail(
      2L,
      "fuller_nonpositive_weights",
      sprintf(
        "min_w=%0.6e; max_w=%0.6e; sigma2=%0.6e",
        suppressWarnings(min(w_j, na.rm = TRUE)),
        suppressWarnings(max(w_j, na.rm = TRUE)),
        sigma2_hat
      )
    ))
  }

  w_inv <- 1 / w_j
  bw_sum <- crossprod(b_mat * sqrt(w_inv))

  if (has_u0) {
    omega_x_sum_w <- matrix(
      c(
        0, 0, 0,
        0, sum(w_inv * s11), sum(w_inv * s12),
        0, sum(w_inv * s12), sum(w_inv * s22)
      ),
      nrow = p,
      byrow = TRUE
    )
  } else {
    omega_x_sum_w <- matrix(c(0, 0, 0, sum(w_inv * s22)), nrow = p, byrow = TRUE)
  }
  omega_sum_w <- matrix(0, nrow = p + 1L, ncol = p + 1L)
  omega_sum_w[1, 1] <- sum(w_inv * omega_y)
  omega_xy_sum_w <- colSums(omega_xy_mat * w_inv)
  omega_sum_w[1, 2:(p + 1L)] <- omega_xy_sum_w
  omega_sum_w[2:(p + 1L), 1] <- omega_xy_sum_w
  omega_sum_w[2:(p + 1L), 2:(p + 1L)] <- omega_x_sum_w

  lambda2_hat <- smallest_det_root(bw_sum, omega_sum_w)
  alpha_step3_scaled <- alpha_step3 / m
  c_correction <- if (!is.na(lambda2_hat) && is.finite(lambda2_hat) && lambda2_hat <= 1 + 1 / m) {
    lambda2_hat - 1 / m - alpha_step3_scaled
  } else {
    1 - alpha_step3_scaled
  }

  s_star <- bw_sum - c_correction * omega_sum_w
  s_x_star <- s_star[2:(p + 1L), 2:(p + 1L), drop = FALSE]
  s_xy_star <- s_star[2:(p + 1L), 1, drop = FALSE]
  sx_star_diag <- fuller_matrix_diagnostics(s_x_star[pred_block_idx, pred_block_idx, drop = FALSE])

  out_fail <- dplyr::mutate(
    out_fail,
    fuller_lambda2 = lambda2_hat,
    fuller_correction_c = c_correction,
    fuller_sx_star_condition = sx_star_diag$condition_number,
    fuller_sx_star_min_eigen = sx_star_diag$min_eigen
  )

  gamma_hat <- tryCatch(as.vector(solve(s_x_star, s_xy_star)), error = function(e) NULL)
  if (is.null(gamma_hat) || any(!is.finite(gamma_hat))) {
    return(fuller_fail(1L, "fuller_gamma_solve_failed"))
  }

  # Standard errors (Fuller, 1987).
  xw_sum <- crossprod(x_mat * sqrt(w_inv))
 
  s_x_observed <- xw_sum 
  sx_observed_diag <- fuller_matrix_diagnostics(s_x_observed[pred_block_idx, pred_block_idx, drop = FALSE])
  sx_star_relative_min_eigen <- fuller_relative_min_eigen(sx_star_diag, sx_observed_diag)
  out_fail <- dplyr::mutate(
    out_fail,
    fuller_sx_star_relative_min_eigen = sx_star_relative_min_eigen,
    fuller_sx_observed_max_eigen = sx_observed_diag$max_eigen
  )
  # VARIANCE-BREAD SWITCH
  # * "modified" (VH primary) inverts the lambda/alpha-modified final Sx*.
  #   This is the second stabilization and is intended to tame the extreme-SE
  #   tail when the corrected predictor moment is nearly singular.
  # * "fuller" uses the unmodified corrected bread from Fuller (1987),
  #   Eq. 3.1.28: sum(w_j x_j x_j') - sum(w_j Omega_xj).
  # Holding `preliminary_moment` fixed, this switch changes Wald SEs only; the
  # final point estimate above always comes from the same modified equation.
  variance_bread_mat <- if (identical(variance_bread, "modified")) {
    s_x_star
  } else {
    xw_sum - omega_x_sum_w
  }
  variance_bread_scaled <- variance_bread_mat / m
  s_x_inv <- tryCatch(solve(variance_bread_scaled), error = function(e) NULL)
  if (is.null(s_x_inv) || any(!is.finite(s_x_inv))) {
    return(fuller_fail(1L, "fuller_sx_star_solve_failed"))
  }

  # Fuller's variance meat uses
  #   tilde_omega_j = omega_xyj - Omega_xj gamma0.
  # The first component is exactly zero because the design-matrix intercept is
  # measured without error. The branches below only expand this expression for
  # dual (u0, u1) versus slope-only (u1) predictors; they are not distinct
  # estimator variants.
  if (has_u0) {
    tilde_mat <- cbind(
      0,
      omega_xy_u0 - (s11 * gamma0_u0 + s12 * gamma0_u1),
      omega_xy_u1 - (s12 * gamma0_u0 + s22 * gamma0_u1)
    )
  } else {
    tilde_mat <- cbind(0, omega_xy_u1 - s22 * gamma0_u1)
  }

  meat_sum <- xw_sum + crossprod(tilde_mat * w_inv) # w_inv is actually squared here
  vcov_gamma <- (s_x_inv %*% meat_sum %*% s_x_inv) / (m^2)
  vcov_gamma <- (vcov_gamma + t(vcov_gamma)) / 2 # trick to ensure symmetry
  target_gamma_idx <- p
  var_u1 <- unname(vcov_gamma[target_gamma_idx, target_gamma_idx])
  if (!is.finite(var_u1) || var_u1 < 0) {
    return(fuller_fail(2L, "fuller_negative_or_nonfinite_variance"))
  }

  
  if (!skip_internal_scaling) {
    # Scale to the one-latent-SD target used elsewhere in stage-2 summaries.
    # The Fuller estimating equations use the step-3 weights w_j, which downweight
    # observations with large composite measurement variance. Using the same
    # weights for the latent-SD scaling tends to be more stable than the raw
    # unweighted corrected covariance when the supplied measurement-error
    # variances are large.
    # https://en.wikipedia.org/wiki/Weighted_arithmetic_mean#
    pred_mat <- if (has_u0) cbind(u0_vec, u1_vec) else matrix(u1_vec, ncol = 1L)
    w_sum <- sum(w_inv)
    w_sq_sum <- sum(w_inv^2)
    denom_eff <- w_sum - w_sq_sum / w_sum
    pred_mean <- colSums(pred_mat * w_inv) / w_sum
    pred_centered <- sweep(pred_mat, 2L, pred_mean, FUN = "-")
    pred_centered_w <- pred_centered * sqrt(w_inv)
  
    if (has_u0) {
      omega_x_sum_pred_w <- matrix(
        c(sum(w_inv * s11), sum(w_inv * s12), sum(w_inv * s12), sum(w_inv * s22)),
        nrow = 2L,
        byrow = TRUE
      )
    } else {
      omega_x_sum_pred_w <- matrix(sum(w_inv * s22), nrow = 1L)
    }
  
    # we have to choose a third alpha since we're now working with centered x
    lambda_scaling_hat <- smallest_det_root(crossprod(pred_centered_w), omega_x_sum_pred_w)
    c_correction_scaling <- if (!is.na(lambda_scaling_hat) && is.finite(lambda_scaling_hat) && lambda_scaling_hat <= 1 + 1 / m) {
      lambda_scaling_hat - 1 / m - alpha_scaling / m
    } else {
      1 - alpha_scaling / m
    }
  
    # Use the standard effective-denominator for a weighted sample covariance so
    # that constant weights reduce exactly to the unweighted (m - 1) denominator.
    sigma_x_observed <- crossprod(pred_centered_w) / denom_eff
    sigma_x_observed <- (sigma_x_observed + t(sigma_x_observed)) / 2
    scaling_observed_diag <- fuller_matrix_diagnostics(sigma_x_observed)
    # Apply the Step-3 correction to the centered predictor covariance.
    sigma_x_hat <- (crossprod(pred_centered_w) - c_correction_scaling * omega_x_sum_pred_w) / denom_eff
    sigma_x_hat <- (sigma_x_hat + t(sigma_x_hat)) / 2
    scaling_diag <- fuller_matrix_diagnostics(sigma_x_hat)
    scaling_relative_min_eigen <- fuller_relative_min_eigen(scaling_diag, scaling_observed_diag)
  
    out_fail <- dplyr::mutate(
      out_fail,
      fuller_lambda_scaling = lambda_scaling_hat,
      fuller_correction_scaling = c_correction_scaling,
      fuller_scaling_condition = scaling_diag$condition_number,
      fuller_scaling_min_eigen = scaling_diag$min_eigen,
      fuller_scaling_observed_max_eigen = scaling_observed_diag$max_eigen,
      fuller_scaling_relative_min_eigen = scaling_relative_min_eigen
    )
    
    target_pred_idx <- ncol(pred_mat)
    if (!is.finite(sigma_x_hat[target_pred_idx, target_pred_idx]) ||
      sigma_x_hat[target_pred_idx, target_pred_idx] <= sqrt(.Machine$double.eps)) {
      return(fuller_fail(2L, "fuller_latent_predictor_var_not_positive"))
    }
    scale_u1 <- sqrt(sigma_x_hat[target_pred_idx, target_pred_idx])
  
    est <- unname(gamma_hat[target_gamma_idx]) * scale_u1
    se <- sqrt(var_u1) * scale_u1
    if (!is.finite(est) || !is.finite(se)) {
      return(fuller_fail(2L, "fuller_nonfinite_estimate_or_se"))
    }
  } else {
    est <- unname(gamma_hat[target_gamma_idx])
    se <- sqrt(var_u1)
    lambda_scaling_hat <- NA
    c_correction_scaling <- NA
    scaling_diag <- data.frame(condition_number = NA, min_eigen = NA)
    scaling_observed_diag <- data.frame(max_eigen = NA)
    scaling_relative_min_eigen <- NA
  }

  tibble::tibble(
    estimate = est,
    se = se,
    ci_low = est - stats::qnorm(0.975) * se,
    ci_high = est + stats::qnorm(0.975) * se,
    status_code = 0L,
    fuller_raw_estimate = unname(gamma_hat[target_gamma_idx]),
    fuller_raw_se = sqrt(var_u1),
    mx_issue_class = "ok",
    mx_issue_detail = "ok",
    fuller_lambda1 = lambda1_hat,
    fuller_lambda2 = lambda2_hat,
    fuller_lambda_scaling = lambda_scaling_hat,
    fuller_sigma2 = sigma2_hat,
    fuller_weight_min = min(w_j),
    fuller_weight_max = max(w_j),
    fuller_correction1 = c_correction1,
    fuller_correction_c = c_correction,
    fuller_correction_scaling = c_correction_scaling,
    fuller_measurement_weight_requested = measurement_weight,
    fuller_measurement_weight_used = measurement_weight,
    fuller_preliminary_moment = preliminary_moment,
    fuller_variance_bread = variance_bread,
    fuller_predictor_outcome_covariance_source = covariance_source,
    fuller_predictor_outcome_covariance_max_abs = max(abs(omega_xy_mat)),
    fuller_alpha_step1_requested = alpha_step1,
    fuller_alpha_step1_used = alpha_step1,
    fuller_alpha_step3_requested = alpha_step3,
    fuller_alpha_step3_used = alpha_step3,
    fuller_alpha_scaling_requested = alpha_scaling,
    fuller_alpha_scaling_used = alpha_scaling,
    fuller_auto_tempered = auto_tempered,
    fuller_sx1_star_condition = sx1_star_diag$condition_number,
    fuller_sx1_star_min_eigen = sx1_star_diag$min_eigen,
    fuller_sx1_observed_max_eigen = sx1_observed_diag$max_eigen,
    fuller_sx1_star_relative_min_eigen = sx1_star_relative_min_eigen,
    fuller_sx_star_condition = sx_star_diag$condition_number,
    fuller_sx_star_min_eigen = sx_star_diag$min_eigen,
    fuller_sx_observed_max_eigen = sx_observed_diag$max_eigen,
    fuller_sx_star_relative_min_eigen = sx_star_relative_min_eigen,
    fuller_scaling_condition = scaling_diag$condition_number,
    fuller_scaling_min_eigen = scaling_diag$min_eigen,
    fuller_scaling_observed_max_eigen = scaling_observed_diag$max_eigen,
    fuller_scaling_relative_min_eigen = scaling_relative_min_eigen,
    fuller_reference_se = NA_real_,
    fuller_se_ratio = NA_real_,
    fuller_auto_guard_pass = NA,
    fuller_auto_guard_reason = NA_character_,
    fuller_auto_guard_score = NA_real_,
    fuller_auto_full_weight_guard_pass = NA,
    fuller_auto_full_weight_guard_reason = NA_character_,
    fuller_auto_full_weight_se_ratio = NA_real_,
    fuller_auto_search_evaluations = NA_integer_,
    fuller_auto_search_nonmonotone = NA
  )
}

score_fuller_auto_candidates <- function(candidate_tbl,
                                         reference_se = NA_real_,
                                         target_condition_number = 1e5,
                                         min_sx_uncorr_eigen = sqrt(.Machine$double.eps),
                                         min_scaling_eigen = sqrt(.Machine$double.eps),
                                         min_sx_uncorr_relative_eigen = 5e-2,
                                         min_scaling_relative_eigen = 5e-2) {
  candidate_tbl <- dplyr::mutate(
    candidate_tbl,
    fuller_max_condition = fuller_row_max_finite(
      fuller_sx_star_condition,
      fuller_scaling_condition
    ),
    fuller_reference_se = reference_se,
    fuller_se_ratio = if (is.finite(reference_se) && reference_se > sqrt(.Machine$double.eps)) {
      se / reference_se
    } else {
      NA_real_
    }
  )

  min_sx_star_eigen <- min_sx_uncorr_eigen
  min_sx_star_relative_eigen <- min_sx_uncorr_relative_eigen
  sx_ok <- is.finite(candidate_tbl$fuller_sx_star_min_eigen) &
    candidate_tbl$fuller_sx_star_min_eigen >= min_sx_star_eigen
  scaling_ok <- is.finite(candidate_tbl$fuller_scaling_min_eigen) &
    candidate_tbl$fuller_scaling_min_eigen >= min_scaling_eigen
  sx_relative_ok <- is.finite(candidate_tbl$fuller_sx_star_relative_min_eigen) &
    candidate_tbl$fuller_sx_star_relative_min_eigen >= min_sx_star_relative_eigen
  scaling_relative_ok <- is.finite(candidate_tbl$fuller_scaling_relative_min_eigen) &
    candidate_tbl$fuller_scaling_relative_min_eigen >= min_scaling_relative_eigen
  condition_ok <- is.finite(candidate_tbl$fuller_max_condition) &
    candidate_tbl$fuller_max_condition <= target_condition_number
  se_ok <- is.finite(candidate_tbl$se) & candidate_tbl$se > 0

  status_ok <- as.integer(candidate_tbl$status_code) == 0L
  reason <- dplyr::case_when(
    !status_ok ~ paste0("status_", candidate_tbl$status_code),
    !se_ok ~ "bad_se",
    !sx_ok ~ "sx_star_eigen_floor",
    !scaling_ok ~ "scaling_eigen_floor",
    !sx_relative_ok ~ "sx_star_relative_eigen_floor",
    !scaling_relative_ok ~ "scaling_relative_eigen_floor",
    !condition_ok ~ "condition_cap",
    TRUE ~ "ok"
  )

  eig_penalty <- fuller_guard_penalty(candidate_tbl$fuller_sx_star_min_eigen, min_sx_star_eigen) +
    fuller_guard_penalty(candidate_tbl$fuller_scaling_min_eigen, min_scaling_eigen) +
    fuller_guard_penalty(candidate_tbl$fuller_sx_star_relative_min_eigen, min_sx_star_relative_eigen) +
    fuller_guard_penalty(candidate_tbl$fuller_scaling_relative_min_eigen, min_scaling_relative_eigen)
  condition_penalty <- ifelse(
    is.finite(candidate_tbl$fuller_max_condition) &
      is.finite(target_condition_number) &
      target_condition_number > 0,
    pmax(0, log(candidate_tbl$fuller_max_condition / target_condition_number)),
    Inf
  )
  status_penalty <- ifelse(status_ok & se_ok, 0, Inf)
  temper_penalty <- pmax(0, 1 - candidate_tbl$fuller_measurement_weight_used) * 0.01

  dplyr::mutate(
    candidate_tbl,
    fuller_auto_guard_pass = status_ok & se_ok & sx_ok & scaling_ok & sx_relative_ok &
      scaling_relative_ok & condition_ok,
    fuller_auto_guard_reason = reason,
    fuller_auto_guard_score = status_penalty + eig_penalty + condition_penalty + temper_penalty
  )
}

fuller_reference_dual_se <- function(stage2_df,
                                     outcome,
                                     predictor_u0,
                                     predictor_u1) {
  ref <- tryCatch(
    fit_observed_dual(
      stage2_df,
      outcome = outcome,
      predictor_u0 = predictor_u0,
      predictor_u1 = predictor_u1
    ),
    error = function(e) NULL
  )
  if (is.null(ref) || nrow(ref) == 0L || !"se" %in% names(ref)) {
    return(NA_real_)
  }

  hc3 <- ref$se[ref$se_type == "hc3" & as.integer(ref$status_code) == 0L]
  if (length(hc3) > 0L && is.finite(hc3[[1]]) && hc3[[1]] > sqrt(.Machine$double.eps)) {
    return(hc3[[1]])
  }

  naive <- ref$se[ref$se_type == "naive" & as.integer(ref$status_code) == 0L]
  if (length(naive) > 0L && is.finite(naive[[1]]) && naive[[1]] > sqrt(.Machine$double.eps)) {
    return(naive[[1]])
  }

  NA_real_
}

#' Fit the corrected, stabilized VH Fuller EIV estimator.
#'
#' @details
#' By default this uses the correct predictor--outcome covariance together with
#' the two finite-sample stabilizations selected for the VH primary method.
#' Set `preliminary_moment = "fuller"` and/or `variance_bread = "fuller"`, or
#' use [fit_fuller_dual_variants()], for the prespecified algebra comparisons.
#'
#' @inheritParams fit_fuller_dual_core
#'
#' @return A one-row tibble using the original Fuller result schema.
fit_fuller_dual <- function(stage2_df,
                            outcome,
                            predictor_u0 = NULL,
                            predictor_u1,
                            meas11 = NULL,
                            meas12 = NULL,
                            meas22,
                            outcome_meas_var = NULL,
                            predictor_outcome_meas_cov_u0 = NULL,
                            predictor_outcome_meas_cov_u1 = NULL,
                            preliminary_moment = c("modified", "fuller"),
                            variance_bread = c("modified", "fuller"),
                            skip_internal_scaling = TRUE) {
  preliminary_moment <- match.arg(preliminary_moment)
  variance_bread <- match.arg(variance_bread)
  out <- fit_fuller_dual_core(
    stage2_df,
    outcome = outcome,
    predictor_u0 = predictor_u0,
    predictor_u1 = predictor_u1,
    meas11 = meas11,
    meas12 = meas12,
    meas22 = meas22,
    outcome_meas_var = outcome_meas_var,
    predictor_outcome_meas_cov_u0 = predictor_outcome_meas_cov_u0,
    predictor_outcome_meas_cov_u1 = predictor_outcome_meas_cov_u1,
    measurement_weight = 1,
    preliminary_moment = preliminary_moment,
    variance_bread = variance_bread,
    auto_tempered = FALSE,
    skip_internal_scaling = skip_internal_scaling
  )
  dplyr::select(out, dplyr::any_of(fuller_dual_result_columns()))
}

#' Add explicit zero predictor--outcome error covariances for VH simulations.
#'
#' @details
#' The corrected Fuller equations allow a nonzero covariance between the
#' measurement error in each predictor score and the measurement error in the
#' Stage-2 outcome. VH Studies 1--4 currently generate these errors from
#' separate processes (or observe the outcome without measurement error), so
#' their covariance matrix is exactly zero. Keeping named zero columns in each
#' Stage-2 data set makes that design assumption auditable and prevents an
#' outcome variance from being mistaken for a predictor--outcome covariance.
#'
#' @param stage2_df Stage-2 data frame.
#' @param include_u0 Whether to add the intercept-like predictor covariance.
#'
#' @return `stage2_df` with the shared Fuller covariance columns added.
add_zero_fuller_predictor_outcome_covariance <- function(stage2_df,
                                                         include_u0 = TRUE) {
  if (isTRUE(include_u0)) {
    stage2_df$fuller_predictor_outcome_meas_cov_u0 <- 0
  }
  stage2_df$fuller_predictor_outcome_meas_cov_u1 <- 0
  stage2_df
}

#' Fit the prespecified Fuller algebra-comparison variants.
#'
#' @details
#' `"stabilized"` is the corrected VH primary estimator: it uses the
#' lambda/alpha-modified preliminary moment and modified variance bread.
#' `"fuller_preliminary"` replaces only the preliminary moment with Fuller's
#' fully subtracted moment. `"fuller_equations"` additionally uses Fuller's
#' fully subtracted variance bread. This nested sequence identifies whether the
#' preliminary or bread stabilization is responsible for taming an SE tail.
#' All variants use the correct predictor--outcome measurement-error covariance.
#'
#' @inheritParams fit_fuller_dual
#' @param variants Character vector containing any of `"stabilized"`,
#' `"fuller_preliminary"`, and `"fuller_equations"`.
#'
#' @return One row per requested variant, identified by `fuller_variant`.
fit_fuller_dual_variants <- function(
    stage2_df,
    outcome,
    predictor_u0 = NULL,
    predictor_u1,
    meas11 = NULL,
    meas12 = NULL,
    meas22,
    outcome_meas_var = NULL,
    predictor_outcome_meas_cov_u0 = NULL,
    predictor_outcome_meas_cov_u1 = NULL,
    variants = c("stabilized", "fuller_preliminary", "fuller_equations"),
    skip_internal_scaling = TRUE) {
  allowed <- c("stabilized", "fuller_preliminary", "fuller_equations")
  if (length(variants) == 0L || any(!variants %in% allowed)) {
    stop("`variants` must contain only: ", paste(allowed, collapse = ", "))
  }
  variants <- unique(variants)
  purrr::map_dfr(variants, function(variant) {
    # These are deliberately nested comparisons rather than a full 2 x 2
    # factorial. The second row changes only the preliminary moment relative to
    # the VH primary method; the third then also changes the variance bread.
    # This one-component-at-a-time sequence identifies which stabilization
    # controls the SE tail while retaining the book equations as an endpoint.
    settings <- switch(
      variant,
      stabilized = c(preliminary = "modified", bread = "modified"),
      fuller_preliminary = c(preliminary = "fuller", bread = "modified"),
      fuller_equations = c(preliminary = "fuller", bread = "fuller")
    )
    out <- fit_fuller_dual(
      stage2_df = stage2_df,
      outcome = outcome,
      predictor_u0 = predictor_u0,
      predictor_u1 = predictor_u1,
      meas11 = meas11,
      meas12 = meas12,
      meas22 = meas22,
      outcome_meas_var = outcome_meas_var,
      predictor_outcome_meas_cov_u0 = predictor_outcome_meas_cov_u0,
      predictor_outcome_meas_cov_u1 = predictor_outcome_meas_cov_u1,
      preliminary_moment = unname(settings[["preliminary"]]),
      variance_bread = unname(settings[["bread"]]),
      skip_internal_scaling = skip_internal_scaling
    )
    dplyr::mutate(out, fuller_variant = variant, .before = 1L)
  })
}

#' Convert an average BLUP measurement model to Fuller's additive-error form.
#'
#' @details
#' Lai's 2S-PAA replaces the row-specific BLUP measurement models
#'
#' `m_i = Lambda_i b_i + delta_i`, `Var(delta_i) = Theta_i`
#'
#' by a common model using the elementwise sample averages `Lambda_bar` and
#' `Theta_bar`. This helper applies the same averaging operation before
#' converting that common model to the identity-loading form required by the
#' Fuller estimator:
#'
#' `c_i = solve(Lambda_bar) m_i`
#'
#' `S_bar = solve(Lambda_bar) Theta_bar solve(Lambda_bar)'`.
#'
#' This operation is deliberately different from first applying each
#' row-specific inverse loading and then averaging the resulting score-error
#' covariance matrices.
#'
#' @param stage2_df Data frame containing the two BLUPs and the row-specific
#'   `lambda` and `theta` entries.
#' @param blup_u0,blup_u1 Column names for the intercept and slope BLUPs.
#'
#' @return A list containing `data`, the input data frame augmented with
#'   `fuller_average_u0`, `fuller_average_u1`, and the three common
#'   `fuller_average_meas*` columns, as well as `lambda_bar`, `theta_bar`, and
#'   `measurement_covariance_bar`.
prepare_fuller_average_measurement <- function(
    stage2_df,
    blup_u0 = "u0_eb",
    blup_u1 = "u1_eb") {
  lambda_cols <- c("lambda11", "lambda12", "lambda21", "lambda22")
  theta_cols <- c("theta11", "theta12", "theta22")
  required_cols <- c(blup_u0, blup_u1, lambda_cols, theta_cols)
  missing_cols <- setdiff(required_cols, names(stage2_df))
  if (length(missing_cols) > 0L) {
    stop(
      "Average-measurement Fuller inputs are missing columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  lambda_means <- colMeans(stage2_df[, lambda_cols, drop = FALSE], na.rm = TRUE)
  theta_means <- colMeans(stage2_df[, theta_cols, drop = FALSE], na.rm = TRUE)
  if (any(!is.finite(lambda_means)) || any(!is.finite(theta_means))) {
    stop("Average-measurement Fuller requires finite mean lambda and theta entries.")
  }

  lambda_bar <- matrix(
    lambda_means[c("lambda11", "lambda12", "lambda21", "lambda22")],
    nrow = 2L,
    byrow = TRUE,
    dimnames = list(c("u0_eb", "u1_eb"), c("u0", "u1"))
  )
  theta_bar <- matrix(
    c(
      theta_means[["theta11"]], theta_means[["theta12"]],
      theta_means[["theta12"]], theta_means[["theta22"]]
    ),
    nrow = 2L,
    byrow = TRUE,
    dimnames = list(c("u0_eb", "u1_eb"), c("u0_eb", "u1_eb"))
  )
  lambda_bar_inv <- tryCatch(
    solve(lambda_bar),
    error = function(e) {
      stop("Mean BLUP loading matrix is not invertible: ", conditionMessage(e))
    }
  )

  blup_matrix <- as.matrix(stage2_df[, c(blup_u0, blup_u1), drop = FALSE])
  storage.mode(blup_matrix) <- "double"
  transformed_scores <- blup_matrix %*% t(lambda_bar_inv)
  measurement_covariance_bar <-
    lambda_bar_inv %*% theta_bar %*% t(lambda_bar_inv)
  measurement_covariance_bar <-
    (measurement_covariance_bar + t(measurement_covariance_bar)) / 2
  if (any(!is.finite(measurement_covariance_bar))) {
    stop("Average-measurement Fuller produced a nonfinite error covariance matrix.")
  }

  out <- stage2_df
  out$fuller_average_u0 <- transformed_scores[, 1L]
  out$fuller_average_u1 <- transformed_scores[, 2L]
  out$fuller_average_meas11 <- measurement_covariance_bar[1L, 1L]
  out$fuller_average_meas12 <- measurement_covariance_bar[1L, 2L]
  out$fuller_average_meas22 <- measurement_covariance_bar[2L, 2L]

  list(
    data = out,
    lambda_bar = lambda_bar,
    theta_bar = theta_bar,
    measurement_covariance_bar = measurement_covariance_bar
  )
}

#' Fit the Fuller analogue of Lai's average-measurement 2S-PAA estimator.
#'
#' @inheritParams prepare_fuller_average_measurement
#' @param outcome Column name for the observed Stage-2 outcome.
#' @param predictor_outcome_meas_cov_u0,predictor_outcome_meas_cov_u1 Optional
#' columns containing the outcome covariance with the intercept- and
#' slope-like BLUP measurement errors. The wrapper premultiplies this
#' covariance vector by `solve(Lambda_bar)` so it is on the transformed Fuller
#' predictor scale. Omitted columns are zero.
#' @param skip_internal_scaling Passed to [fit_fuller_dual()].
#'
#' @return A one-row tibble using the traditional Fuller result schema.
fit_fuller_average_measurement <- function(
    stage2_df,
    outcome,
    blup_u0 = "u0_eb",
    blup_u1 = "u1_eb",
    outcome_meas_var = NULL,
    predictor_outcome_meas_cov_u0 = NULL,
    predictor_outcome_meas_cov_u1 = NULL,
    skip_internal_scaling = TRUE) {
  tryCatch(
    {
      prepared <- prepare_fuller_average_measurement(
        stage2_df,
        blup_u0 = blup_u0,
        blup_u1 = blup_u1
      )
      average_cov_u0 <- NULL
      average_cov_u1 <- NULL
      if (!is.null(predictor_outcome_meas_cov_u0) ||
        !is.null(predictor_outcome_meas_cov_u1)) {
        raw_covariance <- cbind(
          if (!is.null(predictor_outcome_meas_cov_u0)) {
            stage2_df[[predictor_outcome_meas_cov_u0]]
          } else {
            0
          },
          if (!is.null(predictor_outcome_meas_cov_u1)) {
            stage2_df[[predictor_outcome_meas_cov_u1]]
          } else {
            0
          }
        )
        transformed_covariance <- raw_covariance %*%
          t(solve(prepared$lambda_bar))
        average_cov_u0 <- "fuller_average_predictor_outcome_meas_cov_u0"
        average_cov_u1 <- "fuller_average_predictor_outcome_meas_cov_u1"
        prepared$data[[average_cov_u0]] <- transformed_covariance[, 1L]
        prepared$data[[average_cov_u1]] <- transformed_covariance[, 2L]
      }
      fit_fuller_dual(
        prepared$data,
        outcome = outcome,
        predictor_u0 = "fuller_average_u0",
        predictor_u1 = "fuller_average_u1",
        meas11 = "fuller_average_meas11",
        meas12 = "fuller_average_meas12",
        meas22 = "fuller_average_meas22",
        outcome_meas_var = outcome_meas_var,
        predictor_outcome_meas_cov_u0 = average_cov_u0,
        predictor_outcome_meas_cov_u1 = average_cov_u1,
        skip_internal_scaling = skip_internal_scaling
      )
    },
    error = function(e) {
      fallback <- stage2_df
      fallback$fuller_average_u0 <- NA_real_
      fallback$fuller_average_u1 <- NA_real_
      fallback$fuller_average_meas11 <- NA_real_
      fallback$fuller_average_meas12 <- NA_real_
      fallback$fuller_average_meas22 <- NA_real_
      failed_fit <- fit_fuller_dual(
        fallback,
        outcome = outcome,
        predictor_u0 = "fuller_average_u0",
        predictor_u1 = "fuller_average_u1",
        meas11 = "fuller_average_meas11",
        meas12 = "fuller_average_meas12",
        meas22 = "fuller_average_meas22",
        outcome_meas_var = outcome_meas_var,
        predictor_outcome_meas_cov_u0 = predictor_outcome_meas_cov_u0,
        predictor_outcome_meas_cov_u1 = predictor_outcome_meas_cov_u1,
        skip_internal_scaling = skip_internal_scaling
      )
      dplyr::mutate(
        failed_fit,
        status_code = 3L,
        mx_issue_class = "fuller_average_measurement_transform_failed",
        mx_issue_detail = substr(conditionMessage(e), 1L, 500L)
      )
    }
  )
}

#' Fit a Fuller EIV estimator with internal stepdown tempering.
#'
#' @details
#' The stepdown procedure tempers the measurement-error correction when the 
#' corrected predictor covariance matrices approach singularity. It performs a 
#' grid search over a set of `candidate_weights` (ranging from 1 down to 0). 
#' If the full correction (weight = 1) leads to an ill-conditioned matrix 
#' (assessed via minimum eigenvalues relative to the observed data),
#' the algorithm bisects the interval between the largest admissible weight and 
#' the lowest inadmissible weight to find the near-optimal measurement-error
#' correction scaling factor. This acts as a regularizer.
#'
#' @inheritParams fit_fuller_dual_core
#' @param candidate_weights Optional numeric vector of measurement-error
#' correction weights in `[0, 1]` for the initial coarse search. If `NULL`, a
#' uniform grid from 1 to 0 is generated.
#' @param coarse_grid_size Number of grid points used when `candidate_weights`
#' is `NULL`.
#' @param search_tolerance Absolute alpha tolerance for bisection refinement.
#' @param max_refinements Maximum bisection iterations after the coarse search.
#' @param target_condition_number Maximum accepted condition number across
#' Fuller internal corrected/scaling blocks.
#' @param min_sx_uncorr_eigen Absolute minimum eigenvalue floor for the
#' corrected S* predictor block used in the Fuller variance calculation.
#' @param min_scaling_eigen Absolute minimum eigenvalue floor for the corrected
#' scaling covariance.
#' @param min_sx_uncorr_relative_eigen Minimum ratio of the corrected S*
#' predictor block's weakest eigenvalue to the observed block's largest eigenvalue.
#' @param min_scaling_relative_eigen Minimum ratio of the corrected scaling
#' block's weakest eigenvalue to the observed scaling block's largest
#' eigenvalue.
#'
#' @return A one-row tibble with the Fuller core result plus stepdown
#' diagnostics, including the selected measurement-error correction weight.
fit_fuller_dual_stepdown <- function(stage2_df,
                                     outcome,
                                     predictor_u0,
                                     predictor_u1,
                                     meas11,
                                     meas12,
                                     meas22,
                                     outcome_meas_var = NULL,
                                     predictor_outcome_meas_cov_u0 = NULL,
                                     predictor_outcome_meas_cov_u1 = NULL,
                                     candidate_weights = NULL,
                                     coarse_grid_size = 9L,
                                     search_tolerance = 0.005,
                                     max_refinements = 12L,
                                     target_condition_number = 1e5,
                                     min_sx_uncorr_eigen = sqrt(.Machine$double.eps),
                                     min_scaling_eigen = sqrt(.Machine$double.eps),
                                     min_sx_uncorr_relative_eigen = 5e-2,
                                     min_scaling_relative_eigen = 5e-2,
                                     skip_internal_scaling = TRUE) {
  if (is.null(candidate_weights)) {
    coarse_grid_size <- max(3L, as.integer(coarse_grid_size))
    candidate_weights <- seq(1, 0, length.out = coarse_grid_size)
  }
  candidate_weights <- sort(unique(pmin(1, pmax(0, candidate_weights))), decreasing = TRUE)
  search_tolerance <- max(sqrt(.Machine$double.eps), search_tolerance)
  max_refinements <- max(0L, as.integer(max_refinements))

  reference_se <- fuller_reference_dual_se(
    stage2_df,
    outcome = outcome,
    predictor_u0 = predictor_u0,
    predictor_u1 = predictor_u1
  )

  evaluate_weights <- function(weights) {
    weights <- sort(unique(pmin(1, pmax(0, weights))), decreasing = TRUE)
    candidate_results <- lapply(weights, function(weight) {
      fit_fuller_dual_core(
        stage2_df,
        outcome = outcome,
        predictor_u0 = predictor_u0,
        predictor_u1 = predictor_u1,
        meas11 = meas11,
        meas12 = meas12,
        meas22 = meas22,
        outcome_meas_var = outcome_meas_var,
        predictor_outcome_meas_cov_u0 = predictor_outcome_meas_cov_u0,
        predictor_outcome_meas_cov_u1 = predictor_outcome_meas_cov_u1,
        measurement_weight = weight,
        auto_tempered = TRUE, 
        skip_internal_scaling = skip_internal_scaling
      )
    })
    score_fuller_auto_candidates(
      dplyr::bind_rows(candidate_results),
      reference_se = reference_se,
      target_condition_number = target_condition_number,
      min_sx_uncorr_eigen = min_sx_uncorr_eigen,
      min_scaling_eigen = min_scaling_eigen,
      min_sx_uncorr_relative_eigen = min_sx_uncorr_relative_eigen,
      min_scaling_relative_eigen = min_scaling_relative_eigen
    )
  }

  deduplicate_candidates <- function(tbl) {
    tbl <- dplyr::mutate(tbl, .weight_key = round(fuller_measurement_weight_used, 10L))
    tbl <- dplyr::arrange(tbl, dplyr::desc(fuller_measurement_weight_used), fuller_auto_guard_score)
    tbl <- dplyr::group_by(tbl, .weight_key)
    tbl <- dplyr::slice(tbl, 1L)
    tbl <- dplyr::ungroup(tbl)
    dplyr::select(tbl, -.weight_key)
  }

  auto_nonmonotone <- function(tbl) {
    check <- dplyr::arrange(tbl, dplyr::desc(fuller_measurement_weight_used))
    check <- dplyr::mutate(check, pass = as.logical(fuller_auto_guard_pass))
    pass <- ifelse(is.na(check$pass), FALSE, check$pass)
    any(cummax(pass) > 0 & !pass)
  }

  candidate_tbl <- evaluate_weights(candidate_weights)
  full_weight <- max(candidate_weights, na.rm = TRUE)
  full_candidate <- dplyr::filter(candidate_tbl, fuller_measurement_weight_used == full_weight)
  full_candidate <- dplyr::slice(full_candidate, 1L)
  full_guard_pass <- if (nrow(full_candidate) > 0L) full_candidate$fuller_auto_guard_pass[[1]] else NA
  full_guard_reason <- if (nrow(full_candidate) > 0L) full_candidate$fuller_auto_guard_reason[[1]] else NA_character_
  full_se_ratio <- if (nrow(full_candidate) > 0L) full_candidate$fuller_se_ratio[[1]] else NA_real_

  # Only attempt bisection if the full weight failed but at least one admissible weight exists
  needs_refinement <- !isTRUE(full_guard_pass) && any(candidate_tbl$fuller_auto_guard_pass)
  
  if (needs_refinement) {
    admissible <- dplyr::filter(candidate_tbl, fuller_auto_guard_pass)
    lower <- max(admissible$fuller_measurement_weight_used, na.rm = TRUE)
    
    upper_candidates <- dplyr::filter(
      candidate_tbl,
      !fuller_auto_guard_pass,
      fuller_measurement_weight_used > lower
    )
    
    if (nrow(upper_candidates) > 0L) {
      upper <- min(upper_candidates$fuller_measurement_weight_used, na.rm = TRUE)
      
      for (i in seq_len(max_refinements)) {
        if (!is.finite(upper - lower) || (upper - lower) <= search_tolerance) {
          break
        }
        midpoint <- (lower + upper) / 2
        midpoint_tbl <- evaluate_weights(midpoint)
        candidate_tbl <- deduplicate_candidates(dplyr::bind_rows(candidate_tbl, midpoint_tbl))
        midpoint_pass <- isTRUE(midpoint_tbl$fuller_auto_guard_pass[[1]])
        if (midpoint_pass) {
          lower <- midpoint
        } else {
          upper <- midpoint
        }
      }
    }
  }

  candidate_tbl <- deduplicate_candidates(candidate_tbl)
  admissible <- dplyr::filter(candidate_tbl, fuller_auto_guard_pass)
  admissible <- dplyr::arrange(admissible, dplyr::desc(fuller_measurement_weight_used))
  if (nrow(admissible) == 0L) {
    successful <- dplyr::filter(candidate_tbl, as.integer(status_code) == 0L, is.finite(fuller_auto_guard_score))
    successful <- dplyr::arrange(successful, fuller_auto_guard_score, dplyr::desc(fuller_measurement_weight_used))
    chosen_candidate <- if (nrow(successful) > 0L) {
      dplyr::slice(successful, 1L)
    } else {
      evaluate_weights(candidate_weights[[length(candidate_weights)]])
    }
    chosen_candidate$fuller_auto_guard_reason <- paste0("fallback_", chosen_candidate$fuller_auto_guard_reason)
    chosen_guard_pass <- FALSE
  } else {
    chosen_candidate <- dplyr::slice(admissible, 1L)
    chosen_guard_pass <- TRUE
  }

  chosen_weight <- chosen_candidate$fuller_measurement_weight_used[[1]]
  chosen_reason <- chosen_candidate$fuller_auto_guard_reason[[1]]
  chosen_score <- chosen_candidate$fuller_auto_guard_score[[1]]
  chosen_se_ratio <- chosen_candidate$fuller_se_ratio[[1]]
  search_evaluations <- nrow(candidate_tbl)
  search_nonmonotone <- auto_nonmonotone(candidate_tbl)

  out <- fit_fuller_dual_core(
    stage2_df,
    outcome = outcome,
    predictor_u0 = predictor_u0,
    predictor_u1 = predictor_u1,
    meas11 = meas11,
    meas12 = meas12,
    meas22 = meas22,
    outcome_meas_var = outcome_meas_var,
    predictor_outcome_meas_cov_u0 = predictor_outcome_meas_cov_u0,
    predictor_outcome_meas_cov_u1 = predictor_outcome_meas_cov_u1,
    measurement_weight = chosen_weight,
    auto_tempered = TRUE, 
    skip_internal_scaling = skip_internal_scaling
  )
  actual_se_ratio <- if (is.finite(reference_se) && reference_se > sqrt(.Machine$double.eps) &&
    nrow(out) > 0L && is.finite(out$se[[1]])) {
    out$se[[1]] / reference_se
  } else {
    chosen_se_ratio
  }

  dplyr::mutate(
    out,
    fuller_reference_se = reference_se,
    fuller_se_ratio = actual_se_ratio,
    fuller_auto_guard_pass = chosen_guard_pass,
    fuller_auto_guard_reason = chosen_reason,
    fuller_auto_guard_score = chosen_score,
    fuller_auto_full_weight_guard_pass = full_guard_pass,
    fuller_auto_full_weight_guard_reason = full_guard_reason,
    fuller_auto_full_weight_se_ratio = full_se_ratio,
    fuller_auto_search_evaluations = search_evaluations,
    fuller_auto_search_nonmonotone = search_nonmonotone,
    mx_issue_detail = ifelse(
      is.na(mx_issue_detail) | mx_issue_detail == "ok",
      sprintf(
        paste0(
          "auto_weight=%0.4f; guard=%s; target_condition=%0.3e; ",
          "sx_rel_floor=%0.3e; scaling_rel_floor=%0.3e; ",
          "evals=%d; nonmonotone=%s"
        ),
        chosen_weight,
        chosen_reason,
        target_condition_number,
        min_sx_uncorr_relative_eigen,
        min_scaling_relative_eigen,
        search_evaluations,
        search_nonmonotone
      ),
      mx_issue_detail
    )
  )
}

#' Fit a Fuller EIV estimator with alpha stepdown tempering.
#'
#' @details
#' The alpha stepdown procedure searches separately for Step-1 and Step-3
#' correction factors (alpha). It uses a coarse grid followed by bisection to
#' find the smallest alpha that preserves a minimum ratio of corrected to
#' observed predictor-block eigenvalues. Step 1 guards `S1*`; Step 3 guards
#' the final corrected `S*` matrix. Because larger alpha values temper more
#' of the Fuller correction, this keeps the strongest admissible correction.
#'
#' @inheritParams fit_fuller_dual_core
#' @param candidate_alphas Optional numeric vector of alpha candidates for the
#' coarse search. If `NULL`, a uniform grid between the admissible upper and
#' lower bounds is generated.
#' @param coarse_grid_size Number of grid points used when `candidate_alphas`
#' is `NULL`.
#' @param search_tolerance Absolute alpha tolerance for bisection refinement.
#' @param max_refinements Maximum bisection iterations after the coarse search.
#' @param min_sx1_relative_eigen Minimum ratio of the Step-1 corrected
#' predictor-block minimum eigenvalue to the observed block maximum eigenvalue.
#' @param min_sx_star_relative_eigen Minimum ratio of the Step-3 corrected
#' predictor-block minimum eigenvalue to the observed block maximum eigenvalue.
#' @param target_condition_number Maximum accepted condition number for the
#' Step-1 and Step-3 corrected predictor blocks.
#'
#' @return A one-row tibble with the Fuller core result and selected alphas.
fit_fuller_dual_alpha_stepdown <- function(stage2_df,
                                           outcome,
                                           predictor_u0,
                                           predictor_u1,
                                           meas11,
                                           meas12,
                                           meas22,
                                           outcome_meas_var = NULL,
                                           predictor_outcome_meas_cov_u0 = NULL,
                                           predictor_outcome_meas_cov_u1 = NULL,
                                           candidate_alphas = NULL,
                                           coarse_grid_size = 9L,
                                           search_tolerance = 0.005,
                                           max_refinements = 12L,
                                           min_sx1_star_relative_eigen = 5e-2,
                                           min_sx_star_relative_eigen = 5e-2,
                                           min_scaling_relative_eigen = 5e-2,
                                           min_sx1_star_eigen = sqrt(.Machine$double.eps),
                                           min_sx_star_eigen = sqrt(.Machine$double.eps),
                                           min_scaling_eigen = sqrt(.Machine$double.eps),
                                           target_condition_number = 1e5,
                                           skip_internal_scaling = TRUE) {
  if (is.null(candidate_alphas)) {
    coarse_grid_size <- max(3L, as.integer(coarse_grid_size))
  }
  search_tolerance <- max(sqrt(.Machine$double.eps), search_tolerance)
  max_refinements <- max(0L, as.integer(max_refinements))

  has_u0 <- !is.null(predictor_u0)
  p <- if (has_u0) 3L else 2L
  cols_needed <- c(outcome, predictor_u1, meas22)
  if (has_u0) {
    cols_needed <- c(outcome, predictor_u0, predictor_u1, meas11, meas12, meas22)
  }
  if (!is.null(outcome_meas_var)) {
    cols_needed <- c(cols_needed, outcome_meas_var)
  }
  if (!is.null(predictor_outcome_meas_cov_u0)) {
    cols_needed <- c(cols_needed, predictor_outcome_meas_cov_u0)
  }
  if (!is.null(predictor_outcome_meas_cov_u1)) {
    cols_needed <- c(cols_needed, predictor_outcome_meas_cov_u1)
  }
  cols_needed <- unique(cols_needed)
  dat <- stage2_df[, cols_needed, drop = FALSE]
  dat <- dat[stats::complete.cases(dat), , drop = FALSE]
  m <- nrow(dat)
  alpha_lower <- p + 1

  alpha_candidates <- function(alpha_upper) {
    if (!is.finite(alpha_upper) || alpha_upper <= alpha_lower) {
      return(alpha_lower)
    }
    cand <- if (is.null(candidate_alphas)) {
      seq(alpha_lower, alpha_upper, length.out = coarse_grid_size)
    } else {
      candidate_alphas
    }
    cand <- cand[is.finite(cand)]
    if (length(cand) == 0L) {
      cand <- alpha_lower
    }
    cand <- sort(unique(pmin(alpha_upper, pmax(alpha_lower, cand))), decreasing = FALSE)
    if (!(alpha_lower %in% cand)) {
      cand <- c(alpha_lower, cand)
    }
    if (is.null(candidate_alphas) && !(alpha_upper %in% cand)) {
      cand <- c(cand, alpha_upper)
    }
    sort(unique(cand), decreasing = FALSE)
  }

  refine_alpha_choice <- function(candidate_tbl,
                                  alpha_col,
                                  pass_col,
                                  evaluate_fn,
                                  deduplicate_fn) {
    candidate_tbl <- deduplicate_fn(candidate_tbl)
    pass <- as.logical(candidate_tbl[[pass_col]])
    pass <- ifelse(is.na(pass), FALSE, pass)
    alpha_values <- candidate_tbl[[alpha_col]]

    if (!any(pass)) {
      return(list(alpha = alpha_lower, candidates = candidate_tbl))
    }

    upper <- min(alpha_values[pass], na.rm = TRUE)
    lower_candidates <- alpha_values[!pass & alpha_values < upper]

    if (length(lower_candidates) > 0L) {
      lower <- max(lower_candidates, na.rm = TRUE)
      for (i in seq_len(max_refinements)) {
        if (!is.finite(upper - lower) || (upper - lower) <= search_tolerance) {
          break
        }
        midpoint <- (lower + upper) / 2
        midpoint_tbl <- evaluate_fn(midpoint)
        candidate_tbl <- deduplicate_fn(dplyr::bind_rows(candidate_tbl, midpoint_tbl))
        midpoint_pass <- isTRUE(midpoint_tbl[[pass_col]][[1]])
        midpoint_alpha <- midpoint_tbl[[alpha_col]][[1]]
        if (midpoint_pass) {
          upper <- midpoint_alpha
        } else {
          lower <- midpoint_alpha
        }
      }
    }

    candidate_tbl <- deduplicate_fn(candidate_tbl)
    pass <- as.logical(candidate_tbl[[pass_col]])
    pass <- ifelse(is.na(pass), FALSE, pass)
    alpha_values <- candidate_tbl[[alpha_col]]
    list(alpha = min(alpha_values[pass], na.rm = TRUE), candidates = candidate_tbl)
  }

  alpha_upper_from_lambda <- function(lambda_hat) {
    if (!is.finite(lambda_hat) || !is.finite(m) || m <= 0) {
      return(alpha_lower)
    }
    # Larger alpha tempers more of the Fuller measurement-error subtraction.
    # When lambda exceeds the standard Fuller branch threshold, alpha values
    # above m make c = 1 - alpha / m negative, which flips the correction into
    # adding measurement-error covariance. Cap at m so stepdown cannot
    # over-temper past a zero correction.
    if (lambda_hat <= 1 + 1 / m) {
      return(lambda_hat * m - 1)
    }
    m
  }

  score_alpha_step1 <- function(tbl) {
    sx1_star_ok <- is.finite(tbl$fuller_sx1_star_min_eigen) &
      tbl$fuller_sx1_star_min_eigen >= min_sx1_star_eigen
    sx1_star_relative_ok <- is.finite(tbl$fuller_sx1_star_relative_min_eigen) &
      tbl$fuller_sx1_star_relative_min_eigen >= min_sx1_star_relative_eigen
    condition_ok <- is.finite(tbl$fuller_sx1_star_condition) &
      tbl$fuller_sx1_star_condition <= target_condition_number
    status_ok <- as.integer(tbl$status_code) == 0L
    reason <- dplyr::case_when(
      !status_ok ~ paste0("status_", tbl$status_code),
      !sx1_star_ok ~ "sx1_star_eigen_floor",
      !sx1_star_relative_ok ~ "sx1_star_relative_eigen_floor",
      !condition_ok ~ "sx1_star_condition_cap",
      TRUE ~ "ok"
    )
    eig_penalty <- fuller_guard_penalty(tbl$fuller_sx1_star_min_eigen, min_sx1_star_eigen) +
      fuller_guard_penalty(tbl$fuller_sx1_star_relative_min_eigen, min_sx1_star_relative_eigen)
    condition_penalty <- ifelse(
      is.finite(tbl$fuller_sx1_star_condition) &
        is.finite(target_condition_number) &
        target_condition_number > 0,
      pmax(0, log(tbl$fuller_sx1_star_condition / target_condition_number)),
      Inf
    )
    status_penalty <- ifelse(status_ok, 0, Inf)
    dplyr::mutate(
      tbl,
      alpha_step1_guard_pass = status_ok & sx1_star_ok & sx1_star_relative_ok & condition_ok,
      alpha_step1_guard_reason = reason,
      alpha_step1_guard_score = status_penalty + eig_penalty + condition_penalty
    )
  }

  deduplicate_alpha_step1 <- function(tbl) {
    tbl <- dplyr::mutate(tbl, .alpha_key = round(fuller_alpha_step1_used, 10L))
    tbl <- dplyr::arrange(tbl, fuller_alpha_step1_used, dplyr::desc(alpha_step1_guard_pass))
    tbl <- dplyr::group_by(tbl, .alpha_key)
    tbl <- dplyr::slice(tbl, 1L)
    tbl <- dplyr::ungroup(tbl)
    dplyr::select(tbl, -.alpha_key)
  }

  evaluate_alpha_step1 <- function(alphas) {
    alphas <- sort(unique(alphas), decreasing = FALSE)
    candidate_results <- lapply(alphas, function(alpha) {
      fit_fuller_dual_core(
        stage2_df,
        outcome = outcome,
        predictor_u0 = predictor_u0,
        predictor_u1 = predictor_u1,
        meas11 = meas11,
        meas12 = meas12,
        meas22 = meas22,
        outcome_meas_var = outcome_meas_var,
        predictor_outcome_meas_cov_u0 = predictor_outcome_meas_cov_u0,
        predictor_outcome_meas_cov_u1 = predictor_outcome_meas_cov_u1,
        measurement_weight = 1,
        alpha_step1 = alpha,
        alpha_step3 = alpha_lower,
        auto_tempered = TRUE,
        skip_internal_scaling = skip_internal_scaling
      )
    })
    score_alpha_step1(dplyr::bind_rows(candidate_results))
  }

  base_step1 <- fit_fuller_dual_core(
    stage2_df,
    outcome = outcome,
    predictor_u0 = predictor_u0,
    predictor_u1 = predictor_u1,
    meas11 = meas11,
    meas12 = meas12,
    meas22 = meas22,
    outcome_meas_var = outcome_meas_var,
    predictor_outcome_meas_cov_u0 = predictor_outcome_meas_cov_u0,
    predictor_outcome_meas_cov_u1 = predictor_outcome_meas_cov_u1,
    measurement_weight = 1,
    alpha_step1 = alpha_lower,
    alpha_step3 = alpha_lower,
    auto_tempered = TRUE,
    skip_internal_scaling = skip_internal_scaling
  )
  lambda1_hat <- if (nrow(base_step1) > 0L) base_step1$fuller_lambda1[[1]] else NA_real_
  alpha_step1_upper <- alpha_upper_from_lambda(lambda1_hat)
  if (!is.finite(alpha_step1_upper)) {
    alpha_step1_upper <- alpha_lower
  }
  chosen_alpha_step1 <- alpha_lower
  if (alpha_step1_upper > alpha_lower) {
    step1_choice <- refine_alpha_choice(
      evaluate_alpha_step1(alpha_candidates(alpha_step1_upper)),
      alpha_col = "fuller_alpha_step1_used",
      pass_col = "alpha_step1_guard_pass",
      evaluate_fn = evaluate_alpha_step1,
      deduplicate_fn = deduplicate_alpha_step1
    )
    chosen_alpha_step1 <- step1_choice$alpha
  }

  score_alpha_step3 <- function(tbl) {
    sx_star_ok <- is.finite(tbl$fuller_sx_star_min_eigen) &
      tbl$fuller_sx_star_min_eigen >= min_sx_star_eigen
    sx_star_relative_ok <- is.finite(tbl$fuller_sx_star_relative_min_eigen) &
      tbl$fuller_sx_star_relative_min_eigen >= min_sx_star_relative_eigen
    condition_ok <- is.finite(tbl$fuller_sx_star_condition) &
      tbl$fuller_sx_star_condition <= target_condition_number
    status_ok <- as.integer(tbl$status_code) == 0L
    reason <- dplyr::case_when(
      !status_ok ~ paste0("status_", tbl$status_code),
      !sx_star_ok ~ "sx_star_eigen_floor",
      !sx_star_relative_ok ~ "sx_star_relative_eigen_floor",
      !condition_ok ~ "sx_star_condition_cap",
      TRUE ~ "ok"
    )
    eig_penalty <- fuller_guard_penalty(tbl$fuller_sx_star_min_eigen, min_sx_star_eigen) +
      fuller_guard_penalty(tbl$fuller_sx_star_relative_min_eigen, min_sx_star_relative_eigen)
    condition_penalty <- ifelse(
      is.finite(tbl$fuller_sx_star_condition) &
        is.finite(target_condition_number) &
        target_condition_number > 0,
      pmax(0, log(tbl$fuller_sx_star_condition / target_condition_number)),
      Inf
    )
    status_penalty <- ifelse(status_ok, 0, Inf)
    dplyr::mutate(
      tbl,
      alpha_step3_guard_pass = status_ok & sx_star_ok & sx_star_relative_ok & condition_ok,
      alpha_step3_guard_reason = reason,
      alpha_step3_guard_score = status_penalty + eig_penalty + condition_penalty
    )
  }

  deduplicate_alpha_step3 <- function(tbl) {
    tbl <- dplyr::mutate(tbl, .alpha_key = round(fuller_alpha_step3_used, 10L))
    tbl <- dplyr::arrange(tbl, fuller_alpha_step3_used, dplyr::desc(alpha_step3_guard_pass))
    tbl <- dplyr::group_by(tbl, .alpha_key)
    tbl <- dplyr::slice(tbl, 1L)
    tbl <- dplyr::ungroup(tbl)
    dplyr::select(tbl, -.alpha_key)
  }

  evaluate_alpha_step3 <- function(alphas) {
    alphas <- sort(unique(alphas), decreasing = FALSE)
    candidate_results <- lapply(alphas, function(alpha) {
      fit_fuller_dual_core(
        stage2_df,
        outcome = outcome,
        predictor_u0 = predictor_u0,
        predictor_u1 = predictor_u1,
        meas11 = meas11,
        meas12 = meas12,
        meas22 = meas22,
        outcome_meas_var = outcome_meas_var,
        predictor_outcome_meas_cov_u0 = predictor_outcome_meas_cov_u0,
        predictor_outcome_meas_cov_u1 = predictor_outcome_meas_cov_u1,
        measurement_weight = 1,
        alpha_step1 = chosen_alpha_step1,
        alpha_step3 = alpha,
        auto_tempered = TRUE,
        skip_internal_scaling = skip_internal_scaling
      )
    })
    score_alpha_step3(dplyr::bind_rows(candidate_results))
  }

  base_step3 <- fit_fuller_dual_core(
    stage2_df,
    outcome = outcome,
    predictor_u0 = predictor_u0,
    predictor_u1 = predictor_u1,
    meas11 = meas11,
    meas12 = meas12,
    meas22 = meas22,
    outcome_meas_var = outcome_meas_var,
    predictor_outcome_meas_cov_u0 = predictor_outcome_meas_cov_u0,
    predictor_outcome_meas_cov_u1 = predictor_outcome_meas_cov_u1,
    measurement_weight = 1,
    alpha_step1 = chosen_alpha_step1,
    alpha_step3 = alpha_lower,
    auto_tempered = TRUE,
    skip_internal_scaling = skip_internal_scaling
  )
  lambda2_hat <- if (nrow(base_step3) > 0L) base_step3$fuller_lambda2[[1]] else NA_real_
  alpha_step3_upper <- alpha_upper_from_lambda(lambda2_hat)
  if (!is.finite(alpha_step3_upper)) {
    alpha_step3_upper <- alpha_lower
  }
  chosen_alpha_step3 <- alpha_lower
  if (alpha_step3_upper > alpha_lower) {
    step3_choice <- refine_alpha_choice(
      evaluate_alpha_step3(alpha_candidates(alpha_step3_upper)),
      alpha_col = "fuller_alpha_step3_used",
      pass_col = "alpha_step3_guard_pass",
      evaluate_fn = evaluate_alpha_step3,
      deduplicate_fn = deduplicate_alpha_step3
    )
    chosen_alpha_step3 <- step3_choice$alpha
  }
  
  if (!skip_internal_scaling) {
    score_alpha_scaling <- function(tbl) {
      status_ok <- as.integer(tbl$status_code) == 0L
      scaling_ok <- is.finite(tbl$fuller_scaling_min_eigen) &
        tbl$fuller_scaling_min_eigen >= min_scaling_eigen
      scaling_relative_ok <- is.finite(tbl$fuller_scaling_relative_min_eigen) &
        tbl$fuller_scaling_relative_min_eigen >= min_scaling_relative_eigen
      reason <- dplyr::case_when(
        !status_ok ~ paste0("status_", tbl$status_code),
        !scaling_ok ~ "scaling_eigen_floor",
        !scaling_relative_ok ~ "scaling_relative_eigen_floor",
        TRUE ~ "ok"
      )
      eig_penalty <- fuller_guard_penalty(tbl$fuller_scaling_min_eigen, min_scaling_eigen) +
        fuller_guard_penalty(tbl$fuller_scaling_relative_min_eigen, min_scaling_relative_eigen)
      condition_penalty <- ifelse(
        is.finite(tbl$fuller_scaling_condition) &
          is.finite(target_condition_number) &
          target_condition_number > 0,
        pmax(0, log(tbl$fuller_scaling_condition / target_condition_number)),
        Inf
      )
      status_penalty <- ifelse(status_ok, 0, Inf)
      dplyr::mutate(
        tbl,
        alpha_scaling_guard_pass = status_ok & scaling_ok & scaling_relative_ok,
        alpha_scaling_guard_reason = reason,
        alpha_scaling_guard_score = status_penalty + eig_penalty + condition_penalty
      )
    }
  
    deduplicate_alpha_scaling <- function(tbl) {
      tbl <- dplyr::mutate(tbl, .alpha_key = round(fuller_alpha_scaling_used, 10L))
      tbl <- dplyr::arrange(tbl, fuller_alpha_scaling_used, dplyr::desc(alpha_scaling_guard_pass))
      tbl <- dplyr::group_by(tbl, .alpha_key)
      tbl <- dplyr::slice(tbl, 1L)
      tbl <- dplyr::ungroup(tbl)
      dplyr::select(tbl, -.alpha_key)
    }
  
    evaluate_alpha_scaling <- function(alphas) {
      alphas <- sort(unique(alphas), decreasing = FALSE)
      candidate_results <- lapply(alphas, function(alpha_scaling) {
        fit_fuller_dual_core(
          stage2_df,
          outcome = outcome,
          predictor_u0 = predictor_u0,
          predictor_u1 = predictor_u1,
          meas11 = meas11,
          meas12 = meas12,
          meas22 = meas22,
          outcome_meas_var = outcome_meas_var,
          predictor_outcome_meas_cov_u0 = predictor_outcome_meas_cov_u0,
          predictor_outcome_meas_cov_u1 = predictor_outcome_meas_cov_u1,
          measurement_weight = 1,
          alpha_step1 = chosen_alpha_step1,
          alpha_step3 = chosen_alpha_step3,
          alpha_scaling = alpha_scaling,
          auto_tempered = TRUE,
          skip_internal_scaling = skip_internal_scaling
        )
      })
      score_alpha_scaling(dplyr::bind_rows(candidate_results))
    }
  
    base_scaling <- fit_fuller_dual_core(
      stage2_df,
      outcome = outcome,
      predictor_u0 = predictor_u0,
      predictor_u1 = predictor_u1,
      meas11 = meas11,
      meas12 = meas12,
      meas22 = meas22,
      outcome_meas_var = outcome_meas_var,
      predictor_outcome_meas_cov_u0 = predictor_outcome_meas_cov_u0,
      predictor_outcome_meas_cov_u1 = predictor_outcome_meas_cov_u1,
      measurement_weight = 1,
      alpha_step1 = chosen_alpha_step1,
      alpha_step3 = chosen_alpha_step3,
      alpha_scaling = alpha_lower,
      auto_tempered = TRUE,
      skip_internal_scaling = skip_internal_scaling
    )
    # The scaling-alpha upper bound should come from the centered-covariance
    # determinant root, not from `fuller_scaling_relative_min_eigen`. The latter
    # is only an admissibility diagnostic, so using it as lambda can arbitrarily
    # cap the search before a valid scaling alpha is reachable.
    scaling_lambda_hat <- if (nrow(base_scaling) > 0L && "fuller_lambda_scaling" %in% names(base_scaling)) {
      base_scaling$fuller_lambda_scaling[[1]]
    } else {
      NA_real_
    }
    alpha_scaling_upper <- alpha_upper_from_lambda(scaling_lambda_hat)
    if (!is.finite(alpha_scaling_upper)) {
      alpha_scaling_upper <- alpha_lower
    }
    chosen_alpha_scaling <- alpha_lower
    if (alpha_scaling_upper > alpha_lower) {
      scaling_choice <- refine_alpha_choice(
        evaluate_alpha_scaling(alpha_candidates(alpha_scaling_upper)),
        alpha_col = "fuller_alpha_scaling_used",
        pass_col = "alpha_scaling_guard_pass",
        evaluate_fn = evaluate_alpha_scaling,
        deduplicate_fn = deduplicate_alpha_scaling
      )
      chosen_alpha_scaling <- scaling_choice$alpha
    }
  } else {
    chosen_alpha_scaling <- NA
  }

  fit_fuller_dual_core(
    stage2_df,
    outcome = outcome,
    predictor_u0 = predictor_u0,
    predictor_u1 = predictor_u1,
    meas11 = meas11,
    meas12 = meas12,
    meas22 = meas22,
    outcome_meas_var = outcome_meas_var,
    predictor_outcome_meas_cov_u0 = predictor_outcome_meas_cov_u0,
    predictor_outcome_meas_cov_u1 = predictor_outcome_meas_cov_u1,
    measurement_weight = 1,
    alpha_step1 = chosen_alpha_step1,
    alpha_step3 = chosen_alpha_step3,
    alpha_scaling = chosen_alpha_scaling,
    auto_tempered = TRUE,
    skip_internal_scaling = skip_internal_scaling
  )
}

# alias
fit_fuller <- fit_fuller_dual

#' Format stacked-sandwich covariance variants as estimator rows.
#'
#' @details
#' `stacked_sandwich_for_corrected_scores()` returns a fitted stage-2
#' coefficient vector and several covariance matrix variants. Simulation
#' summaries expect one row per method, so this helper converts HC0-HC3
#' covariance entries into method-labelled rows with Wald confidence limits.
#'
#' @param sandwich_out List returned by `stacked_sandwich_for_corrected_scores()`.
#' @param df Degrees of freedom for the `t` critical value.
#' @param alpha_names Character vector naming the entries of
#' `sandwich_out$alpha_hat`. Defaults to `c("(Intercept)", "x")`.
#' @param term Character scalar naming the coefficient to report.
#' @param method_prefix Prefix used to build method names. Variants append
#' `"_hc0"` through `"_hc3"`.
#'
#' @return A tibble with `method`, `estimate`, `se`, `ci_low`, and `ci_high`.
format_stacked_sandwich_rows <- function(sandwich_out,
                                         df,
                                         alpha_names = c("(Intercept)", "x"),
                                         term = "x",
                                         method_prefix = "corrected_full_stacked") {
  variants <- c("hc0", "hc1", "hc2", "hc3")
  x_idx <- which(alpha_names == term)

  if (length(x_idx) != 1L) {
    return(purrr::map_dfr(variants, function(variant) {
      tibble::tibble(
        method = paste0(method_prefix, "_", variant),
        estimate = NA_real_,
        se = NA_real_,
        ci_low = NA_real_,
        ci_high = NA_real_
      )
    }))
  }

  crit <- stats::qt(0.975, df)
  purrr::map_dfr(variants, function(variant) {
    vcov_name <- paste0("vcov_", variant)
    vcov_mat <- sandwich_out[[vcov_name]]
    se <- if (!is.null(vcov_mat) && nrow(vcov_mat) >= x_idx && ncol(vcov_mat) >= x_idx) {
      sqrt(unname(vcov_mat[x_idx, x_idx]))
    } else {
      NA_real_
    }
    est <- unname(sandwich_out$alpha_hat[x_idx])

    tibble::tibble(
      method = paste0(method_prefix, "_", variant),
      estimate = est,
      se = se,
      ci_low = est - crit * se,
      ci_high = est + crit * se
    )
  })
}
