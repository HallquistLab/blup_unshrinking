# Posterior-reliability calibration for random-intercept/random-slope designs.
#
# This module converts design-level targets into the population parameters used
# by `simulate_dataset()`. The calibration has two distinct stages:
#
# 1. Solve for the marginal random-slope variance that produces a requested
#    expected posterior reliability under the planned Z and R matrices.
# 2. Decompose that marginal slope into `gamma * x + u1` for a requested
#    standardized structural coefficient, preserving the marginal
#    intercept-slope covariance.
#
# The target first-stage covariance and the covariance used to draw residual
# random effects are therefore different objects:
#
#   G_marginal = Var(b0, b1), where b1 = gamma * x + u1
#   G_residual = Var(b0, u1)
#
# In particular, the marginal correlation belongs in the reliability
# calculation, while the residual correlation belongs in `simulate_dataset()`.
# Keeping both explicit prevents the structural effect from silently changing
# the intended first-stage covariance structure.

#' Construct the within-cluster random-effect design used by simulate_dataset().
#'
#' @details
#' The simulation generator creates equally spaced time values on `[-1, 1]`,
#' centers them, and by default divides by their sample standard deviation.
#' This helper reproduces that transformation exactly so calibration uses the
#' same random-effect design matrix as data generation.
#'
#' With `standardize = TRUE`, the slope column satisfies
#' `sum(z^2) = n - 1`. Centering also makes the intercept and slope columns
#' orthogonal in balanced iid conditions. The returned matrix is the
#' random-effect design `Z_i`, not the fixed-effect design `X_i`.
#'
#' @param n Integer cluster size, at least 2.
#' @param standardize Logical; standardize the centered time scores to unit
#'   sample standard deviation. This should normally remain `TRUE` to match
#'   `simulate_dataset()`.
#'
#' @return A numeric `n x 2` matrix with columns named `intercept` and `slope`.
make_reliability_time_design <- function(n, standardize = TRUE) {
  n <- as.integer(n[[1]])
  if (!is.finite(n) || n < 2L) {
    stop("`n` must be one integer >= 2.")
  }

  # Match the time-score construction in `simulate_dataset()` rather than using
  # an abstract information approximation based only on cluster size.
  z <- seq(-1, 1, length.out = n)
  z <- z - mean(z)
  if (isTRUE(standardize)) {
    # `scale(..., scale = TRUE)` uses the sample SD, so S_zz equals n - 1.
    z <- as.numeric(scale(z, center = FALSE, scale = TRUE))
  }
  cbind(intercept = 1, slope = z)
}

#' Construct an iid, AR(1), or Toeplitz residual covariance matrix.
#'
#' This standalone constructor mirrors `make_R_matrix()` without requiring the
#' simulation helper module to have been sourced first. It supports the same
#' residual structures used by the simulation machinery:
#'
#' - `"iid"`: identity correlation matrix;
#' - `"ar1"`: correlation `rho^abs(j - k)`;
#' - `"toeplitz"`: user-supplied lag correlations, padded with zero at
#'   unspecified higher lags.
#'
#' The returned matrix is `R_i = sigma^2 C_i`. A Cholesky factorization is used
#' as a final positive-definiteness check because the posterior calculation
#' requires an invertible residual covariance matrix.
#'
#' @param n Integer cluster size.
#' @param sigma Positive residual standard deviation.
#' @param r_spec Residual covariance specification. May be `NULL`, a character
#'   scalar, or a list containing `structure`. AR(1) specifications require
#'   `rho`; Toeplitz specifications require `correlations` or use `rho` as a
#'   backward-compatible alias.
#'
#' @return A positive-definite numeric `n x n` residual covariance matrix.
make_reliability_residual_covariance <- function(n, sigma = 1, r_spec = NULL) {
  n <- as.integer(n[[1]])
  sigma <- as.numeric(sigma[[1]])
  if (!is.finite(n) || n < 1L || !is.finite(sigma) || sigma <= 0) {
    stop("`n` must be positive and `sigma` must be finite and positive.")
  }

  if (is.null(r_spec)) {
    r_spec <- list(structure = "iid")
  } else if (is.character(r_spec) && length(r_spec) == 1L) {
    r_spec <- list(structure = r_spec)
  }
  if (!is.list(r_spec)) {
    stop("`r_spec` must be NULL, a character scalar, or a list.")
  }

  structure_name <- r_spec$structure
  if (is.null(structure_name)) {
    structure_name <- r_spec$type
  }
  if (is.null(structure_name) || length(structure_name) != 1L ||
      !is.character(structure_name)) {
    stop("`r_spec` must include one character `structure` field.")
  }
  structure_name <- tolower(as.character(structure_name[[1]]))

  # Build a correlation matrix first, then apply the residual variance scale.
  # This separation makes the role of sigma in reliability calibration clear.
  cor_mat <- switch(
    structure_name,
    iid = diag(n),
    ar1 = {
      rho <- as.numeric(r_spec$rho[[1]])
      if (!is.finite(rho) || abs(rho) >= 1) {
        stop("AR(1) `r_spec` requires `rho` with absolute value < 1.")
      }
      idx <- seq_len(n)
      outer(idx, idx, function(a, b) rho^abs(a - b))
    },
    toeplitz = {
      correlations <- r_spec$correlations
      if (is.null(correlations)) {
        correlations <- r_spec$rho
      }
      correlations <- as.numeric(correlations)
      if (length(correlations) < 1L || any(!is.finite(correlations)) ||
          any(abs(correlations) >= 1)) {
        stop("Toeplitz `r_spec` requires finite lag correlations with absolute values < 1.")
      }
      lag_cor <- numeric(n)
      lag_cor[[1]] <- 1
      n_lag <- min(n - 1L, length(correlations))
      if (n_lag > 0L) {
        # Unspecified lags remain zero, matching `make_R_matrix()`.
        lag_cor[seq_len(n_lag) + 1L] <- correlations[seq_len(n_lag)]
      }
      stats::toeplitz(lag_cor)
    },
    stop("Unsupported residual covariance structure: ", structure_name)
  )

  R <- sigma^2 * cor_mat
  # Fail during condition construction rather than deep inside a simulation
  # replication if the requested residual structure is inadmissible.
  tryCatch(
    chol(R),
    error = function(e) stop("The requested residual covariance matrix is not positive definite.")
  )
  R
}

#' Construct a bivariate marginal random-effects covariance matrix.
#'
#' @details
#' Given marginal intercept variance `tau0^2`, marginal slope variance
#' `tau1^2`, and correlation `rho`, this function constructs
#'
#' `G = [[tau0^2, rho*tau0*tau1], [rho*tau0*tau1, tau1^2]]`.
#'
#' In posterior-reliability calibration this is the covariance of the total
#' random intercept and total random slope seen by the first-stage model. It is
#' not necessarily the covariance used to draw residual random effects after a
#' structural predictor is added to the slope.
#'
#' @param intercept_variance Positive marginal random-intercept variance.
#' @param slope_variance Positive marginal random-slope variance.
#' @param intercept_slope_correlation Correlation between the marginal random
#'   intercept and marginal random slope; must have absolute value below one.
#'
#' @return A symmetric `2 x 2` covariance matrix with intercept/slope dimnames.
make_random_effect_covariance <- function(
    intercept_variance,
    slope_variance,
    intercept_slope_correlation) {
  intercept_variance <- as.numeric(intercept_variance[[1]])
  slope_variance <- as.numeric(slope_variance[[1]])
  rho <- as.numeric(intercept_slope_correlation[[1]])

  if (!is.finite(intercept_variance) || intercept_variance <= 0 ||
      !is.finite(slope_variance) || slope_variance <= 0) {
    stop("Random-effect variances must be finite and positive.")
  }
  if (!is.finite(rho) || abs(rho) >= 1) {
    stop("`intercept_slope_correlation` must have absolute value < 1.")
  }

  # Convert the scale-free correlation into the covariance needed by G.
  covariance <- rho * sqrt(intercept_variance * slope_variance)
  matrix(
    c(intercept_variance, covariance, covariance, slope_variance),
    nrow = 2L,
    dimnames = list(c("intercept", "slope"), c("intercept", "slope"))
  )
}

#' Compute the conditional posterior covariance for one cluster.
#'
#' @details
#' For the Gaussian mixed model with prior covariance `G`, random-effect design
#' `Z`, and residual covariance `R`, the conditional covariance is
#'
#' `V_post = (G^{-1} + Z' R^{-1} Z)^{-1}`.
#'
#' This covariance depends only on the design and variance components, not on
#' the observed response values. That property allows reliability to be
#' calibrated before any Monte Carlo outcomes are generated.
#'
#' @param G Positive-definite `2 x 2` marginal random-effects covariance matrix.
#' @param Z Numeric `n x 2` random-effect design matrix.
#' @param R Positive-definite `n x n` residual covariance matrix.
#'
#' @return The symmetric `2 x 2` posterior covariance matrix for one cluster.
posterior_random_effect_covariance <- function(G, Z, R) {
  if (!all(dim(G) == c(2L, 2L)) || ncol(Z) != 2L ||
      nrow(R) != nrow(Z) || ncol(R) != nrow(Z)) {
    stop("Dimensions of `G`, `Z`, and `R` are incompatible.")
  }
  # `solve(R, Z)` avoids explicitly constructing R^{-1}; `crossprod` then
  # forms the likelihood information Z'R^{-1}Z.
  likelihood_information <- crossprod(Z, solve(R, Z))
  solve(solve(G) + likelihood_information)
}

#' Compute expected posterior slope reliability across reference clusters.
#'
#' @details
#' Slope reliability for cluster type `i` is defined as
#'
#' `1 - V_post_i[2, 2] / G[2, 2]`.
#'
#' This function averages the posterior-variance ratio before subtracting it
#' from one. When `weights` are supplied, `Z_list` and `R_list` may contain only
#' unique cluster types and the weights record how often each type occurs in a
#' deterministic reference population.
#'
#' The target is an expected design-level reliability. Individual clusters can
#' have different reliabilities when cluster sizes or residual structures vary.
#'
#' @param G Positive-definite `2 x 2` marginal random-effects covariance matrix.
#' @param Z_list Nonempty list of `n_i x 2` random-effect design matrices.
#' @param R_list Nonempty list of matching `n_i x n_i` residual covariance
#'   matrices.
#' @param weights Optional nonnegative frequency weights, one per list element.
#'   Defaults to equal weights.
#'
#' @return A numeric scalar in principle between zero and one.
expected_slope_reliability <- function(G, Z_list, R_list, weights = NULL) {
  if (length(Z_list) == 0L || length(Z_list) != length(R_list)) {
    stop("`Z_list` and `R_list` must be nonempty lists of equal length.")
  }
  if (is.null(weights)) {
    weights <- rep(1, length(Z_list))
  }
  weights <- as.numeric(weights)
  if (length(weights) != length(Z_list) || any(!is.finite(weights)) ||
      any(weights < 0) || sum(weights) <= 0) {
    stop("`weights` must be nonnegative, finite, and match `Z_list`.")
  }

  # Only the slope diagonal is needed for the scalar reliability target, but
  # each posterior is computed from the full 2 x 2 covariance structure.
  posterior_slope_variances <- vapply(
    seq_along(Z_list),
    function(i) posterior_random_effect_covariance(
      G,
      Z_list[[i]],
      R_list[[i]]
    )[2L, 2L],
    numeric(1L)
  )
  # Average the remaining uncertainty fraction, then convert to reliability.
  1 - stats::weighted.mean(
    posterior_slope_variances / G[2L, 2L],
    w = weights
  )
}

#' Compute posterior reliability for a linear random-effect contrast.
#'
#' @details
#' For a contrast `a`, reliability is
#'
#' `1 - E[a' V_post,i a] / (a' G a)`.
#'
#' Using the full posterior covariance keeps the definition valid when the
#' intercept and slope are correlated or when their sampling errors covary.
#'
#' @param G Positive-definite `2 x 2` random-effects covariance matrix.
#' @param Z_list,R_list Reference design and residual covariance lists.
#' @param contrast Length-two numeric contrast.
#' @param weights Optional nonnegative reference-cluster weights.
#'
#' @return Expected posterior reliability for the requested contrast.
expected_contrast_reliability <- function(
    G,
    Z_list,
    R_list,
    contrast,
    weights = NULL) {
  if (length(Z_list) == 0L || length(Z_list) != length(R_list)) {
    stop("`Z_list` and `R_list` must be nonempty lists of equal length.")
  }
  contrast <- as.numeric(contrast)
  if (length(contrast) != 2L || any(!is.finite(contrast))) {
    stop("`contrast` must contain two finite values.")
  }
  if (is.null(weights)) {
    weights <- rep(1, length(Z_list))
  }
  weights <- as.numeric(weights)
  if (length(weights) != length(Z_list) || any(!is.finite(weights)) ||
      any(weights < 0) || sum(weights) <= 0) {
    stop("`weights` must be nonnegative, finite, and match `Z_list`.")
  }

  prior_variance <- drop(t(contrast) %*% G %*% contrast)
  if (!is.finite(prior_variance) || prior_variance <= 0) {
    stop("The requested contrast must have positive finite variance under `G`.")
  }
  posterior_variances <- vapply(
    seq_along(Z_list),
    function(i) {
      V_post <- posterior_random_effect_covariance(
        G,
        Z_list[[i]],
        R_list[[i]]
      )
      drop(t(contrast) %*% V_post %*% contrast)
    },
    numeric(1L)
  )

  1 - stats::weighted.mean(posterior_variances / prior_variance, w = weights)
}

#' Contrast and posterior reliability for slope variation unique of intercept.
#'
#' @details
#' The residualized latent slope is
#'
#' `eta = u1 - G12 / G11 * u0`.
#'
#' Its variance is `G22 - G12^2 / G11`. In a stage-2 regression containing
#' both `u0` and `u1`, reparameterizing the predictors as `(u0, eta)` leaves
#' the coefficient on the focal slope unchanged. Its posterior reliability is
#' therefore aligned with the dual-predictor estimand.
#'
#' @param G Positive-definite `2 x 2` random-effects covariance matrix.
#' @param Z_list,R_list Reference design and residual covariance lists.
#' @param weights Optional nonnegative reference-cluster weights.
#'
#' @return Expected posterior reliability of the residualized slope.
expected_residualized_slope_reliability <- function(
    G,
    Z_list,
    R_list,
    weights = NULL) {
  contrast <- c(-G[1L, 2L] / G[1L, 1L], 1)
  expected_contrast_reliability(
    G = G,
    Z_list = Z_list,
    R_list = R_list,
    contrast = contrast,
    weights = weights
  )
}

#' Calibrate the level-1 residual scale at fixed covariance geometry.
#'
#' @details
#' Unlike `calibrate_slope_variance()`, this calibration holds the entire
#' random-effects covariance matrix fixed and solves for a common residual SD.
#' This permits matched-reliability comparisons without changing the ratio of
#' slope to intercept variance or the condition number of `G`.
#'
#' The supplied `R_shape_list` contains residual covariance matrices at unit
#' scale. Candidate matrices are `sigma^2 * R_shape`. The calibration may
#' target either marginal slope reliability or residualized-slope reliability.
#'
#' @param target_reliability Desired reliability strictly between zero and one.
#' @param G Fixed positive-definite `2 x 2` random-effects covariance matrix.
#' @param Z_list Reference random-effect design matrices.
#' @param R_shape_list Matching unit-scale residual covariance matrices.
#' @param weights Optional nonnegative reference-cluster weights.
#' @param reliability_measure Either `"marginal_slope"` or
#'   `"residualized_slope"`.
#' @param log_sigma_bounds Search interval for `log(sigma)`.
#' @param tolerance Root-finding tolerance.
#'
#' @return Calibration details, including `sigma`, both achieved reliability
#'   definitions, and the scaled residual covariance matrices.
calibrate_residual_scale <- function(
    target_reliability,
    G,
    Z_list,
    R_shape_list,
    weights = NULL,
    reliability_measure = c("marginal_slope", "residualized_slope"),
    log_sigma_bounds = c(-15, 15),
    tolerance = 1e-10) {
  target_reliability <- as.numeric(target_reliability[[1]])
  reliability_measure <- match.arg(reliability_measure)
  if (!is.finite(target_reliability) ||
      target_reliability <= 0 || target_reliability >= 1) {
    stop("`target_reliability` must be strictly between 0 and 1.")
  }
  if (!all(dim(G) == c(2L, 2L)) || any(!is.finite(G)) ||
      min(eigen(G, symmetric = TRUE, only.values = TRUE)$values) <= 0) {
    stop("`G` must be a finite positive-definite 2 x 2 matrix.")
  }
  if (length(Z_list) == 0L || length(Z_list) != length(R_shape_list)) {
    stop("`Z_list` and `R_shape_list` must be nonempty lists of equal length.")
  }

  reliability_at_log_sigma <- function(log_sigma) {
    sigma <- exp(log_sigma)
    R_list <- lapply(R_shape_list, function(R_shape) sigma^2 * R_shape)
    if (identical(reliability_measure, "marginal_slope")) {
      expected_slope_reliability(G, Z_list, R_list, weights = weights)
    } else {
      expected_residualized_slope_reliability(
        G, Z_list, R_list, weights = weights
      )
    }
  }

  attainable <- vapply(log_sigma_bounds, reliability_at_log_sigma, numeric(1L))
  attainable_range <- range(attainable)
  if (target_reliability < attainable_range[[1L]] - tolerance ||
      target_reliability > attainable_range[[2L]] + tolerance) {
    stop(sprintf(
      paste0(
        "Target reliability %.4f is outside the attainable interval ",
        "[%.4f, %.4f] for this fixed covariance geometry."
      ),
      target_reliability, attainable_range[[1L]], attainable_range[[2L]]
    ))
  }

  root <- stats::uniroot(
    function(log_sigma) {
      reliability_at_log_sigma(log_sigma) - target_reliability
    },
    interval = log_sigma_bounds,
    tol = tolerance
  )
  sigma <- exp(root$root)
  R_list <- lapply(R_shape_list, function(R_shape) sigma^2 * R_shape)
  achieved_marginal <- expected_slope_reliability(
    G, Z_list, R_list, weights = weights
  )
  achieved_residualized <- expected_residualized_slope_reliability(
    G, Z_list, R_list, weights = weights
  )

  list(
    target_reliability = target_reliability,
    reliability_measure = reliability_measure,
    achieved_reliability = if (identical(
      reliability_measure, "marginal_slope"
    )) achieved_marginal else achieved_residualized,
    achieved_marginal_slope_reliability = achieved_marginal,
    achieved_residualized_slope_reliability = achieved_residualized,
    sigma = sigma,
    residual_covariances = R_list,
    attainable_reliability = stats::setNames(
      attainable,
      c("lower_sigma_bound", "upper_sigma_bound")
    )
  )
}

#' Calibrate a measurement process while holding covariance shape fixed.
#'
#' @details
#' This is the shared calibration entry point for the amended Vig--Hallquist
#' studies. It fixes the population random-effect covariance geometry---the
#' intercept variance, slope/intercept variance ratio, and intercept--slope
#' correlation---and varies only the common Level-1 residual scale. The target
#' may be marginal posterior slope reliability, posterior reliability of the
#' slope residualized on the intercept, or the time-zero/intercept ICC.
#'
#' The ICC parameterization and posterior-reliability parameterization
#' therefore traverse the same fixed-shape DGM family. They differ in which
#' member of that family is held constant when cluster design changes.
#'
#' @param target_value Target value strictly between zero and one.
#' @param Z_list Reference random-effect design matrices.
#' @param R_shape_list Matching unit-diagonal residual correlation matrices.
#' @param weights Optional nonnegative reference-cluster weights.
#' @param intercept_variance Fixed marginal random-intercept variance.
#' @param slope_intercept_variance_ratio Fixed ratio `G[2,2] / G[1,1]`.
#' @param intercept_slope_correlation Fixed marginal correlation.
#' @param calibration_target One of `"marginal_slope"`,
#'   `"residualized_slope"`, or `"intercept_icc"`.
#' @param tolerance Numerical tolerance passed to posterior-reliability
#'   calibration and used for validation.
#'
#' @return A list containing the fixed `G`, solved `sigma`, scaled residual
#'   covariance matrices, both posterior reliabilities, ICC diagnostics, and
#'   covariance-geometry diagnostics.
calibrate_shape_preserving_measurement <- function(
    target_value,
    Z_list,
    R_shape_list,
    weights = NULL,
    intercept_variance = 0.81,
    slope_intercept_variance_ratio = 1,
    intercept_slope_correlation = 0,
    calibration_target = c(
      "marginal_slope", "residualized_slope", "intercept_icc"
    ),
    tolerance = 1e-10) {
  calibration_target <- match.arg(calibration_target)
  target_value <- as.numeric(target_value[[1L]])
  intercept_variance <- as.numeric(intercept_variance[[1L]])
  slope_intercept_variance_ratio <- as.numeric(
    slope_intercept_variance_ratio[[1L]]
  )
  intercept_slope_correlation <- as.numeric(
    intercept_slope_correlation[[1L]]
  )

  if (!is.finite(target_value) || target_value <= 0 || target_value >= 1) {
    stop("`target_value` must be strictly between zero and one.")
  }
  if (!is.finite(intercept_variance) || intercept_variance <= 0 ||
      !is.finite(slope_intercept_variance_ratio) ||
      slope_intercept_variance_ratio <= 0) {
    stop("The fixed covariance variances and ratio must be positive.")
  }
  if (length(Z_list) == 0L || length(Z_list) != length(R_shape_list)) {
    stop("`Z_list` and `R_shape_list` must be nonempty lists of equal length.")
  }
  unit_diagonal <- vapply(
    R_shape_list,
    function(R_shape) {
      is.matrix(R_shape) && all(is.finite(R_shape)) &&
        max(abs(diag(R_shape) - 1)) <= sqrt(tolerance)
    },
    logical(1L)
  )
  if (!all(unit_diagonal)) {
    stop(
      "`R_shape_list` must contain unit-diagonal residual correlation matrices."
    )
  }

  G <- make_random_effect_covariance(
    intercept_variance = intercept_variance,
    slope_variance = intercept_variance * slope_intercept_variance_ratio,
    intercept_slope_correlation = intercept_slope_correlation
  )

  if (identical(calibration_target, "intercept_icc")) {
    # ICC = tau0^2 / (tau0^2 + sigma^2). Holding G fixed and solving sigma is
    # scale-equivalent to Lai's total-variance-one parameterization.
    sigma <- sqrt(intercept_variance * (1 - target_value) / target_value)
    residual_covariances <- lapply(
      R_shape_list,
      function(R_shape) sigma^2 * R_shape
    )
    achieved_marginal <- expected_slope_reliability(
      G, Z_list, residual_covariances, weights = weights
    )
    achieved_residualized <- expected_residualized_slope_reliability(
      G, Z_list, residual_covariances, weights = weights
    )
  } else {
    residual_calibration <- calibrate_residual_scale(
      target_reliability = target_value,
      G = G,
      Z_list = Z_list,
      R_shape_list = R_shape_list,
      weights = weights,
      reliability_measure = calibration_target,
      tolerance = tolerance
    )
    sigma <- residual_calibration$sigma
    residual_covariances <- residual_calibration$residual_covariances
    achieved_marginal <-
      residual_calibration$achieved_marginal_slope_reliability
    achieved_residualized <-
      residual_calibration$achieved_residualized_slope_reliability
  }

  residualized_slope_variance <-
    G[2L, 2L] - G[1L, 2L]^2 / G[1L, 1L]
  intercept_icc <- G[1L, 1L] / (G[1L, 1L] + sigma^2)
  marginal_slope_icc <- G[2L, 2L] / (G[2L, 2L] + sigma^2)
  conditional_slope_icc <- residualized_slope_variance /
    (residualized_slope_variance + sigma^2)
  achieved_calibration_value <- switch(
    calibration_target,
    marginal_slope = achieved_marginal,
    residualized_slope = achieved_residualized,
    intercept_icc = intercept_icc
  )

  list(
    calibration_version = "shape_preserving_v2",
    calibration_target = calibration_target,
    calibration_target_value = target_value,
    achieved_calibration_value = achieved_calibration_value,
    achieved_marginal_slope_reliability = achieved_marginal,
    achieved_residualized_slope_reliability = achieved_residualized,
    intercept_icc = intercept_icc,
    marginal_slope_icc = marginal_slope_icc,
    conditional_slope_icc = conditional_slope_icc,
    intercept_variance = G[1L, 1L],
    slope_variance = G[2L, 2L],
    residualized_slope_variance = residualized_slope_variance,
    slope_intercept_variance_ratio = G[2L, 2L] / G[1L, 1L],
    intercept_slope_correlation = intercept_slope_correlation,
    G_condition_number = kappa(G, exact = TRUE),
    sigma = sigma,
    G_marginal = G,
    residual_covariances = residual_covariances
  )
}

#' Solve for marginal slope variance from a posterior-reliability target.
#'
#' @details
#' Reliability alone does not identify the full random-effects covariance
#' matrix. This function treats the intercept variance, marginal
#' intercept-slope correlation, `Z_i`, and `R_i` as fixed, and solves only for
#' the marginal slope variance.
#'
#' The root problem is
#'
#' `expected_slope_reliability(G(tau1^2)) - target_reliability = 0`.
#'
#' The search is performed on `log(tau1^2)`. This guarantees a positive
#' candidate variance and provides a stable search over many orders of
#' magnitude. Reliability at the two search bounds is also used to detect
#' infeasible targets before calling `uniroot()`.
#'
#' @param target_reliability Desired expected posterior slope reliability,
#'   strictly between zero and one.
#' @param Z_list Nonempty list of random-effect design matrices.
#' @param R_list Nonempty list of matching residual covariance matrices.
#' @param weights Optional nonnegative frequency weights for unique reference
#'   cluster types.
#' @param intercept_variance Fixed positive marginal intercept variance.
#' @param intercept_slope_correlation Fixed marginal intercept-slope
#'   correlation.
#' @param log_variance_bounds Length-two numeric interval for
#'   `log(slope_variance)`. The defaults cover an intentionally wide range.
#' @param tolerance Numerical tolerance passed to `uniroot()` and used for the
#'   attainable-range check.
#'
#' @return A list containing the target and achieved reliability, calibrated
#'   marginal slope variance, marginal `G` matrix, posterior covariance for each
#'   supplied cluster type, and reliability values at the search bounds.
calibrate_slope_variance <- function(
    target_reliability,
    Z_list,
    R_list,
    weights = NULL,
    intercept_variance = 1,
    intercept_slope_correlation = 0,
    log_variance_bounds = c(-30, 30),
    tolerance = 1e-10) {
  target_reliability <- as.numeric(target_reliability[[1]])
  if (!is.finite(target_reliability) ||
      target_reliability <= 0 || target_reliability >= 1) {
    stop("`target_reliability` must be strictly between 0 and 1.")
  }

  # Work on the log-variance scale so every candidate is strictly positive.
  reliability_at_log_variance <- function(log_slope_variance) {
    G <- make_random_effect_covariance(
      intercept_variance = intercept_variance,
      slope_variance = exp(log_slope_variance),
      intercept_slope_correlation = intercept_slope_correlation
    )
    expected_slope_reliability(G, Z_list, R_list, weights = weights)
  }

  # Correlated intercepts and slopes can impose a nonzero lower bound on slope
  # reliability. Evaluate the requested search interval before root finding so
  # infeasible design cells fail with an interpretable message.
  attainable <- vapply(
    log_variance_bounds,
    reliability_at_log_variance,
    numeric(1L)
  )
  if (target_reliability < attainable[[1L]] - tolerance ||
      target_reliability > attainable[[2L]] + tolerance) {
    stop(sprintf(
      paste0(
        "Target reliability %.4f is outside the attainable interval ",
        "[%.4f, %.4f] for this design and fixed covariance structure."
      ),
      target_reliability, attainable[[1L]], attainable[[2L]]
    ))
  }

  # Reliability is monotone in the marginal slope variance for the supported
  # Gaussian designs, making a one-dimensional bracketed root solve sufficient.
  root <- stats::uniroot(
    function(log_variance) {
      reliability_at_log_variance(log_variance) - target_reliability
    },
    interval = log_variance_bounds,
    tol = tolerance
  )

  # Convert the solved log variance back to the scale used in G and simulation.
  slope_variance <- exp(root$root)
  G <- make_random_effect_covariance(
    intercept_variance = intercept_variance,
    slope_variance = slope_variance,
    intercept_slope_correlation = intercept_slope_correlation
  )

  list(
    target_reliability = target_reliability,
    achieved_reliability = expected_slope_reliability(
      G,
      Z_list,
      R_list,
      weights = weights
    ),
    intercept_variance = as.numeric(intercept_variance[[1]]),
    slope_variance_marginal = slope_variance,
    intercept_slope_correlation_marginal = as.numeric(intercept_slope_correlation[[1]]),
    G_marginal = G,
    # Retain these matrices for diagnostics and independent verification of the
    # calibrated condition.
    posterior_covariances = Map(
      function(Z, R) posterior_random_effect_covariance(G, Z, R),
      Z_list,
      R_list
    ),
    attainable_reliability = stats::setNames(
      attainable,
      c("lower_bound", "upper_bound")
    )
  )
}

#' Decompose a calibrated marginal slope using a standardized coefficient.
#'
#' @details
#' The first-stage model sees a total slope
#'
#' `b1_i = gamma * x_i + u1_i`,
#'
#' where `Var(x) = 1` and `x` is independent of `(b0, u1)`. If
#' `standardized_beta = gamma / SD(b1)`, then
#'
#' `gamma = standardized_beta * SD(b1)`
#'
#' and
#'
#' `Var(u1) = (1 - standardized_beta^2) * Var(b1)`.
#'
#' The marginal covariance `Cov(b0, b1)` is preserved as `Cov(b0, u1)`.
#' Because the residual slope variance is smaller than the marginal slope
#' variance, the residual correlation passed to `simulate_dataset()` generally
#' differs from the marginal correlation used during reliability calibration.
#'
#' @param calibration Output from `calibrate_slope_variance()`.
#' @param standardized_beta Standardized coefficient for the structural
#'   predictor. Its absolute value must be below one.
#'
#' @return A list containing the raw structural coefficient, standardized beta,
#'   residual slope variance and standard deviation, preserved covariance,
#'   residual intercept-slope correlation, and residual `G` matrix used to draw
#'   `(b0, u1)`.
decompose_structural_slope <- function(
    calibration,
    standardized_beta = 0) {
  standardized_beta <- as.numeric(standardized_beta[[1]])
  if (!is.finite(standardized_beta) || abs(standardized_beta) >= 1) {
    stop("`standardized_beta` must be finite with absolute value below one.")
  }

  G_marginal <- calibration$G_marginal
  total_slope_variance <- G_marginal[2L, 2L]
  covariance_01 <- G_marginal[1L, 2L]
  structural_r_squared <- standardized_beta^2
  gamma <- standardized_beta * sqrt(total_slope_variance)
  residual_slope_variance <- total_slope_variance - gamma^2

  # Preserve the marginal intercept-slope covariance. Since x is independent of
  # b0, Cov(b0, b1) = Cov(b0, u1).
  G_residual <- matrix(
    c(
      G_marginal[1L, 1L], covariance_01,
      covariance_01, residual_slope_variance
    ),
    nrow = 2L,
    dimnames = dimnames(G_marginal)
  )
  eigenvalues <- eigen(G_residual, symmetric = TRUE, only.values = TRUE)$values
  if (min(eigenvalues) <= sqrt(.Machine$double.eps)) {
    stop(sprintf(
      paste0(
        "The requested standardized beta %.3f is incompatible with the ",
        "calibrated intercept-slope covariance. Residual G is not positive definite."
      ),
      standardized_beta
    ))
  }

  residual_correlation <- covariance_01 / sqrt(
    G_residual[1L, 1L] * G_residual[2L, 2L]
  )

  list(
    gamma_x_on_slope = gamma,
    standardized_beta = standardized_beta,
    structural_r_squared = structural_r_squared,
    slope_variance_marginal = total_slope_variance,
    slope_variance_residual = residual_slope_variance,
    slope_sd_residual = sqrt(residual_slope_variance),
    intercept_slope_covariance_residual = covariance_01,
    intercept_slope_correlation_residual = residual_correlation,
    G_residual = G_residual
  )
}

#' Calibrate a standardized latent-slope effect between two growth processes.
#'
#' @details
#' The predictor process has marginal random-effect covariance `G_predictor`.
#' The outcome process has a separately reliability-calibrated marginal
#' covariance `G_outcome`. The outcome slope is generated as
#'
#' `b_q1 = theta0 * b_y0 + theta1 * b_y1 + v`.
#'
#' Both structural coefficients are specified on a fully standardized scale.
#' The residual `(b_q0, v)` covariance preserves the calibrated marginal
#' covariance between the outcome intercept and outcome slope. Predictor-process
#' random effects are independent of this residual outcome block.
calibrate_dual_process_effect <- function(
    G_predictor,
    G_outcome,
    standardized_slope_beta,
    structural_target = c("slope_only", "intercept_slope"),
    nuisance_intercept_standardized_beta = 0) {
  structural_target <- match.arg(structural_target)
  standardized_slope_beta <- as.numeric(standardized_slope_beta[[1]])
  nuisance_intercept_standardized_beta <-
    as.numeric(nuisance_intercept_standardized_beta[[1]])

  validate_G <- function(G, name) {
    if (!all(dim(G) == c(2L, 2L)) || any(!is.finite(G)) ||
        min(eigen(G, symmetric = TRUE, only.values = TRUE)$values) <= 0) {
      stop(sprintf("`%s` must be a finite positive-definite 2 x 2 matrix.", name))
    }
  }
  validate_G(G_predictor, "G_predictor")
  validate_G(G_outcome, "G_outcome")
  if (!is.finite(standardized_slope_beta) ||
      !is.finite(nuisance_intercept_standardized_beta)) {
    stop("Standardized structural coefficients must be finite.")
  }

  outcome_slope_sd <- sqrt(G_outcome[2L, 2L])
  predictor_intercept_sd <- sqrt(G_predictor[1L, 1L])
  predictor_slope_sd <- sqrt(G_predictor[2L, 2L])
  theta0_std <- if (identical(structural_target, "slope_only")) {
    0
  } else {
    nuisance_intercept_standardized_beta
  }
  theta0 <- theta0_std * outcome_slope_sd / predictor_intercept_sd
  theta1 <- standardized_slope_beta * outcome_slope_sd / predictor_slope_sd
  coefficients <- c(intercept = theta0, slope = theta1)
  structural_variance <- drop(
    t(coefficients) %*% G_predictor %*% coefficients
  )
  residual_slope_variance <- G_outcome[2L, 2L] - structural_variance
  outcome_covariance <- G_outcome[1L, 2L]
  G_outcome_residual <- matrix(
    c(
      G_outcome[1L, 1L], outcome_covariance,
      outcome_covariance, residual_slope_variance
    ),
    nrow = 2L,
    dimnames = list(c("q_intercept", "q_slope_residual"),
                    c("q_intercept", "q_slope_residual"))
  )

  tolerance <- sqrt(.Machine$double.eps) * max(1, G_outcome[2L, 2L])
  residual_eigenvalues <- eigen(
    G_outcome_residual,
    symmetric = TRUE,
    only.values = TRUE
  )$values
  if (residual_slope_variance <= tolerance ||
      min(residual_eigenvalues) <= tolerance) {
    stop(sprintf(
      paste0(
        "Infeasible dual-process condition: standardized theta0=%.3f and ",
        "theta1=%.3f leave an inadmissible residual outcome-slope covariance."
      ),
      theta0_std, standardized_slope_beta
    ))
  }

  predictor_to_outcome_slope <- drop(G_predictor %*% coefficients)
  G_joint <- matrix(
    0,
    nrow = 4L,
    ncol = 4L,
    dimnames = list(
      c("y_intercept", "y_slope", "q_intercept", "q_slope"),
      c("y_intercept", "y_slope", "q_intercept", "q_slope")
    )
  )
  G_joint[1:2, 1:2] <- G_predictor
  G_joint[3:4, 3:4] <- G_outcome
  G_joint[1:2, 4] <- predictor_to_outcome_slope
  G_joint[4, 1:2] <- predictor_to_outcome_slope
  if (min(eigen(G_joint, symmetric = TRUE, only.values = TRUE)$values) <= tolerance) {
    stop("Calibrated dual-process marginal covariance is not positive definite.")
  }

  focal_predictor_variance <- if (identical(structural_target, "slope_only")) {
    G_predictor[2L, 2L]
  } else {
    G_predictor[2L, 2L] -
      G_predictor[1L, 2L]^2 / G_predictor[1L, 1L]
  }

  list(
    structural_target = structural_target,
    standardized_intercept_beta = theta0_std,
    standardized_slope_beta = standardized_slope_beta,
    theta0_intercept = theta0,
    theta1_slope = theta1,
    structural_variance = structural_variance,
    outcome_slope_residual_variance = residual_slope_variance,
    outcome_residual_correlation = outcome_covariance / sqrt(
      G_outcome_residual[1L, 1L] * G_outcome_residual[2L, 2L]
    ),
    total_structural_r_squared =
      structural_variance / G_outcome[2L, 2L],
    focal_unique_r_squared =
      theta1^2 * focal_predictor_variance / G_outcome[2L, 2L],
    G_outcome_residual = G_outcome_residual,
    G_joint_marginal = G_joint
  )
}

#' Calibrate a focal standardized random-slope effect on an external outcome.
#'
#' @details
#' For BLUPs used as predictors, first-stage reliability determines the
#' marginal random-effect covariance `G`; the structural model does not alter
#' that covariance. The focal effect is therefore calibrated directly as the
#' standardized regression coefficient
#'
#' `beta2_std = beta2 * sqrt(Var(b1)) / sqrt(Var(z))`.
#'
#' This parameterization isolates the slope coefficient when a fixed nuisance
#' intercept coefficient is also present. In contrast, targeting total
#' structural R-squared combines the intercept, slope, and their covariance and
#' cannot represent null or small focal slope effects below the nuisance-effect
#' variance floor.
#'
#' @param G_marginal Positive-definite `2 x 2` marginal random-effect
#'   covariance matrix.
#' @param standardized_slope_beta Target standardized coefficient for the
#'   random slope.
#' @param structural_target Either `"slope_only"` or `"intercept_slope"`.
#' @param nuisance_intercept_beta Raw coefficient for the random intercept in
#'   the dual-predictor target. Ignored for the slope-only target.
#' @param outcome_variance Fixed marginal variance of the external outcome.
#'
#' @return A list with raw coefficients, residual outcome variance, total
#' structural R-squared, and the focal slope's unique variance fraction.
calibrate_blup_predictor_effect <- function(
    G_marginal,
    standardized_slope_beta,
    structural_target = c("slope_only", "intercept_slope"),
    nuisance_intercept_beta = 0,
    outcome_variance = 1) {
  structural_target <- match.arg(structural_target)
  standardized_slope_beta <- as.numeric(standardized_slope_beta[[1]])
  nuisance_intercept_beta <- as.numeric(nuisance_intercept_beta[[1]])
  outcome_variance <- as.numeric(outcome_variance[[1]])

  if (!all(dim(G_marginal) == c(2L, 2L)) ||
      any(!is.finite(G_marginal)) ||
      min(eigen(G_marginal, symmetric = TRUE, only.values = TRUE)$values) <= 0) {
    stop("`G_marginal` must be a finite positive-definite 2 x 2 matrix.")
  }
  if (!is.finite(standardized_slope_beta) ||
      !is.finite(nuisance_intercept_beta) ||
      !is.finite(outcome_variance) || outcome_variance <= 0) {
    stop("Structural coefficients must be finite and `outcome_variance` must be positive.")
  }

  beta1 <- if (identical(structural_target, "slope_only")) {
    0
  } else {
    nuisance_intercept_beta
  }
  beta2 <- standardized_slope_beta *
    sqrt(outcome_variance / G_marginal[2L, 2L])
  coefficients <- c(intercept = beta1, slope = beta2)
  predictor_variance <- drop(t(coefficients) %*% G_marginal %*% coefficients)
  residual_variance <- outcome_variance - predictor_variance

  tolerance <- sqrt(.Machine$double.eps) * max(1, outcome_variance)
  if (residual_variance <= tolerance) {
    stop(sprintf(
      paste0(
        "Infeasible predictor model: explained variance %.6f leaves ",
        "non-positive residual variance for outcome variance %.6f."
      ),
      predictor_variance, outcome_variance
    ))
  }

  focal_predictor_variance <- if (identical(structural_target, "slope_only")) {
    G_marginal[2L, 2L]
  } else {
    G_marginal[2L, 2L] -
      G_marginal[1L, 2L]^2 / G_marginal[1L, 1L]
  }
  focal_unique_r_squared <- beta2^2 * focal_predictor_variance / outcome_variance

  list(
    structural_target = structural_target,
    standardized_slope_beta = standardized_slope_beta,
    beta1_intercept = beta1,
    beta2_slope = beta2,
    predictor_variance = predictor_variance,
    outcome_residual_variance = residual_variance,
    total_structural_r_squared = predictor_variance / outcome_variance,
    focal_unique_r_squared = focal_unique_r_squared
  )
}

#' Build a deterministic cluster-size profile for condition-level calibration.
#'
#' The profile approximates the population distribution of cluster sizes without
#' using replication-specific random draws. Calibrating against it once when the
#' condition grid is built keeps population parameters fixed across replications.
#'
#' @details
#' The returned counts reproduce the cluster-size mechanisms in
#' `simulate_dataset()` using deterministic midpoint quantiles instead of random
#' draws:
#'
#' - `"balanced"` returns one cluster size equal to `mean_n_trial`;
#' - `"unbalanced"` maps uniform midpoint quantiles through the simulation's
#'   `[0.6, 1.4] * mean_n_trial` rule;
#' - `"informative_unbalanced"` maps deterministic rank percentiles through the
#'   same decreasing power function used to relate cluster size to `x`.
#'
#' For unbalanced conditions, increasing `n_reference` gives a finer
#' approximation to the population distribution. These counts are calibration
#' devices only; each Monte Carlo replication still generates its own realized
#' cluster sizes using `simulate_dataset()`.
#'
#' @param mean_n_trial Target mean cluster size used by the simulation.
#' @param balance_mode One of `"balanced"`, `"unbalanced"`, or
#'   `"informative_unbalanced"`.
#' @param min_n_trial Lower bound for random unbalanced counts.
#' @param highly_unbalanced_min_n_trial Lower bound for informative unbalanced
#'   counts.
#' @param highly_unbalanced_power Positive power controlling the informative
#'   imbalance profile.
#' @param n_reference Number of deterministic midpoint quantiles for unbalanced
#'   profiles. Ignored in effect for balanced profiles, which need one type.
#'
#' @return An integer vector of deterministic reference cluster sizes.
reference_trial_counts <- function(
    mean_n_trial,
    balance_mode = "balanced",
    min_n_trial = 2L,
    highly_unbalanced_min_n_trial = 2L,
    highly_unbalanced_power = 3,
    n_reference = 1001L) {
  mean_n_trial <- as.integer(mean_n_trial[[1]])
  min_n_trial <- as.integer(min_n_trial[[1]])
  highly_unbalanced_min_n_trial <- as.integer(highly_unbalanced_min_n_trial[[1]])
  highly_unbalanced_power <- as.numeric(highly_unbalanced_power[[1]])
  n_reference <- as.integer(n_reference[[1]])
  balance_mode <- as.character(balance_mode[[1]])

  if (mean_n_trial < 2L || n_reference < 2L) {
    stop("`mean_n_trial` and `n_reference` must be at least 2.")
  }

  if (identical(balance_mode, "balanced")) {
    return(mean_n_trial)
  }

  # Midpoint quantiles avoid random-number generation and exclude exact 0/1
  # endpoints, making the profile deterministic and stable across runs.
  probabilities <- (seq_len(n_reference) - 0.5) / n_reference
  if (identical(balance_mode, "unbalanced")) {
    # Deterministic analogue of runif(n, 0.6, 1.4) in simulate_dataset().
    return(pmax(
      min_n_trial,
      as.integer(round((0.6 + 0.8 * probabilities) * mean_n_trial))
    ))
  }

  if (identical(balance_mode, "informative_unbalanced")) {
    if (mean_n_trial < highly_unbalanced_min_n_trial) {
      stop("`mean_n_trial` must be >= `highly_unbalanced_min_n_trial`.")
    }
    # The generator gives high-x subjects fewer observations. Reliability
    # calibration needs the resulting marginal count distribution; the
    # structural association with x is retained during actual simulation.
    reverse_percentile <- 1 - probabilities
    amplitude <- (mean_n_trial - highly_unbalanced_min_n_trial) *
      (highly_unbalanced_power + 1)
    return(pmax(
      highly_unbalanced_min_n_trial,
      as.integer(round(
        highly_unbalanced_min_n_trial +
          amplitude * reverse_percentile^highly_unbalanced_power
      ))
    ))
  }

  stop("Unsupported `balance_mode`: ", balance_mode)
}

#' Build deterministic Z and R reference lists for one simulation condition.
#'
#' @details
#' This function converts the deterministic count profile into the design
#' objects needed by `calibrate_slope_variance()`. Repeated cluster sizes are
#' compressed into unique `Z_i`/`R_i` pairs plus integer frequency weights.
#' Compression substantially reduces repeated matrix inversions while producing
#' the same weighted expected reliability as the full reference profile.
#'
#' @param mean_n_trial Target mean cluster size.
#' @param sigma Positive residual standard deviation.
#' @param balance_mode Cluster-size mechanism; see `reference_trial_counts()`.
#' @param min_n_trial Lower bound for random unbalanced counts.
#' @param highly_unbalanced_min_n_trial Lower bound for informative unbalanced
#'   counts.
#' @param highly_unbalanced_power Power controlling informative imbalance.
#' @param r_spec Residual covariance specification passed to
#'   `make_reliability_residual_covariance()`.
#' @param n_reference Number of deterministic reference quantiles.
#'
#' @return A list containing the full reference count profile, unique counts,
#'   frequency weights, matching `Z_list` and `R_list`, and count summaries.
make_reliability_reference_design <- function(
    mean_n_trial,
    sigma,
    balance_mode = "balanced",
    min_n_trial = 2L,
    highly_unbalanced_min_n_trial = 2L,
    highly_unbalanced_power = 3,
    r_spec = NULL,
    n_reference = 1001L) {
  counts <- reference_trial_counts(
    mean_n_trial = mean_n_trial,
    balance_mode = balance_mode,
    min_n_trial = min_n_trial,
    highly_unbalanced_min_n_trial = highly_unbalanced_min_n_trial,
    highly_unbalanced_power = highly_unbalanced_power,
    n_reference = n_reference
  )

  # Store one matrix pair per unique cluster size and carry multiplicities as
  # weights. This is exact for the deterministic reference profile.
  unique_counts <- sort(unique(counts))
  count_weights <- tabulate(match(counts, unique_counts), nbins = length(unique_counts))
  list(
    trial_counts = counts,
    unique_trial_counts = unique_counts,
    count_weights = count_weights,
    Z_list = lapply(unique_counts, make_reliability_time_design),
    R_list = lapply(
      unique_counts,
      make_reliability_residual_covariance,
      sigma = sigma,
      r_spec = r_spec
    ),
    mean_trial_count = mean(counts),
    min_trial_count = min(counts),
    max_trial_count = max(counts)
  )
}

#' Calibrate all simulation inputs for one reliability/beta condition.
#'
#' @details
#' This is the high-level condition-calibration entry point used by simulation
#' design builders. It:
#'
#' 1. constructs a deterministic reference distribution of cluster designs;
#' 2. solves for the marginal slope variance that attains the reliability target;
#' 3. decomposes the marginal slope according to standardized beta;
#' 4. returns the residual slope SD and residual correlation expected by
#'    `simulate_dataset()`.
#'
#' Call this function once when a condition manifest is constructed, not once
#' per Monte Carlo replication. Storing its scalar outputs in the manifest keeps
#' the population parameters fixed while realized subjects, cluster sizes, and
#' outcomes vary across replications.
#'
#' @param target_reliability Desired expected posterior slope reliability.
#' @param standardized_beta Standardized effect of the cluster-level predictor
#'   on the marginal random slope.
#' @param marginal_rho Target correlation between the total random intercept and
#'   total random slope in the first-stage marginal `G`.
#' @param tau0 Positive marginal random-intercept standard deviation.
#' @param mean_n_trial Target mean cluster size.
#' @param sigma Positive residual standard deviation.
#' @param balance_mode Cluster-size mechanism.
#' @param min_n_trial Lower bound for random unbalanced counts.
#' @param highly_unbalanced_min_n_trial Lower bound for informative unbalanced
#'   counts.
#' @param highly_unbalanced_power Power controlling informative imbalance.
#' @param r_spec Residual covariance specification.
#' @param n_reference Number of deterministic reference quantiles used for
#'   unbalanced calibration.
#' @return A list containing:
#'
#' - target and achieved posterior reliability;
#' - standardized beta and its derived structural `R^2`;
#' - marginal and residual intercept-slope correlations;
#' - cell-specific `gamma_x_on_slope`;
#' - marginal and residual slope variances;
#' - residual slope SD (`tau1_residual`);
#' - marginal and residual `G` matrices;
#' - summaries of the deterministic reference count profile.
calibrate_random_slope_condition <- function(
    target_reliability,
    standardized_beta,
    marginal_rho,
    tau0,
    mean_n_trial,
    sigma,
    balance_mode = "balanced",
    min_n_trial = 2L,
    highly_unbalanced_min_n_trial = 2L,
    highly_unbalanced_power = 3,
    r_spec = NULL,
    n_reference = 1001L) {
  # The reference profile is deterministic, so repeated calls for the same
  # condition return identical population parameters.
  reference <- make_reliability_reference_design(
    mean_n_trial = mean_n_trial,
    sigma = sigma,
    balance_mode = balance_mode,
    min_n_trial = min_n_trial,
    highly_unbalanced_min_n_trial = highly_unbalanced_min_n_trial,
    highly_unbalanced_power = highly_unbalanced_power,
    r_spec = r_spec,
    n_reference = n_reference
  )
  # Reliability is defined using the marginal G seen by a Stage-1 model that
  # omits the structural predictor x.
  calibration <- calibrate_slope_variance(
    target_reliability = target_reliability,
    Z_list = reference$Z_list,
    R_list = reference$R_list,
    weights = reference$count_weights,
    intercept_variance = tau0^2,
    intercept_slope_correlation = marginal_rho
  )
  # Convert the marginal slope into the residual random-effect parameters that
  # the existing data generator already accepts.
  structural <- decompose_structural_slope(
    calibration,
    standardized_beta = standardized_beta
  )

  list(
    target_reliability = calibration$target_reliability,
    achieved_reliability = calibration$achieved_reliability,
    structural_r_squared = structural$structural_r_squared,
    standardized_beta = structural$standardized_beta,
    marginal_rho = as.numeric(marginal_rho[[1]]),
    rho_residual = structural$intercept_slope_correlation_residual,
    gamma_x_on_slope = structural$gamma_x_on_slope,
    slope_variance_marginal = structural$slope_variance_marginal,
    slope_variance_residual = structural$slope_variance_residual,
    tau1_residual = structural$slope_sd_residual,
    G_marginal = calibration$G_marginal,
    G_residual = structural$G_residual,
    reference_mean_n_trial = reference$mean_trial_count,
    reference_min_n_trial = reference$min_trial_count,
    reference_max_n_trial = reference$max_trial_count,
    calibration_reference_n = length(reference$trial_counts)
  )
}

#' Calibrate reliability and a standardized BLUP-as-predictor effect.
#'
#' This is the predictor-side counterpart to
#' `calibrate_random_slope_condition()`. It retains the full marginal `G` for
#' data generation and calibrates the external-outcome coefficient without
#' decomposing the random slope.
calibrate_blup_predictor_condition <- function(
    target_reliability,
    standardized_slope_beta,
    structural_target,
    marginal_rho,
    tau0,
    mean_n_trial,
    sigma,
    nuisance_intercept_beta = 0,
    outcome_variance = 1,
    balance_mode = "balanced",
    min_n_trial = 2L,
    highly_unbalanced_min_n_trial = 2L,
    highly_unbalanced_power = 3,
    r_spec = NULL,
    n_reference = 1001L) {
  reference <- make_reliability_reference_design(
    mean_n_trial = mean_n_trial,
    sigma = sigma,
    balance_mode = balance_mode,
    min_n_trial = min_n_trial,
    highly_unbalanced_min_n_trial = highly_unbalanced_min_n_trial,
    highly_unbalanced_power = highly_unbalanced_power,
    r_spec = r_spec,
    n_reference = n_reference
  )
  reliability <- calibrate_slope_variance(
    target_reliability = target_reliability,
    Z_list = reference$Z_list,
    R_list = reference$R_list,
    weights = reference$count_weights,
    intercept_variance = tau0^2,
    intercept_slope_correlation = marginal_rho
  )
  structural <- calibrate_blup_predictor_effect(
    G_marginal = reliability$G_marginal,
    standardized_slope_beta = standardized_slope_beta,
    structural_target = structural_target,
    nuisance_intercept_beta = nuisance_intercept_beta,
    outcome_variance = outcome_variance
  )

  c(
    list(
      target_reliability = reliability$target_reliability,
      achieved_reliability = reliability$achieved_reliability,
      marginal_rho = as.numeric(marginal_rho[[1]]),
      slope_variance_marginal = reliability$slope_variance_marginal,
      tau1_marginal = sqrt(reliability$slope_variance_marginal),
      G_marginal = reliability$G_marginal,
      reference_mean_n_trial = reference$mean_trial_count,
      reference_min_n_trial = reference$min_trial_count,
      reference_max_n_trial = reference$max_trial_count,
      calibration_reference_n = length(reference$trial_counts)
    ),
    structural
  )
}
