#!/usr/bin/env Rscript

# Focused audit simulation for the Fuller EIV implementation.
#
# This file deliberately implements the governing equations independently of
# R/stage2_estimators.R.  Its "current" branch reproduces the production
# implementation at the time of the audit, while the other branches change one
# algebraic choice at a time. The corrected cross-moment branch now reproduces
# production; the recycled branch is retained as a historical bug comparator.
# Keeping the audit implementation separate makes exact agreement with
# production a meaningful check rather than a tautology.

suppressPackageStartupMessages({
  library(data.table)
  library(geigen)
})

smallest_generalized_root <- function(a_mat, b_mat) {
  out <- tryCatch(
    geigen::geigen(a_mat, b_mat, symmetric = FALSE, only.values = TRUE),
    error = function(e) NULL
  )
  if (is.null(out) || is.null(out$values)) {
    return(NA_real_)
  }
  values <- out$values
  if (is.complex(values)) {
    imag_tol <- sqrt(.Machine$double.eps) *
      max(1, suppressWarnings(max(abs(values), na.rm = TRUE)))
    values <- ifelse(abs(Im(values)) <= imag_tol, Re(values), NA_real_)
  }
  values <- values[is.finite(values)]
  if (length(values) == 0L) NA_real_ else min(values)
}

sum_covariance_array <- function(cov_array, weights = NULL) {
  n <- dim(cov_array)[1L]
  if (is.null(weights)) weights <- rep(1, n)
  out <- matrix(0, nrow = dim(cov_array)[2L], ncol = dim(cov_array)[3L])
  for (i in seq_len(n)) out <- out + weights[[i]] * cov_array[i, , ]
  out
}

quadratic_rows <- function(cov_array, coefficient) {
  vapply(
    seq_len(dim(cov_array)[1L]),
    function(i) drop(crossprod(coefficient, cov_array[i, , ] %*% coefficient)),
    numeric(1L)
  )
}

omega_times_coefficient <- function(cov_array, coefficient) {
  t(vapply(
    seq_len(dim(cov_array)[1L]),
    function(i) drop(cov_array[i, , ] %*% coefficient),
    numeric(length(coefficient))
  ))
}

make_full_omega_array <- function(omega_y, omega_xy, omega_x) {
  n <- length(omega_y)
  p <- ncol(omega_xy)
  out <- array(0, dim = c(n, p + 1L, p + 1L))
  for (i in seq_len(n)) {
    out[i, 1L, 1L] <- omega_y[[i]]
    out[i, 1L, 2L:(p + 1L)] <- omega_xy[i, ]
    out[i, 2L:(p + 1L), 1L] <- omega_xy[i, ]
    out[i, 2L:(p + 1L), 2L:(p + 1L)] <- omega_x[i, , ]
  }
  out
}

#' Independently implement one audited Fuller variant.
#'
#' `preliminary = "current_recycled"` reproduces the former Step-1 bug:
#' it applies the lambda/alpha correction to the preliminary moments and
#' subtracts the scalar sum of outcome-error variances from every element of
#' X'y. `current_cross` retains that finite-sample correction but uses the
#' predictor-outcome error covariance. `documented` implements Equation 58 of
#' documentation/BLUP_Project.pdf.
#'
#' `se_bread = "current_final"` uses the final modified S_x* matrix, as the
#' corrected VH primary estimator does. `documented_full` uses Equation 66's fully
#' measurement-error-subtracted S_x matrix.
fit_fuller_audit <- function(y,
                             x,
                             omega_x,
                             omega_y,
                             omega_xy,
                             preliminary = c(
                               "current_recycled", "current_cross", "documented"
                             ),
                             se_bread = c("current_final", "documented_full")) {
  preliminary <- match.arg(preliminary)
  se_bread <- match.arg(se_bread)
  x <- as.matrix(x)
  y <- as.numeric(y)
  omega_y <- as.numeric(omega_y)
  omega_xy <- as.matrix(omega_xy)
  m <- nrow(x)
  p <- ncol(x)
  alpha <- p + 1

  fail <- function(reason) {
    data.table(
      estimate = NA_real_, se = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
      status = reason, lambda1 = NA_real_, lambda2 = NA_real_,
      correction1 = NA_real_, correction2 = NA_real_, sigma2 = NA_real_,
      gamma0_target = NA_real_, bread_min_eigen = NA_real_
    )
  }

  if (length(y) != m || length(omega_y) != m || nrow(omega_xy) != m ||
      ncol(omega_xy) != p || !identical(dim(omega_x), c(m, p, p))) {
    return(fail("dimension_mismatch"))
  }

  omega_full <- make_full_omega_array(omega_y, omega_xy, omega_x)
  omega_full_sum <- sum_covariance_array(omega_full)
  omega_x_sum <- sum_covariance_array(omega_x)
  omega_xy_sum <- colSums(omega_xy)
  b <- cbind(y, x)

  lambda1 <- smallest_generalized_root(crossprod(b), omega_full_sum)
  correction1 <- if (is.finite(lambda1) && lambda1 <= 1 + 1 / m) {
    lambda1 - 1 / m - alpha / m
  } else {
    1 - alpha / m
  }

  if (preliminary == "current_recycled") {
    a0 <- crossprod(x) - correction1 * omega_x_sum
    # Intentional scalar recycling to reproduce R/stage2_estimators.R exactly.
    d0 <- as.vector(crossprod(x, y) - correction1 * sum(omega_y))
  } else if (preliminary == "current_cross") {
    a0 <- crossprod(x) - correction1 * omega_x_sum
    d0 <- as.vector(crossprod(x, y) - correction1 * omega_xy_sum)
  } else {
    a0 <- crossprod(x) - omega_x_sum
    d0 <- as.vector(crossprod(x, y) - omega_xy_sum)
  }

  gamma0 <- tryCatch(as.vector(solve(a0, d0)), error = function(e) NULL)
  if (is.null(gamma0) || any(!is.finite(gamma0))) return(fail("gamma0_failure"))

  composite_me <- omega_y - 2 * rowSums(omega_xy * rep(gamma0, each = m)) +
    quadratic_rows(omega_x, gamma0)
  residual0 <- y - drop(x %*% gamma0)
  sigma2 <- if (is.finite(lambda1) && lambda1 < 1) {
    0
  } else {
    sum(residual0^2) / max(1, m - p) - mean(composite_me)
  }
  sigma2 <- max(0, sigma2)
  weights <- sigma2 + composite_me
  if (any(!is.finite(weights)) || any(weights <= sqrt(.Machine$double.eps))) {
    return(fail("weight_failure"))
  }

  w_inv <- 1 / weights
  b_weighted_sum <- crossprod(b * sqrt(w_inv))
  omega_weighted_sum <- sum_covariance_array(omega_full, w_inv)
  lambda2 <- smallest_generalized_root(b_weighted_sum, omega_weighted_sum)
  correction2 <- if (is.finite(lambda2) && lambda2 <= 1 + 1 / m) {
    lambda2 - 1 / m - alpha / m
  } else {
    1 - alpha / m
  }
  s_star <- b_weighted_sum - correction2 * omega_weighted_sum
  s_x_star <- s_star[2L:(p + 1L), 2L:(p + 1L), drop = FALSE]
  s_xy_star <- s_star[2L:(p + 1L), 1L, drop = FALSE]
  gamma <- tryCatch(as.vector(solve(s_x_star, s_xy_star)), error = function(e) NULL)
  if (is.null(gamma) || any(!is.finite(gamma))) return(fail("gamma_failure"))

  x_weighted_sum <- crossprod(x * sqrt(w_inv))
  omega_x_weighted_sum <- sum_covariance_array(omega_x, w_inv)
  bread <- if (se_bread == "current_final") {
    s_x_star / m
  } else {
    (x_weighted_sum - omega_x_weighted_sum) / m
  }
  bread <- (bread + t(bread)) / 2
  bread_inv <- tryCatch(solve(bread), error = function(e) NULL)
  if (is.null(bread_inv) || any(!is.finite(bread_inv))) return(fail("bread_failure"))

  tilde <- omega_xy - omega_times_coefficient(omega_x, gamma0)
  meat <- x_weighted_sum + crossprod(tilde * w_inv)
  vcov <- (bread_inv %*% meat %*% bread_inv) / m^2
  vcov <- (vcov + t(vcov)) / 2
  variance <- vcov[p, p]
  if (!is.finite(variance) || variance < 0) return(fail("variance_failure"))

  estimate <- gamma[[p]]
  se <- sqrt(variance)
  data.table(
    estimate = estimate,
    se = se,
    ci_low = estimate - qnorm(.975) * se,
    ci_high = estimate + qnorm(.975) * se,
    status = "ok",
    lambda1 = lambda1,
    lambda2 = lambda2,
    correction1 = correction1,
    correction2 = correction2,
    sigma2 = sigma2,
    gamma0_target = gamma0[[p]],
    bread_min_eigen = min(eigen(bread, symmetric = TRUE, only.values = TRUE)$values)
  )
}

audit_variants <- data.table(
  variant = c(
    "current",
    "cross_moment_fixed",
    "documented_preliminary",
    "documented_full"
  ),
  preliminary = c(
    "current_recycled",
    "current_cross",
    "documented",
    "documented"
  ),
  se_bread = c(
    "current_final",
    "current_final",
    "current_final",
    "documented_full"
  )
)

geometry_parameters <- function(geometry) {
  if (geometry == "regular") {
    ratio <- 1
  } else if (geometry == "weak_shape") {
    ratio <- 0.005
  } else {
    stop("Unknown geometry: ", geometry)
  }
  rho <- 0.5
  g <- matrix(
    c(1, rho * sqrt(ratio), rho * sqrt(ratio), ratio),
    nrow = 2L,
    byrow = TRUE
  )
  g_eigenvalues <- eigen(g, symmetric = TRUE, only.values = TRUE)$values
  list(
    g = g,
    ratio = ratio,
    rho = rho,
    kappa = max(g_eigenvalues) / min(g_eigenvalues)
  )
}

measurement_parameters <- function(g, scenario, m) {
  has_predictor_error <- scenario %in% c(
    "predictor_only", "both_homoskedastic", "both_heteroskedastic"
  )
  has_outcome_error <- scenario %in% c(
    "outcome_homoskedastic", "outcome_heteroskedastic",
    "both_homoskedastic", "both_heteroskedastic"
  )
  heterogeneous <- scenario %in% c("outcome_heteroskedastic", "both_heteroskedastic")

  predictor_base <- matrix(0, 2L, 2L)
  if (has_predictor_error) {
    predictor_variance <- c(g[1, 1] * (1 - .80) / .80, g[2, 2] * (1 - .25) / .25)
    predictor_base <- matrix(
      c(
        predictor_variance[[1]], -.20 * sqrt(prod(predictor_variance)),
        -.20 * sqrt(prod(predictor_variance)), predictor_variance[[2]]
      ),
      nrow = 2L,
      byrow = TRUE
    )
  }
  outcome_base <- if (has_outcome_error) 1 else 0
  multiplier <- if (heterogeneous) {
    rep(c(.2, 1, 1.8), length.out = m)
  } else {
    rep(1, m)
  }

  p <- 3L
  omega_x <- array(0, dim = c(m, p, p))
  for (i in seq_len(m)) {
    omega_x[i, 2:3, 2:3] <- multiplier[[i]] * predictor_base
  }
  list(
    omega_x = omega_x,
    omega_y = multiplier * outcome_base,
    omega_xy = matrix(0, nrow = m, ncol = p),
    predictor_base = predictor_base,
    outcome_base = outcome_base,
    multiplier = multiplier
  )
}

simulate_audit_replication <- function(m, geometry, scenario) {
  geometry_info <- geometry_parameters(geometry)
  g <- geometry_info$g
  standardized_beta <- c(.20, .40)
  beta <- standardized_beta / sqrt(diag(g))
  explained <- drop(crossprod(beta, g %*% beta))
  residual_variance <- 1 - explained

  z <- matrix(rnorm(m * 2L), ncol = 2L)
  true_x <- z %*% chol(g)
  true_y <- drop(true_x %*% beta) + rnorm(m, sd = sqrt(residual_variance))
  measurement <- measurement_parameters(g, scenario, m)

  observed_x <- true_x
  for (i in seq_len(m)) {
    pred_cov <- measurement$omega_x[i, 2:3, 2:3, drop = FALSE]
    dim(pred_cov) <- c(2L, 2L)
    if (max(abs(pred_cov)) > 0) {
      observed_x[i, ] <- observed_x[i, ] +
        drop(rnorm(2L) %*% chol(pred_cov))
    }
  }
  observed_y <- true_y + rnorm(m, sd = sqrt(measurement$omega_y))
  x_design <- cbind(1, observed_x)

  list(
    y = observed_y,
    x = x_design,
    omega_x = measurement$omega_x,
    omega_y = measurement$omega_y,
    omega_xy = measurement$omega_xy,
    truth_raw = beta[[2L]],
    reporting_scale = sqrt(g[2L, 2L]),
    g_condition = geometry_info$kappa,
    variance_ratio = geometry_info$ratio
  )
}

run_audit_simulation <- function(n_rep = 1000L, seed = 20260817L, cores = 1L) {
  cells <- CJ(
    scenario = c(
      "no_measurement_error",
      "outcome_homoskedastic",
      "outcome_heteroskedastic",
      "predictor_only",
      "both_homoskedastic",
      "both_heteroskedastic"
    ),
    geometry = c("regular", "weak_shape"),
    m = c(50L, 200L),
    sorted = FALSE
  )
  cells[, cell_id := .I]

  run_cell <- function(cell_id) {
    condition <- cells[cell_id]
    set.seed(seed + 100003L * cell_id)
    output <- vector("list", n_rep * nrow(audit_variants))
    out_index <- 0L
    for (rep_id in seq_len(n_rep)) {
      simulated <- simulate_audit_replication(
        m = condition$m,
        geometry = condition$geometry,
        scenario = condition$scenario
      )
      for (variant_row in seq_len(nrow(audit_variants))) {
        variant <- audit_variants[variant_row]
        fit <- fit_fuller_audit(
          y = simulated$y,
          x = simulated$x,
          omega_x = simulated$omega_x,
          omega_y = simulated$omega_y,
          omega_xy = simulated$omega_xy,
          preliminary = variant$preliminary,
          se_bread = variant$se_bread
        )
        fit[, `:=`(
          scenario = condition$scenario,
          geometry = condition$geometry,
          m = condition$m,
          rep = rep_id,
          variant = variant$variant,
          truth = simulated$truth_raw * simulated$reporting_scale,
          estimate = estimate * simulated$reporting_scale,
          se = se * simulated$reporting_scale,
          ci_low = ci_low * simulated$reporting_scale,
          ci_high = ci_high * simulated$reporting_scale,
          gamma0_target = gamma0_target * simulated$reporting_scale,
          g_condition = simulated$g_condition,
          variance_ratio = simulated$variance_ratio
        )]
        out_index <- out_index + 1L
        output[[out_index]] <- fit
      }
    }
    rbindlist(output, fill = TRUE)
  }

  cores <- max(1L, min(as.integer(cores), nrow(cells)))
  output <- if (cores > 1L && .Platform$OS.type != "windows") {
    parallel::mclapply(seq_len(nrow(cells)), run_cell, mc.cores = cores)
  } else {
    lapply(seq_len(nrow(cells)), run_cell)
  }
  rbindlist(output, fill = TRUE)
}

validate_corrected_against_production <- function(repo_root, n_rep = 100L, seed = 81726L) {
  source(file.path(repo_root, "R", "stats_helpers.R"))
  source(file.path(repo_root, "R", "stage2_estimators.R"))
  set.seed(seed)
  output <- vector("list", n_rep)
  scenarios <- c("outcome_homoskedastic", "both_heteroskedastic")
  geometries <- c("regular", "weak_shape")

  for (i in seq_len(n_rep)) {
    scenario <- scenarios[(i - 1L) %% length(scenarios) + 1L]
    geometry <- geometries[(i - 1L) %% length(geometries) + 1L]
    simulated <- simulate_audit_replication(80L, geometry, scenario)
    dat <- data.frame(
      y = simulated$y,
      x0 = simulated$x[, 2L],
      x1 = simulated$x[, 3L],
      meas11 = simulated$omega_x[, 2L, 2L],
      meas12 = simulated$omega_x[, 2L, 3L],
      meas22 = simulated$omega_x[, 3L, 3L],
      outcome_meas_var = simulated$omega_y,
      outcome_predictor_cov_u0 = simulated$omega_xy[, 2L],
      outcome_predictor_cov_u1 = simulated$omega_xy[, 3L]
    )
    production <- fit_fuller_dual(
      dat,
      outcome = "y",
      predictor_u0 = "x0",
      predictor_u1 = "x1",
      meas11 = "meas11",
      meas12 = "meas12",
      meas22 = "meas22",
      outcome_meas_var = "outcome_meas_var",
      predictor_outcome_meas_cov_u0 = "outcome_predictor_cov_u0",
      predictor_outcome_meas_cov_u1 = "outcome_predictor_cov_u1",
      skip_internal_scaling = TRUE
    )
    independent <- fit_fuller_audit(
      y = simulated$y,
      x = simulated$x,
      omega_x = simulated$omega_x,
      omega_y = simulated$omega_y,
      omega_xy = simulated$omega_xy,
      preliminary = "current_cross",
      se_bread = "current_final"
    )
    output[[i]] <- data.table(
      replication = i,
      scenario = scenario,
      geometry = geometry,
      production_status = production$status_code[[1L]],
      audit_status = independent$status[[1L]],
      estimate_difference = production$estimate[[1L]] - independent$estimate[[1L]],
      se_difference = production$se[[1L]] - independent$se[[1L]],
      lambda1_difference = production$fuller_lambda1[[1L]] - independent$lambda1[[1L]],
      lambda2_difference = production$fuller_lambda2[[1L]] - independent$lambda2[[1L]],
      sigma2_difference = production$fuller_sigma2[[1L]] - independent$sigma2[[1L]]
    )
  }
  rbindlist(output, fill = TRUE)
}

summarize_audit_simulation <- function(results) {
  results[, eligible := status == "ok" & is.finite(estimate) & is.finite(se) & se > 0]
  results[, covered := eligible & ci_low <= truth & ci_high >= truth]
  results[, .(
    n_rep = .N,
    eligibility = mean(eligible),
    mean_estimate = mean(estimate[eligible]),
    truth = unique(truth)[1L],
    bias = mean(estimate[eligible] - truth[eligible]),
    empirical_sd = sd(estimate[eligible]),
    mean_se = mean(se[eligible]),
    median_se = median(se[eligible]),
    mean_se_ratio = mean(se[eligible]) / sd(estimate[eligible]),
    median_se_ratio = median(se[eligible]) / sd(estimate[eligible]),
    coverage = mean(covered[eligible]),
    mean_gamma0 = mean(gamma0_target[eligible]),
    median_bread_min_eigen = median(bread_min_eigen[eligible]),
    mean_sigma2 = mean(sigma2[eligible]),
    g_condition = unique(g_condition)[1L],
    variance_ratio = unique(variance_ratio)[1L]
  ), by = .(scenario, geometry, m, variant)]
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  n_rep <- if (length(args) >= 1L) as.integer(args[[1L]]) else 1000L
  output_dir <- if (length(args) >= 2L) args[[2L]] else
    file.path("outputs", "fuller_algebra_audit")
  cores <- if (length(args) >= 3L) as.integer(args[[3L]]) else
    min(4L, parallel::detectCores(logical = FALSE))
  repo_root <- normalizePath(".", mustWork = TRUE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  message("Running Fuller algebra audit with ", n_rep, " replications per cell")
  validation <- validate_corrected_against_production(repo_root, n_rep = 100L)
  results <- run_audit_simulation(n_rep = n_rep, cores = cores)
  summary <- summarize_audit_simulation(results)

  fwrite(
    validation,
    file.path(output_dir, "fuller_algebra_production_validation.csv")
  )
  fwrite(
    results,
    file.path(output_dir, "fuller_algebra_replication_results.csv.gz")
  )
  fwrite(
    summary,
    file.path(output_dir, "fuller_algebra_summary.csv")
  )
  message("Wrote audit outputs to ", output_dir)
}
