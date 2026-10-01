#!/usr/bin/env Rscript
# Run from the repository root: Rscript --vanilla tests/test_fuller_mlm.R
source('R/fuller_mlm.R')
source('R/stage2_estimators.R')
# The factor-score prototype is not required by this commit. Opt in to its
# additional reference comparisons with --with-fsr when it is available.
with_fsr <- '--with-fsr' %in% commandArgs(trailingOnly = TRUE)
if (with_fsr) source('vh_fsr/R/vh_fuller_q.R')

close <- function(actual, expected, label, tolerance = 1e-9) {
  if (!isTRUE(all.equal(actual, expected, tolerance = tolerance,
                        check.attributes = FALSE))) stop(label)
}
expect_error <- function(expr, pattern) {
  err <- tryCatch(force(expr), error = identity)
  stopifnot(inherits(err, 'error'), grepl(pattern, conditionMessage(err)))
}
compare <- function(a, b) {
  stopifnot(a$converged, b$converged)
  close(a$coefficients, b$coefficients, 'Coefficient/SE/CI mismatch')
  close(a$vcov, b$vcov, 'Full covariance mismatch')
  close(a$preliminary_coefficients, b$preliminary_coefficients, 'Preliminary mismatch')
  close(a$analysis_weights, b$analysis_weights, 'Weights mismatch')
  fields <- c('lambda1', 'lambda2', 'correction1', 'correction3', 'sigma2')
  close(unlist(a$diagnostics[fields]), unlist(b$diagnostics[fields]), 'Diagnostics mismatch')
}

stopifnot(fuller_mlm_resolve_backend('auto') == 'r')
expect_error(fuller_mlm_resolve_backend('rcpp'), 'requires the compiled')
Rcpp::sourceCpp('src/fuller_mlm_kernels.cpp', cacheDir = file.path(tempdir(), 'fuller_mlm_rcpp_cache'))
stopifnot(fuller_mlm_resolve_backend('auto') == 'rcpp')

# Valid heterogeneous, correlated error blocks, independently constructed.
set.seed(241001)
make_data <- function(q, n = 180L) {
  h <- exp(rnorm(n, sd = 0.45))
  joint <- diag(c(0.12, rep(0.16, q))) + 0.025
  errors <- matrix(rnorm(n * (q + 1L)), n) %*% chol(joint) * sqrt(h)
  latent <- matrix(rnorm(n * q), n)
  x <- latent + errors[, -1L, drop = FALSE]
  colnames(x) <- paste0('x', seq_len(q))
  y <- 0.3 + as.vector(latent %*% seq(0.6, -0.25, length.out = q)) +
    rnorm(n, sd = 0.7) + errors[, 1L]
  list(y = y, x = x,
       omega_x = array(rep(joint[-1L, -1L], n) * rep(h, each = q * q), c(q, q, n)),
       omega_y = joint[1L, 1L] * h,
       omega_xy = outer(h, joint[1L, -1L]))
}

for (q in c(1L, 2L, 4L)) {
  dat <- make_data(q)
  r <- do.call(fit_fuller_mlm_variants, c(dat, list(backend = 'r')))
  cpp <- do.call(fit_fuller_mlm_variants, c(dat, list(backend = 'rcpp')))
  for (variant in names(r)) {
    compare(r[[variant]], cpp[[variant]])
    if (with_fsr) {
      d <- r[[variant]]$diagnostics
      ref <- do.call(fit_croon_fuller_q, c(dat, list(
        correction_mode = 'fuller', backend = 'r',
        preliminary_moment = d$preliminary_moment, variance_bread = d$variance_bread)))
      compare(cpp[[variant]], ref)
    }
  }
  # Production exposes only the last predictor's coefficient and SE.
  if (q <= 2L) {
    df <- data.frame(y = dat$y, x1 = dat$x[, q], s22 = dat$omega_x[q, q, ],
                     syy = dat$omega_y, c1 = dat$omega_xy[, q])
    args <- list(stage2_df = df, outcome = 'y', predictor_u1 = 'x1',
                 meas22 = 's22', outcome_meas_var = 'syy',
                 predictor_outcome_meas_cov_u1 = 'c1')
    if (q == 2L) {
      args$stage2_df$x0 <- dat$x[, 1L]
      args$stage2_df$s11 <- dat$omega_x[1L, 1L, ]
      args$stage2_df$s12 <- dat$omega_x[1L, 2L, ]
      args$stage2_df$c0 <- dat$omega_xy[, 1L]
      args$predictor_u0 <- 'x0'
      args$meas11 <- 's11'
      args$meas12 <- 's12'
      args$predictor_outcome_meas_cov_u0 <- 'c0'
    }
    production <- do.call(fit_fuller_dual_variants, args)
    stopifnot(all(production$status_code == 0L))
    for (i in seq_len(nrow(production))) {
      fit <- cpp[[production$fuller_variant[i]]]
      close(fit$coefficients$estimate[q + 1L], production$estimate[i], 'Production estimate')
      close(fit$coefficients$se[q + 1L], production$se[i], 'Production SE')
      fields <- c('lambda1', 'lambda2', 'sigma2', 'correction1', 'correction3')
      legacy <- c('fuller_lambda1', 'fuller_lambda2', 'fuller_sigma2',
                  'fuller_correction1', 'fuller_correction_c')
      close(unlist(fit$diagnostics[fields]), unlist(production[i, legacy]), 'Production steps')
    }
  }
}
cat('R/Rcpp, q = 1/2/4, all variants, and production parity passed.\n')

# Covariance formats, alpha controls, observed covariates, and the OLS limit.
dat <- make_data(3L)
dat$omega_x <- dat$omega_x[, , 1L]
dat$omega_xy <- c(0.01, -0.01, 0)
dat$omega_y <- 0.2
ref <- do.call(fit_fuller_mlm, dat)
arr <- array(rep(dat$omega_x, length(dat$y)), c(3L, 3L, length(dat$y)))
for (layout in list(arr, aperm(arr, c(3L, 1L, 2L)), rep(list(dat$omega_x), length(dat$y)))) {
  dat_layout <- dat
  dat_layout$omega_x <- layout
  for (backend in c('r', 'rcpp')) {
    compare(do.call(fit_fuller_mlm, c(dat_layout, list(backend = backend))), ref)
  }
}
custom <- do.call(fit_fuller_mlm, c(dat, list(alpha_step1 = 7, alpha_step3 = 8)))
custom_r <- do.call(fit_fuller_mlm, c(dat, list(alpha_step1 = 7, alpha_step3 = 8, backend = 'r')))
compare(custom, custom_r)
if (with_fsr) {
  custom_ref <- do.call(fit_vh_fuller_q, c(dat, list(alpha_step1 = 7, alpha_step3 = 8, backend = 'r')))
  compare(custom, custom_ref)
}
small_alpha <- do.call(fit_fuller_mlm, c(dat, list(alpha_step1 = 1, alpha_step3 = 1)))
compare(small_alpha, ref)

observed <- dat
observed$omega_x[3L, ] <- observed$omega_x[, 3L] <- 0
observed$omega_xy[3L] <- 0
compare(do.call(fit_fuller_mlm, c(observed, list(backend = 'r'))),
        do.call(fit_fuller_mlm, c(observed, list(backend = 'rcpp'))))
ols <- lm(dat$y ~ dat$x)
zero <- fit_fuller_mlm(dat$y, dat$x, matrix(0, 3L, 3L))
close(zero$coefficients$estimate, coef(ols), 'OLS coefficients')
close(zero$vcov, vcov(ols), 'OLS covariance')

# Force the small determinant-root branch and zero structural-variance boundary.
weak <- dat
weak$omega_x <- diag(3, 3L)
weak$omega_xy <- NULL
for (backend in c('r', 'rcpp')) {
  fits <- do.call(fit_fuller_mlm_variants, c(weak, list(backend = backend)))
  for (fit in fits) {
    stopifnot(fit$converged, fit$diagnostics$lambda1 < 1, fit$diagnostics$sigma2 == 0)
    reference_fitter <- if (with_fsr) fit_vh_fuller_q else fit_fuller_mlm
    ref_weak <- do.call(reference_fitter, c(weak, list(
      preliminary_moment = fit$diagnostics$preliminary_moment,
      variance_bread = fit$diagnostics$variance_bread, backend = 'r')))
    compare(fit, ref_weak)
  }
}
# Variant selection preserves order and removes duplicates.
selected <- do.call(fit_fuller_mlm_variants, c(dat, list(
  variants = c('fuller_equations', 'stabilized', 'fuller_equations'))))
stopifnot(identical(names(selected), c('fuller_equations', 'stabilized')))
cat('Layouts, alpha controls, observed covariates, OLS, and weak-root branch passed.\n')

for (backend in c('r', 'rcpp')) {
  fit <- function(...) fit_fuller_mlm(dat$y, dat$x, ..., backend = backend)
  expect_error(fit(diag(c(-1, 1, 1))), 'positive semidefinite')
  expect_error(fit(diag(3), omega_xy = c(1, 0, 0)), 'not PSD')
  expect_error(fit(diag(3), omega_xy = c(NA, 0, 0)), 'non-finite')
  expect_error(fit(diag(3), omega_y = -1), 'nonnegative')
  expect_error(fit(diag(3), alpha_step1 = c(1, 2)), 'positive scalars')
  expect_error(fit(diag(3), alpha_step3 = NA), 'positive scalars')
  expect_error(fit(matrix(0, 2, 2)), 'q by q')
  expect_error(fit_fuller_mlm(c(NA, dat$y[-1]), dat$x, diag(3), backend = backend), 'finite')
  # Exactly duplicated error-free columns produce a reported solve failure.
  singular <- fit_fuller_mlm(dat$y, cbind(dat$x[, 1], dat$x[, 1]),
                             matrix(0, 2, 2), backend = backend)
  stopifnot(!singular$converged, singular$status_code == 1L,
            grepl('singular', singular$message))
  no_variance <- fit_fuller_mlm(rep(0, length(dat$y)), dat$x, matrix(0, 3, 3),
                               backend = backend)
  stopifnot(!no_variance$converged, grepl('weights', no_variance$message))
}
expect_error(do.call(fit_fuller_mlm_variants, c(dat, list(variants = NA_character_))), 'variants')
expect_error(do.call(fit_fuller_mlm_variants, c(dat, list(variants = character()))), 'variants')

# Package-like lexical lookup: neither R functions nor exports live globally.
private <- new.env(parent = baseenv())
sys.source('R/fuller_mlm.R', envir = private)
Rcpp::sourceCpp('src/fuller_mlm_kernels.cpp', env = private,
                cacheDir = file.path(tempdir(), 'fuller_mlm_rcpp_cache'))
stopifnot(private$fuller_mlm_resolve_backend('auto') == 'rcpp')
compare(do.call(private$fit_fuller_mlm, dat), do.call(fit_fuller_mlm, dat))
cat('Validation, numerical failures, and private-environment dispatch passed.\n')

if (with_fsr) cat('Optional vh_fsr reference comparisons passed.\n')
