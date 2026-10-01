# Fuller/M-MOM stage-2 regression for unshrunk MLM random-effect proxies.
# Extracted from vh_fsr/R/vh_fuller_q.R; no factor-score or Croon dependencies.
# Row-specific full error covariances are intrinsic to unbalanced MLMs.
# Small matrix operations use R's BLAS/LAPACK; optional C++ kernels handle rows.

#' Fit Fuller/M-MOM regression to unshrunk MLM proxies.
#'
#' Each row represents an independent Stage-2 unit (usually an MLM cluster).
#' The additive measurement model is x_i = xi_i + u_xi and y_i = eta_i + u_yi,
#' with eta_i = intercept + beta' xi_i + residual_i. Supply the error covariance
#' of these proxies after unshrinking, not the posterior covariance of the BLUPs.
#'
#' @param y Numeric outcome vector: observed or already unshrunk.
#' @param x Numeric n by q predictor matrix, without an intercept column.
#'   Unshrink BLUPs before fitting. Observed covariates may be included with
#'   zero rows and columns in their measurement-error covariance blocks.
#' @param omega_x Additive predictor error covariance after unshrinking:
#'   a common q by q matrix, q by q by n array, n by q by q array, or length-n
#'   list of matrices. The canonical q by q by n layout takes precedence.
#'   Predictor order must match the columns of x; names are not used to reorder.
#'   Blocks are symmetrized and checked for positive semidefiniteness, as is
#'   the joint outcome/predictor error block (relative tolerance 1e-8).
#' @param omega_y Outcome error variance: nonnegative scalar or length n.
#' @param omega_xy Predictor/outcome error covariance: length q or n by q.
#'   NULL means zero, appropriate only when these errors are uncorrelated.
#' @param alpha_step1,alpha_step3 Positive scalar Fuller correction constants;
#'   default to q + 2, with smaller supplied values raised to q + 2, as in
#'   fit_fuller_dual_core().
#' @param preliminary_moment "modified" subtracts the Step-1 lambda/alpha
#'   correction times the error moments; "fuller" subtracts the full error
#'   moments. This changes the preliminary coefficients and feasible weights.
#' @param variance_bread "modified" uses the final lambda/alpha-corrected
#'   predictor moment; "fuller" uses the fully subtracted weighted moment.
#'   This changes SEs only; the final point equation is modified in both cases.
#' @param backend "auto" uses registered Rcpp kernels when available, otherwise
#'   R; "r" and "rcpp" require the selected backend. No runtime compilation.
#'
#' @details Estimates are on the input scale. No latent-SD standardization is
#'   performed. Standard errors use Fuller's model-based covariance conditional
#'   on the supplied Stage-1 quantities; Stage-1 uncertainty is not propagated.
#'   Inputs must be finite and aligned; this API does not silently omit rows or
#'   clamp negative input variances. Invalid inputs raise errors; numerical
#'   fitting failures return converged = FALSE, status_code = 1, and a message.
#'   At least max(8, q + 2) complete rows are required. The intercept is added
#'   internally and is measured without error. R and Rcpp use the same equations.
#' @return A fuller_mlm_fit list. On success, status_code is 0 and coefficients
#'   contains term, estimate, se, ci_low, and ci_high (95% normal intervals).
#'   vcov is the full coefficient covariance, including the intercept;
#'   preliminary_coefficients is the Step-1 solution; corrected_design_moment
#'   is the unnormalized final predictor moment. analysis_weights contains
#'   inverse case variances, 1 / (sigma2 + composite measurement variance).
#'   diagnostics records n, q, backend, determinant roots, correction constants,
#'   structural residual variance sigma2, and predictor-moment diagnostics.
#'   c_F aliases correction3. Failed fits have NULL coefficients and vcov.
fit_fuller_mlm <- function(
    y, x, omega_x, omega_y = 0, omega_xy = NULL,
    alpha_step1 = NULL, alpha_step3 = NULL,
    preliminary_moment = c("modified", "fuller"),
    variance_bread = c("modified", "fuller"),
    backend = c("auto", "r", "rcpp")) {
  preliminary_moment <- match.arg(preliminary_moment)
  variance_bread <- match.arg(variance_bread)
  inputs <- fuller_mlm_prepare(
    y, x, omega_x, omega_y, omega_xy, alpha_step1, alpha_step3, backend
  )
  fuller_mlm_fit_prepared(inputs, preliminary_moment, variance_bread, match.call())
}

#' Fit the same three algebra variants as fit_fuller_dual_variants(), for any q.
#'
#' @inheritParams fit_fuller_mlm
#' @param variants Any subset of "stabilized", "fuller_preliminary", and
#'   "fuller_equations". Duplicates are removed, retaining requested order.
#' @details "stabilized" uses modified preliminary and bread moments;
#'   "fuller_preliminary" uses a fully subtracted preliminary moment and modified
#'   bread; "fuller_equations" fully subtracts both. All three retain the final
#'   lambda/alpha-modified point equation and return estimates on the input scale.
#' @return Named list of fuller_mlm_fit objects, one per variant. Inputs and
#'   covariance validation are shared across variants. Inspect converged on
#'   each fit before using its coefficient table.
fit_fuller_mlm_variants <- function(
    y, x, omega_x, omega_y = 0, omega_xy = NULL,
    variants = c("stabilized", "fuller_preliminary", "fuller_equations"),
    alpha_step1 = NULL, alpha_step3 = NULL,
    backend = c("auto", "r", "rcpp")) {
  allowed <- c("stabilized", "fuller_preliminary", "fuller_equations")
  if (!is.character(variants) || length(variants) == 0L ||
      anyNA(variants) || any(!variants %in% allowed)) {
    stop("`variants` must contain only: ", paste(allowed, collapse = ", "))
  }
  variants <- unique(variants)
  inputs <- fuller_mlm_prepare(
    y, x, omega_x, omega_y, omega_xy, alpha_step1, alpha_step3, backend
  )
  call <- match.call()
  fits <- lapply(variants, function(variant) {
    fit <- fuller_mlm_fit_prepared(
      inputs,
      preliminary_moment = if (variant == "stabilized") "modified" else "fuller",
      variance_bread = if (variant == "fuller_equations") "fuller" else "modified",
      call = call
    )
    fit$variant <- variant
    fit
  })
  names(fits) <- variants
  fits
}

print.fuller_mlm_fit <- function(x, ...) {
  cat("Fuller/M-MOM MLM fit\n")
  if (!isTRUE(x$converged)) {
    cat("Fit failed:", x$message, "\n")
  } else {
    print(x$coefficients, row.names = FALSE, ...)
    cat("Model-based Fuller SEs; Stage-1 quantities held fixed.\n")
  }
  invisible(x)
}

fuller_mlm_symmetrize <- function(x) {
  (x + t(x)) / 2
}

fuller_mlm_matrix_diagnostics <- function(x) {
  x <- fuller_mlm_symmetrize(as.matrix(x))
  values <- tryCatch(
    eigen(x, symmetric = TRUE, only.values = TRUE)$values,
    error = function(e) NA_real_
  )
  values <- values[is.finite(values)]
  if (length(values) == 0L) {
    return(list(min_eigen = NA_real_, max_eigen = NA_real_, condition_number = Inf))
  }
  min_eigen <- min(values)
  max_eigen <- max(values)
  condition_number <- if (min_eigen > 0) max_eigen / min_eigen else Inf
  list(
    min_eigen = min_eigen,
    max_eigen = max_eigen,
    condition_number = condition_number
  )
}

fuller_mlm_relative_min_eigen <- function(corrected, observed) {
  if (is.finite(corrected$min_eigen) &&
      is.finite(observed$max_eigen) &&
      observed$max_eigen > sqrt(.Machine$double.eps)) {
    corrected$min_eigen / observed$max_eigen
  } else {
    NA_real_
  }
}

fuller_mlm_smallest_det_root <- function(a, b) {
  if (!requireNamespace("geigen", quietly = TRUE)) {
    stop("The `geigen` package is required for Fuller/M-MOM estimation.")
  }
  # Roots solve det(a - lambda * b) = 0. The error-free intercept makes b
  # singular, so the solver must permit positive-semidefinite b; its symmetric
  # mode requires positive definiteness. Ignore infinite roots and numerical
  # imaginary noise. NA means no usable finite root and selects the alpha-only
  # correction branch in the caller, matching the established R estimator.
  out <- tryCatch(
    geigen::geigen(a, b, symmetric = FALSE, only.values = TRUE),
    error = function(e) NULL
  )
  if (is.null(out) || is.null(out$values)) return(NA_real_)

  values <- out$values
  if (is.complex(values)) {
    scale <- suppressWarnings(max(abs(values), na.rm = TRUE))
    if (!is.finite(scale)) scale <- 1
    tolerance <- sqrt(.Machine$double.eps) * max(1, scale)
    values <- ifelse(abs(Im(values)) <= tolerance, Re(values), NA_real_)
  }
  values <- values[is.finite(values)]
  if (length(values) == 0L) NA_real_ else min(values)
}

fuller_mlm_coerce_omega_x_array <- function(omega_x, n, q) {
  if (is.matrix(omega_x)) {
    if (!identical(dim(omega_x), c(q, q))) {
      stop("A common `omega_x` matrix must be q by q.")
    }
    out <- array(rep(as.numeric(omega_x), n), dim = c(q, q, n))
  } else if (is.array(omega_x) && length(dim(omega_x)) == 3L) {
    dims <- dim(omega_x)
    if (identical(dims, c(q, q, n))) {
      out <- omega_x
    } else if (identical(dims, c(n, q, q))) {
      out <- aperm(omega_x, c(2, 3, 1))
    } else {
      stop("An `omega_x` array must be q by q by n (or n by q by q).")
    }
  } else if (is.list(omega_x)) {
    if (length(omega_x) != n) stop("An `omega_x` list must have length n.")
    out <- array(NA_real_, dim = c(q, q, n))
    for (i in seq_len(n)) {
      this <- as.matrix(omega_x[[i]])
      if (!identical(dim(this), c(q, q))) {
        stop("Every element of an `omega_x` list must be q by q.")
      }
      out[, , i] <- this
    }
  } else {
    stop("`omega_x` must be a common matrix, a three-dimensional array, or a list of matrices.")
  }

  if (any(!is.finite(out))) stop("`omega_x` contains non-finite values.")
  out
}

fuller_mlm_normalize_omega_x <- function(omega_x, n, q, tolerance = 1e-8) {
  out <- fuller_mlm_coerce_omega_x_array(omega_x, n = n, q = q)
  for (i in seq_len(n)) {
    out[, , i] <- fuller_mlm_symmetrize(out[, , i])
    diag_i <- fuller_mlm_matrix_diagnostics(out[, , i])
    scale_i <- max(1, abs(diag_i$max_eigen))
    if (!is.finite(diag_i$min_eigen) || diag_i$min_eigen < -tolerance * scale_i) {
      stop(sprintf("`omega_x` is not positive semidefinite for row %d.", i))
    }
  }
  out
}

fuller_mlm_normalize_omega_xy <- function(omega_xy, n, q) {
  if (is.null(omega_xy)) return(matrix(0, nrow = n, ncol = q))
  if (any(!is.finite(omega_xy))) stop("`omega_xy` contains non-finite values.")
  if (is.vector(omega_xy) && length(omega_xy) == q) {
    return(matrix(rep(as.numeric(omega_xy), each = n), nrow = n, ncol = q))
  }
  omega_xy <- as.matrix(omega_xy)
  if (!identical(dim(omega_xy), c(n, q))) {
    stop("`omega_xy` must be length q or an n by q matrix.")
  }
  if (any(!is.finite(omega_xy))) stop("`omega_xy` contains non-finite values.")
  omega_xy
}

fuller_mlm_weighted_omega_sum <- function(omega_x, weights) {
  q <- dim(omega_x)[1]
  n <- dim(omega_x)[3]
  out <- matrix(0, nrow = q, ncol = q)
  for (i in seq_len(n)) out <- out + weights[[i]] * omega_x[, , i]
  fuller_mlm_symmetrize(out)
}

fuller_mlm_rcpp_function_names <- function() {
  c(
    "fuller_mlm_normalize_omega_x_cpp",
    "fuller_mlm_joint_omega_bad_row_cpp",
    "fuller_mlm_weighted_omega_sum_cpp",
    "fuller_mlm_omega_matvec_cpp",
    "fuller_mlm_tilde_crossproduct_cpp"
  )
}

fuller_mlm_rcpp_available <- function() {
  all(vapply(
    fuller_mlm_rcpp_function_names(),
    exists,
    logical(1),
    mode = "function",
    envir = environment(fuller_mlm_rcpp_available),
    inherits = TRUE
  ))
}

fuller_mlm_resolve_backend <- function(backend = c("auto", "r", "rcpp")) {
  if (length(backend) > 1L) backend <- backend[[1L]]
  backend <- match.arg(tolower(backend), c("auto", "r", "rcpp"))
  if (identical(backend, "auto")) {
    return(if (fuller_mlm_rcpp_available()) "rcpp" else "r")
  }
  if (identical(backend, "rcpp") && !fuller_mlm_rcpp_available()) {
    stop(
      "`backend = \"rcpp\"` requires the compiled Fuller/M-MOM helpers. ",
      "For standalone development, call `Rcpp::sourceCpp(",
      "\"src/fuller_mlm_kernels.cpp\")`; a package build should ",
      "generate and register the Rcpp exports."
    )
  }
  backend
}

fuller_mlm_normalize_omega_x_backend <- function(
    omega_x, n, q, tolerance, backend) {
  if (identical(backend, "r")) {
    return(fuller_mlm_normalize_omega_x(omega_x, n, q, tolerance))
  }
  if (is.matrix(omega_x)) {
    if (!identical(dim(omega_x), c(q, q))) {
      stop("A common `omega_x` matrix must be q by q.")
    }
    if (any(!is.finite(omega_x))) stop("`omega_x` contains non-finite values.")
    storage.mode(omega_x) <- "double"
  } else {
    omega_x <- fuller_mlm_coerce_omega_x_array(omega_x, n = n, q = q)
    storage.mode(omega_x) <- "double"
  }
  normalized <- fuller_mlm_normalize_omega_x_cpp(
    omega_x, n = n, q = q, tolerance = tolerance
  )
  if (normalized$bad_row > 0L) {
    stop(sprintf(
      "`omega_x` is not positive semidefinite for row %d.",
      normalized$bad_row
    ))
  }
  normalized$omega_x
}

fuller_mlm_validate_joint_omega <- function(
    omega_x, omega_y, omega_xy, tolerance, backend) {
  n <- dim(omega_x)[3L]
  if (identical(backend, "rcpp")) {
    bad_row <- fuller_mlm_joint_omega_bad_row_cpp(
      omega_x, omega_y, omega_xy, tolerance = tolerance
    )
    if (bad_row > 0L) {
      stop(sprintf(
        "The joint measurement-error covariance is not PSD for row %d.",
        bad_row
      ))
    }
    return(invisible(NULL))
  }

  for (i in seq_len(n)) {
    full_i <- rbind(
      c(omega_y[[i]], omega_xy[i, ]),
      cbind(omega_xy[i, ], omega_x[, , i])
    )
    diag_i <- fuller_mlm_matrix_diagnostics(full_i)
    scale_i <- max(1, abs(diag_i$max_eigen))
    if (!is.finite(diag_i$min_eigen) ||
        diag_i$min_eigen < -tolerance * scale_i) {
      stop(sprintf(
        "The joint measurement-error covariance is not PSD for row %d.", i
      ))
    }
  }
  invisible(NULL)
}

fuller_mlm_weighted_omega_sum_backend <- function(omega_x, weights, backend) {
  if (identical(backend, "rcpp")) {
    return(fuller_mlm_symmetrize(
      fuller_mlm_weighted_omega_sum_cpp(omega_x, weights)
    ))
  }
  fuller_mlm_weighted_omega_sum(omega_x, weights)
}

fuller_mlm_composite_measurement_variance <- function(
    omega_y, omega_xy, omega_x, gamma_x, backend) {
  n <- length(omega_y)
  if (identical(backend, "rcpp")) {
    omega_gamma <- fuller_mlm_omega_matvec_cpp(omega_x, gamma_x)
    return(
      omega_y - 2 * as.vector(omega_xy %*% gamma_x) +
        as.vector(omega_gamma %*% gamma_x)
    )
  }

  out <- numeric(n)
  for (i in seq_len(n)) {
    out[[i]] <- omega_y[[i]] -
      2 * sum(gamma_x * omega_xy[i, ]) +
      as.numeric(crossprod(gamma_x, omega_x[, , i] %*% gamma_x))
  }
  out
}

fuller_mlm_failure <- function(call, message, diagnostics = list()) {
  structure(
    list(
      call = call,
      converged = FALSE,
      status_code = 1L,
      message = message,
      coefficients = NULL,
      vcov = NULL,
      diagnostics = diagnostics
    ),
    class = "fuller_mlm_fit"
  )
}

# Validate once, including when fitting several algebra variants.
fuller_mlm_prepare <- function(y, x, omega_x, omega_y, omega_xy,
                               alpha_step1, alpha_step3, backend) {
  backend <- fuller_mlm_resolve_backend(backend)
  if (!requireNamespace("geigen", quietly = TRUE)) {
    stop("The `geigen` package is required for Fuller/M-MOM estimation.")
  }
  if (!is.numeric(y) || !is.numeric(as.matrix(x))) {
    stop("`y` and `x` must be numeric.")
  }
  y <- as.numeric(y)
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  n <- length(y)
  q <- ncol(x)
  p <- q + 1L
  if (nrow(x) != n || q < 1L) stop("`x` must be an n by q matrix aligned with `y`.")
  if (n < max(8L, p + 1L)) stop("Fuller/M-MOM requires at least max(8, q + 2) complete rows.")
  if (any(!is.finite(y)) || any(!is.finite(x))) {
    stop("`y` and `x` must be complete and finite; filter rows before fitting.")
  }
  term_names <- colnames(x)
  if (is.null(term_names)) term_names <- paste0("x", seq_len(q))
  if (anyNA(term_names) || any(!nzchar(term_names)) ||
      anyDuplicated(c("(Intercept)", term_names))) {
    stop("Predictor names must be nonempty, unique, and distinct from `(Intercept)`.")
  }
  colnames(x) <- term_names
  design <- cbind(`(Intercept)` = 1, x)
  design_names <- colnames(design)
  predictor_index <- seq.int(2L, p)
  omega_x <- fuller_mlm_normalize_omega_x_backend(
    omega_x, n, q, tolerance = 1e-8, backend = backend
  )
  if (!is.numeric(omega_y) || !length(omega_y) %in% c(1L, n) ||
      any(!is.finite(omega_y)) || any(omega_y < 0)) {
    stop("`omega_y` must be a finite, nonnegative scalar or length-n vector.")
  }
  omega_y <- rep(as.numeric(omega_y), length.out = n)
  if (!is.null(omega_xy) && !is.numeric(omega_xy)) stop("`omega_xy` must be numeric.")
  omega_xy <- fuller_mlm_normalize_omega_xy(omega_xy, n, q)
  fuller_mlm_validate_joint_omega(omega_x, omega_y, omega_xy, 1e-8, backend)
  normalize_alpha <- function(value) {
    if (is.null(value)) return(p + 1)
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) || value <= 0) {
      stop("Fuller alpha constants must be finite, positive scalars.")
    }
    max(value, p + 1)
  }
  alpha_step1 <- normalize_alpha(alpha_step1)
  alpha_step3 <- normalize_alpha(alpha_step3)
  list(y = y, x = x, n = n, q = q, p = p, design = design,
       design_names = design_names, predictor_index = predictor_index,
       omega_x = omega_x, omega_y = omega_y, omega_xy = omega_xy,
       alpha_step1 = alpha_step1, alpha_step3 = alpha_step3, backend = backend)
}

fuller_mlm_fit_prepared <- function(inputs, preliminary_moment, variance_bread, call) {
  # The environment contains only validated inputs, not arbitrary user names.
  list2env(inputs, envir = environment())
  omega_x_sum_q <- fuller_mlm_weighted_omega_sum_backend(
    omega_x, rep(1, n), backend
  )
  omega_x_sum <- matrix(0, nrow = p, ncol = p)
  omega_x_sum[predictor_index, predictor_index] <- omega_x_sum_q
  omega_xy_design <- cbind(0, omega_xy)
  omega_xy_sum <- colSums(omega_xy_design)

  # Step 1: stack (y, intercept, x) and its additive error block. The intercept
  # error row/column stays zero; omega_xy supplies the cross-error correction,
  # independently of the outcome error variance omega_y.
  b_matrix <- cbind(y, design)
  bb_sum <- crossprod(b_matrix)
  omega_sum <- matrix(0, nrow = p + 1L, ncol = p + 1L)
  omega_sum[1, 1] <- sum(omega_y)
  omega_sum[1, 2:(p + 1L)] <- omega_xy_sum
  omega_sum[2:(p + 1L), 1] <- omega_xy_sum
  omega_sum[2:(p + 1L), 2:(p + 1L)] <- omega_x_sum

  lambda1 <- fuller_mlm_smallest_det_root(bb_sum, omega_sum)
  correction1 <- if (is.finite(lambda1) && lambda1 <= 1 + 1 / n) {
    lambda1 - 1 / n - alpha_step1 / n
  } else {
    1 - alpha_step1 / n
  }

  preliminary_correction <- if (identical(preliminary_moment, "fuller")) {
    1
  } else {
    correction1
  }
  a0 <- crossprod(design) - preliminary_correction * omega_x_sum
  b0 <- as.vector(
    crossprod(design, y) -
      preliminary_correction * omega_xy_sum
  )
  gamma0 <- tryCatch(as.vector(solve(a0, b0)), error = function(e) NULL)
  if (is.null(gamma0) || any(!is.finite(gamma0))) {
    return(fuller_mlm_failure(
      call, "The preliminary corrected moment was singular.",
      diagnostics = list(
        computation_backend = backend,
        c_F = correction1
      )
    ))
  }
  names(gamma0) <- design_names
  gamma0_x <- gamma0[predictor_index]

  # Step 2: Var(u_y - gamma0_x' u_x) = omega_y - 2 gamma0_x' omega_xy
  # + gamma0_x' Omega_x gamma0_x. Subtract its mean from the preliminary
  # residual mean square; lambda1 < 1 sets the structural variance to zero.
  composite_measurement_variance <- fuller_mlm_composite_measurement_variance(
    omega_y, omega_xy, omega_x, gamma0_x, backend
  )
  if (any(composite_measurement_variance < -1e-8)) {
    return(fuller_mlm_failure(
      call,
      "The supplied measurement model produced a negative composite error variance.",
      diagnostics = list(computation_backend = backend)
    ))
  }
  composite_measurement_variance <- pmax(composite_measurement_variance, 0)

  residual0 <- y - as.vector(design %*% gamma0)
  sigma2_ols <- sum(residual0^2) / max(1, n - p)
  sigma2 <- if (is.finite(lambda1) && lambda1 < 1) {
    0
  } else {
    sigma2_ols - mean(composite_measurement_variance)
  }
  if (!is.finite(sigma2)) {
    return(fuller_mlm_failure(
      call,
      "The corrected structural residual variance was non-finite.",
      diagnostics = list(computation_backend = backend)
    ))
  }
  sigma2 <- max(0, sigma2)

  case_variance <- sigma2 + composite_measurement_variance
  if (any(!is.finite(case_variance)) ||
      any(case_variance <= sqrt(.Machine$double.eps))) {
    return(fuller_mlm_failure(
      call, "Feasible case weights were nonpositive or non-finite.",
      diagnostics = list(computation_backend = backend)
    ))
  }
  inverse_case_variance <- 1 / case_variance
  final_weights <- inverse_case_variance

  # Step 3: weight each unit by its inverse composite residual variance and
  # recompute the determinant root. Every variant uses this modified final
  # point equation; the preliminary switch affects it through the weights.
  bw_sum <- crossprod(b_matrix * sqrt(final_weights))
  omega_x_sum_w_q <- fuller_mlm_weighted_omega_sum_backend(
    omega_x, final_weights, backend
  )
  omega_x_sum_w <- matrix(0, nrow = p, ncol = p)
  omega_x_sum_w[predictor_index, predictor_index] <- omega_x_sum_w_q
  omega_xy_sum_w <- colSums(omega_xy_design * final_weights)
  omega_sum_w <- matrix(0, nrow = p + 1L, ncol = p + 1L)
  omega_sum_w[1, 1] <- sum(final_weights * omega_y)
  omega_sum_w[1, 2:(p + 1L)] <- omega_xy_sum_w
  omega_sum_w[2:(p + 1L), 1] <- omega_xy_sum_w
  omega_sum_w[2:(p + 1L), 2:(p + 1L)] <- omega_x_sum_w

  lambda2 <- fuller_mlm_smallest_det_root(bw_sum, omega_sum_w)
  correction3 <- if (is.finite(lambda2) && lambda2 <= 1 + 1 / n) {
    lambda2 - 1 / n - alpha_step3 / n
  } else {
    1 - alpha_step3 / n
  }
  s_star <- bw_sum - correction3 * omega_sum_w
  s_x_star <- s_star[2:(p + 1L), 2:(p + 1L), drop = FALSE]
  s_xy_star <- s_star[2:(p + 1L), 1, drop = FALSE]
  gamma <- tryCatch(as.vector(solve(s_x_star, s_xy_star)), error = function(e) NULL)
  if (is.null(gamma) || any(!is.finite(gamma))) {
    return(fuller_mlm_failure(
      call, "The final corrected Fuller/M-MOM moment was singular.",
      diagnostics = list(
        computation_backend = backend,
        c_F = correction3
      )
    ))
  }
  names(gamma) <- design_names

  # Fuller's model-based covariance, conditional on Stage-1 estimates.
  # Changing the bread does not refit gamma. The tilde term uses gamma0_x
  # (the preliminary slopes), and its crossproduct squares the inverse case
  # weights. Dividing the bread by n and the covariance by n^2 cancels the
  # normalization, preserving the established R equations.
  xw_sum <- crossprod(design * sqrt(final_weights))
  variance_bread_matrix <- if (identical(variance_bread, "modified")) {
    s_x_star
  } else {
    xw_sum - omega_x_sum_w
  }
  bread_inverse <- tryCatch(solve(variance_bread_matrix / n), error = function(e) NULL)
  if (is.null(bread_inverse) || any(!is.finite(bread_inverse))) {
    return(fuller_mlm_failure(
      call, "The Fuller variance bread was singular.",
      diagnostics = list(computation_backend = backend)
    ))
  }
  if (identical(backend, "rcpp")) {
    tilde_crossproduct <- fuller_mlm_tilde_crossproduct_cpp(
      omega_x, omega_xy, gamma0_x, inverse_case_variance
    )
  } else {
    tilde <- matrix(0, nrow = n, ncol = p)
    for (i in seq_len(n)) {
      tilde[i, predictor_index] <- omega_xy[i, ] - omega_x[, , i] %*% gamma0_x
    }
    tilde_crossproduct <- crossprod(tilde * inverse_case_variance)
  }
  meat <- xw_sum + tilde_crossproduct
  vcov_gamma <- bread_inverse %*% meat %*% bread_inverse / n^2
  vcov_gamma <- fuller_mlm_symmetrize(vcov_gamma)
  dimnames(vcov_gamma) <- list(design_names, design_names)
  variances <- diag(vcov_gamma)
  if (any(!is.finite(variances)) || any(variances < 0)) {
    return(fuller_mlm_failure(
      call, "The Fuller/M-MOM coefficient variance was negative or non-finite.",
      diagnostics = list(computation_backend = backend)
    ))
  }

  se <- sqrt(variances)
  coefficient_table <- data.frame(
    term = design_names,
    estimate = unname(gamma),
    se = unname(se),
    ci_low = unname(gamma - stats::qnorm(0.975) * se),
    ci_high = unname(gamma + stats::qnorm(0.975) * se),
    row.names = NULL,
    check.names = FALSE
  )

  observed_initial <- fuller_mlm_matrix_diagnostics(
    (crossprod(design) / n)[predictor_index, predictor_index, drop = FALSE]
  )
  corrected_initial <- fuller_mlm_matrix_diagnostics(
    (a0 / n)[predictor_index, predictor_index, drop = FALSE]
  )
  observed_final <- fuller_mlm_matrix_diagnostics(
    xw_sum[predictor_index, predictor_index, drop = FALSE]
  )
  corrected_final <- fuller_mlm_matrix_diagnostics(
    s_x_star[predictor_index, predictor_index, drop = FALSE]
  )

  structure(
    list(
      call = call,
      converged = TRUE,
      status_code = 0L,
      message = "ok",
      coefficients = coefficient_table,
      vcov = vcov_gamma,
      preliminary_coefficients = gamma0,
      corrected_design_moment = s_x_star,
      analysis_weights = final_weights,
      diagnostics = list(
        n = n,
        q = q,
        computation_backend = backend,
        c_F = correction3,
        lambda1 = lambda1,
        lambda2 = lambda2,
        correction1 = correction1,
        correction3 = correction3,
        alpha_step1 = alpha_step1,
        alpha_step3 = alpha_step3,
        sigma2 = sigma2,
        case_variance_min = min(case_variance),
        case_variance_max = max(case_variance),
        preliminary_moment = preliminary_moment,
        variance_bread = variance_bread,
        initial_condition_number = corrected_initial$condition_number,
        initial_relative_min_eigen = fuller_mlm_relative_min_eigen(corrected_initial, observed_initial),
        final_condition_number = corrected_final$condition_number,
        final_relative_min_eigen = fuller_mlm_relative_min_eigen(corrected_final, observed_final),
        predictor_outcome_error_covariance_max_abs = max(abs(omega_xy))
      )
    ),
    class = "fuller_mlm_fit"
  )
}
