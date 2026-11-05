#' Shared utility helpers for MLM BLUP simulations.
#'
#' These helpers keep simulation scripts concise and consistent. They mostly
#' provide NA-stable summaries, compact diagnostic strings, failure-tolerant
#' mixed-model fitting, and matrix regularization for corrected-score
#' estimators.

#' Compute a mean that is stable for empty or all-missing inputs.
#'
#' @details
#' Base `mean(..., na.rm = TRUE)` returns `NaN` when all values are missing.
#' Simulation summaries are easier to combine when such cases are represented
#' as `NA_real_`, so this helper explicitly handles zero-length and all-missing
#' vectors before delegating to `mean()`.
#'
#' @param x Numeric, logical, or otherwise mean-compatible vector.
#'
#' @return
#' A numeric scalar. Returns `NA_real_` when `x` has length zero or contains no
#' non-missing values; otherwise returns `mean(x, na.rm = TRUE)`.
safe_mean <- function(x) {
  # Preserve a single missing numeric value rather than allowing mean() to emit
  # NaN, which is awkward in downstream summary tables and CSV outputs.
  if (length(x) == 0L || all(is.na(x))) {
    return(NA_real_)
  }
  mean(x, na.rm = TRUE)
}

#' Collapse diagnostic messages into a compact single-line string.
#'
#' @details
#' Simulation runs often accumulate warning, convergence, or model-status
#' messages from several sources. This helper drops missing entries, converts
#' the remaining values to character, removes duplicates while preserving their
#' first-observed order, joins them with `" | "`, and normalizes whitespace.
#'
#' @param x Vector of messages, conditions, or values coercible to character.
#'
#' @return
#' A character scalar containing the compacted message, or `NA_character_` when
#' no non-missing message is available.
compact_message <- function(x) {
  # Treat NULL, empty vectors, and all-missing vectors as "no diagnostic detail"
  # so callers can store a clean missing value.
  if (is.null(x) || length(x) == 0L || all(is.na(x))) {
    return(NA_character_)
  }

  # Coerce after omitting missing values so condition messages, numeric status
  # codes, and character warnings can all be summarized by the same helper.
  x <- unique(stats::na.omit(as.character(x)))
  if (length(x) == 0L) {
    return(NA_character_)
  }

  # Keep diagnostics readable in one table cell or CSV field by removing
  # embedded line breaks and repeated whitespace.
  gsub("[[:space:]]+", " ", trimws(paste(x, collapse = " | ")))
}

#' Fit an lme4 mixed model with quiet failure handling.
#'
#' @details
#' Wraps `lme4::lmer()` for simulation loops where occasional convergence
#' warnings, singular fits, or fitting errors are expected. Warnings and
#' messages are suppressed so large simulation grids do not flood the console.
#' Hard errors return `NULL`, allowing callers to record a failed replicate and
#' continue.
#'
#' This helper does not hide successful fits that are singular or converged with
#' warnings; those conditions should be diagnosed separately with helpers such
#' as `get_stage1_diagnostics()`.
#'
#' @param formula Model formula passed to `lme4::lmer()`.
#' @param data Data frame used for model fitting.
#' @param ... Additional arguments passed through to `lme4::lmer()`.
#'
#' @return
#' An `lmerMod` object on successful fitting, or `NULL` if `lmer()` throws an
#' error.
safe_lmer <- function(formula, data, ...) {
  tryCatch(
    # Suppress routine optimizer chatter here; convergence and singularity are
    # captured explicitly downstream when a fit object is available.
    suppressWarnings(suppressMessages(lme4::lmer(formula, data = data, ...))),
    error = function(e) NULL
  )
}

#' Fit an nlme mixed model with quiet failure handling.
#'
#' @details
#' Mirrors `safe_lmer()` for R-side residual covariance models. This wrapper is
#' used when simulations need a Stage-1 fit whose fixed effects, random-effect
#' covariance, and residual covariance all come from the same non-iid likelihood
#' rather than from an `lme4` iid-residual approximation.
#'
#' @param fixed Fixed-effect formula passed to `nlme::lme()`.
#' @param random Random-effect formula passed to `nlme::lme()`.
#' @param data Data frame used for model fitting.
#' @param ... Additional arguments passed through to `nlme::lme()`.
#'
#' @return
#' An `lme` object on successful fitting, or `NULL` if `nlme::lme()` throws an
#' error or if `nlme` is unavailable.
safe_lme <- function(fixed, random, data, ...) {
  if (!requireNamespace("nlme", quietly = TRUE)) {
    return(NULL)
  }
  tryCatch(
    suppressWarnings(suppressMessages(nlme::lme(
      fixed = fixed,
      random = random,
      data = data,
      ...
    ))),
    error = function(e) NULL
  )
}

#' Project a symmetric matrix to positive definite by flooring eigenvalues.
#'
#' @details
#' Corrected-score and errors-in-variables estimators can produce covariance or
#' normal-equation matrices with tiny negative eigenvalues because of sampling
#' noise or numerical roundoff. This helper symmetrizes the input matrix,
#' performs an eigen decomposition, replaces eigenvalues below `min_eigen` with
#' `min_eigen`, reconstructs the matrix, and symmetrizes again to remove
#' floating-point asymmetry.
#'
#' @param mat Numeric square matrix to regularize.
#' @param min_eigen Numeric lower bound for the eigenvalues of the returned
#' matrix. Defaults to `1e-6`.
#'
#' @return
#' A symmetric positive-definite matrix with eigenvalues no smaller than
#' `min_eigen`, up to numerical precision.
project_to_pd <- function(mat, min_eigen = 1e-6) {
  # Force exact symmetry before decomposition; eigen(..., symmetric = TRUE)
  # only consults one triangle of the matrix.
  mat <- (mat + t(mat)) / 2
  eig <- eigen(mat, symmetric = TRUE)

  # Flooring, rather than adding a fixed ridge to the whole matrix, changes
  # only the directions that are non-positive or too close to singular.
  eig$values <- pmax(eig$values, min_eigen)
  out <- eig$vectors %*% diag(eig$values, nrow = nrow(mat)) %*% t(eig$vectors)

  # Reconstruction can introduce tiny asymmetric roundoff; return a matrix that
  # downstream symmetric solvers and diagnostics can safely treat as symmetric.
  (out + t(out)) / 2
}
