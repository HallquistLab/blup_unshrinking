#' Shared Lai-style EB measurement and OpenMx helpers.
#'
#' These helpers support both the standalone Lai replication simulations and
#' the random-slope sandwich coverage simulation. They centralize construction
#' of 2S-PA measurement-model inputs, OpenMx retry handling, OpenMx diagnostic
#' extraction, and the Lai 2S-PA/2S-PAA model wrappers.

#' Warn when a legacy Lai EB-input convenience wrapper is used.
#'
#' @param fun_name Character scalar naming the caller.
#'
#' @return Invisibly returns `NULL`.
warn_legacy_lai_input_helper <- function(fun_name) {
  warning(
    fun_name,
    "() is a legacy Lai EB-input convenience wrapper. Prefer calling ",
    "get_stage1_eb_components() directly, then select_lai_measurement_columns() ",
    "when OpenMx/Lai column names are needed. The direct path keeps EB/BLUPs, ",
    "posterior covariance, Lai lambda/theta, fitted G, and residual R handling ",
    "tied to one explicit Stage-1 fit, which is important for non-diagonal R.",
    call. = FALSE
  )
  invisible(NULL)
}

#' Build the random-effect design matrix for one cluster.
#'
#' @details
#' Lai-style two-stage path analysis treats empirical Bayes random effects as
#' fallible measures of latent cluster effects. The per-cluster measurement
#' model depends on the same random-effect design matrix used in the first-stage
#' mixed model. For random-intercept-only models this is a column of ones. For
#' random-intercept/random-slope models this is `cbind(1, within_var)`.
#'
#' @param cluster_df Data frame containing the level-1 observations for a single
#' cluster.
#' @param within_var Optional character scalar naming the within-cluster slope
#' variable. If `NULL`, an intercept-only design is returned.
#'
#' @return A numeric matrix with one row per observation and one column per
#' random-effect component.
default_re_design <- function(cluster_df, within_var = NULL) {
  if (is.null(within_var)) {
    matrix(1, nrow = nrow(cluster_df), ncol = 1L)
  } else {
    cbind(1, cluster_df[[within_var]])
  }
}

#' Format one cluster's EB measurement quantities as a tibble row.
#'
#' @details
#' `compute_eb_measurement_inputs()` computes EB predictions, posterior
#' covariance, reliability (`lambda`), and unreliability (`theta`) matrices.
#' This helper gives those matrix entries stable column names used by the
#' downstream OpenMx wrappers. Univariate outputs can be prefixed, which is used
#' in Study 3 where the `z` outcome has its own first-stage model.
#'
#' @param id Cluster identifier.
#' @param eb Numeric vector of empirical Bayes random-effect predictions.
#' @param post_vcov Posterior covariance matrix for `eb`.
#' @param lambda Reliability/loading matrix mapping latent effects to EB
#' measures.
#' @param theta Unreliability/residual covariance matrix for EB measures.
#' @param prefix Optional prefix applied to univariate output column names.
#'
#' @return A one-row tibble. Bivariate rows contain `u0_eb`, `u1_eb`,
#' `postvar11`, `postvar12`, `postvar22`, `lambda11` through `lambda22`, and
#' `theta11` through `theta22`. Univariate rows contain prefixed `u0_eb`,
#' `postvar11`, `lambda11`, and `theta11`.
make_eb_output_row <- function(id, eb, post_vcov, lambda, theta, prefix = "") {
  n_re <- length(eb)
  out <- tibble::tibble(id = id)

  if (n_re == 1L) {
    out[[paste0(prefix, "u0_eb")]] <- unname(eb[[1]])
    out[[paste0(prefix, "postvar11")]] <- post_vcov[1, 1]
    out[[paste0(prefix, "lambda11")]] <- lambda[1, 1]
    out[[paste0(prefix, "theta11")]] <- theta[1, 1]
    return(out)
  }

  if (n_re == 2L) {
    out$u0_eb <- unname(eb[[1]])
    out$u1_eb <- unname(eb[[2]])
    out$postvar11 <- post_vcov[1, 1]
    out$postvar12 <- post_vcov[1, 2]
    out$postvar22 <- post_vcov[2, 2]
    out$lambda11 <- lambda[1, 1]
    out$lambda12 <- lambda[1, 2]
    out$lambda21 <- lambda[2, 1]
    out$lambda22 <- lambda[2, 2]
    out$theta11 <- theta[1, 1]
    out$theta12 <- theta[1, 2]
    out$theta22 <- theta[2, 2]
    return(out)
  }

  stop("Only univariate and bivariate random-effect measurement inputs are currently supported.")
}

select_lai_measurement_columns <- function(stage1_components, n_re, prefix = "") {
  if (n_re == 1L) {
    out <- stage1_components[, c("id", "u0_eb", "postvar11", "lambda11", "theta11"), drop = FALSE]
    if (nzchar(prefix)) {
      names(out) <- sub("^(u0_eb|postvar11|lambda11|theta11)$", paste0(prefix, "\\1"), names(out))
    }
    return(out)
  }

  out <- stage1_components[, c(
    "id", "u0_eb", "u1_eb", "postvar11", "postvar12", "postvar22",
    "lambda11", "lambda12", "lambda21", "lambda22", "theta11", "theta12", "theta22"
  ), drop = FALSE]
  if (nzchar(prefix)) {
    names(out)[names(out) != "id"] <- paste0(
      prefix,
      names(out)[names(out) != "id"]
    )
  }
  out
}

#' Order a 2-by-2 loading matrix for OpenMx `unique.bivariate` paths.
#'
#' @details
#' The Stage-1 columns use ordinary row/column subscripts: `lambda12` is row 1,
#' column 2 of the loading matrix, so the natural row-major listing is
#' `lambda11, lambda12, lambda21, lambda22`. In contrast,
#' `mxPath(..., connect = "unique.bivariate")` fills the RAM `A` matrix in
#' column-major order. It therefore requires
#' `lambda11, lambda21, lambda12, lambda22`. Supplying row-major order silently
#' transposes the two off-diagonal loadings.
#'
#' All bivariate 2S-PA and 2S-PAA wrappers use this helper for both fixed values
#' and definition-variable labels so the OpenMx-specific ordering is explicit
#' and cannot drift between model variants.
#'
#' @param prefix Optional column-name prefix, such as `"q_"` for the outcome
#'   process in the dual-process model.
#'
#' @return Character vector of loading-column names in the column-major order
#' required by OpenMx.
openmx_bivariate_loading_columns <- function(prefix = "") {
  paste0(prefix, c("lambda11", "lambda21", "lambda12", "lambda22"))
}

#' Compute Lai-style EB measurement-model inputs from a Stage-1 fit.
#'
#' @details
#' For each cluster this constructs the empirical Bayes random-effect
#' prediction, posterior covariance, and the Lai/OpenMx measurement-model
#' reliability (`lambda`) and unreliability (`theta`) matrices.
#'
#' This wrapper now delegates to `get_stage1_eb_components()`, which computes
#' EB means, posterior covariance, `lambda`, and `theta` from a single
#' consistent Stage-1 source. For iid `lme4` fits, that source is the fitted
#' iid residual covariance. For `nlme::lme` fits, fitted R-side covariance
#' matrices are used. Supplying `R_list` with an `lme4` fit is intentionally
#' rejected by the shared adapter because `lme4` estimates `beta` and `G` under
#' an iid residual likelihood.
#'
#' Lifecycle: this is a legacy convenience wrapper retained for older
#' simulation code. New code should call `get_stage1_eb_components()` directly
#' and then `select_lai_measurement_columns()` if Lai/OpenMx column names are
#' required.
#'
#' @param fit_obj Fitted Stage-1 model supported by
#' `get_stage1_eb_components()`.
#' @param split_dat Named list of cluster-level data frames. Names must include
#' every value in `ordered_ids`.
#' @param ordered_ids Character vector giving the cluster order required in the
#' output.
#' @param within_var Optional within-cluster predictor used for a random-slope
#' design. If `NULL`, an intercept-only design is used.
#' @param group Optional grouping-factor name for supported fit classes.
#' @param prefix Optional prefix for univariate output columns.
#' @param R_list Optional named list of cluster-level residual covariance
#' matrices. Names should match `ordered_ids`. If unnamed, matrices are matched
#' to `ordered_ids` by position.
#' @param .warn_legacy Logical scalar. Internal use only; controls whether the
#' legacy lifecycle warning is emitted.
#'
#' @return A tibble with one row per cluster containing EB predictions,
#' posterior variance entries, and Lai measurement-model `lambda`/`theta`
#' entries.
compute_eb_measurement_inputs <- function(fit_obj,
                                          split_dat,
                                          ordered_ids,
                                          within_var = NULL,
                                          group = NULL,
                                          prefix = "",
                                          R_list = NULL,
                                          .warn_legacy = TRUE) {
  if (isTRUE(.warn_legacy)) {
    warn_legacy_lai_input_helper("compute_eb_measurement_inputs")
  }

  if (any(!ordered_ids %in% names(split_dat))) {
    stop("`split_dat` must contain one named data frame per `ordered_ids` entry.")
  }

  split_dat <- split_dat[ordered_ids]
  data <- dplyr::bind_rows(lapply(seq_along(split_dat), function(i) {
    df_i <- split_dat[[i]]
    df_i[[".lai_cluster_id"]] <- ordered_ids[[i]]
    df_i
  }))
  response_var <- all.vars(stats::formula(fit_obj))[[1]]

  stage1_components <- get_stage1_eb_components(
    fit_obj = fit_obj,
    data = data,
    cluster_var = ".lai_cluster_id",
    outcome_var = response_var,
    within_var = within_var,
    R_list = R_list,
    group = group
  )

  n_re <- if (is.null(within_var)) 1L else 2L
  select_lai_measurement_columns(stage1_components, n_re = n_re, prefix = prefix)
}

#' Compute bivariate EB measurement inputs for random intercept/slope models.
#'
#' @param fit_obj Fitted Stage-1 model supported by
#' `get_stage1_eb_components()` with two random-effect components.
#' @param split_dat Named list of cluster-level data frames.
#' @param ordered_ids Character vector specifying output cluster order.
#' @param within_var Character scalar naming the random-slope predictor.
#' @param R_list Optional residual covariance list passed to
#' `compute_eb_measurement_inputs()`.
#' @param .warn_legacy Logical scalar. Internal use only; controls whether the
#' legacy lifecycle warning is emitted.
#'
#' @return A tibble from `compute_eb_measurement_inputs()` with unprefixed
#' bivariate EB, posterior variance, `lambda`, and `theta` columns.
compute_bivariate_eb_inputs <- function(fit_obj, split_dat, ordered_ids, within_var,
                                        R_list = NULL, .warn_legacy = TRUE) {
  if (isTRUE(.warn_legacy)) {
    warn_legacy_lai_input_helper("compute_bivariate_eb_inputs")
  }

  compute_eb_measurement_inputs(
    fit_obj = fit_obj,
    split_dat = split_dat,
    ordered_ids = ordered_ids,
    within_var = within_var,
    prefix = "",
    R_list = R_list,
    .warn_legacy = FALSE
  )
}

#' Compute univariate EB measurement inputs for random-intercept models.
#'
#' @param fit_obj Fitted Stage-1 model supported by
#' `get_stage1_eb_components()` with a single random-effect component.
#' @param split_dat Named list of cluster-level data frames.
#' @param ordered_ids Character vector specifying output cluster order.
#' @param prefix Prefix for output column names. Defaults to `"z_"` for the
#' Study 3 repeated-`z` model.
#' @param R_list Optional residual covariance list passed to
#' `compute_eb_measurement_inputs()`.
#' @param .warn_legacy Logical scalar. Internal use only; controls whether the
#' legacy lifecycle warning is emitted.
#'
#' @return A tibble from `compute_eb_measurement_inputs()` with prefixed
#' univariate EB, posterior variance, `lambda`, and `theta` columns.
compute_univariate_eb_inputs <- function(fit_obj, split_dat, ordered_ids, prefix = "z_",
                                         R_list = NULL, .warn_legacy = TRUE) {
  if (isTRUE(.warn_legacy)) {
    warn_legacy_lai_input_helper("compute_univariate_eb_inputs")
  }

  compute_eb_measurement_inputs(
    fit_obj = fit_obj,
    split_dat = split_dat,
    ordered_ids = ordered_ids,
    within_var = NULL,
    prefix = prefix,
    R_list = R_list,
    .warn_legacy = FALSE
  )
}

#' Build the stage-2 input data for the structural-slope Lai 2S-PA model.
#'
#' @details
#' This convenience wrapper is used by the random-slope sandwich coverage
#' simulation. It extracts bivariate EB measurement inputs from a first-stage
#' model with random effects `(1 + z | id)` and joins them to an `id_df`
#' containing `id` and the level-2 predictor `x`.
#'
#' @param fit_null Fitted Stage-1 model supported by
#' `get_stage1_eb_components()`.
#' @param split_dat Named list of cluster-level data frames.
#' @param id_df Level-2 data frame containing `id` and `x`.
#' @param R_list Optional residual covariance list passed to
#' `compute_bivariate_eb_inputs()`.
#' @param .warn_legacy Logical scalar. Internal use only; controls whether the
#' legacy lifecycle warning is emitted.
#'
#' @return A tibble suitable for `fit_lai_2spa()`.
compute_lai_2spa_inputs <- function(fit_null, split_dat, id_df, R_list = NULL, .warn_legacy = TRUE) {
  if (isTRUE(.warn_legacy)) {
    warn_legacy_lai_input_helper("compute_lai_2spa_inputs")
  }

  out <- compute_bivariate_eb_inputs(
    fit_obj = fit_null,
    split_dat = split_dat,
    ordered_ids = as.character(id_df$id),
    within_var = "z",
    R_list = R_list,
    .warn_legacy = FALSE
  )
  dplyr::left_join(id_df[, c("id", "x"), drop = FALSE], out, by = "id")
}

#' Run an OpenMx model with quiet retries.
#'
#' @details
#' Simulation grids occasionally produce difficult OpenMx optimization
#' problems. This wrapper calls `OpenMx::mxTryHard()` repeatedly, suppresses
#' routine optimizer chatter, optionally logs warnings, and returns the first
#' status-code-0 fit. If no successful fit is found, it returns the last fit
#' object so callers can still classify the failure.
#'
#' @param mx_mod OpenMx model object.
#' @param max_tries Number of `mxTryHard()` attempts.
#' @param warning_log Optional file path. If supplied, captured warnings are
#' appended with timestamps and attempt numbers.
#'
#' @return An OpenMx fit object, or `NULL` if all attempts fail before producing
#' a fit.
run_mx_safe <- function(mx_mod, max_tries = 5L, warning_log = NULL) {
  last_fit <- NULL
  for (attempt in seq_len(max_tries)) {
    fit_try <- tryCatch(
      local({
        fit_obj <- NULL
        warnings_seen <- character()

        # Capture console output and warnings separately: large simulation
        # batches stay quiet, while warning text remains available for logs and
        # diagnostics when requested.
        withCallingHandlers(
          suppressMessages(
            capture.output(
              fit_obj <- OpenMx::mxTryHard(
                mx_mod,
                silent = TRUE,
                verbose = 0L,
                iterationSummary = FALSE,
                bestInitsOutput = FALSE,
                showInits = FALSE,
                extraTries = 0L,
                intervals = FALSE
              ),
              file = nullfile()
            )
          ),
          warning = function(w) {
            warnings_seen <<- c(warnings_seen, conditionMessage(w))
            invokeRestart("muffleWarning")
          }
        )

        if (length(warnings_seen) > 0L && !is.null(warning_log)) {
          dir.create(dirname(warning_log), recursive = TRUE, showWarnings = FALSE)
          cat(
            sprintf("[%s] mxTryHard attempt %d warnings:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), attempt),
            paste(unique(warnings_seen), collapse = " | "),
            "\n",
            file = warning_log,
            append = TRUE
          )
        }
        fit_obj
      }),
      error = function(e) NULL
    )
    if (!is.null(fit_try)) {
      last_fit <- fit_try
      ok <- tryCatch(
        identical(fit_try$output$status$code, 0L),
        error = function(e) FALSE
      )
      if (ok) {
        return(fit_try)
      }
    }
  }
  last_fit
}

#' Extract an OpenMx standard error and related diagnostics.
#'
#' Normalize a value expected to be one numeric scalar.
#'
#' Failed OpenMx fits can leave algebras as `NULL` or `numeric(0)`. Converting
#' those values without checking their length can silently create zero-row
#' tibbles later in the result pipeline. This helper maps any non-scalar value
#' to an explicit missing value.
mx_numeric_scalar <- function(x, default = NA_real_) {
  value <- suppressWarnings(as.numeric(x))
  if (length(value) != 1L) default else value[[1L]]
}

#' Convert OpenMx audit diagnostics to an exactly one-row tibble.
#'
#' OpenMx can return `NULL`/length-zero algebra values for failed fits. A
#' length-zero list element makes `tibble::as_tibble()` produce zero rows,
#' which previously caused the entire 2S-PA method row to disappear. This
#' normalizer preserves logical boundary fields and converts every other
#' malformed scalar to `NA_real_` before constructing the tibble.
mx_diagnostics_tibble <- function(diagnostics) {
  for (name in names(diagnostics)) {
    value <- diagnostics[[name]]
    if (grepl("_boundary$", name)) {
      diagnostics[[name]] <- if (length(value) == 1L) {
        as.logical(value[[1L]])
      } else {
        NA
      }
    } else {
      diagnostics[[name]] <- mx_numeric_scalar(value)
    }
  }
  out <- tibble::as_tibble(diagnostics)
  if (nrow(out) != 1L) {
    stop("Internal error: OpenMx diagnostics must contain exactly one row.")
  }
  out
}

#' Safely extract an OpenMx algebra standard error.
#'
#' @details
#' `OpenMx::mxSE(..., details = TRUE)` can warn or error when the information
#' matrix is indefinite or the repeated-sampling covariance cannot be computed.
#' In addition, `mxSE()` calls `vcov.MxModel()` internally; older OpenMx
#' versions print an uncaught error when a failed fit has no scalar `fitUnits`.
#' The preflight checks here avoid that noisy call and normalize the failure
#' into fields that `extract_mx_stats()` can classify without interrupting or
#' dropping a simulation row.
#'
#' @param algebra_name Character scalar naming an OpenMx algebra to pass to
#' `mxSE()`.
#' @param mx_fit Fitted OpenMx model.
#'
#' @return A list with `se`, `cov`, `warning`, and `error` entries.
extract_mx_se_details <- function(algebra_name, mx_fit) {
  fit_units <- tryCatch(mx_fit$output$fitUnits, error = function(e) NULL)
  parameter_vcov <- tryCatch(mx_fit$output[["vcov"]], error = function(e) NULL)
  if (length(fit_units) != 1L || is.na(fit_units) ||
      !(fit_units %in% c("-2lnL", "r'Wr", "r'wr"))) {
    return(list(
      se = NA_real_,
      cov = matrix(NA_real_, nrow = 1L, ncol = 1L),
      warning = NA_character_,
      error = "OpenMx fit units were unavailable; algebra SE was not requested."
    ))
  }
  if (is.null(parameter_vcov) || length(parameter_vcov) == 0L) {
    return(list(
      se = NA_real_,
      cov = matrix(NA_real_, nrow = 1L, ncol = 1L),
      warning = NA_character_,
      error = "OpenMx parameter covariance was unavailable; algebra SE was not requested."
    ))
  }
  se_warnings <- character()

  se_details <- withCallingHandlers(
    tryCatch(
      OpenMx::mxSE(algebra_name, mx_fit, forceName = TRUE, silent = TRUE, details = TRUE),
      error = function(e) e
    ),
    warning = function(w) {
      se_warnings <<- c(se_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  if (inherits(se_details, "error")) {
    return(list(
      se = NA_real_,
      cov = matrix(NA_real_, nrow = 1L, ncol = 1L),
      warning = compact_message(se_warnings),
      error = compact_message(conditionMessage(se_details))
    ))
  }

  list(
    se = tryCatch(mx_numeric_scalar(se_details$SE), error = function(e) NA_real_),
    cov = tryCatch(as.matrix(se_details$Cov), error = function(e) matrix(NA_real_, nrow = 1L, ncol = 1L)),
    warning = compact_message(se_warnings),
    error = NA_character_
  )
}

#' Classify OpenMx fit and standard-error problems.
#'
#' @details
#' The simulations need compact categorical failure labels rather than raw
#' OpenMx status messages. This function prioritizes hard fit failures, then
#' non-finite estimates, then standard-error/covariance problems, and finally
#' information-matrix warnings.
#'
#' @param code OpenMx status code.
#' @param status_msg OpenMx status message.
#' @param est Numeric algebra estimate.
#' @param se Numeric standard error.
#' @param info_definite Logical OpenMx information-definiteness flag, or `NA`.
#' @param se_warning Warning text captured from `mxSE()`.
#' @param se_error Error text captured from `mxSE()`.
#' @param se_cov Repeated-sampling covariance matrix returned by `mxSE()`.
#'
#' @return A character scalar issue class such as `"ok"`, `"status_10"`,
#' `"se_vcov_unavailable"`, or `"information_matrix_not_definite"`.
classify_mx_issue <- function(code, status_msg, est, se, info_definite, se_warning, se_error, se_cov) {
  code_int <- if (length(code) == 0L || is.null(code) || is.na(code)) NA_integer_ else as.integer(code)
  cov_diag <- tryCatch(diag(as.matrix(se_cov)), error = function(e) numeric())
  vcov_bad <- length(cov_diag) == 0L || any(!is.finite(cov_diag)) || any(cov_diag < 0)
  se_problem_text <- compact_message(c(se_warning, se_error))

  if (is.na(code_int)) {
    return("mx_fit_null")
  }
  if (identical(code_int, 10L)) {
    if (!is.na(status_msg) && grepl("not positive definite", status_msg, ignore.case = TRUE)) {
      return("non_pd_implied_cov")
    }
    return("status_10")
  }
  if (identical(code_int, 6L)) {
    return("indefinite_hessian_status6")
  }
  if (!identical(code_int, 0L)) {
    return(sprintf("status_%d", code_int))
  }
  if (!is.finite(est)) {
    return("estimate_not_finite")
  }
  if (!is.na(info_definite) && identical(info_definite, FALSE) && !is.finite(se)) {
    return("information_matrix_not_definite")
  }
  if (!is.finite(se)) {
    if ((is.na(se_problem_text) || nzchar(se_problem_text)) &&
      (grepl("vcov", se_problem_text, ignore.case = TRUE) ||
        grepl("standard errors", se_problem_text, ignore.case = TRUE) ||
        grepl("repeated-sampling covariance", se_problem_text, ignore.case = TRUE) ||
        vcov_bad)) {
      return("se_vcov_unavailable")
    }
    return("se_not_finite")
  }
  if (!is.na(info_definite) && identical(info_definite, FALSE)) {
    return("information_matrix_not_definite")
  }

  "ok"
}

#' Extract a standardized estimator row from an OpenMx fit.
#'
#' @details
#' This is the common post-processing path for the Lai OpenMx wrappers. It
#' reads model status, evaluates the requested algebra, obtains the algebra
#' standard error, classifies any problem, and returns a stable one-row tibble
#' with estimate, standard error, confidence limits, and diagnostic metadata.
#'
#' Estimates and SEs are set to `NA` unless the issue class is `"ok"`. This
#' keeps downstream coverage summaries from treating questionable OpenMx
#' solutions as successful estimates while preserving enough metadata to audit
#' the failure mode.
#'
#' @param mx_fit Fitted OpenMx model, or `NULL`.
#' @param algebra_name Character scalar naming the algebra whose estimate and
#' SE should be extracted.
#' @param ci_multiplier Multiplier used for Wald confidence limits. Defaults to
#' normal 1.96-style limits; callers can pass a `t` multiplier when matching a
#' specific simulation convention.
#' @param raw_algebra_name Optional algebra naming the focal structural path
#'   before any reporting-scale multiplication. Its estimate and SE are saved
#'   as `mx_raw_focal_estimate` and `mx_raw_focal_se`.
#' @param latent_covariance_blocks Named list of two-element character vectors.
#'   Each vector gives the fitted latent intercept and slope names whose RAM
#'   `S` block is audited. The default records `u0`/`u1`; dual-process callers
#'   request separate predictor and outcome blocks.
#'
#' Extract fitted latent covariance and unscaled focal-path diagnostics.
mx_latent_covariance_diagnostics <- function(
    mx_fit, raw_algebra_name = NULL,
    latent_covariance_blocks = list(latent = c("u0", "u1")),
    boundary_tolerance = sqrt(.Machine$double.eps)) {
  diagnostics <- list(
    mx_raw_focal_estimate = NA_real_,
    mx_raw_focal_se = NA_real_
  )
  if (!is.null(mx_fit) && !is.null(raw_algebra_name)) {
    diagnostics$mx_raw_focal_estimate <- tryCatch(
      mx_numeric_scalar(OpenMx::mxEvalByName(raw_algebra_name, mx_fit)),
      error = function(e) NA_real_
    )
    diagnostics$mx_raw_focal_se <- tryCatch(
      mx_numeric_scalar(extract_mx_se_details(raw_algebra_name, mx_fit)$se),
      error = function(e) NA_real_
    )
  }

  block_prefix <- function(block_name) {
    if (identical(block_name, "latent") || !nzchar(block_name)) {
      "mx_latent_"
    } else {
      paste0("mx_", block_name, "_latent_")
    }
  }
  for (block_name in names(latent_covariance_blocks)) {
    prefix <- block_prefix(block_name)
    for (metric in c(
      "intercept_variance", "slope_variance",
      "intercept_slope_covariance", "intercept_slope_correlation",
      "covariance_min_eigenvalue", "covariance_condition_number"
    )) {
      diagnostics[[paste0(prefix, metric)]] <- NA_real_
    }
    diagnostics[[paste0(prefix, "covariance_boundary")]] <- NA
  }

  if (is.null(mx_fit) || length(latent_covariance_blocks) == 0L) {
    return(mx_diagnostics_tibble(diagnostics))
  }
  fitted_s <- tryCatch(OpenMx::mxEval(S, mx_fit), error = function(e) NULL)
  if (is.null(fitted_s) || !is.matrix(fitted_s) ||
      any(dim(fitted_s) == 0L)) {
    return(mx_diagnostics_tibble(diagnostics))
  }
  variable_names <- c(mx_fit$manifestVars, mx_fit$latentVars)
  for (block_name in names(latent_covariance_blocks)) {
    latent_names <- latent_covariance_blocks[[block_name]]
    indices <- match(latent_names, variable_names)
    if (length(indices) != 2L || anyNA(indices)) next
    covariance <- as.matrix(fitted_s[indices, indices, drop = FALSE])
    covariance <- (covariance + t(covariance)) / 2
    eigenvalues <- tryCatch(
      eigen(covariance, symmetric = TRUE, only.values = TRUE)$values,
      error = function(e) rep(NA_real_, 2L)
    )
    prefix <- block_prefix(block_name)
    diagnostics[[paste0(prefix, "intercept_variance")]] <- covariance[1L, 1L]
    diagnostics[[paste0(prefix, "slope_variance")]] <- covariance[2L, 2L]
    diagnostics[[paste0(prefix, "intercept_slope_covariance")]] <- covariance[1L, 2L]
    diagnostics[[paste0(prefix, "intercept_slope_correlation")]] <- if (
      all(is.finite(covariance)) && all(diag(covariance) > 0)
    ) {
      covariance[1L, 2L] / sqrt(covariance[1L, 1L] * covariance[2L, 2L])
    } else {
      NA_real_
    }
    diagnostics[[paste0(prefix, "covariance_min_eigenvalue")]] <- if (
      all(is.finite(eigenvalues))
    ) min(eigenvalues) else NA_real_
    diagnostics[[paste0(prefix, "covariance_condition_number")]] <- tryCatch(
      kappa(covariance, exact = TRUE),
      error = function(e) NA_real_
    )
    diagnostics[[paste0(prefix, "covariance_boundary")]] <-
      any(!is.finite(eigenvalues)) || any(eigenvalues <= boundary_tolerance) ||
      any(diag(covariance) <= boundary_tolerance)
  }
  mx_diagnostics_tibble(diagnostics)
}

#' @return A one-row tibble with estimator columns (`estimate`, `se`, `ci_low`,
#' `ci_high`, `status_code`) plus OpenMx diagnostic columns.
extract_mx_stats <- function(
    mx_fit, algebra_name = "xstd_u1", ci_multiplier = stats::qnorm(0.975),
    raw_algebra_name = NULL,
    latent_covariance_blocks = list(latent = c("u0", "u1"))) {
  extended_diagnostics <- mx_latent_covariance_diagnostics(
    mx_fit,
    raw_algebra_name = raw_algebra_name,
    latent_covariance_blocks = latent_covariance_blocks
  )
  if (is.null(mx_fit)) {
    return(dplyr::bind_cols(tibble::tibble(
      estimate = NA_real_, se = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
      status_code = NA_integer_, mx_status_msg = NA_character_, mx_info_definite = NA,
      mx_condition_number = NA_real_, mx_issue_class = "mx_fit_null",
      mx_issue_detail = "OpenMx fit object was NULL."
    ), extended_diagnostics))
  }

  code_value <- suppressWarnings(as.integer(mx_fit$output$status$code))
  code <- if (length(code_value) == 1L) code_value[[1L]] else NA_integer_
  status_msg <- compact_message(mx_fit$output$status$statusMsg)
  info_definite <- if (is.null(mx_fit$output$infoDefinite) || length(mx_fit$output$infoDefinite) == 0L) {
    NA
  } else {
    as.logical(mx_fit$output$infoDefinite[[1L]])
  }
  condition_number <- if (is.null(mx_fit$output$conditionNumber) || length(mx_fit$output$conditionNumber) == 0L) {
    NA_real_
  } else {
    mx_numeric_scalar(mx_fit$output$conditionNumber)
  }

  # Nonzero status codes are treated as estimator failures before attempting
  # algebra extraction or SE computation; both can be misleading after failed
  # optimization.
  if (!identical(code, 0L)) {
    issue_class <- classify_mx_issue(
      code = code,
      status_msg = status_msg,
      est = NA_real_,
      se = NA_real_,
      info_definite = info_definite,
      se_warning = NA_character_,
      se_error = NA_character_,
      se_cov = matrix(NA_real_, nrow = 1L, ncol = 1L)
    )
    issue_detail <- compact_message(c(
      status_msg,
      if (identical(issue_class, "indefinite_hessian_status6")) "OpenMx reported an uncertain solution / indefinite Hessian." else NA_character_
    ))

    return(dplyr::bind_cols(tibble::tibble(
      estimate = NA_real_, se = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
      status_code = code, mx_status_msg = status_msg, mx_info_definite = info_definite,
      mx_condition_number = condition_number, mx_issue_class = issue_class,
      mx_issue_detail = issue_detail
    ), extended_diagnostics))
  }

  est_error <- NA_character_
  est <- tryCatch(
    mx_numeric_scalar(OpenMx::mxEvalByName(algebra_name, mx_fit)),
    error = function(e) {
      est_error <<- compact_message(conditionMessage(e))
      NA_real_
    }
  )
  se_details <- extract_mx_se_details(algebra_name, mx_fit)
  se <- se_details$se

  # Classification uses both OpenMx status fields and mxSE diagnostics so
  # summaries can distinguish convergence failures from information-matrix or
  # covariance-estimation problems.
  issue_class <- classify_mx_issue(
    code = code,
    status_msg = status_msg,
    est = est,
    se = se,
    info_definite = info_definite,
    se_warning = se_details$warning,
    se_error = se_details$error,
    se_cov = se_details$cov
  )
  issue_detail <- compact_message(c(status_msg, est_error, se_details$warning, se_details$error))
  est_out <- if (identical(issue_class, "ok")) est else NA_real_
  se_out <- if (identical(issue_class, "ok")) se else NA_real_

  dplyr::bind_cols(tibble::tibble(
    estimate = est_out,
    se = se_out,
    ci_low = if (is.finite(est_out) && is.finite(se_out)) est_out - ci_multiplier * se_out else NA_real_,
    ci_high = if (is.finite(est_out) && is.finite(se_out)) est_out + ci_multiplier * se_out else NA_real_,
    status_code = code,
    mx_status_msg = status_msg,
    mx_info_definite = info_definite,
    mx_condition_number = condition_number,
    mx_issue_class = issue_class,
    mx_issue_detail = issue_detail
  ), extended_diagnostics)
}

#' Fit a Lai 2S-PA model for `x -> latent random slope`.
#'
#' @details
#' This wrapper is used by the sandwich coverage simulation. The manifest
#' variables are the level-2 predictor `x` and the EB measures `u0_eb` and
#' `u1_eb`; the latent variables are the true random intercept and random
#' slope. Fixed per-row loadings (`lambda`) and residual covariance entries
#' (`theta`) come from `compute_lai_2spa_inputs()`.
#'
#' The reported algebra `gamma_hat = A[5, 1]` is the OpenMx RAM path from `x`
#' to the latent random slope `u1`. Confidence limits use a `t` multiplier with
#' `nrow(stage2_df) - 2` degrees of freedom to match the legacy sandwich
#' simulation convention.
#'
#' @param stage2_df Data frame containing `x`, `u0_eb`, `u1_eb`, `lambda11`,
#' `lambda12`, `lambda21`, `lambda22`, `theta11`, `theta12`, and `theta22`.
#' @param use_average Logical. If `TRUE`, use average `lambda`/`theta` entries
#' as fixed constants; if `FALSE`, use definition variables from each row.
#'
#' @return A one-row tibble from `extract_mx_stats()` for `gamma_hat`.
fit_lai_2spa <- function(stage2_df, use_average = FALSE) {
  loading_cols <- openmx_bivariate_loading_columns()
  if (isTRUE(use_average)) {
    loading_arg <- list(values = colMeans(stage2_df[, loading_cols], na.rm = TRUE))
    theta_arg <- list(values = colMeans(stage2_df[, c("theta11", "theta12", "theta22")], na.rm = TRUE))
  } else {
    loading_arg <- list(labels = paste0("data.", loading_cols))
    theta_arg <- list(labels = paste0("data.", c("theta11", "theta12", "theta22")))
  }

  mx_mod <- OpenMx::mxModel(
    if (isTRUE(use_average)) "lai_2spaa_structural_slope" else "lai_2spa_structural_slope",
    type = "RAM",
    manifestVars = c("x", "u0_eb", "u1_eb"),
    latentVars = c("u0", "u1"),
    OpenMx::mxData(stage2_df, type = "raw"),
    OpenMx::mxPath(from = "x", arrows = 2, free = TRUE, values = stats::var(stage2_df$x), labels = "var_x"),
    OpenMx::mxPath(from = "x", to = "u1", free = TRUE, values = 0, labels = "gamma_x_u1"),

    # EB measures are fallible indicators of the latent random effects. Their
    # loadings are fixed to the per-cluster reliability matrix from the
    # first-stage mixed model.
    do.call(OpenMx::mxPath, c(list(
      from = c("u0", "u1"),
      to = c("u0_eb", "u1_eb"),
      connect = "unique.bivariate",
      free = FALSE
    ), loading_arg)),

    # The EB residual covariance is fixed to theta, the unreliability matrix.
    do.call(OpenMx::mxPath, c(list(
      from = c("u0_eb", "u1_eb"),
      connect = "unique.pairs",
      arrows = 2,
      free = FALSE
    ), theta_arg)),
    OpenMx::mxPath(
      from = c("u0", "u1"),
      connect = "unique.pairs",
      arrows = 2,
      free = TRUE,
      values = c(0.5, 0.1, 0.5),
      labels = c("var_u0", "cov_u0_u1", "var_u1")
    ),
    OpenMx::mxPath(from = "one", to = c("x", "u0", "u1"), free = TRUE, values = c(mean(stage2_df$x), 0, 0)),
    OpenMx::mxAlgebra(A[5, 1], name = "raw_focal_path"),
    OpenMx::mxAlgebra(A[5, 1], name = "gamma_hat")
  )

  extract_mx_stats(
    run_mx_safe(mx_mod),
    algebra_name = "gamma_hat",
    ci_multiplier = stats::qt(0.975, nrow(stage2_df) - 2L),
    raw_algebra_name = "raw_focal_path",
    latent_covariance_blocks = list(latent = c("u0", "u1"))
  )
}

#' Fit Lai 2S-PA/2S-PAA for matched-clustering observed outcomes.
#'
#' @details
#' This is the OpenMx estimator used by Lai replication Studies 1 and 2, where
#' the first-stage random effects and observed outcome `z` share the same
#' clusters. It models EB random-effect measures as fallible indicators of
#' latent `u0` and `u1`, then estimates the effect of latent random slope `u1`
#' on observed `z`.
#'
#' The reported algebra `xstd_u1 = A[3, 5] * sqrt(S[5, 5])` scales the path by
#' the latent slope SD, matching the one-SD estimand used by the stage-2
#' regression helpers.
#'
#' @param stage2_df Data frame with EB measurement inputs and observed outcome
#' `z`.
#' @param use_average Logical. If `TRUE`, fit the 2S-PAA variant using average
#' `lambda`/`theta`; otherwise fit the row-specific 2S-PA variant.
#' @param u0_start Optional starting value for the nuisance latent-intercept
#'   path to `z`. When omitted, the helper uses a compatible value from
#'   `fixed_params` when available.
#' @param reporting_scale Optional fixed scale for the reported slope effect.
#'   When omitted, the fitted latent slope SD is used.
#'
#' @return A one-row tibble from `extract_mx_stats()` for `xstd_u1`.
fit_lai_2spa_observed_outcome <- function(
    stage2_df,
    use_average = FALSE,
    u0_start = NULL,
    reporting_scale = NULL) {
  if (is.null(u0_start)) {
    u0_start <- if (exists("fixed_params", inherits = TRUE) &&
        !is.null(fixed_params$beta_zu0)) {
      fixed_params$beta_zu0
    } else if (exists("fixed_params", inherits = TRUE) &&
        !is.null(fixed_params$beta1z)) {
      fixed_params$beta1z
    } else {
      0
    }
  }
  loading_cols <- openmx_bivariate_loading_columns()
  if (isTRUE(use_average)) {
    loading_arg <- list(values = colMeans(stage2_df[, loading_cols], na.rm = TRUE))
    theta_arg <- list(values = colMeans(stage2_df[, c("theta11", "theta12", "theta22")], na.rm = TRUE))
  } else {
    loading_arg <- list(labels = paste0("data.", loading_cols))
    theta_arg <- list(labels = paste0("data.", c("theta11", "theta12", "theta22")))
  }
  target_algebra <- if (is.null(reporting_scale)) {
    OpenMx::mxAlgebra(A[3, 5] * sqrt(S[5, 5]), name = "xstd_u1")
  } else {
    reporting_scale <- as.numeric(reporting_scale[[1]])
    if (!is.finite(reporting_scale) || reporting_scale <= 0) {
      stop("`reporting_scale` must be finite and positive.")
    }
    eval(substitute(
      OpenMx::mxAlgebra(A[3, 5] * SCALE, name = "xstd_u1"),
      list(SCALE = reporting_scale)
    ))
  }

  mx_mod <- OpenMx::mxModel(
    if (isTRUE(use_average)) "lai_2spaa" else "lai_2spa",
    type = "RAM",
    manifestVars = c("u0_eb", "u1_eb", "z"),
    latentVars = c("u0", "u1"),
    OpenMx::mxData(stage2_df[, c("u0_eb", "u1_eb", "z", "lambda11", "lambda12", "lambda21", "lambda22", "theta11", "theta12", "theta22")], type = "raw"),

    # The u0 path is a nuisance coefficient; the u1 path is the target.
    OpenMx::mxPath(
      from = c("u0", "u1"),
      to = "z",
      free = TRUE,
      values = c(as.numeric(u0_start[[1]]), 0)
    ),
    OpenMx::mxPath(from = "z", arrows = 2, free = TRUE, values = max(1e-3, stats::var(stage2_df$z, na.rm = TRUE))),
    do.call(OpenMx::mxPath, c(list(
      from = c("u0", "u1"), to = c("u0_eb", "u1_eb"), connect = "unique.bivariate", free = FALSE
    ), loading_arg)),
    do.call(OpenMx::mxPath, c(list(
      from = c("u0_eb", "u1_eb"), connect = "unique.pairs", arrows = 2, free = FALSE
    ), theta_arg)),
    OpenMx::mxPath(from = c("u0", "u1"), connect = "unique.pairs", arrows = 2, free = TRUE, values = c(0.5, 0.1, 0.5)),
    OpenMx::mxPath(from = "one", to = c("u0", "u1", "z"), free = TRUE, values = c(0, 0, mean(stage2_df$z, na.rm = TRUE))),
    OpenMx::mxAlgebra(A[3, 5], name = "raw_focal_path"),
    target_algebra
  )

  mx_fit <- run_mx_safe(mx_mod)
  out <- extract_mx_stats(
    mx_fit,
    raw_algebra_name = "raw_focal_path",
    latent_covariance_blocks = list(latent = c("u0", "u1"))
  )

  # Preserve the model-estimated latent slope SD even when the caller requests
  # a fixed reporting scale.  Lai Study 1 multiplied the 2S-PA path by this
  # fitted SD; refreshed simulations also need it to reconstruct that
  # historical reporting convention without refitting the model.
  latent_slope_sd <- tryCatch({
    if (!identical(out$status_code[[1]], 0L)) {
      NA_real_
    } else {
      latent_covariance <- OpenMx::mxEval(S, mx_fit)
      slope_variance <- as.numeric(latent_covariance[5L, 5L])
      if (is.finite(slope_variance) && slope_variance > 0) sqrt(slope_variance) else NA_real_
    }
  }, error = function(e) NA_real_)
  dplyr::mutate(out, lai_fitted_latent_slope_sd = latent_slope_sd)
}

#' Fit Lai 2S-PA/2S-PAA for disparate-clustering Study 3.
#'
#' @details
#' Study 3 has separate first-stage models for the random effects underlying
#' `y` and the repeated observed `z` outcome. This wrapper models `u0_eb` and
#' `u1_eb` as fallible indicators of latent `u0`/`u1`, and `z_u0_eb` as a
#' fallible indicator of `z_lat`. The structural target is the effect of latent
#' random slope `u1` on `z_lat`.
#'
#' The reported algebra `xstd_u1 = A[6, 5] * sqrt(S[5, 5])` matches the
#' standardized latent-slope estimand used elsewhere in the replication.
#'
#' @param stage2_df Data frame containing bivariate EB measurement inputs,
#' prefixed univariate `z_` measurement inputs, and `z_u0_eb`.
#' @param use_average Logical. If `TRUE`, fit the 2S-PAA variant using average
#' loadings/residual variances; otherwise fit the row-specific 2S-PA variant.
#'
#' @return A one-row tibble from `extract_mx_stats()` for `xstd_u1`.
fit_lai_2spa_disparate <- function(stage2_df, use_average = FALSE) {
  loading_cols <- openmx_bivariate_loading_columns()
  if (isTRUE(use_average)) {
    loading_arg <- list(values = colMeans(stage2_df[, loading_cols], na.rm = TRUE))
    theta_arg <- list(values = colMeans(stage2_df[, c("theta11", "theta12", "theta22")], na.rm = TRUE))
    z_loading_arg <- list(values = mean(stage2_df$z_lambda11, na.rm = TRUE))
    z_theta_arg <- list(values = mean(stage2_df$z_theta11, na.rm = TRUE))
  } else {
    loading_arg <- list(labels = paste0("data.", loading_cols))
    theta_arg <- list(labels = paste0("data.", c("theta11", "theta12", "theta22")))
    z_loading_arg <- list(labels = "data.z_lambda11")
    z_theta_arg <- list(labels = "data.z_theta11")
  }

  mx_mod <- OpenMx::mxModel(
    if (isTRUE(use_average)) "lai_2spaa_disparate" else "lai_2spa_disparate",
    type = "RAM",
    manifestVars = c("u0_eb", "u1_eb", "z_u0_eb"),
    latentVars = c("u0", "u1", "z_lat"),
    OpenMx::mxData(stage2_df, type = "raw"),

    # The observed z process is represented by latent z_lat because its EB
    # estimate comes from a separate random-intercept measurement model.
    OpenMx::mxPath(from = c("u0", "u1"), to = "z_lat", free = TRUE, values = c(fixed_params$beta_zu0, 0)),
    OpenMx::mxPath(from = "z_lat", arrows = 2, free = TRUE, values = 1),
    do.call(OpenMx::mxPath, c(list(
      from = c("u0", "u1"), to = c("u0_eb", "u1_eb"), connect = "unique.bivariate", free = FALSE
    ), loading_arg)),
    do.call(OpenMx::mxPath, c(list(
      from = c("u0_eb", "u1_eb"), connect = "unique.pairs", arrows = 2, free = FALSE
    ), theta_arg)),
    do.call(OpenMx::mxPath, c(list(from = "z_lat", to = "z_u0_eb", free = FALSE), z_loading_arg)),
    do.call(OpenMx::mxPath, c(list(from = "z_u0_eb", arrows = 2, free = FALSE), z_theta_arg)),
    OpenMx::mxPath(from = c("u0", "u1"), connect = "unique.pairs", arrows = 2, free = TRUE, values = c(0.5, 0.1, 0.5)),
    OpenMx::mxPath(from = "one", to = c("u0", "u1", "z_lat"), free = TRUE, values = c(0, 0, fixed_params$z_intercept)),
    OpenMx::mxAlgebra(A[6, 5], name = "raw_focal_path"),
    OpenMx::mxAlgebra(A[6, 5] * sqrt(S[5, 5]), name = "xstd_u1")
  )

  extract_mx_stats(
    run_mx_safe(mx_mod),
    raw_algebra_name = "raw_focal_path",
    latent_covariance_blocks = list(latent = c("u0", "u1"))
  )
}

#' Fit Lai 2S-PA/2S-PAA for two random-intercept/random-slope processes.
#'
#' The predictor process supplies latent `y_u0` and `y_u1`; the outcome process
#' supplies latent `q_u0` and `q_u1`. The focal path is `y_u1 -> q_u1`, adjusted
#' for `y_u0`. Both bivariate EB measurement models retain their full loading and
#' residual covariance matrices. When `reporting_scale` is supplied, the raw
#' focal path and its standard error are multiplied by that fixed scale.
fit_lai_2spa_dual_process <- function(
    stage2_df,
    use_average = FALSE,
    theta0_start = 0,
    theta1_start = 0,
    reporting_scale = NULL) {
  y_loading_cols <- openmx_bivariate_loading_columns()
  y_theta_cols <- c("theta11", "theta12", "theta22")
  q_loading_cols <- openmx_bivariate_loading_columns("q_")
  q_theta_cols <- paste0("q_", y_theta_cols)

  if (isTRUE(use_average)) {
    y_loading_arg <- list(values = colMeans(stage2_df[, y_loading_cols], na.rm = TRUE))
    y_theta_arg <- list(values = colMeans(stage2_df[, y_theta_cols], na.rm = TRUE))
    q_loading_arg <- list(values = colMeans(stage2_df[, q_loading_cols], na.rm = TRUE))
    q_theta_arg <- list(values = colMeans(stage2_df[, q_theta_cols], na.rm = TRUE))
  } else {
    y_loading_arg <- list(labels = paste0("data.", y_loading_cols))
    y_theta_arg <- list(labels = paste0("data.", y_theta_cols))
    q_loading_arg <- list(labels = paste0("data.", q_loading_cols))
    q_theta_arg <- list(labels = paste0("data.", q_theta_cols))
  }

  target_algebra <- if (is.null(reporting_scale)) {
    OpenMx::mxAlgebra(
      A[8, 6] * sqrt(S[6, 6] / q1_variance),
      name = "xstd_u1"
    )
  } else {
    reporting_scale <- as.numeric(reporting_scale[[1]])
    if (!is.finite(reporting_scale) || reporting_scale <= 0) {
      stop("`reporting_scale` must be finite and positive.")
    }
    eval(substitute(
      OpenMx::mxAlgebra(A[8, 6] * SCALE, name = "xstd_u1"),
      list(SCALE = reporting_scale)
    ))
  }

  mx_mod <- OpenMx::mxModel(
    if (isTRUE(use_average)) "lai_2spaa_dual_process" else "lai_2spa_dual_process",
    type = "RAM",
    manifestVars = c("u0_eb", "u1_eb", "q_u0_eb", "q_u1_eb"),
    latentVars = c("y_u0", "y_u1", "q_u0", "q_u1"),
    OpenMx::mxData(stage2_df, type = "raw"),
    do.call(OpenMx::mxPath, c(list(
      from = c("y_u0", "y_u1"),
      to = c("u0_eb", "u1_eb"),
      connect = "unique.bivariate",
      free = FALSE
    ), y_loading_arg)),
    do.call(OpenMx::mxPath, c(list(
      from = c("u0_eb", "u1_eb"),
      connect = "unique.pairs",
      arrows = 2,
      free = FALSE
    ), y_theta_arg)),
    do.call(OpenMx::mxPath, c(list(
      from = c("q_u0", "q_u1"),
      to = c("q_u0_eb", "q_u1_eb"),
      connect = "unique.bivariate",
      free = FALSE
    ), q_loading_arg)),
    do.call(OpenMx::mxPath, c(list(
      from = c("q_u0_eb", "q_u1_eb"),
      connect = "unique.pairs",
      arrows = 2,
      free = FALSE
    ), q_theta_arg)),
    OpenMx::mxPath(
      from = c("y_u0", "y_u1"),
      to = "q_u1",
      free = TRUE,
      values = c(theta0_start, theta1_start)
    ),
    OpenMx::mxPath(
      from = c("y_u0", "y_u1"),
      connect = "unique.pairs",
      arrows = 2,
      free = TRUE,
      values = c(0.5, 0.1, 0.5)
    ),
    OpenMx::mxPath(
      from = c("q_u0", "q_u1"),
      connect = "unique.pairs",
      arrows = 2,
      free = TRUE,
      values = c(0.5, 0.1, 0.5)
    ),
    OpenMx::mxPath(
      from = "one",
      to = c("y_u0", "y_u1", "q_u0", "q_u1"),
      free = TRUE,
      values = 0
    ),
    OpenMx::mxAlgebra(
      S[8, 8] +
        A[8, 5]^2 * S[5, 5] +
        A[8, 6]^2 * S[6, 6] +
        2 * A[8, 5] * A[8, 6] * S[5, 6],
      name = "q1_variance"
    ),
    OpenMx::mxAlgebra(A[8, 6], name = "raw_focal_path"),
    target_algebra
  )

  extract_mx_stats(
    run_mx_safe(mx_mod),
    raw_algebra_name = "raw_focal_path",
    latent_covariance_blocks = list(
      predictor = c("y_u0", "y_u1"),
      outcome = c("q_u0", "q_u1")
    )
  )
}
