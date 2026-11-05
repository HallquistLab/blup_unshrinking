mplus_diagnostics_template <- function() {
  tibble::tibble(
    mplus_warning_count = NA_integer_,
    mplus_critical_warning = NA,
    mplus_critical_warning_detail = NA_character_,
    mplus_target_parameter_count = NA_integer_,
    mplus_raw_focal_estimate = NA_real_,
    mplus_raw_focal_se = NA_real_,
    mplus_fitted_latent_slope_variance = NA_real_,
    mplus_fitted_latent_slope_sd = NA_real_,
    mplus_min_reported_variance = NA_real_,
    mplus_boundary_variance_count = NA_integer_,
    mplus_nonpositive_variance_count = NA_integer_,
    mplus_latent_intercept_variance = NA_real_,
    mplus_latent_slope_variance = NA_real_,
    mplus_latent_intercept_slope_covariance = NA_real_,
    mplus_latent_intercept_slope_correlation = NA_real_,
    mplus_latent_covariance_min_eigenvalue = NA_real_,
    mplus_latent_covariance_boundary = NA,
    mplus_predictor_latent_intercept_variance = NA_real_,
    mplus_predictor_latent_slope_variance = NA_real_,
    mplus_predictor_latent_intercept_slope_covariance = NA_real_,
    mplus_predictor_latent_intercept_slope_correlation = NA_real_,
    mplus_predictor_latent_covariance_min_eigenvalue = NA_real_,
    mplus_predictor_latent_covariance_boundary = NA,
    mplus_outcome_latent_intercept_variance = NA_real_,
    mplus_outcome_latent_slope_variance = NA_real_,
    mplus_outcome_latent_intercept_slope_covariance = NA_real_,
    mplus_outcome_latent_intercept_slope_correlation = NA_real_,
    mplus_outcome_latent_covariance_min_eigenvalue = NA_real_,
    mplus_outcome_latent_covariance_boundary = NA
  )
}

mplus_message_lines <- function(messages) {
  if (is.null(messages) || length(messages) == 0L) {
    return(character())
  }
  lines <- unlist(lapply(messages, paste, collapse = " "), use.names = FALSE)
  lines <- trimws(as.character(lines))
  unique(lines[nzchar(lines) & !is.na(lines)])
}

#' Run one Mplus model with a writable Fortran scratch directory.
#'
#' @details
#' Mplus creates a temporary `~model.tst` file before fitting. On managed or
#' sandboxed systems, the inherited `TMPDIR` can point to a readable but
#' non-writable application scratch path. Mplus then exits with internal error
#' `FE1001 9`, and its short error output may contain NUL padding that prevents
#' `MplusAutomation::readModels()` from parsing the actual cause.
#'
#' The directory containing `model_file` is already writable because
#' `mplusModeler()` writes the input, data, and output files there. Following
#' the safe approach used by `submitModels()`, this helper creates a unique
#' private scratch directory beneath it for each Mplus invocation. This matters
#' because Mplus uses the fixed filename `~model.tst`; merely pointing several
#' concurrent models at one writable directory would still permit collisions.
#' The scratch directory is removed and the previous environment value is
#' restored even if fitting errors.
#'
#' @param model An `mplusObject`.
#' @param model_file Path to the generated `.inp` file.
#'
#' @param ... Additional arguments passed to `MplusAutomation::mplusModeler()`.
#'
#' @return The result from `MplusAutomation::mplusModeler()`. If model fitting
#'   throws an R error, a `mplus_modeler_failure` object retaining the error and
#'   warning messages is returned.
run_mplus_modeler_writable_tmp <- function(model, model_file, ...) {
  model_dir <- normalizePath(dirname(model_file), mustWork = TRUE)
  scratch_dir <- tempfile(
    pattern = paste0("mplus-scratch-", Sys.getpid(), "-"),
    tmpdir = model_dir
  )
  if (!dir.create(scratch_dir, mode = "0700", showWarnings = FALSE)) {
    stop("Could not create a private Mplus scratch directory in: ", model_dir)
  }
  old_tmpdir <- Sys.getenv("TMPDIR", unset = NA_character_)
  on.exit({
    if (is.na(old_tmpdir)) {
      Sys.unsetenv("TMPDIR")
    } else {
      Sys.setenv(TMPDIR = old_tmpdir)
    }
    unlink(scratch_dir, recursive = TRUE, force = TRUE)
  }, add = TRUE)
  Sys.setenv(TMPDIR = scratch_dir)

  warning_messages <- character()
  tryCatch(
    withCallingHandlers(
      MplusAutomation::mplusModeler(
        model,
        run = 1L,
        modelout = model_file,
        ...
      ),
      warning = function(w) {
        warning_messages <<- unique(c(warning_messages, conditionMessage(w)))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      messages <- unique(c(conditionMessage(e), warning_messages))
      messages <- messages[!is.na(messages) & nzchar(messages)]
      structure(
        list(
          message = if (length(messages)) {
            paste(messages, collapse = "\n")
          } else {
            paste(class(e), collapse = "/")
          },
          error_message = conditionMessage(e),
          warning_messages = warning_messages,
          condition_class = class(e)
        ),
        class = c("mplus_modeler_failure", "list")
      )
    }
  )
}

#' Identify Mplus warnings that make a returned estimate unsuitable for analysis.
#'
#' Mplus can write a parameter table even when it warns about a nonconverged or
#' inadmissible solution. These prespecified patterns are treated as critical;
#' other warnings are retained for audit but do not automatically exclude a
#' result from performance summaries.
mplus_warning_diagnostics <- function(warnings) {
  warning_lines <- mplus_message_lines(warnings)
  critical_pattern <- paste(
    c(
      "did not converge", "not converge", "did not terminate normally",
      "not positive definite", "non-positive definite", "negative residual variance",
      "singular", "ill-conditioned", "hessian", "standard errors?.*(could not|not)"
    ),
    collapse = "|"
  )
  critical <- grepl(critical_pattern, warning_lines, ignore.case = TRUE, perl = TRUE)
  diagnostics <- mplus_diagnostics_template()
  diagnostics$mplus_warning_count <- as.integer(length(warning_lines))
  diagnostics$mplus_critical_warning <- any(critical)
  diagnostics$mplus_critical_warning_detail <- if (any(critical)) {
      paste(warning_lines[critical], collapse = "\n")
    } else {
      NA_character_
    }
  diagnostics
}

#' Add fitted Mplus latent covariance and boundary diagnostics.
populate_mplus_latent_diagnostics <- function(
    diagnostics, pars, xvar,
    latent_covariance_blocks = list(),
    boundary_tolerance = sqrt(.Machine$double.eps)) {
  param_header <- toupper(as.character(pars$paramHeader))
  param <- toupper(as.character(pars$param))
  estimates <- suppressWarnings(as.numeric(pars$est))
  # Mplus reports exogenous latent variances under `Variances` and variances
  # of endogenous latent variables under `Residual.Variances`. Both are model
  # variance parameters and both can reveal a boundary-adjacent solution.
  variance_parameter <- grepl("VARIANCES$", param_header)

  variance_rows <- which(variance_parameter)
  variance_values <- estimates[variance_rows]
  finite_variances <- variance_values[is.finite(variance_values)]
  diagnostics$mplus_min_reported_variance <- if (length(finite_variances)) {
    min(finite_variances)
  } else {
    NA_real_
  }
  diagnostics$mplus_boundary_variance_count <- as.integer(sum(
    is.finite(variance_values) & variance_values <= boundary_tolerance
  ))
  diagnostics$mplus_nonpositive_variance_count <- as.integer(sum(
    is.finite(variance_values) & variance_values <= 0
  ))

  xvar <- toupper(xvar)
  slope_rows <- which(variance_parameter & param == xvar)
  if (length(slope_rows) == 1L) {
    slope_variance <- estimates[slope_rows]
    diagnostics$mplus_fitted_latent_slope_variance <- slope_variance
    diagnostics$mplus_fitted_latent_slope_sd <- if (
      is.finite(slope_variance) && slope_variance > 0
    ) sqrt(slope_variance) else NA_real_
  }

  for (block_name in names(latent_covariance_blocks)) {
    variables <- toupper(latent_covariance_blocks[[block_name]])
    if (length(variables) != 2L || anyNA(variables)) next
    intercept <- variables[[1L]]
    slope <- variables[[2L]]
    intercept_rows <- which(
      variance_parameter & param == intercept
    )
    block_slope_rows <- which(
      variance_parameter & param == slope
    )
    covariance_rows <- which(
      (param_header == paste0(intercept, ".WITH") & param == slope) |
        (param_header == paste0(slope, ".WITH") & param == intercept)
    )
    if (length(intercept_rows) != 1L || length(block_slope_rows) != 1L ||
        length(covariance_rows) != 1L) next
    covariance <- matrix(
      c(
        estimates[intercept_rows], estimates[covariance_rows],
        estimates[covariance_rows], estimates[block_slope_rows]
      ),
      nrow = 2L,
      byrow = TRUE
    )
    eigenvalues <- tryCatch(
      eigen(covariance, symmetric = TRUE, only.values = TRUE)$values,
      error = function(e) rep(NA_real_, 2L)
    )
    prefix <- if (identical(block_name, "latent") || !nzchar(block_name)) {
      "mplus_latent_"
    } else {
      paste0("mplus_", block_name, "_latent_")
    }
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
    diagnostics[[paste0(prefix, "covariance_boundary")]] <-
      any(!is.finite(eigenvalues)) || any(eigenvalues <= boundary_tolerance) ||
      any(diag(covariance) <= boundary_tolerance)
  }
  diagnostics
}

#' Extract a focal Mplus regression and nonrecoverable fit diagnostics.
#'
#' The reported estimate/SE are multiplied by `reporting_scale`; their raw
#' counterparts are retained separately. Named latent covariance blocks parse
#' both Mplus `Variances` and `Residual.Variances` parameter headers, since an
#' endogenous latent slope is listed under the latter.
extract_mplus_stats <- function(output_file, yvar, xvar, ci_multiplier = stats::qnorm(0.975),
                                reporting_scale = NULL, type = "reg",
                                latent_covariance_blocks = list()) {
  if (type != "reg") {
    stop("`extract_mplus_stats` is only set up for regressions currently")
  }
  
  results <- tryCatch(
    MplusAutomation::readModels(output_file),
    error = function(e) return(NULL)
  )
  if (is.null(results)) {
    return(dplyr::bind_cols(tibble::tibble(
      estimate = NA_real_,
      se = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      status_code = 2L,
      mx_issue_class = "mplus_null_output_file",
      mx_issue_detail = NA_character_
    ), mplus_diagnostics_template()))
  }
  
  errors <- results$errors
  if (length(errors)) {
    parsed_errors <- unlist(lapply(errors, paste, collapse = " "))
    return(dplyr::bind_cols(tibble::tibble(
      estimate = NA_real_,
      se = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      status_code = 1L,
      mx_issue_class = "mplus_error_thrown",
      mx_issue_detail = paste(parsed_errors, collapse = "\n")
    ), mplus_diagnostics_template()))
  } else {
    warning_diagnostics <- mplus_warning_diagnostics(results$warnings)
    
    pars <- results$parameters$unstandardized
    yvar <- toupper(yvar)
    xvar <- toupper(xvar)
    
    scale_u1 <- if (is.null(reporting_scale)) {
      stop("no default reporting scale for mplus models currently")
      # sqrt(with(pars, est[paramHeader == "Variances" & param == xvar]))
    } else {
      reporting_scale <- as.numeric(reporting_scale[[1]])
      if (!is.finite(reporting_scale) || reporting_scale <= 0) {
        stop("`reporting_scale` must be finite and positive.")
      }
      reporting_scale
    }
    
    operator <- switch(
      type,
      reg = ".ON",
      var = ".WITH"
    )
    rows <- which(pars$paramHeader == paste0(yvar, operator) & pars$param == xvar)
    warning_diagnostics$mplus_target_parameter_count <- as.integer(length(rows))
    warning_diagnostics <- populate_mplus_latent_diagnostics(
      warning_diagnostics,
      pars = pars,
      xvar = xvar,
      latent_covariance_blocks = latent_covariance_blocks
    )
    if (length(rows) != 1L) {
      return(dplyr::bind_cols(tibble::tibble(
        estimate = NA_real_,
        se = NA_real_,
        ci_low = NA_real_,
        ci_high = NA_real_,
        status_code = 3L,
        mx_issue_class = "mplus_target_parameter_not_unique",
        mx_issue_detail = sprintf("Expected one target parameter; found %d.", length(rows))
      ), warning_diagnostics))
    }
    raw_est <- suppressWarnings(as.numeric(pars[rows, "est"]))
    raw_se <- suppressWarnings(as.numeric(pars[rows, "se"]))
    warning_diagnostics$mplus_raw_focal_estimate <- raw_est
    warning_diagnostics$mplus_raw_focal_se <- raw_se
    est <- raw_est * scale_u1
    se <- raw_se * scale_u1
    if (!is.finite(est) || !is.finite(se) || se <= 0) {
      return(dplyr::bind_cols(tibble::tibble(
        estimate = NA_real_,
        se = NA_real_,
        ci_low = NA_real_,
        ci_high = NA_real_,
        status_code = 3L,
        mx_issue_class = "mplus_invalid_target_estimate_or_se",
        mx_issue_detail = "The target estimate or standard error was non-finite or non-positive."
      ), warning_diagnostics))
    }
    
    return(dplyr::bind_cols(tibble::tibble(
      estimate = est,
      se = se,
      ci_low = est - ci_multiplier * se,
      ci_high = est + ci_multiplier * se,
      status_code = 0L,
      mx_issue_class = "ok",
      mx_issue_detail = "ok"
    ), warning_diagnostics))
  }
}

fit_mplus_blup_predictor <- function(level1_data, level2_data, join_by = dplyr::join_by("cid_chr" == "id"),
                                     title = "blup_predictor_msem", outcome_variable = "y", slope_only = FALSE,
                                     within_component = "x", between_component = "z", cluster_id = "cid",
                                     random_slope_label = "s", model_file = tempfile(title, fileext = ".inp"),
                                     output_file = gsub(".inp", ".out", model_file), reporting_scale = NULL) {
  
  if (!grepl(".inp$", model_file)) {
    stop("`model_file` must end in .inp")
  }
  if (model_file != gsub(".out", ".inp", model_file)) {
    stop("`model_file` and `output_file` must have the same name")
  }
  if (!(join_by$x %in% names(level1_data) && join_by$y %in% names(level2_data))) {
    stop("invalid join_by statement -- data missing columns")
  }
  
  data <- dplyr::left_join(level1_data, level2_data, by = join_by)
  
  model <- MplusAutomation::mplusObject(
    TITLE = title,
    ANALYSIS = "TYPE = TWOLEVEL RANDOM",
    VARIABLE = glue::glue(
      "WITHIN = {within_component};
       BETWEEN = {between_component};
       CLUSTER IS {cluster_id};"
    ),
    MODEL = glue::glue(
      "%WITHIN%
       {random_slope_label} | {outcome_variable} ON {within_component};
       %BETWEEN%
       {between_component} ON {if (slope_only) '' else outcome_variable} {random_slope_label};
       {outcome_variable} WITH {random_slope_label};"
    ),
    rdata = data
  )
  
  fit <- run_mplus_modeler_writable_tmp(model, model_file)
  
  if (is.null(fit) || inherits(fit, "mplus_modeler_failure")) {
    failure_detail <- if (inherits(fit, "mplus_modeler_failure")) {
      fit$message
    } else {
      "mplusModeler returned NULL without throwing an R error."
    }
    out <- dplyr::bind_cols(tibble::tibble(
      estimate = NA_real_,
      se = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      status_code = 2L,
      mx_issue_class = "mplus_null_fit",
      mx_issue_detail = failure_detail
    ), mplus_diagnostics_template())
  } else {
    out <- extract_mplus_stats(
      output_file,
      xvar = random_slope_label,
      yvar = between_component,
      reporting_scale = reporting_scale,
      latent_covariance_blocks = list(
        latent = c(outcome_variable, random_slope_label)
      )
    )
  }
  
  return(out)
  
}

fit_mplus_dual_process <- function(proc1_data, proc2_data, title = "dual_process_sem", outcome1_var = "y", outcome2_var = "q",
                                   cluster_id = "cid", time_prefix = "t", time_index_var = "trial_index", round_time = 5,
                                   time_value_var = "x", slope_only = FALSE, model_file = tempfile(title, fileext = ".inp"),
                                   output_file = gsub(".inp", ".out", model_file), reporting_scale = NULL) {
  
  if (is.null(time_value_var)) time_value_var <- time_index_var
  
  prepare_process_data <- function(data, outcome, cluster, idx, time_val, time_prefix) {
    
    cols <- c(outcome, cluster, idx, time_val)
    
    if (!all(cols%in% colnames(data))) {
      stop("dual process data missing specified columns: ", 
           setdiff(cols, colnames(data)))
    }
    
    wide <- data %>%
      dplyr::select(!!cluster, !!idx, !!outcome) %>%
      dplyr::arrange(!!cluster, !!idx) %>%
      tidyr::pivot_wider(
        names_from = !!idx, 
        values_from = !!outcome, 
        names_prefix = paste0(time_prefix, outcome)
      )
    
    time_map <- data %>%
      dplyr::select(!!idx, !!time_val) %>%
      dplyr::distinct() %>%
      dplyr::arrange(!!idx)
    
    col_names <- paste0(time_prefix, outcome, time_map[[idx]])
    at_values  <- round(time_map[[time_val]], round_time)
    
    model_string <- paste0(col_names, "@", at_values)
    time_var_range <- paste(col_names[1], col_names[length(col_names)], sep = "-")
    
    return(list(wide = wide, model_string = model_string, time_var_range = time_var_range))
    
  }
  
  y_out <- prepare_process_data(
    data = proc1_data, 
    outcome= outcome1_var, 
    cluster = cluster_id, 
    idx = time_index_var,
    time_val = time_value_var,
    time_prefix = time_prefix
  )
  
  q_out <- prepare_process_data(
    data = proc2_data, 
    outcome= outcome2_var, 
    cluster = cluster_id, 
    idx = time_index_var,
    time_val = time_value_var,
    time_prefix = time_prefix
  )
  
  data <- dplyr::left_join(y_out$wide, q_out$wide, by = cluster_id)
  
  model <- MplusAutomation::mplusObject(
    TITLE = title,
    MODEL = glue::glue(
      "i1 s1 | {paste(y_out$model_string, collapse = '\n')};
       i2 s2 | {paste(q_out$model_string, collapse = '\n')};
       s2 ON {if (slope_only) 's1' else paste0(c('i1', 's1'), collapse = ' ')};
       {paste(y_out$time_var_range)} (1);
       {paste(q_out$time_var_range)} (2);
       i1 WITH s1;
       i2 WITH s2;"
    ),
    rdata = data
  )
  
  fit <- run_mplus_modeler_writable_tmp(model, model_file)
  
  if (is.null(fit) || inherits(fit, "mplus_modeler_failure")) {
    failure_detail <- if (inherits(fit, "mplus_modeler_failure")) {
      fit$message
    } else {
      "mplusModeler returned NULL without throwing an R error."
    }
    out <- dplyr::bind_cols(tibble::tibble(
      estimate = NA_real_,
      se = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      status_code = 2L,
      mx_issue_class = "mplus_null_fit",
      mx_issue_detail = failure_detail
    ), mplus_diagnostics_template())
  } else {
    out <- extract_mplus_stats(
      output_file,
      xvar = "s1",
      yvar = "s2",
      reporting_scale = reporting_scale,
      latent_covariance_blocks = list(
        predictor = c("i1", "s1"),
        outcome = c("i2", "s2")
      )
    )
  }
  
  return(out)
  
}
