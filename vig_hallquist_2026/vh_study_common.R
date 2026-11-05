#' Shared study machinery that is not part of the general score helpers.

# vig_hallquist_truth <- function(condition) {
#   condition$gamma * sqrt(condition$tau1)
# }

simulate_data_blup_as_outcome <- function(condition) {
  sim <- simulate_dataset(
    n_id = condition$num_clus,
    mean_n_trial = condition$mean_clus_size,
    params = list(
      tau0 = fixed_params$tau0,
      rho = condition$rho_residual,
      gamma_x_on_slope = condition$beta1w,
      beta_0 = fixed_params$gamma0_outcome,
      beta_z = fixed_params$gamma1_outcome
    ),
    tau1 = condition$tau1,
    sigma = condition$sigma,
    has_random_slope = TRUE,
    balanced = balance_mode_to_sim_arg(condition$balance_mode),
    min_n_trial = condition$min_clus_size,
    highly_unbalanced_min_n_trial = condition$highly_unbalanced_min_clus_size,
    highly_unbalanced_power = condition$highly_unbalanced_power,
    r_spec = condition_to_r_spec(condition)
  )
  
  list(
    lv1 = tibble::tibble(
      cid = sim$dat$id,
      cid_chr = as.character(sim$dat$id),
      trial_index = sim$dat$trial_index,
      x = sim$dat$z, # level 1 predictor is now "x" across the board
      y = sim$dat$y
    ),
    lv2_true = tibble::tibble(
      id = sim$id_df$id,
      w = sim$id_df$x, # level 2 predictor of slope is now "w" across the board
      z = NA_real_,
      true_u0 = sim$id_df$true_intercept_dev,
      true_u1 = sim$id_df$true_slope_dev
    ), 
    R_list = sim$R_list,
    r_spec = sim$r_spec,
    mean_realized_trials = sim$mean_realized_trials,
    min_realized_trials = sim$min_realized_trials,
    prop_ids_leq_2_trials = sim$prop_ids_leq_2_trials,
    prop_ids_leq_3_trials = sim$prop_ids_leq_3_trials,
    study_structure = "w",
    balance_mode = condition$balance_mode
  )
}

simulate_data_blup_as_predictor <- function(condition) {
  covu <- make_random_effect_covariance(
    intercept_variance = fixed_params$tau0^2,
    slope_variance = condition$tau1^2,
    intercept_slope_correlation = condition$marginal_rho
  )
  
  r_spec <- normalize_r_spec(r_spec = condition_to_r_spec(condition))
  
  if (!is.finite(condition$min_clus_size) || condition$min_clus_size < 1L) {
    stop("`min_clus_size` must be an integer >= 1.")
  }
  
  balanced <- balance_mode_to_sim_arg(condition$balance_mode)
  cluster_sizes <- if (isTRUE(balanced) ||
      (is.character(balanced) && balanced == "balanced")) {
    rep(as.integer(condition$mean_clus_size), as.integer(condition$num_clus))
  } else if (isFALSE(balanced)) {
    pmax(
      as.integer(condition$min_clus_size),
      as.integer(round(
        stats::runif(
          as.integer(condition$num_clus),
          min = 0.6,
          max = 1.4
        ) * as.integer(condition$mean_clus_size)
      ))
    )
  } else {
    stop("Informative cluster-size imbalance is not defined for the BLUP-as-predictor study.")
  }
  cid <- rep(seq_len(as.integer(condition$num_clus)), cluster_sizes)

  x <- unlist(lapply(
    cluster_sizes,
    function(n_i) make_reliability_time_design(n_i)[, "slope"]
  ), use.names = FALSE)

  u <- draw_random_effects(
    n_id = condition$num_clus,
    tau0 = fixed_params$tau0,
    tau1 = condition$tau1,
    rho = condition$marginal_rho
  )
  y <- fixed_params$gamma0_predictor + fixed_params$gamma1_predictor * x + rowSums(cbind(1, x) * u[cid, , drop = FALSE]) +
    draw_level1_residuals(cluster_sizes, sigma = condition$sigma, condition = condition)
  sigma_e <- fixed_params$z_variance -
    drop(t(c(condition$beta1z, condition$beta2z)) %*%
      covu %*% c(condition$beta1z, condition$beta2z))
  if ("outcome_residual_variance" %in% names(condition) &&
      abs(sigma_e - condition$outcome_residual_variance[[1]]) > 1e-10) {
    stop("Stored and recomputed BLUP-as-predictor residual variances disagree.")
  }
  if (!is.finite(sigma_e) || sigma_e <= 0) {
    stop("Calibrated BLUP-as-predictor outcome residual variance must be positive.")
  }
  z <- fixed_params$beta0z + condition$beta1z * u[, 1] + condition$beta2z * u[, 2] +
    stats::rnorm(as.integer(condition$num_clus), sd = sqrt(sigma_e))

  R_list <- stats::setNames(
    lapply(cluster_sizes, make_R_matrix, sigma = condition$sigma, r_spec = r_spec),
    unique(cid)
  )
  
  list(
    lv1 = tibble::tibble(
      cid = cid,
      cid_chr = as.character(cid),
      trial_index = ave(seq_along(cid), cid, FUN = seq_along),
      x = x,
      y = y
    ),
    lv2_true = tibble::tibble(
      id = as.character(seq_len(as.integer(condition$num_clus))),
      w = NA_real_,
      z = z,
      true_u0 = u[, 1],
      true_u1 = u[, 2]
    ), 
    R_list = R_list,
    r_spec = r_spec,
    mean_realized_trials = mean(cluster_sizes),
    min_realized_trials = min(cluster_sizes),
    prop_ids_leq_2_trials = mean(cluster_sizes <= 2L),
    prop_ids_leq_3_trials = mean(cluster_sizes <= 3L),
    study_structure = "z",
    balance_mode = condition$balance_mode
  )
}

simulate_data_dual_blup <- function(condition) {
  n_clus <- as.integer(condition$num_clus[[1]])
  make_sizes <- function(mean_size) {
    balance_mode <- as.character(condition$balance_mode[[1]])
    if (identical(balance_mode, "balanced")) {
      return(rep(as.integer(mean_size), n_clus))
    }
    if (identical(balance_mode, "unbalanced")) {
      return(pmax(
        as.integer(condition$min_clus_size[[1]]),
        as.integer(round(
          stats::runif(n_clus, min = 0.6, max = 1.4) *
            as.integer(mean_size)
        ))
      ))
    }
    stop("Informative cluster-size imbalance is not implemented for Study 3.")
  }

  cluster_sizes_y <- make_sizes(condition$mean_clus_size_y[[1]])
  cluster_sizes_q <- make_sizes(condition$mean_clus_size_q[[1]])
  cid_y <- rep(seq_len(n_clus), cluster_sizes_y)
  cid_q <- rep(seq_len(n_clus), cluster_sizes_q)
  x_y <- unlist(lapply(
    cluster_sizes_y,
    function(n_i) make_reliability_time_design(n_i)[, "slope"]
  ), use.names = FALSE)
  x_q <- unlist(lapply(
    cluster_sizes_q,
    function(n_i) make_reliability_time_design(n_i)[, "slope"]
  ), use.names = FALSE)

  y_effects <- draw_random_effects(
    n_id = n_clus,
    tau0 = fixed_params$tau0,
    tau1 = condition$tau1_y[[1]],
    rho = condition$marginal_rho[[1]]
  )
  q_residual_effects <- draw_random_effects(
    n_id = n_clus,
    tau0 = fixed_params$tau0,
    tau1 = condition$tau1_residual_q[[1]],
    rho = condition$rho_residual_q[[1]]
  )
  q_slope <- condition$theta0[[1]] * y_effects[, 1L] +
    condition$theta1[[1]] * y_effects[, 2L] +
    q_residual_effects[, 2L]

  y <- fixed_params$gamma0_process_y +
    fixed_params$gamma1_process_y * x_y +
    y_effects[cid_y, 1L] +
    y_effects[cid_y, 2L] * x_y +
    draw_level1_residuals(
      cluster_sizes_y,
      sigma = condition$sigma_y[[1]],
      condition = condition
    )
  q <- fixed_params$gamma0_process_q +
    fixed_params$gamma1_process_q * x_q +
    q_residual_effects[cid_q, 1L] +
    q_slope[cid_q] * x_q +
    draw_level1_residuals(
      cluster_sizes_q,
      sigma = condition$sigma_q[[1]],
      condition = condition
    )
  r_spec <- normalize_r_spec(condition_to_r_spec(condition))
  ids <- as.character(seq_len(n_clus))

  list(
    lv1_y = tibble::tibble(
      cid = factor(cid_y),
      cid_chr = as.character(cid_y),
      trial_index = ave(seq_along(cid_y), cid_y, FUN = seq_along),
      x = x_y,
      y = y
    ),
    lv1_q = tibble::tibble(
      cid = factor(cid_q),
      cid_chr = as.character(cid_q),
      trial_index = ave(seq_along(cid_q), cid_q, FUN = seq_along),
      x = x_q,
      q = q
    ),
    lv2_true = tibble::tibble(
      id = ids,
      true_y0 = y_effects[, 1L],
      true_y1 = y_effects[, 2L],
      true_q0 = q_residual_effects[, 1L],
      true_q1 = q_slope
    ),
    R_list_y = stats::setNames(
      lapply(
        cluster_sizes_y,
        make_R_matrix,
        sigma = condition$sigma_y[[1]],
        r_spec = r_spec
      ),
      ids
    ),
    R_list_q = stats::setNames(
      lapply(
        cluster_sizes_q,
        make_R_matrix,
        sigma = condition$sigma_q[[1]],
        r_spec = r_spec
      ),
      ids
    ),
    mean_realized_trials_y = mean(cluster_sizes_y),
    mean_realized_trials_q = mean(cluster_sizes_q),
    min_realized_trials_y = min(cluster_sizes_y),
    min_realized_trials_q = min(cluster_sizes_q),
    study_structure = "dual_process",
    balance_mode = condition$balance_mode
  )
}

condition_to_r_spec <- function(condition) {
  r_structure <- if ("r_structure" %in% names(condition)) {
    as.character(condition$r_structure[[1]])
  } else if ("residual_structure" %in% names(condition)) {
    as.character(condition$residual_structure[[1]])
  } else {
    "iid"
  }

  switch(
    tolower(r_structure),
    iid = list(structure = "iid"),
    ar1 = {
      rho <- if ("r_rho" %in% names(condition)) condition$r_rho[[1]] else condition$residual_rho[[1]]
      list(structure = "ar1", rho = as.numeric(rho))
    },
    toeplitz = {
      correlations <- if ("r_correlations" %in% names(condition)) condition$r_correlations[[1]] else condition$residual_correlations[[1]]
      list(structure = "toeplitz", correlations = as.numeric(correlations))
    },
    stop("Unsupported residual structure: ", r_structure)
  )
}


condition_to_nlme_correlation <- function(condition, cluster_var = "id") {
  r_structure <- if ("r_structure" %in% names(condition)) {
    as.character(condition$r_structure[[1]])
  } else {
    "iid"
  }
  
  switch(
    r_structure,
    iid = NULL,
    ar1 = {
      if (!requireNamespace("nlme", quietly = TRUE)) {
        stop("The `nlme` package is required for AR(1) Stage-1 residual covariance fits.")
      }
      nlme::corAR1(
        form = stats::as.formula(paste0("~trial_index | ", cluster_var))
      )
    },
    stop("Unsupported residual structure for nlme Stage-1 fit: ", r_structure)
  )
}

condition_uses_non_iid_R <- function(condition) {
  !identical(normalize_r_spec(condition_to_r_spec(condition))$structure, "iid")
}


add_trial_index <- function(data, cluster_var, index_var = "trial_index") {
  cluster_ids <- as.character(data[[cluster_var]])
  data[[index_var]] <- ave(seq_len(nrow(data)), cluster_ids, FUN = seq_along)
  data
}

draw_level1_residuals <- function(cluster_sizes, sigma, condition) {
  r_spec <- condition_to_r_spec(condition)
  unlist(lapply(cluster_sizes, function(n_i) {
    draw_residuals_from_R(make_R_matrix(n_i, sigma = sigma, r_spec = r_spec))
  }), use.names = FALSE)
}

# TODO: this is only for BLUPs as predictors right now (see README for models)
fit_stage1 <- function(fixed, random, data, condition, cluster_var) {
  if (condition_uses_non_iid_R(condition)) {
    if (!requireNamespace("nlme", quietly = TRUE)) {
      return(NULL)
    }
    return(safe_lme(
      fixed = fixed,
      random = random,
      data = data,
      correlation = condition_to_nlme_correlation(condition, cluster_var = cluster_var),
      method = "REML",
      control = nlme::lmeControl(returnObject = TRUE, msMaxIter = 100L, opt = "optim")
    ))
  }

  fixed_txt <- paste(deparse(fixed), collapse = " ")
  random_txt <- sub("^~", "", paste(deparse(random), collapse = " "))
  safe_lmer(stats::as.formula(paste0(fixed_txt, " + (", random_txt, ")")), data = data)
}

extract_centered_slope_eb <- function(fit_obj, ordered_ids) {
  if (is.null(fit_obj)) {
    return(tibble::tibble(id = ordered_ids, centered_u1_eb = NA_real_))
  }

  re_df <- tryCatch({
    if (inherits(fit_obj, "lme")) {
      nlme::ranef(fit_obj)
    } else {
      lme4::ranef(fit_obj)[[1]]
    }
  }, error = function(e) NULL)

  if (is.null(re_df)) {
    return(tibble::tibble(id = ordered_ids, centered_u1_eb = NA_real_))
  }

  tibble::tibble(id = rownames(re_df), centered_u1_eb = as.numeric(re_df[[ncol(re_df)]]))
}

combine_dual_stage1_diagnostics <- function(y_diag, q_diag) {
  extreme_by_absolute_value <- function(x, y) {
    values <- c(x, y)
    values <- values[is.finite(values)]
    if (length(values) == 0L) NA_real_ else values[[which.max(abs(values))]]
  }
  max_finite <- function(x, y) {
    values <- c(x, y)
    values <- values[is.finite(values)]
    if (length(values) == 0L) NA_real_ else max(values)
  }

  y_problem <- isTRUE(y_diag$stage1_singular_problem[[1]])
  q_problem <- isTRUE(q_diag$stage1_singular_problem[[1]])
  problem_parts <- c(
    if (y_problem) paste0("y:", y_diag$stage1_problem_detail[[1]]) else NA_character_,
    if (q_problem) paste0("q:", q_diag$stage1_problem_detail[[1]]) else NA_character_
  )

  dplyr::bind_cols(
    tibble::tibble(
      stage1_singular_problem = y_problem || q_problem,
      stage1_problem_detail = if (y_problem || q_problem) {
        compact_message(problem_parts)
      } else {
        "ok"
      },
      stage1_lmer_singular =
        isTRUE(y_diag$stage1_lmer_singular[[1]]) ||
        isTRUE(q_diag$stage1_lmer_singular[[1]]),
      stage1_re_corr = extreme_by_absolute_value(
        y_diag$stage1_re_corr[[1]],
        q_diag$stage1_re_corr[[1]]
      ),
      stage1_eb_corr = extreme_by_absolute_value(
        y_diag$stage1_eb_corr[[1]],
        q_diag$stage1_eb_corr[[1]]
      ),
      stage1_design_kappa = max_finite(
        y_diag$stage1_design_kappa[[1]],
        q_diag$stage1_design_kappa[[1]]
      )
    ),
    y_diag %>% dplyr::rename_with(~ paste0("stage1_y_", sub("^stage1_", "", .x))),
    q_diag %>% dplyr::rename_with(~ paste0("stage1_q_", sub("^stage1_", "", .x)))
  )
}

# study method vectors in study specific files
study_methods_for_condition <- function(condition) {
  study_key <- as.character(condition$study[[1]])
  switch(
    study_key,
    study0 = study1_methods(),
    study1 = study1_methods(),
    study2 = study2_methods(),
    study3 = study3_methods(),
    study4 = study4_methods(),
    study5 = study5_methods(),
    study1v2 = study1_methods(),
    study2v2 = study2_methods(),
    study3v2 = study3_methods(),
    study4v2 = study4_methods(),
    iccbridge = study2_methods(),
    stop("Unsupported study key: ", study_key)
  )
}

make_failed_result <- function(condition, methods, truth) {
  dplyr::bind_cols(
    tibble::tibble(
      study = condition$study,
      method = methods,
      estimate = NA_real_,
      se = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      status_code = NA_integer_,
      truth = truth
    ),
    empty_stage1_diagnostics()[rep(1L, length(methods)), , drop = FALSE]
  )
}

add_study_result_context <- function(results, condition, truth) {
  dplyr::mutate(results, study = condition$study, truth = truth)
}

#' Summarize fitted Stage-1 measurement quality within one replication.
#'
#' These quantities are deliberately saved before the cluster-level score data
#' are discarded. They audit the full 2-by-2 Lai loading/residual matrices,
#' posterior uncertainty, likelihood-only corrected-score uncertainty, and the
#' empirical relationship between simulated random slopes and their estimated
#' scores. Posterior reliabilities are computed from the fitted `G`, rather
#' than copied from the population calibration manifest.
summarize_stage1_measurement_diagnostics <- function(
    stage1_scores, G_hat, suffix = NULL, data_prefix = "",
    true_slope_col = NULL, blup_slope_col = paste0(data_prefix, "u1_eb"),
    corrected_slope_col = if (nzchar(data_prefix)) {
      paste0(data_prefix, "corrected_slope")
    } else {
      "corrected_slope_full"
    }) {
  result_prefix <- paste0(
    "stage1",
    if (!is.null(suffix) && nzchar(suffix)) paste0("_", suffix) else "",
    "_"
  )
  named_value <- function(name, value) {
    stats::setNames(list(value), paste0(result_prefix, name))
  }
  finite_values <- function(column) {
    if (!(column %in% names(stage1_scores))) return(numeric())
    values <- suppressWarnings(as.numeric(stage1_scores[[column]]))
    values[is.finite(values)]
  }
  finite_mean <- function(x) if (length(x)) mean(x) else NA_real_
  finite_sd <- function(x) if (length(x) > 1L) stats::sd(x) else NA_real_
  finite_min <- function(x) if (length(x)) min(x) else NA_real_
  finite_max <- function(x) if (length(x)) max(x) else NA_real_

  diagnostics <- list()
  append_summary <- function(label, column) {
    values <- finite_values(column)
    diagnostics <<- c(
      diagnostics,
      named_value(paste0(label, "_mean"), finite_mean(values)),
      named_value(paste0(label, "_min"), finite_min(values)),
      named_value(paste0(label, "_max"), finite_max(values))
    )
  }

  for (entry in c("11", "12", "21", "22")) {
    append_summary(
      paste0("lambda", entry),
      paste0(data_prefix, "lambda", entry)
    )
  }
  for (entry in c("11", "12", "22")) {
    append_summary(
      paste0("theta", entry),
      paste0(data_prefix, "theta", entry)
    )
    append_summary(
      paste0("posterior_variance", entry),
      paste0(data_prefix, "postvar", entry)
    )
    corrected_error_column <- paste0(data_prefix, "ols_var", entry)
    if (!(corrected_error_column %in% names(stage1_scores))) {
      corrected_error_column <- paste0(data_prefix, "mle_var", entry)
    }
    append_summary(
      paste0("corrected_score_error_covariance", entry),
      corrected_error_column
    )
  }

  marginal_reliability <- residualized_reliability <- rep(
    NA_real_, nrow(stage1_scores)
  )
  if (is.matrix(G_hat) && all(dim(G_hat) >= 2L)) {
    post11 <- suppressWarnings(as.numeric(
      stage1_scores[[paste0(data_prefix, "postvar11")]]
    ))
    post12 <- suppressWarnings(as.numeric(
      stage1_scores[[paste0(data_prefix, "postvar12")]]
    ))
    post22 <- suppressWarnings(as.numeric(
      stage1_scores[[paste0(data_prefix, "postvar22")]]
    ))
    if (is.finite(G_hat[2L, 2L]) && G_hat[2L, 2L] > 0) {
      marginal_reliability <- 1 - post22 / G_hat[2L, 2L]
    }
    if (is.finite(G_hat[1L, 1L]) && G_hat[1L, 1L] > 0) {
      residualization <- G_hat[1L, 2L] / G_hat[1L, 1L]
      residualized_variance <- G_hat[2L, 2L] -
        G_hat[1L, 2L]^2 / G_hat[1L, 1L]
      residualized_postvar <- post22 - 2 * residualization * post12 +
        residualization^2 * post11
      if (is.finite(residualized_variance) && residualized_variance > 0) {
        residualized_reliability <-
          1 - residualized_postvar / residualized_variance
      }
    }
  }
  for (metric in c("marginal_slope", "residualized_slope")) {
    values <- if (identical(metric, "marginal_slope")) {
      marginal_reliability
    } else {
      residualized_reliability
    }
    values <- values[is.finite(values)]
    diagnostics <- c(
      diagnostics,
      named_value(
        paste0("fitted_", metric, "_posterior_reliability_mean"),
        finite_mean(values)
      ),
      named_value(
        paste0("fitted_", metric, "_posterior_reliability_sd"),
        finite_sd(values)
      ),
      named_value(
        paste0("fitted_", metric, "_posterior_reliability_min"),
        finite_min(values)
      ),
      named_value(
        paste0("fitted_", metric, "_posterior_reliability_max"),
        finite_max(values)
      )
    )
  }

  score_diagnostics <- function(score_column, label) {
    if (is.null(true_slope_col) ||
        !(true_slope_col %in% names(stage1_scores)) ||
        !(score_column %in% names(stage1_scores))) {
      return(c(
        named_value(paste0(label, "_slope_bias"), NA_real_),
        named_value(paste0(label, "_slope_rmse"), NA_real_)
      ))
    }
    truth <- suppressWarnings(as.numeric(stage1_scores[[true_slope_col]]))
    score <- suppressWarnings(as.numeric(stage1_scores[[score_column]]))
    keep <- is.finite(truth) & is.finite(score)
    error <- score[keep] - truth[keep]
    c(
      named_value(
        paste0(label, "_slope_bias"),
        if (length(error)) mean(error) else NA_real_
      ),
      named_value(
        paste0(label, "_slope_rmse"),
        if (length(error)) sqrt(mean(error^2)) else NA_real_
      )
    )
  }
  diagnostics <- c(
    diagnostics,
    score_diagnostics(blup_slope_col, "blup"),
    score_diagnostics(corrected_slope_col, "corrected_score")
  )

  true_blup_correlation <- NA_real_
  if (!is.null(true_slope_col) &&
      all(c(true_slope_col, blup_slope_col) %in% names(stage1_scores))) {
    truth <- suppressWarnings(as.numeric(stage1_scores[[true_slope_col]]))
    blup <- suppressWarnings(as.numeric(stage1_scores[[blup_slope_col]]))
    keep <- is.finite(truth) & is.finite(blup)
    if (sum(keep) >= 3L && stats::sd(truth[keep]) > 0 &&
        stats::sd(blup[keep]) > 0) {
      true_blup_correlation <- stats::cor(truth[keep], blup[keep])
    }
  }
  diagnostics <- c(
    diagnostics,
    named_value("true_blup_slope_correlation", true_blup_correlation),
    named_value(
      "true_blup_slope_r_squared",
      if (is.finite(true_blup_correlation)) true_blup_correlation^2 else NA_real_
    )
  )

  tibble::as_tibble(diagnostics)
}

add_stage1_estimates <- function(
    results, fit_obj, data, cluster_var, within_var = NULL, R_list = NULL,
    group = NULL, suffix = NULL, stage1_scores = NULL,
    data_prefix = "", true_slope_col = NULL,
    blup_slope_col = paste0(data_prefix, "u1_eb"),
    corrected_slope_col = if (nzchar(data_prefix)) {
      paste0(data_prefix, "corrected_slope")
    } else {
      "corrected_slope_full"
    }) {
  if (is.null(fit_obj)) {
    return(results)
  }
  stage1_components <- extract_stage1_components(
    fit_obj = fit_obj,
    data = data,
    cluster_var = cluster_var,
    within_var = within_var,
    R_list = R_list,
    group = group
  )
  G_hat <- stage1_components$G_hat
  if (is.null(G_hat)) {
    return(results)
  } else {
    result_prefix <- paste0(
      "stage1",
      if (!is.null(suffix) && nzchar(suffix)) paste0("_", suffix) else "",
      "_"
    )
    g_eigenvalues <- tryCatch(
      eigen((G_hat + t(G_hat)) / 2, symmetric = TRUE, only.values = TRUE)$values,
      error = function(e) rep(NA_real_, nrow(G_hat))
    )
    g_diagnostics <- tibble::tibble(
      intercept_variance = G_hat[1, 1],
      slope_variance = G_hat[2, 2],
      intercept_slope_covariance = G_hat[1, 2],
      intercept_slope_correlation = if (
        all(is.finite(G_hat[1:2, 1:2])) &&
          G_hat[1, 1] > 0 && G_hat[2, 2] > 0
      ) {
        G_hat[1, 2] / sqrt(G_hat[1, 1] * G_hat[2, 2])
      } else {
        NA_real_
      },
      latent_covariance_min_eigenvalue = if (any(is.finite(g_eigenvalues))) {
        min(g_eigenvalues, na.rm = TRUE)
      } else {
        NA_real_
      },
      latent_covariance_boundary = any(!is.finite(g_eigenvalues)) ||
        any(g_eigenvalues <= sqrt(.Machine$double.eps), na.rm = TRUE)
    )
    names(g_diagnostics) <- paste0(result_prefix, names(g_diagnostics))
    if (!is.null(stage1_scores)) {
      measurement_diagnostics <- summarize_stage1_measurement_diagnostics(
        stage1_scores = stage1_scores,
        G_hat = G_hat,
        suffix = suffix,
        data_prefix = data_prefix,
        true_slope_col = true_slope_col,
        blup_slope_col = blup_slope_col,
        corrected_slope_col = corrected_slope_col
      )
      g_diagnostics <- dplyr::bind_cols(
        g_diagnostics,
        measurement_diagnostics
      )
    }
    dplyr::bind_cols(results, g_diagnostics)
  }
}
