#' Vig-Hallquist (2026) simulation Study 4: heterogeneous cluster information.

study4_methods <- function() {
  c(
    "oracle_dual",
    "naive_dual_blup",
    "closed_form_dual",
    "fuller_closed_form",
    "fuller_book_preliminary_closed_form",
    "fuller_book_closed_form",
    "fuller_average_measurement",
    "fuller_alpha_stepdown_closed_form",
    "lai_2spa",
    "lai_2spaa",
    "msem"
  )
}

add_study4_method_roles <- function(results) {
  dplyr::mutate(
    results,
    method_role = dplyr::case_when(
      method %in% c("oracle_dual", "msem") ~ "benchmark",
      method == "closed_form_dual" ~ "unshrinking_only_diagnostic",
      method == "naive_dual_blup" ~ "blup_projection_diagnostic",
      method %in% c("fuller_average_measurement", "lai_2spaa") ~
        "average_measurement_approximation",
      TRUE ~ "row_specific_primary"
    )
  )
}

#' Apply the shared VH point/interval eligibility contract to Study 4.
add_study4_analysis_eligibility <- function(results) {
  add_vh_analysis_eligibility(results)
}

#' Allocate one Study 4 information profile to clusters.
#'
#' Counts follow the population profile as closely as possible. Cluster-size
#' labels are then permuted independently of all subsequently drawn latent
#' variables and residuals.
make_study4_cluster_sizes <- function(num_clus, information_profile) {
  num_clus <- as.integer(num_clus[[1]])
  if (!is.finite(num_clus) || num_clus < 2L) {
    stop("`num_clus` must be an integer of at least 2.")
  }
  spec <- study4_profile_spec(information_profile)
  expected_counts <- num_clus * spec$weights
  counts <- floor(expected_counts)
  remainder <- num_clus - sum(counts)
  if (remainder > 0L) {
    fractional <- expected_counts - counts
    if (sum(fractional) <= sqrt(.Machine$double.eps)) {
      chosen <- sample(seq_along(counts), size = remainder, replace = TRUE)
    } else {
      chosen <- sample(
        seq_along(counts),
        size = remainder,
        replace = FALSE,
        prob = fractional
      )
    }
    counts <- counts + tabulate(chosen, nbins = length(counts))
  }
  sample(rep(spec$sizes, counts), size = num_clus, replace = FALSE)
}

#' Simulate the Study 4 heterogeneous-information predictor design.
simulate_data_study4 <- function(condition) {
  n_clus <- as.integer(condition$num_clus[[1]])
  sigma <- as.numeric(condition$sigma[[1]])
  profile <- study4_profile_spec(condition$information_profile[[1]])
  cluster_sizes <- make_study4_cluster_sizes(
    n_clus,
    condition$information_profile[[1]]
  )
  ids <- as.character(seq_len(n_clus))
  Z_list <- lapply(
    cluster_sizes,
    study4_time_design,
    information_matched = profile$information_matched
  )
  R_list <- stats::setNames(
    lapply(
      cluster_sizes,
      make_R_matrix,
      sigma = sigma,
      r_spec = list(structure = "iid")
    ),
    ids
  )

  u <- draw_random_effects(
    n_id = n_clus,
    tau0 = fixed_params$tau0,
    tau1 = condition$tau1[[1]],
    rho = condition$marginal_rho[[1]]
  )
  G_marginal <- make_random_effect_covariance(
    intercept_variance = fixed_params$tau0^2,
    slope_variance = condition$tau1[[1]]^2,
    intercept_slope_correlation = condition$marginal_rho[[1]]
  )
  cluster_reliability <- vapply(
    seq_len(n_clus),
    function(i) {
      V_i <- posterior_random_effect_covariance(
        G_marginal,
        Z_list[[i]],
        R_list[[i]]
      )
      1 - V_i[2L, 2L] / G_marginal[2L, 2L]
    },
    numeric(1L)
  )

  lv1 <- purrr::map_dfr(seq_len(n_clus), function(i) {
    x_i <- Z_list[[i]][, "slope"]
    e_i <- draw_residuals_from_R(R_list[[i]])
    tibble::tibble(
      cid = factor(rep(i, cluster_sizes[[i]]), levels = seq_len(n_clus)),
      cid_chr = ids[[i]],
      trial_index = seq_len(cluster_sizes[[i]]),
      x = x_i,
      y = fixed_params$gamma0_predictor +
        fixed_params$gamma1_predictor * x_i +
        u[i, 1L] + u[i, 2L] * x_i + e_i
    )
  })

  sigma_z2 <- as.numeric(condition$outcome_residual_variance[[1]])
  if (!is.finite(sigma_z2) || sigma_z2 <= 0) {
    stop("Study 4 outcome residual variance must be finite and positive.")
  }
  z <- fixed_params$beta0z +
    condition$beta1z[[1]] * u[, 1L] +
    condition$beta2z[[1]] * u[, 2L] +
    stats::rnorm(n_clus, sd = sqrt(sigma_z2))
  lv2_true <- tibble::tibble(
    id = ids,
    w = NA_real_,
    z = z,
    true_u0 = u[, 1L],
    true_u1 = u[, 2L],
    cluster_size = cluster_sizes,
    cluster_reliability = cluster_reliability
  )

  list(
    lv1 = lv1,
    lv2_true = lv2_true,
    R_list = R_list,
    r_spec = list(structure = "iid"),
    cluster_sizes = cluster_sizes,
    cluster_reliability = cluster_reliability,
    mean_realized_trials = mean(cluster_sizes),
    min_realized_trials = min(cluster_sizes),
    max_realized_trials = max(cluster_sizes),
    prop_ids_leq_3_trials = mean(cluster_sizes <= 3L),
    study_structure = "z",
    balance_mode = "heterogeneous_profile"
  )
}

simulate_study4 <- function(condition) {
  simulate_data_study4(condition)
}

study4_measurement_diagnostics <- function(sim, stage2_df) {
  safe_values <- function(x) {
    as.numeric(x[is.finite(x)])
  }
  safe_mean <- function(x) {
    x <- safe_values(x)
    if (length(x) == 0L) NA_real_ else mean(x)
  }
  safe_sd <- function(x) {
    x <- safe_values(x)
    if (length(x) < 2L) NA_real_ else stats::sd(x)
  }
  safe_range_value <- function(x, which = c("min", "max")) {
    which <- match.arg(which)
    x <- safe_values(x)
    if (length(x) == 0L) {
      return(NA_real_)
    }
    if (identical(which, "min")) min(x) else max(x)
  }
  safe_iqr <- function(x) {
    x <- safe_values(x)
    if (length(x) == 0L) {
      return(NA_real_)
    }
    quartiles <- study4_weighted_quantile(
      x,
      rep(1, length(x)),
      c(0.25, 0.75)
    )
    diff(quartiles)
  }
  matrices_from_columns <- function(data, columns, symmetric = FALSE) {
    if (!all(columns %in% names(data))) {
      return(list())
    }
    lapply(seq_len(nrow(data)), function(i) {
      values <- as.numeric(unlist(
        data[i, columns, drop = FALSE],
        use.names = FALSE
      ))
      if (isTRUE(symmetric)) {
        matrix(
          c(values[[1L]], values[[2L]], values[[2L]], values[[3L]]),
          nrow = 2L,
          byrow = TRUE
        )
      } else {
        matrix(values, nrow = 2L, byrow = TRUE)
      }
    })
  }
  score_group_summary <- function(score, cluster_size) {
    keep <- is.finite(score) & is.finite(stage2_df$true_u1) &
      stage2_df$cluster_size == cluster_size
    errors <- score[keep] - stage2_df$true_u1[keep]
    c(
      bias = if (length(errors) == 0L) NA_real_ else mean(errors),
      rmse = if (length(errors) == 0L) NA_real_ else sqrt(mean(errors^2))
    )
  }

  lambda_matrices <- matrices_from_columns(
    stage2_df,
    c("lambda11", "lambda12", "lambda21", "lambda22")
  )
  theta_matrices <- matrices_from_columns(
    stage2_df,
    c("theta11", "theta12", "theta22"),
    symmetric = TRUE
  )
  ols_covariances <- matrices_from_columns(
    stage2_df,
    c("ols_var11", "ols_var12", "ols_var22"),
    symmetric = TRUE
  )
  average_measurement_slope <- tryCatch(
    prepare_fuller_average_measurement(stage2_df)$data$fuller_average_u1,
    error = function(e) rep(NA_real_, nrow(stage2_df))
  )
  small_size <- min(stage2_df$cluster_size)
  large_size <- max(stage2_df$cluster_size)
  blup_small <- score_group_summary(stage2_df$u1_eb, small_size)
  blup_large <- score_group_summary(stage2_df$u1_eb, large_size)
  corrected_small <- score_group_summary(
    stage2_df$corrected_slope_full,
    small_size
  )
  corrected_large <- score_group_summary(
    stage2_df$corrected_slope_full,
    large_size
  )
  average_small <- score_group_summary(average_measurement_slope, small_size)
  average_large <- score_group_summary(average_measurement_slope, large_size)

  tibble::tibble(
    mean_realized_trials = mean(sim$cluster_sizes),
    min_realized_trials = min(sim$cluster_sizes),
    max_realized_trials = max(sim$cluster_sizes),
    prop_ids_leq_3_trials = mean(sim$cluster_sizes <= 3L),
    realized_reliability_mean = mean(sim$cluster_reliability),
    realized_reliability_sd = safe_sd(sim$cluster_reliability),
    realized_reliability_iqr = safe_iqr(sim$cluster_reliability),
    realized_reliability_min = safe_range_value(
      sim$cluster_reliability,
      "min"
    ),
    realized_reliability_max = safe_range_value(
      sim$cluster_reliability,
      "max"
    ),
    lambda22_mean = safe_mean(stage2_df$lambda22),
    lambda22_sd = safe_sd(stage2_df$lambda22),
    lambda22_min = safe_range_value(stage2_df$lambda22, "min"),
    lambda22_max = safe_range_value(stage2_df$lambda22, "max"),
    lambda_matrix_frobenius_rms_dispersion =
      study4_matrix_rms_dispersion(lambda_matrices),
    theta22_mean = safe_mean(stage2_df$theta22),
    theta22_sd = safe_sd(stage2_df$theta22),
    theta22_min = safe_range_value(stage2_df$theta22, "min"),
    theta22_max = safe_range_value(stage2_df$theta22, "max"),
    theta_matrix_frobenius_rms_dispersion =
      study4_matrix_rms_dispersion(theta_matrices),
    ols_var22_mean = safe_mean(stage2_df$ols_var22),
    ols_var22_sd = safe_sd(stage2_df$ols_var22),
    ols_var22_min = safe_range_value(stage2_df$ols_var22, "min"),
    ols_var22_max = safe_range_value(stage2_df$ols_var22, "max"),
    ols_cov_matrix_frobenius_rms_dispersion =
      study4_matrix_rms_dispersion(ols_covariances),
    score_small_cluster_size = small_size,
    score_large_cluster_size = large_size,
    blup_slope_bias_small = blup_small[["bias"]],
    blup_slope_rmse_small = blup_small[["rmse"]],
    blup_slope_bias_large = blup_large[["bias"]],
    blup_slope_rmse_large = blup_large[["rmse"]],
    corrected_slope_bias_small = corrected_small[["bias"]],
    corrected_slope_rmse_small = corrected_small[["rmse"]],
    corrected_slope_bias_large = corrected_large[["bias"]],
    corrected_slope_rmse_large = corrected_large[["rmse"]],
    average_measurement_slope_bias_small = average_small[["bias"]],
    average_measurement_slope_rmse_small = average_small[["rmse"]],
    average_measurement_slope_bias_large = average_large[["bias"]],
    average_measurement_slope_rmse_large = average_large[["rmse"]]
  )
}

run_study4_rep <- function(condition) {
  truth <- as.numeric(condition$standardized_beta_target[[1]])
  latent_slope_sd <- as.numeric(condition$tau1[[1]])
  sim <- simulate_study4(condition)
  fit_y <- fit_stage1(
    y ~ x,
    random = ~x | cid,
    data = sim$lv1,
    condition = condition,
    cluster_var = "cid"
  )
  if (is.null(fit_y)) {
    return(
      make_failed_result(condition, study4_methods(), truth) %>%
        add_study4_method_roles() %>%
        add_study4_analysis_eligibility()
    )
  }

  stage1_y <- tryCatch(
    get_stage1_eb_components(
      fit_obj = fit_y,
      data = sim$lv1,
      cluster_var = "cid",
      outcome_var = "y",
      within_var = "x"
    ),
    error = function(e) NULL
  )
  corrected_y <- tryCatch(
    get_closed_form_corrected_scores(
      fit_obj = fit_y,
      data = sim$lv1,
      cluster_var = "cid",
      outcome_var = "y",
      within_var = "x"
    ),
    error = function(e) NULL
  )
  if (is.null(stage1_y) || is.null(corrected_y)) {
    return(
      make_failed_result(condition, study4_methods(), truth) %>%
        add_study4_method_roles() %>%
        add_study4_analysis_eligibility()
    )
  }

  stage2_df <- sim$lv2_true %>%
    dplyr::left_join(stage1_y, by = "id") %>%
    dplyr::left_join(
      corrected_y %>%
        dplyr::select(
          id,
          corrected_intercept_full,
          corrected_slope_full,
          ols_var11,
          ols_var12,
          ols_var22
      ),
      by = "id"
    ) %>%
    add_zero_fuller_predictor_outcome_covariance()
  stage1_diag <- get_stage1_diagnostics(fit_y, stage2_df)
  measurement_diag <- study4_measurement_diagnostics(sim, stage2_df)

  results <- dplyr::bind_rows(
    fit_observed_dual(
      stage2_df,
      outcome = "z",
      predictor_u0 = "true_u0",
      predictor_u1 = "true_u1",
      reporting_scale = latent_slope_sd
    ) %>%
      dplyr::filter(se_type == "naive") %>%
      dplyr::transmute(
        method = "oracle_dual",
        estimate, se, ci_low, ci_high, status_code
      ),
    fit_observed_dual(
      stage2_df,
      outcome = "z",
      predictor_u0 = "u0_eb",
      predictor_u1 = "u1_eb",
      reporting_scale = latent_slope_sd
    ) %>%
      dplyr::filter(se_type == "naive") %>%
      dplyr::transmute(
        method = "naive_dual_blup",
        estimate, se, ci_low, ci_high, status_code
      ),
    fit_observed_dual(
      stage2_df,
      outcome = "z",
      predictor_u0 = "corrected_intercept_full",
      predictor_u1 = "corrected_slope_full",
      reporting_scale = latent_slope_sd
    ) %>%
      dplyr::filter(se_type == "naive") %>%
      dplyr::transmute(
        method = "closed_form_dual",
        estimate, se, ci_low, ci_high, status_code
      ),
    fit_fuller_dual_variants(
      stage2_df,
      outcome = "z",
      predictor_u0 = "corrected_intercept_full",
      predictor_u1 = "corrected_slope_full",
      meas11 = "ols_var11",
      meas12 = "ols_var12",
      meas22 = "ols_var22",
      predictor_outcome_meas_cov_u0 =
        "fuller_predictor_outcome_meas_cov_u0",
      predictor_outcome_meas_cov_u1 =
        "fuller_predictor_outcome_meas_cov_u1"
    ) %>%
      rescale_fuller_to_population_sd(latent_slope_sd) %>%
      dplyr::mutate(
        method = unname(c(
          stabilized = "fuller_closed_form",
          fuller_preliminary = "fuller_book_preliminary_closed_form",
          fuller_equations = "fuller_book_closed_form"
        )[fuller_variant])
      ) %>%
      dplyr::select(method, dplyr::everything()),
    fit_fuller_average_measurement(
      stage2_df,
      outcome = "z",
      blup_u0 = "u0_eb",
      blup_u1 = "u1_eb",
      predictor_outcome_meas_cov_u0 =
        "fuller_predictor_outcome_meas_cov_u0",
      predictor_outcome_meas_cov_u1 =
        "fuller_predictor_outcome_meas_cov_u1"
    ) %>%
      rescale_fuller_to_population_sd(latent_slope_sd) %>%
      dplyr::mutate(method = "fuller_average_measurement") %>%
      dplyr::select(method, dplyr::everything()),
    fit_fuller_dual_alpha_stepdown(
      stage2_df,
      outcome = "z",
      predictor_u0 = "corrected_intercept_full",
      predictor_u1 = "corrected_slope_full",
      meas11 = "ols_var11",
      meas12 = "ols_var12",
      meas22 = "ols_var22",
      predictor_outcome_meas_cov_u0 =
        "fuller_predictor_outcome_meas_cov_u0",
      predictor_outcome_meas_cov_u1 =
        "fuller_predictor_outcome_meas_cov_u1"
    ) %>%
      rescale_fuller_to_population_sd(latent_slope_sd) %>%
      dplyr::mutate(method = "fuller_alpha_stepdown_closed_form") %>%
      dplyr::select(method, dplyr::everything()),
    fit_lai_2spa_observed_outcome(
      stage2_df,
      use_average = FALSE,
      u0_start = condition$beta1z[[1]],
      reporting_scale = latent_slope_sd
    ) %>%
      dplyr::mutate(method = "lai_2spa") %>%
      dplyr::select(method, dplyr::everything()),
    fit_lai_2spa_observed_outcome(
      stage2_df,
      use_average = TRUE,
      u0_start = condition$beta1z[[1]],
      reporting_scale = latent_slope_sd
    ) %>%
      dplyr::mutate(method = "lai_2spaa") %>%
      dplyr::select(method, dplyr::everything()),
    fit_mplus_blup_predictor(
      level1_data = sim$lv1,
      level2_data = sim$lv2_true,
      outcome_variable = "y",
      within_component = "x",
      between_component = "z",
      cluster_id = "cid",
      reporting_scale = latent_slope_sd
    ) %>%
      dplyr::mutate(method = "msem") %>%
      dplyr::select(method, dplyr::everything())
  ) %>%
    add_study4_method_roles() %>%
    add_study4_analysis_eligibility()

  results <- results %>%
    add_stage1_estimates(
      fit_obj = fit_y,
      data = sim$lv1,
      cluster_var = "cid",
      within_var = "x",
      stage1_scores = stage2_df,
      true_slope_col = "true_u1"
    )

  dplyr::bind_cols(
    results,
    stage1_diag[rep(1L, nrow(results)), , drop = FALSE],
    measurement_diag[rep(1L, nrow(results)), , drop = FALSE]
  ) %>%
    add_study_result_context(condition, truth)
}
