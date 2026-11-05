#' Study-specific condition grids for Vig-Hallquist (2026)

fixed_params <- list(
  x_mean = 0,
  x_variance = 1,
  z_variance = 1,
  tau0 = 0.9,
  beta0z = 1.5,
  beta1z = 0.4,
  theta0_standardized = 0.4,
  gamma0_outcome = 1.0,
  gamma1_outcome = 0.6,
  gamma0_predictor = 0.0,
  gamma1_predictor = 0.5,
  gamma0_process_y = 0.0,
  gamma1_process_y = 0.5,
  gamma0_process_q = 0.5,
  gamma1_process_q = 0.4
)

cluster_balance_modes <- function() {
  c("balanced", "unbalanced", "informative_unbalanced")
}

balance_mode_to_sim_arg <- function(balance_mode) {
  switch(
    as.character(balance_mode[[1]]),
    balanced = "balanced",
    unbalanced = FALSE,
    informative_unbalanced = "highly_unbalanced",
    stop("Unsupported balance mode: ", balance_mode)
  )
}

residual_structures <- function() {
  tibble::tribble(
    ~r_structure, ~r_rho,
    "iid", NA_real_,
    "ar1", 0.3,
    "ar1", 0.6
  )
}

#' Add fixed posterior-reliability calibration parameters to a condition grid.
#'
#' Calibration is performed once while the grid is constructed. The resulting
#' population parameters are stored in the manifest and reused unchanged across
#' Monte Carlo replications.
#' 
#' NOTE: BLUP as outcome study called cluster size `n_trial` so there is some
#' mismatch in this function in terms of naming to be aware of.
calibrate_reliability_design <- function(
    design,
    tau0 = 0.9,
    calibration_reference_n = 1001L) {
  required <- c(
    "target_reliability", "standardized_beta_target", "marginal_rho",
    "mean_clus_size", "sigma", "balance_mode", "min_clus_size",
    "highly_unbalanced_min_clus_size", "highly_unbalanced_power",
    "r_structure", "r_rho"
  )
  missing <- setdiff(required, names(design))
  if (length(missing) > 0L) {
    stop("Reliability design is missing columns: ", paste(missing, collapse = ", "))
  }
  if (!exists("calibrate_random_slope_condition", mode = "function")) {
    stop(
      "`calibrate_random_slope_condition()` is unavailable. ",
      "Source `R/reliability_calibration.R` before building this grid."
    )
  }
  
  calibrated_rows <- lapply(seq_len(nrow(design)), function(i) {
    condition <- design[i, , drop = FALSE]
    calibrated <- calibrate_random_slope_condition(
      target_reliability = condition$target_reliability[[1]],
      standardized_beta = condition$standardized_beta_target[[1]],
      marginal_rho = condition$marginal_rho[[1]],
      tau0 = tau0,
      mean_n_trial = condition$mean_clus_size[[1]],
      sigma = condition$sigma[[1]],
      balance_mode = condition$balance_mode[[1]],
      min_n_trial = condition$min_clus_size[[1]],
      highly_unbalanced_min_n_trial = condition$highly_unbalanced_min_clus_size[[1]],
      highly_unbalanced_power = condition$highly_unbalanced_power[[1]],
      r_spec = condition_to_r_spec(condition),
      n_reference = calibration_reference_n
    )
    
    data.frame(
      calibration_tau0 = tau0,
      achieved_reliability = calibrated$achieved_reliability,
      standardized_beta = calibrated$standardized_beta,
      structural_r2 = calibrated$structural_r_squared,
      gamma_x_on_slope = calibrated$gamma_x_on_slope,
      slope_variance_marginal = calibrated$slope_variance_marginal,
      slope_variance_residual = calibrated$slope_variance_residual,
      tau1_residual = calibrated$tau1_residual,
      rho_residual = calibrated$rho_residual,
      reference_mean_clus_size = calibrated$reference_mean_n_trial,
      reference_min_clus_size = calibrated$reference_min_n_trial,
      reference_max_clus_size = calibrated$reference_max_n_trial,
      calibration_reference_n = calibrated$calibration_reference_n
    )
  })
  
  calibrated <- dplyr::bind_rows(calibrated_rows)
  dplyr::bind_cols(design, calibrated) %>%
    # `rho` remains the planned marginal correlation for reporting. The runner
    # explicitly passes `rho_residual` to simulate_dataset().
    dplyr::mutate(
      rho = .data$marginal_rho,
      tau1 = .data$tau1_residual
    )
}

# this is a smoke version of study 1 for testing
make_study0_design <- function() {
  tidyr::tibble(
    study = "study0",
    num_clus = 30L,
    mean_clus_size = 5L,
    target_reliability = 0.5,
    marginal_rho = 0.0,
    standardized_beta_target = 0.6,
    balance_mode = "balanced",
    min_clus_size = 2L,
    highly_unbalanced_min_clus_size = 2L,
    highly_unbalanced_power = 3,
    r_structure = "iid",
    r_rho = NA_real_,
    sigma = 1.0,
    study_label = "BLUP as Outcome (Smoke)",
    study_structure = "w"
  ) %>%
    calibrate_reliability_design(
      tau0 = fixed_params$tau0,
      calibration_reference_n = 1001L
    ) %>%
    rename(
      beta1w = gamma_x_on_slope # level 2 structural slope is now "w"
    )
}

make_study1_design <- function() {
  tidyr::crossing(
    study = "study1",
    num_clus = c(30L, 50L, 100L, 150L, 300L),
    mean_clus_size = c(3L, 5L, 10L, 25L),
    target_reliability = c(0.25, 0.5, 0.8),
    marginal_rho = c(-0.50, 0, 0.50),
    standardized_beta_target = c(0, 0.2, 0.4, 0.6),
    balance_mode = "balanced",
    min_clus_size = 2L,
    highly_unbalanced_min_clus_size = 2L,
    highly_unbalanced_power = 3,
    r_structure = "iid",
    r_rho = NA_real_,
    sigma = 1.0,
    study_label = "BLUP as Outcome",
    study_structure = "w"
    ) %>%
    calibrate_reliability_design(
      tau0 = fixed_params$tau0,
      calibration_reference_n = 1001L
    ) %>%
    rename(
      beta1w = gamma_x_on_slope # level 2 structural slope is now "w"
    )
}

make_study2_design <- function() {
  first_stage_grid <- tidyr::crossing(
    mean_clus_size = c(3L, 5L, 10L, 25L),
    target_reliability = c(0.25, 0.5, 0.8),
    marginal_rho = c(-0.50, 0, 0.50),
    balance_mode = "balanced",
    min_clus_size = 2L,
    highly_unbalanced_min_clus_size = 2L,
    highly_unbalanced_power = 3,
    r_structure = "iid",
    r_rho = NA_real_,
    sigma = 1.0
  )

  calibrated_first_stage <- lapply(seq_len(nrow(first_stage_grid)), function(i) {
    condition <- first_stage_grid[i, , drop = FALSE]
    reference <- make_reliability_reference_design(
      mean_n_trial = condition$mean_clus_size[[1]],
      sigma = condition$sigma[[1]],
      balance_mode = condition$balance_mode[[1]],
      min_n_trial = condition$min_clus_size[[1]],
      highly_unbalanced_min_n_trial = condition$highly_unbalanced_min_clus_size[[1]],
      highly_unbalanced_power = condition$highly_unbalanced_power[[1]],
      r_spec = condition_to_r_spec(condition),
      n_reference = 1001L
    )
    reliability <- calibrate_slope_variance(
      target_reliability = condition$target_reliability[[1]],
      Z_list = reference$Z_list,
      R_list = reference$R_list,
      weights = reference$count_weights,
      intercept_variance = fixed_params$tau0^2,
      intercept_slope_correlation = condition$marginal_rho[[1]]
    )
    dplyr::bind_cols(
      condition,
      tibble::tibble(
        calibration_tau0 = fixed_params$tau0,
        achieved_reliability = reliability$achieved_reliability,
        slope_variance_marginal = reliability$slope_variance_marginal,
        tau1 = sqrt(reliability$slope_variance_marginal),
        rho = condition$marginal_rho[[1]],
        G_marginal = list(reliability$G_marginal),
        reference_mean_clus_size = reference$mean_trial_count,
        reference_min_clus_size = reference$min_trial_count,
        reference_max_clus_size = reference$max_trial_count,
        calibration_reference_n = length(reference$trial_counts)
      )
    )
  }) %>%
    dplyr::bind_rows()

  design <- tidyr::crossing(
    calibrated_first_stage,
    study = "study2",
    num_clus = c(30L, 50L, 100L, 150L, 300L),
    standardized_beta_target = c(0, 0.2, 0.4, 0.6),
    structural_target = "intercept_slope", # only one structural target for now
    study_label = "BLUP as Predictor",
    study_structure = "z"
  )

  structural_rows <- lapply(seq_len(nrow(design)), function(i) {
    condition <- design[i, , drop = FALSE]
    structural <- calibrate_blup_predictor_effect(
      G_marginal = condition$G_marginal[[1]],
      standardized_slope_beta = condition$standardized_beta_target[[1]],
      structural_target = condition$structural_target[[1]],
      nuisance_intercept_beta = fixed_params$beta1z,
      outcome_variance = fixed_params$z_variance
    )
    tibble::tibble(
      standardized_beta = structural$standardized_slope_beta,
      beta1z = structural$beta1_intercept,
      beta2z = structural$beta2_slope,
      structural_r2 = structural$total_structural_r_squared,
      focal_unique_r2 = structural$focal_unique_r_squared,
      outcome_residual_variance = structural$outcome_residual_variance
    )
  }) %>%
    dplyr::bind_rows()

  dplyr::bind_cols(
    dplyr::select(design, -G_marginal),
    structural_rows
  )
}

make_study3_design <- function() {
  first_stage_grid <- tidyr::crossing(
    mean_clus_size = c(3L, 5L, 10L),
    target_reliability = c(0.25, 0.8),
    marginal_rho = c(0, 0.5),
    balance_mode = "balanced",
    min_clus_size = 2L,
    highly_unbalanced_min_clus_size = 2L,
    highly_unbalanced_power = 3,
    r_structure = "iid",
    r_rho = NA_real_,
    sigma = 1.0
  )

  calibrated_first_stage <- lapply(seq_len(nrow(first_stage_grid)), function(i) {
    condition <- first_stage_grid[i, , drop = FALSE]
    reference <- make_reliability_reference_design(
      mean_n_trial = condition$mean_clus_size[[1]],
      sigma = condition$sigma[[1]],
      balance_mode = condition$balance_mode[[1]],
      min_n_trial = condition$min_clus_size[[1]],
      highly_unbalanced_min_n_trial =
        condition$highly_unbalanced_min_clus_size[[1]],
      highly_unbalanced_power = condition$highly_unbalanced_power[[1]],
      r_spec = condition_to_r_spec(condition),
      n_reference = 1001L
    )
    reliability <- calibrate_slope_variance(
      target_reliability = condition$target_reliability[[1]],
      Z_list = reference$Z_list,
      R_list = reference$R_list,
      weights = reference$count_weights,
      intercept_variance = fixed_params$tau0^2,
      intercept_slope_correlation = condition$marginal_rho[[1]]
    )
    dplyr::bind_cols(
      condition,
      tibble::tibble(
        achieved_reliability = reliability$achieved_reliability,
        slope_variance_marginal = reliability$slope_variance_marginal,
        tau1 = sqrt(reliability$slope_variance_marginal),
        G_marginal = list(reliability$G_marginal),
        reference_mean_clus_size = reference$mean_trial_count,
        reference_min_clus_size = reference$min_trial_count,
        reference_max_clus_size = reference$max_trial_count,
        calibration_reference_n = length(reference$trial_counts)
      )
    )
  }) %>%
    dplyr::bind_rows()

  design <- tidyr::crossing(
    study = "study3",
    num_clus = c(50L, 100L, 150L, 300L),
    mean_clus_size_y = c(3L, 5L, 10L),
    mean_clus_size_q = c(3L, 5L, 10L),
    target_reliability_y = c(0.25, 0.8),
    target_reliability_q = c(0.25, 0.8),
    marginal_rho = c(0, 0.5),
    standardized_beta_target = c(0, 0.2, 0.5),
    structural_target = "intercept_slope", # only one structural target for now
    study_label = "BLUP as Predictor and Outcome",
    study_structure = "dual_process"
  )

  y_calibration <- calibrated_first_stage %>%
    dplyr::transmute(
      mean_clus_size_y = mean_clus_size,
      target_reliability_y = target_reliability,
      marginal_rho,
      achieved_reliability_y = achieved_reliability,
      slope_variance_marginal_y = slope_variance_marginal,
      tau1_y = tau1,
      G_y_marginal = G_marginal,
      reference_mean_clus_size_y = reference_mean_clus_size,
      reference_min_clus_size_y = reference_min_clus_size,
      reference_max_clus_size_y = reference_max_clus_size
    )
  q_calibration <- calibrated_first_stage %>%
    dplyr::transmute(
      mean_clus_size_q = mean_clus_size,
      target_reliability_q = target_reliability,
      marginal_rho,
      achieved_reliability_q = achieved_reliability,
      slope_variance_marginal_q = slope_variance_marginal,
      tau1_q = tau1,
      G_q_marginal = G_marginal,
      reference_mean_clus_size_q = reference_mean_clus_size,
      reference_min_clus_size_q = reference_min_clus_size,
      reference_max_clus_size_q = reference_max_clus_size
    )

  design <- design %>%
    dplyr::left_join(
      y_calibration,
      by = c("mean_clus_size_y", "target_reliability_y", "marginal_rho")
    ) %>%
    dplyr::left_join(
      q_calibration,
      by = c("mean_clus_size_q", "target_reliability_q", "marginal_rho")
    )

  structural_rows <- lapply(seq_len(nrow(design)), function(i) {
    condition <- design[i, , drop = FALSE]
    structural <- calibrate_dual_process_effect(
      G_predictor = condition$G_y_marginal[[1]],
      G_outcome = condition$G_q_marginal[[1]],
      standardized_slope_beta = condition$standardized_beta_target[[1]],
      structural_target = condition$structural_target[[1]],
      nuisance_intercept_standardized_beta =
        fixed_params$theta0_standardized
    )
    tibble::tibble(
      standardized_beta = structural$standardized_slope_beta,
      standardized_theta0 = structural$standardized_intercept_beta,
      theta0 = structural$theta0_intercept,
      theta1 = structural$theta1_slope,
      structural_r2 = structural$total_structural_r_squared,
      focal_unique_r2 = structural$focal_unique_r_squared,
      slope_variance_residual_q =
        structural$outcome_slope_residual_variance,
      tau1_residual_q =
        sqrt(structural$outcome_slope_residual_variance),
      rho_residual_q = structural$outcome_residual_correlation
    )
  }) %>%
    dplyr::bind_rows()

  dplyr::bind_cols(
    dplyr::select(design, -G_y_marginal, -G_q_marginal),
    structural_rows
  ) %>%
    dplyr::mutate(
      calibration_tau0 = fixed_params$tau0,
      sigma_y = 1.0,
      sigma_q = 1.0,
      balance_mode = "balanced",
      min_clus_size = 2L,
      highly_unbalanced_min_clus_size = 2L,
      highly_unbalanced_power = 3,
      r_structure = "iid",
      r_rho = NA_real_
    )
}

#' Return the cluster-information specification for one Study 4 profile.
#'
#' Primary profiles all have mean cluster size 10. The information-matched
#' profile is a focused falsification control: its cluster sizes vary, but its
#' slope columns are rescaled during calibration and generation so that
#' `sum(x^2) = 9` for every cluster.
study4_profile_spec <- function(information_profile) {
  information_profile <- as.character(information_profile[[1]])
  switch(
    information_profile,
    homogeneous = list(
      sizes = 10L,
      weights = 1,
      information_matched = FALSE
    ),
    moderate = list(
      sizes = c(5L, 15L),
      weights = c(0.5, 0.5),
      information_matched = FALSE
    ),
    severe = list(
      sizes = c(3L, 17L),
      weights = c(0.5, 0.5),
      information_matched = FALSE
    ),
    severe_information_matched = list(
      sizes = c(3L, 17L),
      weights = c(0.5, 0.5),
      information_matched = TRUE
    ),
    stop("Unsupported Study 4 information profile: ", information_profile)
  )
}

#' Construct the random-effect design for one Study 4 cluster type.
study4_time_design <- function(cluster_size, information_matched = FALSE) {
  Z <- make_reliability_time_design(cluster_size)
  if (isTRUE(information_matched)) {
    Z[, "slope"] <- Z[, "slope"] * sqrt(9 / (cluster_size - 1))
  }
  Z
}

#' Compute weighted quantiles for a discrete Study 4 profile.
study4_weighted_quantile <- function(x, weights, probs) {
  keep <- is.finite(x) & is.finite(weights) & weights > 0
  x <- as.numeric(x[keep])
  weights <- as.numeric(weights[keep])
  if (length(x) == 0L) {
    return(rep(NA_real_, length(probs)))
  }
  ord <- order(x)
  x <- x[ord]
  cumulative_weight <- cumsum(weights[ord] / sum(weights))
  vapply(
    probs,
    function(prob) {
      x[[which(cumulative_weight >= min(1, max(0, prob)))[[1L]]]]
    },
    numeric(1L)
  )
}

#' Weighted RMS Frobenius distance of matrices from their elementwise mean.
study4_matrix_rms_dispersion <- function(matrices, weights = NULL) {
  if (is.null(weights)) {
    weights <- rep(1, length(matrices))
  }
  valid <- vapply(
    matrices,
    function(mat) is.matrix(mat) && all(dim(mat) == c(2L, 2L)) && all(is.finite(mat)),
    logical(1L)
  ) & is.finite(weights) & weights > 0
  if (!any(valid)) {
    return(NA_real_)
  }
  matrices <- matrices[valid]
  weights <- weights[valid] / sum(weights[valid])
  mean_matrix <- Reduce(
    `+`,
    Map(function(mat, weight) mat * weight, matrices, weights)
  )
  sqrt(sum(vapply(
    seq_along(matrices),
    function(i) weights[[i]] * sum((matrices[[i]] - mean_matrix)^2),
    numeric(1L)
  )))
}

#' Summarize population or fitted Study 4 measurement matrices.
study4_measurement_matrix_summary <- function(
    lambda_matrices,
    theta_matrices,
    score_error_covariances,
    weights = NULL) {
  if (is.null(weights)) {
    weights <- rep(1, length(lambda_matrices))
  }
  weights <- weights / sum(weights)
  lambda22 <- vapply(lambda_matrices, function(mat) mat[2L, 2L], numeric(1L))
  theta22 <- vapply(theta_matrices, function(mat) mat[2L, 2L], numeric(1L))
  ols_var22 <- vapply(
    score_error_covariances,
    function(mat) mat[2L, 2L],
    numeric(1L)
  )

  list(
    lambda22_mean = stats::weighted.mean(lambda22, weights),
    lambda22_min = min(lambda22),
    lambda22_max = max(lambda22),
    lambda_matrix_frobenius_rms_dispersion =
      study4_matrix_rms_dispersion(lambda_matrices, weights),
    theta22_mean = stats::weighted.mean(theta22, weights),
    theta22_min = min(theta22),
    theta22_max = max(theta22),
    theta_matrix_frobenius_rms_dispersion =
      study4_matrix_rms_dispersion(theta_matrices, weights),
    ols_var22_mean = stats::weighted.mean(ols_var22, weights),
    ols_var22_min = min(ols_var22),
    ols_var22_max = max(ols_var22),
    ols_cov_matrix_frobenius_rms_dispersion =
      study4_matrix_rms_dispersion(score_error_covariances, weights)
  )
}

#' Calibrate a Study 4 first-stage information profile.
calibrate_study4_profile <- function(
    information_profile,
    target_reliability,
    marginal_rho,
    sigma = 1.0) {
  spec <- study4_profile_spec(information_profile)
  Z_list <- lapply(
    spec$sizes,
    study4_time_design,
    information_matched = spec$information_matched
  )
  R_list <- lapply(
    spec$sizes,
    make_reliability_residual_covariance,
    sigma = sigma,
    r_spec = list(structure = "iid")
  )
  calibration <- calibrate_slope_variance(
    target_reliability = target_reliability,
    Z_list = Z_list,
    R_list = R_list,
    weights = spec$weights,
    intercept_variance = fixed_params$tau0^2,
    intercept_slope_correlation = marginal_rho
  )
  type_reliability <- vapply(
    calibration$posterior_covariances,
    function(V_i) {
      1 - V_i[2L, 2L] / calibration$G_marginal[2L, 2L]
    },
    numeric(1L)
  )
  weighted_mean <- stats::weighted.mean(type_reliability, spec$weights)
  weighted_variance <- stats::weighted.mean(
    (type_reliability - weighted_mean)^2,
    spec$weights
  )
  reliability_quartiles <- study4_weighted_quantile(
    type_reliability,
    spec$weights,
    c(0.25, 0.75)
  )
  population_measurement <- Map(function(Z_i, R_i) {
    sigma_i <- Z_i %*% calibration$G_marginal %*% t(Z_i) + R_i
    A_i <- t(solve(sigma_i, Z_i %*% calibration$G_marginal))
    R_inv_Z <- solve(R_i, Z_i)
    list(
      lambda = A_i %*% Z_i,
      theta = A_i %*% R_i %*% t(A_i),
      score_error_covariance = solve(crossprod(Z_i, R_inv_Z))
    )
  }, Z_list, R_list)
  measurement_diagnostics <- study4_measurement_matrix_summary(
    lambda_matrices = lapply(population_measurement, `[[`, "lambda"),
    theta_matrices = lapply(population_measurement, `[[`, "theta"),
    score_error_covariances = lapply(
      population_measurement,
      `[[`,
      "score_error_covariance"
    ),
    weights = spec$weights
  )

  list(
    calibration = calibration,
    profile = spec,
    type_reliability = type_reliability,
    reliability_mean = weighted_mean,
    reliability_sd = sqrt(weighted_variance),
    reliability_iqr = diff(reliability_quartiles),
    reliability_min = min(type_reliability),
    reliability_max = max(type_reliability),
    measurement_diagnostics = measurement_diagnostics
  )
}

make_study4_design <- function() {
  primary <- tidyr::crossing(
    study = "study4",
    num_clus = c(50L, 150L, 300L),
    target_reliability = c(0.25, 0.5, 0.8),
    information_profile = c("homogeneous", "moderate", "severe"),
    marginal_rho = c(0, 0.5),
    standardized_beta_target = c(0, 0.4),
    structural_target = "intercept_slope",
    is_falsification_control = FALSE
  )
  information_matched_control <- tidyr::crossing(
    study = "study4",
    num_clus = c(50L, 150L, 300L),
    target_reliability = 0.5,
    information_profile = "severe_information_matched",
    marginal_rho = 0,
    standardized_beta_target = 0.4,
    structural_target = "intercept_slope",
    is_falsification_control = TRUE
  )
  design <- dplyr::bind_rows(primary, information_matched_control)

  calibration_rows <- lapply(seq_len(nrow(design)), function(i) {
    condition <- design[i, , drop = FALSE]
    calibrated <- calibrate_study4_profile(
      information_profile = condition$information_profile[[1]],
      target_reliability = condition$target_reliability[[1]],
      marginal_rho = condition$marginal_rho[[1]],
      sigma = 1.0
    )
    cal <- calibrated$calibration
    spec <- calibrated$profile
    reliability <- calibrated$type_reliability
    measurement <- calibrated$measurement_diagnostics
    structural <- calibrate_blup_predictor_effect(
      G_marginal = cal$G_marginal,
      standardized_slope_beta = condition$standardized_beta_target[[1]],
      structural_target = condition$structural_target[[1]],
      nuisance_intercept_beta = fixed_params$beta1z,
      outcome_variance = fixed_params$z_variance
    )

    tibble::tibble(
      mean_clus_size = stats::weighted.mean(spec$sizes, spec$weights),
      profile_min_clus_size = min(spec$sizes),
      profile_max_clus_size = max(spec$sizes),
      profile_small_clus_size = min(spec$sizes),
      profile_large_clus_size = max(spec$sizes),
      profile_small_weight = spec$weights[[1L]],
      profile_large_weight = if (length(spec$weights) > 1L) {
        spec$weights[[length(spec$weights)]]
      } else {
        0
      },
      information_matched = spec$information_matched,
      calibration_tau0 = fixed_params$tau0,
      achieved_reliability = cal$achieved_reliability,
      reliability_sd = calibrated$reliability_sd,
      reliability_iqr = calibrated$reliability_iqr,
      reliability_min = calibrated$reliability_min,
      reliability_max = calibrated$reliability_max,
      reliability_small = reliability[[1L]],
      reliability_large = reliability[[length(reliability)]],
      population_lambda22_mean = measurement$lambda22_mean,
      population_lambda22_min = measurement$lambda22_min,
      population_lambda22_max = measurement$lambda22_max,
      population_lambda_matrix_frobenius_rms_dispersion =
        measurement$lambda_matrix_frobenius_rms_dispersion,
      population_theta22_mean = measurement$theta22_mean,
      population_theta22_min = measurement$theta22_min,
      population_theta22_max = measurement$theta22_max,
      population_theta_matrix_frobenius_rms_dispersion =
        measurement$theta_matrix_frobenius_rms_dispersion,
      population_ols_var22_mean = measurement$ols_var22_mean,
      population_ols_var22_min = measurement$ols_var22_min,
      population_ols_var22_max = measurement$ols_var22_max,
      population_ols_cov_matrix_frobenius_rms_dispersion =
        measurement$ols_cov_matrix_frobenius_rms_dispersion,
      slope_variance_marginal = cal$slope_variance_marginal,
      tau1 = sqrt(cal$slope_variance_marginal),
      rho = condition$marginal_rho[[1]],
      standardized_beta = structural$standardized_slope_beta,
      beta1z = structural$beta1_intercept,
      beta2z = structural$beta2_slope,
      structural_r2 = structural$total_structural_r_squared,
      focal_unique_r2 = structural$focal_unique_r_squared,
      outcome_residual_variance = structural$outcome_residual_variance
    )
  }) %>%
    dplyr::bind_rows()

  dplyr::bind_cols(design, calibration_rows) %>%
    dplyr::mutate(
      balance_mode = "heterogeneous_profile",
      min_clus_size = .data$profile_min_clus_size,
      highly_unbalanced_min_clus_size = .data$profile_min_clus_size,
      highly_unbalanced_power = NA_real_,
      r_structure = "iid",
      r_rho = NA_real_,
      sigma = 1.0,
      study_label = "Heterogeneous Cluster Information",
      study_structure = "z"
    )
}

# -------------------------------------------------------------------------
# Amended (v2) fixed-covariance-shape designs
# -------------------------------------------------------------------------

#' Calibrate a standard Study 1--3 first-stage grid at fixed covariance shape.
#'
#' This helper is deliberately measurement-only. Structural coefficients are
#' added by the study-specific wrappers after the target marginal covariance
#' and Level-1 residual scale have been fixed.
calibrate_shape_preserving_first_stage_grid <- function(
    first_stage_grid,
    slope_intercept_variance_ratio = 1,
    calibration_target = "marginal_slope",
    calibration_reference_n = 1001L) {
  required <- c(
    "mean_clus_size", "target_reliability", "marginal_rho",
    "balance_mode", "min_clus_size", "highly_unbalanced_min_clus_size",
    "highly_unbalanced_power", "r_structure", "r_rho"
  )
  missing <- setdiff(required, names(first_stage_grid))
  if (length(missing) > 0L) {
    stop(
      "Shape-preserving first-stage grid is missing columns: ",
      paste(missing, collapse = ", ")
    )
  }

  lapply(seq_len(nrow(first_stage_grid)), function(i) {
    condition <- first_stage_grid[i, , drop = FALSE]
    reference <- make_reliability_reference_design(
      mean_n_trial = condition$mean_clus_size[[1]],
      sigma = 1,
      balance_mode = condition$balance_mode[[1]],
      min_n_trial = condition$min_clus_size[[1]],
      highly_unbalanced_min_n_trial =
        condition$highly_unbalanced_min_clus_size[[1]],
      highly_unbalanced_power = condition$highly_unbalanced_power[[1]],
      r_spec = condition_to_r_spec(condition),
      n_reference = calibration_reference_n
    )
    calibration <- calibrate_shape_preserving_measurement(
      target_value = condition$target_reliability[[1]],
      Z_list = reference$Z_list,
      R_shape_list = reference$R_list,
      weights = reference$count_weights,
      intercept_variance = fixed_params$tau0^2,
      slope_intercept_variance_ratio = slope_intercept_variance_ratio,
      intercept_slope_correlation = condition$marginal_rho[[1]],
      calibration_target = calibration_target
    )

    dplyr::bind_cols(
      condition,
      tibble::tibble(
        study_version = "v2",
        calibration_version = calibration$calibration_version,
        calibration_arm = paste0(
          "shape_preserving_", calibration$calibration_target
        ),
        calibration_metric = calibration$calibration_target,
        calibration_target = calibration$calibration_target,
        calibration_target_value = calibration$calibration_target_value,
        achieved_calibration_value =
          calibration$achieved_calibration_value,
        covariance_shape_fixed = TRUE,
        calibration_tau0 = fixed_params$tau0,
        calibration_tau0_sq = calibration$intercept_variance,
        calibration_tau1_sq = calibration$slope_variance,
        achieved_reliability =
          calibration$achieved_marginal_slope_reliability,
        achieved_partial_reliability =
          calibration$achieved_residualized_slope_reliability,
        intercept_icc = calibration$intercept_icc,
        marginal_slope_icc = calibration$marginal_slope_icc,
        conditional_slope_icc = calibration$conditional_slope_icc,
        slope_variance_marginal = calibration$slope_variance,
        residualized_slope_variance =
          calibration$residualized_slope_variance,
        slope_intercept_variance_ratio =
          calibration$slope_intercept_variance_ratio,
        G_condition_number = calibration$G_condition_number,
        sigma = calibration$sigma,
        tau1_marginal = sqrt(calibration$slope_variance),
        rho = condition$marginal_rho[[1]],
        G_marginal = list(calibration$G_marginal),
        reference_mean_clus_size = reference$mean_trial_count,
        reference_min_clus_size = reference$min_trial_count,
        reference_max_clus_size = reference$max_trial_count,
        calibration_reference_n = length(reference$trial_counts)
      )
    )
  }) %>%
    dplyr::bind_rows()
}

#' Build amended Study 1 (fixed covariance shape; marginal reliability target).
make_study1_v2_design <- function(
    slope_intercept_variance_ratio = 1) {
  first_stage <- tidyr::crossing(
    mean_clus_size = c(3L, 5L, 10L, 25L),
    target_reliability = c(0.25, 0.5, 0.8),
    marginal_rho = c(-0.50, 0, 0.50),
    balance_mode = "balanced",
    min_clus_size = 2L,
    highly_unbalanced_min_clus_size = 2L,
    highly_unbalanced_power = 3,
    r_structure = "iid",
    r_rho = NA_real_
  ) %>%
    calibrate_shape_preserving_first_stage_grid(
      slope_intercept_variance_ratio = slope_intercept_variance_ratio
    )

  design <- tidyr::crossing(
    first_stage,
    study = "study1v2",
    num_clus = c(30L, 50L, 100L, 150L, 300L),
    standardized_beta_target = c(0, 0.2, 0.4, 0.6),
    study_label = "BLUP as Outcome (v2 fixed shape)",
    study_structure = "w"
  )
  structural_rows <- lapply(seq_len(nrow(design)), function(i) {
    condition <- design[i, , drop = FALSE]
    structural <- decompose_structural_slope(
      list(G_marginal = condition$G_marginal[[1]]),
      standardized_beta = condition$standardized_beta_target[[1]]
    )
    tibble::tibble(
      standardized_beta = structural$standardized_beta,
      structural_r2 = structural$structural_r_squared,
      beta1w = structural$gamma_x_on_slope,
      slope_variance_residual = structural$slope_variance_residual,
      tau1 = structural$slope_sd_residual,
      rho_residual = structural$intercept_slope_correlation_residual,
      structural_residual_G_min_eigen = min(eigen(
        structural$G_residual,
        symmetric = TRUE,
        only.values = TRUE
      )$values)
    )
  }) %>%
    dplyr::bind_rows()

  dplyr::bind_cols(
    dplyr::select(design, -G_marginal, -tau1_marginal),
    structural_rows
  )
}

#' Build amended Study 2 (fixed covariance shape; marginal reliability target).
make_study2_v2_design <- function(
    slope_intercept_variance_ratio = 1) {
  first_stage <- tidyr::crossing(
    mean_clus_size = c(3L, 5L, 10L, 25L),
    target_reliability = c(0.25, 0.5, 0.8),
    marginal_rho = c(-0.50, 0, 0.50),
    balance_mode = "balanced",
    min_clus_size = 2L,
    highly_unbalanced_min_clus_size = 2L,
    highly_unbalanced_power = 3,
    r_structure = "iid",
    r_rho = NA_real_
  ) %>%
    calibrate_shape_preserving_first_stage_grid(
      slope_intercept_variance_ratio = slope_intercept_variance_ratio
    )

  design <- tidyr::crossing(
    first_stage,
    study = "study2v2",
    num_clus = c(30L, 50L, 100L, 150L, 300L),
    standardized_beta_target = c(0, 0.2, 0.4, 0.6),
    structural_target = "intercept_slope",
    study_label = "BLUP as Predictor (v2 fixed shape)",
    study_structure = "z"
  )
  structural_rows <- lapply(seq_len(nrow(design)), function(i) {
    condition <- design[i, , drop = FALSE]
    structural <- calibrate_blup_predictor_effect(
      G_marginal = condition$G_marginal[[1]],
      standardized_slope_beta = condition$standardized_beta_target[[1]],
      structural_target = condition$structural_target[[1]],
      nuisance_intercept_beta = fixed_params$beta1z,
      outcome_variance = fixed_params$z_variance
    )
    tibble::tibble(
      standardized_beta = structural$standardized_slope_beta,
      beta1z = structural$beta1_intercept,
      beta2z = structural$beta2_slope,
      structural_r2 = structural$total_structural_r_squared,
      focal_unique_r2 = structural$focal_unique_r_squared,
      outcome_residual_variance = structural$outcome_residual_variance
    )
  }) %>%
    dplyr::bind_rows()

  dplyr::bind_cols(
    dplyr::select(design, -G_marginal),
    structural_rows
  ) %>%
    dplyr::mutate(tau1 = .data$tau1_marginal) %>%
    dplyr::select(-tau1_marginal)
}

#' Build amended Study 3 using the common calibration independently twice.
#'
#' The reliability and effect arguments retain the published v2 grid as their
#' defaults.  They are exposed so focused sensitivity modules can reuse the
#' exact Study-3 DGM without adding new levels to the primary factorial design.
make_study3_v2_design <- function(
    slope_intercept_variance_ratio = 1,
    target_reliabilities = c(0.25, 0.8),
    standardized_beta_targets = c(0, 0.2, 0.5)) {
  first_stage <- tidyr::crossing(
    mean_clus_size = c(3L, 5L, 10L),
    target_reliability = target_reliabilities,
    marginal_rho = c(0, 0.5),
    balance_mode = "balanced",
    min_clus_size = 2L,
    highly_unbalanced_min_clus_size = 2L,
    highly_unbalanced_power = 3,
    r_structure = "iid",
    r_rho = NA_real_
  ) %>%
    calibrate_shape_preserving_first_stage_grid(
      slope_intercept_variance_ratio = slope_intercept_variance_ratio
    )

  y_calibration <- first_stage %>%
    dplyr::transmute(
      mean_clus_size_y = mean_clus_size,
      target_reliability_y = target_reliability,
      marginal_rho,
      sigma_y = sigma,
      achieved_reliability_y = achieved_reliability,
      achieved_partial_reliability_y = achieved_partial_reliability,
      intercept_icc_y = intercept_icc,
      marginal_slope_icc_y = marginal_slope_icc,
      conditional_slope_icc_y = conditional_slope_icc,
      G_condition_number_y = G_condition_number,
      slope_variance_marginal_y = slope_variance_marginal,
      residualized_slope_variance_y = residualized_slope_variance,
      tau1_y = tau1_marginal,
      G_y_marginal = G_marginal,
      reference_mean_clus_size_y = reference_mean_clus_size,
      reference_min_clus_size_y = reference_min_clus_size,
      reference_max_clus_size_y = reference_max_clus_size
    )
  q_calibration <- first_stage %>%
    dplyr::transmute(
      mean_clus_size_q = mean_clus_size,
      target_reliability_q = target_reliability,
      marginal_rho,
      sigma_q = sigma,
      achieved_reliability_q = achieved_reliability,
      achieved_partial_reliability_q = achieved_partial_reliability,
      intercept_icc_q = intercept_icc,
      marginal_slope_icc_q = marginal_slope_icc,
      conditional_slope_icc_q = conditional_slope_icc,
      G_condition_number_q = G_condition_number,
      slope_variance_marginal_q = slope_variance_marginal,
      residualized_slope_variance_q = residualized_slope_variance,
      tau1_q = tau1_marginal,
      G_q_marginal = G_marginal,
      reference_mean_clus_size_q = reference_mean_clus_size,
      reference_min_clus_size_q = reference_min_clus_size,
      reference_max_clus_size_q = reference_max_clus_size
    )

  design <- tidyr::crossing(
    study = "study3v2",
    num_clus = c(50L, 100L, 150L, 300L),
    mean_clus_size_y = c(3L, 5L, 10L),
    mean_clus_size_q = c(3L, 5L, 10L),
    target_reliability_y = target_reliabilities,
    target_reliability_q = target_reliabilities,
    marginal_rho = c(0, 0.5),
    standardized_beta_target = standardized_beta_targets,
    structural_target = "intercept_slope",
    study_label = "BLUP as Predictor and Outcome (v2 fixed shape)",
    study_structure = "dual_process"
  ) %>%
    dplyr::left_join(
      y_calibration,
      by = c("mean_clus_size_y", "target_reliability_y", "marginal_rho")
    ) %>%
    dplyr::left_join(
      q_calibration,
      by = c("mean_clus_size_q", "target_reliability_q", "marginal_rho")
    )

  structural_rows <- lapply(seq_len(nrow(design)), function(i) {
    condition <- design[i, , drop = FALSE]
    structural <- calibrate_dual_process_effect(
      G_predictor = condition$G_y_marginal[[1]],
      G_outcome = condition$G_q_marginal[[1]],
      standardized_slope_beta = condition$standardized_beta_target[[1]],
      structural_target = condition$structural_target[[1]],
      nuisance_intercept_standardized_beta =
        fixed_params$theta0_standardized
    )
    tibble::tibble(
      standardized_beta = structural$standardized_slope_beta,
      standardized_theta0 = structural$standardized_intercept_beta,
      theta0 = structural$theta0_intercept,
      theta1 = structural$theta1_slope,
      structural_r2 = structural$total_structural_r_squared,
      focal_unique_r2 = structural$focal_unique_r_squared,
      slope_variance_residual_q =
        structural$outcome_slope_residual_variance,
      tau1_residual_q = sqrt(structural$outcome_slope_residual_variance),
      rho_residual_q = structural$outcome_residual_correlation,
      structural_residual_G_min_eigen = min(eigen(
        structural$G_outcome_residual,
        symmetric = TRUE,
        only.values = TRUE
      )$values),
      structural_joint_G_min_eigen = min(eigen(
        structural$G_joint_marginal,
        symmetric = TRUE,
        only.values = TRUE
      )$values)
    )
  }) %>%
    dplyr::bind_rows()

  dplyr::bind_cols(
    dplyr::select(design, -G_y_marginal, -G_q_marginal),
    structural_rows
  ) %>%
    dplyr::mutate(
      study_version = "v2",
      calibration_version = "shape_preserving_v2",
      calibration_arm = "shape_preserving_marginal",
      calibration_metric = "marginal_slope",
      calibration_target = "marginal_slope",
      covariance_shape_fixed = TRUE,
      calibration_tau0 = fixed_params$tau0,
      calibration_tau0_sq = fixed_params$tau0^2,
      calibration_tau1_sq = fixed_params$tau0^2 *
        slope_intercept_variance_ratio,
      slope_intercept_variance_ratio = slope_intercept_variance_ratio,
      G_condition_number = .data$G_condition_number_y,
      balance_mode = "balanced",
      min_clus_size = 2L,
      highly_unbalanced_min_clus_size = 2L,
      highly_unbalanced_power = 3,
      r_structure = "iid",
      r_rho = NA_real_
    )
}

#' Calibrate one amended Study 4 information profile at fixed G.
calibrate_study4_profile_v2 <- function(
    information_profile,
    target_reliability,
    marginal_rho,
    slope_intercept_variance_ratio = 1) {
  spec <- study4_profile_spec(information_profile)
  Z_list <- lapply(
    spec$sizes,
    study4_time_design,
    information_matched = spec$information_matched
  )
  R_shape_list <- lapply(spec$sizes, diag)
  calibration <- calibrate_shape_preserving_measurement(
    target_value = target_reliability,
    Z_list = Z_list,
    R_shape_list = R_shape_list,
    weights = spec$weights,
    intercept_variance = fixed_params$tau0^2,
    slope_intercept_variance_ratio = slope_intercept_variance_ratio,
    intercept_slope_correlation = marginal_rho,
    calibration_target = "marginal_slope"
  )
  type_reliability <- vapply(seq_along(Z_list), function(i) {
    V_i <- posterior_random_effect_covariance(
      calibration$G_marginal,
      Z_list[[i]],
      calibration$residual_covariances[[i]]
    )
    1 - V_i[2L, 2L] / calibration$G_marginal[2L, 2L]
  }, numeric(1L))
  type_partial_reliability <- vapply(seq_along(Z_list), function(i) {
    expected_residualized_slope_reliability(
      calibration$G_marginal,
      list(Z_list[[i]]),
      list(calibration$residual_covariances[[i]])
    )
  }, numeric(1L))
  weighted_mean <- stats::weighted.mean(type_reliability, spec$weights)
  weighted_variance <- stats::weighted.mean(
    (type_reliability - weighted_mean)^2,
    spec$weights
  )
  reliability_quartiles <- study4_weighted_quantile(
    type_reliability,
    spec$weights,
    c(0.25, 0.75)
  )
  population_measurement <- Map(function(Z_i, R_i) {
    sigma_i <- Z_i %*% calibration$G_marginal %*% t(Z_i) + R_i
    A_i <- t(solve(sigma_i, Z_i %*% calibration$G_marginal))
    R_inv_Z <- solve(R_i, Z_i)
    list(
      lambda = A_i %*% Z_i,
      theta = A_i %*% R_i %*% t(A_i),
      score_error_covariance = solve(crossprod(Z_i, R_inv_Z))
    )
  }, Z_list, calibration$residual_covariances)
  measurement_diagnostics <- study4_measurement_matrix_summary(
    lambda_matrices = lapply(population_measurement, `[[`, "lambda"),
    theta_matrices = lapply(population_measurement, `[[`, "theta"),
    score_error_covariances = lapply(
      population_measurement,
      `[[`,
      "score_error_covariance"
    ),
    weights = spec$weights
  )

  list(
    calibration = calibration,
    profile = spec,
    type_reliability = type_reliability,
    type_partial_reliability = type_partial_reliability,
    reliability_mean = weighted_mean,
    reliability_sd = sqrt(weighted_variance),
    reliability_iqr = diff(reliability_quartiles),
    reliability_min = min(type_reliability),
    reliability_max = max(type_reliability),
    partial_reliability_mean = stats::weighted.mean(
      type_partial_reliability,
      spec$weights
    ),
    partial_reliability_min = min(type_partial_reliability),
    partial_reliability_max = max(type_partial_reliability),
    measurement_diagnostics = measurement_diagnostics
  )
}

#' Build amended Study 4 (fixed covariance shape across information profiles).
make_study4_v2_design <- function(
    slope_intercept_variance_ratio = 1) {
  primary <- tidyr::crossing(
    study = "study4v2",
    num_clus = c(50L, 150L, 300L),
    target_reliability = c(0.25, 0.5, 0.8),
    information_profile = c("homogeneous", "moderate", "severe"),
    marginal_rho = c(0, 0.5),
    standardized_beta_target = c(0, 0.4),
    structural_target = "intercept_slope",
    is_falsification_control = FALSE
  )
  information_matched_control <- tidyr::crossing(
    study = "study4v2",
    num_clus = c(50L, 150L, 300L),
    target_reliability = 0.5,
    information_profile = "severe_information_matched",
    marginal_rho = 0,
    standardized_beta_target = 0.4,
    structural_target = "intercept_slope",
    is_falsification_control = TRUE
  )
  design <- dplyr::bind_rows(primary, information_matched_control)

  calibration_rows <- lapply(seq_len(nrow(design)), function(i) {
    condition <- design[i, , drop = FALSE]
    calibrated <- calibrate_study4_profile_v2(
      information_profile = condition$information_profile[[1]],
      target_reliability = condition$target_reliability[[1]],
      marginal_rho = condition$marginal_rho[[1]],
      slope_intercept_variance_ratio = slope_intercept_variance_ratio
    )
    cal <- calibrated$calibration
    spec <- calibrated$profile
    reliability <- calibrated$type_reliability
    partial_reliability <- calibrated$type_partial_reliability
    measurement <- calibrated$measurement_diagnostics
    structural <- calibrate_blup_predictor_effect(
      G_marginal = cal$G_marginal,
      standardized_slope_beta = condition$standardized_beta_target[[1]],
      structural_target = condition$structural_target[[1]],
      nuisance_intercept_beta = fixed_params$beta1z,
      outcome_variance = fixed_params$z_variance
    )

    tibble::tibble(
      study_version = "v2",
      calibration_version = cal$calibration_version,
      calibration_arm = "shape_preserving_marginal",
      calibration_metric = cal$calibration_target,
      calibration_target = cal$calibration_target,
      calibration_target_value = cal$calibration_target_value,
      achieved_calibration_value = cal$achieved_calibration_value,
      covariance_shape_fixed = TRUE,
      mean_clus_size = stats::weighted.mean(spec$sizes, spec$weights),
      profile_min_clus_size = min(spec$sizes),
      profile_max_clus_size = max(spec$sizes),
      profile_small_clus_size = min(spec$sizes),
      profile_large_clus_size = max(spec$sizes),
      profile_small_weight = spec$weights[[1L]],
      profile_large_weight = if (length(spec$weights) > 1L) {
        spec$weights[[length(spec$weights)]]
      } else {
        0
      },
      information_matched = spec$information_matched,
      calibration_tau0 = fixed_params$tau0,
      calibration_tau0_sq = cal$intercept_variance,
      calibration_tau1_sq = cal$slope_variance,
      achieved_reliability = cal$achieved_marginal_slope_reliability,
      achieved_partial_reliability =
        cal$achieved_residualized_slope_reliability,
      intercept_icc = cal$intercept_icc,
      marginal_slope_icc = cal$marginal_slope_icc,
      conditional_slope_icc = cal$conditional_slope_icc,
      reliability_sd = calibrated$reliability_sd,
      reliability_iqr = calibrated$reliability_iqr,
      reliability_min = calibrated$reliability_min,
      reliability_max = calibrated$reliability_max,
      reliability_small = reliability[[1L]],
      reliability_large = reliability[[length(reliability)]],
      partial_reliability_mean = calibrated$partial_reliability_mean,
      partial_reliability_min = calibrated$partial_reliability_min,
      partial_reliability_max = calibrated$partial_reliability_max,
      partial_reliability_small = partial_reliability[[1L]],
      partial_reliability_large =
        partial_reliability[[length(partial_reliability)]],
      population_lambda22_mean = measurement$lambda22_mean,
      population_lambda22_min = measurement$lambda22_min,
      population_lambda22_max = measurement$lambda22_max,
      population_lambda_matrix_frobenius_rms_dispersion =
        measurement$lambda_matrix_frobenius_rms_dispersion,
      population_theta22_mean = measurement$theta22_mean,
      population_theta22_min = measurement$theta22_min,
      population_theta22_max = measurement$theta22_max,
      population_theta_matrix_frobenius_rms_dispersion =
        measurement$theta_matrix_frobenius_rms_dispersion,
      population_ols_var22_mean = measurement$ols_var22_mean,
      population_ols_var22_min = measurement$ols_var22_min,
      population_ols_var22_max = measurement$ols_var22_max,
      population_ols_cov_matrix_frobenius_rms_dispersion =
        measurement$ols_cov_matrix_frobenius_rms_dispersion,
      slope_variance_marginal = cal$slope_variance,
      residualized_slope_variance = cal$residualized_slope_variance,
      slope_intercept_variance_ratio =
        cal$slope_intercept_variance_ratio,
      G_condition_number = cal$G_condition_number,
      sigma = cal$sigma,
      tau1 = sqrt(cal$slope_variance),
      rho = condition$marginal_rho[[1]],
      standardized_beta = structural$standardized_slope_beta,
      beta1z = structural$beta1_intercept,
      beta2z = structural$beta2_slope,
      structural_r2 = structural$total_structural_r_squared,
      focal_unique_r2 = structural$focal_unique_r_squared,
      outcome_residual_variance = structural$outcome_residual_variance
    )
  }) %>%
    dplyr::bind_rows()

  dplyr::bind_cols(design, calibration_rows) %>%
    dplyr::mutate(
      balance_mode = "heterogeneous_profile",
      min_clus_size = .data$profile_min_clus_size,
      highly_unbalanced_min_clus_size = .data$profile_min_clus_size,
      highly_unbalanced_power = NA_real_,
      r_structure = "iid",
      r_rho = NA_real_,
      study_label = "Heterogeneous Cluster Information (v2 fixed shape)",
      study_structure = "z"
    )
}

#' Crosswalk between ICC inputs and posterior slope reliability.
#'
#' The first family contains the amended VH posterior-reliability targets under
#' sample-SD time scaling. The second contains Lai-style ICC inputs under the
#' RMS time scaling used in the Study 1 supplement. Both hold the random-effect
#' variance ratio fixed within a row, so ICC and posterior reliability become
#' transparent re-expressions of the same regular covariance family.
make_icc_posterior_reliability_crosswalk <- function(
    vh_cluster_sizes = c(3L, 5L, 10L, 25L),
    lai_cluster_sizes = c(3L, 10L, 25L),
    reliability_targets = c(0.25, 0.5, 0.8),
    icc_targets = c(0.05, 0.2, 0.5),
    marginal_rhos = c(-0.5, 0, 0.5),
    variance_ratios = c(0.5, 1, 2)) {
  make_time_design <- function(n, scaling) {
    if (identical(scaling, "vh_sample_sd")) {
      return(make_reliability_time_design(n))
    }
    x <- seq(-1, 1, length.out = n)
    x <- x / sqrt(mean(x^2))
    cbind(intercept = 1, slope = x)
  }
  grids <- dplyr::bind_rows(
    tidyr::crossing(
      design_family = "VH v2 posterior reliability target",
      time_standardization = "vh_sample_sd",
      mean_clus_size = as.integer(vh_cluster_sizes),
      marginal_rho = marginal_rhos,
      slope_intercept_variance_ratio = variance_ratios,
      calibration_metric = "marginal_slope",
      calibration_target_value = reliability_targets
    ),
    tidyr::crossing(
      design_family = "Lai-style ICC reference",
      time_standardization = "lai_rms",
      mean_clus_size = as.integer(lai_cluster_sizes),
      marginal_rho = marginal_rhos,
      slope_intercept_variance_ratio = variance_ratios,
      calibration_metric = "intercept_icc",
      calibration_target_value = icc_targets
    )
  )

  lapply(seq_len(nrow(grids)), function(i) {
    row <- grids[i, , drop = FALSE]
    n_i <- row$mean_clus_size[[1]]
    calibration <- calibrate_shape_preserving_measurement(
      target_value = row$calibration_target_value[[1]],
      Z_list = list(make_time_design(
        n_i,
        row$time_standardization[[1]]
      )),
      R_shape_list = list(diag(n_i)),
      intercept_variance = fixed_params$tau0^2,
      slope_intercept_variance_ratio =
        row$slope_intercept_variance_ratio[[1]],
      intercept_slope_correlation = row$marginal_rho[[1]],
      calibration_target = row$calibration_metric[[1]]
    )
    dplyr::bind_cols(
      row,
      tibble::tibble(
        sigma = calibration$sigma,
        intercept_icc = calibration$intercept_icc,
        marginal_slope_icc = calibration$marginal_slope_icc,
        conditional_slope_icc = calibration$conditional_slope_icc,
        achieved_marginal_reliability =
          calibration$achieved_marginal_slope_reliability,
        achieved_residualized_reliability =
          calibration$achieved_residualized_slope_reliability,
        G_condition_number = calibration$G_condition_number
      )
    )
  }) %>%
    dplyr::bind_rows()
}

#' Build a runnable ICC-versus-posterior-reliability bridge manifest.
#'
#' ICC anchors are chosen separately within each rho/reliability level so the
#' two arms are identical at the reference cluster size m=10. The posterior-R
#' arm recalibrates sigma at every m, whereas the ICC arm holds that reference
#' ICC fixed across m. Differences away from m=10 therefore isolate the choice
#' of what is held constant as cluster information changes.
make_icc_reliability_bridge_design <- function(
    reference_cluster_size = 10L,
    slope_intercept_variance_ratio = 1) {
  target_grid <- tidyr::crossing(
    posterior_reliability_anchor = c(0.25, 0.5, 0.8),
    marginal_rho = c(-0.5, 0, 0.5)
  )
  anchor_rows <- lapply(seq_len(nrow(target_grid)), function(i) {
    row <- target_grid[i, , drop = FALSE]
    anchor <- calibrate_shape_preserving_measurement(
      target_value = row$posterior_reliability_anchor[[1]],
      Z_list = list(make_reliability_time_design(reference_cluster_size)),
      R_shape_list = list(diag(reference_cluster_size)),
      intercept_variance = fixed_params$tau0^2,
      slope_intercept_variance_ratio =
        slope_intercept_variance_ratio,
      intercept_slope_correlation = row$marginal_rho[[1]],
      calibration_target = "marginal_slope"
    )
    dplyr::bind_cols(
      row,
      tibble::tibble(icc_anchor = anchor$intercept_icc)
    )
  }) %>%
    dplyr::bind_rows()

  first_stage_grid <- tidyr::crossing(
    calibration_arm = c(
      "posterior_reliability_targeted", "icc_anchored_at_m10"
    ),
    mean_clus_size = c(3L, 10L, 25L),
    anchor_rows,
    balance_mode = "balanced",
    min_clus_size = 2L,
    highly_unbalanced_min_clus_size = 2L,
    highly_unbalanced_power = 3,
    r_structure = "iid",
    r_rho = NA_real_
  )
  calibrated <- lapply(seq_len(nrow(first_stage_grid)), function(i) {
    condition <- first_stage_grid[i, , drop = FALSE]
    n_i <- condition$mean_clus_size[[1]]
    target_metric <- if (identical(
      condition$calibration_arm[[1]],
      "posterior_reliability_targeted"
    )) "marginal_slope" else "intercept_icc"
    target_value <- if (identical(target_metric, "marginal_slope")) {
      condition$posterior_reliability_anchor[[1]]
    } else {
      condition$icc_anchor[[1]]
    }
    cal <- calibrate_shape_preserving_measurement(
      target_value = target_value,
      Z_list = list(make_reliability_time_design(n_i)),
      R_shape_list = list(diag(n_i)),
      intercept_variance = fixed_params$tau0^2,
      slope_intercept_variance_ratio =
        slope_intercept_variance_ratio,
      intercept_slope_correlation = condition$marginal_rho[[1]],
      calibration_target = target_metric
    )
    dplyr::bind_cols(
      condition,
      tibble::tibble(
        study_version = "bridge",
        calibration_version = "shape_preserving_v2",
        calibration_metric = target_metric,
        calibration_target = target_metric,
        calibration_target_value = target_value,
        achieved_calibration_value = cal$achieved_calibration_value,
        covariance_shape_fixed = TRUE,
        target_reliability =
          condition$posterior_reliability_anchor[[1]],
        achieved_reliability =
          cal$achieved_marginal_slope_reliability,
        achieved_partial_reliability =
          cal$achieved_residualized_slope_reliability,
        intercept_icc = cal$intercept_icc,
        marginal_slope_icc = cal$marginal_slope_icc,
        conditional_slope_icc = cal$conditional_slope_icc,
        calibration_tau0 = fixed_params$tau0,
        calibration_tau0_sq = cal$intercept_variance,
        calibration_tau1_sq = cal$slope_variance,
        sigma = cal$sigma,
        slope_variance_marginal = cal$slope_variance,
        residualized_slope_variance = cal$residualized_slope_variance,
        slope_intercept_variance_ratio =
          cal$slope_intercept_variance_ratio,
        G_condition_number = cal$G_condition_number,
        tau1 = sqrt(cal$slope_variance),
        rho = condition$marginal_rho[[1]],
        G_marginal = list(cal$G_marginal)
      )
    )
  }) %>%
    dplyr::bind_rows()

  design <- tidyr::crossing(
    calibrated,
    study = "iccbridge",
    num_clus = 100L,
    standardized_beta_target = c(0, 0.4),
    structural_target = "intercept_slope",
    study_label = "ICC versus posterior reliability bridge",
    study_structure = "z"
  )
  structural_rows <- lapply(seq_len(nrow(design)), function(i) {
    condition <- design[i, , drop = FALSE]
    structural <- calibrate_blup_predictor_effect(
      G_marginal = condition$G_marginal[[1]],
      standardized_slope_beta = condition$standardized_beta_target[[1]],
      structural_target = condition$structural_target[[1]],
      nuisance_intercept_beta = fixed_params$beta1z,
      outcome_variance = fixed_params$z_variance
    )
    tibble::tibble(
      standardized_beta = structural$standardized_slope_beta,
      beta1z = structural$beta1_intercept,
      beta2z = structural$beta2_slope,
      structural_r2 = structural$total_structural_r_squared,
      focal_unique_r2 = structural$focal_unique_r_squared,
      outcome_residual_variance = structural$outcome_residual_variance
    )
  }) %>%
    dplyr::bind_rows()

  bridge_design <- dplyr::bind_cols(
    dplyr::select(design, -G_marginal),
    structural_rows
  )

  # Pair the two calibration arms before condition IDs are assigned. The
  # paired arms share every design factor except `calibration_arm`; using this
  # stable group in the replication seed gives them common random numbers.
  # At m = 10 the two calibrated DGMs are identical, so this creates an exact
  # replication-by-replication anchor check.
  bridge_design %>%
    dplyr::mutate(
      bridge_pair_label = sprintf(
        "m%02d_rho%+.2f_rel%.2f_beta%.2f",
        as.integer(mean_clus_size),
        as.numeric(marginal_rho),
        as.numeric(posterior_reliability_anchor),
        as.numeric(standardized_beta_target)
      ),
      bridge_pair_id = as.integer(factor(
        bridge_pair_label,
        levels = sort(unique(bridge_pair_label))
      )),
      # Keep this namespace separate from ordinary condition IDs while
      # remaining comfortably inside R's integer seed range.
      simulation_seed_group = 10000L + bridge_pair_id
    )
}

#' Build Study 5: a matched-reliability calibration bridge.
#'
#' @details
#' All three arms target reliability .25, but they reach it through different
#' one-dimensional calibrations:
#'
#' - `current_g22` reproduces Study 2 by fixing the intercept variance
#'   `tau0^2 = .81` and `sigma = 1`, then solving only for the slope variance
#'   `tau1^2` using marginal slope reliability. `current_g22` is retained as
#'   the old code label for this legacy `solve-G_22` arm;
#' - `shape_preserving_marginal` fixes `tau1^2/tau0^2 = 1` and solves for
#'   `sigma` using the same marginal reliability;
#' - `shape_preserving_partial` fixes the same covariance shape and solves for
#'   `sigma` using reliability of the slope residualized on the intercept.
#'
#' Holding the estimator bundle, structural effect, sample size, and target
#' information constant makes the arm contrasts diagnostic of covariance
#' geometry versus the reliability definition.
make_study5_design <- function() {
  target <- 0.25
  shape_variance <- fixed_params$tau0^2
  first_stage_grid <- tidyr::crossing(
    calibration_arm = c(
      "current_g22",
      "shape_preserving_marginal",
      "shape_preserving_partial"
    ),
    mean_clus_size = c(10L, 25L),
    marginal_rho = c(0, 0.5),
    target_calibration_reliability = target,
    target_reliability = target,
    balance_mode = "balanced",
    min_clus_size = 2L,
    highly_unbalanced_min_clus_size = 2L,
    highly_unbalanced_power = 3,
    r_structure = "iid",
    r_rho = NA_real_
  )

  calibrated_first_stage <- lapply(seq_len(nrow(first_stage_grid)), function(i) {
    condition <- first_stage_grid[i, , drop = FALSE]
    reference <- make_reliability_reference_design(
      mean_n_trial = condition$mean_clus_size[[1]],
      sigma = 1,
      balance_mode = condition$balance_mode[[1]],
      min_n_trial = condition$min_clus_size[[1]],
      highly_unbalanced_min_n_trial =
        condition$highly_unbalanced_min_clus_size[[1]],
      highly_unbalanced_power = condition$highly_unbalanced_power[[1]],
      r_spec = condition_to_r_spec(condition),
      n_reference = 1001L
    )

    arm <- condition$calibration_arm[[1]]
    if (identical(arm, "current_g22")) {
      calibration_metric <- "marginal_slope"
      sigma <- 1
      marginal_calibration <- calibrate_slope_variance(
        target_reliability = target,
        Z_list = reference$Z_list,
        R_list = reference$R_list,
        weights = reference$count_weights,
        intercept_variance = shape_variance,
        intercept_slope_correlation = condition$marginal_rho[[1]]
      )
      G <- marginal_calibration$G_marginal
      achieved_marginal <- marginal_calibration$achieved_reliability
      achieved_partial <- expected_residualized_slope_reliability(
        G,
        reference$Z_list,
        reference$R_list,
        weights = reference$count_weights
      )
    } else {
      calibration_metric <- if (identical(
        arm, "shape_preserving_marginal"
      )) "marginal_slope" else "residualized_slope"
      shape_calibration <- calibrate_shape_preserving_measurement(
        target_value = target,
        Z_list = reference$Z_list,
        R_shape_list = reference$R_list,
        weights = reference$count_weights,
        intercept_variance = shape_variance,
        slope_intercept_variance_ratio = 1,
        intercept_slope_correlation = condition$marginal_rho[[1]],
        calibration_target = calibration_metric
      )
      G <- shape_calibration$G_marginal
      sigma <- shape_calibration$sigma
      achieved_marginal <-
        shape_calibration$achieved_marginal_slope_reliability
      achieved_partial <-
        shape_calibration$achieved_residualized_slope_reliability
    }

    residualized_variance <-
      G[2L, 2L] - G[1L, 2L]^2 / G[1L, 1L]
    achieved_calibration <- if (identical(
      calibration_metric, "marginal_slope"
    )) achieved_marginal else achieved_partial

    dplyr::bind_cols(
      condition,
      tibble::tibble(
        calibration_metric = calibration_metric,
        calibration_target = calibration_metric,
        calibration_target_value = target,
        calibration_version = if (identical(arm, "current_g22")) {
          "legacy_solve_tau1_sq"
        } else {
          "shape_preserving_v2"
        },
        covariance_shape_fixed = !identical(arm, "current_g22"),
        calibration_tau0 = fixed_params$tau0,
        sigma = sigma,
        achieved_calibration_reliability = achieved_calibration,
        achieved_reliability = achieved_marginal,
        achieved_partial_reliability = achieved_partial,
        slope_variance_marginal = G[2L, 2L],
        residualized_slope_variance = residualized_variance,
        slope_intercept_variance_ratio = G[2L, 2L] / G[1L, 1L],
        G_condition_number = kappa(G, exact = TRUE),
        intercept_icc = G[1L, 1L] / (G[1L, 1L] + sigma^2),
        marginal_slope_icc = G[2L, 2L] / (G[2L, 2L] + sigma^2),
        conditional_slope_icc =
          residualized_variance / (residualized_variance + sigma^2),
        tau1 = sqrt(G[2L, 2L]),
        rho = condition$marginal_rho[[1]],
        G_marginal = list(G),
        reference_mean_clus_size = reference$mean_trial_count,
        reference_min_clus_size = reference$min_trial_count,
        reference_max_clus_size = reference$max_trial_count,
        calibration_reference_n = length(reference$trial_counts)
      )
    )
  }) %>%
    dplyr::bind_rows()

  design <- tidyr::crossing(
    calibrated_first_stage,
    study = "study5",
    num_clus = 100L,
    standardized_beta_target = c(0, 0.2, 0.4),
    structural_target = "intercept_slope",
    study_label = "Matched Reliability Calibration Bridge",
    study_structure = "z"
  )

  structural_rows <- lapply(seq_len(nrow(design)), function(i) {
    condition <- design[i, , drop = FALSE]
    structural <- calibrate_blup_predictor_effect(
      G_marginal = condition$G_marginal[[1]],
      standardized_slope_beta = condition$standardized_beta_target[[1]],
      structural_target = condition$structural_target[[1]],
      nuisance_intercept_beta = fixed_params$beta1z,
      outcome_variance = fixed_params$z_variance
    )
    tibble::tibble(
      standardized_beta = structural$standardized_slope_beta,
      beta1z = structural$beta1_intercept,
      beta2z = structural$beta2_slope,
      structural_r2 = structural$total_structural_r_squared,
      focal_unique_r2 = structural$focal_unique_r_squared,
      outcome_residual_variance = structural$outcome_residual_variance
    )
  }) %>%
    dplyr::bind_rows()

  dplyr::bind_cols(
    dplyr::select(design, -G_marginal),
    structural_rows
  )
}

#' Build the compact covariance-shape robustness bridge.
#'
#' This module deliberately does not cross the variance ratio with every
#' factor in Studies 1--4.  Instead, it reuses five pre-specified sentinel
#' DGMs: one moderate-information cell from each amended study and one
#' low-information Study-2 stress cell.  Within a sentinel, only the marginal
#' slope/intercept variance ratio and the null/non-null structural effect vary.
#' The Level-1 residual scale is recalibrated at every ratio, so marginal
#' posterior slope reliability remains fixed and the contrast isolates
#' covariance shape rather than measurement-information level.
#'
#' Ratio-one conditions are intentionally rerun.  All three ratio arms within
#' a DGM/effect pair share `simulation_seed_group`, providing common random
#' numbers for direct replication-level comparisons without touching the
#' completed primary-v2 output directories.
make_variance_ratio_sensitivity_design <- function(
    variance_ratios = c(0.5, 1, 2),
    standardized_beta_targets = c(0, 0.4),
    simulation_seed_group_offset = 12000L) {
  variance_ratios <- sort(unique(as.numeric(variance_ratios)))
  standardized_beta_targets <- sort(unique(as.numeric(
    standardized_beta_targets
  )))
  if (
    length(variance_ratios) == 0L ||
      any(!is.finite(variance_ratios)) ||
      any(variance_ratios <= 0)
  ) {
    stop("`variance_ratios` must contain positive finite values.")
  }
  if (
    length(standardized_beta_targets) == 0L ||
      any(!is.finite(standardized_beta_targets))
  ) {
    stop("`standardized_beta_targets` must contain finite values.")
  }
  if (
    length(simulation_seed_group_offset) != 1L ||
      !is.finite(simulation_seed_group_offset) ||
      simulation_seed_group_offset < 1L
  ) {
    stop("`simulation_seed_group_offset` must be a positive finite scalar.")
  }

  block_levels <- c(
    "study1_outcome_moderate",
    "study2_predictor_moderate",
    "study3_dual_moderate",
    "study4_heterogeneous_moderate",
    "study2_predictor_low_reliability"
  )
  block_labels <- c(
    study1_outcome_moderate =
      "Study 1: slope outcome, moderate information",
    study2_predictor_moderate =
      "Study 2: slope predictor, moderate information",
    study3_dual_moderate =
      "Study 3: slope predictor and outcome, moderate information",
    study4_heterogeneous_moderate =
      "Study 4: severe heterogeneous information, moderate mean reliability",
    study2_predictor_low_reliability =
      "Study 2: slope predictor, low-information stress"
  )

  tag_block <- function(data, block) {
    dplyr::mutate(
      data,
      sensitivity_block = block,
      sensitivity_block_label = unname(block_labels[[block]])
    )
  }

  build_ratio <- function(variance_ratio) {
    study1 <- make_study1_v2_design(
      slope_intercept_variance_ratio = variance_ratio
    ) %>%
      dplyr::filter(
        .data$num_clus == 50L,
        .data$mean_clus_size == 3L,
        .data$target_reliability == 0.5,
        .data$marginal_rho == 0.5,
        .data$standardized_beta_target %in% standardized_beta_targets
      ) %>%
      tag_block("study1_outcome_moderate")

    study2_base <- make_study2_v2_design(
      slope_intercept_variance_ratio = variance_ratio
    )
    study2_moderate <- study2_base %>%
      dplyr::filter(
        .data$num_clus == 50L,
        .data$mean_clus_size == 3L,
        .data$target_reliability == 0.5,
        .data$marginal_rho == -0.5,
        .data$standardized_beta_target %in% standardized_beta_targets
      ) %>%
      tag_block("study2_predictor_moderate")
    study2_low <- study2_base %>%
      dplyr::filter(
        .data$num_clus == 50L,
        .data$mean_clus_size == 3L,
        .data$target_reliability == 0.25,
        .data$marginal_rho == -0.5,
        .data$standardized_beta_target %in% standardized_beta_targets
      ) %>%
      tag_block("study2_predictor_low_reliability")

    study3 <- make_study3_v2_design(
      slope_intercept_variance_ratio = variance_ratio,
      target_reliabilities = 0.5,
      standardized_beta_targets = standardized_beta_targets
    ) %>%
      dplyr::filter(
        .data$num_clus == 50L,
        .data$mean_clus_size_y == 3L,
        .data$mean_clus_size_q == 3L,
        .data$target_reliability_y == 0.5,
        .data$target_reliability_q == 0.5,
        .data$marginal_rho == 0.5
      ) %>%
      dplyr::mutate(
        mean_clus_size = .data$mean_clus_size_y,
        target_reliability = .data$target_reliability_y,
        achieved_reliability = .data$achieved_reliability_y,
        achieved_partial_reliability =
          .data$achieved_partial_reliability_y,
        intercept_icc = .data$intercept_icc_y,
        marginal_slope_icc = .data$marginal_slope_icc_y,
        conditional_slope_icc = .data$conditional_slope_icc_y,
        sigma = .data$sigma_y
      ) %>%
      tag_block("study3_dual_moderate")

    study4 <- make_study4_v2_design(
      slope_intercept_variance_ratio = variance_ratio
    ) %>%
      dplyr::filter(
        .data$num_clus == 50L,
        .data$target_reliability == 0.5,
        .data$information_profile == "severe",
        .data$marginal_rho == 0.5,
        !.data$is_falsification_control,
        .data$standardized_beta_target %in% standardized_beta_targets
      ) %>%
      tag_block("study4_heterogeneous_moderate")

    dplyr::bind_rows(
      study1,
      study2_moderate,
      study3,
      study4,
      study2_low
    ) %>%
      dplyr::mutate(
        simulation_module = "variance_ratio_sensitivity_v1",
        variance_ratio_arm = paste0(
          "tau1_sq_over_tau0_sq_",
          format(variance_ratio, trim = TRUE, scientific = FALSE)
        ),
        sensitivity_reference_ratio = 1,
        sensitivity_is_low_reliability_stress =
          .data$sensitivity_block ==
            "study2_predictor_low_reliability"
      )
  }

  design <- lapply(variance_ratios, build_ratio) %>%
    dplyr::bind_rows()
  expected_per_ratio <- length(block_levels) *
    length(standardized_beta_targets)
  observed_per_ratio <- design %>%
    dplyr::count(.data$slope_intercept_variance_ratio, name = "n")
  if (
    nrow(design) != length(variance_ratios) * expected_per_ratio ||
      any(observed_per_ratio$n != expected_per_ratio)
  ) {
    stop(
      "Variance-ratio sensitivity sentinel selection is incomplete: expected ",
      expected_per_ratio, " rows per ratio."
    )
  }

  pair_index <- design %>%
    dplyr::distinct(
      .data$sensitivity_block,
      .data$standardized_beta_target
    ) %>%
    dplyr::mutate(
      block_order = match(.data$sensitivity_block, block_levels)
    ) %>%
    dplyr::arrange(
      .data$block_order,
      .data$standardized_beta_target
    ) %>%
    dplyr::mutate(
      variance_ratio_pair_id = dplyr::row_number(),
      variance_ratio_pair_label = paste0(
        .data$sensitivity_block,
        "__beta_",
        format(
          .data$standardized_beta_target,
          trim = TRUE,
          scientific = FALSE
        )
      )
    ) %>%
    dplyr::select(-"block_order")

  design %>%
    dplyr::left_join(
      pair_index,
      by = c("sensitivity_block", "standardized_beta_target")
    ) %>%
    dplyr::mutate(
      sensitivity_block = factor(
        .data$sensitivity_block,
        levels = block_levels
      ),
      simulation_seed_group = as.integer(
        simulation_seed_group_offset + .data$variance_ratio_pair_id
      )
    ) %>%
    dplyr::arrange(
      .data$sensitivity_block,
      .data$standardized_beta_target,
      .data$slope_intercept_variance_ratio
    )
}

# Explicit legacy aliases make it possible to call the original one-coordinate
# solve directly even after the amended wrappers become the preferred designs.
make_study1_legacy_design <- make_study1_design
make_study2_legacy_design <- make_study2_design
make_study3_legacy_design <- make_study3_design
make_study4_legacy_design <- make_study4_design

study_condition_counts <- function() {
  c(
    study1 = 5L * 4L * 3L * 3L * 4L,
    study2 = 4L * 3L * 3L * 5L * 4L,
    study3 = 4L * 3L * 3L * 2L * 2L * 2L * 3L,
    study4 = (3L * 3L * 3L * 2L * 2L) + 3L,
    study0 = 1L,
    study5 = 3L * 2L * 2L * 3L,
    study1v2 = 5L * 4L * 3L * 3L * 4L,
    study2v2 = 4L * 3L * 3L * 5L * 4L,
    study3v2 = 4L * 3L * 3L * 2L * 2L * 2L * 3L,
    study4v2 = (3L * 3L * 3L * 2L * 2L) + 3L,
    iccbridge = 2L * 3L * 3L * 3L * 2L,
    ratiosensitivity = 5L * 3L * 2L
  )
}

select_design <- function(study_arg = "all", max_conditions = NA_integer_) {
  study_arg <- tolower(trimws(as.character(study_arg[[1]])))
  legacy_primary <- paste0("study", 1:5)
  amended_primary <- paste0("study", 1:4, "v2")
  requested <- switch(
    study_arg,
    all = legacy_primary,
    legacy = legacy_primary,
    alllegacy = legacy_primary,
    v2 = amended_primary,
    amended = amended_primary,
    allv2 = amended_primary,
    allversions = c(legacy_primary, amended_primary, "iccbridge"),
    all_versions = c(legacy_primary, amended_primary, "iccbridge"),
    {
      tokens <- trimws(unlist(strsplit(study_arg, ",", fixed = TRUE)))
      vapply(tokens, function(token) {
        token <- sub("^study", "", token)
        switch(
          token,
          `0` = "study0",
          `1` = "study1",
          `2` = "study2",
          `3` = "study3",
          `4` = "study4",
          `5` = "study5",
          `1v2` = "study1v2",
          `2v2` = "study2v2",
          `3v2` = "study3v2",
          `4v2` = "study4v2",
          bridge = "iccbridge",
          iccbridge = "iccbridge",
          sensitivity = "ratiosensitivity",
          ratiosensitivity = "ratiosensitivity",
          ratio_sensitivity = "ratiosensitivity",
          variance_ratio = "ratiosensitivity",
          variance_ratio_sensitivity = "ratiosensitivity",
          shapebridge = "ratiosensitivity",
          NA_character_
        )
      }, character(1L))
    }
  )
  if (length(requested) == 0L || anyNA(requested)) {
    stop(paste0(
      "No study conditions selected. Use `all` for the original Studies 1--5, ",
      "`allv2` for amended Studies 1--4, `1`--`5`, `1v2`--`4v2`, ",
      "`iccbridge`, `ratiosensitivity`, or a comma-separated combination."
    ))
  }

  builders <- list(
    study0 = make_study0_design,
    study1 = make_study1_legacy_design,
    study2 = make_study2_legacy_design,
    study3 = make_study3_legacy_design,
    study4 = make_study4_legacy_design,
    study5 = make_study5_design,
    study1v2 = make_study1_v2_design,
    study2v2 = make_study2_v2_design,
    study3v2 = make_study3_v2_design,
    study4v2 = make_study4_v2_design,
    iccbridge = make_icc_reliability_bridge_design,
    ratiosensitivity = make_variance_ratio_sensitivity_design
  )
  counts <- study_condition_counts()
  offsets <- c(0L, cumsum(counts)[-length(counts)])
  names(offsets) <- names(counts)
  built <- lapply(requested, function(study_name) {
    study_design <- builders[[study_name]]()
    if (is.null(study_design) || nrow(study_design) == 0L) {
      return(tibble::tibble())
    }
    if (nrow(study_design) != counts[[study_name]]) {
      stop(
        "Canonical condition count is out of sync for ", study_name,
        ": expected ", counts[[study_name]], ", got ", nrow(study_design), "."
      )
    }
    if (!("study_version" %in% names(study_design))) {
      study_design$study_version <- if (study_name %in% paste0("study", 1:4)) {
        "legacy"
      } else if (identical(study_name, "study5")) {
        "bridge"
      } else {
        "legacy_smoke"
      }
    }
    if (!("calibration_version" %in% names(study_design))) {
      study_design$calibration_version <- if (
        study_name %in% paste0("study", 1:4)
      ) {
        "legacy_solve_tau1_sq"
      } else if (identical(study_name, "study5")) {
        "study5_calibration_bridge"
      } else {
        "legacy_solve_tau1_sq"
      }
    }
    study_design %>%
      dplyr::mutate(
        condition_id = offsets[[study_name]] + seq_len(dplyr::n())
      ) %>%
      dplyr::relocate(condition_id)
  })
  design <- dplyr::bind_rows(built)

  if (!is.na(max_conditions) && max_conditions > 0L) {
    design <- design %>% dplyr::slice_head(n = max_conditions)
  }

  if (nrow(design) == 0L) {
    stop("No study conditions selected.")
  }

  design
}
