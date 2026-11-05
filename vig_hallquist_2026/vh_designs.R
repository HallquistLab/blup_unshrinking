# Generated v2-only release source.

fixed_params <- list(x_mean = 0, x_variance = 1, z_variance = 1, tau0 = 0.9, beta0z = 1.5, beta1z = 0.4, theta0_standardized = 0.4, gamma0_outcome = 1, gamma1_outcome = 0.6, gamma0_predictor = 0, gamma1_predictor = 0.5, gamma0_process_y = 0, gamma1_process_y = 0.5, gamma0_process_q = 0.5, gamma1_process_q = 0.4)

cluster_balance_modes <- function() {
    c("balanced", "unbalanced", "informative_unbalanced")
}

balance_mode_to_sim_arg <- function(balance_mode) {
    switch(as.character(balance_mode[[1]]), balanced = "balanced", unbalanced = FALSE, informative_unbalanced = "highly_unbalanced", stop("Unsupported balance mode: ", balance_mode))
}

residual_structures <- function() {
    tibble::tribble(~r_structure, ~r_rho, "iid", NA_real_, "ar1", 0.3, "ar1", 0.6)
}

calibrate_reliability_design <- function(design, tau0 = 0.9, calibration_reference_n = 1001L) {
    required <- c("target_reliability", "standardized_beta_target", "marginal_rho", "mean_clus_size", "sigma", "balance_mode", "min_clus_size", "highly_unbalanced_min_clus_size", "highly_unbalanced_power", "r_structure", "r_rho")
    missing <- setdiff(required, names(design))
    if (length(missing) > 0L) {
        stop("Reliability design is missing columns: ", paste(missing, collapse = ", "))
    }
    if (!exists("calibrate_random_slope_condition", mode = "function")) {
        stop("`calibrate_random_slope_condition()` is unavailable. ", "Source `R/reliability_calibration.R` before building this grid.")
    }
    calibrated_rows <- lapply(seq_len(nrow(design)), function(i) {
        condition <- design[i, , drop = FALSE]
        calibrated <- calibrate_random_slope_condition(target_reliability = condition$target_reliability[[1]], standardized_beta = condition$standardized_beta_target[[1]], marginal_rho = condition$marginal_rho[[1]], tau0 = tau0, mean_n_trial = condition$mean_clus_size[[1]], sigma = condition$sigma[[1]], balance_mode = condition$balance_mode[[1]], min_n_trial = condition$min_clus_size[[1]], highly_unbalanced_min_n_trial = condition$highly_unbalanced_min_clus_size[[1]], highly_unbalanced_power = condition$highly_unbalanced_power[[1]], 
            r_spec = condition_to_r_spec(condition), n_reference = calibration_reference_n)
        data.frame(calibration_tau0 = tau0, achieved_reliability = calibrated$achieved_reliability, standardized_beta = calibrated$standardized_beta, structural_r2 = calibrated$structural_r_squared, gamma_x_on_slope = calibrated$gamma_x_on_slope, slope_variance_marginal = calibrated$slope_variance_marginal, slope_variance_residual = calibrated$slope_variance_residual, tau1_residual = calibrated$tau1_residual, rho_residual = calibrated$rho_residual, reference_mean_clus_size = calibrated$reference_mean_n_trial, 
            reference_min_clus_size = calibrated$reference_min_n_trial, reference_max_clus_size = calibrated$reference_max_n_trial, calibration_reference_n = calibrated$calibration_reference_n)
    })
    calibrated <- dplyr::bind_rows(calibrated_rows)
    dplyr::bind_cols(design, calibrated) %>% dplyr::mutate(rho = .data$marginal_rho, tau1 = .data$tau1_residual)
}

study4_profile_spec <- function(information_profile) {
    information_profile <- as.character(information_profile[[1]])
    switch(information_profile, homogeneous = list(sizes = 10L, weights = 1, information_matched = FALSE), moderate = list(sizes = c(5L, 15L), weights = c(0.5, 0.5), information_matched = FALSE), severe = list(sizes = c(3L, 17L), weights = c(0.5, 0.5), information_matched = FALSE), severe_information_matched = list(sizes = c(3L, 17L), weights = c(0.5, 0.5), information_matched = TRUE), stop("Unsupported Study 4 information profile: ", information_profile))
}

study4_time_design <- function(cluster_size, information_matched = FALSE) {
    Z <- make_reliability_time_design(cluster_size)
    if (isTRUE(information_matched)) {
        Z[, "slope"] <- Z[, "slope"] * sqrt(9/(cluster_size - 1))
    }
    Z
}

study4_weighted_quantile <- function(x, weights, probs) {
    keep <- is.finite(x) & is.finite(weights) & weights > 0
    x <- as.numeric(x[keep])
    weights <- as.numeric(weights[keep])
    if (length(x) == 0L) {
        return(rep(NA_real_, length(probs)))
    }
    ord <- order(x)
    x <- x[ord]
    cumulative_weight <- cumsum(weights[ord]/sum(weights))
    vapply(probs, function(prob) {
        x[[which(cumulative_weight >= min(1, max(0, prob)))[[1L]]]]
    }, numeric(1L))
}

study4_matrix_rms_dispersion <- function(matrices, weights = NULL) {
    if (is.null(weights)) {
        weights <- rep(1, length(matrices))
    }
    valid <- vapply(matrices, function(mat) is.matrix(mat) && all(dim(mat) == c(2L, 2L)) && all(is.finite(mat)), logical(1L)) & is.finite(weights) & weights > 0
    if (!any(valid)) {
        return(NA_real_)
    }
    matrices <- matrices[valid]
    weights <- weights[valid]/sum(weights[valid])
    mean_matrix <- Reduce(`+`, Map(function(mat, weight) mat * weight, matrices, weights))
    sqrt(sum(vapply(seq_along(matrices), function(i) weights[[i]] * sum((matrices[[i]] - mean_matrix)^2), numeric(1L))))
}

study4_measurement_matrix_summary <- function(lambda_matrices, theta_matrices, score_error_covariances, weights = NULL) {
    if (is.null(weights)) {
        weights <- rep(1, length(lambda_matrices))
    }
    weights <- weights/sum(weights)
    lambda22 <- vapply(lambda_matrices, function(mat) mat[2L, 2L], numeric(1L))
    theta22 <- vapply(theta_matrices, function(mat) mat[2L, 2L], numeric(1L))
    ols_var22 <- vapply(score_error_covariances, function(mat) mat[2L, 2L], numeric(1L))
    list(lambda22_mean = stats::weighted.mean(lambda22, weights), lambda22_min = min(lambda22), lambda22_max = max(lambda22), lambda_matrix_frobenius_rms_dispersion = study4_matrix_rms_dispersion(lambda_matrices, weights), theta22_mean = stats::weighted.mean(theta22, weights), theta22_min = min(theta22), theta22_max = max(theta22), theta_matrix_frobenius_rms_dispersion = study4_matrix_rms_dispersion(theta_matrices, weights), ols_var22_mean = stats::weighted.mean(ols_var22, weights), ols_var22_min = min(ols_var22), 
        ols_var22_max = max(ols_var22), ols_cov_matrix_frobenius_rms_dispersion = study4_matrix_rms_dispersion(score_error_covariances, weights))
}

calibrate_study4_profile <- function(information_profile, target_reliability, marginal_rho, sigma = 1) {
    spec <- study4_profile_spec(information_profile)
    Z_list <- lapply(spec$sizes, study4_time_design, information_matched = spec$information_matched)
    R_list <- lapply(spec$sizes, make_reliability_residual_covariance, sigma = sigma, r_spec = list(structure = "iid"))
    calibration <- calibrate_slope_variance(target_reliability = target_reliability, Z_list = Z_list, R_list = R_list, weights = spec$weights, intercept_variance = fixed_params$tau0^2, intercept_slope_correlation = marginal_rho)
    type_reliability <- vapply(calibration$posterior_covariances, function(V_i) {
        1 - V_i[2L, 2L]/calibration$G_marginal[2L, 2L]
    }, numeric(1L))
    weighted_mean <- stats::weighted.mean(type_reliability, spec$weights)
    weighted_variance <- stats::weighted.mean((type_reliability - weighted_mean)^2, spec$weights)
    reliability_quartiles <- study4_weighted_quantile(type_reliability, spec$weights, c(0.25, 0.75))
    population_measurement <- Map(function(Z_i, R_i) {
        sigma_i <- Z_i %*% calibration$G_marginal %*% t(Z_i) + R_i
        A_i <- t(solve(sigma_i, Z_i %*% calibration$G_marginal))
        R_inv_Z <- solve(R_i, Z_i)
        list(lambda = A_i %*% Z_i, theta = A_i %*% R_i %*% t(A_i), score_error_covariance = solve(crossprod(Z_i, R_inv_Z)))
    }, Z_list, R_list)
    measurement_diagnostics <- study4_measurement_matrix_summary(lambda_matrices = lapply(population_measurement, `[[`, "lambda"), theta_matrices = lapply(population_measurement, `[[`, "theta"), score_error_covariances = lapply(population_measurement, `[[`, "score_error_covariance"), weights = spec$weights)
    list(calibration = calibration, profile = spec, type_reliability = type_reliability, reliability_mean = weighted_mean, reliability_sd = sqrt(weighted_variance), reliability_iqr = diff(reliability_quartiles), reliability_min = min(type_reliability), reliability_max = max(type_reliability), measurement_diagnostics = measurement_diagnostics)
}

calibrate_shape_preserving_first_stage_grid <- function(first_stage_grid, slope_intercept_variance_ratio = 1, calibration_target = "marginal_slope", calibration_reference_n = 1001L) {
    required <- c("mean_clus_size", "target_reliability", "marginal_rho", "balance_mode", "min_clus_size", "highly_unbalanced_min_clus_size", "highly_unbalanced_power", "r_structure", "r_rho")
    missing <- setdiff(required, names(first_stage_grid))
    if (length(missing) > 0L) {
        stop("Shape-preserving first-stage grid is missing columns: ", paste(missing, collapse = ", "))
    }
    lapply(seq_len(nrow(first_stage_grid)), function(i) {
        condition <- first_stage_grid[i, , drop = FALSE]
        reference <- make_reliability_reference_design(mean_n_trial = condition$mean_clus_size[[1]], sigma = 1, balance_mode = condition$balance_mode[[1]], min_n_trial = condition$min_clus_size[[1]], highly_unbalanced_min_n_trial = condition$highly_unbalanced_min_clus_size[[1]], highly_unbalanced_power = condition$highly_unbalanced_power[[1]], r_spec = condition_to_r_spec(condition), n_reference = calibration_reference_n)
        calibration <- calibrate_shape_preserving_measurement(target_value = condition$target_reliability[[1]], Z_list = reference$Z_list, R_shape_list = reference$R_list, weights = reference$count_weights, intercept_variance = fixed_params$tau0^2, slope_intercept_variance_ratio = slope_intercept_variance_ratio, intercept_slope_correlation = condition$marginal_rho[[1]], calibration_target = calibration_target)
        dplyr::bind_cols(condition, tibble::tibble(study_version = "v2", calibration_version = calibration$calibration_version, calibration_arm = paste0("shape_preserving_", calibration$calibration_target), calibration_metric = calibration$calibration_target, calibration_target = calibration$calibration_target, calibration_target_value = calibration$calibration_target_value, achieved_calibration_value = calibration$achieved_calibration_value, covariance_shape_fixed = TRUE, calibration_tau0 = fixed_params$tau0, 
            calibration_tau0_sq = calibration$intercept_variance, calibration_tau1_sq = calibration$slope_variance, achieved_reliability = calibration$achieved_marginal_slope_reliability, achieved_partial_reliability = calibration$achieved_residualized_slope_reliability, intercept_icc = calibration$intercept_icc, marginal_slope_icc = calibration$marginal_slope_icc, conditional_slope_icc = calibration$conditional_slope_icc, slope_variance_marginal = calibration$slope_variance, residualized_slope_variance = calibration$residualized_slope_variance, 
            slope_intercept_variance_ratio = calibration$slope_intercept_variance_ratio, G_condition_number = calibration$G_condition_number, sigma = calibration$sigma, tau1_marginal = sqrt(calibration$slope_variance), rho = condition$marginal_rho[[1]], G_marginal = list(calibration$G_marginal), reference_mean_clus_size = reference$mean_trial_count, reference_min_clus_size = reference$min_trial_count, reference_max_clus_size = reference$max_trial_count, calibration_reference_n = length(reference$trial_counts)))
    }) %>% dplyr::bind_rows()
}

make_study1_v2_design <- function(slope_intercept_variance_ratio = 1) {
    first_stage <- tidyr::crossing(mean_clus_size = c(3L, 5L, 10L, 25L), target_reliability = c(0.25, 0.5, 0.8), marginal_rho = c(-0.5, 0, 0.5), balance_mode = "balanced", min_clus_size = 2L, highly_unbalanced_min_clus_size = 2L, highly_unbalanced_power = 3, r_structure = "iid", r_rho = NA_real_) %>% calibrate_shape_preserving_first_stage_grid(slope_intercept_variance_ratio = slope_intercept_variance_ratio)
    design <- tidyr::crossing(first_stage, study = "study1v2", num_clus = c(30L, 50L, 100L, 150L, 300L), standardized_beta_target = c(0, 0.2, 0.4, 0.6), study_label = "BLUP as Outcome (v2 fixed shape)", study_structure = "w")
    structural_rows <- lapply(seq_len(nrow(design)), function(i) {
        condition <- design[i, , drop = FALSE]
        structural <- decompose_structural_slope(list(G_marginal = condition$G_marginal[[1]]), standardized_beta = condition$standardized_beta_target[[1]])
        tibble::tibble(standardized_beta = structural$standardized_beta, structural_r2 = structural$structural_r_squared, beta1w = structural$gamma_x_on_slope, slope_variance_residual = structural$slope_variance_residual, tau1 = structural$slope_sd_residual, rho_residual = structural$intercept_slope_correlation_residual, structural_residual_G_min_eigen = min(eigen(structural$G_residual, symmetric = TRUE, only.values = TRUE)$values))
    }) %>% dplyr::bind_rows()
    dplyr::bind_cols(dplyr::select(design, -G_marginal, -tau1_marginal), structural_rows)
}

make_study2_v2_design <- function(slope_intercept_variance_ratio = 1) {
    first_stage <- tidyr::crossing(mean_clus_size = c(3L, 5L, 10L, 25L), target_reliability = c(0.25, 0.5, 0.8), marginal_rho = c(-0.5, 0, 0.5), balance_mode = "balanced", min_clus_size = 2L, highly_unbalanced_min_clus_size = 2L, highly_unbalanced_power = 3, r_structure = "iid", r_rho = NA_real_) %>% calibrate_shape_preserving_first_stage_grid(slope_intercept_variance_ratio = slope_intercept_variance_ratio)
    design <- tidyr::crossing(first_stage, study = "study2v2", num_clus = c(30L, 50L, 100L, 150L, 300L), standardized_beta_target = c(0, 0.2, 0.4, 0.6), structural_target = "intercept_slope", study_label = "BLUP as Predictor (v2 fixed shape)", study_structure = "z")
    structural_rows <- lapply(seq_len(nrow(design)), function(i) {
        condition <- design[i, , drop = FALSE]
        structural <- calibrate_blup_predictor_effect(G_marginal = condition$G_marginal[[1]], standardized_slope_beta = condition$standardized_beta_target[[1]], structural_target = condition$structural_target[[1]], nuisance_intercept_beta = fixed_params$beta1z, outcome_variance = fixed_params$z_variance)
        tibble::tibble(standardized_beta = structural$standardized_slope_beta, beta1z = structural$beta1_intercept, beta2z = structural$beta2_slope, structural_r2 = structural$total_structural_r_squared, focal_unique_r2 = structural$focal_unique_r_squared, outcome_residual_variance = structural$outcome_residual_variance)
    }) %>% dplyr::bind_rows()
    dplyr::bind_cols(dplyr::select(design, -G_marginal), structural_rows) %>% dplyr::mutate(tau1 = .data$tau1_marginal) %>% dplyr::select(-tau1_marginal)
}

make_study3_v2_design <- function(slope_intercept_variance_ratio = 1, target_reliabilities = c(0.25, 0.8), standardized_beta_targets = c(0, 0.2, 0.5)) {
    first_stage <- tidyr::crossing(mean_clus_size = c(3L, 5L, 10L), target_reliability = target_reliabilities, marginal_rho = c(0, 0.5), balance_mode = "balanced", min_clus_size = 2L, highly_unbalanced_min_clus_size = 2L, highly_unbalanced_power = 3, r_structure = "iid", r_rho = NA_real_) %>% calibrate_shape_preserving_first_stage_grid(slope_intercept_variance_ratio = slope_intercept_variance_ratio)
    y_calibration <- first_stage %>% dplyr::transmute(mean_clus_size_y = mean_clus_size, target_reliability_y = target_reliability, marginal_rho, sigma_y = sigma, achieved_reliability_y = achieved_reliability, achieved_partial_reliability_y = achieved_partial_reliability, intercept_icc_y = intercept_icc, marginal_slope_icc_y = marginal_slope_icc, conditional_slope_icc_y = conditional_slope_icc, G_condition_number_y = G_condition_number, slope_variance_marginal_y = slope_variance_marginal, residualized_slope_variance_y = residualized_slope_variance, 
        tau1_y = tau1_marginal, G_y_marginal = G_marginal, reference_mean_clus_size_y = reference_mean_clus_size, reference_min_clus_size_y = reference_min_clus_size, reference_max_clus_size_y = reference_max_clus_size)
    q_calibration <- first_stage %>% dplyr::transmute(mean_clus_size_q = mean_clus_size, target_reliability_q = target_reliability, marginal_rho, sigma_q = sigma, achieved_reliability_q = achieved_reliability, achieved_partial_reliability_q = achieved_partial_reliability, intercept_icc_q = intercept_icc, marginal_slope_icc_q = marginal_slope_icc, conditional_slope_icc_q = conditional_slope_icc, G_condition_number_q = G_condition_number, slope_variance_marginal_q = slope_variance_marginal, residualized_slope_variance_q = residualized_slope_variance, 
        tau1_q = tau1_marginal, G_q_marginal = G_marginal, reference_mean_clus_size_q = reference_mean_clus_size, reference_min_clus_size_q = reference_min_clus_size, reference_max_clus_size_q = reference_max_clus_size)
    design <- tidyr::crossing(study = "study3v2", num_clus = c(50L, 100L, 150L, 300L), mean_clus_size_y = c(3L, 5L, 10L), mean_clus_size_q = c(3L, 5L, 10L), target_reliability_y = target_reliabilities, target_reliability_q = target_reliabilities, marginal_rho = c(0, 0.5), standardized_beta_target = standardized_beta_targets, structural_target = "intercept_slope", study_label = "BLUP as Predictor and Outcome (v2 fixed shape)", study_structure = "dual_process") %>% dplyr::left_join(y_calibration, by = c("mean_clus_size_y", 
        "target_reliability_y", "marginal_rho")) %>% dplyr::left_join(q_calibration, by = c("mean_clus_size_q", "target_reliability_q", "marginal_rho"))
    structural_rows <- lapply(seq_len(nrow(design)), function(i) {
        condition <- design[i, , drop = FALSE]
        structural <- calibrate_dual_process_effect(G_predictor = condition$G_y_marginal[[1]], G_outcome = condition$G_q_marginal[[1]], standardized_slope_beta = condition$standardized_beta_target[[1]], structural_target = condition$structural_target[[1]], nuisance_intercept_standardized_beta = fixed_params$theta0_standardized)
        tibble::tibble(standardized_beta = structural$standardized_slope_beta, standardized_theta0 = structural$standardized_intercept_beta, theta0 = structural$theta0_intercept, theta1 = structural$theta1_slope, structural_r2 = structural$total_structural_r_squared, focal_unique_r2 = structural$focal_unique_r_squared, slope_variance_residual_q = structural$outcome_slope_residual_variance, tau1_residual_q = sqrt(structural$outcome_slope_residual_variance), rho_residual_q = structural$outcome_residual_correlation, 
            structural_residual_G_min_eigen = min(eigen(structural$G_outcome_residual, symmetric = TRUE, only.values = TRUE)$values), structural_joint_G_min_eigen = min(eigen(structural$G_joint_marginal, symmetric = TRUE, only.values = TRUE)$values))
    }) %>% dplyr::bind_rows()
    dplyr::bind_cols(dplyr::select(design, -G_y_marginal, -G_q_marginal), structural_rows) %>% dplyr::mutate(study_version = "v2", calibration_version = "shape_preserving_v2", calibration_arm = "shape_preserving_marginal", calibration_metric = "marginal_slope", calibration_target = "marginal_slope", covariance_shape_fixed = TRUE, calibration_tau0 = fixed_params$tau0, calibration_tau0_sq = fixed_params$tau0^2, calibration_tau1_sq = fixed_params$tau0^2 * slope_intercept_variance_ratio, slope_intercept_variance_ratio = slope_intercept_variance_ratio, 
        G_condition_number = .data$G_condition_number_y, balance_mode = "balanced", min_clus_size = 2L, highly_unbalanced_min_clus_size = 2L, highly_unbalanced_power = 3, r_structure = "iid", r_rho = NA_real_)
}

calibrate_study4_profile_v2 <- function(information_profile, target_reliability, marginal_rho, slope_intercept_variance_ratio = 1) {
    spec <- study4_profile_spec(information_profile)
    Z_list <- lapply(spec$sizes, study4_time_design, information_matched = spec$information_matched)
    R_shape_list <- lapply(spec$sizes, diag)
    calibration <- calibrate_shape_preserving_measurement(target_value = target_reliability, Z_list = Z_list, R_shape_list = R_shape_list, weights = spec$weights, intercept_variance = fixed_params$tau0^2, slope_intercept_variance_ratio = slope_intercept_variance_ratio, intercept_slope_correlation = marginal_rho, calibration_target = "marginal_slope")
    type_reliability <- vapply(seq_along(Z_list), function(i) {
        V_i <- posterior_random_effect_covariance(calibration$G_marginal, Z_list[[i]], calibration$residual_covariances[[i]])
        1 - V_i[2L, 2L]/calibration$G_marginal[2L, 2L]
    }, numeric(1L))
    type_partial_reliability <- vapply(seq_along(Z_list), function(i) {
        expected_residualized_slope_reliability(calibration$G_marginal, list(Z_list[[i]]), list(calibration$residual_covariances[[i]]))
    }, numeric(1L))
    weighted_mean <- stats::weighted.mean(type_reliability, spec$weights)
    weighted_variance <- stats::weighted.mean((type_reliability - weighted_mean)^2, spec$weights)
    reliability_quartiles <- study4_weighted_quantile(type_reliability, spec$weights, c(0.25, 0.75))
    population_measurement <- Map(function(Z_i, R_i) {
        sigma_i <- Z_i %*% calibration$G_marginal %*% t(Z_i) + R_i
        A_i <- t(solve(sigma_i, Z_i %*% calibration$G_marginal))
        R_inv_Z <- solve(R_i, Z_i)
        list(lambda = A_i %*% Z_i, theta = A_i %*% R_i %*% t(A_i), score_error_covariance = solve(crossprod(Z_i, R_inv_Z)))
    }, Z_list, calibration$residual_covariances)
    measurement_diagnostics <- study4_measurement_matrix_summary(lambda_matrices = lapply(population_measurement, `[[`, "lambda"), theta_matrices = lapply(population_measurement, `[[`, "theta"), score_error_covariances = lapply(population_measurement, `[[`, "score_error_covariance"), weights = spec$weights)
    list(calibration = calibration, profile = spec, type_reliability = type_reliability, type_partial_reliability = type_partial_reliability, reliability_mean = weighted_mean, reliability_sd = sqrt(weighted_variance), reliability_iqr = diff(reliability_quartiles), reliability_min = min(type_reliability), reliability_max = max(type_reliability), partial_reliability_mean = stats::weighted.mean(type_partial_reliability, spec$weights), partial_reliability_min = min(type_partial_reliability), partial_reliability_max = max(type_partial_reliability), 
        measurement_diagnostics = measurement_diagnostics)
}

make_study4_v2_design <- function(slope_intercept_variance_ratio = 1) {
    primary <- tidyr::crossing(study = "study4v2", num_clus = c(50L, 150L, 300L), target_reliability = c(0.25, 0.5, 0.8), information_profile = c("homogeneous", "moderate", "severe"), marginal_rho = c(0, 0.5), standardized_beta_target = c(0, 0.4), structural_target = "intercept_slope", is_falsification_control = FALSE)
    information_matched_control <- tidyr::crossing(study = "study4v2", num_clus = c(50L, 150L, 300L), target_reliability = 0.5, information_profile = "severe_information_matched", marginal_rho = 0, standardized_beta_target = 0.4, structural_target = "intercept_slope", is_falsification_control = TRUE)
    design <- dplyr::bind_rows(primary, information_matched_control)
    calibration_rows <- lapply(seq_len(nrow(design)), function(i) {
        condition <- design[i, , drop = FALSE]
        calibrated <- calibrate_study4_profile_v2(information_profile = condition$information_profile[[1]], target_reliability = condition$target_reliability[[1]], marginal_rho = condition$marginal_rho[[1]], slope_intercept_variance_ratio = slope_intercept_variance_ratio)
        cal <- calibrated$calibration
        spec <- calibrated$profile
        reliability <- calibrated$type_reliability
        partial_reliability <- calibrated$type_partial_reliability
        measurement <- calibrated$measurement_diagnostics
        structural <- calibrate_blup_predictor_effect(G_marginal = cal$G_marginal, standardized_slope_beta = condition$standardized_beta_target[[1]], structural_target = condition$structural_target[[1]], nuisance_intercept_beta = fixed_params$beta1z, outcome_variance = fixed_params$z_variance)
        tibble::tibble(study_version = "v2", calibration_version = cal$calibration_version, calibration_arm = "shape_preserving_marginal", calibration_metric = cal$calibration_target, calibration_target = cal$calibration_target, calibration_target_value = cal$calibration_target_value, achieved_calibration_value = cal$achieved_calibration_value, covariance_shape_fixed = TRUE, mean_clus_size = stats::weighted.mean(spec$sizes, spec$weights), profile_min_clus_size = min(spec$sizes), profile_max_clus_size = max(spec$sizes), 
            profile_small_clus_size = min(spec$sizes), profile_large_clus_size = max(spec$sizes), profile_small_weight = spec$weights[[1L]], profile_large_weight = if (length(spec$weights) > 1L) {
                spec$weights[[length(spec$weights)]]
            }
            else {
                0
            }, information_matched = spec$information_matched, calibration_tau0 = fixed_params$tau0, calibration_tau0_sq = cal$intercept_variance, calibration_tau1_sq = cal$slope_variance, achieved_reliability = cal$achieved_marginal_slope_reliability, achieved_partial_reliability = cal$achieved_residualized_slope_reliability, intercept_icc = cal$intercept_icc, marginal_slope_icc = cal$marginal_slope_icc, conditional_slope_icc = cal$conditional_slope_icc, reliability_sd = calibrated$reliability_sd, 
            reliability_iqr = calibrated$reliability_iqr, reliability_min = calibrated$reliability_min, reliability_max = calibrated$reliability_max, reliability_small = reliability[[1L]], reliability_large = reliability[[length(reliability)]], partial_reliability_mean = calibrated$partial_reliability_mean, partial_reliability_min = calibrated$partial_reliability_min, partial_reliability_max = calibrated$partial_reliability_max, partial_reliability_small = partial_reliability[[1L]], partial_reliability_large = partial_reliability[[length(partial_reliability)]], 
            population_lambda22_mean = measurement$lambda22_mean, population_lambda22_min = measurement$lambda22_min, population_lambda22_max = measurement$lambda22_max, population_lambda_matrix_frobenius_rms_dispersion = measurement$lambda_matrix_frobenius_rms_dispersion, population_theta22_mean = measurement$theta22_mean, population_theta22_min = measurement$theta22_min, population_theta22_max = measurement$theta22_max, population_theta_matrix_frobenius_rms_dispersion = measurement$theta_matrix_frobenius_rms_dispersion, 
            population_ols_var22_mean = measurement$ols_var22_mean, population_ols_var22_min = measurement$ols_var22_min, population_ols_var22_max = measurement$ols_var22_max, population_ols_cov_matrix_frobenius_rms_dispersion = measurement$ols_cov_matrix_frobenius_rms_dispersion, slope_variance_marginal = cal$slope_variance, residualized_slope_variance = cal$residualized_slope_variance, slope_intercept_variance_ratio = cal$slope_intercept_variance_ratio, G_condition_number = cal$G_condition_number, 
            sigma = cal$sigma, tau1 = sqrt(cal$slope_variance), rho = condition$marginal_rho[[1]], standardized_beta = structural$standardized_slope_beta, beta1z = structural$beta1_intercept, beta2z = structural$beta2_slope, structural_r2 = structural$total_structural_r_squared, focal_unique_r2 = structural$focal_unique_r_squared, outcome_residual_variance = structural$outcome_residual_variance)
    }) %>% dplyr::bind_rows()
    dplyr::bind_cols(design, calibration_rows) %>% dplyr::mutate(balance_mode = "heterogeneous_profile", min_clus_size = .data$profile_min_clus_size, highly_unbalanced_min_clus_size = .data$profile_min_clus_size, highly_unbalanced_power = NA_real_, r_structure = "iid", r_rho = NA_real_, study_label = "Heterogeneous Cluster Information (v2 fixed shape)", study_structure = "z")
}

study_condition_counts <- function() {
  c(study1v2 = 720L, study2v2 = 720L, study3v2 = 864L, study4v2 = 111L)
}

select_design <- function(study_arg = "allv2", max_conditions = NA_integer_) {
  study_arg <- tolower(trimws(as.character(study_arg[[1]])))
  amended_primary <- paste0("study", 1:4, "v2")
  requested <- switch(study_arg, v2 = amended_primary, amended = amended_primary, allv2 = amended_primary, {
    tokens <- trimws(unlist(strsplit(study_arg, ",", fixed = TRUE)))
    vapply(tokens, function(token) {
      token <- sub("^study", "", token)
      switch(token, `1v2` = "study1v2", `2v2` = "study2v2", `3v2` = "study3v2", `4v2` = "study4v2", NA_character_)
    }, character(1L))
  })
  if (length(requested) == 0L || anyNA(requested)) stop("Only amended v2 Studies 1--4 are included in this release; use `allv2`, `v2`, or `1v2`--`4v2`.")
  builders <- list(study1v2 = make_study1_v2_design, study2v2 = make_study2_v2_design, study3v2 = make_study3_v2_design, study4v2 = make_study4_v2_design)
  counts <- study_condition_counts()
  offsets <- c(0L, cumsum(counts)[-length(counts)])
  names(offsets) <- names(counts)
  built <- lapply(requested, function(study_name) {
    study_design <- builders[[study_name]]()
    if (nrow(study_design) != counts[[study_name]]) stop("Canonical condition count is out of sync for ", study_name)
    study_design %>% dplyr::mutate(condition_id = offsets[[study_name]] + seq_len(dplyr::n())) %>% dplyr::relocate(condition_id)
  })
  design <- dplyr::bind_rows(built)
  if (!is.na(max_conditions) && max_conditions > 0L) design <- design %>% dplyr::slice_head(n = max_conditions)
  if (nrow(design) == 0L) stop("No study conditions selected.")
  design
}
