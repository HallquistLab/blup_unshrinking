#' Vig-Hallquist (2026) simulation study 2: BLUP as predictor
#' 
study2_methods <- function() {
  c(
    "oracle_dual",
    "naive_dual_blup",
    "closed_form_dual",
    "fuller_closed_form",
    "fuller_book_preliminary_closed_form",
    "fuller_book_closed_form",
    "fuller_alpha_stepdown_closed_form",
    "lai_2spa",
    "msem"
  )
}

add_study2_method_roles <- function(results) {
  dplyr::mutate(
    results,
    method_role = dplyr::case_when(
      method == "oracle_dual" ~ "benchmark",
      grepl("_slope", method) ~ "slope_only_diagnostic",
      TRUE ~ "primary_dual"
    )
  )
}

#' Apply the shared VH point/interval eligibility contract to Study 2.
#'
#' The explicit point and interval fields distinguish a usable estimate from a
#' usable Wald interval. Legacy `analysis_*` columns remain point-eligibility
#' aliases; raw estimates are never altered.
add_study2_analysis_eligibility <- function(results) {
  add_vh_analysis_eligibility(results)
}

simulate_study2 <- function(condition) {
  simulate_data_blup_as_predictor(condition)
}

run_study2_rep <- function(condition) {
  truth <- as.numeric(condition$standardized_beta_target[[1]])
  latent_slope_sd <- as.numeric(condition$tau1[[1]])
  sim <- simulate_study2(condition)
  fit_y <- fit_stage1(
    y ~ x,
    random = ~x | cid,
    data = sim$lv1,
    condition = condition,
    cluster_var = "cid"
  )
  if (is.null(fit_y)) {
    return(
      make_failed_result(condition, study2_methods(), truth) %>%
        dplyr::mutate(
          analysis_eligible = FALSE,
          analysis_exclusion_reason = "estimation_unavailable"
        ) %>%
        add_study2_analysis_eligibility() %>%
        add_study2_method_roles()
    )
  }

  ordered_ids <- sim$lv2_true$id
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
      make_failed_result(condition, study2_methods(), truth) %>%
        dplyr::mutate(
          analysis_eligible = FALSE,
          analysis_exclusion_reason = "estimation_unavailable"
        ) %>%
        add_study2_analysis_eligibility() %>%
        add_study2_method_roles()
    )
  }

  centered_dat <- dplyr::mutate(
    sim$lv1,
    yc = y - ave(y, cid),
    xc = x - ave(x, cid)
  )
  fit_centered <- fit_stage1(
    yc ~ 0 + xc,
    random = ~0 + xc | cid,
    data = centered_dat,
    condition = condition,
    cluster_var = "cid"
  )
  centered_u1 <- extract_centered_slope_eb(
    fit_centered,
    ordered_ids = ordered_ids
  )

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
    dplyr::left_join(centered_u1, by = "id") %>%
    add_zero_fuller_predictor_outcome_covariance()
  stage1_diag <- get_stage1_diagnostics(fit_y, stage2_df)

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
  )

  results <- results %>%
    add_stage1_estimates(
      fit_obj = fit_y,
      data = sim$lv1,
      cluster_var = "cid",
      within_var = "x",
      stage1_scores = stage2_df,
      true_slope_col = "true_u1"
    ) %>%
    add_study2_analysis_eligibility() %>%
    add_study2_method_roles()

  dplyr::bind_cols(
    results,
    stage1_diag[rep(1L, nrow(results)), , drop = FALSE]
  ) %>%
    add_study_result_context(condition, truth)
}
