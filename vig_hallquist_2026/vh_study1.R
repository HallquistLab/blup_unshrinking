#' Vig-Hallquist (2026) simulation study 1: BLUP as outcome

study1_methods <- function() {
  c(
    "oracle",
    "naive_blup",
    "closed_form",
    "fuller_closed_form",
    "fuller_book_preliminary_closed_form",
    "fuller_book_closed_form",
    "fuller_alpha_stepdown_closed_form",
    "single_subject_ols",
    "lai_2spa",
    "direct_mlm"
  )
}

#' Apply the shared VH point/interval eligibility contract to Study 1.
add_study1_analysis_eligibility <- function(results) {
  add_vh_analysis_eligibility(results)
}

simulate_study1 <- function(condition) {
  simulate_data_blup_as_outcome(condition)
}

run_study1_rep <- function(condition) {
  truth <- as.numeric(condition$standardized_beta_target[[1]])
  outcome_scale <- 1 / sqrt(
    as.numeric(condition$slope_variance_marginal[[1]])
  )
  sim <- simulate_study1(condition)
  fit_y <- fit_stage1(
    y ~ x,
    random = ~x | cid,
    data = sim$lv1,
    condition = condition,
    cluster_var = "cid"
  )
  if (is.null(fit_y)) {
    return(
      make_failed_result(condition, study1_methods(), truth) %>%
        add_study1_analysis_eligibility()
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
      make_failed_result(condition, study1_methods(), truth) %>%
        add_study1_analysis_eligibility()
    )
  }

  individual_slopes <- sim$lv1 %>%
    dplyr::group_by(cid_chr) %>%
    dplyr::group_modify(function(.x, .y) {
      fit <- tryCatch(stats::lm(y ~ x, data = .x), error = function(e) NULL)
      tibble::tibble(
        individual_slope = if (is.null(fit)) {
          NA_real_
        } else {
          unname(stats::coef(fit)[["x"]])
        }
      )
    }) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(id = cid_chr, individual_slope)

  stage2_df <- sim$lv2_true %>%
    dplyr::mutate(x = w) %>%
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
    dplyr::left_join(individual_slopes, by = "id") %>%
    dplyr::mutate(zero = 0) %>%
    add_zero_fuller_predictor_outcome_covariance(include_u0 = FALSE)
  stage1_diag <- get_stage1_diagnostics(fit_y, stage2_df)

  direct_data <- sim$lv1 %>%
    dplyr::left_join(
      sim$lv2_true %>% dplyr::select(cid_chr = id, w),
      by = "cid_chr"
    )
  direct_fit <- safe_lmer(
    y ~ w + x + w:x + (1 + x | cid),
    data = direct_data
  )
  direct_row <- if (is.null(direct_fit)) {
    tibble::tibble(
      method = "direct_mlm",
      estimate = NA_real_,
      se = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      status_code = NA_integer_
    )
  } else {
    coef_tab <- summary(direct_fit)$coefficients
    term <- if ("w:x" %in% rownames(coef_tab)) "w:x" else "x:w"
    estimate <- unname(coef_tab[term, "Estimate"])
    se <- unname(coef_tab[term, "Std. Error"])
    tibble::tibble(
      method = "direct_mlm",
      estimate = estimate,
      se = se,
      ci_low = estimate - stats::qnorm(0.975) * se,
      ci_high = estimate + stats::qnorm(0.975) * se,
      status_code = 0L
    )
  }

  results <- dplyr::bind_rows(
    fit_observed_single(
      stage2_df,
      outcome = "true_u1",
      predictor = "w",
      reporting_scale = 1
    ) %>%
      dplyr::filter(se_type == "naive") %>%
      dplyr::transmute(
        method = "oracle", estimate, se, ci_low, ci_high, status_code
      ),
    fit_observed_single(
      stage2_df,
      outcome = "u1_eb",
      predictor = "w",
      reporting_scale = 1
    ) %>%
      dplyr::filter(se_type == "naive") %>%
      dplyr::transmute(
        method = "naive_blup", estimate, se, ci_low, ci_high, status_code
      ),
    fit_observed_single(
      stage2_df,
      outcome = "corrected_slope_full",
      predictor = "w",
      reporting_scale = 1
    ) %>%
      dplyr::filter(se_type == "naive") %>%
      dplyr::transmute(
        method = "closed_form", estimate, se, ci_low, ci_high, status_code
      ),
    fit_fuller_dual_variants(
      stage2_df,
      outcome = "corrected_slope_full",
      predictor_u1 = "w",
      meas22 = "zero",
      outcome_meas_var = "ols_var22",
      predictor_outcome_meas_cov_u1 =
        "fuller_predictor_outcome_meas_cov_u1"
    ) %>%
      rescale_fuller_to_population_sd(1) %>%
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
      outcome = "corrected_slope_full",
      predictor_u0 = NULL,
      predictor_u1 = "w",
      meas11 = NULL,
      meas12 = NULL,
      meas22 = "zero",
      outcome_meas_var = "ols_var22",
      predictor_outcome_meas_cov_u1 =
        "fuller_predictor_outcome_meas_cov_u1"
    ) %>%
      rescale_fuller_to_population_sd(1) %>%
      dplyr::mutate(method = "fuller_alpha_stepdown_closed_form") %>%
      dplyr::select(method, dplyr::everything()),
    fit_observed_single(
      stage2_df,
      outcome = "individual_slope",
      predictor = "w",
      reporting_scale = 1
    ) %>%
      dplyr::filter(se_type == "naive") %>%
      dplyr::transmute(
        method = "single_subject_ols",
        estimate, se, ci_low, ci_high, status_code
      ),
    fit_lai_2spa(stage2_df, use_average = FALSE) %>%
      dplyr::mutate(method = "lai_2spa") %>%
      dplyr::select(method, dplyr::everything()),
    direct_row
  ) %>%
    dplyr::mutate(
      dplyr::across(
        dplyr::any_of(c("estimate", "se", "ci_low", "ci_high")),
        ~ .x * outcome_scale
      )
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
    add_study1_analysis_eligibility()

  dplyr::bind_cols(
    results,
    stage1_diag[rep(1L, nrow(results)), , drop = FALSE]
  ) %>%
    add_study_result_context(condition, truth)
}
