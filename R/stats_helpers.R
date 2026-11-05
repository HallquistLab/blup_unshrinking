#' ---
#' title: "Statistical Extraction Helpers"
#' description: "Functions to extract point estimates and standard errors from fitted models."
#' ---

#' Extract statistics from a linear model (lm) object for a specific term.
#'
#' This helper extracts the estimate, standard error, and 95% confidence intervals
#' for a given predictor from an `lm` object. It gracefully handles missing terms
#' by returning NAs instead of throwing an error.
#'
#' @param fit An `lm` object resulting from a linear regression.
#' @param term The string name of the predictor term to extract (e.g., "x", "true_slope_dev").
#' @param use_t Logical; if TRUE, uses the t-distribution for critical values (based on residual df).
#'              If FALSE, uses the normal approximation (z-distribution).
#' @return A 1-row `tibble` containing: `estimate`, `se`, `ci_low`, `ci_high`.
extract_lm_stats <- function(fit, term = "x", use_t = FALSE) {
  coef_tab <- summary(fit)$coefficients
  
  # Return NAs if the requested term is not in the model
  if (!(term %in% rownames(coef_tab))) {
    return(tibble(estimate = NA_real_, se = NA_real_, ci_low = NA_real_, ci_high = NA_real_))
  }
  
  est <- unname(coef_tab[term, "Estimate"])
  se <- unname(coef_tab[term, "Std. Error"])
  
  # Determine the critical value for the 95% CI
  if (use_t) {
    df <- fit$df.residual
    crit <- qt(0.975, df)
  } else {
    crit <- qnorm(0.975)
  }
  
  tibble(
    estimate = est,
    se = se,
    ci_low = est - crit * se,
    ci_high = est + crit * se
  )
}

#' Extract statistics from a linear mixed model (lmer) object for a specific term.
#'
#' This helper extracts the estimate, standard error, and 95% confidence intervals
#' for a given fixed effect from an `lme4::lmer` object. It includes special logic
#' to handle interaction terms whose names might be rearranged by R's formula parser
#' (e.g., matching "x:z" even if the model summary shows "z:x").
#'
#' @param fit An `lmerMod` object.
#' @param term The string name of the predictor to extract. Supports regex for interactions.
#' @param use_t Logical; if TRUE, uses the t-distribution for CIs. Mixed models don't have
#'              an unambiguous degree of freedom, so an approximation is used if `df` is NULL.
#' @param df Explicit degrees of freedom to use if `use_t = TRUE`. If NULL, defaults to an
#'           approximation: `(number of level-2 clusters) - (number of fixed effects)`.
#' @return A 1-row `tibble` containing: `estimate`, `se`, `ci_low`, `ci_high`.
extract_lmer_stats <- function(fit, term = "x:z", use_t = FALSE, df = NULL) {
  coef_tab <- coef(summary(fit))
  
  # 1. Try exact match first
  if (term %in% rownames(coef_tab)) {
    term_name <- term
  } else {
    # 2. Try regex match (useful for cross-level interactions like x:z or z:x)
    term_name <- grep(term, rownames(coef_tab), value = TRUE)
    
    # 3. Fallback specifically for the x:z interaction if R flipped the order to z:x
    if (length(term_name) == 0 && term == "x:z") {
      term_name <- grep("z:x", rownames(coef_tab), value = TRUE)
    }
  }

  # If we couldn't uniquely identify the term, return NAs
  if (length(term_name) != 1L) {
    return(tibble(estimate = NA_real_, se = NA_real_, ci_low = NA_real_, ci_high = NA_real_))
  }

  est <- unname(coef_tab[term_name, "Estimate"])
  se <- unname(coef_tab[term_name, "Std. Error"])
  
  # Determine the critical value for the 95% CI
  if (use_t) {
    if (is.null(df)) {
      # Fallback approximation for degrees of freedom: N_clusters - N_fixed_effects
      # This is conservative but roughly appropriate for level-2 predictors
      n_id <- length(unique(fit@frame$id))
      df <- n_id - length(lme4::fixef(fit))
    }
    crit <- qt(0.975, df)
  } else {
    crit <- qnorm(0.975) # Normal approximation (standard Wald CI)
  }
  
  tibble(
    estimate = est,
    se = se,
    ci_low = est - crit * se,
    ci_high = est + crit * se
  )
}
