#' BLUP and corrected-score helpers.
#'
#' These helpers extract empirical Bayes random-effect predictions from
#' `lme4` models and construct corrected cluster scores used in downstream
#' stage-2 regressions. The main correction implemented here removes the
#' Gaussian prior contribution from the posterior BLUP mean using the posterior
#' and prior precision matrices.

#' Return a fallback value when an object is `NULL`.
#'
#' @details
#' This small infix helper is used for optional arguments where `NULL` means
#' "use the default discovered from the fitted model." Unlike `%in%`-style
#' vectorized helpers, it tests only whether the left-hand side object itself is
#' `NULL`.
#'
#' @param x Candidate value.
#' @param y Fallback value used when `x` is `NULL`.
#'
#' @return
#' `x` when it is not `NULL`; otherwise `y`.
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

#' Convert random-effect names into stable column-name fragments.
#'
#' @details
#' `lme4` random-effect names can contain punctuation such as `(Intercept)` or
#' operators from model formulas. This helper maps those names to simple
#' alphanumeric-plus-underscore fragments so score columns can be named
#' predictably, for example `blup_intercept` and `corrected_slope`.
#'
#' @param x Character vector of raw random-effect names, usually from
#' `colnames(lme4::ranef(fit)[[group]])`.
#'
#' @return
#' A character vector with `(Intercept)` mapped to `intercept`, non-alphanumeric
#' runs replaced by underscores, leading/trailing underscores stripped, and any
#' empty result replaced by `"re"`.
sanitize_re_name <- function(x) {
  # Preserve the conventional intercept label before removing punctuation.
  x <- gsub("^\\(Intercept\\)$", "intercept", x)
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  ifelse(nzchar(x), x, "re")
}

#' Warn when a legacy iid-only lme4 score extractor is used.
#'
#' @param fun_name Character scalar naming the caller.
#'
#' @return Invisibly returns `NULL`.
warn_legacy_iid_score_extractor <- function(fun_name) {
  warning(
    fun_name,
    "() is a legacy iid-only lme4 workflow. Prefer get_stage1_eb_components(), ",
    "which computes EB/BLUPs, posterior covariance, diagonal and matrix ",
    "corrections, and Lai measurement inputs from one consistent Stage-1 fit. ",
    "That is required for non-diagonal residual covariance because lme4::ranef(..., condVar = TRUE) ",
    "and lme4::VarCorr() are tied to lme4's iid-residual likelihood.",
    call. = FALSE
  )
  invisible(NULL)
}

#' Normalize a cluster-level residual covariance list.
#'
#' @param R_list Optional list of residual covariance matrices.
#' @param cluster_ids Character vector of cluster identifiers in data order.
#'
#' @return
#' `NULL` or a named `R_list` ordered to `cluster_ids`.
normalize_R_list <- function(R_list, cluster_ids) {
  if (is.null(R_list)) {
    return(NULL)
  }
  if (!is.list(R_list)) {
    stop("`R_list` must be NULL or a list of cluster-level residual covariance matrices.")
  }
  if (is.null(names(R_list))) {
    if (length(R_list) != length(cluster_ids)) {
      stop("Unnamed `R_list` must have one matrix per cluster in `data` order.")
    }
    names(R_list) <- cluster_ids
  }
  R_list <- R_list[cluster_ids]
  if (any(vapply(R_list, is.null, logical(1)))) {
    stop("`R_list` does not contain one matrix per cluster in `data`.")
  }
  R_list
}

#' Coerce an nlme VarCov object to a plain numeric matrix.
#'
#' @param x Object returned by `nlme::getVarCov()`.
#'
#' @return Numeric matrix.
as_plain_vcov_matrix <- function(x) {
  x_unclass <- unclass(x)
  if (is.list(x_unclass)) {
    x_unclass <- x_unclass[[1L]]
  }
  as.matrix(x_unclass)
}

#' Extract common Stage-1 ingredients for GLS-aware score correction.
#'
#' @details
#' `get_gls_corrected_scores()` needs only fixed effects, the random-effect
#' covariance `G`, cluster-level residual covariance matrices `R_i`, and stable
#' random-effect names. This adapter keeps those ingredients independent of the
#' package used to fit Stage 1.
#'
#' `lme4` objects are accepted only for iid residual covariance. Passing an
#' explicit non-NULL `R_list` with a `merMod` object is refused because the
#' `lme4` fixed-effect and random-effect covariance estimates come from the
#' iid-residual likelihood, not from the supplied `R_i` likelihood.
#'
#' `nlme::lme` objects can carry fitted R-side correlation structures. When
#' `R_list` is `NULL`, their fitted conditional residual covariance matrices
#' are extracted with `nlme::getVarCov(..., type = "conditional")`. A supplied
#' `R_list` overrides the fitted residual covariance list, which is useful for
#' known-`R_i` simulations.
#'
#' @param fit_obj Fitted Stage-1 model object.
#' @param data Original long-format data.
#' @param cluster_var Character scalar naming the grouping column.
#' @param within_var Optional within-cluster random-slope predictor.
#' @param R_list Optional cluster residual covariance list.
#' @param group Optional grouping-factor name.
#'
#' @return
#' A list with `beta_hat`, `G_hat`, `R_list`, and `re_names_raw`.
extract_stage1_components <- function(fit_obj, data, cluster_var, within_var = NULL,
                                      R_list = NULL, group = NULL) {
  UseMethod("extract_stage1_components")
}

#' @export
extract_stage1_components.merMod <- function(fit_obj, data, cluster_var, within_var = NULL,
                                             R_list = NULL, group = NULL) {
  cluster_ids <- unique(as.character(data[[cluster_var]]))
  if (!is.null(R_list)) {
    stop(
      "`R_list` cannot be supplied with an `lme4` merMod object. ",
      "`lme4` estimates beta and G under iid residual covariance; use an ",
      "`nlme::lme` fit or another R-aware Stage-1 fit instead."
    )
  }

  beta_hat <- lme4::fixef(fit_obj)
  re_list <- lme4::ranef(fit_obj)
  group <- group %||% names(re_list)[[1]]
  g_hat <- as.matrix(lme4::VarCorr(fit_obj)[[group]])
  sigma2_hat <- stats::sigma(fit_obj)^2
  split_dat <- split(data, as.character(data[[cluster_var]]), drop = TRUE)[cluster_ids]
  R_iid <- stats::setNames(
    lapply(split_dat, function(df_i) sigma2_hat * diag(nrow(df_i))),
    cluster_ids
  )

  list(
    beta_hat = beta_hat,
    G_hat = g_hat,
    R_list = R_iid,
    re_names_raw = colnames(re_list[[group]])
  )
}

#' @export
extract_stage1_components.lme <- function(fit_obj, data, cluster_var, within_var = NULL,
                                          R_list = NULL, group = NULL) {
  if (!requireNamespace("nlme", quietly = TRUE)) {
    stop("The `nlme` package is required to extract Stage-1 components from `lme` objects.")
  }

  cluster_ids <- unique(as.character(data[[cluster_var]]))
  beta_hat <- nlme::fixef(fit_obj)
  g_hat <- as_plain_vcov_matrix(nlme::getVarCov(fit_obj, type = "random.effects"))

  if (is.null(R_list)) {
    R_list <- stats::setNames(lapply(cluster_ids, function(cluster_id) {
      as_plain_vcov_matrix(nlme::getVarCov(
        fit_obj,
        individuals = cluster_id,
        type = "conditional"
      ))
    }), cluster_ids)
  } else {
    R_list <- normalize_R_list(R_list, cluster_ids)
  }

  re_names_raw <- colnames(g_hat)
  if (is.null(re_names_raw)) {
    re_names_raw <- if (is.null(within_var)) "(Intercept)" else c("(Intercept)", within_var)
  }

  list(
    beta_hat = beta_hat,
    G_hat = g_hat,
    R_list = R_list,
    re_names_raw = re_names_raw
  )
}

#' @export
extract_stage1_components.default <- function(fit_obj, data, cluster_var, within_var = NULL,
                                             R_list = NULL, group = NULL) {
  stop(
    "Unsupported Stage-1 fit object class for GLS-corrected scores: ",
    paste(class(fit_obj), collapse = ", "),
    ". Supported classes are `lme4` merMod and `nlme` lme."
  )
}

#' Extract fixed effects from a supported Stage-1 model.
#'
#' @param fit_obj Fitted Stage-1 model.
#'
#' @return Named numeric fixed-effect vector.
stage1_fixef <- function(fit_obj) {
  if (inherits(fit_obj, "lme")) {
    if (!requireNamespace("nlme", quietly = TRUE)) {
      stop("The `nlme` package is required to extract fixed effects from `lme` objects.")
    }
    return(nlme::fixef(fit_obj))
  }
  lme4::fixef(fit_obj)
}

#' Format manual EB/posterior components as a Stage-2 row.
#'
#' @param cluster_id Cluster identifier.
#' @param eb EB/BLUP vector.
#' @param post_vcov Posterior covariance matrix.
#' @param lambda Lai reliability/loading matrix.
#' @param theta Lai EB measurement residual covariance matrix.
#' @param mle Direct MLE/Bartlett score vector.
#' @param mle_vcov MLE/Bartlett score covariance matrix.
#' @param corrected Matrix-corrected MLE/Bartlett score vector.
#' @param corrected_vars Diagonal of the matrix-corrected score covariance.
#' @param diag_corrected Diagonal-only corrected score vector.
#' @param diag_vars Diagonal-only corrected score variances.
#' @param re_names Sanitized random-effect names.
#'
#' @return One-row tibble with BLUP, posterior, Lai, and corrected-score columns.
format_stage1_eb_row <- function(cluster_id, eb, post_vcov, lambda, theta,
                                 mle, mle_vcov, corrected, corrected_vars,
                                 diag_corrected, diag_vars, re_names) {
  n_re <- length(eb)
  out <- tibble::tibble(id = cluster_id)

  for (j in seq_len(n_re)) {
    out[[paste0("blup_", re_names[[j]])]] <- eb[[j]]
    out[[paste0("mle_", re_names[[j]])]] <- mle[[j]]
    out[[paste0("mle_", re_names[[j]], "_var")]] <- mle_vcov[j, j]
    out[[paste0("corrected_", re_names[[j]])]] <- corrected[[j]]
    out[[paste0("corrected_", re_names[[j]], "_var")]] <- corrected_vars[[j]]
    out[[paste0("corrected_", re_names[[j]], "_diag")]] <- diag_corrected[[j]]
    out[[paste0("corrected_", re_names[[j]], "_diag_var")]] <- diag_vars[[j]]
  }

  if (n_re == 1L) {
    out$u0_eb <- eb[[1]]
    out$postvar11 <- post_vcov[1, 1]
    out$lambda11 <- lambda[1, 1]
    out$theta11 <- theta[1, 1]
    out$mle_var11 <- mle_vcov[1, 1]
  } else if (n_re == 2L) {
    out$u0_eb <- eb[[1]]
    out$u1_eb <- eb[[2]]
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
    out$mle_var11 <- mle_vcov[1, 1]
    out$mle_var12 <- mle_vcov[1, 2]
    out$mle_var22 <- mle_vcov[2, 2]
  } else {
    stop("Only univariate and bivariate random-effect components are currently supported.")
  }

  out
}

#' Compute R-aware EB means, posterior covariance, and score corrections.
#'
#' @details
#' This is the canonical Stage-1 ingredient extractor for downstream BLUP,
#' corrected-score, Fuller, and Lai/OpenMx estimators. It intentionally computes
#' EB means and posterior covariance from the Gaussian conditioning equations
#' rather than relying on package-specific random-effect extractors:
#'
#' `Sigma_i = Z_i G Z_i' + R_i`
#'
#' `A_i = G Z_i' Sigma_i^{-1}`
#'
#' `b_Ei = A_i (y_i - X_i beta)`
#'
#' `V_post_i = (G^{-1} + Z_i' R_i^{-1} Z_i)^{-1}`
#'
#' The same pieces also give Lai's reliability and unreliability matrices:
#'
#' `lambda_i = A_i Z_i`
#'
#' `theta_i = A_i R_i A_i'`
#'
#' For `merMod` fits, `R_list` must be `NULL`; the adapter supplies the fitted
#' iid residual covariance `sigma^2 I`. For `nlme::lme` fits, `R_list = NULL`
#' uses the fitted conditional residual covariance matrices from the object.
#'
#' @param fit_obj Fitted Stage-1 model supported by `extract_stage1_components()`.
#' @param data Original long-format data frame.
#' @param cluster_var Character scalar naming the cluster/grouping column.
#' @param outcome_var Character scalar naming the outcome column.
#' @param within_var Optional character scalar naming the random-slope
#' predictor. If `NULL`, an intercept-only random-effect design is used.
#' @param R_list Optional named list of cluster-level residual covariance
#' matrices.
#' @param group Optional grouping-factor name for `merMod` objects.
#'
#' @return A tibble with one row per cluster containing EB/BLUP columns
#' (`u0_eb`, `u1_eb`, `blup_*`), posterior variances (`postvar*`), Lai
#' measurement matrices (`lambda*`, `theta*`), direct MLE/Bartlett scores
#' (`mle_*`), matrix-corrected scores (`corrected_*`), and diagonal-only
#' corrected scores (`corrected_*_diag`).
get_stage1_eb_components <- function(fit_obj, data, cluster_var, outcome_var, within_var = NULL,
                                     R_list = NULL, group = NULL) {
  cluster_ids <- unique(as.character(data[[cluster_var]]))
  split_dat <- split(data, as.character(data[[cluster_var]]), drop = TRUE)[cluster_ids]
  n_re <- if (is.null(within_var)) 1L else 2L

  stage1 <- extract_stage1_components(
    fit_obj = fit_obj,
    data = data,
    cluster_var = cluster_var,
    within_var = within_var,
    R_list = R_list,
    group = group
  )
  beta_hat <- stage1$beta_hat
  g_hat <- stage1$G_hat
  R_list <- normalize_R_list(stage1$R_list, cluster_ids)

  re_names_raw <- stage1$re_names_raw
  if (is.null(re_names_raw) || length(re_names_raw) != n_re) {
    re_names_raw <- if (is.null(within_var)) "(Intercept)" else c("(Intercept)", within_var)
  }
  re_names <- sanitize_re_name(re_names_raw)

  if (nrow(g_hat) != n_re || ncol(g_hat) != n_re) {
    stop("The fitted random-effect covariance dimension does not match `within_var`.")
  }
  purrr::map_dfr(cluster_ids, function(cluster_id) {
    df_i <- split_dat[[cluster_id]]
    if (is.null(within_var)) {
      z_mat <- matrix(1, nrow = nrow(df_i), ncol = 1L)
      x_mat <- z_mat
      beta_vec <- beta_hat[[1]]
    } else {
      z_vec <- df_i[[within_var]]
      z_mat <- cbind(1, z_vec)
      x_mat <- z_mat
      beta_vec <- c(beta_hat[[1]], beta_hat[[within_var]])
    }

    R_i <- as.matrix(R_list[[cluster_id]])
    pieces <- tryCatch({
      if (!is.matrix(R_i) || nrow(R_i) != nrow(df_i) || ncol(R_i) != nrow(df_i)) {
        stop("R_i has incompatible dimensions.")
      }

      resid_i <- df_i[[outcome_var]] - as.numeric(x_mat %*% beta_vec)
      sigma_y_i <- z_mat %*% g_hat %*% t(z_mat) + R_i
      sigma_y_inv <- solve(sigma_y_i)
      a_i <- g_hat %*% t(z_mat) %*% sigma_y_inv
      R_inv_Z <- solve(R_i, z_mat)
      R_inv_resid <- solve(R_i, resid_i)
      info_like <- crossprod(z_mat, R_inv_Z)
      # Use the covariance identity rather than (G^{-1} + Z'R^{-1}Z)^{-1}
      # so singular but positive-semidefinite Stage-1 random-effect
      # covariance estimates still produce EB means and Lai measurement inputs.
      post_vcov <- g_hat - a_i %*% z_mat %*% g_hat
      post_vcov <- (post_vcov + t(post_vcov)) / 2
      mle_vcov <- solve(info_like)
      eb <- as.numeric(a_i %*% resid_i)
      mle <- as.numeric(mle_vcov %*% crossprod(z_mat, R_inv_resid))
      corrected <- unweight_random_effects(
        post_mean = eb,
        post_vcov = post_vcov,
        prior_mean = rep(0, n_re),
        prior_vcov = g_hat,
        return_var = TRUE
      )

      diag_corrected <- diag_vars <- rep(NA_real_, n_re)
      for (j in seq_len(n_re)) {
        denom <- (1 / post_vcov[j, j]) - (1 / g_hat[j, j])
        if (is.finite(denom) && denom > sqrt(.Machine$double.eps)) {
          diag_vars[[j]] <- 1 / denom
          diag_corrected[[j]] <- diag_vars[[j]] * eb[[j]] / post_vcov[j, j]
        }
      }

      list(
        eb = eb,
        post_vcov = post_vcov,
        lambda = a_i %*% z_mat,
        theta = a_i %*% R_i %*% t(a_i),
        mle = mle,
        mle_vcov = mle_vcov,
        corrected = corrected$scores,
        corrected_vars = corrected$vars,
        diag_corrected = diag_corrected,
        diag_vars = diag_vars
      )
    }, error = function(e) {
      list(
        eb = rep(NA_real_, n_re),
        post_vcov = matrix(NA_real_, nrow = n_re, ncol = n_re),
        lambda = matrix(NA_real_, nrow = n_re, ncol = n_re),
        theta = matrix(NA_real_, nrow = n_re, ncol = n_re),
        mle = rep(NA_real_, n_re),
        mle_vcov = matrix(NA_real_, nrow = n_re, ncol = n_re),
        corrected = rep(NA_real_, n_re),
        corrected_vars = rep(NA_real_, n_re),
        diag_corrected = rep(NA_real_, n_re),
        diag_vars = rep(NA_real_, n_re)
      )
    })

    format_stage1_eb_row(
      cluster_id = cluster_id,
      eb = pieces$eb,
      post_vcov = pieces$post_vcov,
      lambda = pieces$lambda,
      theta = pieces$theta,
      mle = pieces$mle,
      mle_vcov = pieces$mle_vcov,
      corrected = pieces$corrected,
      corrected_vars = pieces$corrected_vars,
      diag_corrected = pieces$diag_corrected,
      diag_vars = pieces$diag_vars,
      re_names = re_names
    )
  })
}

#' Unweight random effects from the prior using matrix inversion.
#'
#' @details
#' In the Gaussian mixed-model setup, each empirical Bayes random-effect
#' prediction is a posterior mean that combines the cluster likelihood with the
#' population random-effect prior. If
#' `V_post` is the conditional posterior covariance, `V_prior` is the
#' population random-effect covariance, `m_post` is the BLUP/posterior mean, and
#' `m_prior` is the prior mean, this helper solves
#'
#' `m_like = (V_post^{-1} - V_prior^{-1})^{-1}
#'           (V_post^{-1} m_post - V_prior^{-1} m_prior)`.
#'
#' The result is a likelihood-only or "unweighted" score estimate under this
#' normal-normal algebra. Singular or numerically invalid systems return an
#' all-`NA_real_` vector of the same length as `post_mean`, allowing simulation
#' loops to continue without special error handling.
#'
#' @param post_mean Numeric vector of posterior means, typically one row of
#' `lme4::ranef(fit, condVar = TRUE)[[group]]`.
#' @param post_vcov Numeric posterior covariance matrix for the same cluster.
#' @param prior_mean Numeric prior mean vector for the random effects, usually
#' zeros in `lme4` fits.
#' @param prior_vcov Numeric prior covariance matrix from
#' `lme4::VarCorr(fit)[[group]]`.
#' @param return_var Logical scalar indicating whether to return the variance of
#' the unweighted score. If `TRUE`, the function returns a list with `scores
#' = corrected score` and `var = (V_post^{-1} - V_prior^{-1})^{-1}`. If `FALSE`,
#' only the score is returned.
#'
#' @return
#' Numeric vector of corrected likelihood-only random-effect scores, or an
#' `NA_real_` vector if the precision subtraction or solve fails.
unweight_random_effects <- function(post_mean, post_vcov, prior_mean, prior_vcov, return_var = FALSE) {
  out <- tryCatch({
    # Subtract the prior precision contribution from the posterior precision
    # and solve for the score implied by the cluster likelihood alone.
    solve(
      solve(post_vcov) - solve(prior_vcov),
      solve(post_vcov, post_mean) - solve(prior_vcov, prior_mean)
    )
  }, error = function(e) {
    rep(NA_real_, length(post_mean))
  })
  scores <- as.numeric(out)

  if (return_var) {
    out_var <- tryCatch({
      solve(solve(post_vcov) - solve(prior_vcov))
    }, error = function(e) {
      matrix(NA_real_, nrow = length(post_mean), ncol = length(post_mean))
    })
    return(list(scores = scores, vars = diag(out_var)))
  } else {
    return(scores)
  }
}

#' Extract legacy iid-only EB/BLUPs and Vig-style corrected scores.
#'
#' @details
#' This helper is retained for legacy iid `lme4` workflows. New simulation code
#' should prefer `get_stage1_eb_components()`, which supports both iid `merMod`
#' fits and R-aware `nlme::lme` fits through the Stage-1 adapter layer. The
#' newer helper also returns EB/BLUPs, posterior covariance, diagonal and matrix
#' corrections, and Lai measurement inputs from one consistent Stage-1 fit.
#'
#' Calls `lme4::ranef(fit_null, condVar = TRUE)` to obtain cluster-level BLUPs
#' and their conditional posterior covariance matrices, then applies
#' `unweight_random_effects()` cluster by cluster. The returned tibble uses
#' model-derived columns of the form `blup_<effect>` and `corrected_<effect>`,
#' where `<effect>` is sanitized by `sanitize_re_name()`. For example,
#' `(Intercept)` becomes `blup_intercept` / `corrected_intercept`, while a
#' random slope for `z` becomes `blup_z` / `corrected_z`.
#'
#' Because `lme4` assumes iid level-1 residual covariance, this function should
#' not be used when data are generated or modeled with non-diagonal residual
#' covariance such as AR(1) or Toeplitz structures.
#'
#' @param fit_null Fitted `lme4` mixed model with random effects available via
#' `ranef(..., condVar = TRUE)`.
#' @param group Optional character scalar naming the grouping factor to extract.
#' If `NULL`, the first grouping factor in the `ranef()` list is used.
#' @return
#' A tibble with one row per group level. The `id` column contains group-level
#' row names from `ranef()`, followed by BLUP and corrected-score columns.
get_corrected_scores <- function(fit_null, group = NULL) {
  warn_legacy_iid_score_extractor("get_corrected_scores")

  re_list <- lme4::ranef(fit_null, condVar = TRUE)
  group <- group %||% names(re_list)[[1]]
  re_df <- re_list[[group]]

  # `postVar` is a p x p x n array of conditional posterior covariance
  # matrices aligned with the rows of the random-effects data frame.
  post_var_arr <- attr(re_df, "postVar")
  prior_vcov <- as.matrix(lme4::VarCorr(fit_null)[[group]])
  n_re <- ncol(re_df)
  prior_mean <- rep(0, n_re)
  re_names <- sanitize_re_name(colnames(re_df))

  purrr::map_dfr(seq_len(nrow(re_df)), function(i) {
    post_mean <- as.numeric(re_df[i, ])
    post_vcov <- post_var_arr[, , i]

    # For one-dimensional random-effect structures, subsetting can drop the
    # posterior covariance to a scalar; restore matrix shape for the algebra.
    if (!is.matrix(post_vcov)) {
      post_vcov <- matrix(post_vcov, nrow = 1L, ncol = 1L)
    }
    corrected <- unweight_random_effects(
      post_mean = post_mean,
      post_vcov = post_vcov,
      prior_mean = prior_mean,
      prior_vcov = prior_vcov,
      return_var = TRUE
    )

    out <- tibble::tibble(id = rownames(re_df)[[i]])

    # Always emit formula-derived names so this helper also works for models
    # whose random effects are not simply intercept/slope.
    for (j in seq_len(n_re)) {
      out[[paste0("blup_", re_names[[j]])]] <- post_mean[[j]]
      out[[paste0("corrected_", re_names[[j]])]] <- corrected$scores[[j]]
      out[[paste0("corrected_", re_names[[j]], "_var")]] <- corrected$vars[[j]]
    }

    out
  })
}

#' Extract legacy iid-only diagonal corrected random-effect scores.
#'
#' @details
#' This helper is retained for legacy iid `lme4` workflows. New simulation code
#' should prefer `get_stage1_eb_components()` and use its
#' `corrected_<effect>_diag` and `corrected_<effect>_diag_var` columns. The
#' newer helper computes diagonal-only corrections from the same Stage-1
#' ingredients used for EB/BLUPs, posterior covariance, full matrix correction,
#' and Lai measurement inputs.
#'
#' This helper is intentionally separate from `get_corrected_scores()` because
#' the diagonal-only correction is a diagnostic comparison, not the full
#' corrected-score estimator. It applies the scalar prior-unweighting formula to
#' each random-effect component independently and ignores posterior/prior
#' covariance among random effects.
#'
#' For a component `j`, the returned value is
#' `(V_post[j,j]^{-1} - V_prior[j,j]^{-1})^{-1}
#'  V_post[j,j]^{-1} m_post[j]`.
#'
#' Because `lme4` assumes iid level-1 residual covariance, this function should
#' not be used when data are generated or modeled with non-diagonal residual
#' covariance such as AR(1) or Toeplitz structures.
#'
#' @param fit_null Fitted `lme4` mixed model with random effects available via
#' `ranef(..., condVar = TRUE)`.
#' @param group Optional character scalar naming the grouping factor to extract.
#' If `NULL`, the first grouping factor in the `ranef()` list is used.
#'
#' @return
#' A tibble with one row per group level, an `id` column, and columns named
#' `corrected_<effect>_diag` for each random-effect component.
get_diagonal_corrected_scores <- function(fit_null, group = NULL) {
  warn_legacy_iid_score_extractor("get_diagonal_corrected_scores")

  re_list <- lme4::ranef(fit_null, condVar = TRUE)
  group <- group %||% names(re_list)[[1]]
  re_df <- re_list[[group]]
  post_var_arr <- attr(re_df, "postVar")
  prior_vcov <- as.matrix(lme4::VarCorr(fit_null)[[group]])
  n_re <- ncol(re_df)
  re_names <- sanitize_re_name(colnames(re_df))

  purrr::map_dfr(seq_len(nrow(re_df)), function(i) {
    post_mean <- as.numeric(re_df[i, ])
    post_vcov <- post_var_arr[, , i]
    if (!is.matrix(post_vcov)) {
      post_vcov <- matrix(post_vcov, nrow = 1L, ncol = 1L)
    }

    out <- tibble::tibble(id = rownames(re_df)[[i]])
    for (j in seq_len(n_re)) {
      diag_corrected <- tryCatch({
        denom <- (1 / post_vcov[j, j]) - (1 / prior_vcov[j, j])
        if (!is.finite(denom) || denom <= sqrt(.Machine$double.eps)) {
          NA_real_
        } else {
          (1 / denom) * (post_mean[[j]] / post_vcov[j, j])
        }
      }, error = function(e) {
        NA_real_
      })
      diag_var <- tryCatch({
        denom <- (1 / post_vcov[j, j]) - (1 / prior_vcov[j, j])
        if (!is.finite(denom) || denom <= sqrt(.Machine$double.eps)) {
          NA_real_
        } else {
          1 / denom
        }
      }, error = function(e) {
        NA_real_
      })

      out[[paste0("corrected_", re_names[[j]], "_diag")]] <- diag_corrected
      out[[paste0("corrected_", re_names[[j]], "_diag_var")]] <- diag_var
    }

    out
  })
}

#' Extract GLS-aware EB/BLUPs, posterior covariance, and MLE scores.
#'
#' @details
#' This helper is the residual-covariance-aware counterpart to
#' `get_corrected_scores()`. It does not use `lme4::ranef(..., condVar = TRUE)`
#' for the cluster predictions, because those conditional modes and posterior
#' covariance matrices assume the residual covariance structure fitted by
#' the original model object. Instead, it rebuilds the Gaussian score
#' ingredients directly from common Stage-1 components extracted by
#' `extract_stage1_components()`.
#'
#' For cluster `i`, let `r_i = y_i - X_i beta_hat`, `Z_i` be the random-effect
#' design, `G` be the fitted random-effect covariance, `R_i` be the residual
#' covariance, and `Sigma_i = Z_i G Z_i' + R_i` be the marginal covariance.
#' The EB/BLUP scoring matrix in Vig and Hallquist (V&H) equation 29 is
#'
#' `A_Ei = G Z_i' Sigma_i^{-1}`,
#'
#' giving the EB/BLUP prediction
#'
#' `b_Ei = A_Ei r_i = G Z_i' Sigma_i^{-1} r_i`.
#'
#' The conditional posterior covariance of `b_i | y_i` is
#'
#' `V_post_i = (G^{-1} + Z_i' R_i^{-1} Z_i)^{-1}`.
#'
#' The MLE/Bartlett scoring matrix in V&H equation 30 is
#'
#' `A_Mi = (Z_i' R_i^{-1} Z_i)^{-1} Z_i' R_i^{-1}`,
#'
#' giving the likelihood-only score
#'
#' `b_Mi = A_Mi r_i`.
#'
#' The V&H equation 33 gives the same MLE score as a conversion from
#' the EB/BLUP score:
#'
#' `b_Mi = ((G Z_i' R_i^{-1} Z_i)^{-1} + I) b_Ei`.
#'
#' This implementation uses the equivalent posterior-minus-prior precision
#' identity already used by `unweight_random_effects()`:
#'
#' `b_Mi = (V_post_i^{-1} - G^{-1})^{-1} V_post_i^{-1} b_Ei`.
#'
#' The returned `corrected_*` columns are these matrix-converted MLE/Bartlett
#' scores. The returned `mle_*` columns are the direct GLS scores from equation
#' 30, included for diagnostics; the two should agree up to numerical
#' tolerance when `G` and `R_i` are nonsingular. The `corrected_*_var` columns
#' are the diagonal elements of `(Z_i' R_i^{-1} Z_i)^{-1}`, matching the
#' conditional variance expression in V&H equation 38.
#'
#' @param fit_obj Fitted Stage-1 model. Supported classes are `lme4` `merMod`
#'   objects for iid residual covariance and `nlme` `lme` objects for fitted
#'   R-side residual covariance structures. Supplying non-NULL `R_list` with a
#'   `merMod` object is an error because `lme4` estimates `beta_hat` and `G_hat`
#'   under an iid-residual likelihood.
#' @param data Original long-format data frame.
#' @param cluster_var Character scalar naming the cluster/grouping column.
#' @param outcome_var Character scalar naming the outcome column.
#' @param within_var Optional character scalar naming the random-slope
#'   predictor. If `NULL`, an intercept-only random-effect design is used.
#' @param R_list Optional named list of cluster-level residual covariance
#'   matrices. Names should match `cluster_var` values. If unnamed, matrices are
#'   matched to clusters in data order. For `nlme` fits, `NULL` means use the
#'   fitted conditional residual covariance matrices from the object.
#' @param group Optional grouping-factor name for extracting `G_hat` and
#'   random-effect names from `merMod` objects. Defaults to the first grouping
#'   factor.
#'
#' @return
#' A tibble with one row per cluster. It includes EB/BLUP columns
#' `blup_<effect>`, posterior covariance columns `postvar*`, direct
#' MLE/Bartlett columns `mle_<effect>`, matrix-converted score columns
#' `corrected_<effect>`, and measurement-variance columns
#' `corrected_<effect>_var`.
get_gls_corrected_scores <- function(fit_obj, data, cluster_var, outcome_var, within_var = NULL,
                                     R_list = NULL, group = NULL) {
  cluster_ids <- unique(as.character(data[[cluster_var]]))
  split_dat <- split(data, as.character(data[[cluster_var]]), drop = TRUE)[cluster_ids]
  n_re <- if (is.null(within_var)) 1L else 2L

  stage1 <- extract_stage1_components(
    fit_obj = fit_obj,
    data = data,
    cluster_var = cluster_var,
    within_var = within_var,
    R_list = R_list,
    group = group
  )
  beta_hat <- stage1$beta_hat
  g_hat <- stage1$G_hat
  R_list <- normalize_R_list(stage1$R_list, cluster_ids)

  re_names_raw <- stage1$re_names_raw
  if (is.null(re_names_raw) || length(re_names_raw) != n_re) {
    re_names_raw <- if (is.null(within_var)) "(Intercept)" else c("(Intercept)", within_var)
  }
  re_names <- sanitize_re_name(re_names_raw)

  if (nrow(g_hat) != n_re || ncol(g_hat) != n_re) {
    stop("The fitted random-effect covariance dimension does not match `within_var`.")
  }

  prior_mean <- rep(0, n_re)
  g_inv <- tryCatch(solve(g_hat), error = function(e) NULL)

  purrr::map_dfr(cluster_ids, function(cluster_id) {
    df_i <- split_dat[[cluster_id]]

    if (is.null(within_var)) {
      z_mat <- matrix(1, nrow = nrow(df_i), ncol = 1L)
      x_mat <- z_mat
      beta_vec <- beta_hat[[1]]
    } else {
      z_vec <- df_i[[within_var]]
      z_mat <- cbind(1, z_vec)
      x_mat <- z_mat
      beta_vec <- c(beta_hat[[1]], beta_hat[[within_var]])
    }

    resid_i <- df_i[[outcome_var]] - as.numeric(x_mat %*% beta_vec)
    R_i <- as.matrix(R_list[[cluster_id]])

    pieces <- tryCatch({
      if (is.null(g_inv)) {
        stop("G is singular.")
      }
      if (!is.matrix(R_i) || nrow(R_i) != nrow(df_i) || ncol(R_i) != nrow(df_i)) {
        stop("R_i has incompatible dimensions.")
      }

      # EB/BLUP from V&H equations 29 and 31:
      # b_Ei = G Z_i' Sigma_i^-1 r_i, Sigma_i = Z_i G Z_i' + R_i.
      sigma_y_i <- z_mat %*% g_hat %*% t(z_mat) + R_i
      sigma_inv_resid <- solve(sigma_y_i, resid_i)
      eb <- as.numeric(g_hat %*% crossprod(z_mat, sigma_inv_resid))

      # Posterior covariance and MLE/Bartlett information. Avoid explicitly
      # forming R_i^-1; solve(R_i, Z_i) and solve(R_i, r_i) provide the needed
      # products for Z_i' R_i^-1 Z_i and Z_i' R_i^-1 r_i.
      R_inv_Z <- solve(R_i, z_mat)
      R_inv_resid <- solve(R_i, resid_i)
      info_like <- crossprod(z_mat, R_inv_Z)
      post_vcov <- solve(g_inv + info_like)

      # Direct MLE/Bartlett score from V&H equation 30.
      mle_vcov <- solve(info_like)
      mle <- as.numeric(mle_vcov %*% crossprod(z_mat, R_inv_resid))

      # Matrix conversion from EB/BLUP to MLE/Bartlett. This is the precision subtraction form of V&H equation 33.
      corrected <- unweight_random_effects(
        post_mean = eb,
        post_vcov = post_vcov,
        prior_mean = prior_mean,
        prior_vcov = g_hat,
        return_var = TRUE
      )

      list(
        eb = eb,
        post_vcov = post_vcov,
        mle = mle,
        mle_vcov = mle_vcov,
        corrected = corrected$scores,
        corrected_vars = corrected$vars
      )
    }, error = function(e) {
      list(
        eb = rep(NA_real_, n_re),
        post_vcov = matrix(NA_real_, nrow = n_re, ncol = n_re),
        mle = rep(NA_real_, n_re),
        mle_vcov = matrix(NA_real_, nrow = n_re, ncol = n_re),
        corrected = rep(NA_real_, n_re),
        corrected_vars = rep(NA_real_, n_re)
      )
    })

    out <- tibble::tibble(id = cluster_id)
    for (j in seq_len(n_re)) {
      out[[paste0("blup_", re_names[[j]])]] <- pieces$eb[[j]]
      out[[paste0("mle_", re_names[[j]])]] <- pieces$mle[[j]]
      out[[paste0("corrected_", re_names[[j]])]] <- pieces$corrected[[j]]
      out[[paste0("corrected_", re_names[[j]], "_var")]] <- pieces$corrected_vars[[j]]
      out[[paste0("mle_", re_names[[j]], "_var")]] <- pieces$mle_vcov[j, j]
    }

    if (n_re == 1L) {
      out$postvar11 <- pieces$post_vcov[1, 1]
      out$mle_var11 <- pieces$mle_vcov[1, 1]
    } else {
      out$postvar11 <- pieces$post_vcov[1, 1]
      out$postvar12 <- pieces$post_vcov[1, 2]
      out$postvar22 <- pieces$post_vcov[2, 2]
      out$mle_var11 <- pieces$mle_vcov[1, 1]
      out$mle_var12 <- pieces$mle_vcov[1, 2]
      out$mle_var22 <- pieces$mle_vcov[2, 2]
    }

    out
  })
}

#' Compute likelihood-only cluster scores directly from residualized cluster GLS.
#'
#' @details
#' This helper computes likelihood-only random-effect scores without using
#' posterior BLUP covariance matrices. For cluster `i`, let `r_i = y_i - X_i
#' beta` be the fixed-effect residual, `Z_i` be the random-effect design, and
#' `R_i` be the level-1 residual covariance matrix. The MLE/Bartlett score is
#' the GLS coefficient from regressing `r_i` on `Z_i`:
#'
#' `b_Mi = (Z_i' R_i^{-1} Z_i)^{-1} Z_i' R_i^{-1} (y_i - X_i beta)`.
#'
#' Its conditional sampling covariance is `(Z_i' R_i^{-1} Z_i)^{-1}`. When
#' `R_list` is `NULL`, `R_i` defaults to `sigma^2 I`, the `sigma^2` factor
#' cancels from the coefficient equation, and the score reduces to the original
#' OLS closed form `(Z_i'Z_i)^{-1} Z_i'r_i`; the covariance becomes
#' `sigma^2 (Z_i'Z_i)^{-1}`.
#'
#' The iid path is performance-specialized because it is the common simulation
#' case. It batches clusters by identical `Z_i`, solves `(Z_i'Z_i)^{-1}` once
#' per unique design, and multiplies that projection into all residual vectors
#' for the design group at once. The explicit-`R_list` path keeps the direct GLS
#' formula cluster by cluster because arbitrary non-diagonal `R_i` generally
#' differs across clusters.
#'
#' For a random-intercept-only model, the returned score is the GLS intercept of
#' the fixed-effect residuals. For a random-intercept/random-slope model, the
#' returned scores are the GLS intercept and within-cluster slope from residuals
#' on `(1, within_var)` inside each cluster.
#'
#' The function is useful for balanced comparisons where the likelihood-only
#' cluster score has a direct closed form and should not depend on the
#' prior-unweighting matrix algebra used by `get_corrected_scores()`. The
#' returned `ols_var*` columns are retained as compatibility aliases for the
#' measurement-error estimators; under non-diagonal `R_i` they contain GLS/MLE
#' score variances.
#'
#' @param fit_obj Fitted `lme4` model whose fixed effects are used to residualize
#' the outcome.
#' @param data Original long-format data frame used for the mixed model.
#' @param cluster_var Character scalar naming the cluster/grouping column.
#' @param outcome_var Character scalar naming the outcome column.
#' @param within_var Optional character scalar naming the within-cluster
#' predictor. If `NULL`, an intercept-only cluster score is computed.
#' @param R_list Optional named list of cluster-level residual covariance
#' matrices. Names should match `cluster_var` values. If unnamed, matrices are
#' matched to clusters in the order they appear in `data`. If `NULL`, iid
#' residual covariance `sigma(fit_obj)^2 I` is used for each cluster.
#'
#' @return
#' A tibble with one row per cluster and an `id` column. Intercept-only calls
#' return `corrected` and `corrected_intercept_full`; calls with `within_var`
#' return `corrected_intercept_full` and `corrected_slope_full`. Singular or
#' numerically invalid cluster systems return `NA_real_` values for that
#' cluster's corrected scores and variances.
get_closed_form_corrected_scores <- function(fit_obj, data, cluster_var, outcome_var, within_var = NULL,
                                             R_list = NULL) {
  cluster_ids <- unique(as.character(data[[cluster_var]]))
  split_dat <- split(data, as.character(data[[cluster_var]]), drop = TRUE)[cluster_ids]
  beta_hat <- stage1_fixef(fit_obj)

  if (is.null(R_list) && inherits(fit_obj, "lme")) {
    R_list <- extract_stage1_components(
      fit_obj = fit_obj,
      data = data,
      cluster_var = cluster_var,
      within_var = within_var
    )$R_list
  }

  if (!is.null(R_list) && !is.list(R_list)) {
    stop("`R_list` must be NULL or a list of cluster-level residual covariance matrices.")
  }
  if (!is.null(R_list) && is.null(names(R_list))) {
    if (length(R_list) != length(cluster_ids)) {
      stop("Unnamed `R_list` must have one matrix per cluster in `data` order.")
    }
    names(R_list) <- cluster_ids
  }

  if (is.null(R_list)) {
    sigma2_hat <- stats::sigma(fit_obj)^2

    # Fast iid branch: R_i = sigma^2 I. The scalar sigma^2 does not affect the
    # coefficient estimate, so the score is OLS on fixed-effect residuals. We
    # still use sigma^2 below for the score covariance.
    id_vec <- as.character(data[[cluster_var]])

    if (is.null(within_var)) {
      # Random intercept only:
      #   Z_i = 1_i
      #   r_i = y_i - beta_0
      resid_all <- data[[outcome_var]] - beta_hat[[1]]
      resid_by_id <- split(resid_all, id_vec, drop = TRUE)[cluster_ids]
      z_by_id <- lapply(resid_by_id, function(resid_i) matrix(1, nrow = length(resid_i), ncol = 1L))
      n_re <- 1L
    } else {
      # Random intercept and slope:
      #   Z_i = X_i = [1, z_i]
      #   r_i = y_i - beta_0 - beta_z z_i
      # This assumes the fixed-effect and random-effect designs match, as in
      # the current random-slope simulations.
      z_all <- data[[within_var]]
      resid_all <- data[[outcome_var]] - beta_hat[[1]] - beta_hat[[within_var]] * z_all
      resid_by_id <- split(resid_all, id_vec, drop = TRUE)[cluster_ids]
      z_vec_by_id <- split(z_all, id_vec, drop = TRUE)[cluster_ids]
      z_by_id <- lapply(z_vec_by_id, function(z_vec) cbind(1, z_vec))
      n_re <- 2L
    }

    # Bucket clusters by exact Z design. Balanced simulations usually have one
    # bucket; unbalanced simulations often have one bucket per realized trial
    # count because z is generated deterministically from the count.
    design_buckets <- new.env(parent = emptyenv())
    for (i in seq_along(cluster_ids)) {
      z_mat <- z_by_id[[i]]
      z_values <- as.numeric(z_mat)
      # First partition by dimensions to avoid comparing vectors of different
      # lengths, then compare the full numeric design for exact reuse.
      bucket_key <- paste(nrow(z_mat), ncol(z_mat), sep = "|")
      bucket <- if (exists(bucket_key, envir = design_buckets, inherits = FALSE)) {
        get(bucket_key, envir = design_buckets, inherits = FALSE)
      } else {
        list()
      }

      matched <- FALSE
      if (length(bucket) > 0L) {
        for (j in seq_along(bucket)) {
          if (identical(z_values, bucket[[j]]$z_values)) {
            bucket[[j]]$indices <- c(bucket[[j]]$indices, i)
            matched <- TRUE
            break
          }
        }
      }
      if (!matched) {
        bucket[[length(bucket) + 1L]] <- list(
          z_values = z_values,
          z_mat = z_mat,
          indices = i
        )
      }
      assign(bucket_key, bucket, envir = design_buckets)
    }

    corrected <- matrix(NA_real_, nrow = n_re, ncol = length(cluster_ids))
    vcov_arr <- array(NA_real_, dim = c(n_re, n_re, length(cluster_ids)))

    # For each unique design:
    #   P_Z = (Z'Z)^-1 Z'
    #   B_group = P_Z R_group
    # where R_group is the matrix of residual vectors for all clusters sharing Z.
    for (bucket_key in ls(design_buckets, all.names = TRUE)) {
      bucket <- get(bucket_key, envir = design_buckets, inherits = FALSE)
      for (entry in bucket) {
        idx <- entry$indices
        z_mat <- entry$z_mat
        ztz_inv <- tryCatch(solve(crossprod(z_mat)), error = function(e) NULL)
        if (is.null(ztz_inv)) {
          next
        }

        resid_mat <- do.call(cbind, resid_by_id[idx])
        corrected[, idx] <- ztz_inv %*% crossprod(z_mat, resid_mat)
        vcov_i <- ztz_inv * sigma2_hat
        # Every cluster in this design group has the same iid score covariance.
        for (cluster_pos in idx) {
          vcov_arr[, , cluster_pos] <- vcov_i
        }
      }
    }

    out <- tibble::tibble(id = cluster_ids)
    if (n_re == 1L) {
      out$corrected <- corrected[1, ]
      out$corrected_intercept_full <- corrected[1, ]
      out$gls_var11 <- vcov_arr[1, 1, ]
      out$ols_var11 <- vcov_arr[1, 1, ]
    } else {
      out$corrected_intercept_full <- corrected[1, ]
      out$corrected_slope_full <- corrected[2, ]
      out$gls_var11 <- vcov_arr[1, 1, ]
      out$gls_var12 <- vcov_arr[1, 2, ]
      out$gls_var22 <- vcov_arr[2, 2, ]
      out$ols_var11 <- vcov_arr[1, 1, ]
      out$ols_var12 <- vcov_arr[1, 2, ]
      out$ols_var22 <- vcov_arr[2, 2, ]
    }
    return(out)
  }

  # General GLS branch for explicit R_i. This is the path used for AR(1),
  # Toeplitz, or any other non-diagonal residual covariance supplied by the
  # simulator. It intentionally mirrors the formula in the roxygen block.
  purrr::map_dfr(cluster_ids, function(cluster_id) {
    df_i <- split_dat[[cluster_id]]

    if (is.null(within_var)) {
      # Random-intercept-only case: the cluster score is the intercept from a
      # regression of fixed-effect residuals on a column of ones.
      z_mat <- matrix(1, nrow = nrow(df_i), ncol = 1L)
      x_mat <- z_mat
      beta_vec <- beta_hat[[1]]
    } else {
      # Random intercept and slope case: use the same fixed-effect design for
      # residualization, then solve the within-cluster random-effect design.
      z_vec <- df_i[[within_var]]
      z_mat <- cbind(1, z_vec)
      x_mat <- z_mat
      beta_vec <- c(beta_hat[[1]], beta_hat[[within_var]])
    }

    # Residualize using fixed effects only, then estimate the cluster-specific
    # likelihood-only coefficients from the random-effect design. With the
    # default iid R_i this is exactly the original OLS closed form; with
    # non-diagonal R_i it is the MLE/Bartlett GLS score.
    resid_i <- df_i[[outcome_var]] - as.numeric(x_mat %*% beta_vec)

    R_i <- as.matrix(R_list[[cluster_id]])

    gls <- tryCatch({
      if (!is.matrix(R_i) || nrow(R_i) != nrow(df_i) || ncol(R_i) != nrow(df_i)) {
        stop("R_i has incompatible dimensions.")
      }
      # Avoid forming R_i^-1 explicitly. solve(R_i, Z_i) and solve(R_i, r_i)
      # supply R_i^-1 Z_i and R_i^-1 r_i for the GLS normal equations.
      R_inv_Z <- solve(R_i, z_mat)
      R_inv_resid <- solve(R_i, resid_i)
      z_Rinv_z <- crossprod(z_mat, R_inv_Z)
      gls_vcov <- solve(z_Rinv_z)
      list(
        corrected = as.numeric(gls_vcov %*% crossprod(z_mat, R_inv_resid)),
        vcov = gls_vcov
      )
    }, error = function(e) {
      list(
        corrected = rep(NA_real_, ncol(z_mat)),
        vcov = matrix(NA_real_, nrow = ncol(z_mat), ncol = ncol(z_mat))
      )
    })

    out <- tibble::tibble(id = cluster_id)
    if (ncol(z_mat) == 1L) {
      out$corrected <- gls$corrected[[1]]
      out$corrected_intercept_full <- gls$corrected[[1]]
      out$gls_var11 <- gls$vcov[1, 1]
      out$ols_var11 <- gls$vcov[1, 1]
    } else {
      out$corrected_intercept_full <- gls$corrected[[1]]
      out$corrected_slope_full <- gls$corrected[[2]]
      out$gls_var11 <- gls$vcov[1, 1]
      out$gls_var12 <- gls$vcov[1, 2]
      out$gls_var22 <- gls$vcov[2, 2]
      out$ols_var11 <- gls$vcov[1, 1]
      out$ols_var12 <- gls$vcov[1, 2]
      out$ols_var22 <- gls$vcov[2, 2]
    }
    out
  })
}



#' Compute various cluster scores from the same Stage-1 ingredients.
run_scoring_methods <- function(fit_obj, data, cluster_var, outcome_var, within_var = NULL,
                                     R_list = NULL, group = NULL) {
  cluster_ids <- unique(as.character(data[[cluster_var]]))
  split_dat <- split(data, as.character(data[[cluster_var]]), drop = TRUE)[cluster_ids]
  n_re <- if (is.null(within_var)) 1L else 2L

  stage1 <- extract_stage1_components(
    fit_obj = fit_obj,
    data = data,
    cluster_var = cluster_var,
    within_var = within_var,
    R_list = R_list,
    group = group
  )
  beta_hat <- stage1$beta_hat
  g_hat <- stage1$G_hat
  R_list <- normalize_R_list(stage1$R_list, cluster_ids)

  re_names_raw <- stage1$re_names_raw
  if (is.null(re_names_raw) || length(re_names_raw) != n_re) {
    re_names_raw <- if (is.null(within_var)) "(Intercept)" else c("(Intercept)", within_var)
  }
  re_names <- sanitize_re_name(re_names_raw)

  if (nrow(g_hat) != n_re || ncol(g_hat) != n_re) {
    stop("The fitted random-effect covariance dimension does not match `within_var`.")
  }
  purrr::map_dfr(cluster_ids, function(cluster_id) {
    df_i <- split_dat[[cluster_id]]
    if (is.null(within_var)) {
      z_mat <- matrix(1, nrow = nrow(df_i), ncol = 1L)
      x_mat <- z_mat
      beta_vec <- beta_hat[[1]]
    } else {
      z_vec <- df_i[[within_var]]
      z_mat <- cbind(1, z_vec)
      x_mat <- z_mat
      beta_vec <- c(beta_hat[[1]], beta_hat[[within_var]])
    }

    R_i <- as.matrix(R_list[[cluster_id]])

    try_scoring_methods(
      z_mat = z_mat,
      x_mat = x_mat,
      y_vec = df_i[[outcome_var]],
      beta_vec = beta_vec,
      g_hat = g_hat,
      R_i = R_i,
      n_re = n_re
    )
  })
}

#' Try multiple scoring methods on the same cluster data, returning a list of
#' results with error handling for each method.
try_scoring_methods <- function(
  z_mat,
  x_mat,
  y_vec,
  beta_vec,
  g_hat,
  R_i,
  n_re
) {
  methods <- c("closed_form", "unweight_full", "unweight_diag", "g_r_conversion", "g_sigma_conversion")

  # if (!is.matrix(R_i) || nrow(R_i) != nrow(df_i) || ncol(R_i) != nrow(df_i)) {
  #       stop("R_i has incompatible dimensions.")
  #     }

  dplyr::bind_rows(
    try_scoring_closed_form(
      z_mat = z_mat,
      x_mat = x_mat,
      y_vec = y_vec,
      beta_vec = beta_vec,
      R_i = R_i,
      n_re = n_re
    ),
    try_scoring_unweight_full(
      z_mat = z_mat,
      x_mat = x_mat,
      y_vec = y_vec,
      beta_vec = beta_vec,
      g_hat = g_hat,
      R_i = R_i,
      n_re = n_re
    ),
    try_scoring_unweight_diag(
      z_mat = z_mat,
      x_mat = x_mat,
      y_vec = y_vec,
      beta_vec = beta_vec,
      g_hat = g_hat,
      R_i = R_i,
      n_re = n_re
    ),
    try_scoring_g_r_conversion(
      z_mat = z_mat,
      x_mat = x_mat,
      y_vec = y_vec,
      beta_vec = beta_vec,
      g_hat = g_hat,
      R_i = R_i,
      n_re = n_re
    ),
    try_scoring_g_sigma_conversion(
      z_mat = z_mat,
      x_mat = x_mat,
      y_vec = y_vec,
      beta_vec = beta_vec,
      g_hat = g_hat,
      R_i = R_i,
      n_re = n_re
    )
  )
}

try_scoring_closed_form <- function(
  z_mat,
  x_mat,
  y_vec,
  beta_vec,
  R_i,
  n_re
) {

  resid_i <- y_vec - as.numeric(x_mat %*% beta_vec)
  
  result <- tryCatch({

    failed_step <- "r inverse z"
    R_inv_Z <- solve(R_i, z_mat)

    failed_step <- "r inverse resid"
    R_inv_resid <- solve(R_i, resid_i)

    failed_step <- "z crossprod r inverse z"
    z_Rinv_z <- crossprod(z_mat, R_inv_Z)

    failed_step <- "z_r_z inverse"
    gls_vcov <- solve(z_Rinv_z)

    failed_step <- "closed form score calculation"
    gls_score <- as.numeric(gls_vcov %*% crossprod(z_mat, R_inv_resid))

    dplyr::tibble(
      method = "closed_form",
      success = TRUE,
      err_step = "none",
      err_msg = "none",
      score1 = gls_score[1],
      score2 = if (n_re > 1L) gls_score[2] else NA_real_
    )
  }, error = function(e) {
    dplyr::tibble(
      method = "closed_form",
      success = FALSE,
      err_step = failed_step,
      err_msg = e$message,
      score1 = NA_real_,
      score2 = NA_real_
    )
  })

  return(result)
}

try_scoring_unweight_full <- function(
  z_mat,
  x_mat,
  y_vec,
  beta_vec,
  g_hat,
  R_i,
  n_re
) {
  
  resid_i <- y_vec - as.numeric(x_mat %*% beta_vec)
  sigma_y_i <- z_mat %*% g_hat %*% t(z_mat) + R_i

  result <- tryCatch({
    failed_step <- "sigma y inverse"
    sigma_y_inv <- solve(sigma_y_i)

    failed_step <- "shrinkage matrix calculation"
    a_i <- g_hat %*% t(z_mat) %*% sigma_y_inv

    failed_step <- "eb score calculation"
    eb_score <- as.numeric(a_i %*% resid_i)

    failed_step <- "posterior covariance calculation"
    post_vcov <- g_hat - a_i %*% z_mat %*% g_hat
    post_vcov <- (post_vcov + t(post_vcov)) / 2

    failed_step <- "uw score calculation"
    uw <- solve(
      solve(post_vcov) - solve(g_hat),
      solve(post_vcov, eb_score)
    )
    
    dplyr::tibble(
      method = "unweight_full",
      success = TRUE,
      err_step = "none",
      err_msg = "none",
      score1 = uw[1],
      score2 = if (n_re > 1L) uw[2] else NA_real_
    )
  }, error = function(e) {
    dplyr::tibble(
      method = "unweight_full",
      success = FALSE,
      err_step = failed_step,
      err_msg = e$message,
      score1 = NA_real_,
      score2 = NA_real_
    )
  })

  return(result)
}

try_scoring_unweight_diag <- function(
  z_mat,
  x_mat,
  y_vec,
  beta_vec,
  g_hat,
  R_i,
  n_re
) {
  
  resid_i <- y_vec - as.numeric(x_mat %*% beta_vec)
  sigma_y_i <- z_mat %*% g_hat %*% t(z_mat) + R_i

  result <- tryCatch({
    failed_step <- "sigma y inverse"
    sigma_y_inv <- solve(sigma_y_i)

    failed_step <- "shrinkage matrix calculation"
    a_i <- g_hat %*% t(z_mat) %*% sigma_y_inv

    failed_step <- "eb score calculation"
    eb_score <- as.numeric(a_i %*% resid_i)

    failed_step <- "posterior covariance calculation"
    post_vcov <- g_hat - a_i %*% z_mat %*% g_hat
    post_vcov <- (post_vcov + t(post_vcov)) / 2

    diag_corrected <- rep(NA_real_, n_re)
    for (j in seq_len(n_re)) {
      failed_step <- paste("unweight diag calculation for component", j)
      denom <- (1 / post_vcov[j, j]) - (1 / g_hat[j, j])
      if (is.finite(denom) && denom > sqrt(.Machine$double.eps)) {
        failed_step <- paste("unweight diag score calculation for component", j)
        diag_corrected[[j]] <- (1 / denom) * eb_score[[j]] / post_vcov[j, j]
      }
    }
    
    dplyr::tibble(
      method = "unweight_diag",
      success = TRUE,
      err_step = "none",
      err_msg = "none",
      score1 = diag_corrected[1],
      score2 = if (n_re > 1L) diag_corrected[2] else NA_real_
    )
  }, error = function(e) {
    dplyr::tibble(
      method = "unweight_diag",
      success = FALSE,
      err_step = failed_step,
      err_msg = e$message,
      score1 = NA_real_,
      score2 = NA_real_
    )
  })

  return(result)
}

try_scoring_g_r_conversion <- function(
  z_mat,
  x_mat,
  y_vec,
  beta_vec,
  g_hat,
  R_i,
  n_re
) {
  
  resid_i <- y_vec - as.numeric(x_mat %*% beta_vec)
  sigma_y_i <- z_mat %*% g_hat %*% t(z_mat) + R_i

  result <- tryCatch({
    failed_step <- "g_z_r_inv_z calculation"
    gzrz <- g_hat %*% t(z_mat) %*% solve(R_i, z_mat)

    # failed_step <- "g_z_r_inv_z inverse"
    # gzrz_inv <- solve(gzrz)
    # g_r_conversion <- gzrz_inv + diag(n_re)

    # # regular EB calculations
    # failed_step <- "sigma y inverse"
    # sigma_y_inv <- solve(sigma_y_i)

    # failed_step <- "shrinkage matrix calculation"
    # a_i <- g_hat %*% t(z_mat) %*% sigma_y_inv

    failed_step <- "eb score calculation"
    eb_score <- g_hat %*% t(z_mat) %*% solve(sigma_y_i, resid_i)

    # final conversion step
    failed_step <- "g_r conversion score calculation"
    g_r_score <- as.numeric(solve(gzrz, eb_score) + eb_score)

    dplyr::tibble(
      method = "g_r_conversion",
      success = TRUE,
      err_step = "none",
      err_msg = "none",
      score1 = g_r_score[1],
      score2 = if (n_re > 1L) g_r_score[2] else NA_real_
    )

    }, error = function(e) {
      dplyr::tibble(
        method = "g_r_conversion",
        success = FALSE,
        err_step = failed_step,
        err_msg = e$message,
        score1 = NA_real_,
        score2 = NA_real_
      )
  })

  return(result)

}

try_scoring_g_sigma_conversion <- function(
  z_mat,
  x_mat,
  y_vec,
  beta_vec,
  g_hat,
  R_i,
  n_re
) {
  
  resid_i <- y_vec - as.numeric(x_mat %*% beta_vec)
  sigma_y_i <- z_mat %*% g_hat %*% t(z_mat) + R_i

  result <- tryCatch({

    # failed_step <- "sigma y inverse"
    # sigma_y_inv <- solve(sigma_y_i)

    failed_step <- "g_z_sigma_inv_z calculation"
    gzsz <- g_hat %*% t(z_mat) %*% solve(sigma_y_i, z_mat)

    # failed_step <- "g_z_sigma_inv_z inverse"
    # g_sigma_conversion <- solve(gzsz)

    # regular EB calculations
    failed_step <- "eb score calculation"
    eb_score <- g_hat %*% t(z_mat) %*% solve(sigma_y_i, resid_i)

    # final conversion step
    failed_step <- "g_sigma conversion score calculation"
    g_sigma_score <- as.numeric(solve(gzsz, eb_score))

    dplyr::tibble(
      method = "g_sigma_conversion",
      success = TRUE,
      err_step = "none",
      err_msg = "none",
      score1 = g_sigma_score[1],
      score2 = if (n_re > 1L) g_sigma_score[2] else NA_real_
    )

    }, error = function(e) {
      dplyr::tibble(
        method = "g_sigma_conversion",
        success = FALSE,
        err_step = failed_step,
        err_msg = e$message,
        score1 = NA_real_,
        score2 = NA_real_
      )
  })

  return(result)

}