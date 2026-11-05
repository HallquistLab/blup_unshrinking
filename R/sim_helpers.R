#' ---
#' title: "Simulation Helpers for MLM BLUP Correction"
#' description: "Core data generation functions for univariate and bivariate multilevel models."
#' ---

#' Draw random effects for intercept and slope from a bivariate normal distribution.
#'
#' This function generates subject-specific random deviations (u_0i, u_1i) from a 
#' zero-mean bivariate normal distribution with specified standard deviations and correlation.
#' It uses the Cholesky decomposition of the covariance matrix.
#'
#' @param n_id The number of subjects (clusters).
#' @param tau0 The standard deviation of the random intercept (tau_0).
#' @param tau1 The standard deviation of the random slope (tau_1).
#' @param rho The correlation between the random intercept and random slope.
#' @return A matrix of size `n_id` x 2 containing random effects (u0, u1).
draw_random_effects <- function(n_id, tau0, tau1, rho) {
  # Construct the 2x2 covariance matrix G
  sd_vec <- c(tau0, tau1)
  cor_mat <- matrix(c(1, rho, rho, 1), nrow = 2)
  cov_mat <- diag(sd_vec) %*% cor_mat %*% diag(sd_vec)
  
  # Generate standard normal draws and correlate them using Cholesky factor
  z_draws <- matrix(stats::rnorm(n_id * 2L), ncol = 2)
  re <- z_draws %*% chol(cov_mat)
  colnames(re) <- c("u0", "u1")
  re
}

#' Standardize a residual covariance specification.
#'
#' @details
#' `simulate_dataset()` uses this small specification object to generate the
#' within-cluster residual covariance matrix `R_i`. The default `"iid"`
#' structure preserves the original simulation behavior.
#'
#' @param r_spec Residual covariance specification. Supported forms are `NULL`,
#'   a character scalar (`"iid"`, `"ar1"`, or `"toeplitz"`), or a list with a
#'   `structure` or `type` field plus structure-specific parameters.
#'
#' @return A normalized list with a lower-case `structure` field.
normalize_r_spec <- function(r_spec = NULL) {
  if (is.null(r_spec)) {
    r_spec <- list(structure = "iid")
  }
  if (is.character(r_spec) && length(r_spec) == 1L) {
    r_spec <- list(structure = r_spec)
  }
  if (!is.list(r_spec)) {
    stop("`r_spec` must be NULL, a character scalar, or a list.")
  }

  structure <- r_spec$structure
  if (is.null(structure)) {
    structure <- r_spec$type
  }
  if (is.null(structure) || length(structure) != 1L || !is.character(structure)) {
    stop("`r_spec` must include a character `structure` field.")
  }
  r_spec$structure <- tolower(structure)
  r_spec$type <- NULL

  if (!(r_spec$structure %in% c("iid", "ar1", "toeplitz"))) {
    stop("Unsupported `r_spec$structure`: ", r_spec$structure)
  }
  if (identical(r_spec$structure, "ar1")) {
    rho <- r_spec$rho
    if (is.null(rho) || length(rho) != 1L || !is.finite(rho) || abs(rho) >= 1) {
      stop("AR(1) `r_spec` requires finite scalar `rho` with abs(rho) < 1.")
    }
    r_spec$rho <- as.numeric(rho)
  }
  if (identical(r_spec$structure, "toeplitz")) {
    correlations <- r_spec$correlations
    if (is.null(correlations)) {
      correlations <- r_spec$rho
    }
    if (is.null(correlations) || length(correlations) < 1L ||
      any(!is.finite(correlations)) || any(abs(correlations) >= 1)) {
      stop("Toeplitz `r_spec` requires finite lag `correlations` with absolute values < 1.")
    }
    r_spec$correlations <- as.numeric(correlations)
    r_spec$rho <- NULL
  }

  r_spec
}

#' Construct a cluster-level residual covariance matrix.
#'
#' @details
#' This factory returns `R_i = sigma^2 C_i`, where `C_i` is an iid, AR(1), or
#' Toeplitz residual correlation matrix. Toeplitz lag correlations shorter than
#' `n - 1` are padded with zeros, making finite-band structures easy to specify.
#'
#' @param n Cluster size.
#' @param sigma Residual standard deviation.
#' @param r_spec Normalized residual covariance specification from
#'   `normalize_r_spec()`.
#'
#' @return A positive-definite `n x n` covariance matrix.
make_R_matrix <- function(n, sigma, r_spec = NULL) {
  n <- as.integer(n[[1]])
  sigma <- as.numeric(sigma[[1]])
  r_spec <- normalize_r_spec(r_spec)

  if (!is.finite(n) || n < 1L) {
    stop("`n` must be a positive integer.")
  }
  if (!is.finite(sigma) || sigma <= 0) {
    stop("`sigma` must be positive.")
  }

  cor_mat <- switch(
    r_spec$structure,
    iid = diag(n),
    ar1 = {
      idx <- seq_len(n)
      outer(idx, idx, function(a, b) r_spec$rho^abs(a - b))
    },
    toeplitz = {
      lag_cor <- numeric(n)
      lag_cor[[1]] <- 1
      n_lag <- min(n - 1L, length(r_spec$correlations))
      if (n_lag > 0L) {
        lag_cor[seq_len(n_lag) + 1L] <- r_spec$correlations[seq_len(n_lag)]
      }
      stats::toeplitz(lag_cor)
    }
  )

  r_mat <- (sigma^2) * cor_mat
  chol_ok <- tryCatch({
    chol(r_mat)
    TRUE
  }, error = function(e) FALSE)
  if (!chol_ok) {
    stop("The requested residual covariance matrix is not positive definite.")
  }

  r_mat
}

#' Draw residuals from a cluster-level residual covariance matrix.
#'
#' @param r_mat Positive-definite residual covariance matrix.
#'
#' @return Numeric residual vector with covariance `r_mat`.
draw_residuals_from_R <- function(r_mat) {
  chol_r <- chol(r_mat)
  drop(t(chol_r) %*% stats::rnorm(nrow(r_mat)))
}

#' Unified dataset simulation for random intercept and random slope MLMs.
#'
#' This function generates hierarchical data according to a specific Data Generating Process (DGP).
#' It supports two primary modes:
#' 
#' 1. Random Intercept Only (`has_random_slope = FALSE`):
#'    Level 1: y_{ij} = \beta_0 + \eta_i + e_{ij},  where e_{ij} ~ N(0, \sigma^2)
#'    Level 2: \eta_i = \beta_1 x_i + u_{0i},       where u_{0i} ~ N(0, \tau^2)
#'    (where \eta_i is the true latent intercept for subject i)
#'
#' 2. Random Slope (`has_random_slope = TRUE`):
#'    Level 1: y_{ij} = \beta_0 + u_{0i} + (\beta_z + u_{1i}*) z_{ij} + e_{ij},  where e_{ij} ~ N(0, \sigma^2)
#'    Level 2: u_{1i}* = \gamma x_i + u_{1i},  where (u_{0i}, u_{1i}) ~ MVN(0, G)
#'    (where u_{1i}* is the true latent slope for subject i, and x_i predicts this slope)
#'
#' @param n_id Number of subjects.
#' @param mean_n_trial Target mean trials per subject (j index).
#' @param params List of population parameters: 
#'        - beta_0 (fixed intercept)
#'        - beta_z (fixed slope for level-1 predictor z)
#'        - gamma_x_on_slope (effect of level-2 predictor x on the random slope u_1i)
#'        - beta_1 (effect of level-2 predictor x on the random intercept, used if has_random_slope=FALSE)
#'        - tau0 / tau (standard deviation of the random intercept)
#'        - tau1 (standard deviation of the random slope)
#'        - rho (correlation between intercept and slope)
#'        - sigma (residual standard error)
#' @param tau1 Standard deviation of the random slope (if relevant, can override params$tau1).
#' @param sigma Residual standard error (can override params$sigma).
#' @param has_random_slope Logical; if TRUE, generates a random-slope design (bivariate).
#' @param balanced Logical; if TRUE, all subjects have exactly `mean_n_trial` trials. If FALSE, trials are drawn uniformly from [0.6*mean, 1.4*mean].
#' @param min_n_trial Minimum trials per subject in the jittered unbalanced design.
#' @param highly_unbalanced_min_n_trial Minimum trials per subject in the highly unbalanced design.
#' @param highly_unbalanced_power Power used to map the rank-based imbalance profile to trial counts.
#' @param r_spec Residual covariance specification for cluster-level `R_i`
#'   matrices. Defaults to iid residuals. Supported structures are `"iid"`,
#'   `list(structure = "ar1", rho = ...)`, and
#'   `list(structure = "toeplitz", correlations = c(...))`.
#' @return A list containing:
#'         - `dat`: The full long-format data frame.
#'         - `id_df`: Subject-level data containing true latent variables and predictors.
#'         - `R_list`: Named list of generated cluster-level residual covariance matrices.
#'         - `r_spec`: Normalized residual covariance specification used for generation.
#'         - `mean_realized_trials`: The empirical mean number of trials generated per subject.
#'         - `min_realized_trials`: The minimum realized trial count in the generated sample.
#'         - `prop_ids_leq_2_trials`: Proportion of subjects with 2 or fewer trials.
#'         - `prop_ids_leq_3_trials`: Proportion of subjects with 3 or fewer trials.
simulate_dataset <- function(n_id, mean_n_trial, params, tau1 = NULL, sigma = NULL,
                             has_random_slope = TRUE, balanced = FALSE,
                             min_n_trial = 6L, highly_unbalanced_min_n_trial = 4L,
                             highly_unbalanced_power = 3, r_spec = NULL) {
  n_id <- as.integer(n_id[[1]])
  mean_n_trial <- as.integer(mean_n_trial[[1]])
  sigma <- if (!is.null(sigma)) as.numeric(sigma[[1]]) else params$sigma
  r_spec <- normalize_r_spec(r_spec = r_spec)
  min_n_trial <- as.integer(min_n_trial[[1]])
  highly_unbalanced_min_n_trial <- as.integer(highly_unbalanced_min_n_trial[[1]])
  highly_unbalanced_power <- as.numeric(highly_unbalanced_power[[1]])

  if (!is.finite(min_n_trial) || min_n_trial < 1L) {
    stop("`min_n_trial` must be an integer >= 1.")
  }
  if (!is.finite(highly_unbalanced_min_n_trial) || highly_unbalanced_min_n_trial < 1L) {
    stop("`highly_unbalanced_min_n_trial` must be an integer >= 1.")
  }
  if (!is.finite(highly_unbalanced_power) || highly_unbalanced_power <= 0) {
    stop("`highly_unbalanced_power` must be positive.")
  }

  # Subjects' level-2 covariate x (standardized)
  x_subj <- scale(stats::rnorm(n_id), center = TRUE, scale = TRUE)[, 1]

  if (has_random_slope) {
    tau1 <- if (!is.null(tau1)) as.numeric(tau1[[1]]) else params$tau1
    
    # Draw correlated random effects (u_0i, u_1i)
    re <- draw_random_effects(n_id = n_id, tau0 = params$tau0, tau1 = tau1, rho = params$rho)
    
    # Intercept deviations are just the drawn u_0i
    intercept_dev_subj <- re[, "u0"]
    
    # The true latent slope (u_1i*) is predicted by x_i with coefficient gamma, plus unexplained variance u_1i
    slope_dev_subj <- params$gamma_x_on_slope * x_subj + re[, "u1"]
  } else {
    # Random intercept only: independent draws for u_0i
    intercept_dev_subj <- stats::rnorm(n_id, mean = 0, sd = params$tau)
    slope_dev_subj <- rep(0, n_id) # No random slope present
  }

  # Determine number of trials (j) for each subject (i)
  if (is.character(balanced) && balanced == "highly_unbalanced") {
    # Generate extreme imbalance: number of trials is correlated with x
    # This induces heteroscedasticity in the estimation error of the random effects
    if (mean_n_trial < highly_unbalanced_min_n_trial) {
      stop("`mean_n_trial` must be >= `highly_unbalanced_min_n_trial`.")
    }
    x_percentile <- rank(x_subj) / n_id
    rev_percentile <- 1 - x_percentile # High x gets low trials (higher variance)
    amplitude <- (mean_n_trial - highly_unbalanced_min_n_trial) * (highly_unbalanced_power + 1)
    base_trials <- highly_unbalanced_min_n_trial + amplitude * (rev_percentile ^ highly_unbalanced_power)
    trial_counts <- pmax(highly_unbalanced_min_n_trial, as.integer(round(base_trials)))
  } else if (isTRUE(balanced) || (is.character(balanced) && balanced == "balanced")) {
    trial_counts <- rep(mean_n_trial, n_id)
  } else {
    # Unbalanced: jitter trials +/- 40% around the mean with a configurable lower bound.
    trial_counts <- pmax(min_n_trial, as.integer(round(stats::runif(n_id, min = 0.6, max = 1.4) * mean_n_trial)))
  }

  cluster_ids <- as.character(seq_len(n_id))
  R_list <- stats::setNames(
    lapply(trial_counts, make_R_matrix, sigma = sigma, r_spec = r_spec),
    cluster_ids
  )

  # Generate level-1 observations
  dat <- purrr::map_dfr(seq_len(n_id), function(i) {
    # Within-subject covariate z_{ij} (centered/scaled to range roughly [-1, 1])
    # If has_random_slope is FALSE, z is not used in the generating model except perhaps as a placeholder
    raw_z <- seq(-1, 1, length.out = trial_counts[[i]])
    if (length(raw_z) == 1L) {
      z_i <- 0
    } else {
      z_scaled <- scale(raw_z, center = TRUE, scale = TRUE)[, 1]
      if (anyNA(z_scaled)) {
        z_i <- raw_z - mean(raw_z)
      } else {
        z_i <- z_scaled
      }
    }

    e_i <- draw_residuals_from_R(R_list[[i]])
    if (has_random_slope) {
      # y_{ij} = \beta_0 + u_{0i} + (\beta_z + u_{1i}*) z_{ij} + e_{ij}
      y_i <- params$beta_0 + intercept_dev_subj[[i]] + (params$beta_z + slope_dev_subj[[i]]) * z_i + e_i
    } else {
      # In the intercept-only case, params$beta_1 is the effect of x on the latent intercept
      # \eta_i = \beta_1 x_i + u_{0i}
      eta_i <- params$beta_1 * x_subj[[i]] + intercept_dev_subj[[i]]
      
      # y_{ij} = \beta_0 + \eta_i + e_{ij}
      y_i <- params$beta_0 + eta_i + e_i
    }

    tibble(
      id = as.character(i),
      trial_index = seq_len(trial_counts[[i]]),
      x = x_subj[[i]],
      z = z_i,
      true_intercept_dev = intercept_dev_subj[[i]],
      true_slope_dev = slope_dev_subj[[i]],
      eta = if (!has_random_slope) params$beta_1 * x_subj[[i]] + intercept_dev_subj[[i]] else NA_real_,
      y = y_i
    )
  })

  # Extract true latent values at the subject level (Level 2)
  id_df <- dat %>%
    distinct(id, x, true_intercept_dev, true_slope_dev, eta)

  list(
    dat = mutate(dat, id = factor(id)),
    id_df = id_df,
    R_list = R_list,
    r_spec = r_spec,
    mean_realized_trials = mean(trial_counts),
    min_realized_trials = min(trial_counts),
    prop_ids_leq_2_trials = mean(trial_counts <= 2L),
    prop_ids_leq_3_trials = mean(trial_counts <= 3L)
  )
}
