#' ---
#' title: "Derivative Backends"
#' description: "Provides interchangeable backends for calculating gradients and Hessians."
#' ---

#' Create a derivative backend function list.
#'
#' This factory function returns a list containing standardized `gradient` and 
#' `hessian` functions, allowing the simulation scripts to easily swap between
#' different methods for calculating derivatives (e.g., finite differences, 
#' numDeriv, or analytical solutions).
#'
#' @param method Character string specifying the backend ("handcoded", "numDeriv", etc.).
#' @param grad_eps Step size for numerical gradients.
#' @param hess_eps Step size for numerical Hessians.
#' @param numderiv_method Method to pass to `numDeriv` if selected.
#' @return A list with `name`, `gradient`, and `hessian` components.
make_derivative_backend <- function(method = c("handcoded", "numDeriv", "merDeriv", "tmb", "analytical"),
                                    grad_eps = 1e-6,
                                    hess_eps = 1e-4,
                                    numderiv_method = "Richardson") {

  method <- match.arg(method)

  finite_diff_gradient <- function(fn, x, eps = grad_eps) {
    x <- as.numeric(x)
    step <- eps * pmax(1, abs(x))
    out <- numeric(length(x))

    for (j in seq_along(x)) {
      x_plus <- x
      x_minus <- x
      x_plus[j] <- x_plus[j] + step[j]
      x_minus[j] <- x_minus[j] - step[j]
      out[j] <- (fn(x_plus) - fn(x_minus)) / (2 * step[j])
    }

    out
  }

  finite_diff_hessian <- function(fn, x, eps = hess_eps) {
    x <- as.numeric(x)
    step <- eps * pmax(1, abs(x))
    p <- length(x)
    out <- matrix(0, p, p)
    f0 <- fn(x)

    for (j in seq_len(p)) {
      x_plus <- x
      x_minus <- x
      x_plus[j] <- x_plus[j] + step[j]
      x_minus[j] <- x_minus[j] - step[j]
      out[j, j] <- (fn(x_plus) - 2 * f0 + fn(x_minus)) / (step[j]^2)
    }

    if (p > 1L) {
      for (j in seq_len(p - 1L)) {
        for (k in (j + 1L):p) {
          x_pp <- x
          x_pm <- x
          x_mp <- x
          x_mm <- x
          x_pp[c(j, k)] <- x_pp[c(j, k)] + c(step[j], step[k])
          x_pm[j] <- x_pm[j] + step[j]
          x_pm[k] <- x_pm[k] - step[k]
          x_mp[j] <- x_mp[j] - step[j]
          x_mp[k] <- x_mp[k] + step[k]
          x_mm[c(j, k)] <- x_mm[c(j, k)] - c(step[j], step[k])
          val <- (fn(x_pp) - fn(x_pm) - fn(x_mp) + fn(x_mm)) / (4 * step[j] * step[k])
          out[j, k] <- val
          out[k, j] <- val
        }
      }
    }

    out
  }

  if (identical(method, "handcoded")) {
    return(list(
      name = "handcoded",
      gradient = function(fn, x) finite_diff_gradient(fn = fn, x = x),
      hessian = function(fn, x) finite_diff_hessian(fn = fn, x = x)
    ))
  }

  if (identical(method, "analytical")) {
    return(list(
      name = "analytical",
      gradient = function(...) stop("analytical backend uses get_exact_stage1_matrices"),
      hessian = function(...) stop("analytical backend uses get_exact_stage1_matrices")
    ))
  }

  if (identical(method, "merDeriv")) {
    return(list(
      name = "merDeriv",
      gradient = function(fn, x) finite_diff_gradient(fn = fn, x = x),
      hessian = function(fn, x) finite_diff_hessian(fn = fn, x = x)
    ))
  }

  if (identical(method, "tmb")) {
    return(list(
      name = "tmb",
      gradient = function(fn, x) finite_diff_gradient(fn = fn, x = x),
      hessian = function(fn, x) finite_diff_hessian(fn = fn, x = x)
    ))
  }

  if (!requireNamespace("numDeriv", quietly = TRUE)) {
    stop("The `numDeriv` package is required for the `numDeriv` backend.")
  }

  list(
    name = "numDeriv",
    gradient = function(fn, x) {
      as.numeric(numDeriv::grad(
        func = function(z) fn(as.numeric(z)),
        x = as.numeric(x),
        method = numderiv_method
      ))
    },
    hessian = function(fn, x) {
      out <- numDeriv::hessian(
        func = function(z) fn(as.numeric(z)),
        x = as.numeric(x),
        method = numderiv_method
      )
      matrix(as.numeric(out), nrow = length(x), ncol = length(x))
    }
  )
}

#' Compute exact analytical score matrix and expected Hessian for stage-1 MLMs.
#'
#' @details
#' This backend bypasses numerical differentiation by directly evaluating the
#' analytical gradients (scores) and the expected information matrix (Hessian)
#' for a Gaussian linear mixed model with a random intercept and random slope.
#'
#' @param cluster_objects A list of cluster-level data structures, each containing
#' `y` (outcome), `X` (fixed-effect design matrix), `Z` (random-effect design matrix),
#' and `n` (cluster size).
#' @param psi_var A numeric vector of length 6 containing the current parameter
#' estimates: `beta0`, `beta1`, `var(u0)`, `cov(u0, u1)`, `var(u1)`, and `sigma^2`.
#'
#' @return A list with two components:
#' \describe{
#'   \item{s1_mat}{An `N \times 6` matrix of cluster-specific score contributions.}
#'   \item{a11_hat}{A `6 \times 6` matrix representing the average expected Hessian across all clusters.}
#' }
get_exact_stage1_matrices <- function(cluster_objects, psi_var) {
  beta <- psi_var[1:2]
  G <- matrix(c(psi_var[3], psi_var[4], psi_var[4], psi_var[5]), 2, 2)
  sig2 <- psi_var[6]

  # Structural derivatives of the 2x2 random-effect covariance matrix G with
  # respect to its three unique elements: var(u0), cov(u0, u1), var(u1).
  dG_list <- list(
    matrix(c(1, 0, 0, 0), 2, 2),
    matrix(c(0, 1, 1, 0), 2, 2),
    matrix(c(0, 0, 0, 1), 2, 2)
  )

  N <- length(cluster_objects)
  s1_mat <- matrix(0, nrow = N, ncol = 6)
  A11 <- matrix(0, nrow = 6, ncol = 6)

  for (i in seq_len(N)) {
    obj <- cluster_objects[[i]]
    Z <- obj$Z
    X <- obj$X
    r <- obj$y - drop(X %*% beta)

    # V is the marginal covariance matrix of y_i: Z G Z' + sigma^2 I
    V <- Z %*% G %*% t(Z) + sig2 * diag(obj$n)
    chol_V <- chol(V)
    Vinv <- chol2inv(chol_V)

    # Precompute common matrix-vector and matrix-matrix products to avoid
    # redundant computation in the score and Hessian equations below.
    Vinv_r <- drop(Vinv %*% r)
    Zt_Vinv_Z <- t(Z) %*% Vinv %*% Z
    Zt_Vinv_r <- t(Z) %*% Vinv_r

    # beta score: Gradient of the log-likelihood with respect to fixed effects.
    s_beta <- drop(t(X) %*% Vinv_r)

    # theta scores: Gradients with respect to the variance components.
    s_theta <- numeric(4)
    for (k in 1:3) {
      dG <- dG_list[[k]]
      tr_term <- sum(Zt_Vinv_Z * dG)
      quad_term <- sum(Zt_Vinv_r * drop(dG %*% Zt_Vinv_r))
      s_theta[k] <- -0.5 * (tr_term - quad_term)
    }
    # sig2 score: Gradient with respect to the residual variance.
    s_theta[4] <- -0.5 * (sum(diag(Vinv)) - sum(Vinv_r * Vinv_r))

    s1_mat[i, ] <- c(s_beta, s_theta)

    # --- Hessian blocks ---
    # The expected Hessian takes a block-diagonal-like structure where the
    # cross-derivatives between fixed effects and variance components are zero.
    A_bb <- t(X) %*% Vinv %*% X

    A_b_theta <- matrix(0, 2, 4)
    Xt_Vinv_Z <- t(X) %*% Vinv %*% Z
    for (k in 1:3) {
      dG <- dG_list[[k]]
      A_b_theta[, k] <- Xt_Vinv_Z %*% dG %*% Zt_Vinv_r
    }
    A_b_theta[, 4] <- t(X) %*% Vinv %*% Vinv_r

    A_tt <- matrix(0, 4, 4)
    # G components with G components
    for (k in 1:3) {
      for (l in k:3) {
        tr_term <- sum(diag(Zt_Vinv_Z %*% dG_list[[k]] %*% Zt_Vinv_Z %*% dG_list[[l]]))
        quad_term <- sum(Zt_Vinv_r * drop(dG_list[[k]] %*% Zt_Vinv_Z %*% dG_list[[l]] %*% Zt_Vinv_r))
        val <- -0.5 * (tr_term - 2 * quad_term)
        A_tt[k, l] <- val
        A_tt[l, k] <- val
      }
    }
    # cross G with sig2
    Zt_Vinv2_Z <- t(Z) %*% Vinv %*% Vinv %*% Z
    Zt_Vinv2_r <- t(Z) %*% Vinv %*% Vinv_r
    for (k in 1:3) {
      tr_term <- sum(Zt_Vinv2_Z * dG_list[[k]])
      quad_term <- sum(Zt_Vinv_r * drop(dG_list[[k]] %*% Zt_Vinv2_r))
      val <- -0.5 * (tr_term - 2 * quad_term)
      A_tt[k, 4] <- val
      A_tt[4, k] <- val
    }
    # sig2 with sig2
    tr_term <- sum(Vinv * Vinv)
    quad_term <- sum(Vinv_r * drop(Vinv %*% Vinv_r))
    A_tt[4, 4] <- -0.5 * (tr_term - 2 * quad_term)

    A11_i <- rbind(
      cbind(A_bb, A_b_theta),
      cbind(t(A_b_theta), A_tt)
    )
    A11 <- A11 + A11_i
  }

  A11 <- A11 / N

  list(s1_mat = s1_mat, a11_hat = A11)
}
