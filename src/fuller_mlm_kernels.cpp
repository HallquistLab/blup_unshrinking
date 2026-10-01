// Fast row-wise kernels for the Fuller/M-MOM MLM estimator.
//
// The estimator's small dense matrix algebra remains in R, where crossprod(),
// solve(), eigen(), and geigen::geigen() already use compiled BLAS/LAPACK.
// These kernels only replace the R loops over subject-specific q by q
// measurement-error covariance matrices.
//
// These exported kernels are internal implementation helpers. The R fitter
// validates finite inputs and supplies canonical q by q by n arrays in R's
// column-major order: row + q * column + q * q * subject (zero-based indices).
// q excludes the error-free regression intercept. No kernel unshrinks BLUPs.

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include <algorithm>
#include <cmath>
#include <vector>

using namespace Rcpp;

namespace {

void check_omega_dimensions(
    const IntegerVector& dimensions,
    int expected_n,
    int expected_q) {
  if (dimensions.size() != 3 ||
      dimensions[0] != expected_q ||
      dimensions[1] != expected_q ||
      dimensions[2] != expected_n) {
    stop("Internal error: omega_x must be a q by q by n array.");
  }
}

bool matrix_is_psd(
    const arma::mat& matrix,
    double tolerance) {
  arma::vec eigenvalues;
  const bool succeeded = arma::eig_sym(eigenvalues, matrix);
  if (!succeeded || eigenvalues.n_elem == 0) return false;
  const double min_eigen = eigenvalues.min();
  const double max_eigen = eigenvalues.max();
  const double scale = std::max(1.0, std::abs(max_eigen));
  return std::isfinite(min_eigen) && min_eigen >= -tolerance * scale;
}

}  // namespace

// Symmetrize a common q by q matrix or canonical q by q by n array, then
// return the canonical array and bad_row (0 if valid, otherwise 1-based).
// The caller must reject bad_row > 0 before using the partially filled array.
// [[Rcpp::export]]
List fuller_mlm_normalize_omega_x_cpp(
    NumericVector omega_x,
    int n,
    int q,
    double tolerance = 1e-8) {
  const IntegerVector dimensions = omega_x.attr("dim");
  const bool common_matrix = dimensions.size() == 2;
  if (common_matrix) {
    if (dimensions[0] != q || dimensions[1] != q) {
      stop("Internal error: common omega_x must be a q by q matrix.");
    }
  } else {
    check_omega_dimensions(dimensions, n, q);
  }

  NumericVector normalized = common_matrix ?
    NumericVector(static_cast<R_xlen_t>(q) * q * n) : clone(omega_x);
  arma::mat omega_i(q, q);
  int bad_row = 0;

  const int matrices_to_validate = common_matrix ? 1 : n;
  for (int i = 0; i < matrices_to_validate; ++i) {
    const R_xlen_t source_offset = common_matrix ?
      0 : static_cast<R_xlen_t>(i) * q * q;
    for (int column = 0; column < q; ++column) {
      for (int row = 0; row < q; ++row) {
        const double value = 0.5 * (
          omega_x[source_offset + row + q * column] +
          omega_x[source_offset + column + q * row]
        );
        omega_i(row, column) = value;
        if (!common_matrix) {
          normalized[source_offset + row + q * column] = value;
        }
      }
    }
    if (!matrix_is_psd(omega_i, tolerance)) {
      bad_row = i + 1;
      break;
    }
  }

  if (common_matrix && bad_row == 0) {
    for (int i = 0; i < n; ++i) {
      const R_xlen_t destination_offset = static_cast<R_xlen_t>(i) * q * q;
      for (int column = 0; column < q; ++column) {
        for (int row = 0; row < q; ++row) {
          normalized[destination_offset + row + q * column] =
            omega_i(row, column);
        }
      }
    }
  }

  normalized.attr("dim") = IntegerVector::create(q, q, n);
  return List::create(
    _["omega_x"] = normalized,
    _["bad_row"] = bad_row
  );
}

// Validate the full outcome/predictor measurement-error block for every row.
// Return 0 if all blocks are PSD within tolerance, or the first bad row (1-based).
// [[Rcpp::export]]
int fuller_mlm_joint_omega_bad_row_cpp(
    const NumericVector& omega_x,
    const NumericVector& omega_y,
    const NumericMatrix& omega_xy,
    double tolerance = 1e-8) {
  const IntegerVector dimensions = omega_x.attr("dim");
  if (dimensions.size() != 3) {
    stop("Internal error: omega_x must be a q by q by n array.");
  }
  const int q = dimensions[0];
  const int n = dimensions[2];
  check_omega_dimensions(dimensions, n, q);
  if (omega_y.size() != n || omega_xy.nrow() != n || omega_xy.ncol() != q) {
    stop("Internal error: joint measurement-error inputs are not aligned.");
  }

  arma::mat joint(q + 1, q + 1, arma::fill::zeros);
  for (int i = 0; i < n; ++i) {
    joint.zeros();
    joint(0, 0) = omega_y[i];
    const R_xlen_t offset = static_cast<R_xlen_t>(i) * q * q;
    for (int column = 0; column < q; ++column) {
      joint(0, column + 1) = omega_xy(i, column);
      joint(column + 1, 0) = omega_xy(i, column);
      for (int row = 0; row < q; ++row) {
        joint(row + 1, column + 1) =
          omega_x[offset + row + q * column];
      }
    }
    if (!matrix_is_psd(joint, tolerance)) return i + 1;
  }
  return 0;
}

// Sum_i weights[i] * omega_x[, , i].
// [[Rcpp::export]]
NumericMatrix fuller_mlm_weighted_omega_sum_cpp(
    const NumericVector& omega_x,
    const NumericVector& weights) {
  const IntegerVector dimensions = omega_x.attr("dim");
  if (dimensions.size() != 3) {
    stop("Internal error: omega_x must be a q by q by n array.");
  }
  const int q = dimensions[0];
  const int n = dimensions[2];
  check_omega_dimensions(dimensions, n, q);
  if (weights.size() != n) {
    stop("Internal error: omega_x and weights are not aligned.");
  }

  NumericMatrix result(q, q);
  for (int i = 0; i < n; ++i) {
    const double weight = weights[i];
    const R_xlen_t offset = static_cast<R_xlen_t>(i) * q * q;
    for (int column = 0; column < q; ++column) {
      for (int row = 0; row < q; ++row) {
        result(row, column) +=
          weight * omega_x[offset + row + q * column];
      }
    }
  }
  return result;
}

// Return the n by q matrix whose ith row is omega_x[, , i] %*% beta.
// [[Rcpp::export]]
NumericMatrix fuller_mlm_omega_matvec_cpp(
    const NumericVector& omega_x,
    const NumericVector& beta) {
  const IntegerVector dimensions = omega_x.attr("dim");
  if (dimensions.size() != 3) {
    stop("Internal error: omega_x must be a q by q by n array.");
  }
  const int q = dimensions[0];
  const int n = dimensions[2];
  check_omega_dimensions(dimensions, n, q);
  if (beta.size() != q) {
    stop("Internal error: omega_x and beta are not conformable.");
  }

  NumericMatrix result(n, q);
  for (int i = 0; i < n; ++i) {
    const R_xlen_t offset = static_cast<R_xlen_t>(i) * q * q;
    for (int row = 0; row < q; ++row) {
      double value = 0.0;
      for (int column = 0; column < q; ++column) {
        value += omega_x[offset + row + q * column] * beta[column];
      }
      result(i, row) = value;
    }
  }
  return result;
}

// For t_i = (0, omega_xy[i, ] - Omega_x[i] * gamma0_x), return the p by p
// sum of inverse_case_variance[i]^2 * t_i * t_i', with p = q + 1. The leading
// row/column is zero because the intercept has no measurement error. gamma0_x
// contains preliminary slopes, not the final coefficient estimates.
// [[Rcpp::export]]
NumericMatrix fuller_mlm_tilde_crossproduct_cpp(
    const NumericVector& omega_x,
    const NumericMatrix& omega_xy,
    const NumericVector& gamma0_x,
    const NumericVector& inverse_case_variance) {
  const int n = omega_xy.nrow();
  const int q = omega_xy.ncol();
  const IntegerVector dimensions = omega_x.attr("dim");
  check_omega_dimensions(dimensions, n, q);
  if (gamma0_x.size() != q || inverse_case_variance.size() != n) {
    stop("Internal error: variance-kernel inputs are not aligned.");
  }
  NumericMatrix result(q + 1, q + 1);
  std::vector<double> tilde(q);
  for (int i = 0; i < n; ++i) {
    const R_xlen_t offset = static_cast<R_xlen_t>(i) * q * q;
    for (int row = 0; row < q; ++row) {
      double value = omega_xy(i, row);
      for (int column = 0; column < q; ++column) {
        value -= omega_x[offset + row + q * column] * gamma0_x[column];
      }
      tilde[row] = value * inverse_case_variance[i];
    }
    for (int column = 0; column < q; ++column) {
      for (int row = 0; row < q; ++row) {
        result(row + 1, column + 1) += tilde[row] * tilde[column];
      }
    }
  }
  return result;
}
