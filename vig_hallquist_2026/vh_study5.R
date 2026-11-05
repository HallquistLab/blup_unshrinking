# Study 5: matched-reliability calibration bridge.
#
# The data model and estimator bundle are deliberately identical to Study 2.
# Study 5 changes only how the first-stage population parameters are calibrated,
# so any arm differences can be attributed to reliability definition and
# covariance geometry rather than estimator drift.

study5_methods <- function() {
  study2_methods()
}

simulate_study5 <- function(condition) {
  simulate_data_blup_as_predictor(condition)
}

run_study5_rep <- function(condition) {
  # `run_study2_rep()` obtains its data through `simulate_study2()`, which is
  # the same shared BLUP-as-predictor generator used above. Passing a Study 5
  # condition preserves the Study 5 context in every result row.
  run_study2_rep(condition)
}
