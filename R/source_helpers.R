# Central source order for repository-level helper modules.
#
# Simulation entry points source helpers into their own execution environment
# instead of using an installed package namespace. Keep the order here as the
# single place new scripts and tests should consult.

project_helper_source_order <- function() {
  c(
    "core_utils.R",
    "simulation_runner_helpers.R",
    "derivative_backends.R",
    "tmb_stage1_helpers.R",
    "stacked_sandwich_helpers.R",
    "sim_helpers.R",
    "reliability_calibration.R",
    "stats_helpers.R",
    "stage2_estimators.R",
    "fuller_mlm.R",
    "blup_helpers.R",
    "sim_diagnostics.R",
    "lai_openmx_helpers.R",
    "mplus_helpers.R"
  )
}

source_project_helpers <- function(repo_root = ".",
                                   helper_files = project_helper_source_order(),
                                   envir = parent.frame()) {
  helper_paths <- file.path(repo_root, "R", helper_files)
  missing <- helper_paths[!file.exists(helper_paths)]

  if (length(missing) > 0L) {
    stop("Missing project helper file(s): ", paste(missing, collapse = ", "))
  }

  for (helper_path in helper_paths) {
    source(helper_path, local = envir)
  }

  invisible(helper_paths)
}
