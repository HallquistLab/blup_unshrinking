#!/usr/bin/env Rscript

#' ---
#' title: "REPAIR: Random-effect prediction adjustment for inference in regression"
#' description: "Legacy and amended simulation studies including: (1) RE as outcome,
#'               (2) RE as predictor, (3) RE as predictor and outcome,
#'               (4) heterogeneous cluster information and reliability,
#'               (5) matched reliability calibration bridge, and an
#'               ICC-versus-posterior-reliability bridge"
#' ---

# TODO: check these exports
suppressPackageStartupMessages({
  library(data.table)
  library(lme4)
  library(MASS)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(foreach)
  library(doParallel)
  library(OpenMx)
  library(glmnet)
  library(sandwich)
  library(geigen)
  library(MplusAutomation)
  library(glue)
})

current_script_path <- function() {
  frame_files <- vapply(sys.frames(), function(frame) {
    if (!is.null(frame$ofile)) frame$ofile else NA_character_
  }, character(1))
  frame_files <- stats::na.omit(frame_files)
  if (length(frame_files) > 0L) {
    return(normalizePath(frame_files[[length(frame_files)]], mustWork = FALSE))
  }

  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0L) {
    return(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE))
  }

  NA_character_
}

find_repo_root <- function() {
  script_path <- current_script_path()
  script_dir <- if (!is.na(script_path)) dirname(script_path) else getwd()
  candidates <- unique(normalizePath(c(
    getwd(),
    file.path(getwd(), ".."),
    script_dir,
    file.path(script_dir, ".."),
    file.path(script_dir, "..", "..")
  ), mustWork = FALSE))

  for (candidate in candidates) {
    # these are just example files that should exist in the repo root
    if (file.exists(file.path(candidate, "R", "core_utils.R")) &&
      file.exists(file.path(candidate, "R", "lai_openmx_helpers.R"))) {
      return(candidate)
    }
  }

  stop("Could not locate repository root containing shared R/ helpers and vig_hallquist_2026/ scripts.")
}

repo_root <- find_repo_root()
vig_hallquist_dir <- file.path(repo_root, "vig_hallquist_2026")

source(file.path(repo_root, "R", "source_helpers.R"), local = TRUE)
source_project_helpers(repo_root)

source(file.path(vig_hallquist_dir, "vh_study_common.R"), local = TRUE)
source(file.path(vig_hallquist_dir, "vh_designs.R"), local = TRUE)
source(file.path(vig_hallquist_dir, "vh_study1.R"), local = TRUE)
source(file.path(vig_hallquist_dir, "vh_study2.R"), local = TRUE)
for (study_file in c("vh_study3.R", "vh_study4.R", "vh_study5.R")) {
  study_path <- file.path(vig_hallquist_dir, study_file)
  if (file.exists(study_path)) {
    source(study_path, local = TRUE)
  }
}
source(file.path(vig_hallquist_dir, "vh_runner.R"), local = TRUE)

Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
OpenMx::mxOption(NULL, "Number of Threads", 1L)

args <- commandArgs(trailingOnly = TRUE)

if (sys.nframe() == 0L) {
  n_sim <- if (length(args) >= 1L) as.integer(args[[1]]) else 100L
  study_arg <- if (length(args) >= 2L) args[[2]] else "all"
  out_dir <- if (length(args) >= 3L) args[[3]] else file.path(vig_hallquist_dir, "outputs", "vig_hallquist")
  n_cores <- if (length(args) >= 4L) as.integer(args[[4]]) else 1L
  max_conditions <- if (length(args) >= 5L) parse_optional_integer_arg(args[[5]]) else NA_integer_
  chunk_index <- if (length(args) >= 6L) parse_optional_integer_arg(args[[6]]) else NA_integer_
  chunk_size <- if (length(args) >= 7L) parse_optional_integer_arg(args[[7]]) else NA_integer_
  resume_existing <- if (length(args) >= 8L) parse_logical_arg(args[[8]], default = TRUE) else TRUE
  max_aggregate_replication_rows <- as.numeric(Sys.getenv(
    "VIG_HALLQUIST_MAX_AGGREGATE_ROWS",
    "2000000"
  ))

  run_simulation(
    n_sim = n_sim,
    study_arg = study_arg,
    out_dir = out_dir,
    n_cores = n_cores,
    max_conditions = max_conditions,
    chunk_index = chunk_index,
    chunk_size = chunk_size,
    resume_existing = resume_existing,
    max_aggregate_replication_rows = max_aggregate_replication_rows
  )
}
