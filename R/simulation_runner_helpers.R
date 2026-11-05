#' Shared condition-runner utilities for script-oriented simulations.
#'
#' These helpers keep chunking, progress files, and atomic CSV writes
#' consistent across the larger simulation drivers without imposing a package
#' structure on the repository.

parse_optional_integer_arg <- function(x) {
  if (is.null(x) || length(x) == 0L || is.na(x) || x %in% c("", "NA", "NaN", "NULL", "null")) {
    return(NA_integer_)
  }
  as.integer(x)
}

parse_logical_arg <- function(x, default = TRUE) {
  if (is.null(x) || length(x) == 0L || is.na(x) || x %in% c("", "NA", "NaN", "NULL", "null")) {
    return(default)
  }
  as.logical(as.integer(x))
}

slice_condition_chunk <- function(condition_grid, chunk_index = NA_integer_, chunk_size = NA_integer_) {
  if (is.na(chunk_index) && is.na(chunk_size)) {
    attr(condition_grid, "chunk_meta") <- list(
      chunk_index = NA_integer_,
      chunk_size = NA_integer_,
      condition_start = min(condition_grid$condition_id),
      condition_end = max(condition_grid$condition_id),
      n_conditions = nrow(condition_grid)
    )
    return(condition_grid)
  }

  if (is.na(chunk_index) || is.na(chunk_size)) {
    stop("`chunk_index` and `chunk_size` must be supplied together.")
  }
  if (chunk_index < 1L || chunk_size < 1L) {
    stop("`chunk_index` and `chunk_size` must be positive integers.")
  }

  start_idx <- ((chunk_index - 1L) * chunk_size) + 1L
  end_idx <- min(nrow(condition_grid), chunk_index * chunk_size)
  if (start_idx > nrow(condition_grid)) {
    stop("Requested chunk starts after the end of the selected condition grid.")
  }

  out <- condition_grid %>%
    dplyr::slice(start_idx:end_idx)

  attr(out, "chunk_meta") <- list(
    chunk_index = chunk_index,
    chunk_size = chunk_size,
    condition_start = min(out$condition_id),
    condition_end = max(out$condition_id),
    n_conditions = nrow(out)
  )

  out
}

make_chunk_label <- function(chunk_meta) {
  if (is.null(chunk_meta) || is.na(chunk_meta$chunk_index)) {
    "full_selection"
  } else {
    sprintf(
      "chunk_%03d_conditions_%04d_%04d",
      chunk_meta$chunk_index,
      chunk_meta$condition_start,
      chunk_meta$condition_end
    )
  }
}

write_csv_atomic <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp_path <- paste0(path, ".tmp.", Sys.getpid())
  utils::write.csv(x, file = tmp_path, row.names = FALSE)
  if (!file.rename(tmp_path, path)) {
    file.copy(tmp_path, path, overwrite = TRUE)
    unlink(tmp_path)
  }
  invisible(path)
}

write_progress_row <- function(progress_path, row_df) {
  dir.create(dirname(progress_path), recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(progress_path)) {
    write_csv_atomic(row_df, progress_path)
  } else {
    utils::write.table(
      row_df,
      file = progress_path,
      sep = ",",
      row.names = FALSE,
      col.names = FALSE,
      append = TRUE
    )
  }
  invisible(progress_path)
}

condition_file_paths <- function(out_dir, condition_id, prefix = "condition", compressed_replications = FALSE) {
  condition_dir <- file.path(out_dir, "conditions")
  dir.create(condition_dir, recursive = TRUE, showWarnings = FALSE)

  rep_ext <- if (isTRUE(compressed_replications)) "csv.gz" else "csv"
  list(
    replications = file.path(condition_dir, sprintf("%s_%04d_replication_results.%s", prefix, condition_id, rep_ext)),
    summary = file.path(condition_dir, sprintf("%s_%04d_summary.csv", prefix, condition_id))
  )
}

# Default numeric whitelist: columns we should attempt to coerce to numeric on read.
# This helps with extremely large numbers in previous runs; current runs have a 
# screen for this issue
default_numeric_whitelist <- c("stage1_design_kappa")

# Try to coerce one or more columns to numeric if all non-missing entries parse as numeric.
coerce_numeric_if_possible <- function(df, cols) {
  cols <- intersect(as.character(cols), names(df))
  for (col in cols) {
    x <- df[[col]]
    if (is.numeric(x)) {
      next
    }
    if (all(is.na(x) | suppressWarnings(!is.na(as.numeric(x))))) {
      df[[col]] <- as.numeric(x)
    }
  }
  return(df)
}

read_condition_results_file <- function(path, numeric_whitelist = default_numeric_whitelist) {
  if (requireNamespace("data.table", quietly = TRUE)) {
    out <- tibble::as_tibble(data.table::fread(path))
  } else {
    out <- utils::read.csv(path, stringsAsFactors = FALSE)
  }

  for (col in numeric_whitelist) {
    out <- coerce_numeric_if_possible(out, default_numeric_whitelist)
  }

  out
}

load_completed_condition_results <- function(condition_grid,
                                             out_dir,
                                             path_fun = condition_file_paths,
                                             ...) {
  files <- vapply(
    condition_grid$condition_id,
    function(condition_id) {
      paths <- path_fun(out_dir, condition_id, ...)
      gz_path <- paths$replications
      csv_path <- sub("\\.gz$", "", gz_path)
      if (file.exists(gz_path)) {
        gz_path
      } else if (file.exists(csv_path)) {
        csv_path
      } else {
        NA_character_
      }
    },
    character(1)
  )

  existing <- stats::na.omit(files)
  if (length(existing) == 0L) {
    return(tibble::tibble())
  }

  purrr::map_dfr(existing, read_condition_results_file)
}
