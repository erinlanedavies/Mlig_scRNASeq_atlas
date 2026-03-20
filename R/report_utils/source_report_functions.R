# R/source_report_functions.R

source_report_functions <- function(base_dir = NULL,
                                    rel_dir = "R/report_utils",
                                    files = c("io.R", "qc.R", "metadata.R")) {

  if (is.null(base_dir) || !nzchar(base_dir)) {
    base_dir <- getwd()
  }

  base_dir <- normalizePath(base_dir, mustWork = TRUE)
  utils_dir <- file.path(base_dir, rel_dir)

  if (!dir.exists(utils_dir)) {
    stop("Functions directory not found: ", utils_dir)
  }

  for (f in files) {
    f_path <- file.path(utils_dir, f)
    if (!file.exists(f_path)) {
      stop("Required functions file not found: ", f_path)
    }
    source(f_path)
  }

  invisible(TRUE)
}

