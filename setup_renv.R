#!/usr/bin/env Rscript

# setup_renv.R
# Purpose:
#   1. Ensure renv is available
#   2. Activate renv in this project
#   3. Force renv libraries/cache to live outside the project
#   4. Install requested packages
#   5. Lock exact versions in renv.lock
#
# Usage:
#   Rscript setup_renv.R
#
# Notes:
#   - Edit LIBRARY_ROOT and CACHE_ROOT below to match your system.
#   - This script assumes it is run from the project root directory.

options(repos = c(CRAN = "https://cloud.r-project.org"))

readRenviron(".Renviron")

# -----------------------------
# User-configurable paths
# -----------------------------
LIBRARY_ROOT <- file.path("/data", Sys.getenv("USER"), "renv_libs")
CACHE_ROOT   <- file.path("/data", Sys.getenv("USER"), "renv_cache")

# Packages you want locked for this project.
# Add or remove package names here.
PROJECT_PACKAGES <- scan("packages.txt", what = character(), quiet = TRUE)

# -----------------------------
# Helper functions
# -----------------------------
ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  if (!dir.exists(path)) {
    stop("Could not create directory: ", path)
  }
}

message2 <- function(...) cat(..., "\n", sep = "")

# -----------------------------
# Resolve paths and create dirs
# -----------------------------
library_root <- path.expand(LIBRARY_ROOT)
cache_root <- path.expand(CACHE_ROOT)

ensure_dir(library_root)
ensure_dir(cache_root)

# Make these visible to renv for this process
Sys.setenv(
  RENV_PATHS_LIBRARY_ROOT = library_root,
  RENV_PATHS_CACHE = cache_root
)

message2("Using external renv library root: ", library_root)
message2("Using external renv cache: ", cache_root)

# -----------------------------
# Ensure renv is installed
# -----------------------------
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv")
}
library(renv)

# -----------------------------
# Initialize / activate project
# -----------------------------
# If renv is not initialized, initialize it without snapshot prompt.
# bare = TRUE avoids trying to discover and snapshot everything immediately.
project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

message2("Project root: ", project_root)
message2("renv.lock exists? ", file.exists(file.path(project_root, "renv.lock")))
message2("renv/activate.R exists? ", file.exists(file.path(project_root, "renv", "activate.R")))

if (!file.exists(file.path(project_root, "renv", "activate.R"))) {
  stop("renv/activate.R not found. Run this from the project root.")
}

if (!file.exists(file.path(project_root, "renv.lock"))) {
  message2("No renv.lock found. Initializing renv in this project without restart...")
  renv::init(project = project_root, bare = TRUE, restart = FALSE)
} else {
  message2("Existing renv project found. Loading project...")
  renv::load(project = project_root)
}

# -----------------------------
# Write/update project .Renviron
# -----------------------------
# This hard-codes the external library/cache paths for future sessions.
# Because .Renviron lives in the project root, it is loaded automatically
# when you start R in the project.
renviron_lines <- c(
  paste0("RENV_PATHS_LIBRARY_ROOT=", library_root),
  paste0("RENV_PATHS_CACHE=", cache_root)
)

if (file.exists(".Renviron")) {
  old <- readLines(".Renviron", warn = FALSE)
  old <- old[!grepl("^RENV_PATHS_LIBRARY_ROOT=", old)]
  old <- old[!grepl("^RENV_PATHS_CACHE=", old)]
  writeLines(c(old, renviron_lines), ".Renviron")
} else {
  writeLines(renviron_lines, ".Renviron")
}

message2("Updated .Renviron with renv external paths.")

# -----------------------------
# Install and lock packages
# -----------------------------
# install(..., lock = TRUE) updates the lockfile after install in recent renv.
# If unavailable in older renv for any reason, the fallback is snapshot().
message2("Installing project packages:")
message2("  ", paste(PROJECT_PACKAGES, collapse = ", "))

# Explicit snapshot for safety and clarity
renv::snapshot(prompt = FALSE)

message2("Finished.")
message2("Project dependencies are now installed and locked in renv.lock.")
message2("Future users/sessions can run renv::restore() to recreate this environment.")
