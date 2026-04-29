#!/usr/bin/env Rscript

options(repos = c(CRAN = "https://cloud.r-project.org"))

LIBRARY_ROOT <- Sys.getenv(
  "RENV_PATHS_LIBRARY_ROOT",
  file.path("/data", Sys.getenv("USER"), "renv_libs")
)

CACHE_ROOT <- Sys.getenv(
  "RENV_PATHS_CACHE",
  file.path("/data", Sys.getenv("USER"), "renv_cache")
)

PACKAGES_FILE <- "packages.txt"
UPDATE_ALL_DEPENDENCIES <- TRUE

ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  if (!dir.exists(path)) {
    stop("Could not create directory: ", path)
  }
}

message2 <- function(...) cat(..., "\n", sep = "")

read_package_file <- function(path) {
  if (!file.exists(path)) {
    stop("Package file not found: ", path)
  }

  pkgs <- readLines(path, warn = FALSE)
  pkgs <- sub("#.*$", "", pkgs)
  pkgs <- trimws(pkgs)
  pkgs <- pkgs[nzchar(pkgs)]

  unique(pkgs)
}

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

library_root <- path.expand(LIBRARY_ROOT)
cache_root <- path.expand(CACHE_ROOT)

ensure_dir(library_root)
ensure_dir(cache_root)

Sys.setenv(
  RENV_PATHS_LIBRARY_ROOT = library_root,
  RENV_PATHS_CACHE = cache_root
)

if (nzchar(Sys.getenv("RETICULATE_PYTHON"))) {
  message2("Using RETICULATE_PYTHON: ", Sys.getenv("RETICULATE_PYTHON"))
}

if (nzchar(Sys.getenv("LD_LIBRARY_PATH"))) {
  message2("Using LD_LIBRARY_PATH: ", Sys.getenv("LD_LIBRARY_PATH"))
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

options(repos = BiocManager::repositories())

if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv", repos = "https://cloud.r-project.org")
}

library(renv)

options(repos = BiocManager::repositories())

if (!file.exists(file.path(project_root, "renv", "activate.R"))) {
  message2("Initializing renv...")
  renv::init(project = project_root, bare = TRUE, restart = FALSE)
} else {
  message2("Loading existing renv project...")
  renv::load(project = project_root)
}

options(repos = BiocManager::repositories())

renv::settings$snapshot.type("explicit", project = project_root)

renviron_lines <- c(
  paste0("RENV_PATHS_LIBRARY_ROOT=", library_root),
  paste0("RENV_PATHS_CACHE=", cache_root),
  paste0("RETICULATE_PYTHON=", Sys.getenv("RETICULATE_PYTHON")),
  paste0("RETICULATE_PYENV_ROOT=", Sys.getenv("RETICULATE_PYENV_ROOT")),
  paste0("BASILISK_DIR=", Sys.getenv("BASILISK_DIR"))
)

renviron_lines <- renviron_lines[nzchar(sub("^[^=]+=", "", renviron_lines))]

renviron_path <- file.path(project_root, ".Renviron")

if (file.exists(renviron_path)) {
  old <- readLines(renviron_path, warn = FALSE)
  old <- old[!grepl("^RENV_PATHS_LIBRARY_ROOT=", old)]
  old <- old[!grepl("^RENV_PATHS_CACHE=", old)]
  old <- old[!grepl("^RETICULATE_PYTHON=", old)]
  old <- old[!grepl("^RETICULATE_PYENV_ROOT=", old)]
  old <- old[!grepl("^BASILISK_DIR=", old)]
  writeLines(c(old, renviron_lines), renviron_path)
} else {
  writeLines(renviron_lines, renviron_path)
}

PROJECT_PACKAGES <- read_package_file(PACKAGES_FILE)

message2("Installing packages:")
message2(paste(PROJECT_PACKAGES, collapse = ", "))

renv::install(PROJECT_PACKAGES)

options(repos = BiocManager::repositories())

if (isTRUE(UPDATE_ALL_DEPENDENCIES)) {
  message2("Updating all out-of-date packages...")
  renv::update(prompt = FALSE)
} else {
  message2("Updating packages listed in packages.txt...")
  renv::update(packages = PROJECT_PACKAGES, prompt = FALSE)
}

options(repos = BiocManager::repositories())

message2("Snapshotting lockfile...")

renv::snapshot(
  project = project_root,
  packages = PROJECT_PACKAGES,
  prompt = FALSE
)

message2("Checking Bioconductor validity...")
BiocManager::valid()

message2("Done.")
