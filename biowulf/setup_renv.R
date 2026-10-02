#!/usr/bin/env Rscript

# ==============================================================================
# setup_renv.R
#
# Reproduce the R environment recorded in renv.lock.
#
# This script:
#   1. validates that it is being run from the project
#   2. activates renv
#   3. restores packages exactly from renv.lock
#   4. checks that the resulting environment is consistent
#
# It intentionally does NOT:
#   - update packages
#   - snapshot the environment
#   - modify renv.lock
# ==============================================================================

options(
  repos = c(CRAN = "https://cloud.r-project.org")
)

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

lockfile <- file.path(project_root, "renv.lock")
activate <- file.path(project_root, "renv", "activate.R")

if (!file.exists(lockfile)) {
  stop("renv.lock not found: ", lockfile)
}

if (!file.exists(activate)) {
  stop("renv activation script not found: ", activate)
}


# ------------------------------------------------------------------------------
# Activate renv
# ------------------------------------------------------------------------------

cat("\n=== Activating renv ===\n")

source(activate)

cat("renv version: ", as.character(packageVersion("renv")), "\n", sep = "")
cat("R version:    ", R.version.string, "\n", sep = "")

cat("\nLibrary paths:\n")
print(.libPaths())


# ------------------------------------------------------------------------------
# Restore
# ------------------------------------------------------------------------------

cat("\n=== Restoring R environment ===\n\n")

renv::restore(
  project = project_root,
  lockfile = lockfile,
  prompt = FALSE
)


# ------------------------------------------------------------------------------
# Validate
# ------------------------------------------------------------------------------

cat("\n=== Checking renv status ===\n\n")

status <- renv::status(
  project = project_root
)

if (!status$synchronized) {
  stop(
    "\nrenv restore completed, but the project is not synchronized ",
    "with renv.lock."
  )
}


# ------------------------------------------------------------------------------
# Finished
# ------------------------------------------------------------------------------

cat("\n============================================================\n")
cat("R environment restored successfully.\n")
cat("============================================================\n")
