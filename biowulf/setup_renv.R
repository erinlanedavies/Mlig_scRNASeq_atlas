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

cat("\n=== Validating restored R environment ===\n\n")

lock <- renv::lockfile_read(lockfile)

locked <- vapply(
    lock$Packages,
    function(pkg) pkg$Version,
    character(1)
)

installed <- installed.packages()
installed_versions <- installed[, "Version"]

missing <- setdiff(names(locked), names(installed_versions))

common <- intersect(names(locked), names(installed_versions))

wrong_version <- common[
    installed_versions[common] != locked[common]
]

if (length(missing)) {
    cat("\nMissing locked packages:\n")
    print(missing)
}

if (length(wrong_version)) {
    cat("\nPackages with incorrect versions:\n")

    for (pkg in wrong_version) {
        cat(
            sprintf(
                "  %-25s installed: %-15s expected: %s\n",
                pkg,
                installed_versions[[pkg]],
                locked[[pkg]]
            )
        )
    }
}

if (length(missing) || length(wrong_version)) {
    stop(
        "\nR environment restoration failed validation: ",
        "the installed library does not satisfy renv.lock."
    )
}

extra <- setdiff(names(installed_versions), names(locked))

if (length(extra)) {
    cat(
        "\nNote: additional packages were installed during dependency ",
        "resolution:\n"
    )
    cat("  ", paste(extra, collapse = ", "), "\n", sep = "")
}

cat(
    "\nAll packages recorded in renv.lock are installed ",
    "at the expected versions.\n"
)

cat("\n============================================================\n")
cat("R environment restored successfully.\n")
cat("============================================================\n")

# ------------------------------------------------------------------------------
# Finished
# ------------------------------------------------------------------------------

cat("\n============================================================\n")
cat("R environment restored successfully.\n")
cat("============================================================\n")
