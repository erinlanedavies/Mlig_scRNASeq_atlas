# Managing R packages with `renv`

Run these commands from the **Mlig_scRNASeq_atlas** project root. Commit `renv.lock` and `renv/activate.R`, but **not** the project package library.

## Load the project environment

Open the `.Rproj` file or start R in the project root. If `renv` has not activated automatically:

```r
source("renv/activate.R")
renv::status()
```

To install the package versions recorded in `renv.lock` (for example, on another machine):

```r
renv::restore()
```

`restore()` changes the local package library to match the lockfile; it does **not** update the lockfile.

## Add or update a package

```r
renv::install("clustree")                    # CRAN
renv::install("bioc::EnhancedVolcano")        # Bioconductor
# renv::install("owner/repository")           # GitHub
```

Reference the package in a project `.R` or `.Rmd` file, e.g. `library(EnhancedVolcano)` or `EnhancedVolcano::EnhancedVolcano(...)`. This project uses dependency discovery, so installed packages without detectable references may be omitted from a normal snapshot.

## Save package versions to `renv.lock`

```r
renv::dependencies(path = ".")   # Check detected dependencies and parsing warnings
renv::snapshot()                 # Review proposed changes before confirming
renv::status()                   # Check library / lockfile / project consistency
```

Fix any source-file parsing errors reported by `dependencies()`. Check the proposed snapshot for unexpected removals or additions before accepting it. `snapshot()` writes the lockfile; it does not save analysis objects or install packages.

If a package must be recorded even though it is not detected, use `renv::record()` cautiously, then inspect the resulting lockfile entry and run `renv::status()`. For Bioconductor packages, verify that the entry has `Source: Bioconductor` and the correct repository information. Avoid using `snapshot(packages = ...)` to modify the main lockfile without reviewing the diff: it may replace its package set rather than append to it.

## Commit and reproduce

```bash
git diff -- renv.lock
git add renv.lock
git commit -m "Update R dependencies"
```

On another machine, open the project and run `renv::restore()`.

## Troubleshooting

- **Installed and recorded, but `used = n`:** Run `renv::dependencies(path = ".")`. Check that the package is referenced in a scanned `.R`/`.Rmd` file and that the file parses successfully. Do not remove a needed package solely because it is marked unused.
- **`renv` version mismatch:** Compare `packageVersion("renv")` with `renv::lockfile_read()$Packages$renv$Version`. Install the intended version if necessary, then restart R.
- **Unexpected lockfile changes:** Inspect `git diff -- renv.lock` before committing. Restore the previous file from Git only if you intend to discard all uncommitted lockfile changes.

