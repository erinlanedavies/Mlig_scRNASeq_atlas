#!/usr/bin/env bash

# ==============================================================================
# run_r.sh
#
# Run R scripts and R Markdown reports for the Mlig_scRNASeq_atlas project
# on NIH Biowulf.
#
# The launcher configures the runtime environment:
#   - Biowulf modules
#   - project renv library/cache
#   - Python interpreter used by reticulate
#   - writable cache location
#
# Runtime paths have production defaults but can be overridden with environment
# variables. This allows the same launcher to be used for clean-room testing.
#
# Production:
#   ./biowulf/run_r.sh scripts/KnitReports.R
#
# Clean-room example:
#   ENV_PREFIX=/data/$USER/conda_envs/mlig_scrna_runtime_fresh_test \
#   RENV_PATHS_LIBRARY_ROOT=/vf/users/$USER/mlig_fresh_test/renv_libs \
#   RENV_PATHS_CACHE=/vf/users/$USER/mlig_fresh_test/renv_cache \
#   ./biowulf/run_r.sh scripts/KnitReports.R
#
# ==============================================================================

set -euo pipefail


# ------------------------------------------------------------------------------
# Project location
# ------------------------------------------------------------------------------

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$PROJECT_ROOT"


# ------------------------------------------------------------------------------
# Runtime environment locations
# ------------------------------------------------------------------------------

# Each value can be overridden by an existing environment variable.
# Otherwise, use the normal production location.

ENV_PREFIX="${ENV_PREFIX:-/data/${USER}/conda_envs/mlig_scrna_runtime}"
RENV_PATHS_LIBRARY_ROOT="${RENV_PATHS_LIBRARY_ROOT:-/vf/users/${USER}/renv_libs}"
RENV_PATHS_CACHE="${RENV_PATHS_CACHE:-/vf/users/${USER}/renv_cache}"
XDG_CACHE_HOME="${XDG_CACHE_HOME:-/data/${USER}/.cache}"

export RENV_PATHS_LIBRARY_ROOT
export RENV_PATHS_CACHE
export XDG_CACHE_HOME


# ------------------------------------------------------------------------------
# Validate command-line arguments
# ------------------------------------------------------------------------------

if [[ $# -ne 1 ]]; then
    echo "Usage:"
    echo "  $0 path/to/script.R"
    echo "  $0 path/to/report.Rmd"
    exit 1
fi

TARGET="$1"

if [[ ! -f "$TARGET" ]]; then
    echo "ERROR: File not found: $TARGET" >&2
    exit 1
fi


# ------------------------------------------------------------------------------
# Load Biowulf modules
# ------------------------------------------------------------------------------

# Start from a predictable module environment.
module purge

# Compiler and native libraries.
module load gcc/11.3.0
module load hdf5/1.12.2
module load netcdf/4.9.0_gcc-11.3.0
module load openmpi/5.0.5/gcc-11.3.0

# Libraries required by R packages.
module load pcre2/10.40_gcc-11.3.0
module load libtiff/4.6.0_gcc-11.3.0
module load libwebp/1.6.0-gcc-11.3.0

# R.
module load R/4.5.2

# Document-generation tools.
module load pandoc/2.18
module load tex/2024


# ------------------------------------------------------------------------------
# Validate R
# ------------------------------------------------------------------------------

if ! command -v Rscript >/dev/null 2>&1; then
    echo "ERROR: Rscript was not found after loading the R module." >&2
    exit 1
fi


# ------------------------------------------------------------------------------
# Native-library discovery
# ------------------------------------------------------------------------------

# Preserve PKG_CONFIG_PATH and LD_LIBRARY_PATH established by the Biowulf
# modules. Do not add paths from the Conda environment.

PKG_CONFIG_BIN="$(command -v pkg-config || true)"

if [[ "$PKG_CONFIG_BIN" != "/usr/bin/pkg-config" ]]; then
    echo "ERROR: Expected /usr/bin/pkg-config but found:" >&2
    echo "  ${PKG_CONFIG_BIN:-not found}" >&2
    exit 1
fi

if ! pkg-config --exists libwebp; then
    echo "ERROR: pkg-config cannot locate libwebp." >&2
    exit 1
fi

if ! pkg-config --exists libwebpmux; then
    echo "ERROR: pkg-config cannot locate libwebpmux." >&2
    exit 1
fi


# ------------------------------------------------------------------------------
# Python / reticulate
# ------------------------------------------------------------------------------

# Do not activate the Conda environment. Point reticulate directly to its
# Python executable instead.

PYTHON="${ENV_PREFIX}/bin/python"

if [[ ! -x "$PYTHON" ]]; then
    echo "ERROR: Python environment not found:" >&2
    echo "  $ENV_PREFIX" >&2
    echo >&2
    echo "Create the Python environment before running the pipeline." >&2
    exit 1
fi

export RETICULATE_PYTHON="$PYTHON"


# ------------------------------------------------------------------------------
# Cache
# ------------------------------------------------------------------------------

mkdir -p "$XDG_CACHE_HOME"


# ------------------------------------------------------------------------------
# Temporary files
# ------------------------------------------------------------------------------

# Biowulf normally sets TMPDIR=/lscratch/$SLURM_JOB_ID for allocated jobs.
# Preserve that setting.

if [[ -n "${TMPDIR:-}" ]]; then
    mkdir -p "$TMPDIR"
fi


# ------------------------------------------------------------------------------
# Runtime information
# ------------------------------------------------------------------------------

# Print the important environment choices so that job logs document exactly
# which reproducible environments were used.

echo
echo "=== Mlig_scRNASeq_atlas runtime ==="
echo "Project:            $PROJECT_ROOT"
echo "R:                  $(command -v R)"
echo "Python:             $RETICULATE_PYTHON"
echo "renv library root:  $RENV_PATHS_LIBRARY_ROOT"
echo "renv cache:         $RENV_PATHS_CACHE"
echo "Cache:              $XDG_CACHE_HOME"
echo "Target:             $TARGET"
echo "=================================="
echo


# ------------------------------------------------------------------------------
# Run target
# ------------------------------------------------------------------------------

case "$TARGET" in

    *.R)
        echo "Running R script: $TARGET"

        # --vanilla skips .Rprofile, so explicitly activate the project's
        # renv environment before sourcing the script.
        Rscript --vanilla -e \
            "renv::load(project = '$PROJECT_ROOT'); source('$TARGET', chdir = FALSE)"
        ;;

    *.Rmd)
        echo "Rendering R Markdown: $TARGET"

        # render_rmd.R explicitly activates the project's renv environment.
        Rscript --vanilla \
            "${PROJECT_ROOT}/biowulf/render_rmd.R" \
            "$TARGET"
        ;;

    *)
        echo "ERROR: Unsupported file type: $TARGET" >&2
        echo "Only .R and .Rmd files are supported." >&2
        exit 1
        ;;

esac
