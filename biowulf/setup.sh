#!/usr/bin/env bash

# ==============================================================================
# setup.sh
#
# One-command environment setup for Mlig_scRNASeq_atlas on NIH Biowulf.
#
# This script:
#   1. determines the project root
#   2. creates the Python environment
#   3. loads the required Biowulf modules
#   4. configures persistent renv locations
#   5. restores the R environment from renv.lock
#   6. validates R and Python
#
# Usage:
#
#   ./biowulf/setup.sh
#
# ==============================================================================

set -euo pipefail


# ------------------------------------------------------------------------------
# Project
# ------------------------------------------------------------------------------

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$PROJECT_ROOT"


# ------------------------------------------------------------------------------
# Environment locations
# ------------------------------------------------------------------------------

ENV_PREFIX="${ENV_PREFIX:-/data/${USER}/conda_envs/mlig_scrna_runtime}"
RENV_PATHS_LIBRARY_ROOT="${RENV_PATHS_LIBRARY_ROOT:-/vf/users/${USER}/renv_libs}"
RENV_PATHS_CACHE="${RENV_PATHS_CACHE:-/vf/users/${USER}/renv_cache}"
XDG_CACHE_HOME="${XDG_CACHE_HOME:-/data/${USER}/.cache}"

export ENV_PREFIX
export RENV_PATHS_LIBRARY_ROOT
export RENV_PATHS_CACHE
export XDG_CACHE_HOME

# ------------------------------------------------------------------------------
# Installation log
# ------------------------------------------------------------------------------

LOG_DIR="${PROJECT_ROOT}/logs"
mkdir -p "$LOG_DIR"

LOG_FILE="${LOG_DIR}/setup_$(date '+%Y%m%d_%H%M%S').log"

# Send all subsequent stdout and stderr both to the terminal and to the log.
exec > >(tee -a "$LOG_FILE") 2>&1

echo "Installation log:"
echo "  $LOG_FILE"
echo

# ------------------------------------------------------------------------------
# Introduction
# ------------------------------------------------------------------------------

echo
echo "============================================================"
echo "Mlig_scRNASeq_atlas environment setup"
echo "============================================================"
echo
echo "Project:"
echo "  $PROJECT_ROOT"
echo
echo "Python environment:"
echo "  $ENV_PREFIX"
echo
echo "renv library root:"
echo "  $RENV_PATHS_LIBRARY_ROOT"
echo
echo "renv cache:"
echo "  $RENV_PATHS_CACHE"
echo


# ------------------------------------------------------------------------------
# reticulate Python environment
# ------------------------------------------------------------------------------

echo "=== Setting up Python environment ==="
echo

if [[ -x "${ENV_PREFIX}/bin/python" ]]; then

    echo "Python environment already exists:"
    echo "  $ENV_PREFIX"
    echo
    echo "Skipping Conda environment creation."

else

    ./environments/reticulate/setup.sh

fi


# ------------------------------------------------------------------------------
# Biowulf modules
# ------------------------------------------------------------------------------

echo
echo "=== Loading Biowulf modules ==="
echo

module purge

module load gcc/11.3.0
module load hdf5/1.12.2
module load netcdf/4.9.0_gcc-11.3.0
module load openmpi/5.0.5/gcc-11.3.0

module load pcre2/10.40_gcc-11.3.0
module load libtiff/4.6.0_gcc-11.3.0
module load libwebp/1.6.0-gcc-11.3.0

module load R/4.5.2

module load pandoc/2.18
module load tex/2024


# ------------------------------------------------------------------------------
# Validate native environment
# ------------------------------------------------------------------------------

if ! command -v Rscript >/dev/null 2>&1; then
    echo "ERROR: Rscript not found." >&2
    exit 1
fi

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
# Cache directories
# ------------------------------------------------------------------------------

mkdir -p "$RENV_PATHS_LIBRARY_ROOT"
mkdir -p "$RENV_PATHS_CACHE"
mkdir -p "$XDG_CACHE_HOME"


# ------------------------------------------------------------------------------
# Python / reticulate
# ------------------------------------------------------------------------------

PYTHON="${ENV_PREFIX}/bin/python"

if [[ ! -x "$PYTHON" ]]; then
    echo "ERROR: Python executable not found after setup:" >&2
    echo "  $PYTHON" >&2
    exit 1
fi

export RETICULATE_PYTHON="$PYTHON"


# ------------------------------------------------------------------------------
# Restore R environment
# ------------------------------------------------------------------------------

echo
echo "=== Restoring R environment ==="
echo

Rscript --vanilla biowulf/setup_renv.R


# ------------------------------------------------------------------------------
# Validate Python
# ------------------------------------------------------------------------------

echo
echo "=== Validating Python environment ==="
echo

"$PYTHON" - <<'PY'
import sys
import anndata
import h5py
import numpy
import pandas
import scipy

print("Python:", sys.version.split()[0])
print("anndata:", anndata.__version__)
print("h5py:", h5py.__version__)
print("numpy:", numpy.__version__)
print("pandas:", pandas.__version__)
print("scipy:", scipy.__version__)
PY


# ------------------------------------------------------------------------------
# Validate R + reticulate
# ------------------------------------------------------------------------------

echo
echo "=== Validating R environment ==="
echo

Rscript --vanilla - <<'RS'
renv::load()

cat("R:", R.version.string, "\n")
cat("renv:", as.character(packageVersion("renv")), "\n")

packages <- c(
    "Seurat",
    "EnhancedVolcano",
    "DoubletFinder",
    "harmony",
    "leidenbase",
    "scater",
    "sctransform",
    "igraph",
    "ragg",
    "haven"
)

cat("\nKey R packages:\n")

for (pkg in packages) {

    if (!requireNamespace(pkg, quietly = TRUE)) {
        stop("Package unavailable: ", pkg)
    }

    cat(
        sprintf(
            "  %-20s %s\n",
            pkg,
            as.character(packageVersion(pkg))
        )
    )
}


# ------------------------------------------------------------------------------
# reticulate
# ------------------------------------------------------------------------------

cat("\nreticulate:\n")

if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Package unavailable: reticulate")
}

config <- reticulate::py_config()

cat("  Python:", config$python, "\n")
cat("  Version:", as.character(config$version), "\n")


# ------------------------------------------------------------------------------
# Required Python imports through reticulate
# ------------------------------------------------------------------------------

for (module in c("anndata", "h5py", "numpy", "pandas", "scipy")) {

    if (!reticulate::py_module_available(module)) {
        stop("Python module unavailable through reticulate: ", module)
    }
}

cat("\nPython modules are accessible through reticulate.\n")
RS


# ------------------------------------------------------------------------------
# Finished
# ------------------------------------------------------------------------------

echo
echo "============================================================"
echo "Environment setup completed successfully."
echo "============================================================"
echo
echo "Run the project with:"
echo
echo "  ./biowulf/run_r.sh scripts/KnitReports.R"
echo
