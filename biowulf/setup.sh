#!/usr/bin/env bash

# ==============================================================================
# setup.sh
#
# One-command environment setup for Mlig_scRNASeq_atlas on NIH Biowulf.
#
# This script:
#   1. determines the project root
#   2. creates the Python environment used by R through reticulate
#   3. creates the Python environment used by SAMap
#   4. removes inherited Conda state before native R package compilation
#   5. loads the required Biowulf modules
#   6. configures persistent renv locations
#   7. restores the R environment from renv.lock
#   8. validates R, reticulate, and SAMap
#
# Normal usage:
#
#   ./biowulf/setup.sh
#
# Environment locations/names can be overridden for clean-room testing:
#
#   ENV_PREFIX=/data/${USER}/conda_envs/mlig_scrna_runtime_test \
#   SAMAP_ENV_NAME=mlig_samap_test \
#   RENV_PATHS_LIBRARY_ROOT=/vf/users/${USER}/renv_libs_test \
#   RENV_PATHS_CACHE=/vf/users/${USER}/renv_cache_test \
#   XDG_CACHE_HOME=/data/${USER}/.cache_mlig_test \
#     ./biowulf/setup.sh
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
SAMAP_ENV_NAME="${SAMAP_ENV_NAME:-mlig_samap}"

RENV_PATHS_LIBRARY_ROOT="${RENV_PATHS_LIBRARY_ROOT:-/vf/users/${USER}/renv_libs}"
RENV_PATHS_CACHE="${RENV_PATHS_CACHE:-/vf/users/${USER}/renv_cache}"
XDG_CACHE_HOME="${XDG_CACHE_HOME:-/data/${USER}/.cache}"

export ENV_PREFIX
export SAMAP_ENV_NAME
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
echo "reticulate Python environment:"
echo "  $ENV_PREFIX"
echo
echo "SAMap Python environment:"
echo "  $SAMAP_ENV_NAME"
echo
echo "renv library root:"
echo "  $RENV_PATHS_LIBRARY_ROOT"
echo
echo "renv cache:"
echo "  $RENV_PATHS_CACHE"
echo


# ------------------------------------------------------------------------------
# Validate Conda
# ------------------------------------------------------------------------------

if ! command -v conda >/dev/null 2>&1; then
    echo "ERROR: conda was not found in PATH." >&2
    echo
    echo "Initialize Conda before running this script." >&2
    exit 1
fi

echo "Conda:"
conda --version
echo


# ------------------------------------------------------------------------------
# reticulate Python environment
# ------------------------------------------------------------------------------

echo "=== Setting up reticulate Python environment ==="
echo

if [[ -x "${ENV_PREFIX}/bin/python" ]]; then

    echo "reticulate Python environment already exists:"
    echo "  $ENV_PREFIX"
    echo
    echo "Skipping reticulate environment creation."

else

    ./environments/reticulate/setup.sh

fi


# ------------------------------------------------------------------------------
# SAMap Python environment
# ------------------------------------------------------------------------------

echo
echo "=== Setting up SAMap Python environment ==="
echo

if conda env list | awk '{print $1}' | grep -Fxq "${SAMAP_ENV_NAME}"; then

    echo "SAMap environment already exists:"
    echo "  ${SAMAP_ENV_NAME}"
    echo
    echo "Skipping SAMap environment creation."

else

    ./environments/samap/setup.sh

fi


# ------------------------------------------------------------------------------
# Locate SAMap environment
# ------------------------------------------------------------------------------

SAMAP_ENV_PREFIX="$(
    conda env list |
        awk -v env="${SAMAP_ENV_NAME}" '$1 == env {print $NF}'
)"

if [[ -z "${SAMAP_ENV_PREFIX}" ]]; then
    echo "ERROR: Could not determine SAMap environment prefix." >&2
    exit 1
fi

SAMAP_PYTHON="${SAMAP_ENV_PREFIX}/bin/python"

if [[ ! -x "${SAMAP_PYTHON}" ]]; then
    echo "ERROR: SAMap Python executable not found:" >&2
    echo "  ${SAMAP_PYTHON}" >&2
    exit 1
fi


# ------------------------------------------------------------------------------
# Remove inherited Conda environment
#
# R packages must be compiled against the Biowulf/system native libraries,
# not libraries from an interactive Conda environment.
#
# In particular, an auto-activated Conda base environment can place its
# xml2-config on PATH. This can cause R packages such as igraph to link
# against Conda's libxml2 / ICU stack rather than the system libraries.
#
# Both Python environments have already been created above, so Conda activation
# is no longer required for the R build.
# ------------------------------------------------------------------------------

if [[ -n "${CONDA_PREFIX:-}" ]]; then

    ACTIVE_CONDA_PREFIX="${CONDA_PREFIX}"

    echo
    echo "=== Removing inherited Conda environment from R build ==="
    echo
    echo "Detected active Conda environment:"
    echo "  ${ACTIVE_CONDA_PREFIX}"
    echo

    PATH="$(
        printf '%s\n' "$PATH" |
            tr ':' '\n' |
            grep -v "^${ACTIVE_CONDA_PREFIX}/bin$" |
            paste -sd:
    )"

    export PATH

    unset CONDA_PREFIX
    unset CONDA_DEFAULT_ENV
    unset CONDA_PROMPT_MODIFIER
    unset CONDA_SHLVL

    echo "Inherited Conda environment removed from R build PATH."
    echo

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

echo
echo "=== Validating native build environment ==="
echo

if ! command -v Rscript >/dev/null 2>&1; then
    echo "ERROR: Rscript not found." >&2
    exit 1
fi


# ------------------------------------------------------------------------------
# pkg-config
# ------------------------------------------------------------------------------

PKG_CONFIG_BIN="$(command -v pkg-config || true)"

if [[ "$PKG_CONFIG_BIN" != "/usr/bin/pkg-config" ]]; then
    echo "ERROR: Expected /usr/bin/pkg-config but found:" >&2
    echo "  ${PKG_CONFIG_BIN:-not found}" >&2
    exit 1
fi

echo "pkg-config:"
echo "  ${PKG_CONFIG_BIN}"


# ------------------------------------------------------------------------------
# xml2-config
#
# This check prevents Conda's libxml2 / ICU libraries from contaminating
# compilation of R packages such as igraph.
# ------------------------------------------------------------------------------

XML2_CONFIG_BIN="$(command -v xml2-config || true)"

if [[ "$XML2_CONFIG_BIN" != "/usr/bin/xml2-config" ]]; then
    echo "ERROR: Expected /usr/bin/xml2-config but found:" >&2
    echo "  ${XML2_CONFIG_BIN:-not found}" >&2
    echo >&2
    echo "A non-system xml2-config can cause R packages to link against" >&2
    echo "incompatible libxml2 / ICU libraries." >&2
    exit 1
fi

echo "xml2-config:"
echo "  ${XML2_CONFIG_BIN}"


# ------------------------------------------------------------------------------
# libxml2
# ------------------------------------------------------------------------------

LIBXML2_PREFIX="$(pkg-config --variable=prefix libxml-2.0)"

if [[ "$LIBXML2_PREFIX" != "/usr" ]]; then
    echo "ERROR: Expected system libxml2 prefix /usr but found:" >&2
    echo "  ${LIBXML2_PREFIX}" >&2
    exit 1
fi

echo "libxml2 prefix:"
echo "  ${LIBXML2_PREFIX}"


# ------------------------------------------------------------------------------
# libwebp
# ------------------------------------------------------------------------------

if ! pkg-config --exists libwebp; then
    echo "ERROR: pkg-config cannot locate libwebp." >&2
    exit 1
fi

if ! pkg-config --exists libwebpmux; then
    echo "ERROR: pkg-config cannot locate libwebpmux." >&2
    exit 1
fi

echo "libwebp:"
echo "  $(pkg-config --modversion libwebp)"

echo "libwebpmux:"
echo "  $(pkg-config --modversion libwebpmux)"

echo
echo "Native build environment validated successfully."


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
# Validate reticulate Python environment
# ------------------------------------------------------------------------------

echo
echo "=== Validating reticulate Python environment ==="
echo

"$PYTHON" - <<'PY'
import sys
import anndata
import h5py
import numpy
import pandas
import scipy

print("Python:", sys.version.split()[0])
print("Executable:", sys.executable)
print()
print("Python packages:")
print("  anndata:", anndata.__version__)
print("  h5py:   ", h5py.__version__)
print("  numpy:  ", numpy.__version__)
print("  pandas: ", pandas.__version__)
print("  scipy:  ", scipy.__version__)
PY

echo
echo "Checking reticulate Python package dependencies..."
echo

"$PYTHON" -m pip check


# ------------------------------------------------------------------------------
# Validate SAMap environment
#
# This is performed even when the environment already existed and creation was
# skipped above.
# ------------------------------------------------------------------------------

echo
echo "=== Validating SAMap environment ==="
echo

"${SAMAP_PYTHON}" -m pip check

echo

"${SAMAP_PYTHON}" - <<'PY'
import sys
from importlib.metadata import version

packages = [
    "sc-samap",
    "samap-extension",
    "scanpy",
    "anndata",
    "numpy",
    "pandas",
    "scipy",
]

print("Python:", sys.version.split()[0])
print("Executable:", sys.executable)
print()

print("Python packages:")

for package in packages:
    print(f"  {package:20s} {version(package)}")

print()

from samap import SAMAP
from samap.analysis import (
    get_mapping_scores,
    GenePairFinder,
    CellTypeTriangles,
    sankey_plot,
)
from samap_extension import plotting

print("SAMAP class:", SAMAP)
print("SAMAP module:", SAMAP.__module__)
print("samap-extension:", plotting.__file__)

print()
print("SAMap environment validation successful.")
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
echo "Installed environments:"
echo
echo "  R:"
echo "    ${RENV_PATHS_LIBRARY_ROOT}"
echo
echo "  reticulate Python:"
echo "    ${ENV_PREFIX}"
echo
echo "  SAMap:"
echo "    ${SAMAP_ENV_NAME}"
echo "    ${SAMAP_ENV_PREFIX}"
echo
echo "Run the R pipeline with:"
echo
echo "  ./biowulf/run_r.sh scripts/KnitReports.R"
echo
echo "Activate the SAMap environment with:"
echo
echo "  conda activate ${SAMAP_ENV_NAME}"
echo
