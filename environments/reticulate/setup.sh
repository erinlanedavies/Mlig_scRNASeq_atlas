#!/usr/bin/env bash

set -euo pipefail

# ==============================================================================
# Install the Python environment used by R through reticulate
#
# This script creates the Conda environment required by the R/Seurat pipeline
# for Python interoperability through the reticulate package.
#
# It does NOT install R packages. R dependencies are managed separately by renv.
#
# Usage:
#
#   ./environments/reticulate/setup.sh
#
# To install into a different location:
#
#   ENV_PREFIX=/path/to/environment \
#     ./environments/reticulate/setup.sh
#
# ==============================================================================


# ------------------------------------------------------------------------------
# Resolve repository paths
# ------------------------------------------------------------------------------

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

ENV_PREFIX="${ENV_PREFIX:-/data/${USER}/conda_envs/mlig_scrna_runtime}"

ENV_FILE="${SCRIPT_DIR}/environment.yml"
LOCK_FILE="${SCRIPT_DIR}/conda-linux-64.lock.txt"

# ------------------------------------------------------------------------------
# Basic information
# ------------------------------------------------------------------------------

echo "============================================================"
echo "M. lignano scRNA-seq atlas"
echo "Reticulate Python environment setup"
echo "============================================================"
echo
echo "Environment files: ${SCRIPT_DIR}"
echo "Environment prefix: ${ENV_PREFIX}"
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
# Validate environment files
# ------------------------------------------------------------------------------

if [[ ! -f "${ENV_FILE}" ]]; then
    echo "ERROR: Environment specification not found:" >&2
    echo "  ${ENV_FILE}" >&2
    exit 1
fi

if [[ ! -f "${LOCK_FILE}" ]]; then
    echo "ERROR: Conda lock file not found:" >&2
    echo "  ${LOCK_FILE}" >&2
    exit 1
fi


# ------------------------------------------------------------------------------
# Refuse to overwrite an existing environment
# ------------------------------------------------------------------------------

if [[ -d "${ENV_PREFIX}" ]]; then
    echo "ERROR: Environment already exists:" >&2
    echo "  ${ENV_PREFIX}" >&2
    echo
    echo "Remove it explicitly before reinstalling." >&2
    exit 1
fi


# ------------------------------------------------------------------------------
# Create environment
#
# The explicit Linux-64 lock file contains the exact Conda package artifacts
# from the validated environment and is therefore preferred for reproducible
# installation on Biowulf.
# ------------------------------------------------------------------------------

echo "Creating Conda environment from validated lock file..."
echo

conda create \
    --prefix "${ENV_PREFIX}" \
    --file "${LOCK_FILE}" \
    --yes

echo
echo "Conda environment created."
echo


# ------------------------------------------------------------------------------
# Validate Python
# ------------------------------------------------------------------------------

PYTHON="${ENV_PREFIX}/bin/python"

if [[ ! -x "${PYTHON}" ]]; then
    echo "ERROR: Python executable was not created:" >&2
    echo "  ${PYTHON}" >&2
    exit 1
fi

echo "Python:"
"${PYTHON}" --version
echo


# ------------------------------------------------------------------------------
# Validate required Python modules
# ------------------------------------------------------------------------------

echo "Validating Python environment..."
echo

"${PYTHON}" - <<'PY'
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


# ------------------------------------------------------------------------------
# Check package consistency
# ------------------------------------------------------------------------------

echo
echo "Checking Python package dependencies..."
echo

"${PYTHON}" -m pip check


# ------------------------------------------------------------------------------
# Final status
# ------------------------------------------------------------------------------

echo
echo "============================================================"
echo "Reticulate Python environment installed successfully."
echo "============================================================"
echo
echo "Environment:"
echo "  ${ENV_PREFIX}"
echo
echo "Python:"
echo "  ${PYTHON}"
echo
echo "The R pipeline should access this environment through:"
echo
echo "  RETICULATE_PYTHON=${PYTHON}"
echo
