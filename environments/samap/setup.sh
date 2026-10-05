#!/usr/bin/env bash

set -euo pipefail

# ==============================================================================
# setup.sh
#
# Install the Python environment used for SAMap analyses.
#
# This environment is independent of:
#   - the Python environment used by R through reticulate
#   - the R environment managed by renv
#
# Usage:
#
#   ./environments/samap/setup.sh
#
# To create an environment with a different name:
#
#   SAMAP_ENV_NAME=mlig_samap_test \
#     ./environments/samap/setup.sh
#
# ==============================================================================


# ------------------------------------------------------------------------------
# Resolve repository paths
# ------------------------------------------------------------------------------

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

SAMAP_ENV_NAME="${SAMAP_ENV_NAME:-mlig_samap}"

ENV_FILE="${SCRIPT_DIR}/environment.yml"
LOCK_FILE="${SCRIPT_DIR}/environment.lock.yml"


# ------------------------------------------------------------------------------
# Basic information
# ------------------------------------------------------------------------------

echo "============================================================"
echo "M. lignano SAMap environment setup"
echo "============================================================"
echo
echo "Project root:"
echo "  ${PROJECT_ROOT}"
echo
echo "Environment:"
echo "  ${SAMAP_ENV_NAME}"
echo
echo "Environment specification:"
echo "  ${ENV_FILE}"
echo
echo "Environment lock:"
echo "  ${LOCK_FILE}"
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

if conda env list | awk '{print $1}' | grep -Fxq "${SAMAP_ENV_NAME}"; then
    echo "ERROR: Conda environment '${SAMAP_ENV_NAME}' already exists." >&2
    echo
    echo "Remove it explicitly before reinstalling:"
    echo
    echo "  conda env remove -n ${SAMAP_ENV_NAME}"
    echo
    exit 1
fi


# ------------------------------------------------------------------------------
# Create environment
# ------------------------------------------------------------------------------

echo "Creating SAMap environment..."
echo

conda env create \
    --name "${SAMAP_ENV_NAME}" \
    --file "${LOCK_FILE}"

echo
echo "SAMap environment created."
echo


# ------------------------------------------------------------------------------
# Locate environment
# ------------------------------------------------------------------------------

SAMAP_ENV_PREFIX="$(
    conda env list |
        awk -v env="${SAMAP_ENV_NAME}" '$1 == env {print $NF}'
)"

if [[ -z "${SAMAP_ENV_PREFIX}" ]]; then
    echo "ERROR: Could not determine SAMap environment prefix." >&2
    exit 1
fi

PYTHON="${SAMAP_ENV_PREFIX}/bin/python"

if [[ ! -x "${PYTHON}" ]]; then
    echo "ERROR: Python executable was not found:" >&2
    echo "  ${PYTHON}" >&2
    exit 1
fi

echo "Environment prefix:"
echo "  ${SAMAP_ENV_PREFIX}"
echo
echo "Python:"
"${PYTHON}" --version
echo


# ------------------------------------------------------------------------------
# Validate package consistency
# ------------------------------------------------------------------------------

echo "Checking Python package dependencies..."
echo

"${PYTHON}" -m pip check


# ------------------------------------------------------------------------------
# Validate SAMap
# ------------------------------------------------------------------------------

echo
echo "Validating SAMap installation..."
echo

"${PYTHON}" - <<'PY'
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
# Final status
# ------------------------------------------------------------------------------

echo
echo "============================================================"
echo "SAMap environment setup completed successfully."
echo "============================================================"
echo
echo "Environment:"
echo "  ${SAMAP_ENV_NAME}"
echo
echo "Environment prefix:"
echo "  ${SAMAP_ENV_PREFIX}"
echo
echo "Python:"
echo "  ${PYTHON}"
echo
echo "Activate it with:"
echo
echo "  conda activate ${SAMAP_ENV_NAME}"
echo
