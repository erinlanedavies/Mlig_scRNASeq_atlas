#!/usr/bin/env bash
set -euo pipefail

# Create the Python environment used by the M. lignano scRNA-seq atlas
# on NIH Biowulf.
#
# The environment is created from the explicit linux-64 Conda lock file,
# which records the exact package builds used for the validated runtime.

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

ENV_PREFIX="${ENV_PREFIX:-/data/${USER}/conda_envs/mlig_scrna_runtime}"
LOCK_FILE="${PROJECT_ROOT}/biowulf/conda-linux-64.lock.txt"

if [ ! -f "${LOCK_FILE}" ]; then
    echo "ERROR: Conda lock file not found:" >&2
    echo "  ${LOCK_FILE}" >&2
    exit 1
fi

# Locate Conda.
if command -v conda >/dev/null 2>&1; then
    CONDA_BASE="$(conda info --base)"
elif [ -d "${HOME}/miniforge3" ]; then
    CONDA_BASE="${HOME}/miniforge3"
elif [ -d "${HOME}/miniconda3" ]; then
    CONDA_BASE="${HOME}/miniconda3"
elif [ -d "/data/${USER}/miniforge3" ]; then
    CONDA_BASE="/data/${USER}/miniforge3"
elif [ -d "/data/${USER}/miniconda3" ]; then
    CONDA_BASE="/data/${USER}/miniconda3"
else
    echo "ERROR: Conda installation not found." >&2
    exit 1
fi

source "${CONDA_BASE}/etc/profile.d/conda.sh"

# Do not overwrite an existing environment.
if [ -e "${ENV_PREFIX}" ]; then
    echo "ERROR: Conda environment already exists:" >&2
    echo "  ${ENV_PREFIX}" >&2
    echo >&2
    echo "Remove or rename it before creating a new environment." >&2
    exit 1
fi

echo "Creating Conda environment:"
echo "  ${ENV_PREFIX}"
echo
echo "Using explicit lock:"
echo "  ${LOCK_FILE}"
echo

conda create \
    --yes \
    --prefix "${ENV_PREFIX}" \
    --file "${LOCK_FILE}"

echo
echo "Validating Python environment..."

"${ENV_PREFIX}/bin/python" - <<'PY'
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

echo
echo "Conda environment created successfully."
echo "Path: ${ENV_PREFIX}"
