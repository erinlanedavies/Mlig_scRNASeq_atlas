#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

ENV_NAME="${ENV_NAME:-mlig_samap}"

echo "============================================================"
echo "M. lignano SAMap environment setup"
echo "============================================================"
echo
echo "Project root: ${PROJECT_ROOT}"
echo "Environment:  ${ENV_NAME}"
echo

# ------------------------------------------------------------
# Check Conda
# ------------------------------------------------------------

if ! command -v conda >/dev/null 2>&1; then
    echo "ERROR: conda was not found in PATH." >&2
    exit 1
fi

echo "Conda:"
conda --version
echo

# ------------------------------------------------------------
# Refuse to overwrite an existing environment
# ------------------------------------------------------------

if conda env list | awk '{print $1}' | grep -Fxq "${ENV_NAME}"; then
    echo "ERROR: Conda environment '${ENV_NAME}' already exists." >&2
    echo
    echo "Remove it first with:"
    echo
    echo "  conda env remove -n ${ENV_NAME}"
    echo
    exit 1
fi

# ------------------------------------------------------------
# Create environment from the validated lock
# ------------------------------------------------------------

echo "Creating SAMap environment..."
echo

conda env create \
    --name "${ENV_NAME}" \
    --file "${SCRIPT_DIR}/environment.lock.yml"

echo
echo "Environment created."
echo

# ------------------------------------------------------------
# Locate its Python executable
# ------------------------------------------------------------

ENV_PREFIX="$(conda env list | awk -v env="${ENV_NAME}" '$1 == env {print $NF}')"

if [[ -z "${ENV_PREFIX}" ]]; then
    echo "ERROR: Could not determine environment prefix." >&2
    exit 1
fi

PYTHON="${ENV_PREFIX}/bin/python"

if [[ ! -x "${PYTHON}" ]]; then
    echo "ERROR: Python not found at ${PYTHON}" >&2
    exit 1
fi

# ------------------------------------------------------------
# Validate package consistency
# ------------------------------------------------------------

echo "Checking Python package dependencies..."
"${PYTHON}" -m pip check

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

print("Python:", sys.version)
print("Executable:", sys.executable)
print()

for package in packages:
    print(f"{package:20s} {version(package)}")

print()

from samap import SAMAP
from samap.analysis import (
    get_mapping_scores,
    GenePairFinder,
    CellTypeTriangles,
    sankey_plot,
)
from samap_extension import plotting

print("SAMAP:", SAMAP)
print("SAMAP module:", SAMAP.__module__)
print("samap-extension:", plotting.__file__)
print()
print("SAMap environment validation successful.")
PY

echo
echo "============================================================"
echo "SAMap environment setup completed successfully."
echo "============================================================"
echo
echo "Activate it with:"
echo
echo "  conda activate ${ENV_NAME}"
echo
