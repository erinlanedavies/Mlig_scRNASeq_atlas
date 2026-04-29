#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$PROJECT_ROOT"

# -----------------------------
# Biowulf modules
# -----------------------------
module purge

module load gcc/11.3.0
module load hdf5/1.12.2
module load netcdf/4.9.0
module load openmpi/5.0.5
module load pandoc/2.18
module load tex/2024
module load pcre2/10.40
module load R/4.5.2

# -----------------------------
# Locate conda / mamba / micromamba
# -----------------------------
CONDA_BASE=""

if command -v conda >/dev/null 2>&1; then
  CONDA_BASE="$(conda info --base)"
elif [ -d "$HOME/miniforge3" ]; then
  CONDA_BASE="$HOME/miniforge3"
elif [ -d "$HOME/miniconda3" ]; then
  CONDA_BASE="$HOME/miniconda3"
elif [ -d "/data/${USER}/miniforge3" ]; then
  CONDA_BASE="/data/${USER}/miniforge3"
elif [ -d "/data/${USER}/miniconda3" ]; then
  CONDA_BASE="/data/${USER}/miniconda3"
else
  echo "ERROR: conda was not found." >&2
  echo "Load conda first, or edit CONDA_BASE in biowulf/run_r.sh." >&2
  exit 1
fi

source "${CONDA_BASE}/etc/profile.d/conda.sh"

# -----------------------------
# Conda runtime env
# -----------------------------
ENV_PREFIX="/data/${USER}/conda_envs/mlig_scrna_runtime"

if [ ! -x "${ENV_PREFIX}/bin/python" ]; then
  conda create -y -p "${ENV_PREFIX}" \
    -c conda-forge \
    python=3.10
fi

conda activate "${ENV_PREFIX}"

# Ensure libraries needed by R packages compiled from source are present
conda install -y -p "${ENV_PREFIX}" \
  -c conda-forge \
  anndata \
  h5py \
  numpy \
  pandas \
  scipy \
  libiconv \
  icu \
  libwebp \
  pkg-config \
  freetype \
  libpng \
  jpeg \
  libjpeg-turbo \
  libtiff \
  zlib \
  bzip2 \
  harfbuzz \
  fribidi \
  fontconfig

export CXX11="g++"
export CXX11STD="-std=gnu++14"
export CXX11FLAGS="-O2 -march=haswell -mtune=generic"
export PKG_CXXFLAGS="-std=gnu++14"

export PATH="${CONDA_PREFIX}/bin:${PATH}"
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH:-}"
export PKG_CONFIG_PATH="${CONDA_PREFIX}/lib/pkgconfig:${CONDA_PREFIX}/share/pkgconfig:${PKG_CONFIG_PATH:-}"
export CPATH="${CONDA_PREFIX}/include:${CONDA_PREFIX}/include/freetype2:${CPATH:-}"
export C_INCLUDE_PATH="${CONDA_PREFIX}/include:${CONDA_PREFIX}/include/freetype2:${C_INCLUDE_PATH:-}"
export CPLUS_INCLUDE_PATH="${CONDA_PREFIX}/include:${CONDA_PREFIX}/include/freetype2:${CPLUS_INCLUDE_PATH:-}"
export LIBRARY_PATH="${CONDA_PREFIX}/lib:${LIBRARY_PATH:-}"

export INCLUDE_DIR="${CONDA_PREFIX}/include:${CONDA_PREFIX}/include/freetype2"
export LIB_DIR="${CONDA_PREFIX}/lib"

export RETICULATE_PYTHON="${CONDA_PREFIX}/bin/python"
export RETICULATE_PYENV_ROOT="/data/${USER}/reticulate_pyenv"
export BASILISK_DIR="/data/${USER}/basilisk_envs"

# -----------------------------
# renv external paths
# -----------------------------
export RENV_PATHS_LIBRARY_ROOT="/data/${USER}/renv_libs"
export RENV_PATHS_CACHE="/data/${USER}/renv_cache"

mkdir -p "$RENV_PATHS_LIBRARY_ROOT" "$RENV_PATHS_CACHE"
mkdir -p "$RETICULATE_PYENV_ROOT" "$BASILISK_DIR"

# -----------------------------
# Dispatch
# -----------------------------
TARGET="${1:-setup_renv.R}"

case "$TARGET" in
  *.R)
    Rscript --vanilla "$TARGET"
    ;;
  *.Rmd)
    Rscript --vanilla biowulf/render_rmd.R "$TARGET"

    ;;
  *)
    echo "Unsupported target: $TARGET" >&2
    echo "Use an .R or .Rmd file." >&2
    exit 1
    ;;
esac
