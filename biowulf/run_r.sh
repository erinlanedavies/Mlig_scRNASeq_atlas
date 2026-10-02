#!/usr/bin/env bash

set -euo pipefail


# ============================================================
# Project configuration
# ============================================================

# Determine the repository root from the location of this script.
# This allows run_r.sh to be called from any working directory.
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

cd "$PROJECT_ROOT"


# Pre-created Conda environment containing only the Python runtime
# required by reticulate.
ENV_PREFIX="/data/${USER}/conda_envs/mlig_scrna_runtime"

# Locations for large R/Python runtime files.
RENV_ROOT="/vf/users/${USER}"
DATA_ROOT="/data/${USER}"


# ============================================================
# Biowulf modules
# ============================================================

# Start from a clean module environment.
module purge

# Compiler and native libraries used by R packages.
module load gcc/11.3.0
module load hdf5/1.12.2
module load netcdf/4.9.0
module load openmpi/5.0.5

# R Markdown dependencies.
module load pandoc/2.18
module load tex/2024

# R dependencies.
module load pcre2/10.40
module load R/4.5.2

# Required when building packages such as ragg from source.
module load libwebp/1.6.0-gcc-11.3.0


# ============================================================
# Python / Conda runtime
# ============================================================

# IMPORTANT:
#
# Do NOT activate the Conda environment here.
#
# Activating Conda places its bin directory at the beginning of
# PATH. That can cause R package configure scripts to use Conda's
# pkg-config, headers, and native libraries instead of the
# Biowulf/system toolchain.
#
# R therefore uses the Biowulf native build environment, while
# reticulate is pointed directly at the Python interpreter in the
# Conda environment.

if [ ! -x "${ENV_PREFIX}/bin/python" ]; then
    echo "ERROR: Python environment not found:" >&2
    echo "  ${ENV_PREFIX}" >&2
    echo >&2
    echo "Create it first with:" >&2
    echo "  biowulf/setup_conda.sh" >&2
    exit 1
fi

export RETICULATE_PYTHON="${ENV_PREFIX}/bin/python"


# Locations used by reticulate/basilisk if additional environments
# or runtime files are created.
export RETICULATE_PYENV_ROOT="${DATA_ROOT}/reticulate_pyenv"
export BASILISK_DIR="${DATA_ROOT}/basilisk_envs"

mkdir -p \
    "$RETICULATE_PYENV_ROOT" \
    "$BASILISK_DIR"


# ============================================================
# R compilation environment
# ============================================================

# Some packages in the project require at least C++14.
#
# In particular, this is required for compatibility between
# packages such as presto and the version of RcppArmadillo in the
# locked environment.

export CXX11="g++"
export CXX11STD="-std=gnu++14"
export CXX11FLAGS="-O2 -march=haswell -mtune=generic"
export PKG_CXXFLAGS="-std=gnu++14"


# IMPORTANT:
#
# Do NOT point any of the following variables at the Conda
# environment:
#
#   LD_LIBRARY_PATH
#   PKG_CONFIG_PATH
#   CPATH
#   C_INCLUDE_PATH
#   CPLUS_INCLUDE_PATH
#   LIBRARY_PATH
#   INCLUDE_DIR
#   LIB_DIR
#
# The Biowulf modules loaded above establish the native build
# environment. In particular, the libwebp module modifies
# PKG_CONFIG_PATH so that the system pkg-config can locate
# libwebp and libwebpmux.


# ============================================================
# Validate native build environment
# ============================================================

PKG_CONFIG="$(command -v pkg-config || true)"

if [ "$PKG_CONFIG" != "/usr/bin/pkg-config" ]; then
    echo "ERROR: unexpected pkg-config executable:" >&2
    echo "  ${PKG_CONFIG:-not found}" >&2
    echo >&2
    echo "Expected:" >&2
    echo "  /usr/bin/pkg-config" >&2
    echo >&2
    echo "A Conda environment may still be active." >&2
    echo "Run 'conda deactivate' and try again." >&2
    exit 1
fi


# Verify that the WebP libraries required by ragg can be found
# through the Biowulf/system pkg-config environment.
if ! pkg-config --exists libwebp libwebpmux; then
    echo "ERROR: pkg-config cannot find libwebp/libwebpmux." >&2
    echo "Check the Biowulf libwebp module." >&2
    exit 1
fi


# ============================================================
# renv paths
# ============================================================

# Keep the project-specific R library and renv package cache
# outside the project/home filesystem.
export RENV_PATHS_LIBRARY_ROOT="${RENV_ROOT}/renv_libs"
export RENV_PATHS_CACHE="${RENV_ROOT}/renv_cache"

mkdir -p \
    "$RENV_PATHS_LIBRARY_ROOT" \
    "$RENV_PATHS_CACHE"


# ============================================================
# General cache
# ============================================================

# Prevent R and related software from filling ~/.cache.
export XDG_CACHE_HOME="${DATA_ROOT}/.cache"

mkdir -p "$XDG_CACHE_HOME"


# ============================================================
# renv staging
# ============================================================

# renv uses project-local renv/staging while installing packages.
# The project filesystem has a relatively limited quota, so
# redirect staging to /data.

RENV_STAGING="${DATA_ROOT}/renv_staging/Mlig_scRNASeq_atlas"

mkdir -p "$RENV_STAGING"


# Refuse to replace a real renv/staging directory automatically.
# This protects against accidentally deleting package installation
# files or other data.
if [ -e "${PROJECT_ROOT}/renv/staging" ] &&
   [ ! -L "${PROJECT_ROOT}/renv/staging" ]; then

    echo "ERROR: ${PROJECT_ROOT}/renv/staging exists and is not a symlink." >&2
    echo >&2
    echo "If no renv installation is currently running, remove it manually:" >&2
    echo "  rm -rf ${PROJECT_ROOT}/renv/staging" >&2
    echo >&2
    echo "Then run this script again." >&2
    exit 1
fi


# Refresh the symlink so that it always points at the expected
# /data location.
if [ -L "${PROJECT_ROOT}/renv/staging" ]; then
    rm "${PROJECT_ROOT}/renv/staging"
fi

ln -s "$RENV_STAGING" "${PROJECT_ROOT}/renv/staging"


# ============================================================
# Target validation
# ============================================================

if [ "$#" -lt 1 ]; then
    echo "Usage:" >&2
    echo "  biowulf/run_r.sh <script.R>" >&2
    echo "  biowulf/run_r.sh <report.Rmd>" >&2
    exit 1
fi

TARGET="$1"

if [ ! -f "$TARGET" ]; then
    echo "ERROR: target not found:" >&2
    echo "  $TARGET" >&2
    exit 1
fi


# ============================================================
# Execute target
# ============================================================

case "$TARGET" in

    *.R)

        # --vanilla prevents user/site startup configuration from
        # silently changing the runtime.
        #
        # Because --vanilla also prevents the project's .Rprofile
        # from activating renv, activate the project explicitly
        # before sourcing the requested script.

        Rscript --vanilla -e \
            "renv::load(project = '$PROJECT_ROOT'); source('$TARGET', chdir = FALSE)"
        ;;


    *.Rmd)

        # render_rmd.R explicitly activates renv/activate.R before
        # calling rmarkdown::render().
        #
        # Pass the requested Rmd filename as its single trailing
        # command-line argument.

        Rscript --vanilla \
            biowulf/render_rmd.R \
            "$TARGET"
        ;;


    *)

        echo "ERROR: unsupported target:" >&2
        echo "  $TARGET" >&2
        echo >&2
        echo "Expected an .R or .Rmd file." >&2
        exit 1
        ;;

esac
