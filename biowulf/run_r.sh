#!/bin/bash
set -e

module purge
conda deactivate 2>/dev/null || true
module load R

PREFIX="/vf/users/$USER/r-syslibs"
PROJECT_DIR="/spin1/home/linux/pereiralobof2/Projects/Erin/M_lignano/Mlig_scRNASeq_atlas"

export PATH="${PREFIX}/bin:${PATH}"
export PKG_CONFIG_PATH="${PREFIX}/lib/pkgconfig:${PKG_CONFIG_PATH}"
export LD_LIBRARY_PATH="${PREFIX}/lib:${LD_LIBRARY_PATH}"

cd "${PROJECT_DIR}"

echo "Using R: $(which R)"
echo "PREFIX=${PREFIX}"
echo "PKG_CONFIG_PATH=${PKG_CONFIG_PATH}"
echo "LD_LIBRARY_PATH=${LD_LIBRARY_PATH}"

if [ -z "$1" ]; then
    R
else
    Rscript "$1"
fi
