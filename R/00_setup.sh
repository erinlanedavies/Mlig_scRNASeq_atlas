#!/usr/bin/env bash
set -euo pipefail

# 1) conda
mamba create -n Mlig_scRNASeq_atlas --file conda-spec.txt -y || true
mamba activate Mlig_scRNASeq_atlas

# 2) renv
R -q -e "renv::restore(prompt = FALSE)"
