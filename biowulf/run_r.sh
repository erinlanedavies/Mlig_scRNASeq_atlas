#!/usr/bin/env bash
set -euo pipefail

# ---- Start an interactive sesstion ---------------------------
start_sinteractive -m 150g -c 10 -t 36:00:00

# ---- Conda layer (tools / python) ----------------------------
conda activate Mlig_scRNASeq_atlas

# ---- R layer -------------------------------------------------
module load rstudio-server

# ---- R packages (renv) ---------------------------------------
R -q -e "renv::restore(prompt = FALSE)"

