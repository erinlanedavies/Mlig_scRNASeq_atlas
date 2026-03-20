# Environment Setup and Reproducibility Guide

## Purpose

This document describes how to fully reproduce the computational environment required for the **Mlig_scRNASeq_atlas** project.

The setup combines:

- `renv` for R package version locking
- External storage for R libraries (due to limited home space)
- A dedicated system library prefix (via conda) for native dependencies required by compiled R packages

This design ensures:

- Reproducibility across sessions and machines
- Isolation from system-level inconsistencies on Biowulf
- Stability for packages with compiled dependencies (e.g., `igraph`, `ragg`, `leidenbase`)

---

## Project Assumptions

### Paths

Project root:
/spin1/home/linux/pereiralobof2/Projects/Erin/M_lignano/Mlig_scRNASeq_atlas

renv library root:
/vf/users/pereiralobof2/renv_libs

renv cache:
/vf/users/pereiralobof2/renv_cache

System libraries prefix:
/vf/users/pereiralobof2/r-syslibs

### Environment

- Cluster: NIH Biowulf
- R: module-loaded (`module load R`)
- OS: Rocky Linux 8.x

---

## Why This Setup Exists

Several issues required a controlled environment:

- `rlang` version mismatch when loading `Seurat`
- GitHub dependency (`DoubletFinder`) required authentication
- Bioconductor packages (`scater`, `speckle`) not installed via CRAN
- `igraph` and `leidenbase` failed due to missing ICU 75
- `ragg` failed due to missing image libraries

Final solution:

- Use `renv` for R
- Use a dedicated conda prefix for native libraries
- Launch all R sessions through a controlled script

---

## One-Time Machine Setup

```bash
conda create -y -p /vf/users/$USER/r-syslibs \
  -c conda-forge \
  icu=75 \
  libjpeg-turbo \
  freetype \
  libpng \
  libtiff \
  libwebp \
  pkg-config
```

---

## Launcher Script

All R sessions must be launched with:

```bash
./biowulf/run_r.sh
```

---

## renv Workflow

First-time:

```r
renv::restore()
```

Add package:

```r
renv::install("leidenbase")
renv::snapshot(packages = "leidenbase", prompt = FALSE)
```

---

## Troubleshooting

ICU error:

```
libicui18n.so.75: cannot open shared object file
```

Fix: ensure LD_LIBRARY_PATH includes syslibs.

---

## Summary

Always run:

```bash
./biowulf/run_r.sh
```
