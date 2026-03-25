# M. lignano scRNA-seq Atlas Pipeline

A modular, reproducible, and publication-ready pipeline for single-cell RNA-seq analysis of *Macrostomum lignano* using R/Seurat.

Designed for large-scale atlas construction, integration, and downstream biological interpretation.

---

## Repository Structure

```
.
├── biowulf/                 # HPC execution script (NIH Biowulf)
├── config/                  # Configuration files (optional/extendable)
├── data/
│   └── metadata/            # Metadata tables (colors, markers, annotations)
├── docs/                    # Documentation (environment setup)
├── environment.yml          # Conda environment definition
├── R/                       # Core R functions (modularized)
│   └── report_utils/        # Thematic helper modules
├── Rmarkdown/               # Main analysis pipeline (Rmd reports)
├── scripts/                 # Utility scripts (batch knitting)
├── renv/                    # R environment management
└── README.md
```

---

## Pipeline Overview

The workflow is organized into sequential RMarkdown steps:

| Step | File | Description |
|------|------|------------|
| 00 | `00_CreateSeuratFromCounts.Rmd` | Build Seurat object from raw counts |
| 01 | `01_ComputeQCMetrics.Rmd` | Compute QC metrics |
| 02 | `02_FilterPlotQCData.Rmd` | Filter cells and generate QC plots |
| 03 | `03_ComputeSCT.Rmd` | Normalize (RNA + SCTransform) |
| 04 | `04_IntegrateSeuratObjects.Rmd` | Dataset integration |
| 05 | `05_ComputeModuleScores.Rmd` | Module scoring and visualization |

---

## Key Features

### Reproducibility
- Fully parameterized RMarkdown workflow
- `renv`-based dependency locking
- Conda environment support (`environment.yml`)

### Modular Architecture
Functions organized by theme:
- `qc.R`, `qc_plot_helpers.R`
- `seurat.R`, `integration_sct.R`
- `module_score_helpers.R`
- `metadata.R`, `io.R`

### Caching Strategy
- Intermediate `.rds` files reused automatically
- Avoids recomputation of:
  - normalization
  - integration
  - module scores

### Publication-Ready Outputs
- Figures saved as:
  - PDF (vector)
  - PNG (high DPI)
- Consistent themes and palettes

### Metadata

External tables define all data inputs (`data/metadata/`).

#### Data ingestion

- Sample-level information (`sample_metadata.tsv`)
- Curated list of cell cycle genes (`cell_cycle_genes.tsv`)
- Curated list of marker genes for major cell types from the literature (`literature_marker_genes.tsv`)

#### Visualization

- sample colors (`sample_color_codes.tsv`)
- cell type colors (`cell_type_color_codes_level_1.tsv`)
- stage colors (`stage_color_codes.tsv`)

---

## Inputs

Located in `data/metadata/`:
- `sample_metadata.tsv`
- `cell_cycle_genes.tsv`
- `literature_marker_genes.tsv`
- color tables (sample, stage, cell type)

---

## Outputs

Stored in used-defined directory.
- Seurat objects (`.rds`)
- Figures (UMAPs, dot plots, feature plots, clustree)
- HTML reports

---

## Usage

Run pipeline through batch run to ensure all dependencies are loaded from renv.lock and conda environment:

```bash
./biowulf/run_r.sh Rscript scripts/KnitReports.R
```

---

## Environment Setup

Please refer to docs/environment_setup.md for instructions.

---

## Important Notes

### Assay Usage
- RNA assay: normalization, module scoring, differential gene expression
- SCT/integrated assay: integration and clustering

### Color/label Consistency
All plots use external metadata tables for reproducibility.

---

## Design Principles

- Reproducibility first
- Modular and extensible
- Publication-ready by default
- Defensive programming (checks, warnings)

---

## Author

Francisco Pereira Lobo (pereiralobof2@nih.gov, franciscolobo@gmail.com)
Bioinformatics Staff Scientist
Cancer and Developmental Biology Laboratory
National Cancer Institute
National Institutes of Health

