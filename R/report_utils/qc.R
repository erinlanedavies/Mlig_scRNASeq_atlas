# R/report_utils/qc.R

# Compute percentage of mitochondrial RNA
# Priority:
# 1) If mito_genes (vector) is non-empty, use it as a list of gene IDs.
# 2) Else, if mito_species is either human or mouse, use species-specific regular expressions.
# 3) Else, skip with a warning.
compute_percent_mito <- function(obj,
                                 mito_genes,
                                 mito_species,
                                 mito_colname = "percent.mito") {
  mito_genes <- unique(mito_genes)

  if (!is.null(mito_genes) && length(mito_genes) > 0) {
    present_mito <- intersect(mito_genes, rownames(obj))
    if (length(present_mito) == 0) {
      warning("None of the provided mito_genes found in rownames; skipping mito calculation.")
      return(obj)
    }
    obj[[mito_colname]] <- Seurat::PercentageFeatureSet(obj, features = present_mito)
    return(obj)
  }

  pattern <- NULL
  species_lower <- tolower(mito_species)

  if (species_lower %in% c("human", "homo sapiens", "homo_sapiens")) {
    pattern <- "^MT-"
  } else if (species_lower %in% c("mouse", "mus musculus", "mus_musculus")) {
    pattern <- "^mt-"
  } else {
    warning(
      "mito_species '", mito_species,
      "' not recognized and mito_genes not provided; skipping mito calculation."
    )
    return(obj)
  }

  obj[[mito_colname]] <- Seurat::PercentageFeatureSet(obj, pattern = pattern)
  return(obj)
}


# Run DoubletFinder-based doublet detection with internal preprocessing.
# Steps:
# - SCTransform
# - PCA, UMAP, neighbors, clusters
# - paramSweep + pK selection
# - doubletFinder
# - clean metadata columns
# - optional DietSeurat

DoubletDetection <- function(srt_obj,
                             params = list(),
                             ...) {

  args <- params

  srt_obj_tmp <- srt_obj
  srt_obj_tmp <- Seurat::SCTransform(srt_obj_tmp)
  srt_obj_tmp <- Seurat::RunPCA(srt_obj_tmp)
  srt_obj_tmp <- Seurat::RunUMAP(srt_obj_tmp, dims = 1:40)
  srt_obj_tmp <- Seurat::FindNeighbors(srt_obj_tmp, reduction = "pca", dims = 1:40)
  srt_obj_tmp <- Seurat::FindClusters(srt_obj_tmp, verbose = FALSE)

  # Parameter sweep to choose best pK
  sweep.res.list <- DoubletFinder::paramSweep(srt_obj_tmp, PCs = 1:40, sct = TRUE)
  sweep.stats    <- DoubletFinder::summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn          <- DoubletFinder::find.pK(sweep.stats)

  pK <- as.numeric(as.character(bcmvn$pK))
  BCmetric <- bcmvn$BCmetric
  pK_choose <- pK[which.max(BCmetric)]

  # Expected doublets (7.5 percent heuristic)
  nExp_poi <- round(0.075 * nrow(srt_obj_tmp@meta.data))

  srt_obj_tmp <- DoubletFinder::doubletFinder(
    srt_obj_tmp,
    PCs = 1:40,
    pN = 0.25,
    pK = pK_choose,
    nExp = nExp_poi,
    sct = TRUE
  )

  # Clean DF columns and keep only final classification
  md_cols <- colnames(srt_obj_tmp@meta.data)

  pann_col <- md_cols[grepl("^pANN_", md_cols)]
  if (length(pann_col) > 0) {
    srt_obj_tmp@meta.data[, pann_col] <- NULL
  }

  df_col <- md_cols[grepl("^DF.classifications_", md_cols)]
  if (length(df_col) == 1) {
    srt_obj_tmp[["DoubletStatus"]] <- srt_obj_tmp@meta.data[, df_col]
    srt_obj_tmp@meta.data[, df_col] <- NULL
  }

  # Optional slimming of object
  if (isTRUE(args$diet_seurat)) {
    Seurat::DefaultAssay(srt_obj_tmp) <- "RNA"
    srt_obj_tmp <- Seurat::DietSeurat(srt_obj_tmp, assays = "RNA", dimreducs = NULL)

    drop_cols <- c("nCount_SCT", "nFeature_SCT", "SCT_snn_res.0.8")
    for (cc in drop_cols) {
      if (cc %in% colnames(srt_obj_tmp@meta.data)) {
        srt_obj_tmp@meta.data[[cc]] <- NULL
      }
    }
  }

  return(srt_obj_tmp)
}

# Compute outlier flags using scater::isOutlier for standard QC metrics.
# Expects:
# params list with:
#   QC_proxy_vars
#   nmad_lower
#   nmad_upper
#   log

QCReportFilter <- function(srt_obj,
                           params = list(),
                           ...) {

  options(future.globals.maxSize = 2 * 1024^4)

  nmads <- list(
    percent.mito = c(params$nmad_lower[1], params$nmad_upper[1]),
    nCount_RNA   = c(params$nmad_lower[2], params$nmad_upper[2]),
    nFeature_RNA = c(params$nmad_lower[3], params$nmad_upper[3])
  )

  qc.nCount_lower <- scater::isOutlier(
    srt_obj$nCount_RNA,
    nmads = nmads$nCount_RNA[1],
    log   = params$log,
    type  = "lower"
  )

  qc.nCount_upper <- scater::isOutlier(
    srt_obj$nCount_RNA,
    nmads = nmads$nCount_RNA[2],
    log   = params$log,
    type  = "higher"
  )

  qc.nFeature_lower <- scater::isOutlier(
    srt_obj$nFeature_RNA,
    nmads = nmads$nFeature_RNA[1],
    log   = params$log,
    type  = "lower"
  )

  qc.nFeature_upper <- scater::isOutlier(
    srt_obj$nFeature_RNA,
    nmads = nmads$nFeature_RNA[2],
    log   = params$log,
    type  = "higher"
  )

  qc.mito_lower <- scater::isOutlier(
    srt_obj$percent.mito,
    nmads = nmads$percent.mito[1],
    log   = params$log,
    type  = "lower"
  )

  qc.mito_upper <- scater::isOutlier(
    srt_obj$percent.mito,
    nmads = nmads$percent.mito[2],
    log   = params$log,
    type  = "higher"
  )

  discard <- qc.nCount_lower | qc.nCount_upper |
             qc.nFeature_lower | qc.nFeature_upper |
             qc.mito_lower | qc.mito_upper

  srt_obj$isOutlier <- discard
  return(srt_obj)
}

