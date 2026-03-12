# R/report_utils/seurat.R

# adds cell cycle scores
score_cell_cycle <- function(srt_obj,
                             s_features,
                             g2m_features,
                             assay = "RNA",
                             normalize_if_needed = TRUE,
                             verbose = FALSE) {
  if (!inherits(srt_obj, "Seurat")) {
    stop("score_cell_cycle() expects a Seurat object; got: ", paste(class(srt_obj), collapse = ", "))
  }
  if (!assay %in% Seurat::Assays(srt_obj)) {
    stop("Assay not present in object: ", assay)
  }

  # Clean and intersect features with the object
  s_features   <- unique(stats::na.omit(as.character(s_features)))
  g2m_features <- unique(stats::na.omit(as.character(g2m_features)))

  s_features   <- intersect(s_features, rownames(srt_obj))
  g2m_features <- intersect(g2m_features, rownames(srt_obj))

  if (length(s_features) < 10 || length(g2m_features) < 10) {
    stop("Too few cell cycle genes found in the object after intersect(). ",
         "S: ", length(s_features), ", G2M: ", length(g2m_features), ".")
  }

  Seurat::DefaultAssay(srt_obj) <- assay

  # Ensure RNA@data exists (CellCycleScoring typically uses the data slot)
  if (isTRUE(normalize_if_needed) && assay == "RNA") {
    data_mat <- Seurat::GetAssayData(srt_obj, assay = "RNA", slot = "data")
    if (ncol(data_mat) == 0 || nrow(data_mat) == 0) {
      srt_obj <- Seurat::NormalizeData(srt_obj, assay = "RNA", verbose = verbose)
    }
  }

  Seurat::CellCycleScoring(
    object = srt_obj,
    s.features = s_features,
    g2m.features = g2m_features
  )
}

# split a seurat object into a list
split_seurat <- function(srt_obj, split_by) {
  if (!split_by %in% colnames(srt_obj@meta.data)) {
    stop("split_by column not found in Seurat meta.data: ", split_by)
  }
  Seurat::SplitObject(srt_obj, split.by = split_by)
}


#sct_transform_list <- function(obj_list,
#                               method = "glmGamPoi",
#                               vars_to_regress = c("percent.mito", "S.Score", "G2M.Score"),
#                               return_only_var_genes = FALSE) {
#
#  lapply(
#    X = obj_list,
#    FUN = Seurat::SCTransform,
#    method = method,
#    vars.to.regress = vars_to_regress,
#    return.only.var.genes = return_only_var_genes
#  )
#}

sct_transform_list <- function(obj_list,
                               method = "glmGamPoi",
                               vars_to_regress = c("percent.mito", "S.Score", "G2M.Score"),
                               return_only_var_genes = FALSE,
                               verbose = TRUE) {

  if (!is.list(obj_list) || length(obj_list) == 0) {
    stop("obj_list must be a non-empty list of Seurat objects.")
  }

  n <- length(obj_list)
  nms <- names(obj_list)
  if (is.null(nms) || any(!nzchar(nms))) {
    nms <- paste0("sample_", seq_len(n))
  }

  out <- vector("list", n)
  names(out) <- nms

  for (i in seq_len(n)) {
    sid <- nms[i]
    if (isTRUE(verbose)) {
      message(sprintf("[SCTransform] Processing %s (%d/%d)", sid, i, n))
    }

    out[[i]] <- Seurat::SCTransform(
      object = obj_list[[i]],
      method = method,
      vars.to.regress = vars_to_regress,
      return.only.var.genes = return_only_var_genes,
      verbose = FALSE
    )
  }

  out
}


# run SCTransform in a list of seurat objects
#sct_transform_list <- function(obj_list,
#                               method = "glmGamPoi",
#                               vars_to_regress = c("percent.mito", "S.Score", "G2M.Score"),
#                               return_only_var_genes = FALSE,
#                               parallel = TRUE) {

#  if (!is.list(obj_list) || length(obj_list) == 0) {
#    stop("obj_list must be a non-empty list of Seurat objects.")
#  }

#  FUN <- function(x) {
#    Seurat::SCTransform(
#      x,
#      method = method,
#      vars.to.regress = vars_to_regress,
#      return.only.var.genes = return_only_var_genes
#    )
#  }

#  if (isTRUE(parallel)) {
#    if (!requireNamespace("future.apply", quietly = TRUE)) {
#      stop("Package 'future.apply' is required for parallel=TRUE.")
#    }
#    res <- future.apply::future_lapply(obj_list, FUN = FUN)
#  } else {
#    res <- lapply(obj_list, FUN)
#  }

#  names(res) <- names(obj_list)
#  res
#}


# select global features to integrate a list of seurat objects
select_global_features <- function(sample_list, nfeatures = 3000, assay = "RNA") {
  Seurat::SelectIntegrationFeatures(object.list = sample_list, nfeatures = nfeatures, assay = assay)
}

# integrate biological replicates within a list of seurat objects using SCT assay
integrate_replicates_within_group <- function(sample_list,
                                              group_obj,
                                              group_col,
                                              sample_id_col,
                                              global_features,
                                              dims = 1:40,
                                              reduction = "rpca") {

  if (!group_col %in% colnames(group_obj@meta.data)) {
    stop("group_col not found in meta.data: ", group_col)
  }
  if (!sample_id_col %in% colnames(group_obj@meta.data)) {
    stop("sample_id_col not found in meta.data: ", sample_id_col)
  }

  groups <- Seurat::SplitObject(group_obj, split.by = group_col)

  integrated_groups <- lapply(names(groups), function(gname) {
    grp <- groups[[gname]]

    samps <- unique(as.character(grp@meta.data[[sample_id_col]]))
    subs  <- sample_list[samps]

    if (length(subs) == 0) stop("No samples found for group: ", gname)

    subs <- Seurat::PrepSCTIntegration(
      object.list = subs,
      anchor.features = global_features,
      verbose = FALSE
    )

    subs <- lapply(
      subs,
      Seurat::RunPCA,
      assay = "SCT",
      features = global_features,
      verbose = FALSE
    )

    anchors <- Seurat::FindIntegrationAnchors(
      object.list = subs,
      normalization.method = "SCT",
      anchor.features = global_features,
      reduction = reduction,
      dims = dims
    )

    integ <- Seurat::IntegrateData(
      anchorset = anchors,
      normalization.method = "SCT",
      dims = dims
    )

    Seurat::DefaultAssay(integ) <- "integrated"
    integ
  })

  names(integrated_groups) <- names(groups)
  integrated_groups
}

diet_integrated_objects <- function(integrated_list,
                                    assays = "integrated",
                                    dimreducs = "pca") {
  lapply(integrated_list, function(obj) {
    Seurat::DietSeurat(obj, assays = assays, dimreducs = dimreducs)
  })
}

# merge a list of integrated seurat objects
merge_integrated_list <- function(obj_list, project = "Integration") {
  if (length(obj_list) == 0) stop("obj_list is empty.")
  if (length(obj_list) == 1) return(obj_list[[1]])

  Seurat::merge(
    x = obj_list[[1]],
    y = obj_list[-1],
    add.cell.ids = names(obj_list),
    project = project
  )
}


# run downstream analyses for an integrated seurat object
run_integrated_reductions <- function(srt_obj, dims_pca = 1:30, verbose = FALSE) {
  Seurat::DefaultAssay(srt_obj) <- "integrated"
  srt_obj <- Seurat::FindVariableFeatures(srt_obj)
  srt_obj <- Seurat::ScaleData(srt_obj, verbose = verbose)
  srt_obj <- Seurat::RunPCA(srt_obj, verbose = verbose)
  srt_obj <- Seurat::RunUMAP(srt_obj, dims = dims_pca)
  srt_obj
}

# find neighbors and clusters for a seurat object
run_graph_clustering <- function(srt_obj,
                                 dims_neighbors = 1:40,
                                 resolutions = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2, 3, 4),
                                 algorithm = 4,
                                 seed = 42,
                                 verbose = FALSE) {

  srt_obj <- Seurat::FindNeighbors(
    srt_obj,
    reduction = "pca",
    dims = dims_neighbors,
    verbose = verbose
  )

  srt_obj <- Seurat::FindClusters(
    srt_obj,
    resolution = resolutions,
    algorithm = algorithm,
    random.seed = seed
  )

  srt_obj
}

add_module_scores <- function(srt, gene_lists, assay = DefaultAssay(srt), slot = "data") {
  stopifnot(is.list(gene_lists), length(gene_lists) > 0)
  if (is.null(names(gene_lists)) || any(names(gene_lists) == "")) {
    names(gene_lists) <- paste0("Module", seq_along(gene_lists))
  }
  present <- rownames(GetAssayData(srt, assay = assay, slot = slot))
  for (nm in names(gene_lists)) {
    genes <- unique(gene_lists[[nm]])
    genes <- genes[genes %in% present]           # keep only genes in the assay
    #    if (length(genes) < 3L) {
    #      warning(nm, ": <3 genes found in assay ", assay, "; skipped")
    #      next
    #    }
    srt <- AddModuleScore(srt, features = list(genes), name = nm, assay = assay)
    colnames(srt@meta.data)[ncol(srt@meta.data)] <- nm  # rename "<name>1" → "<name>"
  }
  srt
}
