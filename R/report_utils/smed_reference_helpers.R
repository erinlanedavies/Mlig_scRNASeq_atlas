# smed_reference_helpers.R
# Helper functions for building an S. mediterranea Seurat reference atlas

library(Seurat)
library(dplyr)
library(stringr)
library(ggplot2)
library(scCustomize)

get_smed_cell_cycle_features <- function() {
  S_genes <- c(
    "dd_Smed_v4_5688_0_1",
    "dd_Smed_v4_8206_0_1",
    "dd_Smed_v4_1651_0_1",
    "dd_Smed_v4_15389_0_1",
    "dd_Smed_v4_11558_0_1",
    "dd_Smed_v4_5663_0_1",
    "dd_Smed_v4_4379_0_1",
    "dd_Smed_v4_8626_0_1",
    "dd_Smed_v4_16942_0_1",
    "dd_Smed_v4_14547_0_1",
    "dd_Smed_v4_17862_0_1",
    "dd_Smed_v4_5764_0_1",
    "dd_Smed_v4_1260_0_1",
    "dd_Smed_v4_5956_0_1",
    "dd_Smed_v4_1568_0_1",
    "dd_Smed_v4_15465_0_1",
    "dd_Smed_v4_4341_0_1",
    "dd_Smed_v4_9628_0_1",
    "dd_Smed_v4_11434_0_1",
    "dd_Smed_v4_20567_0_1",
    "dd_Smed_v4_18778_0_1",
    "dd_Smed_v4_12580_0_1"
  )

  G2M_genes <- c(
    "dd_Smed_v4_6668_0_1",
    "dd_Smed_v4_5865_0_1",
    "dd_Smed_v4_13972_0_1",
    "dd_Smed_v4_88923_0_1",
    "dd_Smed_v4_14243_0_1",
    "dd_Smed_v4_2016_0_1",
    "dd_Smed_v4_14261_0_1",
    "dd_Smed_v4_970_0_1",
    "dd_Smed_v4_15869_0_1",
    "dd_Smed_v4_17887_0_1",
    "dd_Smed_v4_7553_0_1",
    "dd_Smed_v4_9585_0_1",
    "dd_Smed_v4_12391_0_1"
  )

  list(
    s_features = stringr::str_replace_all(S_genes, "_", "-"),
    g2m_features = stringr::str_replace_all(G2M_genes, "_", "-")
  )
}

read_smed_dge_matrix <- function(path) {
  if (!file.exists(path)) {
    stop("DGE matrix file not found: ", path)
  }

  read.table(path, header = TRUE, check.names = FALSE)
}

read_smed_dge_matrices <- function(principal_path,
                                   brain_path,
                                   sexual_path = NULL,
                                   include_sexual = FALSE) {
  out <- list(
    Principal = read_smed_dge_matrix(principal_path),
    Brain = read_smed_dge_matrix(brain_path)
  )

  if (isTRUE(include_sexual)) {
    if (is.null(sexual_path)) {
      stop("sexual_path must be provided when include_sexual = TRUE.")
    }
    out$Sexual <- read_smed_dge_matrix(sexual_path)
  }

  out
}

create_smed_seurat_object <- function(counts,
                                      sample_id,
                                      min.cells = 3,
                                      min.features = 200) {
  obj <- Seurat::CreateSeuratObject(
    counts = counts,
    min.cells = min.cells,
    min.features = min.features,
    project = sample_id
  )

  obj$sample_ID <- sample_id
  obj
}

create_smed_seurat_list <- function(dge_list,
                                    min.cells = 3,
                                    min.features = 200) {
  if (!is.list(dge_list) || length(dge_list) == 0) {
    stop("dge_list must be a non-empty list.")
  }

  out <- lapply(names(dge_list), function(nm) {
    create_smed_seurat_object(
      counts = dge_list[[nm]],
      sample_id = nm,
      min.cells = min.cells,
      min.features = min.features
    )
  })

  names(out) <- names(dge_list)
  out
}

normalize_and_score_smed_cell_cycle <- function(obj,
                                                s_features,
                                                g2m_features,
                                                assay = "RNA") {
  if (!inherits(obj, "Seurat")) {
    stop("obj must be a Seurat object.")
  }

  if (!assay %in% Seurat::Assays(obj)) {
    stop("Assay not found in object: ", assay)
  }

  Seurat::DefaultAssay(obj) <- assay
  obj <- Seurat::NormalizeData(obj, assay = assay, verbose = FALSE)
  obj <- Seurat::CellCycleScoring(
    object = obj,
    s.features = s_features,
    g2m.features = g2m_features
  )

  obj
}

merge_smed_reference_objects <- function(obj_list,
                                         project = "scRNA_Seq_S_mediterranea") {
  if (!is.list(obj_list) || length(obj_list) == 0) {
    stop("obj_list must be a non-empty list.")
  }

  if (length(obj_list) == 1) {
    return(obj_list[[1]])
  }

  merge(
    x = obj_list[[1]],
    y = obj_list[-1],
    add.cell.ids = names(obj_list),
    project = project
  )
}

# run_smed_reference_workflow <- function(obj,
#                                         seed = 42,
#                                         npcs = 100,
#                                         integration_dims = 1:40,
#                                         clustering_resolution = 0.6,
#                                         k_param = 20) {
#   if (!inherits(obj, "Seurat")) {
#     stop("obj must be a Seurat object.")
#   }
# 
#   obj <- Seurat::SCTransform(
#     object = obj,
#     vst.flavor = "v2",
#     vars.to.regress = c("G2M.Score", "S.Score"),
#     seed.use = seed,
#     verbose = FALSE
#   )
# 
#   obj <- Seurat::RunPCA(
#     object = obj,
#     npcs = npcs,
#     verbose = TRUE,
#     seed.use = seed
#   )
# 
#   obj <- Seurat::IntegrateLayers(
#     object = obj,
#     method = Seurat::CCAIntegration,
#     orig.reduction = "pca",
#     new.reduction = "integrated.cca",
#     normalization.method = "SCT",
#     assay = "SCT",
#     verbose = FALSE
#   )
# 
#   obj <- Seurat::RunUMAP(
#     object = obj,
#     reduction = "integrated.cca",
#     reduction.name = "umap.integrated.cca",
#     dims = integration_dims,
#     return.model = TRUE,
#     seed.use = seed
#   )
# 
#   obj <- Seurat::FindNeighbors(
#     object = obj,
#     reduction = "integrated.cca",
#     dims = integration_dims,
#     k.param = k_param,
#     return.neighbor = TRUE,
#     nn.method = "annoy",
#     annoy.metric = "cosine",
#     verbose = FALSE
#   )
# 
#   obj <- Seurat::FindNeighbors(
#     object = obj,
#     reduction = "integrated.cca",
#     dims = integration_dims,
#     k.param = k_param,
#     compute.SNN = TRUE,
#     nn.method = "annoy",
#     annoy.metric = "cosine",
#     verbose = FALSE
#   )
# 
#   obj <- Seurat::FindClusters(
#     object = obj,
#     resolution = clustering_resolution,
#     algorithm = 4,
#     method = "igraph",
#     random.seed = seed
#   )
# 
#   obj
# }

run_smed_reference_workflow <- function(obj,
                                        seed = 42,
                                        npcs = 100,
                                        integration_dims = 1:40,
                                        clustering_resolution = 0.6,
                                        k_param = 20) {
  if (!inherits(obj, "Seurat")) {
    stop("obj must be a Seurat object.")
  }
  
  obj <- Seurat::SCTransform(
    object = obj,
    vst.flavor = "v2",
    vars.to.regress = c("G2M.Score", "S.Score"),
    seed.use = seed,
    verbose = FALSE
  )
  
  obj <- Seurat::RunPCA(
    object = obj,
    npcs = npcs,
    verbose = TRUE,
    seed.use = seed
  )
  
  # obj <- Seurat::IntegrateLayers(
  #   object = obj,
  #   method = Seurat::CCAIntegration,
  #   orig.reduction = "pca",
  #   new.reduction = "integrated.cca",
  #   normalization.method = "SCT",
  #   assay = "SCT",
  #   verbose = FALSE
  # )
  
  obj <- Seurat::RunUMAP(
    object = obj,
    dims = integration_dims,
    return.model = TRUE,
    seed.use = seed
  )
  
  obj <- Seurat::FindNeighbors(
    object = obj,
    dims = integration_dims,
    k.param = k_param,
    return.neighbor = TRUE,
    nn.method = "annoy",
    annoy.metric = "cosine",
    verbose = FALSE
  )
  
  obj <- Seurat::FindNeighbors(
    object = obj,
    dims = integration_dims,
    k.param = k_param,
    compute.SNN = TRUE,
    nn.method = "annoy",
    annoy.metric = "cosine",
    verbose = FALSE
  )
  
  obj <- Seurat::FindClusters(
    object = obj,
    resolution = clustering_resolution,
    algorithm = 4,
    method = "igraph",
    random.seed = seed
  )
  
  obj
}


read_smed_cell_annotation_table <- function(path) {
  if (!file.exists(path)) {
    stop("Cell annotation file not found: ", path)
  }

  read.table(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
}

strip_cell_prefix <- function(x, regex = "^[^_]*_") {
  sub(regex, "", x)
}

attach_smed_cell_labels <- function(obj,
                                    cell_metadata,
                                    cell_id_col = "Cell.ID",
                                    metadata_cols,
                                    strip_prefix_regex = "^[^_]*_") {
  if (!inherits(obj, "Seurat")) {
    stop("obj must be a Seurat object.")
  }
  
  if (!is.data.frame(cell_metadata)) {
    stop("cell_metadata must be a data.frame.")
  }
  
  if (!cell_id_col %in% colnames(cell_metadata)) {
    stop("cell_id_col not found in cell_metadata: ", cell_id_col)
  }
  
  if (missing(metadata_cols) || length(metadata_cols) == 0) {
    stop("metadata_cols must contain at least one column name.")
  }
  
  if (!all(metadata_cols %in% colnames(cell_metadata))) {
    missing_cols <- setdiff(metadata_cols, colnames(cell_metadata))
    stop(
      "The following metadata columns were not found in cell_metadata: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  # Rename cells in the Seurat object to match the external annotation table
  new_names <- strip_cell_prefix(Seurat::Cells(obj), regex = strip_prefix_regex)
  obj <- Seurat::RenameCells(obj, new.names = new_names)
  
  # Keep only the requested metadata columns plus the cell ID column
  meta_df <- cell_metadata[, c(cell_id_col, metadata_cols), drop = FALSE]
  
  # Remove duplicated cell IDs, keeping the first occurrence
  if (anyDuplicated(meta_df[[cell_id_col]]) > 0) {
    warning("Duplicated cell IDs found in cell_metadata. Keeping the first occurrence for each cell ID.")
    meta_df <- meta_df[!duplicated(meta_df[[cell_id_col]]), , drop = FALSE]
  }
  
  rownames(meta_df) <- as.character(meta_df[[cell_id_col]])
  meta_df[[cell_id_col]] <- NULL
  
  # Reorder to match Seurat cell order
  matched_meta <- meta_df[Seurat::Cells(obj), , drop = FALSE]
  
  # Add all requested metadata columns, preserving original column names
  obj <- Seurat::AddMetaData(obj, metadata = matched_meta)
  
  obj
}
plot_smed_marker_feature <- function(obj,
                                     feature,
                                     title = NULL,
                                     aspect_ratio = 1,
                                     pt.size = 0.1) {
  p <- scCustomize::FeaturePlot_scCustom(
    seurat_object = obj,
    features = feature,
    aspect_ratio = aspect_ratio,
    pt.size = pt.size
  )

  if (!is.null(title) && nzchar(title)) {
    p <- p + ggplot2::ggtitle(title)
  }

  p
}

plot_smed_marker_dotplot <- function(obj,
                                     features,
                                     group_by) {
  scCustomize::Clustered_DotPlot(
    seurat_object = obj,
    features = features,
    group.by = group_by
  )
}

save_smed_reference_plot <- function(plot_obj,
                                     file_stub,
                                     out_dir,
                                     width = 10,
                                     height = 8,
                                     dpi = 600) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  pdf_file <- file.path(out_dir, paste0(file_stub, ".pdf"))
  png_file <- file.path(out_dir, paste0(file_stub, ".png"))
  
  # Handle Clustered_DotPlot() when plot_km_elbow = TRUE
  if (is.list(plot_obj) && !inherits(plot_obj, "ggplot")) {
    # Prefer the first ComplexHeatmap-like element if present
    ht_idx <- which(vapply(
      plot_obj,
      function(x) {
        inherits(x, "Heatmap") ||
          inherits(x, "HeatmapList") ||
          inherits(x, "AdditiveUnit")
      },
      logical(1)
    ))
    
    if (length(ht_idx) > 0) {
      plot_obj <- plot_obj[[ht_idx[1]]]
    }
  }
  
  is_complex_heatmap <-
    inherits(plot_obj, "Heatmap") ||
    inherits(plot_obj, "HeatmapList") ||
    inherits(plot_obj, "AdditiveUnit")
  
  is_ggplot_like <-
    inherits(plot_obj, "gg") ||
    inherits(plot_obj, "ggplot") ||
    inherits(plot_obj, "patchwork")
  
  if (is_complex_heatmap) {
    grDevices::pdf(pdf_file, width = width, height = height)
    ComplexHeatmap::draw(plot_obj)
    grDevices::dev.off()
    
    grDevices::png(
      filename = png_file,
      width = width,
      height = height,
      units = "in",
      res = dpi,
      bg = "white"
    )
    ComplexHeatmap::draw(plot_obj)
    grDevices::dev.off()
    
  } else if (is_ggplot_like) {
    ggplot2::ggsave(
      filename = pdf_file,
      plot = plot_obj,
      width = width,
      height = height,
      units = "in",
      dpi = dpi,
      bg = "white"
    )
    
    ggplot2::ggsave(
      filename = png_file,
      plot = plot_obj,
      width = width,
      height = height,
      units = "in",
      dpi = dpi,
      bg = "white"
    )
  } else {
    stop(
      "Unsupported plot object class: ",
      paste(class(plot_obj), collapse = ", ")
    )
  }
}
read_marker_dotplot_table <- function(path,
                                      celltype_col = "CellType",
                                      gene_col = "GeneID",
                                      sanitize_gene_ids = TRUE,
                                      unique_rows = TRUE) {
  if (!file.exists(path)) {
    stop("Marker dot plot table not found: ", path)
  }
  
  df <- read.table(
    file = path,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  # Validate required columns
  if (!all(c(celltype_col, gene_col) %in% colnames(df))) {
    stop(
      "Marker table must contain columns: ",
      celltype_col, " and ", gene_col
    )
  }
  
  # Clean
  df[[celltype_col]] <- as.character(df[[celltype_col]])
  df[[gene_col]] <- as.character(df[[gene_col]])
  
  df <- df[
    !is.na(df[[celltype_col]]) &
      !is.na(df[[gene_col]]) &
      nzchar(df[[celltype_col]]) &
      nzchar(df[[gene_col]]),
    ,
    drop = FALSE
  ]
  
  # Optional normalization (important for Smed IDs)
  if (isTRUE(sanitize_gene_ids)) {
    df[[gene_col]] <- gsub("_", "-", df[[gene_col]])
  }
  
  if (isTRUE(unique_rows)) {
    df <- unique(df)
  }
  
  rownames(df) <- NULL
  df
}

build_marker_dotplot_features <- function(df,
                                          celltype_col = "CellType",
                                          gene_col = "GeneID") {
  split_list <- split(df[[gene_col]], df[[celltype_col]])
  
  split_list <- lapply(split_list, function(x) unique(as.character(x)))
  
  list(
    features = unique(unlist(split_list)),
    grouped_features = split_list
  )
}

build_marker_feature_plot_table <- function(df,
                                            celltype_col = "CellType",
                                            gene_col = "GeneID") {
  if (!is.data.frame(df)) {
    stop("df must be a data.frame.")
  }
  
  if (!all(c(celltype_col, gene_col) %in% colnames(df))) {
    stop("df must contain columns: ", celltype_col, " and ", gene_col)
  }
  
  out <- df[, c(celltype_col, gene_col), drop = FALSE]
  colnames(out) <- c("CellType", "GeneID")
  
  out$CellType <- as.character(out$CellType)
  out$GeneID <- as.character(out$GeneID)
  
  out <- out[
    !is.na(out$CellType) &
      !is.na(out$GeneID) &
      nzchar(out$CellType) &
      nzchar(out$GeneID),
    ,
    drop = FALSE
  ]
  
  out <- unique(out)
  rownames(out) <- NULL
  out
}