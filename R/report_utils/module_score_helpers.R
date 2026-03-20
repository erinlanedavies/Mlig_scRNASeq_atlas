# module_score_helpers.R
# Helper functions for module score computation and visualization

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(stringr)

read_marker_table <- function(path,
                              gene_col = "Gene_ID",
                              celltype_col = "Cell_tissue_type") {
  if (!file.exists(path)) {
    stop("Marker gene table not found: ", path)
  }

  df <- read.table(
    file = path,
    sep = "\t",
    header = TRUE,
    quote = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (!all(c(gene_col, celltype_col) %in% colnames(df))) {
    stop("Marker table must contain columns: ", gene_col, " and ", celltype_col)
  }

  df[[gene_col]] <- gsub("_", "-", df[[gene_col]])
  df
}

default_celltype_name_map <- function() {
  c(
    "anchor_cell_of_adhesive_organ" = "anchor_cell_of_adhesive_organ",
    "cement_gland_cell" = "cement_gland_cell",
    "epidermal_secretory_cell" = "epidermal_secretory",
    "female_antrum" = "female_antrum",
    "female_germline_ovary" = "female_germline_ovary",
    "gut" = "gut",
    "gut_specialized_cell" = "gut_specialized",
    "male_germline_testis" = "male_germline_testis",
    "muscle" = "muscle",
    "nervous_system" = "nervous_system",
    "prostate" = "prostate",
    "rhabdite-containing_cell" = "rhabdite",
    "secretory_cell_of_adhesive_organ" = "secretory_cell_of_adhesive_organ",
    "stem_cells_progenitors" = "progenitor",
    "epidermal" = "epidermal"
  )
}

build_marker_gene_list <- function(marker_df,
                                   gene_col = "Gene_ID",
                                   celltype_col = "Cell_tissue_type",
                                   name_map = default_celltype_name_map()) {
  split_genes <- split(marker_df[[gene_col]], marker_df[[celltype_col]])
  split_genes <- lapply(split_genes, function(x) unique(stats::na.omit(as.character(x))))

  out <- list()
  for (nm in names(split_genes)) {
    if (nm %in% names(name_map)) {
      out[[name_map[[nm]]]] <- split_genes[[nm]]
    }
  }

  out
}

summarize_marker_table <- function(marker_df,
                                   gene_col = "Gene_ID",
                                   celltype_col = "Cell_tissue_type") {
  marker_df %>%
    dplyr::group_by(.data[[celltype_col]]) %>%
    dplyr::summarise(
      n_marker_genes = dplyr::n_distinct(.data[[gene_col]]),
      .groups = "drop"
    ) %>%
    dplyr::arrange(desc(.data$n_marker_genes), .data[[celltype_col]])
}

module_score_overlap_table <- function(obj,
                                       gene_lists,
                                       object_name = "object") {
  data.frame(
    object = object_name,
    module = names(gene_lists),
    n_markers_total = vapply(gene_lists, length, integer(1)),
    n_markers_present = vapply(gene_lists, function(g) sum(g %in% rownames(obj)), integer(1)),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::mutate(
      pct_present = round(100 * .data$n_markers_present / .data$n_markers_total, 2)
    )
}

add_module_scores_named <- function(obj,
                                    gene_lists,
                                    assay = "RNA",
                                    prefix = "MS_",
                                    nbin = 24,
                                    ctrl = 100,
                                    seed = 1) {
  if (!inherits(obj, "Seurat")) {
    stop("obj must be a Seurat object.")
  }

  if (!assay %in% Assays(obj)) {
    stop("Assay not found in object: ", assay)
  }

  DefaultAssay(obj) <- assay

  cleaned <- lapply(gene_lists, function(g) unique(g[g %in% rownames(obj)]))
  cleaned <- cleaned[vapply(cleaned, length, integer(1)) > 0]

  if (length(cleaned) == 0) {
    stop("None of the gene lists overlap with the object features.")
  }

  set.seed(seed)
  obj <- Seurat::AddModuleScore(
    object = obj,
    features = unname(cleaned),
    assay = assay,
    name = prefix,
    nbin = nbin,
    ctrl = ctrl,
    search = FALSE
  )

  new_cols <- paste0(prefix, seq_along(cleaned))
  final_cols <- paste0(prefix, names(cleaned))
  colnames(obj@meta.data)[match(new_cols, colnames(obj@meta.data))] <- final_cols

  attr(obj, "module_score_columns") <- final_cols
  obj
}

find_cluster_column <- function(obj, resolution = 0.6) {
  res_str <- as.character(resolution)
  pattern <- paste0("_snn_res\\.", gsub("\\.", "\\\\.", res_str), "$")
  hits <- grep(pattern, colnames(obj@meta.data), value = TRUE)
  if (length(hits) > 0) return(hits[1])

  fallback <- grep("_snn_res\\.", colnames(obj@meta.data), value = TRUE)
  if (length(fallback) > 0) return(fallback[1])

  return(NA_character_)
}

plot_module_feature_grid <- function(obj,
                                     module_cols,
                                     reduction = "umap",
                                     ncol = 5,
                                     pt.size = 0.05,
                                     min.cutoff = "q05",
                                     max.cutoff = "q95",
                                     order = TRUE,
                                     cols = c("lightgrey", "darkred")) {
  if (!inherits(obj, "Seurat")) {
    stop("obj must be a Seurat object.")
  }

  if (!"RNA" %in% Assays(obj)) {
    stop("RNA assay not found in object.")
  }

  DefaultAssay(obj) <- "RNA"

  plots <- Seurat::FeaturePlot(
    object = obj,
    features = module_cols,
    reduction = reduction,
    combine = FALSE,
    pt.size = pt.size,
    min.cutoff = min.cutoff,
    max.cutoff = max.cutoff,
    order = order,
    cols = cols
  )

  plots <- lapply(seq_along(plots), function(i) {
    plots[[i]] +
      ggplot2::ggtitle(module_cols[i]) +
      ggplot2::coord_fixed() +
      ggplot2::theme_classic(base_size = 11) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
        axis.title = ggplot2::element_blank(),
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        aspect.ratio = 1
      )
  })

  patchwork::wrap_plots(plots, ncol = ncol)
}

plot_module_dotplot <- function(obj,
                                module_cols,
                                group.by,
                                cols = c("lightgrey", "darkred"),
                                dot.scale = 6) {
  DefaultAssay(obj) <- "RNA"

  Seurat::DotPlot(
    object = obj,
    features = module_cols,
    group.by = group.by,
    cols = cols,
    dot.scale = dot.scale,
    assay = "RNA"
  ) +
    Seurat::RotatedAxis() +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      axis.title = ggplot2::element_text(face = "bold"),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
      axis.text.y = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5)
    )
}

save_publication_plot <- function(plot_obj,
                                  file_stub,
                                  out_dir,
                                  width = 12,
                                  height = 8,
                                  dpi = 600) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  ggplot2::ggsave(
    filename = file.path(out_dir, paste0(file_stub, ".pdf")),
    plot = plot_obj,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    bg = "white"
  )

  ggplot2::ggsave(
    filename = file.path(out_dir, paste0(file_stub, ".png")),
    plot = plot_obj,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    bg = "white"
  )
}

