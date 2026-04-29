# go_enrichment_helpers.R
# Helper functions for FindAllMarkers-based GO enrichment using topGO
# and gprofiler-style plotting.

library(dplyr)
library(stringr)
library(ggplot2)
library(forcats)
library(topGO)

read_gene2go_map <- function(path,
                             gene_col = 1,
                             go_col = 2,
                             go_sep = ", ") {
  if (!file.exists(path)) {
    stop("gene2GO map not found: ", path)
  }

  geneID2GOs <- read.table(
    path,
    sep = "\t",
    header = FALSE,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = ""
  )

  if (ncol(geneID2GOs) < max(gene_col, go_col)) {
    stop("gene2GO map has fewer columns than expected.")
  }

  gene_ids <- as.character(geneID2GOs[[gene_col]])
  go_strings <- as.character(geneID2GOs[[go_col]])

  keep <- !is.na(gene_ids) & !is.na(go_strings) &
    nzchar(gene_ids) & nzchar(go_strings)

  gene_ids <- gene_ids[keep]
  go_strings <- go_strings[keep]

  gene2GO <- strsplit(go_strings, go_sep, fixed = TRUE)
  gene2GO <- lapply(gene2GO, function(x) {
    x <- trimws(x)
    unique(x[nzchar(x)])
  })

  names(gene2GO) <- gene_ids
  gene2GO
}

read_marker_table_for_go <- function(path,
                                     cluster_col = "cluster",
                                     gene_col = "gene",
                                     padj_col = "p_val_adj",
                                     logfc_col = "avg_log2FC") {
  if (!file.exists(path)) {
    stop("Marker table not found: ", path)
  }

  df <- read.table(
    path,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    quote = "",
    comment.char = ""
  )

  validate_marker_table_for_go(
    markers = df,
    cluster_col = cluster_col,
    gene_col = gene_col,
    padj_col = padj_col,
    logfc_col = logfc_col
  )

  df[[cluster_col]] <- as.character(df[[cluster_col]])
  df[[gene_col]] <- as.character(df[[gene_col]])
  df[[padj_col]] <- as.numeric(df[[padj_col]])
  df[[logfc_col]] <- as.numeric(df[[logfc_col]])

  df
}

validate_marker_table_for_go <- function(markers,
                                         cluster_col = "cluster",
                                         gene_col = "gene",
                                         padj_col = "p_val_adj",
                                         logfc_col = "avg_log2FC") {
  if (!is.data.frame(markers)) {
    stop("markers must be a data.frame.")
  }

  required_cols <- c(cluster_col, gene_col, padj_col, logfc_col)
  missing_cols <- setdiff(required_cols, colnames(markers))

  if (length(missing_cols) > 0) {
    stop(
      "Marker table is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  invisible(TRUE)
}
compute_or_load_findallmarkers <- function(obj,
                                           markers_file,
                                           compute_markers = TRUE,
                                           reuse_existing_markers = TRUE,
                                           assay = "RNA",
                                           ident_col,
                                           only_pos = TRUE,
                                           logfc_threshold = 0,
                                           min_pct = 0.1,
                                           test_use = "wilcox",
                                           cluster_col = "cluster",
                                           gene_col = "gene",
                                           padj_col = "p_val_adj",
                                           logfc_col = "avg_log2FC") {
  if (!inherits(obj, "Seurat")) {
    stop("obj must be a Seurat object.")
  }
  
  if (missing(ident_col) || is.null(ident_col) || !nzchar(ident_col)) {
    stop("ident_col must be provided.")
  }
  
  if (!ident_col %in% colnames(obj@meta.data)) {
    stop(
      "ident_col not found in Seurat metadata: ", ident_col,
      "\nAvailable columns:\n",
      paste(colnames(obj@meta.data), collapse = ", ")
    )
  }
  
  if (!assay %in% names(obj@assays)) {
    stop(
      "assay not found in Seurat object: ", assay,
      "\nAvailable assays:\n",
      paste(names(obj@assays), collapse = ", ")
    )
  }
  
  should_reuse <-
    isTRUE(reuse_existing_markers) &&
    !isTRUE(compute_markers) &&
    file.exists(markers_file)
  
  if (should_reuse) {
    message("Reading existing FindAllMarkers table: ", markers_file)
    
    return(
      read_marker_table_for_go(
        path = markers_file,
        cluster_col = cluster_col,
        gene_col = gene_col,
        padj_col = padj_col,
        logfc_col = logfc_col
      )
    )
  }
  
  message("Running Seurat::FindAllMarkers().")
  message("Assay: ", assay)
  message("Identity column: ", ident_col)
  
  Seurat::DefaultAssay(obj) <- assay
  Seurat::Idents(obj) <- obj@meta.data[[ident_col]]
  
  markers <- Seurat::FindAllMarkers(
    object = obj,
    only.pos = only_pos,
    logfc.threshold = logfc_threshold,
    min.pct = min_pct,
    test.use = test_use
  )
  
  if (!is.data.frame(markers)) {
    stop("FindAllMarkers() did not return a data.frame.")
  }
  
  if (nrow(markers) == 0) {
    stop(
      "FindAllMarkers() returned zero marker rows. ",
      "Try lowering findallmarkers_logfc_threshold or findallmarkers_min_pct, ",
      "or check that marker_ident_col defines valid identities."
    )
  }
  
  message("FindAllMarkers returned columns:")
  message("  ", paste(colnames(markers), collapse = ", "))
  
  # Seurat can return different fold-change column names depending on version/settings.
  if (!logfc_col %in% colnames(markers)) {
    fc_candidates <- c(
      "avg_log2FC",
      "avg_logFC",
      "avg_log2fc",
      "avg_logfc"
    )
    
    fc_hit <- fc_candidates[fc_candidates %in% colnames(markers)]
    
    if (length(fc_hit) == 0) {
      stop(
        "No recognized log fold-change column found in FindAllMarkers output.\n",
        "Expected one of: ", paste(fc_candidates, collapse = ", "), "\n",
        "Observed columns: ", paste(colnames(markers), collapse = ", ")
      )
    }
    
    message("Renaming fold-change column '", fc_hit[1], "' to '", logfc_col, "'.")
    colnames(markers)[colnames(markers) == fc_hit[1]] <- logfc_col
  }
  
  required_cols <- c(cluster_col, gene_col, padj_col, logfc_col)
  missing_cols <- setdiff(required_cols, colnames(markers))
  
  if (length(missing_cols) > 0) {
    stop(
      "FindAllMarkers output is missing required columns: ",
      paste(missing_cols, collapse = ", "),
      "\nObserved columns: ",
      paste(colnames(markers), collapse = ", ")
    )
  }
  
  markers[[cluster_col]] <- as.character(markers[[cluster_col]])
  markers[[gene_col]] <- as.character(markers[[gene_col]])
  markers[[padj_col]] <- as.numeric(markers[[padj_col]])
  markers[[logfc_col]] <- as.numeric(markers[[logfc_col]])
  
  dir.create(dirname(markers_file), recursive = TRUE, showWarnings = FALSE)
  
  write.table(
    markers,
    file = markers_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  message("FindAllMarkers table saved to: ", markers_file)
  
  markers
}
summarize_marker_selection_for_go <- function(markers,
                                              cluster_col = "cluster",
                                              gene_col = "gene",
                                              padj_col = "p_val_adj",
                                              logfc_col = "avg_log2FC",
                                              padj_threshold = 1e-5,
                                              logfc_threshold = 2) {
  validate_marker_table_for_go(markers, cluster_col, gene_col, padj_col, logfc_col)

  markers %>%
    dplyr::group_by(.data[[cluster_col]]) %>%
    dplyr::summarise(
      n_marker_rows = dplyr::n(),
      n_unique_genes = dplyr::n_distinct(.data[[gene_col]]),
      n_selected_genes = dplyr::n_distinct(
        .data[[gene_col]][
          !is.na(.data[[padj_col]]) &
            !is.na(.data[[logfc_col]]) &
            .data[[padj_col]] < padj_threshold &
            .data[[logfc_col]] > logfc_threshold
        ]
      ),
      .groups = "drop"
    ) %>%
    dplyr::arrange(.data[[cluster_col]])
}

get_selected_marker_genes_for_go <- function(markers,
                                             cluster_col = "cluster",
                                             gene_col = "gene",
                                             padj_col = "p_val_adj",
                                             logfc_col = "avg_log2FC",
                                             padj_threshold = 1e-5,
                                             logfc_threshold = 2,
                                             all_genes = NULL) {
  validate_marker_table_for_go(markers, cluster_col, gene_col, padj_col, logfc_col)

  out <- markers %>%
    dplyr::filter(
      !is.na(.data[[padj_col]]),
      !is.na(.data[[logfc_col]]),
      .data[[padj_col]] < padj_threshold,
      .data[[logfc_col]] > logfc_threshold
    ) %>%
    dplyr::transmute(
      ClusterID = as.character(.data[[cluster_col]]),
      gene = as.character(.data[[gene_col]]),
      p_val_adj = .data[[padj_col]],
      avg_log2FC = .data[[logfc_col]]
    ) %>%
    dplyr::distinct(.data$ClusterID, .data$gene, .keep_all = TRUE)

  if (!is.null(all_genes)) {
    out <- out %>%
      dplyr::mutate(in_topgo_universe = .data$gene %in% all_genes)
  }

  out %>%
    dplyr::arrange(.data$ClusterID, .data$p_val_adj, dplyr::desc(.data$avg_log2FC))
}

run_topgo_for_gene_set <- function(selected_genes,
                                   all_genes,
                                   gene2GO,
                                   ontology = "BP",
                                   node_size = 3) {
  selected_genes <- unique(as.character(selected_genes))
  selected_genes <- selected_genes[selected_genes %in% all_genes]

  geneList <- ifelse(all_genes %in% selected_genes, 1L, 0L)
  names(geneList) <- all_genes

  GOdata <- methods::new(
    "topGOdata",
    ontology = ontology,
    allGenes = geneList,
    geneSelectionFun = function(x) x == 1,
    annot = topGO::annFUN.gene2GO,
    nodeSize = node_size,
    gene2GO = gene2GO
  )

  resultFisher <- topGO::runTest(
    GOdata,
    algorithm = "weight01",
    statistic = "fisher"
  )

  results <- topGO::GenTable(
    GOdata,
    Fisher = resultFisher,
    topNodes = length(GOdata@graph@nodes)
  )

  results$GO.ID <- as.character(results$GO.ID)
  results$Term <- as.character(results$Term)
  results$Fisher_raw <- as.character(results$Fisher)
  results$Fisher <- suppressWarnings(
    as.numeric(gsub("^<\\s*", "", results$Fisher_raw))
  )
  results$Fisher_FDR <- stats::p.adjust(results$Fisher, method = "BH")

  list(
    GOdata = GOdata,
    resultFisher = resultFisher,
    results = results
  )
}

run_topgo_enrichment_by_cluster <- function(markers,
                                            all_genes,
                                            gene2GO,
                                            cluster_col = "cluster",
                                            gene_col = "gene",
                                            padj_col = "p_val_adj",
                                            logfc_col = "avg_log2FC",
                                            padj_threshold = 1e-5,
                                            logfc_threshold = 2,
                                            ontology = "BP",
                                            node_size = 3) {
  validate_marker_table_for_go(markers, cluster_col, gene_col, padj_col, logfc_col)

  clusters <- sort(unique(as.character(markers[[cluster_col]])))

  results_list <- vector("list", length(clusters))
  names(results_list) <- clusters

  for (cluster_id in clusters) {
    message("Running topGO for cluster: ", cluster_id)

    selected_genes <- markers[[gene_col]][
      markers[[cluster_col]] == cluster_id &
        !is.na(markers[[padj_col]]) &
        !is.na(markers[[logfc_col]]) &
        markers[[padj_col]] < padj_threshold &
        markers[[logfc_col]] > logfc_threshold
    ]

    selected_genes <- unique(as.character(selected_genes))
    selected_genes <- selected_genes[selected_genes %in% all_genes]

    if (length(selected_genes) == 0) {
      warning("No selected genes in universe for cluster: ", cluster_id)
      results_list[[cluster_id]] <- data.frame(
        GO.ID = character(),
        Term = character(),
        Annotated = integer(),
        Significant = integer(),
        Expected = numeric(),
        Fisher = numeric(),
        Fisher_raw = character(),
        Fisher_FDR = numeric(),
        ClusterID = character(),
        n_selected_genes = integer(),
        stringsAsFactors = FALSE
      )
      next
    }

    res <- run_topgo_for_gene_set(
      selected_genes = selected_genes,
      all_genes = all_genes,
      gene2GO = gene2GO,
      ontology = ontology,
      node_size = node_size
    )

    tbl <- res$results
    tbl$ClusterID <- cluster_id
    tbl$n_selected_genes <- length(selected_genes)

    results_list[[cluster_id]] <- tbl
  }

  dplyr::bind_rows(results_list)
}

go_source_from_ontology <- function(ontology) {
  ontology <- toupper(ontology)

  if (ontology == "BP") return("GO:BP")
  if (ontology == "MF") return("GO:MF")
  if (ontology == "CC") return("GO:CC")

  paste0("GO:", ontology)
}

topgo_to_gostplot_table <- function(results,
                                    ontology = "BP",
                                    max_terms_per_cluster = 20) {
  if (nrow(results) == 0) {
    return(data.frame())
  }

  source <- go_source_from_ontology(ontology)

  results %>%
    dplyr::filter(!is.na(.data$Fisher_FDR)) %>%
    dplyr::group_by(.data$ClusterID) %>%
    dplyr::arrange(.data$Fisher_FDR, .by_group = TRUE) %>%
    dplyr::slice_head(n = max_terms_per_cluster) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      query = as.character(.data$ClusterID),
      source = source,
      term_id = as.character(.data$GO.ID),
      term_name = as.character(.data$Term),
      p_value = as.numeric(.data$Fisher_FDR),
      negative_log10_p_value = -log10(pmax(.data$p_value, .Machine$double.xmin)),
      term_size = suppressWarnings(as.integer(.data$Annotated)),
      intersection_size = suppressWarnings(as.integer(.data$Significant)),
      effective_domain_size = NA_integer_,
      precision = NA_real_,
      recall = NA_real_
    ) %>%
    dplyr::select(
      query,
      source,
      term_id,
      term_name,
      p_value,
      negative_log10_p_value,
      term_size,
      intersection_size,
      effective_domain_size,
      precision,
      recall,
      ClusterID,
      GO.ID,
      Term,
      Fisher,
      Fisher_FDR,
      Annotated,
      Significant,
      Expected,
      n_selected_genes
    )
}

plot_go_enrichment_dotplot <- function(gost_style_results,
                                       top_n_terms_per_cluster = 20,
                                       label_width = 70) {
  if (nrow(gost_style_results) == 0) {
    return(
      ggplot2::ggplot() +
        ggplot2::theme_void() +
        ggplot2::ggtitle("No enriched GO terms passed the selected threshold")
    )
  }

  plot_df <- gost_style_results %>%
    dplyr::group_by(.data$ClusterID) %>%
    dplyr::arrange(.data$p_value, .by_group = TRUE) %>%
    dplyr::slice_head(n = top_n_terms_per_cluster) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      term_label = paste0(.data$term_name, " (", .data$term_id, ")"),
      term_label = stringr::str_wrap(.data$term_label, width = label_width),
      ClusterID = as.factor(.data$ClusterID),
      source = as.factor(.data$source),
      term_label = forcats::fct_reorder(.data$term_label, .data$negative_log10_p_value)
    )

  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = .data$ClusterID,
      y = .data$term_label,
      size = .data$intersection_size,
      color = .data$negative_log10_p_value
    )
  ) +
    ggplot2::geom_point(alpha = 0.9) +
    ggplot2::facet_grid(
      rows = ggplot2::vars(.data$source),
      scales = "free_y",
      space = "free_y"
    ) +
    ggplot2::scale_size_continuous(name = "Significant genes") +
    ggplot2::scale_color_viridis_c(name = "-log10(FDR)") +
    ggplot2::labs(
      x = "Cluster / cell type",
      y = "GO term",
      title = "GO enrichment from FindAllMarkers genes"
    ) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
      axis.text.y = ggplot2::element_text(size = 8),
      axis.title = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      strip.text.y = ggplot2::element_text(face = "bold"),
      legend.position = "right"
    )
}

