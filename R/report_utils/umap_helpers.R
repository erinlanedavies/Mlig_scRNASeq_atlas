## Themes / plot modifiers

theme_no_axes <- ggplot2::theme(
  axis.title   = ggplot2::element_blank(),
  axis.text    = ggplot2::element_blank(),
  axis.ticks   = ggplot2::element_blank(),
  axis.line    = ggplot2::element_blank(),
  panel.grid   = ggplot2::element_blank(),
  panel.border = ggplot2::element_blank()
)

clean_umap_theme <- function(legend_pos = "none",
                             title = FALSE,
                             title_size = 16) {
  base <- ggplot2::theme_void() +
    ggplot2::theme(
      legend.position = legend_pos,
      legend.title    = ggplot2::element_blank(),
      plot.margin     = ggplot2::margin(10, 10, 10, 10)
    )
  
  if (isTRUE(title)) {
    base + ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = title_size, face = "bold"),
      plot.title.position = "plot"
    )
  } else {
    # Hard-disable titles (prevents accidental blank title space)
    base + ggplot2::theme(
      plot.title = ggplot2::element_blank(),
      plot.title.position = "plot"
    )
  }
}


add_umap_compass <- function(p,
                             anchor_frac = c(0.05, 0.05),
                             arrow_frac  = 0.096,
                             arrow_lwd   = 0.6,
                             text_size   = 3.0,
                             x_label = "UMAP 1",
                             y_label = "UMAP 2",
                             hjust_x = -0.02) {
  
  b  <- ggplot2::ggplot_build(p)
  pp <- b$layout$panel_params[[1]]
  xr <- pp$x.range
  yr <- pp$y.range
  
  arrow_len <- arrow_frac * min(diff(xr), diff(yr))
  x0 <- xr[1] + anchor_frac[1] * diff(xr)
  y0 <- yr[1] + anchor_frac[2] * diff(yr)
  
  p +
    ggplot2::annotate(
      "segment",
      x = x0, xend = x0 + arrow_len, y = y0, yend = y0,
      linewidth = arrow_lwd,
      arrow = ggplot2::arrow(length = grid::unit(2.0, "mm"), type = "closed")
    ) +
    ggplot2::annotate(
      "segment",
      x = x0, xend = x0, y = y0, yend = y0 + arrow_len,
      linewidth = arrow_lwd,
      arrow = ggplot2::arrow(length = grid::unit(2.0, "mm"), type = "closed")
    ) +
    ggplot2::annotate(
      "text",
      x = x0 + arrow_len, y = y0,
      label = x_label,
      hjust = hjust_x, vjust = 0.5,
      size = text_size
    ) +
    ggplot2::annotate(
      "text",
      x = x0, y = y0 + arrow_len,
      label = y_label,
      hjust = 0.5, vjust = -0.25,
      size = text_size
    )
}



add_halo_labels <- function(obj, p, reduction, cluster_col = NULL, fun = median, size = 3) {
  
  umap_mat <- Seurat::Embeddings(obj, reduction)[, 1:2, drop = FALSE]
  
  clust <- if (is.null(cluster_col)) {
    as.character(Seurat::Idents(obj))
  } else {
    as.character(obj@meta.data[[cluster_col]])
  }
  
  umap_df <- data.frame(
    x = umap_mat[, 1],
    y = umap_mat[, 2],
    cluster = clust
  )
  
  centers <- stats::aggregate(cbind(x, y) ~ cluster, data = umap_df, FUN = fun)
  
  p +
    shadowtext::geom_shadowtext(
      data = centers,
      ggplot2::aes(x = x, y = y, label = cluster),
      colour = "black",
      bg.colour = "white",
      bg.r = 0.25,
      size = size,
      fontface = "bold"
    )
}

add_halo_text_colored <- function(obj, p, reduction,
                                  cluster_col,
                                  color_map,
                                  label_nudges = NULL,
                                  fun = median,
                                  size = 5,
                                  bg_r = 0.16) {
  
  umap_mat <- Seurat::Embeddings(obj, reduction)[, 1:2, drop = FALSE]
  clust <- as.character(obj@meta.data[[cluster_col]])
  
  df <- data.frame(
    x = umap_mat[, 1],
    y = umap_mat[, 2],
    cluster = clust
  )
  
  centers <- stats::aggregate(cbind(x, y) ~ cluster, data = df, FUN = fun)
  
  if (!is.null(label_nudges)) {
    centers <- dplyr::left_join(centers, label_nudges, by = "cluster")
    centers$nudge_x[is.na(centers$nudge_x)] <- 0
    centers$nudge_y[is.na(centers$nudge_y)] <- 0
  } else {
    centers$nudge_x <- 0
    centers$nudge_y <- 0
  }
  
  centers$x <- centers$x + centers$nudge_x
  centers$y <- centers$y + centers$nudge_y
  
  centers$lab_color <- unname(color_map[as.character(centers$cluster)])
  centers$lab_color[is.na(centers$lab_color)] <- "black"
  
  p +
    shadowtext::geom_shadowtext(
      data = centers,
      ggplot2::aes(x = x, y = y, label = cluster),
      inherit.aes = FALSE,
      color = centers$lab_color,
      bg.colour = "white",
      bg.r = bg_r,
      size = size,
      fontface = "bold"
    )
}


add_halo_text_colored_nudged <- function(
  obj,
  p,
  reduction,
  cluster_col,
  color_map,
  label_nudges = NULL,
  fun = median,
  size = 5,
  bg_r = 0.18
) {
  # UMAP coordinates
  um <- Seurat::Embeddings(obj, reduction)[, 1:2, drop = FALSE]
  cl <- as.character(obj@meta.data[[cluster_col]])

  df <- data.frame(
    x = um[, 1],
    y = um[, 2],
    cluster = cl
  )

  centers <- aggregate(cbind(x, y) ~ cluster, data = df, FUN = fun)

  # Join nudges (default 0,0)
  if (!is.null(label_nudges)) {
    centers <- dplyr::left_join(centers, label_nudges, by = "cluster")
    centers$nudge_x[is.na(centers$nudge_x)] <- 0
    centers$nudge_y[is.na(centers$nudge_y)] <- 0
  } else {
    centers$nudge_x <- 0
    centers$nudge_y <- 0
  }

  centers$x <- centers$x + centers$nudge_x
  centers$y <- centers$y + centers$nudge_y

  centers$lab_color <- unname(color_map[centers$cluster])
  centers$lab_color[is.na(centers$lab_color)] <- "black"

  p +
    shadowtext::geom_shadowtext(
      data = centers,
      aes(x = x, y = y, label = cluster),
      inherit.aes = FALSE,
      colour = centers$lab_color,
      bg.colour = "white",
      bg.r = bg_r,
      size = size,
      fontface = "bold"
    )
}

