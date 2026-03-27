theme_dotplot_atlas <- function(n_clusters,
                                cluster_text_size = 14,   # ⬆ bigger cluster labels
                                module_text_size  = 12,
                                module_lineheight = 0.92,
                                cluster_angle      = 45,
                                cluster_spacing_lines = NULL,
                                module_spacing_lines  = 0.20,
                                legend_key_cm = 0.35) {
  
  if (is.null(cluster_spacing_lines)) {
    cluster_spacing_lines <- dotplot_cluster_spacing(n_clusters)
  }
  
  ggplot2::theme(
    # no title (as per your rule)
    plot.title = ggplot2::element_blank(),
    
    # CLUSTER LABELS (x-axis after coord_flip)
    axis.text.x = ggplot2::element_text(
      size  = cluster_text_size,
      angle = cluster_angle,
      hjust = 1,
      vjust = 1,
      face  = "bold"
    ),
    
    # MODULE SCORE LABELS (y-axis)
    axis.text.y = ggplot2::element_text(
      size       = module_text_size,
      lineheight = module_lineheight,
      face       = "bold"
    ),
    
    panel.spacing.x = grid::unit(cluster_spacing_lines, "lines"),
    panel.spacing.y = grid::unit(module_spacing_lines,  "lines"),
    
    panel.grid.major.x = ggplot2::element_line(
      linewidth = 0.25,
      colour = "grey90"
    ),
    panel.grid.major.y = ggplot2::element_blank(),
    panel.grid.minor   = ggplot2::element_blank(),
    
    legend.key.height = grid::unit(legend_key_cm, "cm"),
    legend.key.width  = grid::unit(legend_key_cm, "cm"),
    
    plot.margin = grid::unit(c(0.8, 1.0, 0.8, 0.8), "cm")
  )
}

dotplot_cluster_spacing <- function(n_clusters,
                                    min_lines = 0.25,
                                    max_lines = 2.8) {
  x <- log10(max(n_clusters, 2))
  sp <- 0.55 + 1.35 * x
  sp <- max(min(sp, max_lines), min_lines)
  sp
}
