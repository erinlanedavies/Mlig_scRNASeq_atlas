
# qc_plot_helpers.R
# Helper functions for QC plotting and saving figures

library(Seurat)
library(ggplot2)
library(scales)

qc_theme <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      strip.background = element_rect(fill = "grey95", color = "black"),
      strip.text = element_text(face = "bold"),
      legend.position = "right"
    )
}

make_qc_violin_plot <- function(srt_obj, features, sample_id_col = NULL, pt.size = 0) {

	  plot_list <- VlnPlot(
			           object = srt_obj,
				       features = features,
				       group.by = sample_id_col,
				           pt.size = pt.size,
				           combine = FALSE
					     )

  plot_list <- lapply(plot_list, function(p) {
			          p +
					        qc_theme() +
						      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
					        })

    patchwork::wrap_plots(plot_list, ncol = length(features))
}

make_feature_scatter <- function(srt_obj,
                                 feature1,
                                 feature2,
                                 title = NULL,
                                 pt.size = 0.05,
                                 alpha = 0.25) {

  df <- FetchData(srt_obj, vars = c(feature1, feature2))
  colnames(df) <- c("x", "y")

  ggplot(df, aes(x = x, y = y)) +
    geom_point(size = pt.size, alpha = alpha) +
    scale_x_continuous(labels = comma) +
    scale_y_continuous(labels = comma) +
    labs(title = title, x = feature1, y = feature2) +
    qc_theme()
}

save_qc_plot <- function(plot_obj, filename, out_dir, width = 11, height = 4.5) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(
    filename = file.path(out_dir, filename),
    plot = plot_obj,
    width = width,
    height = height,
    units = "in",
    dpi = 300,
    bg = "white"
  )
}

