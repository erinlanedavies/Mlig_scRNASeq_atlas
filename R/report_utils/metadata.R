# R/report_utils/metadata.R

# Add sample-level metadata to a Seurat object.
# metadata_row is a one-row data.frame (or tibble) with columns specified in metadata_cols.
add_sample_metadata <- function(obj, metadata_row, sample_id_col, metadata_cols) {
  sample_id <- metadata_row[[sample_id_col]]
  obj$sample_id <- as.character(sample_id)

  for (col in metadata_cols) {
    if (!col %in% colnames(metadata_row)) {
      warning("Metadata column '", col, "' not found in sample table; skipping.")
      next
    }
    value <- metadata_row[[col]]
    obj[[col]] <- as.character(value)
  }

  return(obj)
}

