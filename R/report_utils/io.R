# R/report_utils/io.R

# Resolve a full path given base_dir (if any) and a path from the table
resolve_path <- function(base_dir, path_from_table) {
  if (!is.null(base_dir) && nzchar(base_dir)) {
    return(file.path(base_dir, path_from_table))
  } else {
    return(path_from_table)
  }
}

# Load a count matrix given a path and an input type
# - For CellBender: expects an h5 file readable by scCustomize::Read_CellBender_h5_Mat
# - For CellRanger: accepts either a folder (Read10X) or h5 file (Read10X_h5)
load_counts <- function(path, input_type) {
  if (!file.exists(path) && !dir.exists(path)) {
    stop("Counts path does not exist: ", path)
  }

  if (input_type == "CellBender") {
    if (!grepl("\\.h5$", path, ignore.case = TRUE)) {
      stop("For CellBender input_type, path should point to a .h5 file: ", path)
    }
    mat <- scCustomize::Read_CellBender_h5_Mat(file_name = path)

  } else if (input_type == "CellRanger") {
    if (dir.exists(path)) {
      mat <- Seurat::Read10X(data.dir = path)
    } else if (grepl("\\.h5$", path, ignore.case = TRUE)) {
      mat <- Seurat::Read10X_h5(filename = path)
    } else {
      stop("For CellRanger input_type, path should be a folder or a .h5 file: ", path)
    }

  } else {
    stop("Unknown input_type: ", input_type, ". Use 'CellBender' or 'CellRanger'.")
  }

  return(mat)
}


# read seurat RDS
read_seurat_rds <- function(path) {
  if (!file.exists(path)) stop("Seurat RDS not found: ", path)
  readRDS(path)
}


# save RDS files
save_rds <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(object, file = path)
  path
}


#read table with cell cycle genes
read_cell_cycle_table <- function(path) {
  if (!file.exists(path)) stop("Cell cycle gene table not found: ", path)

  df <- read.table(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

  # Filter table as in your Rmd
  df <- df[df$Mlig_geneID != "N.A." &
             df$Obs != "low-qual" &
             df$Obs != "low-qual?", , drop = FALSE]

  # Replace underscores with dashes (base R)
  df$Mlig_geneID <- gsub("_", "-", df$Mlig_geneID)

  df
}

#create output directory, if needed, while preventing overwriting
make_out_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}


# read cell cycle table

get_cell_cycle_features <- function(cell_cycle_df,
                                    gene_col = "Mlig_geneID",
                                    state_col = "CellCycleState") {

  if (!all(c(gene_col, state_col) %in% colnames(cell_cycle_df))) {
    stop("cell_cycle_df must contain columns: ", gene_col, " and ", state_col)
  }

  s_feat   <- cell_cycle_df[[gene_col]][cell_cycle_df[[state_col]] == "S"]
  g2m_feat <- cell_cycle_df[[gene_col]][cell_cycle_df[[state_col]] == "G2"]

  list(
    s_features = unique(stats::na.omit(s_feat)),
    g2m_features = unique(stats::na.omit(g2m_feat))
  )
}

