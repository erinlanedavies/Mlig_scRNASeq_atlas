split_seurat <- function(srt_obj, split_by) {
  if (!split_by %in% colnames(srt_obj@meta.data)) {
    stop("split_by column not found in Seurat meta.data: ", split_by)
  }
  Seurat::SplitObject(srt_obj, split.by = split_by)
}

sct_transform_list <- function(obj_list,
                               method = "glmGamPoi",
                               vars_to_regress = c("percent.mito", "S.Score", "G2M.Score"),
                               return_only_var_genes = FALSE,
                               verbose = TRUE,
                               ...) {

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
      verbose = FALSE,
      ...
    )
  }
  out
}
