library(Seurat)
library(scCustomize)
library(stringr)

#rm(list=ls())

## Functions

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


# select global features to integrate a list of seurat objects
select_global_features <- function(sample_list, nfeatures = 3000, assay = "RNA") {
  Seurat::SelectIntegrationFeatures(object.list = sample_list, nfeatures = nfeatures, assay = assay)
}

#R
run_branch <- function(obj, reduction, dims, prefix,
                       resolution = 0.5, seed = 1, algorithm = 4) {
  
  # Store BOTH nn and snn graphs with explicit names
  nn_name  <- paste0(prefix, "_nn")
  snn_name <- paste0(prefix, "_snn")
  
  # 1) Neighbors (creates nn + snn graphs)
  obj <- Seurat::FindNeighbors(
    object     = obj,
    reduction  = reduction,
    dims       = dims,
    graph.name = c(nn_name, snn_name),
    verbose    = FALSE
  )
  
  # 2) Clustering (one metadata column per resolution)
  if (length(resolution) == 1) {
    cluster_col <- paste0(prefix, "_clusters_res.", resolution)
    obj <- Seurat::FindClusters(
      object       = obj,
      graph.name   = snn_name,
      resolution   = resolution,
      algorithm    = algorithm,
      cluster.name = cluster_col,
      verbose      = FALSE
    )
  } else {
    for (r in resolution) {
      cluster_col <- paste0(prefix, "_clusters_res.", r)
      obj <- Seurat::FindClusters(
        object       = obj,
        graph.name   = snn_name,
        resolution   = r,
        algorithm    = algorithm,
        cluster.name = cluster_col,
        verbose      = FALSE
      )
    }
  }
  
  # 3) UMAP (only depends on reduction+dims, not on clustering resolution)
  umap_name <- paste0(prefix, "_umap")
  obj <- Seurat::RunUMAP(
    object         = obj,
    reduction      = reduction,
    dims           = dims,
    reduction.name = umap_name,
    seed.use       = seed,
    verbose        = FALSE
  )
  
  obj
}

#' Downsample Seurat to equal cells per cluster (v4/v5 safe)
#' @param srt        Seurat object
#' @param group_by   metadata column for groups (e.g., "seurat_clusters")
#' @param n          target cells per group; if NULL uses the minimum group size
#' @param seed       RNG seed
#' @param return_ids if TRUE return cell barcodes instead of a Seurat object
#' @return           Seurat object or character vector of barcodes
downsample_equal_per_cluster <- function(srt, group_by = "seurat_clusters",
                                         n = NULL, seed = 1L, return_ids = FALSE) {
  stopifnot(group_by %in% colnames(srt@meta.data))
  grp <- srt@meta.data[[group_by]]
  if (!is.factor(grp)) grp <- factor(grp)
  tab <- table(grp)
  if (length(tab) < 2L) stop("Need at least 2 groups to downsample.")
  
  target <- if (is.null(n)) min(tab) else as.integer(n)
  if (any(tab < target)) {
    warning(sprintf("Some groups have < %d cells; they will be kept at full size.", target))
  }
  
  set.seed(seed)
  cells_by_grp <- split(rownames(srt@meta.data), grp)
  keep <- unlist(lapply(names(cells_by_grp), function(k) {
    ids <- cells_by_grp[[k]]
    sample(ids, size = min(length(ids), target), replace = FALSE)
  }), use.names = FALSE)
  
  if (isTRUE(return_ids)) return(keep)
  
  # Version-agnostic subset
  if ("Subset" %in% getNamespaceExports("Seurat")) {
    Seurat::Subset(srt, cells = keep)        # Seurat v5
  } else {
    subset(srt, cells = keep)        # Seurat v4
  }
}


set.seed(42)

setwd("/data/Mlig_scRNA_Seq/Francisco/Mlig_scRNASeq_atlas/")

## Data to add module scores


#load list of curated marker genes
marker_genes <- read.table("~/Projects/Erin/M_lignano/Mlig_scRNASeq_atlas/data/metadata/literature_marker_genes.tsv", sep = "\t", header = TRUE, quote = "")

marker_genes$Gene_ID <- str_replace(marker_genes$Gene_ID, pattern = "_", replacement = "-")
anchor_cell_of_adhesive_organ <- marker_genes$Gene_ID[marker_genes$Cell_tissue_type == "anchor_cell_of_adhesive_organ"]
cement_gland_cell <- marker_genes$Gene_ID[marker_genes$Cell_tissue_type == "cement_gland_cell"]
epidermal_secretory_cell <- marker_genes$Gene_ID[marker_genes$Cell_tissue_type == "epidermal_secretory_cell"]
female_antrum <- marker_genes$Gene_ID[marker_genes$Cell_tissue_type == "female_antrum"]
female_germline_ovary <- marker_genes$Gene_ID[marker_genes$Cell_tissue_type == "female_germline_ovary"]
gut <- marker_genes$Gene_ID[marker_genes$Cell_tissue_type == "gut"]
gut_specialized_cell <- marker_genes$Gene_ID[marker_genes$Cell_tissue_type == "gut_specialized_cell"]
male_germline_testis <- marker_genes$Gene_ID[marker_genes$Cell_tissue_type == "male_germline_testis"]
muscle <- marker_genes$Gene_ID[marker_genes$Cell_tissue_type == "muscle"]
nervous_system <- marker_genes$Gene_ID[marker_genes$Cell_tissue_type == "nervous_system"]
prostate <- marker_genes$Gene_ID[marker_genes$Cell_tissue_type == "prostate"]
rhabdite <- marker_genes$Gene_ID[marker_genes$Cell_tissue_type == "rhabdite-containing_cell"]
secretory_cell_of_adhesive_organ <- marker_genes$Gene_ID[marker_genes$Cell_tissue_type == "secretory_cell_of_adhesive_organ"]
stem_cells_progenitors <- marker_genes$Gene_ID[marker_genes$Cell_tissue_type == "stem_cells_progenitors"]
epidermal <- marker_genes$Gene_ID[marker_genes$Cell_tissue_type == "epidermal"]

marker_genes_list <- list(anchor_cell_of_adhesive_organ,
                          cement_gland_cell,
                          epidermal_secretory_cell,
                          female_antrum,
                          female_germline_ovary,
                          gut,
                          gut_specialized_cell,
                          male_germline_testis,
                          muscle,
                          nervous_system,
                          prostate,
                          rhabdite,
                          secretory_cell_of_adhesive_organ,
                          stem_cells_progenitors,
                          epidermal)

names(marker_genes_list) <- c("anchor_cell_of_adhesive_organ",
                              "cement_gland_cell",
                              "epidermal_secretory",
                              "female_antrum",
                              "female_germline_ovary",
                              "gut",
                              "gut_specialized",
                              "male_germline_testis",
                              "muscle",
                              "nervous_system",
                              "prostate",
                              "rhabdite",
                              "secretory_cell_of_adhesive_organ",
                              "progenitor",
                              "epidermal")




# Import list of srt objects (filtered, SCTransformed)
sample_list <- readRDS(file = "/data/Mlig_scRNA_Seq/Francisco/Mlig_scRNASeq_atlas/results/objs/sct_list_v2_tmp.rds")

srt_obj <- readRDS(file = "/data/Mlig_scRNA_Seq/Francisco/Mlig_scRNASeq_atlas/results/objs/srt_obj_norm_scaled_RNA_plus_SCT_v2.rds")

#srt_obj <- add_module_scores(srt_obj, gene_lists = marker_genes_list)

sample_list_all <- sample_list

srt_obj_all <- srt_obj

rm(sample_list)
rm(srt_obj)

gc()

# #srt_obj_all$Stage <- factor(srt_obj_all$Stage, levels = c("0-1_days", "3_days", "7_days", "Adult"))
# 
# # srt_obj_all$sample_id <- factor(srt_obj_all$sample_id, levels = c("17_0_1d_Hatchlings",
# #                                                                   "18_0_1d_Hatchling",
# #                                                                   "19_0_1d_Hatchling",
# #                                                                   "11_unsorted_3dph",
# #                                                                   "12_unsorted_3dph",
# #                                                                   "16_3d_juveniles",
# #                                                                   "8_unsorted_Hatchling_7d",
# #                                                                   "9_unsorted_Hatchling_7d",
# #                                                                   "15_7d_juveniles",
# #                                                                   "2_Adult_Unsorted",
# #                                                                   "2_Adult_Unsorted_DASH",
# #                                                                   "7_Sorted_Adult"))
# # 

## Process each stage individually
# 
# # Adult
# sample_IDs_adult <- c("2_Adult_Unsorted", "2_Adult_Unsorted_DASH", "7_Sorted_Adult")
# sample_list_adult <- sample_list_all[names(sample_list_all) %in% sample_IDs_adult]
# srt_obj_adult <- subset(srt_obj_all, sample_id %in% sample_IDs_adult)
# 
# sample_list_adult <- lapply(sample_list_adult, function(obj) {
#   DefaultAssay(obj) <- "SCT"
#   obj
# })
# 
# ## 6000 features from adults

#variable_features_SCT_adult <- SelectIntegrationFeatures(sample_list_adult, nfeatures = 6000)

saveRDS(variable_features_SCT_adult, file = "/data/Mlig_scRNA_Seq/Francisco/Mlig_scRNASeq_atlas/results/objs/variable_features_SCT_adult.rds")

variable_features_SCT_adult <- readRDS(variable_features_SCT_adult, file = "/data/Mlig_scRNA_Seq/Francisco/Mlig_scRNASeq_atlas/results/objs/variable_features_SCT_adult.rds")

variable_features_SCT_all <- SelectIntegrationFeatures(sample_list_all, nfeatures = 2000)

saveRDS(variable_features_SCT_all, file = "/data/Mlig_scRNA_Seq/Francisco/Mlig_scRNASeq_atlas/results/objs/variable_features_SCT_all.rds")

sample_list_all <- PrepSCTIntegration(
  sample_list_all,
  anchor.features = variable_features_SCT_all
)

 
sample_list_all <- lapply(sample_list_all, function(x) {
  x <- RunPCA(x, assay = "SCT", features = variable_features_SCT_all)
})
 

gc()

integration_anchors <- FindIntegrationAnchors(object.list = sample_list_all,
                                              normalization.method = "SCT",
                                              anchor.features = variable_features_SCT_all,
                                              reduction = "rpca")

# # this command creates an 'integrated' data assay
srt_integrated <- IntegrateData(anchorset = integration_anchors)
# 
# # specify that we will perform downstream analysis on the corrected data note that the
# # original unmodified data still resides in the 'RNA' assay
DefaultAssay(srt_integrated) <- "integrated"

gc()

# # Run the standard workflow for visualization and clustering
srt_integrated <- ScaleData(srt_integrated, verbose = FALSE)
srt_integrated <- RunPCA(srt_integrated, npcs = 40, verbose = FALSE)
srt_integrated <- RunUMAP(srt_integrated, reduction = "pca", dims = 1:40)
srt_integrated <- FindNeighbors(srt_integrated, reduction = "pca", dims = 1:40)
srt_integrated <- FindClusters(srt_integrated, resolution = 0.5)
# 
# DefaultAssay(srt_integrated) <- "RNA"

#srt_integrated <- readRDS(file = "/data/Mlig_scRNA_Seq/Francisco/Mlig_scRNASeq_atlas/results/objs/srt_obj_norm_scaled_RNA_plus_SCT_v2_integrated.rds")

srt_integrated <- JoinLayers(srt_integrated, assay = "RNA", layers = c("counts", "data", "scale.data"))

srt_integrated <- add_module_scores(srt_integrated, gene_lists = marker_genes_list)

FeaturePlot_scCustom(srt_integrated, features = names(marker_genes_list))

saveRDS(srt_integrated, file = "/data/Mlig_scRNA_Seq/Francisco/Mlig_scRNASeq_atlas/results/objs/srt_obj_norm_scaled_RNA_plus_SCT_v2_integrated.rds")



integration_anchors <- FindIntegrationAnchors(object.list = sample_list_all,
                                              normalization.method = "SCT",
                                              anchor.features = variable_features_SCT_all,
                                              reduction = "cca")

gc()

saveRDS(integration_anchors, file = "/data/Mlig_scRNA_Seq/Francisco/Mlig_scRNASeq_atlas/results/objs/integration_anchors_cca.rds")


# # this command creates an 'integrated' data assay
srt_integrated <- IntegrateData(anchorset = integration_anchors)
# 
# # specify that we will perform downstream analysis on the corrected data note that the
# # original unmodified data still resides in the 'RNA' assay
DefaultAssay(srt_integrated) <- "integrated"

gc()

# # Run the standard workflow for visualization and clustering
srt_integrated <- ScaleData(srt_integrated, verbose = FALSE)
srt_integrated <- RunPCA(srt_integrated, npcs = 40, verbose = FALSE)
srt_integrated <- RunUMAP(srt_integrated, reduction = "pca", dims = 1:40)
srt_integrated <- FindNeighbors(srt_integrated, reduction = "pca", dims = 1:40)
srt_integrated <- FindClusters(srt_integrated, resolution = 0.5)


# sample_list_adult <- lapply(sample_list_adult, function(obj) {
#   DefaultAssay(obj) <- "RNA"
#   obj
# })

#sample_list_adult <- lapply(sample_list_adult, function(obj) {
#  obj <- FindVariableFeatures(obj, nfeatures = 10000)
#  obj
#})

variable_features_RNA_adult <- SelectIntegrationFeatures(sample_list_adult, nfeatures = 6000, assay = rep("RNA", length(sample_list_adult)))

srt_obj_adult <- RunPCA(srt_obj_adult,
                        assay = "SCT",
                        features = variable_features_SCT_adult,
                        reduction.name = "pca_sct",
                        reduction.key = "pca.sct_",
                        seed.use = 42,
                        verbose  = TRUE)

#srt_obj_adult <- split(srt_obj_adult, f = srt_obj_adult$sample_id)

ElbowPlot_scCustom(srt_obj_adult, reduction = "pca_sct")

# 20 and 40 dimensions

DefaultAssay(srt_obj_adult) <- "SCT"

srt_obj_adult <- IntegrateLayers(object = srt_obj_adult,
                                 method = RPCAIntegration,
                                 features = variable_features_SCT_adult,
                                 assay = "SCT",                      # explicit, avoids picking RNA
                                 normalization.method = "SCT",
                                 orig.reduction = "pca_sct",
                                 new.reduction = "integrated.sct.rpca",
                                 verbose = TRUE)

srt_obj_adult <- run_branch(srt_obj_adult,
                            reduction = "integrated.sct.rpca",
                            dims = 1:40,
                            prefix = "integrated.sct.rpca_40",
                            resolution = c(0.6, 1.2, 2, 4)
)

srt_obj_adult <- run_branch(srt_obj_adult,
                            reduction = "integrated.sct.rpca",
                            dims = 1:20,
                            prefix = "integrated.sct.rpca_20",
                            resolution = c(0.6, 1.2, 2, 4)
)

gc()

FeaturePlot_scCustom(srt_obj_adult, features = names(marker_genes_list), aspect_ratio = 1, reduction = "integrated.sct.rpca_40_umap")

FeaturePlot_scCustom(srt_obj_adult, features = names(marker_genes_list), aspect_ratio = 1, reduction = "integrated.sct.rpca_20_umap")

saveRDS(srt_obj_adult, file = "/data/Mlig_scRNA_Seq/Francisco/Mlig_scRNASeq_atlas/results/objs/srt_obj_adult.rds")


## 6000 features from adults only, all stages

srt_obj_all <- RunPCA(srt_obj_all,
                  assay = "SCT",
                  features = variable_features_SCT_adult,
                  reduction.name = "pca_sct",
                  reduction.key = "pca.sct_",
                  seed.use = 42,
                  verbose  = TRUE)


ElbowPlot_scCustom(srt_obj_all, reduction = "pca_sct")

# 15 and 40 dimensions

srt_obj_all <- IntegrateLayers(object = srt_obj_all,
                           method = RPCAIntegration,
                           features = variable_features_SCT_adult,
                           assay = "SCT",                      # explicit, avoids picking RNA
                           normalization.method = "SCT",
                           orig.reduction = "pca_sct",
                           new.reduction = "integrated.sct.rpca",
                           verbose = TRUE)



# RPCA, 15 PCs
srt_obj_all <- run_branch(srt_obj_all,
                          reduction = "integrated.sct.rpca",
                          dims = 1:15,
                          prefix = "integrated.sct.rpca_15",
                          resolution = c(0.6, 1.2, 2, 4)
)

gc()

# RPCA, 40 PCs
srt_obj_all <- run_branch(srt_obj_all,
                          reduction = "integrated.sct.rpca",
                          dims = 1:40,
                          prefix = "integrated.sct.rpca_40",
                          resolution = c(0.6, 1.2, 2, 4)
)

gc()


FeaturePlot_scCustom(srt_obj_all, features = names(marker_genes_list), aspect_ratio = 1, reduction = "integrated.sct.rpca_40_umap")

FeaturePlot_scCustom(srt_obj_all, features = names(marker_genes_list), aspect_ratio = 1, reduction = "integrated.sct.rpca_20_umap")






## 2000 features from adults

# Adult
sample_IDs_adult <- c("2_Adult_Unsorted", "2_Adult_Unsorted_DASH", "7_Sorted_Adult")
sample_list_adult <- sample_list_all[names(sample_list_all) %in% sample_IDs_adult]
srt_obj_adult <- subset(srt_obj, sample_id %in% sample_IDs_adult)

DefaultAssay(srt_obj_adult) <- "SCT"

variable_features_SCT_adult <- SelectIntegrationFeatures(sample_list_adult, nfeatures = 2000)

srt_obj_adult <- RunPCA(srt_obj_adult,
                        assay = "SCT",
                        features = variable_features_SCT_adult,
                        reduction.name = "pca_sct",
                        reduction.key = "pca.sct_",
                        seed.use = 42,
                        verbose  = TRUE)

#srt_obj_adult <- split(srt_obj_adult, f = srt_obj_adult$sample_id)

ElbowPlot_scCustom(srt_obj_adult, reduction = "pca_sct")

# 20 and 40 dimensions

DefaultAssay(srt_obj_adult) <- "SCT"

srt_obj_adult <- IntegrateLayers(object = srt_obj_adult,
                                 method = RPCAIntegration,
                                 features = variable_features_SCT_adult,
                                 assay = "SCT",                      # explicit, avoids picking RNA
                                 normalization.method = "SCT",
                                 orig.reduction = "pca_sct",
                                 new.reduction = "integrated.sct.rpca",
                                 verbose = TRUE)

srt_obj_adult <- run_branch(srt_obj_adult,
                            reduction = "integrated.sct.rpca",
                            dims = 1:40,
                            prefix = "integrated.sct.rpca_40",
                            resolution = c(0.6, 1.2, 2, 4)
)

srt_obj_adult <- run_branch(srt_obj_adult,
                            reduction = "integrated.sct.rpca",
                            dims = 1:20,
                            prefix = "integrated.sct.rpca_20",
                            resolution = c(0.6, 1.2, 2, 4)
)

gc()

FeaturePlot_scCustom(srt_obj_adult, features = names(marker_genes_list), aspect_ratio = 1, reduction = "integrated.sct.rpca_40_umap")

FeaturePlot_scCustom(srt_obj_adult, features = names(marker_genes_list), aspect_ratio = 1, reduction = "integrated.sct.rpca_20_umap")

saveRDS(srt_obj_adult, file = "/data/Mlig_scRNA_Seq/Francisco/Mlig_scRNASeq_atlas/results/objs/srt_obj_adult.rds")


## 2000 features from adults only, all stages

srt_obj_all <- srt_obj

srt_obj_all <- RunPCA(srt_obj_all,
                      assay = "SCT",
                      features = variable_features_SCT_adult,
                      reduction.name = "pca_sct",
                      reduction.key = "pca.sct_",
                      seed.use = 42,
                      verbose  = TRUE)

ElbowPlot_scCustom(srt_obj_all, reduction = "pca_sct")

# 20 and 40 dimensions

srt_obj_all <- IntegrateLayers(object = srt_obj_all,
                               method = RPCAIntegration,
                               features = variable_features_SCT_adult,
                               assay = "SCT",                      # explicit, avoids picking RNA
                               normalization.method = "SCT",
                               orig.reduction = "pca_sct",
                               new.reduction = "integrated.sct.rpca",
                               verbose = TRUE)

# RPCA, 20 PCs
srt_obj_all <- run_branch(srt_obj_all,
                          reduction = "integrated.sct.rpca",
                          dims = 1:20,
                          prefix = "integrated.sct.rpca_20",
                          resolution = c(0.6, 1.2, 2, 4)
)

gc()

# RPCA, 40 PCs
srt_obj_all <- run_branch(srt_obj_all,
                          reduction = "integrated.sct.rpca",
                          dims = 1:40,
                          prefix = "integrated.sct.rpca_40",
                          resolution = c(0.6, 1.2, 2, 4)
)

gc()


FeaturePlot_scCustom(srt_obj_all, features = names(marker_genes_list), aspect_ratio = 1, reduction = "integrated.sct.rpca_40_umap")

FeaturePlot_scCustom(srt_obj_all, features = names(marker_genes_list), aspect_ratio = 1, reduction = "integrated.sct.rpca_20_umap")





# 
# srt_obj2 <- RunUMAP(srt_obj, dims = 1:40, reduction = "integrated.rna.rpca", n.neighbors = 50L, min.dist = 0.01)
# srt_obj3 <- RunUMAP(srt_obj, dims = 1:40, reduction = "integrated.rna.rpca", n.neighbors = 5L, min.dist = 0.01)
# srt_obj4 <- RunUMAP(srt_obj, dims = 1:40, reduction = "integrated.sct.rpca", n.neighbors = 50L, min.dist = 0.01)
# srt_obj5 <- RunUMAP(srt_obj, dims = 1:40, reduction = "integrated.sct.rpca", n.neighbors = 5L, min.dist = 0.01)
# 
# srt_obj2 <- RunUMAP(srt_obj, dims = 1:15, reduction = "integrated.rna.rpca", n.neighbors = 50L, min.dist = 0.01)
# srt_obj3 <- RunUMAP(srt_obj, dims = 1:15, reduction = "integrated.rna.rpca", n.neighbors = 5L, min.dist = 0.01)
# srt_obj4 <- RunUMAP(srt_obj, dims = 1:15, reduction = "integrated.sct.rpca", n.neighbors = 50L, min.dist = 0.01)
# srt_obj5 <- RunUMAP(srt_obj, dims = 1:15, reduction = "integrated.sct.rpca", n.neighbors = 5L, min.dist = 0.01)

# srt_obj5 <- RunUMAP(srt_obj, dims = 1:40, reduction = "integrated.sct.rpca", n.neighbors = 5L, min.dist = 0.01)
# srt_obj6 <- RunUMAP(srt_obj, dims = 1:40, reduction = "integrated.sct.harmony", n.neighbors = 5L, min.dist = 0.01)


gc()

saveRDS(srt_obj, file = "/data/Mlig_scRNA_Seq/Francisco/Mlig_scRNASeq_atlas/results/objs/srt_obj_norm_scaled_RNA_plus_SCT_plus_integration.rds")

srt_obj <- readRDS(file = "/data/Mlig_scRNA_Seq/Francisco/Mlig_scRNASeq_atlas/results/objs/srt_obj_norm_scaled_RNA_plus_SCT_plus_integration.rds")


## Epidermal markers

# 0-1 days
sample_IDs_0_1_days <- c("17_0_1d_Hatchlings", "18_0_1d_Hatchling", "19_0_1d_Hatchling")
sample_list_0_1_days <- sample_list_all[names(sample_list_all) %in% sample_IDs_0_1_days]
srt_obj_0_1_days <- subset(srt_obj_all, sample_id %in% sample_IDs_0_1_days)
variable_features_SCT_0_1_days <- SelectIntegrationFeatures(sample_list_0_1_days, nfeatures = 2000, assay = rep("SCT", length(sample_list_0_1_days)))
variable_features_RNA_0_1_days <- SelectIntegrationFeatures(sample_list_0_1_days, nfeatures = 2000, assay = rep("RNA", length(sample_list_0_1_days)))

# 3 days
sample_IDs_3_days <- c("11_unsorted_3dph", "12_unsorted_3dph", "16_3d_juveniles")
sample_list_3_days <- sample_list_all[names(sample_list_all) %in% sample_IDs_3_days]
srt_obj_3_days <- subset(srt_obj_all, sample_id %in% sample_IDs_3_days)
variable_features_SCT_3_days <- SelectIntegrationFeatures(sample_list_3_days, nfeatures = 2000, assay = rep("SCT", length(sample_list_3_days)))
variable_features_RNA_3_days <- SelectIntegrationFeatures(sample_list_3_days, nfeatures = 2000, assay = rep("RNA", length(sample_list_3_days)))

# 7 days
sample_IDs_7_days <- c("8_unsorted_Hatchling_7d", "9_unsorted_Hatchling_7d", "15_7d_juveniles")
sample_list_7_days <- sample_list_all[names(sample_list_all) %in% sample_IDs_7_days]
srt_obj_7_days <- subset(srt_obj_all, sample_id %in% sample_IDs_7_days)
variable_features_SCT_7_days <- SelectIntegrationFeatures(sample_list_7_days, nfeatures = 2000, assay = rep("SCT", length(sample_list_7_days)))
variable_features_RNA_7_days <- SelectIntegrationFeatures(sample_list_7_days, nfeatures = 2000, assay = rep("RNA", length(sample_list_7_days)))



### Upset plots of features for the SCT and RNA assays
library(ComplexHeatmap)

lt_SCT = list("0_1_days" = variable_features_SCT_0_1_days,
              "3_days"  = variable_features_SCT_3_days,
              "7_days"  = variable_features_SCT_7_days,
              "Adult" = variable_features_SCT_adult
)

#mt_SCT <- list_to_matrix(lt_SCT)
mt_SCT = make_comb_mat(lt_SCT)
UpSet(mt_SCT)


lt_RNA = list("0_1_days" = variable_features_RNA_0_1_days,
              "3_days"  = variable_features_RNA_3_days,
              "7_days"  = variable_features_RNA_7_days,
              "Adult" = variable_features_RNA_adult
)

#mt_SCT <- list_to_matrix(lt_SCT)
mt_RNA = make_comb_mat(lt_RNA)
UpSet(mt_RNA)


## Both sets combined
lt_all = list("0_1_days_SCT" = variable_features_SCT_0_1_days,
              "3_days_SCT"  = variable_features_SCT_3_days,
              "7_days_SCT"  = variable_features_SCT_7_days,
              "Adult_SCT" = variable_features_SCT_adult,
              "0_1_days_RNA" = variable_features_RNA_0_1_days,
              "3_days_RNA"  = variable_features_RNA_3_days,
              "7_days_RNA"  = variable_features_RNA_7_days,
              "Adult_RNA" = variable_features_RNA_adult
)

## Both sets combined
lt_all = list("0_1_days_SCT" = variable_features_SCT_0_1_days,
              "0_1_days_RNA" = variable_features_RNA_0_1_days,
              "3_days_SCT"  = variable_features_SCT_3_days,
              "3_days_RNA"  = variable_features_RNA_3_days,
              "7_days_SCT"  = variable_features_SCT_7_days,
              "7_days_RNA"  = variable_features_RNA_7_days,
              "Adult_SCT" = variable_features_SCT_adult,
              "Adult_RNA" = variable_features_RNA_adult
              
)


#mt_SCT <- list_to_matrix(lt_SCT)
mt_all = make_comb_mat(lt_all)
UpSet(mt_all)


## Integrate using unique features selected using the stage-aware strategy

# SCT

variable_features_SCT_stage_aware <- unique(c(variable_features_SCT_0_1_days,
                                              variable_features_SCT_3_days,
                                              variable_features_SCT_7_days,
                                              variable_features_SCT_adult))


DefaultAssay(srt_obj) <- "RNA"

