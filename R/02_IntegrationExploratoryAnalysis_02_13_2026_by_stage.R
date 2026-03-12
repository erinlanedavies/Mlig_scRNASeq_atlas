library(Seurat)
library(scCustomize)

rm(list=ls())

## Functions

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

set.seed(42)

setwd("/data/Mlig_scRNA_Seq/Francisco/Mlig_scRNASeq_atlas/")

# Import list of srt objects (filtered, SCTransformed)
sample_list <- readRDS(file = "/data/Mlig_scRNA_Seq/Francisco/Mlig_scRNASeq_atlas/results/objs/sct_list_tmp.rds")

srt_obj <- readRDS(file = "/data/Mlig_scRNA_Seq/Francisco/Mlig_scRNASeq_atlas/results/objs/srt_obj_norm_scaled_RNA_plus_SCT.rds")

sample_list_all <- sample_list

srt_obj_all <- srt_obj

srt_obj_all$Stage <- factor(srt_obj_all$Stage, levels = c("0-1_days", "3_days", "7_days", "Adult"))

srt_obj_all$sample_id <- factor(srt_obj_all$sample_id, levels = c("17_0_1d_Hatchlings",
                                                                  "18_0_1d_Hatchling",
                                                                  "19_0_1d_Hatchling",
                                                                  "11_unsorted_3dph",
                                                                  "12_unsorted_3dph",
                                                                  "16_3d_juveniles",
                                                                  "8_unsorted_Hatchling_7d",
                                                                  "9_unsorted_Hatchling_7d",
                                                                  "15_7d_juveniles",
                                                                  "2_Adult_Unsorted",
                                                                  "2_Adult_Unsorted_DASH",
                                                                  "7_Sorted_Adult"))


## Process each stage individually

# Adult
sample_IDs_adult <- c("2_Adult_Unsorted", "2_Adult_Unsorted_DASH", "7_Sorted_Adult")
sample_list_adult <- sample_list_all[names(sample_list_all) %in% sample_IDs_adult]
srt_obj_adult <- subset(srt_obj_all, sample_id %in% sample_IDs_adult)
variable_features_SCT_adult <- SelectIntegrationFeatures(sample_list_adult, nfeatures = 2000, assay = rep("SCT", length(sample_list_adult)))
variable_features_RNA_adult <- SelectIntegrationFeatures(sample_list_adult, nfeatures = 2000, assay = rep("RNA", length(sample_list_adult)))

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

variable_features_SCT_stage_aware <- 




## SCT analysis
## Even though we'll use the SCT data, the split needs to be done in the RNA assay, see
## https://github.com/satijalab/seurat/issues/8361 
DefaultAssay(srt_obj) <- "RNA"
variable_features_SCT <- select_global_features(sample_list, nfeatures = 3000, assay = rep("SCT", length(sample_list)))

# Join layers to compute global PCA values
srt_obj <- JoinLayers(srt_obj)

srt_obj <- RunPCA(srt_obj,
                  assay = "SCT",
                  features = variable_features_SCT,
                  reduction.name = "pca_sct",
                  reduction.key = "pca.sct_",
                  seed.use = 42,
                  verbose  = TRUE)

# Split by sample ID to integrate
srt_obj <- split(srt_obj, f = srt_obj$sample_id)

DefaultAssay(srt_obj) <- "SCT"


srt_obj <- IntegrateLayers(object = srt_obj,
                           method = RPCAIntegration,
                           features = variable_features_SCT,
                           assay = "SCT",                      # explicit, avoids picking RNA
                           normalization.method = "SCT",
                           orig.reduction = "pca_sct",
                           new.reduction = "integrated.sct.rpca",
                           verbose = TRUE)

srt_obj <- IntegrateLayers(object = srt_obj,
                           method = HarmonyIntegration,
                           assay = "SCT",                      # explicit, avoids picking RNA
                           features = variable_features_SCT,
                           normalization.method = "SCT",
                           orig.reduction = "pca_sct",
                           new.reduction = "integrated.sct.harmony",
                           verbose = TRUE)


## RNA analysis
DefaultAssay(srt_obj) <- "RNA"
srt_obj <- JoinLayers(srt_obj)

variable_features_RNA <- VariableFeatures(FindVariableFeatures(srt_obj, assay = "RNA", nfeatures = 3000))

srt_obj <- RunPCA(srt_obj,
                  assay = "RNA",
                  features = variable_features_RNA,
                  reduction.name = "pca_rna",
                  reduction.key = "pca.rna_",
                  seed.use = 42,
                  verbose  = TRUE)

# Split for integration
srt_obj <- split(srt_obj, f = srt_obj$sample_id)

srt_obj <- IntegrateLayers(object = srt_obj,
                           assay = "RNA",
#                           layers = c("counts", "data", "scale.data"),
                           method = RPCAIntegration,
                           features = variable_features_RNA,
                           normalization.method = "LogNormalize",
                           orig.reduction = "pca_rna",
                           new.reduction = "integrated.rna.rpca",
                           verbose = TRUE)

srt_obj <- IntegrateLayers(object = srt_obj,
                           method = HarmonyIntegration,
                           features = variable_features_SCT,
                           normalization.method = "LogNormalize",
                           orig.reduction = "pca_rna",
                           new.reduction = "integrated.rna.harmony",
                           verbose = TRUE)

srt_obj <- JoinLayers(srt_obj, assay = "RNA", layers = c("counts", "data", "scale.data"))



## SCT downstream processing

# RPCA, 15 PCs
srt_obj <- run_branch(srt_obj,
                      reduction = "integrated.sct.rpca",
                      dims = 1:15,
                      prefix = "integrated.sct.rpca_15",
                      resolution = c(0.6, 1.2, 2, 4)
)

gc()

# RPCA, 40 PCs
srt_obj <- run_branch(srt_obj,
                      reduction = "integrated.sct.rpca",
                      dims = 1:40,
                      prefix = "integrated.sct.rpca_40",
                      resolution = c(0.6, 1.2, 2, 4)
)

gc()

# Harmony, 15 PCs
srt_obj <- run_branch(srt_obj,
                      reduction = "integrated.sct.harmony",
                      dims = 1:15,
                      prefix = "integrated.sct.harmony_15",
                      resolution = c(0.6, 1.2, 2, 4)
)

gc()

# Harmony, 40 PCs
srt_obj <- run_branch(srt_obj,
                      reduction = "integrated.sct.harmony",
                      dims = 1:40,
                      prefix = "integrated.sct.harmony_40",
                      resolution = c(0.6, 1.2, 2, 4)
)

gc()

# No integration, 15 PCs (PCA space)
srt_obj <- run_branch(srt_obj,
                      reduction = "pca_sct",
                      dims = 1:15,
                      prefix = "pca_sct_15",
                      resolution = c(0.6, 1.2, 2, 4)
)

gc()

# No integration, 40 PCs
srt_obj <- run_branch(srt_obj,
                      reduction = "pca_sct",
                      dims = 1:40,
                      prefix = "pca_sct_40",
                      resolution = c(0.6, 1.2, 2, 4)
)

gc()

## RNA downstream processing

# RPCA, 15 PCs
srt_obj <- run_branch(srt_obj,
                      reduction = "integrated.rna.rpca",
                      dims = 1:15,
                      prefix = "integrated.rna.rpca_15",
                      resolution = c(0.6, 1.2, 2, 4)
)

gc()

# RPCA, 40 PCs
srt_obj <- run_branch(srt_obj,
                      reduction = "integrated.rna.rpca",
                      dims = 1:40,
                      prefix = "integrated.rna.rpca_40",
                      resolution = c(0.6, 1.2, 2, 4)
)


gc()

# Harmony, 15 PCs
srt_obj <- run_branch(srt_obj,
                      reduction = "integrated.rna.harmony",
                      dims = 1:15,
                      prefix = "integrated.rna.harmony_15",
                      resolution = c(0.6, 1.2, 2, 4)
)

gc()

# Harmony, 40 PCs
srt_obj <- run_branch(srt_obj,
                      reduction = "integrated.rna.harmony",
                      dims = 1:40,
                      prefix = "integrated.rna.harmony_40",
                      resolution = c(0.6, 1.2, 2, 4)
)

gc()

# No integration, 15 PCs (PCA space)
srt_obj <- run_branch(srt_obj,
                      reduction = "pca_rna",
                      dims = 1:15,
                      prefix = "pca_rna_15",
                      resolution = c(0.6, 1.2, 2, 4)
)

gc()

# No integration, 40 PCs
srt_obj <- run_branch(srt_obj,
                      reduction = "pca_sct",
                      dims = 1:40,
                      prefix = "pca_rna_40",
                      resolution = c(0.6, 1.2, 2, 4)
)
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

epidermal_markers <- c("Mlig455-003023", "Mlig455-006415", "Mlig455-009402", "Mlig455-015811", "Mlig455-050841", "Mlig455-024763", "Mlig455-032600", "Mlig455-040753")
