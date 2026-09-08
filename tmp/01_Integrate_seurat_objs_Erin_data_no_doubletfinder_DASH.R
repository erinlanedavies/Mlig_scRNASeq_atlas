## 02_Integrate_seurat_objs_Eugene.R
#Author: Francisco Lobo, 05/2024

# Takes as input Seurat objects, integrates them and prints results

#loading libraries
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(dplyr)
library(harmony)
library(speckle)
library(sctransform)
library(scCustomize)
library(parallel)
library(clustree)

##clean workspace
rm(list=ls())
set.seed(42)

#set basedir
basedir <- "/data/Mlig_scRNA_Seq/Francisco/results/"
setwd(basedir)

out_dir <- paste0(basedir, "images/raw_pdf/feat/")

srt_obj <- readRDS("")

# cell_cycle_df <- read.table(file = "../data/metadata/cell_cycle_genes.tsv", sep = "\t", header = TRUE)
#
# cell_cycle_df <- cell_cycle_df[cell_cycle_df$Mlig_geneID != "N.A." & cell_cycle_df$Obs != "low-qual" & cell_cycle_df$Obs != "low-qual?",]
#
# cell_cycle_df$Mlig_geneID <- str_replace(cell_cycle_df$Mlig_geneID, pattern = "_", replacement = "-")
#
# s_feat <- cell_cycle_df$Mlig_geneID[cell_cycle_df$CellCycleState == "S"]
# g2m_feat <- cell_cycle_df$Mlig_geneID[cell_cycle_df$CellCycleState == "G2"]
#
# srt_obj <- CellCycleScoring(srt_obj, s.features = s_feat, g2m.features = g2m_feat)
#
#   sample_list <- SplitObject(srt_obj, split.by="sample_ID")
#
# gc()
#
# options(future.globals.maxSize = 50000 * 1024^2)
#
# # sample_list <- lapply(X = sample_list,
# #                    FUN = SCTransform,
# #                    method = "glmGamPoi",
# #                    vars.to.regress = c("percent.mito", "S.Score", "G2M.Score"),
# #                    return.only.var.genes = FALSE)
#
# gc()
#

srt_obj <- merge(x = sample_list[1],
                 y = sample_list[-1],
                 add.cell.ids = c("Obj1",
                                  "Obj2"
                                  )
                 )

# global_features <- SelectIntegrationFeatures(object.list = sample_list, nfeatures = 3000)
#
# # 3. Split by timepoint
# time_list <- SplitObject(All.sct, split.by = "sample_type")
#
# options(future.globals.maxSize = 50000 * 1024^2)
#
# # 4. Per timepoint, integrate *only* its replicates
# integrated_timepoints <- lapply(names(time_list), function(tp) {
#   options(future.globals.maxSize = 50000 * 1024^2)
#   print(tp)
#   # which samples belong to this time
#   samps <- unique(time_list[[tp]]$sample_ID)
#   subs  <- sample_list[samps]
#
#   # (Optional diagnostic) see if any globals are still missing
#   missing <- lapply(subs, function(x) {
#     setdiff(global_features, rownames(x[["SCT"]]@scale.data))
#   })
#   if (any(lengths(missing) > 0)) {
#     message("Still missing features in SCT slot for time=", tp)
#     print(missing)
#   }
#
#   # Prep, PCA, RPCA‐anchors, Integrate
#   subs <- PrepSCTIntegration(
#     object.list     = subs,
#     anchor.features = global_features,
#     verbose         = FALSE
#   )
#
#   subs <- lapply(subs, RunPCA,
#                  assay    = "SCT",
#                  features = global_features,
#                  verbose  = FALSE)
#
#   anchors <- FindIntegrationAnchors(
#     object.list          = subs,
#     normalization.method = "SCT",
#     anchor.features      = global_features,
#     reduction            = "rpca",
#     dims                 = 1:40
#   )
#
#   integ <- IntegrateData(
#     anchorset            = anchors,
#     normalization.method = "SCT",
#     dims                 = 1:40
#   )
#
#   DefaultAssay(integ) <- "integrated"
#   integ
# })
#
# names(integrated_timepoints) <- names(time_list)
#
# cleaned_timepoints <- lapply(integrated_timepoints, function(obj) {
#   DietSeurat(
#     obj,
#     assays    = "integrated",   # keep only the integrated assay
#     dimreducs = "pca"           # optionally keep the PCA reduction
#   )
# })
#
# # 2. Merge them in one go
# seurat_integrated <- merge(
#   x            = cleaned_timepoints[[1]],
#   y            = cleaned_timepoints[-1],
# #  add.cell.ids = names(cleaned_timepoints),
#   project      = "ReplicateIntegration"
# )
#
# # 3. Re‐run scaling, PCA & UMAP on the merged “integrated” assay
# DefaultAssay(seurat_integrated) <- "integrated"
# seurat_integrated <- FindVariableFeatures(seurat_integrated)
# seurat_integrated <- ScaleData(seurat_integrated, verbose = FALSE)
# seurat_integrated <- RunPCA(seurat_integrated, verbose = FALSE)
# seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:30)
#
# #
# # # 5. Merge all time‐point–integrated objects in one shot
# # seurat_integrated <- merge(
# #   x            = integrated_timepoints[[1]],
# #   y            = integrated_timepoints[-1],
# # #  add.cell.ids = names(integrated_timepoints),
# #   project      = "Timepoint_Replicate_Integration"
# # )
#
# seurat_object_7 <- FindNeighbors(seurat_integrated,
#                                  reduction = "pca",
#                                  dims = 1:40,
#                                  k.param = 20,
#                                  return.neighbor = TRUE,
#                                  nn.method = "annoy",
#                                  annoy.metric = "cosine",
#                                  verbose = FALSE)
#
# seurat_object_8 <- FindNeighbors(seurat_object_7,
#                                  reduction = "pca",
#                                  dims = 1:50,
#                                  k.param = 20,
#                                  compute.SNN = TRUE,
#                                  nn.method = "annoy",
#                                  annoy.metric = "cosine",
#                                  verbose = FALSE)
#
# seurat_integrated_final <- FindClusters(seurat_object_8,
#                                        resolution = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2, 3, 4),
# #                                       resolution = c(0.8),
#                                        algorithm = 4,
#                                        method = "igraph",
#                                        random.seed = 42)
#

#saveRDS(seurat_integrated_final, file = "objs/seurat_integrated_assay_DASH.rds")

#seurat_integrated_final <- readRDS(file = "objs/seurat_integrated_assay_DASH.rds")
#
# tmp <- DimPlot_SCplus2(seurat_integrated_final, split.by = "sample_ID", group.by = "integrated_snn_res.0.2", reduction = "pca")
#
# tmp2 <- DimPlot_SCplus2(seurat_integrated_final, split.by = "sample_ID", group.by = "sample_ID", reduction = "umap")
#
# tmp3 <- DimPlot_SCplus2(seurat_integrated_final, split.by = "sample_type", group.by = "sample_type", reduction = "umap")
#
# plot_by_sample <- plot_grid(tmp2$main_plot, tmp2$facet_plot, rel_widths = c(1, 0.8), ncol = 2)
#
# plot_by_stage <- plot_grid(tmp3$main_plot, tmp3$facet_plot, rel_widths = c(1, 0.8), ncol = 2)
#
# gc()
#
# orig_cells <- Cells(All.sct)
#
# new_cells <- Cells(seurat_integrated_final)
#
# # After merging/integration, collapse the RNA assay layers:
# All.sct <- JoinLayers(All.sct, assay = "RNA")
# gc()
#
# # If your SCT assay also has multiple layers:
# #All.sct[["SCT"]] <- JoinLayers(All.sct, assay = "SCT")
#
# #rna_counts <- GetAssayData(All.sct_new, assay = "RNA", slot = "counts")[, orig_cells]
#
# rna_assay <- GetAssay(object = All.sct, assay = "RNA")
#
# seurat_integrated_final[["RNA"]] <- rna_assay
#
# sct_assay <- GetAssay(object = All.sct, assay = "SCT")
#
# seurat_integrated_final[["SCT"]] <- sct_assay
#
# saveRDS(seurat_integrated_final, file = "objs/seurat_integrated_all_assays_DASH.rds")

seurat_integrated_final <- readRDS(file = "objs/seurat_integrated_all_assays_DASH.rds")

cols_to_keep <- c("orig.ident",
                  "nCount_RNA",
                  "nFeature_RNA",
                  "sample_ID",
                  "sample_type",
                  "percent.mito",
                  "nCount_SCT",
                  "nFeature_SCT",
                  "S.Score",
                  "G2M.Score",
                  "Phase",
                  "nCount_integrated",
                  "nFeature_integrated",
                  "integrated_snn_res.0.2",
                  "integrated_snn_res.0.4",
                  "integrated_snn_res.0.6",
                  "integrated_snn_res.0.8",
                  "integrated_snn_res.1",
                  "integrated_snn_res.1.2",
                  "integrated_snn_res.1.4",
                  "integrated_snn_res.1.6",
                  "integrated_snn_res.1.8",
                  "integrated_snn_res.2",
                  "integrated_snn_res.3",
                  "integrated_snn_res.4")

seurat_integrated_final@meta.data <- seurat_integrated_final@meta.data[, cols_to_keep, drop = FALSE]

Assays(seurat_integrated_final)

DefaultAssay(seurat_integrated_final) <- "RNA"

gc()

seurat_integrated_final$sample_ID <- sub("^S", "", seurat_integrated_final$sample_ID)

## Add module scores to see if results look OK

## Helper function to work with a list of features and remove the annoying "1" that AddModuleScore adds at the end of module names

