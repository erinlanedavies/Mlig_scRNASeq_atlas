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

source("~/bin/R_functions/single_cell.R")

# #load objs
# 
# S17_0_1d_Hatchlings_filt <- readRDS(file = "/data/Mlig_scRNA_Seq/Francisco/results/objs/S17_0_1d_Hatchlings_filt.rds")
# S18_0_1d_Hatchling_filt <- readRDS(file = "/data/Mlig_scRNA_Seq/Francisco/results/objs/S18_0_1d_Hatchling_filt.rds")
# S19_0_1d_Hatchling_filt <- readRDS(file = "/data/Mlig_scRNA_Seq/Francisco/results/objs/S19_0_1d_Hatchling_filt.rds")
# S11_unsorted_3dph_filt <- readRDS(file = "/data/Mlig_scRNA_Seq/Francisco/results/objs/S11_unsorted_3dph_filt.rds")
# S12_unsorted_3dph_filt <- readRDS(file = "/data/Mlig_scRNA_Seq/Francisco/results/objs/S12_unsorted_3dph_filt.rds")
# S16_3d_juveniles_filt <- readRDS(file = "/data/Mlig_scRNA_Seq/Francisco/results/objs/S16_3d_juveniles_filt.rds")
# S8_unsorted_Hatchling_7d_filt <- readRDS(file = "/data/Mlig_scRNA_Seq/Francisco/results/objs/S8_unsorted_Hatchling_7d_filt.rds")
# S15_7d_juveniles_filt <- readRDS(file = "/data/Mlig_scRNA_Seq/Francisco/results/objs/S15_7d_juveniles_filt.rds")
# S9_unsorted_Hatchling_7d_filt <- readRDS(file = "/data/Mlig_scRNA_Seq/Francisco/results/objs/S9_unsorted_Hatchling_7d_filt.rds")
# S1_Adult_Unsorted_filt <- readRDS(file = "/data/Mlig_scRNA_Seq/Francisco/results/objs/S1_Adult_Unsorted_filt.rds")
# S2_Adult_Unsorted_DASH_filt <- readRDS(file = "/data/Mlig_scRNA_Seq/Francisco/results/objs/S2_Adult_Unsorted_DASH_filt.rds")
# S2_Adult_Unsorted_filt <- readRDS(file = "/data/Mlig_scRNA_Seq/Francisco/results/objs/S2_Adult_Unsorted_filt.rds")
# S7_Sorted_Adult_filt <- readRDS(file = "/data/Mlig_scRNA_Seq/Francisco/results/objs/S7_Sorted_Adult_filt.rds")
# 
# #plot doublet detection
# p1 <- DimPlot_scCustom(S17_0_1d_Hatchlings_filt,
#                        aspect_ratio = c(1),
#                        group.by = "DF.classifications_0.25_0.01_511",
#                        label = FALSE,
#                        pt.size = 0.5,
#                        colors_use = DiscretePalette_scCustomize(num_colors = 2, palette = "ditto_seq")) +
#   ggtitle("S17_0_1d_Hatchlings") + NoLegend()
# 
# p2 <- DimPlot_scCustom(S18_0_1d_Hatchling_filt,
#                        aspect_ratio = c(1),
#                        group.by = "DF.classifications_0.25_0.3_627",
#                        label = FALSE,
#                        pt.size = 0.5,
#                        colors_use = DiscretePalette_scCustomize(num_colors = 2, palette = "ditto_seq")) +
#   ggtitle("S18_0_1d_Hatchling") + NoLegend()
# 
# p3 <- DimPlot_scCustom(S19_0_1d_Hatchling_filt,
#                        aspect_ratio = c(1),
#                        group.by = "DF.classifications_0.25_0.27_504",
#                        label = FALSE,
#                        pt.size = 0.5,
#                        colors_use = DiscretePalette_scCustomize(num_colors = 2, palette = "ditto_seq")) +
#   ggtitle("S19_0_1d_Hatchling") + NoLegend()
# 
# p4 <- DimPlot_scCustom(S11_unsorted_3dph_filt,
#                        aspect_ratio = c(1),
#                        group.by = "DF.classifications_0.25_0.26_545",
#                        label = FALSE,
#                        pt.size = 0.5,
#                        colors_use = DiscretePalette_scCustomize(num_colors = 2, palette = "ditto_seq")) +
#   ggtitle("S11_unsorted_3dph") + NoLegend()
# 
# p5 <- DimPlot_scCustom(S12_unsorted_3dph_filt,
#                        aspect_ratio = c(1),
#                        group.by = "DF.classifications_0.25_0.02_527",
#                        label = FALSE,
#                        pt.size = 0.5,
#                        colors_use = DiscretePalette_scCustomize(num_colors = 2, palette = "ditto_seq")) +
#   ggtitle("S12_unsorted_3dph") + NoLegend()
# 
# p6 <- DimPlot_scCustom(S16_3d_juveniles_filt,
#                        aspect_ratio = c(1),
#                        group.by = "DF.classifications_0.25_0.02_466",
#                        label = FALSE,
#                        pt.size = 0.5,
#                        colors_use = DiscretePalette_scCustomize(num_colors = 2, palette = "ditto_seq")) +
#   ggtitle("S16_3d_juveniles") + NoLegend()
# 
# p7 <- DimPlot_scCustom(S8_unsorted_Hatchling_7d_filt,
#                        aspect_ratio = c(1),
#                        group.by = "DF.classifications_0.25_0.21_685",
#                        label = FALSE,
#                        pt.size = 0.5,
#                        colors_use = DiscretePalette_scCustomize(num_colors = 2, palette = "ditto_seq")) +
#   ggtitle("S8_unsorted_Hatchling_7d") + NoLegend()
# 
# p8 <- DimPlot_scCustom(S15_7d_juveniles_filt,
#                        aspect_ratio = c(1),
#                        group.by = "DF.classifications_0.25_0.03_608",
#                        label = FALSE,
#                        pt.size = 0.5,
#                        colors_use = DiscretePalette_scCustomize(num_colors = 2, palette = "ditto_seq")) +
#   ggtitle("S15_7d_juveniles") + NoLegend()
# 
# p9 <- DimPlot_scCustom(S9_unsorted_Hatchling_7d_filt,
#                        aspect_ratio = c(1),
#                        group.by = "DF.classifications_0.25_0.1_555",
#                        label = FALSE,
#                        pt.size = 0.5,
#                        colors_use = DiscretePalette_scCustomize(num_colors = 2, palette = "ditto_seq")) +
#   ggtitle("S9_unsorted_Hatchling_7d") + NoLegend()
# 
# p10 <- DimPlot_scCustom(S1_Adult_Unsorted_filt,
#                        aspect_ratio = c(1),
#                        group.by = "DF.classifications_0.25_0.005_287",
#                        label = FALSE,
#                        pt.size = 0.5,
#                        colors_use = DiscretePalette_scCustomize(num_colors = 2, palette = "ditto_seq")) +
#   ggtitle("S1_Adult_Unsorted") + NoLegend()
# 
# p11 <- DimPlot_scCustom(S2_Adult_Unsorted_DASH_filt,
#                        aspect_ratio = c(1),
#                        group.by = "DF.classifications_0.25_0.005_283",
#                        label = FALSE,
#                        pt.size = 0.5,
#                        colors_use = DiscretePalette_scCustomize(num_colors = 2, palette = "ditto_seq")) +
#   ggtitle("S2_Adult_Unsorted_DASH") + NoLegend()
# 
# p12 <- DimPlot_scCustom(S7_Sorted_Adult_filt,
#                        aspect_ratio = c(1),
#                        group.by = "DF.classifications_0.25_0.005_405",
#                        label = FALSE,
#                        pt.size = 0.5,
#                        colors_use = DiscretePalette_scCustomize(num_colors = 2, palette = "ditto_seq")) +
#   ggtitle("S7_Sorted_Adult") + NoLegend()
# 
# outfile_path <- paste0(out_dir, "DoubletFinder_round_two.png")
# png(filename = outfile_path, units = "cm", width = 40, height = 40, res = 200)
# plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12)
# dev.off()
# 
# #creating a dataset of filtered data
# All_filtered <- merge(S17_0_1d_Hatchlings_filt,
#                       y =c(S18_0_1d_Hatchling_filt,
#                            S19_0_1d_Hatchling_filt,
#                            S11_unsorted_3dph_filt,
#                            S12_unsorted_3dph_filt,
#                            S16_3d_juveniles_filt,
#                            S8_unsorted_Hatchling_7d_filt,
#                            S15_7d_juveniles_filt,
#                            S9_unsorted_Hatchling_7d_filt,
#                            S1_Adult_Unsorted_filt,
#                            S2_Adult_Unsorted_DASH_filt,
#                            S2_Adult_Unsorted_filt,
#                            S7_Sorted_Adult_filt),
#                       add.cell.ids = c("S17_0_1d_Hatchlings_filt",
#                                        "S18_0_1d_Hatchling_filt",
#                                        "S19_0_1d_Hatchling_filt",
#                                        "S11_unsorted_3dph_filt",
#                                        "S12_unsorted_3dph_filt",
#                                        "S16_3d_juveniles_filt",
#                                        "S8_unsorted_Hatchling_7d_filt",
#                                        "S15_7d_juveniles_filt",
#                                        "S9_unsorted_Hatchling_7d_filt",
#                                        "S1_Adult_Unsorted_filt",
#                                        "S2_Adult_Unsorted_DASH_filt",
#                                        "S2_Adult_Unsorted_filt",
#                                        "S7_Sorted_Adult_filt"),
#                       project = "scRNA_Seq_M_lignano")
# 
# VlnPlot(All_filtered, group.by = "orig.ident", features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, alpha = 0, raster = FALSE)


# #Doublet removal
# 
# tmp <- subset(S17_0_1d_Hatchlings_filt, subset = `DF.classifications_0.25_0.01_511` == "Singlet")
# DefaultAssay(tmp) <- "RNA"
# tmp <- DietSeurat(tmp, assays = "RNA", dimreducs = NULL)
# S17_0_1d_Hatchlings_filt <- tmp
# 
# tmp <- subset(S18_0_1d_Hatchling_filt, subset = `DF.classifications_0.25_0.3_627` == "Singlet")
# DefaultAssay(tmp) <- "RNA"
# tmp <- DietSeurat(tmp, assays = "RNA", dimreducs = NULL)
# S18_0_1d_Hatchling_filt <- tmp
# 
# tmp <- subset(S19_0_1d_Hatchling_filt, subset = `DF.classifications_0.25_0.27_504` == "Singlet")
# DefaultAssay(tmp) <- "RNA"
# tmp <- DietSeurat(tmp, assays = "RNA", dimreducs = NULL)
# S19_0_1d_Hatchling_filt <- tmp
# 
# tmp <- subset(S11_unsorted_3dph_filt, subset = `DF.classifications_0.25_0.26_545` == "Singlet")
# DefaultAssay(tmp) <- "RNA"
# tmp <- DietSeurat(tmp, assays = "RNA", dimreducs = NULL)
# S11_unsorted_3dph_filt <- tmp
# 
# tmp <- subset(S12_unsorted_3dph_filt, subset = `DF.classifications_0.25_0.02_527` == "Singlet")
# DefaultAssay(tmp) <- "RNA"
# tmp <- DietSeurat(tmp, assays = "RNA", dimreducs = NULL)
# S12_unsorted_3dph_filt <- tmp
# 
# tmp <- subset(S16_3d_juveniles_filt, subset = `DF.classifications_0.25_0.02_466` == "Singlet")
# DefaultAssay(tmp) <- "RNA"
# tmp <- DietSeurat(tmp, assays = "RNA", dimreducs = NULL)
# S16_3d_juveniles_filt <- tmp
# 
# tmp <- subset(S8_unsorted_Hatchling_7d_filt, subset = `DF.classifications_0.25_0.21_685` == "Singlet")
# DefaultAssay(tmp) <- "RNA"
# tmp <- DietSeurat(tmp, assays = "RNA", dimreducs = NULL)
# S8_unsorted_Hatchling_7d_filt <- tmp
# 
# tmp <- subset(S15_7d_juveniles_filt, subset = `DF.classifications_0.25_0.03_608` == "Singlet")
# DefaultAssay(tmp) <- "RNA"
# tmp <- DietSeurat(tmp, assays = "RNA", dimreducs = NULL)
# S15_7d_juveniles_filt <- tmp
# 
# tmp <- subset(S9_unsorted_Hatchling_7d_filt, subset = `DF.classifications_0.25_0.1_555` == "Singlet")
# DefaultAssay(tmp) <- "RNA"
# tmp <- DietSeurat(tmp, assays = "RNA", dimreducs = NULL)
# S9_unsorted_Hatchling_7d_filt <- tmp
# 
# tmp <- subset(S1_Adult_Unsorted_filt, subset = `DF.classifications_0.25_0.005_287` == "Singlet")
# DefaultAssay(tmp) <- "RNA"
# tmp <- DietSeurat(tmp, assays = "RNA", dimreducs = NULL)
# S1_Adult_Unsorted_filt <- tmp
# 
# tmp <- subset(S2_Adult_Unsorted_DASH_filt, subset = `DF.classifications_0.25_0.005_283` == "Singlet")
# DefaultAssay(tmp) <- "RNA"
# tmp <- DietSeurat(tmp, assays = "RNA", dimreducs = NULL)
# S2_Adult_Unsorted_DASH_filt <- tmp
# 
# tmp <- subset(S7_Sorted_Adult_filt, subset = `DF.classifications_0.25_0.005_405` == "Singlet")
# DefaultAssay(tmp) <- "RNA"
# tmp <- DietSeurat(tmp, assays = "RNA", dimreducs = NULL)
# S7_Sorted_Adult_filt <- tmp
# 
# gc()
# 
# cell_cycle_df <- read.table(file = "/data/Mlig_scRNA_Seq/Francisco/data/cell_cycle_genes.tsv", sep = "\t", header = TRUE)
# 
# cell_cycle_df <- cell_cycle_df[cell_cycle_df$Mlig_geneID != "N.A." & cell_cycle_df$Obs != "low-qual" & cell_cycle_df$Obs != "low-qual?",]
# 
# cell_cycle_df$Mlig_geneID <- str_replace(cell_cycle_df$Mlig_geneID, pattern = "_", replacement = "-")
# 
# s_feat <- cell_cycle_df$Mlig_geneID[cell_cycle_df$CellCycleState == "S"]
# g2m_feat <- cell_cycle_df$Mlig_geneID[cell_cycle_df$CellCycleState == "G2"]
# 
# All_filtered <- CellCycleScoring(All_filtered, s.features = s_feat, g2m.features = g2m_feat)
# 
# if (variable_to_integrate == "sample_type") {
#   All.list <- SplitObject(All_filtered, split.by="sample_type")
# }
# 
# if (variable_to_integrate == "sample_ID") {
#   All.list <- SplitObject(All_filtered, split.by="sample_ID")
# }
# 
# rm(All_filtered)
# 
# gc()
# 
# options(future.globals.maxSize = 50000 * 1024^2)
# 
# # All.list <- lapply(X = All.list, 
# #                    FUN = SCTransform,
# #                    method = "glmGamPoi",
# #                    vars.to.regress = c("percent.mito", "S.Score", "G2M.Score"),
# #                    return.only.var.genes = FALSE)
# 
# All.list <- lapply(X = All.list, 
#                    FUN = SCTransform,
#                    method = "glmGamPoi",
#                    vars.to.regress = c("percent.mito"),
#                    return.only.var.genes = FALSE)
# 
# gc()
# 
# var.features <- SelectIntegrationFeatures(object.list = All.list, nfeatures = 3000)
# http://localhost:33733/auth-sign-in?user=pereiralobof2&password=tY9JNTo7vTd3C3DriPiYuVCP
# 
# #adding variable features
# 
# DefaultAssay(All_filtered) <- "RNA"
# 
# All.no_sct <- DietSeurat(All_filtered, assays = c("RNA"))
# 
# All.no_sct[["RNA"]] <- JoinLayers(All.no_sct[["RNA"]], layers = c("counts"))
# 
# All.no_sct <- NormalizeData(All.no_sct, normalization.method = "LogNormalize", scale.factor = 10000)
# 
# All.no_sct <- ScaleData(All.no_sct)#, vars.to.regress = "percent.mito")
# 
# All.no_sct <- FindVariableFeatures(All.no_sct)
# 
# gc()
# 
# #run PCA
# All.no_sct <- RunPCA(All.no_sct, npcs = 100, verbose = TRUE)
# 
# seurat_object_6 <- RunUMAP(All.no_sct,
#                            reduction = "pca",
#                            reduction.name = "umap.pca",
#                            dims = 1:40,
#                            return.model = TRUE,
# #                           n.neighbors = 30,
# #                           min.dist = 0.1,
# #                           spread = 5,
#                            seed.use = 42)
# 
# seurat_object_6 <- FindNeighbors(seurat_object_6,
#                                  reduction = "pca",
#                                  dims = 1:40,
#                                  k.param = 10,
#                                  return.neighbor = TRUE,
#                                  nn.method = "annoy",
#                                  annoy.metric = "cosine",
#                                  verbose = FALSE)
# 
# seurat_object_6 <- FindNeighbors(seurat_object_6,
#                                  reduction = "pca",
#                                  dims = 1:50,
#                                  k.param = 10,
#                                  compute.SNN = TRUE,
#                                  nn.method = "annoy",
#                                  annoy.metric = "cosine",
#                                  verbose = FALSE)
# 
# seurat_object_6 <- FindClusters(seurat_object_6,
#                                        resolution = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2, 3, 4),
# #                                       resolution = c(0.8),
#                                        algorithm = 4,
#                                        method = "igraph", 
#                                        random.seed = 42)
# 
# 
# gc()
# 
# 
# tmp4 <- DimPlot_SCplus2(seurat_object_6, split.by = "sample_ID", group.by = "RNA_snn_res.0.8", reduction = "pca")
# 
# tmp5 <- DimPlot_SCplus2(seurat_object_6, split.by = "sample_ID", group.by = "RNA_snn_res.0.8", reduction = "umap.pca")
# 
# tmp6 <- DimPlot_SCplus2(seurat_object_6, split.by = "sample_type", group.by = "sample_type", reduction = "umap.pca")
# 
# tmp4$facet_plot
# tmp5$facet_plot
# tmp6$facet_plot
# 
# saveRDS(seurat_object_6, file = "objs/seurat_no_integrated_DASH_noDASH.rds")
# 
# seurat_object_SCT_only$sample_type <- factor(seurat_object_SCT_only$sample_type, levels = c("01D", "3D", "7D", "Ad"))
# 
# seurat_object_SCT_only$sample_ID <- factor(seurat_object_SCT_only$sample_ID, levels = c("17_0_1d_Hatchlings", "18_0_1d_Hatchling", "19_0_1d_Hatchling", "S11_unsorted_3dph", "S12_unsorted_3dph", "S16_3d_juveniles", "8_unsorted_Hatchling_7d", "9_unsorted_Hatchling_7d", "15_7d_juveniles", "1_Adult_Unsorted", "2_Adult_Unsorted_DASH", "7_Sorted_Adult"))
# 
# SCT_plot_by_sample <- DimPlot_clusters_4(seurat_object_SCT_only, metadata_col = "sample_ID", color_factor = "sample_ID", reduction = "umap.pca")
# 
# SCT_plot_by_sample_final <- plot_grid(SCT_plot_by_sample$main_plot, SCT_plot_by_sample$facet_plot, rel_widths = c(1, 0.8), ncol = 2)
# 
# png("images/UMAP_clusters_by_sample_SCT_only_no_doubletfinder.png", width = 40, height = 40, units = "cm", res = 200)
# plot(SCT_plot_by_sample_final)
# dev.off()
# 
# SCT_plot_by_stage <- DimPlot_clusters_4(seurat_object_SCT_only, metadata_col = "sample_type", color_factor = "sample_type", reduction = "umap.pca")
# 
# SCT_plot_by_stage_final <- plot_grid(SCT_plot_by_stage$main_plot, SCT_plot_by_stage$facet_plot, rel_widths = c(1, 0.8), ncol = 2)
# 
# png("images/UMAP_clusters_by_stage_SCT_only_no_doubletfinder.png", width = 40, height = 40, units = "cm", res = 200)
# plot(SCT_plot_by_stage_final)
# dev.off()
# 
# SCT_plot_by_cluster <- DimPlot_clusters_4(seurat_object_SCT_only, metadata_col = "SCT_snn_res.0.8", color_factor = "SCT_snn_res.0.8", reduction = "umap.pca")
# 
# SCT_plot_by_cluster_final <- plot_grid(SCT_plot_by_cluster$main_plot, SCT_plot_by_cluster$facet_plot, rel_widths = c(1, 0.8), ncol = 2)
# 
# png("images/UMAP_clusters_by_cluster_SCT_only_no_doubletfinder.png", width = 40, height = 40, units = "cm", res = 200)
# plot(SCT_plot_by_cluster_final)
# dev.off()
# 
# 
# 
# #run harmony
# #All.sct <- RunHarmony(All.sct, assay.use="SCT", group.by.vars = "origin", plot_convergence = TRUE)
# seurat_object_5 <- RunHarmony(All.sct, assay.use="SCT", group.by.vars = c("sample_type"), plot_convergence = TRUE)
# #All.sct <- RunHarmony(All.sct, assay.use="SCT", group.by.vars = c("sample_ID"), plot_convergence = TRUE)
# 
# gc()
# 
# seurat_object_6 <- RunUMAP(seurat_object_5,
#                            reduction = "harmony",
#                            reduction.name = "umap.harmony",
#                            dims = 1:40,
#                            return.model = TRUE,
#                            # n.neighbors = 30,
#                            # min.dist = 0.01,
#                            # spread = 5,
#                            seed.use = 42)
# 
# seurat_object_7 <- FindNeighbors(seurat_object_6,
#                                  reduction = "harmony",
#                                  dims = 1:40,
#                                  k.param = 20,
#                                  return.neighbor = TRUE,
#                                  nn.method = "annoy",
#                                  annoy.metric = "cosine",
#                                  verbose = FALSE)
# 
# seurat_object_8 <- FindNeighbors(seurat_object_7,
#                                  reduction = "harmony",
#                                  dims = 1:50,
#                                  k.param = 20,
#                                  compute.SNN = TRUE,
#                                  nn.method = "annoy",
#                                  annoy.metric = "cosine",
#                                  verbose = FALSE)
# 
# seurat_object_SCT_harmony <- FindClusters(seurat_object_8,
#                                        #                                resolution = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2),
#                                        resolution = c(0.8),
#                                        algorithm = 4,
#                                        method = "igraph", 
#                                        random.seed = 42)
# 
# saveRDS(seurat_object_SCT_harmony, file = "objs/srt_obj_SCT_harmony_no_doubletfinder.rds")
# 
# seurat_object_SCT_harmony$sample_type <- factor(seurat_object_SCT_harmony$sample_type, levels = c("01D", "3D", "7D", "Ad"))
# 
# seurat_object_SCT_harmony$sample_ID <- factor(seurat_object_SCT_harmony$sample_ID, levels = c("17_0_1d_Hatchlings", "18_0_1d_Hatchling", "19_0_1d_Hatchling", "S11_unsorted_3dph", "S12_unsorted_3dph", "S16_3d_juveniles", "8_unsorted_Hatchling_7d", "9_unsorted_Hatchling_7d", "15_7d_juveniles", "1_Adult_Unsorted", "2_Adult_Unsorted_DASH", "7_Sorted_Adult"))
# 
# SCT_plot_by_sample <- DimPlot_clusters_4(seurat_object_SCT_harmony, metadata_col = "sample_ID", color_factor = "sample_ID", reduction = "umap.harmony")
# 
# SCT_plot_by_sample_final <- plot_grid(SCT_plot_by_sample$main_plot, SCT_plot_by_sample$facet_plot, rel_widths = c(1, 0.8), ncol = 2)
# 
# png("images/UMAP_clusters_by_sample_SCT_harmony_no_doubletfinder.png", width = 40, height = 40, units = "cm", res = 200)
# plot(SCT_plot_by_sample_final)
# dev.off()
# 
# SCT_plot_by_stage <- DimPlot_clusters_4(seurat_object_SCT_harmony, metadata_col = "sample_type", color_factor = "sample_type", reduction = "umap.harmony")
# 
# SCT_plot_by_stage_final <- plot_grid(SCT_plot_by_stage$main_plot, SCT_plot_by_stage$facet_plot, rel_widths = c(1, 0.8), ncol = 2)
# 
# png("images/UMAP_clusters_by_stage_SCT_harmony_no_doubletfinder.png", width = 40, height = 40, units = "cm", res = 200)
# plot(SCT_plot_by_stage_final)
# dev.off()
# 
# SCT_plot_by_cluster <- DimPlot_clusters_4(seurat_object_SCT_harmony, metadata_col = "SCT_snn_res.0.8", color_factor = "SCT_snn_res.0.8", reduction = "umap.harmony")
# 
# SCT_plot_by_cluster_final <- plot_grid(SCT_plot_by_cluster$main_plot, SCT_plot_by_cluster$facet_plot, rel_widths = c(1, 0.8), ncol = 2)
# 
# png("images/UMAP_clusters_by_cluster_SCT_harmony_no_doubletfinder.png", width = 40, height = 40, units = "cm", res = 200)
# plot(SCT_plot_by_cluster_final)
# dev.off()
# 
# 
# 
# ##rPCA plus SCT
# 
# # 1. Split by sample and SCTransform *with* full residuals
# sample_list <- All.list
# 
# # 2. Pick your global features
# global_features <- var.features
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

#load list of curated marker genes
marker_genes <- read.table("/data/Mlig_scRNA_Seq/Francisco/data/M_lignano_marker_genes.txt", sep = "\t", header = TRUE, quote = "")

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
neural_cell_in_dorsal_tail <- marker_genes$Gene_ID[marker_genes$Cell_tissue_type == "neural_cell_in_dorsal_tail"]
neural_cell_in_tail_paired <- marker_genes$Gene_ID[marker_genes$Cell_tissue_type == "neural_cell_in_tail_(paired)"]
prostate <- marker_genes$Gene_ID[marker_genes$Cell_tissue_type == "prostate"]
rhabdite <- marker_genes$Gene_ID[marker_genes$Cell_tissue_type == "rhabdite-containing_cell"]
secretory_cell_of_adhesive_organ <- marker_genes$Gene_ID[marker_genes$Cell_tissue_type == "secretory_cell_of_adhesive_organ"]
stem_cells_progenitors <- marker_genes$Gene_ID[marker_genes$Cell_tissue_type == "stem_cells_progenitors"]

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
                          neural_cell_in_dorsal_tail,
                          neural_cell_in_tail_paired,
                          prostate,
                          rhabdite,
                          secretory_cell_of_adhesive_organ,
                          stem_cells_progenitors)

names(marker_genes_list) <- c("anchor_cell_of_adhesive_organ", "cement_gland_cell", "epidermal_secretory", "female_antrum", "female_germline_ovary", "gut", "gut_specialized", "male_germline_testis", "muscle", "nervous_system", "neural_cell_in_dorsal_tail", "neural_cell_in_tail_paired", "prostate", "rhabdite", "secretory_cell_of_adhesive_organ", "progenitor")

seurat_integrated_final <- NormalizeData(seurat_integrated_final, normalization.method = "LogNormalize", scale.factor = 10000)

seurat_integrated_final <- add_module_scores(seurat_integrated_final, marker_genes_list, assay = "RNA")

seurat_integrated_final$sample_type <- factor(seurat_integrated_final$sample_type, levels = c("01D", "3D", "7D", "Ad"))

seurat_integrated_final$sample_ID <- factor(seurat_integrated_final$sample_ID, levels = c("17_0_1d_Hatchlings",
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

saveRDS(seurat_integrated_final, file = "../results/objs/seurat_integrated_final_plus_module_scores.rds")

# #FeaturePlot_scCustom(All.sct, features = "male_germline_testis_score1")
# 
# mod_cols <- c("anchor_cell_of_adhesive_organ", "cement_gland_cell", "epidermal_secretory", "female_antrum", "female_germline_ovary", "gut", "gut_specialized", "male_germline_testis", "muscle", "nervous_system", "neural_cell_in_dorsal_tail", "neural_cell_in_tail_paired", "prostate", "rhabdite", "secretory_cell_of_adhesive_organ", "progenitor")
# seurat_integrated_final$anchor_cell_of_adhesive_organ <- seurat_integrated_final$anchor_cell_of_adhesive_organ_score1
# #seurat_integrated_final$anchor_cell_of_adhesive_organ <- (seurat_integrated_final$anchor_cell_of_adhesive_organ - mean(seurat_integrated_final$anchor_cell_of_adhesive_organ)) / sd(seurat_integrated_final$anchor_cell_of_adhesive_organ)
# #delete old metadata value
# seurat_integrated_final$anchor_cell_of_adhesive_organ_score1 <- NULL
# 
# 
# seurat_integrated_final$cement_gland_cell <- seurat_integrated_final$cement_gland_cell_score1
# #seurat_integrated_final$cement_gland_cell <- (seurat_integrated_final$cement_gland_cell - mean(seurat_integrated_final$cement_gland_cell)) / sd(seurat_integrated_final$cement_gland_cell)
# seurat_integrated_final$cement_gland_cell_score1 <- NULL
# 
# seurat_integrated_final$epidermal_secretory <- seurat_integrated_final$epidermal_secretory_score1
# #seurat_integrated_final$epidermal_secretory <- (seurat_integrated_final$epidermal_secretory - mean(seurat_integrated_final$epidermal_secretory)) / sd(seurat_integrated_final$epidermal_secretory)
# seurat_integrated_final$epidermal_secretory_score1 <- NULL
# 
# seurat_integrated_final$female_antrum <- seurat_integrated_final$female_antrum_score1
# #seurat_integrated_final$female_antrum <- (seurat_integrated_final$female_antrum - mean(seurat_integrated_final$female_antrum)) / sd(seurat_integrated_final$female_antrum)
# seurat_integrated_final$female_antrum_score1 <- NULL
# 
# seurat_integrated_final$female_germline_ovary <- seurat_integrated_final$female_germline_ovary_score1
# #seurat_integrated_final$female_germline_ovary <- (seurat_integrated_final$female_germline_ovary - mean(seurat_integrated_final$female_germline_ovary)) / sd(seurat_integrated_final$female_germline_ovary)
# seurat_integrated_final$female_germline_ovary_score1 <- NULL
# 
# seurat_integrated_final$gut <- seurat_integrated_final$gut_score1
# #seurat_integrated_final$gut <- (seurat_integrated_final$gut - mean(seurat_integrated_final$gut)) / sd(seurat_integrated_final$gut)
# seurat_integrated_final$gut_score1 <- NULL
# 
# seurat_integrated_final$gut_specialized <- seurat_integrated_final$gut_specialized_score1
# #seurat_integrated_final$gut_specialized <- (seurat_integrated_final$gut_specialized - mean(seurat_integrated_final$gut_specialized)) / sd(seurat_integrated_final$gut_specialized)
# seurat_integrated_final$gut_specialized_score1 <- NULL
# 
# seurat_integrated_final$male_germline_testis <- seurat_integrated_final$male_germline_testis_score1
# #seurat_integrated_final$male_germline_testis <- (seurat_integrated_final$male_germline_testis - mean(seurat_integrated_final$male_germline_testis)) / sd(seurat_integrated_final$male_germline_testis)
# seurat_integrated_final$male_germline_testis_score1 <- NULL
# 
# seurat_integrated_final$muscle <- seurat_integrated_final$muscle_score1
# #seurat_integrated_final$muscle <- (seurat_integrated_final$muscle - mean(seurat_integrated_final$muscle)) / sd(seurat_integrated_final$muscle)
# seurat_integrated_final$muscle_score1 <- NULL
# 
# seurat_integrated_final$nervous_system <- seurat_integrated_final$nervous_system_score1
# #seurat_integrated_final$nervous_system <- (seurat_integrated_final$nervous_system - mean(seurat_integrated_final$nervous_system)) / sd(seurat_integrated_final$nervous_system)
# seurat_integrated_final$nervous_system_score1 <- NULL
# 
# seurat_integrated_final$neural_cell_in_dorsal_tail <- seurat_integrated_final$neural_cell_in_dorsal_tail_score1
# #seurat_integrated_final$neural_cell_in_dorsal_tail <- (seurat_integrated_final$neural_cell_in_dorsal_tail - mean(seurat_integrated_final$neural_cell_in_dorsal_tail)) / sd(seurat_integrated_final$neural_cell_in_dorsal_tail)
# seurat_integrated_final$neural_cell_in_dorsal_tail_score1 <- NULL
# 
# seurat_integrated_final$neural_cell_in_tail_paired <- seurat_integrated_final$neural_cell_in_tail_paired_score1
# #seurat_integrated_final$neural_cell_in_tail_paired <- (seurat_integrated_final$neural_cell_in_tail_paired - mean(seurat_integrated_final$neural_cell_in_tail_paired)) / sd(seurat_integrated_final$neural_cell_in_tail_paired)
# seurat_integrated_final$neural_cell_in_tail_paired_score1 <- NULL
# 
# seurat_integrated_final$prostate <- seurat_integrated_final$prostate_score1
# #seurat_integrated_final$prostate <- (seurat_integrated_final$prostate - mean(seurat_integrated_final$prostate)) / sd(seurat_integrated_final$prostate)
# seurat_integrated_final$prostate_score1 <- NULL
# 
# seurat_integrated_final$rhabdite <- seurat_integrated_final$rhabdite_score1
# #seurat_integrated_final$rhabdite <- (seurat_integrated_final$rhabdite - mean(seurat_integrated_final$rhabdite)) / sd(seurat_integrated_final$rhabdite)
# seurat_integrated_final$rhabdite_score1 <- NULL
# 
# seurat_integrated_final$secretory_cell_of_adhesive_organ <- seurat_integrated_final$secretory_cell_of_adhesive_organ_score1
# #seurat_integrated_final$secretory_cell_of_adhesive_organ <- (seurat_integrated_final$secretory_cell_of_adhesive_organ - mean(seurat_integrated_final$secretory_cell_of_adhesive_organ)) / sd(seurat_integrated_final$secretory_cell_of_adhesive_organ)
# seurat_integrated_final$secretory_cell_of_adhesive_organ_score1 <- NULL
# 
# seurat_integrated_final$progenitor <- seurat_integrated_final$progenitor_score1
# #seurat_integrated_final$progenitor <- (seurat_integrated_final$progenitor - mean(seurat_integrated_final$progenitor)) / sd(seurat_integrated_final$progenitor)
# seurat_integrated_final$progenitor_score1 <- NULL
# 
# 
# #defining cell types to be considered
# features <- c("cement_gland_cell",
#               "epidermal_secretory",
#               "female_germline_ovary",
#               "female_antrum",
#               "gut",
#               "gut_specialized",
#               "male_germline_testis",
#               "muscle",
#               "nervous_system",
#               "prostate",
#               "rhabdite",
#               "progenitor",
#               "anchor_cell_of_adhesive_organ",
#               "secretory_cell_of_adhesive_organ",
#               "neural_cell_in_dorsal_tail",
#               "neural_cell_in_tail_paired")
# 
# #plot
# png(file = "images/FeaturePlot_scores_manual_curation_integration_biological_replicates_NODASH.png", width = 40, height = 40, units = "cm", res = 600)
# FeaturePlot_scCustom(seurat_integrated_final, features = features, aspect_ratio = 1, num_columns = 4)
# dev.off()

Idents(seurat_integrated_final) <- "integrated_snn_res.1.2"

DimPlot_scCustom(seurat_integrated_final, aspect_ratio = 1)
# 
# DotPlot_scCustom(seurat_integrated_final, features = features, x_lab_rotate = TRUE) +  theme(
#   plot.margin = margin(t = 1, r = 1, b = 1, l = 3, unit = "cm") # Increased left margin
# )
# 
# Idents(seurat_integrated_final) <- "integrated_snn_res.0.4"
# 
# DimPlot_scCustom(seurat_integrated_final, aspect_ratio = 1)
# 
# DotPlot_scCustom(seurat_integrated_final, features = features, x_lab_rotate = TRUE) +  theme(
#   plot.margin = margin(t = 1, r = 1, b = 1, l = 3, unit = "cm") # Increased left margin
# )
# 
# for (resolution in resolutions) {
#   print(resolution)
#   curr_res <- paste0("integrated_snn_res.", resolution)
#   Idents(All.sct) <- curr_res
#   outfile <- paste0("results/images/DotPlot_clusters_", resolution, "_SCT_only_no_doubletfinder_integration_biological_replicates.png")
#   #  SCT_plot_by_group <- DimPlot_SCplus(All.sct, split.by = curr_res, group.by = curr_res, reduction = "umap")
#   
#   p <-   
#     DotPlot_scCustom(All.sct, features = c("anchor_cell_of_adhesive_organ",
#                                            "cement_gland_cell",
#                                            "epidermal_secretory",
#                                            "female_antrum",
#                                            "female_germline_ovary",
#                                            "gut",
#                                            "gut_specialized",
#                                            "male_germline_testis",
#                                            "muscle",
#                                            "nervous_system",
#                                            "neural_cell_in_dorsal_tail",
#                                            "neural_cell_in_tail_paired",
#                                            "prostate",
#                                            "rhabdite",
#                                            "secretory_cell_of_adhesive_organ",
#                                            "progenitor"),
#                      x_lab_rotate = TRUE) + theme(
#                        plot.margin = unit(c(1, 1, 1, 3), "cm")
#                      )
#   
#   
#   
#   png(outfile, width = 35, height = 45, units = "cm", res = 200)
#   print(p)
#   dev.off()
# }
# 
# # Remove a few cell types as discussed with Matt & Erin
# features <- c(#"anchor_cell_of_adhesive_organ",
#   "cement_gland_cell",
#   "epidermal_secretory",
#   #"female_antrum",
#   "female_germline_ovary",
#   "gut",
#   "gut_specialized",
#   "male_germline_testis",
#   "muscle",
#   "nervous_system",
#   #              "neural_cell_in_dorsal_tail",
#   #              "neural_cell_in_tail_paired",
#   "prostate",
#   "rhabdite",
#   #              "secretory_cell_of_adhesive_organ",
#   "progenitor")
# 
# 
# # Ensure all features are in the meta.data
# missing_features <- setdiff(features, colnames(All.sct@meta.data))
# 
# if (length(missing_features) > 0) {
#   stop(paste("The following features are not found in meta.data:", paste(missing_features, collapse = ", ")))
# }
# 
# # Calculate the average values for each feature by cluster
# avg_features_per_cluster <- All.sct@meta.data %>%
#   dplyr::group_by(integrated_snn_res.4) %>%
#   #  dplyr::summarise(across(all_of(features), ~ mean(.x, na.rm = FALSE), .names = "{col}"))
#   dplyr::summarise(across(all_of(features), ~ median(.x, na.rm = FALSE), .names = "{col}"))
# 
# # # Calculate the median values for each feature by cluster
# # avg_features_per_cluster <- All.sct@meta.data %>%
# #   dplyr::group_by(seurat_clusters) %>%
# #   dplyr::summarise(across(all_of(features), ~ median(.x, na.rm = FALSE), .names = "{col}"))
# 
# avg_features_per_cluster <- as.data.frame(avg_features_per_cluster)
# 
# rownames(avg_features_per_cluster) <- avg_features_per_cluster$integrated_snn_res.4
# 
# avg_features_per_cluster <- avg_features_per_cluster[,2:length(colnames(avg_features_per_cluster))]
# 
# cluster_to_max_feature <- vec <- character(length = length(rownames(avg_features_per_cluster)))
# 
# names(cluster_to_max_feature) <- rownames(avg_features_per_cluster)
# 
# #iterate through features (module scores) and clusters to select, for each cluster, its ID based on average module score
# for (i in rownames(avg_features_per_cluster)) {
#   max_val <- -10000
#   cluster_ID <- "NA"
#   for (i2 in colnames(avg_features_per_cluster)) {
#     if (avg_features_per_cluster[i,i2] > max_val) {
#       if (avg_features_per_cluster[i,i2] > 0) {
#         #      if (avg_features_per_cluster[i,i2] > 0.01) { #cutoff to associate
#         cluster_ID <- i2
#         max_val <- avg_features_per_cluster[i,i2]
#       }
#     }
#   }
#   cluster_to_max_feature[i] <- cluster_ID
#   print(paste0("Cluster ID: ", i, " - annotation: ", cluster_ID, " - value:", max_val))
# }
# 
# 
# #create a new metadata layer with manual curation cluster IDs
# All.sct$cluster_ID_manual_curation <- unname(cluster_to_max_feature[All.sct$integrated_snn_res.4])
# 
# #plot manual curation
# png(filename = "results/images/UMAP_manual_annotation_stringent_integration_biological_replicates.png", units = "cm", width = 40, height = 40, res = 200)
# 
# DimPlot_scCustom(All.sct, group.by = "cluster_ID_manual_curation",
#                  aspect_ratio = 1, label = TRUE, label.box = TRUE, repel = TRUE,
#                  colors_use = DiscretePalette_scCustomize(num_colors = length(table(All.sct$cluster_ID_manual_curation)), palette = "ditto_seq"))
# 
# dev.off()
# 
# 


### Adding SAMap - S. mediterranea

#read the output of SAMap
tmp <- read.csv("/data/Mlig_scRNA_Seq/Francisco/data/cross_species_annotation/cluster_annotation_three_spp_DASH.tsv", header = TRUE, row.names = 1)

#remove M. lignano duplicates from columns
pattern <- "^ml_"
tmp <- tmp[, !grepl(pattern, colnames(tmp))]

#remove S. mansoni from columns
pattern <- "^sc_"
tmp <- tmp[,!grepl(pattern, colnames(tmp))]

#remove S. mediterranea duplicates from rows
pattern <- "^sm_"
tmp <- tmp[!grepl(pattern, rownames(tmp)),]

#remove S. mansoni from rows
pattern <- "^sc_"
tmp <- tmp[!grepl(pattern, rownames(tmp)),]

#remove species two-letter codes from cluster IDs
rownames(tmp) <- gsub(pattern = "ml_", replacement = "", rownames(tmp))
colnames(tmp) <- gsub(pattern = "sm_", replacement = "", colnames(tmp))

avg_features_per_cluster <- tmp

#rm(tmp)

cluster_to_max_feature <- vec <- character(length = length(rownames(avg_features_per_cluster)))

names(cluster_to_max_feature) <- sort(rownames(avg_features_per_cluster))

#iterate through features (module scores) and clusters to select, for each cluster, its ID based on average SAMap values
for (i in rownames(avg_features_per_cluster)) {
  max_val <- -10000
  cluster_ID <- "NA"
  for (i2 in colnames(avg_features_per_cluster)) {
    if (avg_features_per_cluster[i,i2] > max_val) {
      #      if (avg_features_per_cluster[i,i2] > 0) {
      if (avg_features_per_cluster[i,i2] > 0) { #cutoff
        cluster_ID <- i2
        max_val <- avg_features_per_cluster[i,i2]
      }
    }
  }
  cluster_to_max_feature[i] <- cluster_ID
  print(paste0("Cluster ID: ", i, " - annotation: ", cluster_ID, " - score: ", max_val))
}

#check if everything worked all right
#unique(seurat_integrated_final$seurat_clusters)

#create a new metadata layer with cross-species cluster IDs
seurat_integrated_final$cluster_ID_cross_species_annotation_sm_relaxed <- unname(cluster_to_max_feature[as.character(seurat_integrated_final$integrated_snn_res.1.2)])

# cell_clusters <- Idents(seurat_integrated_final)
# cell_level_data <- avg_features_per_cluster[as.character(cell_clusters),]
# colnames(cell_level_data) <- paste0("Score_cross_spp_", colnames(cell_level_data))
# 
# seurat_integrated_final <- AddMetaData(seurat_integrated_final, metadata = cell_level_data)

DimPlot_scCustom(seurat_integrated_final, group.by = "cluster_ID_cross_species_annotation_sm")



cluster_to_max_feature <- vec <- character(length = length(rownames(avg_features_per_cluster)))

names(cluster_to_max_feature) <- sort(rownames(avg_features_per_cluster))

#iterate through features (module scores) and clusters to select, for each cluster, its ID based on average SAMap values
for (i in rownames(avg_features_per_cluster)) {
  max_val <- -10000
  cluster_ID <- "NA"
  for (i2 in colnames(avg_features_per_cluster)) {
    if (avg_features_per_cluster[i,i2] > max_val) {
      #      if (avg_features_per_cluster[i,i2] > 0) {
      if (avg_features_per_cluster[i,i2] > 0.3) { #cutoff
        cluster_ID <- i2
        max_val <- avg_features_per_cluster[i,i2]
      }
    }
  }
  cluster_to_max_feature[i] <- cluster_ID
  print(paste0("Cluster ID: ", i, " - annotation: ", cluster_ID, " - score: ", max_val))
}

#check if everything worked all right
#unique(seurat_integrated_final$seurat_clusters)

#create a new metadata layer with cross-species cluster IDs
seurat_integrated_final$cluster_ID_cross_species_annotation_sm_stringent <- unname(cluster_to_max_feature[as.character(seurat_integrated_final$integrated_snn_res.1.2)])


### Adding SAMap - S. mansoni

#read the output of SAMap
tmp <- read.csv("/data/Mlig_scRNA_Seq/Francisco/data/cross_species_annotation/cluster_annotation_three_spp_DASH.tsv", header = TRUE, row.names = 1)

#remove M. lignano duplicates from columns
pattern <- "^ml_"
tmp <- tmp[, !grepl(pattern, colnames(tmp))]

#remove S. mediterranea from columns
pattern <- "^sm_"
tmp <- tmp[,!grepl(pattern, colnames(tmp))]

#remove S. mediterranea from rows
pattern <- "^sm_"
tmp <- tmp[!grepl(pattern, rownames(tmp)),]

#remove S. mansoni from rows
pattern <- "^sc_"
tmp <- tmp[!grepl(pattern, rownames(tmp)),]

#remove species two-letter codes from cluster IDs
rownames(tmp) <- gsub(pattern = "ml_", replacement = "", rownames(tmp))
colnames(tmp) <- gsub(pattern = "sc_", replacement = "", colnames(tmp))

avg_features_per_cluster <- tmp

#rm(tmp)

cluster_to_max_feature <- vec <- character(length = length(rownames(avg_features_per_cluster)))

names(cluster_to_max_feature) <- sort(rownames(avg_features_per_cluster))

#iterate through features (module scores) and clusters to select, for each cluster, its ID based on average SAMap values
for (i in rownames(avg_features_per_cluster)) {
  max_val <- -10000
  cluster_ID <- "NA"
  for (i2 in colnames(avg_features_per_cluster)) {
    if (avg_features_per_cluster[i,i2] > max_val) {
      #      if (avg_features_per_cluster[i,i2] > 0) {
      if (avg_features_per_cluster[i,i2] > 0) { #cutoff
        cluster_ID <- i2
        max_val <- avg_features_per_cluster[i,i2]
      }
    }
  }
  cluster_to_max_feature[i] <- cluster_ID
  print(paste0("Cluster ID: ", i, " - annotation: ", cluster_ID, " - score: ", max_val))
}

#check if everything worked all right
#unique(seurat_integrated_final$seurat_clusters)

#create a new metadata layer with cross-species cluster IDs
seurat_integrated_final$cluster_ID_cross_species_annotation_sc_relaxed <- unname(cluster_to_max_feature[as.character(seurat_integrated_final$integrated_snn_res.1.2)])



cluster_to_max_feature <- vec <- character(length = length(rownames(avg_features_per_cluster)))

names(cluster_to_max_feature) <- sort(rownames(avg_features_per_cluster))

#iterate through features (module scores) and clusters to select, for each cluster, its ID based on average SAMap values
for (i in rownames(avg_features_per_cluster)) {
  max_val <- -10000
  cluster_ID <- "NA"
  for (i2 in colnames(avg_features_per_cluster)) {
    if (avg_features_per_cluster[i,i2] > max_val) {
      #      if (avg_features_per_cluster[i,i2] > 0) {
      if (avg_features_per_cluster[i,i2] > 0.3) { #cutoff
        cluster_ID <- i2
        max_val <- avg_features_per_cluster[i,i2]
      }
    }
  }
  cluster_to_max_feature[i] <- cluster_ID
  print(paste0("Cluster ID: ", i, " - annotation: ", cluster_ID, " - score: ", max_val))
}

#check if everything worked all right
#unique(seurat_integrated_final$seurat_clusters)

#create a new metadata layer with cross-species cluster IDs
seurat_integrated_final$cluster_ID_cross_species_annotation_sc_stringent <- unname(cluster_to_max_feature[as.character(seurat_integrated_final$integrated_snn_res.1.2)])


# cell_clusters <- Idents(seurat_integrated_final)
# cell_level_data <- avg_features_per_cluster[as.character(cell_clusters),]
# colnames(cell_level_data) <- paste0("Score_cross_spp_", colnames(cell_level_data))
# 
# seurat_integrated_final <- AddMetaData(seurat_integrated_final, metadata = cell_level_data)

p_sc_s <- DimPlot_scCustom(seurat_integrated_final, group.by = "cluster_ID_cross_species_annotation_sc_stringent", aspect_ratio = 1, label = TRUE, repel = TRUE)
p_sc_r <- DimPlot_scCustom(seurat_integrated_final, group.by = "cluster_ID_cross_species_annotation_sc_relaxed", aspect_ratio = 1, label = TRUE, repel = TRUE, )
p_sm_s <- DimPlot_scCustom(seurat_integrated_final, group.by = "cluster_ID_cross_species_annotation_sm_stringent", aspect_ratio = 1, label = TRUE, repel = TRUE)
p_sm_r <- DimPlot_scCustom(seurat_integrated_final, group.by = "cluster_ID_cross_species_annotation_sm_relaxed", aspect_ratio = 1, label = TRUE, repel = TRUE)

pdf("figures/DimPlot_cross_spp_ann_S_mansoni_stringent.pdf", width = 10, height = 10)
p_sc_s
dev.off()

pdf("figures/DimPlot_cross_spp_ann_S_mansoni_relaxed.pdf", width = 10, height = 10)
p_sc_r
dev.off()

pdf("figures/DimPlot_cross_spp_ann_S_mediterranea_stringent.pdf", width = 10, height = 10)
p_sm_s
dev.off()

pdf("figures/DimPlot_cross_spp_ann_S_mediterranea_relaxed.pdf", width = 10, height = 10)
p_sm_r
dev.off()

colnames(seurat_integrated_final@meta.data)

saveRDS(seurat_integrated_final, file = "../results/objs/seurat_integrated_final_plus_module_scores_and_cross_spp_annots.rds")


### Comparing DASH versus non-DASH

# SCT
DefaultAssay(seurat_integrated_final) <- "SCT"
seurat_integrated_final <- PrepSCTFindMarkers(seurat_integrated_final)
DEGS_Ad_non_DASH_vs_DASH_SCT <- FindMarkers(seurat_integrated_final, ident.1 = "2_Adult_Unsorted", ident.2 = "2_Adult_Unsorted_DASH")


# RNA
DefaultAssay(seurat_integrated_final) <- "RNA"
seurat_integrated_final <- NormalizeData(seurat_integrated_final, normalization.method = "LogNormalize", scale.factor = 10000)
DEGS_Ad_non_DASH_vs_DASH_RNA <- FindMarkers(seurat_integrated_final, ident.1 = "2_Adult_Unsorted", ident.2 = "2_Adult_Unsorted_DASH")




## Generating GO enrichment analysis for Sam's talk

library(topGO)

seurat_integrated_final <- PrepSCTFindMarkers(seurat_integrated_final)

Idents(seurat_integrated_final) <- "integrated_snn_res.1.2"

marker.genes <- FindAllMarkers(seurat_integrated_final, assay = "SCT")

all.genes <- rownames(seurat_integrated_final@assays$SCT)

geneID2GOs <- read.table("/data/Mlig_scRNA_Seq/Francisco/data/geneid2go.map", sep = "\t", header = FALSE)
gene2GO <- strsplit(geneID2GOs$V2, ", ",fixed=T)
names(gene2GO) <- geneID2GOs$V1

results_enrichment <- data.frame(GO.ID=character(),
                                 Term=character(),
                                 Annotated=integer(),
                                 Significant=integer(),
                                 Expected =double(),
                                 Fisher=double(),
                                 Fisher_FDR=double(),
                                 ClusterID=character()
)

for (cluster in unique(Idents(seurat_integrated_final))) {
  print(cluster)
  genes <- marker.genes$gene[marker.genes$cluster == cluster & marker.genes$p_val_adj < 0.00001 & marker.genes$avg_log2FC > 2]
  # define geneList as 1 if gene is in expressed.genes, 0 otherwise
  geneList <- ifelse(all.genes %in% genes, 1, 0)
  names(geneList) <- all.genes
  # Create topGOdata object
  GOdata <- new("topGOdata",
                ontology = "BP",
                allGenes = geneList,
                geneSelectionFun = function(x)(x == 1),
                annot = annFUN.gene2GO,
                nodeSize = 3,
                gene2GO = gene2GO)
  
  # Test for enrichment using Fisher's Exact Test
  resultFisher <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
  results <- GenTable(GOdata, Fisher = resultFisher, topNodes = length(GOdata@graph@nodes), numChar = 1000)
  results$Fisher <- ifelse(results$Fisher == "< 1e-30", 1e-30, as.numeric(results$Fisher))
  results$Fisher_FDR <- p.adjust(results$Fisher)
  results$ClusterID <- cluster
  print(results[results$Fisher_FDR < 0.2,])
  results_enrichment <- rbind(results_enrichment, results)
}

openxlsx::write.xlsx(x = results_enrichment, file = "../results/GO_enrichment_analysis_Sam_talk.xlsx")

showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 20, useInfo = 'all')

