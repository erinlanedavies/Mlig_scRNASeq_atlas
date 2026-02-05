## 00_QC_analysis.R
#Author: Francisco Lobo, 04/2024

# Takes as input the outputs of Starsolo and creates Seurat objects

#installing libraries
#BiocManager::install("speckle")

#loading libraries
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(dplyr)
library(harmony)
library(scCustomize)
library(speckle)
library(scater)
library(DoubletFinder)


##clean workspace
rm(list=ls())

#set basedir
basedir <- "/data/Mlig_scRNA_Seq/Francisco/data/STARsolo_output_Erin/"
setwd(basedir)

## Read Cellbender output (matrix of gene expression per cell)

#0-1 day

results_path <- paste0(basedir, "17_0_1d_Hatchlings_learning_rate/Solo.out/Gene/cellbender_filtered.h5")
S17_0_1d_Hatchlings_counts <- Read_CellBender_h5_Mat(file_name = results_path)

results_path <- paste0(basedir, "18_0_1d_Hatchling_learning_rate/Solo.out/Gene/cellbender_filtered.h5")
S18_0_1d_Hatchling_counts <- Read_CellBender_h5_Mat(file_name = results_path)

results_path <- paste0(basedir, "19_0_1d_Hatchling_learning_rate/Solo.out/Gene/cellbender_filtered.h5")
S19_0_1d_Hatchling_counts <- Read_CellBender_h5_Mat(file_name = results_path)


#3 days
results_path <- paste0(basedir, "11_unsorted_3dph_learning_rate/Solo.out/Gene/cellbender_filtered.h5")
S11_unsorted_3dph_counts <- Read_CellBender_h5_Mat(file_name = results_path)

results_path <- paste0(basedir, "12_unsorted_3dph_learning_rate/Solo.out/Gene/cellbender_filtered.h5")
S12_unsorted_3dph_counts <- Read_CellBender_h5_Mat(file_name = results_path)

results_path <- paste0(basedir, "16_3d_juveniles_learning_rate/Solo.out/Gene/cellbender_filtered.h5")
S16_3d_juveniles_counts <- Read_CellBender_h5_Mat(file_name = results_path)


#7 days
results_path <- paste0(basedir, "8_unsorted_Hatchling_7d/Solo.out/Gene/cellbender_filtered.h5")
S8_unsorted_Hatchling_7d_counts <- Read_CellBender_h5_Mat(file_name = results_path)

results_path <- paste0(basedir, "15_7d_juveniles/Solo.out/Gene/cellbender_filtered.h5")
S15_7d_juveniles_counts <- Read_CellBender_h5_Mat(file_name = results_path)

results_path <- paste0(basedir, "9_unsorted_Hatchling_7d_learning_rate/Solo.out/Gene/cellbender_filtered.h5")
S9_unsorted_Hatchling_7d_counts <- Read_CellBender_h5_Mat(file_name = results_path)


#Adult
#results_path <- paste0(basedir, "1_Adult_Unsorted/Solo.out/Gene/cellbender_filtered.h5")
#S1_Adult_Unsorted_counts <- Read_CellBender_h5_Mat(file_name = results_path)

results_path <- paste0(basedir, "2_Adult_Unsorted/Solo.out/Gene/cellbender_filtered.h5")
S2_Adult_Unsorted_counts <- Read_CellBender_h5_Mat(file_name = results_path)

results_path <- paste0(basedir, "2_Adult_Unsorted_DASH/Solo.out/Gene/cellbender_filtered.h5")
S2_Adult_Unsorted_DASH_counts <- Read_CellBender_h5_Mat(file_name = results_path)

results_path <- paste0(basedir, "7_Sorted_Adult/Solo.out/Gene/cellbender_filtered.h5")
S7_Sorted_Adult_counts <- Read_CellBender_h5_Mat(file_name = results_path)


######################################

## Create Seurat objects

## 0-1 day

#From counts to objs
S17_0_1d_Hatchlings <- CreateSeuratObject(counts = S17_0_1d_Hatchlings_counts, min.cells=3, min.features = 200, project = "01D_sample1")
S18_0_1d_Hatchling <- CreateSeuratObject(counts = S18_0_1d_Hatchling_counts, min.cells=3, min.features = 200, project = "01D_sample2")
S19_0_1d_Hatchling <- CreateSeuratObject(counts = S19_0_1d_Hatchling_counts, min.cells=3, min.features = 200, project = "01D_sample3")

#metadata
S17_0_1d_Hatchlings$sample_ID <- "17_0_1d_Hatchlings"
S18_0_1d_Hatchling$sample_ID <- "18_0_1d_Hatchling"
S19_0_1d_Hatchling$sample_ID <- "19_0_1d_Hatchling"

S17_0_1d_Hatchlings$sample_type <- "01D"
S18_0_1d_Hatchling$sample_type <- "01D"
S19_0_1d_Hatchling$sample_type <- "01D"

S17_0_1d_Hatchlings$origin <- "Erin"
S18_0_1d_Hatchling$origin <- "Erin"
S19_0_1d_Hatchling$origin <- "Erin"


## 3 days

#From counts to objs
S11_unsorted_3dph <- CreateSeuratObject(counts = S11_unsorted_3dph_counts, min.cells=3, min.features = 200, project = "3D_sample1")
S12_unsorted_3dph <- CreateSeuratObject(counts = S12_unsorted_3dph_counts, min.cells=3, min.features = 200, project = "3D_sample2")
S16_3d_juveniles <- CreateSeuratObject(counts = S16_3d_juveniles_counts, min.cells=3, min.features = 200, project = "3D_sample3")

#metadata
S11_unsorted_3dph$sample_ID <- "11_unsorted_3dph"
S12_unsorted_3dph$sample_ID <- "12_unsorted_3dph"
S16_3d_juveniles$sample_ID <- "16_3d_juveniles"

S11_unsorted_3dph$sample_type <- "3D"
S12_unsorted_3dph$sample_type <- "3D"
S16_3d_juveniles$sample_type <- "3D"

S11_unsorted_3dph$origin <- "Erin"
S12_unsorted_3dph$origin <- "Erin"
S16_3d_juveniles$origin <- "Erin"


## 7 days

results_path <- paste0(basedir, "8_unsorted_Hatchling_7d/Solo.out/Gene/cellbender_filtered.h5")
S8_unsorted_Hatchling_7d_counts <- Read_CellBender_h5_Mat(file_name = results_path)

results_path <- paste0(basedir, "15_7d_juveniles/Solo.out/Gene/cellbender_filtered.h5")
S15_7d_juveniles_counts <- Read_CellBender_h5_Mat(file_name = results_path)

results_path <- paste0(basedir, "9_unsorted_Hatchling_7d_learning_rate/Solo.out/Gene/cellbender_filtered.h5")
S9_unsorted_Hatchling_7d_counts <- Read_CellBender_h5_Mat(file_name = results_path)


#From counts to objs
S8_unsorted_Hatchling_7d <- CreateSeuratObject(counts = S8_unsorted_Hatchling_7d_counts, min.cells=3, min.features = 200, project = "7D_sample1")
S15_7d_juveniles <- CreateSeuratObject(counts = S15_7d_juveniles_counts, min.cells=3, min.features = 200, project = "7D_sample2")
S9_unsorted_Hatchling_7d <- CreateSeuratObject(counts = S9_unsorted_Hatchling_7d_counts, min.cells=3, min.features = 200, project = "7D_sample3")

#metadata
S8_unsorted_Hatchling_7d$sample_ID <- "8_unsorted_Hatchling_7d"
S15_7d_juveniles$sample_ID <- "15_7d_juveniles"
S9_unsorted_Hatchling_7d$sample_ID <- "9_unsorted_Hatchling_7d"

S8_unsorted_Hatchling_7d$sample_type <- "7D"
S15_7d_juveniles$sample_type <- "7D"
S9_unsorted_Hatchling_7d$sample_type <- "7D"

S8_unsorted_Hatchling_7d$origin <- "Erin"
S15_7d_juveniles$origin <- "Erin"
S9_unsorted_Hatchling_7d$origin <- "Erin"

## Adult

#From counts to objs
S1_Adult_Unsorted <- CreateSeuratObject(counts = S1_Adult_Unsorted_counts, min.cells=3, min.features = 200, project = "Ad_sample1")
S2_Adult_Unsorted <- CreateSeuratObject(counts = S2_Adult_Unsorted_counts, min.cells=3, min.features = 200, project = "Ad_sample2")
S2_Adult_Unsorted_DASH <- CreateSeuratObject(counts = S2_Adult_Unsorted_DASH_counts, min.cells=3, min.features = 200, project = "Ad_sample3")
S7_Sorted_Adult <- CreateSeuratObject(counts = S7_Sorted_Adult_counts, min.cells=3, min.features = 200, project = "Ad_sample1noDASH")

#metadata
S1_Adult_Unsorted$sample_ID <- "1_Adult_Unsorted"
S2_Adult_Unsorted$sample_ID <- "2_Adult_Unsorted"
S2_Adult_Unsorted_DASH$sample_ID <- "2_Adult_Unsorted_DASH"
S7_Sorted_Adult$sample_ID <- "7_Sorted_Adult"

S1_Adult_Unsorted$sample_type <- "Ad"
S2_Adult_Unsorted$sample_type <- "Ad"
S2_Adult_Unsorted_DASH$sample_type <- "Ad"
S7_Sorted_Adult$sample_type <- "Ad"

S1_Adult_Unsorted$origin <- "Erin"
S2_Adult_Unsorted$origin <- "Erin"
S2_Adult_Unsorted_DASH$origin <- "Erin"
S7_Sorted_Adult$origin <- "Erin"


#Computing the %mitochondrial RNA

# /data/Mlig_scRNA_Seq/projectOutputs/mligDavies/STARsolo/Mlig_4_5.STAR_2_7_9a/

# STAR --runThreadN 8 --runMode genomeGenerate --genomeFastaFiles /data/paezbaenal2/resources/mligDavies/fastas/Mlig_4_5.fa --sjdbGTFfile /data/paezbaenal2/resources/mligDavies/gtfs/merged_mito_rRNA_Mlig_RNA_4_5_v5.coregenes.gtf --genomeDir /data/paezbaenal2/projectOutputs/mligDavies/STARsolo/Mlig_4_5.STAR_2_7_9a

mito_prot_coding_genes <- c("Mlig999001.1",
                            "Mlig999002.1",
                            "Mlig999012.1",
                            "Mlig999013.1",
                            "Mlig999014.1",
                            "Mlig999015.1",
                            "Mlig999016.1",
                            "Mlig999020.1",
                            "Mlig999023.1",
                            "Mlig999025.1",
                            "Mlig999027.1",
                            "Mlig999031.1",
                            "Mlig999036.1")

S17_0_1d_Hatchlings[["percent.mito"]] <- PercentageFeatureSet(S17_0_1d_Hatchlings, features = mito_prot_coding_genes)
S18_0_1d_Hatchling[["percent.mito"]] <- PercentageFeatureSet(S18_0_1d_Hatchling, features = mito_prot_coding_genes)
S19_0_1d_Hatchling[["percent.mito"]] <- PercentageFeatureSet(S19_0_1d_Hatchling, features = mito_prot_coding_genes)
S11_unsorted_3dph[["percent.mito"]] <- PercentageFeatureSet(S11_unsorted_3dph, features = mito_prot_coding_genes)
S12_unsorted_3dph[["percent.mito"]] <- PercentageFeatureSet(S12_unsorted_3dph, features = mito_prot_coding_genes)
S16_3d_juveniles[["percent.mito"]] <- PercentageFeatureSet(S16_3d_juveniles, features = mito_prot_coding_genes)
S8_unsorted_Hatchling_7d[["percent.mito"]] <- PercentageFeatureSet(S8_unsorted_Hatchling_7d, features = mito_prot_coding_genes)
S15_7d_juveniles[["percent.mito"]] <- PercentageFeatureSet(S15_7d_juveniles, features = mito_prot_coding_genes)
S9_unsorted_Hatchling_7d[["percent.mito"]] <- PercentageFeatureSet(S9_unsorted_Hatchling_7d, features = mito_prot_coding_genes)
S1_Adult_Unsorted[["percent.mito"]] <- PercentageFeatureSet(S1_Adult_Unsorted, features = mito_prot_coding_genes)
S2_Adult_Unsorted[["percent.mito"]] <- PercentageFeatureSet(S2_Adult_Unsorted, features = mito_prot_coding_genes)
S2_Adult_Unsorted_DASH[["percent.mito"]] <- PercentageFeatureSet(S2_Adult_Unsorted_DASH, features = mito_prot_coding_genes)
S7_Sorted_Adult[["percent.mito"]] <- PercentageFeatureSet(S7_Sorted_Adult, features = mito_prot_coding_genes)


#filtering & doublet detection

#S17_0_1d_Hatchlings
#filtering
qc.lib <- scater::isOutlier(S17_0_1d_Hatchlings$nCount_RNA, nmads=2, log=TRUE, type="both")
qc.nexprs <- isOutlier(S17_0_1d_Hatchlings$nFeature_RNA, log=TRUE, nmads=2, type="both")
qc.mito <- isOutlier(S17_0_1d_Hatchlings$percent.mito, nmads=1, type="higher")
discard <- qc.lib | qc.nexprs  | qc.mito
S17_0_1d_Hatchlings$discard<-discard
S17_0_1d_Hatchlings_filt <-subset(S17_0_1d_Hatchlings, subset = discard!=TRUE)

#doublet detection
tmp <- SCTransform(S17_0_1d_Hatchlings_filt, vars.to.regress = "percent.mito")
tmp <- RunPCA(tmp)
tmp <- RunUMAP(tmp, dims = 1:40)
tmp <- FindNeighbors(tmp, reduction="pca", dims = 1:40)
tmp <- FindClusters(object = tmp, verbose = FALSE)

sweep.res.list <- paramSweep(tmp, PCs = 1:40, sct = TRUE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)

pK=as.numeric(as.character(bcmvn$pK))
BCmetric=bcmvn$BCmetric
pK_choose = pK[which(BCmetric %in% max(BCmetric))]

par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
plot(x = pK, y = BCmetric, pch = 16,type="b",
     col = "blue",lty=1)
abline(v=pK_choose,lwd=2,col='red',lty=2)
title("The BCmvn distributions")
text(pK_choose,max(BCmetric),as.character(pK_choose),pos = 4,col = "red")

nExp_poi <- round(0.075*nrow(tmp@meta.data))

S17_0_1d_Hatchlings_filt <- doubletFinder(tmp, PCs = 1:40, pN = 0.25, pK = pK_choose, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)

saveRDS(object = S17_0_1d_Hatchlings_filt, file = "/data/Mlig_scRNA_Seq/Francisco/results/objs/S17_0_1d_Hatchlings_filt.rds")


#S18_0_1d_Hatchling
#filtering
qc.lib <- scater::isOutlier(S18_0_1d_Hatchling$nCount_RNA, nmads=2, log=TRUE, type="both")
qc.nexprs <- isOutlier(S18_0_1d_Hatchling$nFeature_RNA, log=TRUE, nmads=2, type="both")
qc.mito <- isOutlier(S18_0_1d_Hatchling$percent.mito, nmads=1, type="higher")
discard <- qc.lib | qc.nexprs  | qc.mito
S18_0_1d_Hatchling$discard<-discard
S18_0_1d_Hatchling_filt <-subset(S18_0_1d_Hatchling, subset= discard!=TRUE)

#doublet detection
tmp <- SCTransform(S18_0_1d_Hatchling_filt, vars.to.regress = "percent.mito")
tmp <- RunPCA(tmp)
tmp <- RunUMAP(tmp, dims = 1:40)
tmp <- FindNeighbors(tmp, reduction="pca", dims = 1:40)
tmp <- FindClusters(object = tmp, verbose = FALSE)

sweep.res.list <- paramSweep(tmp, PCs = 1:40, sct = TRUE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)

pK=as.numeric(as.character(bcmvn$pK))
BCmetric=bcmvn$BCmetric
pK_choose = pK[which(BCmetric %in% max(BCmetric))]

par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
plot(x = pK, y = BCmetric, pch = 16,type="b",
     col = "blue",lty=1)
abline(v=pK_choose,lwd=2,col='red',lty=2)
title("The BCmvn distributions")
text(pK_choose,max(BCmetric),as.character(pK_choose),pos = 4,col = "red")

nExp_poi <- round(0.075*nrow(tmp@meta.data))

S18_0_1d_Hatchling_filt <- doubletFinder(tmp, PCs = 1:40, pN = 0.25, pK = pK_choose, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)

saveRDS(object = S18_0_1d_Hatchling_filt, file = "/data/Mlig_scRNA_Seq/Francisco/results/objs/S18_0_1d_Hatchling_filt.rds")


#S19_0_1d_Hatchling
#filtering
qc.lib <- scater::isOutlier(S19_0_1d_Hatchling$nCount_RNA, nmads=2, log=TRUE, type="both")
qc.nexprs <- isOutlier(S19_0_1d_Hatchling$nFeature_RNA, log=TRUE, nmads=2, type="both")
qc.mito <- isOutlier(S19_0_1d_Hatchling$percent.mito, nmads=1, type="higher")
discard <- qc.lib | qc.nexprs  | qc.mito
S19_0_1d_Hatchling$discard<-discard
S19_0_1d_Hatchling_filt <-subset(S19_0_1d_Hatchling, subset= discard!=TRUE)

#doublet detection
tmp <- SCTransform(S19_0_1d_Hatchling_filt, vars.to.regress = "percent.mito")
tmp <- RunPCA(tmp)
tmp <- RunUMAP(tmp, dims = 1:40)
tmp <- FindNeighbors(tmp, reduction="pca", dims = 1:40)
tmp <- FindClusters(object = tmp, verbose = FALSE)

sweep.res.list <- paramSweep(tmp, PCs = 1:40, sct = TRUE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)

pK=as.numeric(as.character(bcmvn$pK))
BCmetric=bcmvn$BCmetric
pK_choose = pK[which(BCmetric %in% max(BCmetric))]

par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
plot(x = pK, y = BCmetric, pch = 16,type="b",
     col = "blue",lty=1)
abline(v=pK_choose,lwd=2,col='red',lty=2)
title("The BCmvn distributions")
text(pK_choose,max(BCmetric),as.character(pK_choose),pos = 4,col = "red")

nExp_poi <- round(0.075*nrow(tmp@meta.data))

S19_0_1d_Hatchling_filt <- doubletFinder(tmp, PCs = 1:40, pN = 0.25, pK = pK_choose, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)

saveRDS(object = S19_0_1d_Hatchling_filt, file = "/data/Mlig_scRNA_Seq/Francisco/results/objs/S19_0_1d_Hatchling_filt.rds")


#S11_unsorted_3dph
#filtering
qc.lib <- scater::isOutlier(S11_unsorted_3dph$nCount_RNA, nmads=2, log=TRUE, type="both")
qc.nexprs <- isOutlier(S11_unsorted_3dph$nFeature_RNA, log=TRUE, nmads=2, type="both")
qc.mito <- isOutlier(S11_unsorted_3dph$percent.mito, nmads=1, type="higher")
discard <- qc.lib | qc.nexprs  | qc.mito
S11_unsorted_3dph$discard<-discard
S11_unsorted_3dph_filt <-subset(S11_unsorted_3dph, subset= discard!=TRUE)

#doublet detection
tmp <- SCTransform(S11_unsorted_3dph_filt, vars.to.regress = "percent.mito")
tmp <- RunPCA(tmp)
tmp <- RunUMAP(tmp, dims = 1:40)
tmp <- FindNeighbors(tmp, reduction="pca", dims = 1:40)
tmp <- FindClusters(object = tmp, verbose = FALSE)

sweep.res.list <- paramSweep(tmp, PCs = 1:40, sct = TRUE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)

pK=as.numeric(as.character(bcmvn$pK))
BCmetric=bcmvn$BCmetric
pK_choose = pK[which(BCmetric %in% max(BCmetric))]

par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
plot(x = pK, y = BCmetric, pch = 16,type="b",
     col = "blue",lty=1)
abline(v=pK_choose,lwd=2,col='red',lty=2)
title("The BCmvn distributions")
text(pK_choose,max(BCmetric),as.character(pK_choose),pos = 4,col = "red")

nExp_poi <- round(0.075*nrow(tmp@meta.data))

S11_unsorted_3dph_filt <- doubletFinder(tmp, PCs = 1:40, pN = 0.25, pK = pK_choose, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)

saveRDS(object = S11_unsorted_3dph_filt, file = "/data/Mlig_scRNA_Seq/Francisco/results/objs/S11_unsorted_3dph_filt.rds")


#S12_unsorted_3dph
#filtering
qc.lib <- scater::isOutlier(S12_unsorted_3dph$nCount_RNA, nmads=2, log=TRUE, type="both")
qc.nexprs <- isOutlier(S12_unsorted_3dph$nFeature_RNA, log=TRUE, nmads=2, type="both")
qc.mito <- isOutlier(S12_unsorted_3dph$percent.mito, nmads=1, type="higher")
discard <- qc.lib | qc.nexprs  | qc.mito
S12_unsorted_3dph$discard<-discard
S12_unsorted_3dph_filt <-subset(S12_unsorted_3dph, subset= discard!=TRUE)

#doublet detection
tmp <- SCTransform(S12_unsorted_3dph_filt, vars.to.regress = "percent.mito")
tmp <- RunPCA(tmp)
tmp <- RunUMAP(tmp, dims = 1:40)
tmp <- FindNeighbors(tmp, reduction="pca", dims = 1:40)
tmp <- FindClusters(object = tmp, verbose = FALSE)

sweep.res.list <- paramSweep(tmp, PCs = 1:40, sct = TRUE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)

pK=as.numeric(as.character(bcmvn$pK))
BCmetric=bcmvn$BCmetric
pK_choose = pK[which(BCmetric %in% max(BCmetric))]

par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
plot(x = pK, y = BCmetric, pch = 16,type="b",
     col = "blue",lty=1)
abline(v=pK_choose,lwd=2,col='red',lty=2)
title("The BCmvn distributions")
text(pK_choose,max(BCmetric),as.character(pK_choose),pos = 4,col = "red")

nExp_poi <- round(0.075*nrow(tmp@meta.data))

S12_unsorted_3dph_filt <- doubletFinder(tmp, PCs = 1:40, pN = 0.25, pK = pK_choose, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)

saveRDS(object = S12_unsorted_3dph_filt, file = "/data/Mlig_scRNA_Seq/Francisco/results/objs/S12_unsorted_3dph_filt.rds")


#S16_3d_juveniles
#filtering
qc.lib <- scater::isOutlier(S16_3d_juveniles$nCount_RNA, nmads=2, log=TRUE, type="both")
qc.nexprs <- isOutlier(S16_3d_juveniles$nFeature_RNA, log=TRUE, nmads=2, type="both")
qc.mito <- isOutlier(S16_3d_juveniles$percent.mito, nmads=1, type="higher")
discard <- qc.lib | qc.nexprs  | qc.mito
S16_3d_juveniles$discard<-discard
S16_3d_juveniles_filt <-subset(S16_3d_juveniles, subset= discard!=TRUE)

#doublet detection
tmp <- SCTransform(S16_3d_juveniles_filt, vars.to.regress = "percent.mito")
tmp <- RunPCA(tmp)
tmp <- RunUMAP(tmp, dims = 1:40)
tmp <- FindNeighbors(tmp, reduction="pca", dims = 1:40)
tmp <- FindClusters(object = tmp, verbose = FALSE)

sweep.res.list <- paramSweep(tmp, PCs = 1:40, sct = TRUE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)

pK=as.numeric(as.character(bcmvn$pK))
BCmetric=bcmvn$BCmetric
pK_choose = pK[which(BCmetric %in% max(BCmetric))]

par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
plot(x = pK, y = BCmetric, pch = 16,type="b",
     col = "blue",lty=1)
abline(v=pK_choose,lwd=2,col='red',lty=2)
title("The BCmvn distributions")
text(pK_choose,max(BCmetric),as.character(pK_choose),pos = 4,col = "red")

nExp_poi <- round(0.075*nrow(tmp@meta.data))

S16_3d_juveniles_filt <- doubletFinder(tmp, PCs = 1:40, pN = 0.25, pK = pK_choose, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)

saveRDS(object = S16_3d_juveniles_filt, file = "/data/Mlig_scRNA_Seq/Francisco/results/objs/S16_3d_juveniles_filt.rds")


#S8_unsorted_Hatchling_7d
#filtering
qc.lib <- scater::isOutlier(S8_unsorted_Hatchling_7d$nCount_RNA, nmads=2, log=TRUE, type="both")
qc.nexprs <- isOutlier(S8_unsorted_Hatchling_7d$nFeature_RNA, log=TRUE, nmads=2, type="both")
qc.mito <- isOutlier(S8_unsorted_Hatchling_7d$percent.mito, nmads=1, type="higher")
discard <- qc.lib | qc.nexprs  | qc.mito
S8_unsorted_Hatchling_7d$discard<-discard
S8_unsorted_Hatchling_7d_filt <-subset(S8_unsorted_Hatchling_7d, subset= discard!=TRUE)

#doublet detection
tmp <- SCTransform(S8_unsorted_Hatchling_7d_filt, vars.to.regress = "percent.mito")
tmp <- RunPCA(tmp)
tmp <- RunUMAP(tmp, dims = 1:40)
tmp <- FindNeighbors(tmp, reduction="pca", dims = 1:40)
tmp <- FindClusters(object = tmp, verbose = FALSE)

sweep.res.list <- paramSweep(tmp, PCs = 1:40, sct = TRUE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)

pK=as.numeric(as.character(bcmvn$pK))
BCmetric=bcmvn$BCmetric
pK_choose = pK[which(BCmetric %in% max(BCmetric))]

par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
plot(x = pK, y = BCmetric, pch = 16,type="b",
     col = "blue",lty=1)
abline(v=pK_choose,lwd=2,col='red',lty=2)
title("The BCmvn distributions")
text(pK_choose,max(BCmetric),as.character(pK_choose),pos = 4,col = "red")

nExp_poi <- round(0.075*nrow(tmp@meta.data))

S8_unsorted_Hatchling_7d_filt <- doubletFinder(tmp, PCs = 1:40, pN = 0.25, pK = pK_choose, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)

saveRDS(object = S8_unsorted_Hatchling_7d_filt, file = "/data/Mlig_scRNA_Seq/Francisco/results/objs/S8_unsorted_Hatchling_7d_filt.rds")


#S15_7d_juveniles
#filtering
qc.lib <- scater::isOutlier(S15_7d_juveniles$nCount_RNA, nmads=2, log=TRUE, type="both")
qc.nexprs <- isOutlier(S15_7d_juveniles$nFeature_RNA, log=TRUE, nmads=2, type="both")
qc.mito <- isOutlier(S15_7d_juveniles$percent.mito, nmads=1, type="higher")
discard <- qc.lib | qc.nexprs  | qc.mito
S15_7d_juveniles$discard<-discard
S15_7d_juveniles_filt <-subset(S15_7d_juveniles, subset= discard!=TRUE)

#doublet detection
tmp <- SCTransform(S15_7d_juveniles_filt, vars.to.regress = "percent.mito")
tmp <- RunPCA(tmp)
tmp <- RunUMAP(tmp, dims = 1:40)
tmp <- FindNeighbors(tmp, reduction="pca", dims = 1:40)
tmp <- FindClusters(object = tmp, verbose = FALSE)

sweep.res.list <- paramSweep(tmp, PCs = 1:40, sct = TRUE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)

pK=as.numeric(as.character(bcmvn$pK))
BCmetric=bcmvn$BCmetric
pK_choose = pK[which(BCmetric %in% max(BCmetric))]

par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
plot(x = pK, y = BCmetric, pch = 16,type="b",
     col = "blue",lty=1)
abline(v=pK_choose,lwd=2,col='red',lty=2)
title("The BCmvn distributions")
text(pK_choose,max(BCmetric),as.character(pK_choose),pos = 4,col = "red")

nExp_poi <- round(0.075*nrow(tmp@meta.data))

S15_7d_juveniles_filt <- doubletFinder(tmp, PCs = 1:40, pN = 0.25, pK = pK_choose, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)

saveRDS(object = S15_7d_juveniles_filt, file = "/data/Mlig_scRNA_Seq/Francisco/results/objs/S15_7d_juveniles_filt.rds")


#S9_unsorted_Hatchling_7d
#filtering
qc.lib <- scater::isOutlier(S9_unsorted_Hatchling_7d$nCount_RNA, nmads=2, log=TRUE, type="both")
qc.nexprs <- isOutlier(S9_unsorted_Hatchling_7d$nFeature_RNA, log=TRUE, nmads=2, type="both")
qc.mito <- isOutlier(S9_unsorted_Hatchling_7d$percent.mito, nmads=1, type="higher")
discard <- qc.lib | qc.nexprs  | qc.mito
S9_unsorted_Hatchling_7d$discard<-discard
S9_unsorted_Hatchling_7d_filt <-subset(S9_unsorted_Hatchling_7d, subset= discard!=TRUE)

#doublet detection
tmp <- SCTransform(S9_unsorted_Hatchling_7d_filt, vars.to.regress = "percent.mito")
tmp <- RunPCA(tmp)
tmp <- RunUMAP(tmp, dims = 1:40)
tmp <- FindNeighbors(tmp, reduction="pca", dims = 1:40)
tmp <- FindClusters(object = tmp, verbose = FALSE)

sweep.res.list <- paramSweep(tmp, PCs = 1:40, sct = TRUE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)

pK=as.numeric(as.character(bcmvn$pK))
BCmetric=bcmvn$BCmetric
pK_choose = pK[which(BCmetric %in% max(BCmetric))]

par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
plot(x = pK, y = BCmetric, pch = 16,type="b",
     col = "blue",lty=1)
abline(v=pK_choose,lwd=2,col='red',lty=2)
title("The BCmvn distributions")
text(pK_choose,max(BCmetric),as.character(pK_choose),pos = 4,col = "red")

nExp_poi <- round(0.075*nrow(tmp@meta.data))

S9_unsorted_Hatchling_7d_filt <- doubletFinder(tmp, PCs = 1:40, pN = 0.25, pK = pK_choose, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)

saveRDS(object = S9_unsorted_Hatchling_7d_filt, file = "/data/Mlig_scRNA_Seq/Francisco/results/objs/S9_unsorted_Hatchling_7d_filt.rds")


#S1_Adult_Unsorted
#filtering
qc.lib <- scater::isOutlier(S1_Adult_Unsorted$nCount_RNA, nmads=2, log=TRUE, type="both")
qc.nexprs <- isOutlier(S1_Adult_Unsorted$nFeature_RNA, log=TRUE, nmads=2, type="both")
qc.mito <- isOutlier(S1_Adult_Unsorted$percent.mito, nmads=1, type="higher")
discard <- qc.lib | qc.nexprs  | qc.mito
S1_Adult_Unsorted$discard<-discard
S1_Adult_Unsorted_filt <-subset(S1_Adult_Unsorted, subset= discard!=TRUE)

#doublet detection
tmp <- SCTransform(S1_Adult_Unsorted_filt, vars.to.regress = "percent.mito")
tmp <- RunPCA(tmp)
tmp <- RunUMAP(tmp, dims = 1:40)
tmp <- FindNeighbors(tmp, reduction="pca", dims = 1:40)
tmp <- FindClusters(object = tmp, verbose = FALSE)

sweep.res.list <- paramSweep(tmp, PCs = 1:40, sct = TRUE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)

pK=as.numeric(as.character(bcmvn$pK))
BCmetric=bcmvn$BCmetric
pK_choose = pK[which(BCmetric %in% max(BCmetric))]

par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
plot(x = pK, y = BCmetric, pch = 16,type="b",
     col = "blue",lty=1)
abline(v=pK_choose,lwd=2,col='red',lty=2)
title("The BCmvn distributions")
text(pK_choose,max(BCmetric),as.character(pK_choose),pos = 4,col = "red")

nExp_poi <- round(0.075*nrow(tmp@meta.data))

S1_Adult_Unsorted_filt <- doubletFinder(tmp, PCs = 1:40, pN = 0.25, pK = pK_choose, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)

saveRDS(object = S1_Adult_Unsorted_filt, file = "/data/Mlig_scRNA_Seq/Francisco/results/objs/S1_Adult_Unsorted_filt.rds")


#S2_Adult_Unsorted_DASH
#filtering
qc.lib <- scater::isOutlier(S2_Adult_Unsorted_DASH$nCount_RNA, nmads=2, log=TRUE, type="both")
qc.nexprs <- isOutlier(S2_Adult_Unsorted_DASH$nFeature_RNA, log=TRUE, nmads=2, type="both")
qc.mito <- isOutlier(S2_Adult_Unsorted_DASH$percent.mito, nmads=1, type="higher")
discard <- qc.lib | qc.nexprs  | qc.mito
S2_Adult_Unsorted_DASH$discard<-discard
S2_Adult_Unsorted_DASH_filt <-subset(S2_Adult_Unsorted_DASH, subset= discard!=TRUE)

#doublet detection
tmp <- SCTransform(S2_Adult_Unsorted_DASH_filt, vars.to.regress = "percent.mito")
tmp <- RunPCA(tmp)
tmp <- RunUMAP(tmp, dims = 1:40)
tmp <- FindNeighbors(tmp, reduction="pca", dims = 1:40)
tmp <- FindClusters(object = tmp, verbose = FALSE)

sweep.res.list <- paramSweep(tmp, PCs = 1:40, sct = TRUE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)

pK=as.numeric(as.character(bcmvn$pK))
BCmetric=bcmvn$BCmetric
pK_choose = pK[which(BCmetric %in% max(BCmetric))]

par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
plot(x = pK, y = BCmetric, pch = 16,type="b",
     col = "blue",lty=1)
abline(v=pK_choose,lwd=2,col='red',lty=2)
title("The BCmvn distributions")
text(pK_choose,max(BCmetric),as.character(pK_choose),pos = 4,col = "red")

nExp_poi <- round(0.075*nrow(tmp@meta.data))

S2_Adult_Unsorted_DASH_filt <- doubletFinder(tmp, PCs = 1:40, pN = 0.25, pK = pK_choose, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)

saveRDS(object = S2_Adult_Unsorted_DASH_filt, file = "/data/Mlig_scRNA_Seq/Francisco/results/objs/S2_Adult_Unsorted_DASH_filt.rds")



#S2_Adult_Unsorted
#filtering
qc.lib <- scater::isOutlier(S2_Adult_Unsorted$nCount_RNA, nmads=2, log=TRUE, type="both")
qc.nexprs <- isOutlier(S2_Adult_Unsorted$nFeature_RNA, log=TRUE, nmads=2, type="both")
qc.mito <- isOutlier(S2_Adult_Unsorted$percent.mito, nmads=1, type="higher")
discard <- qc.lib | qc.nexprs  | qc.mito
S2_Adult_Unsorted$discard<-discard
S2_Adult_Unsorted_filt <-subset(S2_Adult_Unsorted, subset= discard!=TRUE)

#doublet detection
tmp <- SCTransform(S2_Adult_Unsorted_filt, vars.to.regress = "percent.mito")
tmp <- RunPCA(tmp)
tmp <- RunUMAP(tmp, dims = 1:40)
tmp <- FindNeighbors(tmp, reduction="pca", dims = 1:40)
tmp <- FindClusters(object = tmp, verbose = FALSE)

sweep.res.list <- paramSweep(tmp, PCs = 1:40, sct = TRUE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)

pK=as.numeric(as.character(bcmvn$pK))
BCmetric=bcmvn$BCmetric
pK_choose = pK[which(BCmetric %in% max(BCmetric))]

par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
plot(x = pK, y = BCmetric, pch = 16,type="b",
     col = "blue",lty=1)
abline(v=pK_choose,lwd=2,col='red',lty=2)
title("The BCmvn distributions")
text(pK_choose,max(BCmetric),as.character(pK_choose),pos = 4,col = "red")

nExp_poi <- round(0.075*nrow(tmp@meta.data))

S2_Adult_Unsorted_filt <- doubletFinder(tmp, PCs = 1:40, pN = 0.25, pK = pK_choose, nExp = nExp_poi, reuse.pANN = NULL, sct = TRUE)

saveRDS(object = S2_Adult_Unsorted_filt, file = "/data/Mlig_scRNA_Seq/Francisco/results/objs/S2_Adult_Unsorted_filt.rds")




#S7_Sorted_Adult
#filtering
qc.lib <- scater::isOutlier(S7_Sorted_Adult$nCount_RNA, nmads=2, log=TRUE, type="both")
qc.nexprs <- isOutlier(S7_Sorted_Adult$nFeature_RNA, log=TRUE, nmads=2, type="both")
qc.mito <- isOutlier(S7_Sorted_Adult$percent.mito, nmads=1, type="higher")
discard <- qc.lib | qc.nexprs  | qc.mito
S7_Sorted_Adult$discard<-discard
S7_Sorted_Adult_filt <-subset(S7_Sorted_Adult, subset= discard!=TRUE)

#doublet detection
tmp <- SCTransform(S7_Sorted_Adult_filt, vars.to.regress = "percent.mito")
tmp <- RunPCA(tmp)
tmp <- RunUMAP(tmp, dims = 1:40)
tmp <- FindNeighbors(tmp, reduction="pca", dims = 1:40)
tmp <- FindClusters(object = tmp, verbose = FALSE)

sweep.res.list <- paramSweep(tmp, PCs = 1:40, sct = TRUE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)

pK=as.numeric(as.character(bcmvn$pK))
BCmetric=bcmvn$BCmetric
pK_choose = pK[which(BCmetric %in% max(BCmetric))]

par(mar=c(5,4,4,8)+1,cex.main=1.2,font.main=2)
plot(x = pK, y = BCmetric, pch = 16,type="b",
     col = "blue",lty=1)
abline(v=pK_choose,lwd=2,col='red',lty=2)
title("The BCmvn distributions")
text(pK_choose,max(BCmetric),as.character(pK_choose),pos = 4,col = "red")

nExp_poi <- round(0.075*nrow(tmp@meta.data))

S7_Sorted_Adult_filt <- doubletFinder(tmp, PCs = 1:40, pN = 0.25, pK = pK_choose, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)

saveRDS(object = S7_Sorted_Adult_filt, file = "/data/Mlig_scRNA_Seq/Francisco/results/objs/S7_Sorted_Adult_filt.rds")



#cleaning up
rm(S17_0_1d_Hatchlings_counts, S18_0_1d_Hatchling_counts, S19_0_1d_Hatchling_counts, S11_unsorted_3dph_counts, S12_unsorted_3dph_counts, S16_3d_juveniles_counts, S8_unsorted_Hatchling_7d_counts, S15_7d_juveniles_counts, S9_unsorted_Hatchling_7d_counts, S1_Adult_Unsorted_counts, S2_Adult_Unsorted_DASH_counts, S7_Sorted_Adult_counts)

gc()


#S2_Adult_Unsorted_DASH_filt <- readRDS(file = "/data/Mlig_scRNA_Seq/Francisco/results/objs/S2_Adult_Unsorted_DASH_filt.rds")
#S17_0_1d_Hatchlings_filt <- readRDS(file = "/data/Mlig_scRNA_Seq/Francisco/results/objs/S17_0_1d_Hatchlings_filt.rds")
#S18_0_1d_Hatchling_filt <- readRDS(file = "/data/Mlig_scRNA_Seq/Francisco/results/objs/S18_0_1d_Hatchling_filt.rds")
#S19_0_1d_Hatchling_filt <- readRDS(file = "/data/Mlig_scRNA_Seq/Francisco/results/objs/S19_0_1d_Hatchling_filt.rds")
#S11_unsorted_3dph_filt <- readRDS(file = "/data/Mlig_scRNA_Seq/Francisco/results/objs/S11_unsorted_3dph_filt.rds")
#S12_unsorted_3dph_filt <- readRDS(file = "/data/Mlig_scRNA_Seq/Francisco/results/objs/S12_unsorted_3dph_filt.rds")
#S16_3d_juveniles_filt <- readRDS(file = "/data/Mlig_scRNA_Seq/Francisco/results/objs/S16_3d_juveniles_filt.rds")
#S8_unsorted_Hatchling_7d_filt <- readRDS(file = "/data/Mlig_scRNA_Seq/Francisco/results/objs/S8_unsorted_Hatchling_7d_filt.rds")
#S15_7d_juveniles_filt <- readRDS(file = "/data/Mlig_scRNA_Seq/Francisco/results/objs/S15_7d_juveniles_filt.rds")
#S9_unsorted_Hatchling_7d_filt <- readRDS(file = "/data/Mlig_scRNA_Seq/Francisco/results/objs/S9_unsorted_Hatchling_7d_filt.rds")
#S1_Adult_Unsorted_filt <- readRDS(file = "/data/Mlig_scRNA_Seq/Francisco/results/objs/S1_Adult_Unsorted_filt.rds")
#S2_Adult_Unsorted_DASH_filt <- readRDS(file = "/data/Mlig_scRNA_Seq/Francisco/results/objs/S2_Adult_Unsorted_DASH_filt.rds")
#S7_Sorted_Adult_filt <- readRDS(file = "/data/Mlig_scRNA_Seq/Francisco/results/objs/S7_Sorted_Adult_filt.rds")

#merging
#All <- merge(S17_0_1d_Hatchlings, y =c(S18_0_1d_Hatchling, S19_0_1d_Hatchling, S11_unsorted_3dph, S12_unsorted_3dph, S16_3d_juveniles, S8_unsorted_Hatchling_7d, S15_7d_juveniles, S9_unsorted_Hatchling_7d, S1_Adult_Unsorted, S2_Adult_Unsorted_DASH, S7_Sorted_Adult), add.cell.ids = c("S17_0_1d_Hatchlings", "S18_0_1d_Hatchling", "S19_0_1d_Hatchling", "S11_unsorted_3dph", "S12_unsorted_3dph", "S16_3d_juveniles", "S8_unsorted_Hatchling_7d", "S15_7d_juveniles", "S9_unsorted_Hatchling_7d", "S1_Adult_Unsorted",  "S2_Adult_Unsorted_DASH", "S7_Sorted_Adult"), project = "scRNA_Seq_M_lignano")

All_filtered <- merge(S17_0_1d_Hatchlings_filt, y =c(S18_0_1d_Hatchling_filt, S19_0_1d_Hatchling_filt, S11_unsorted_3dph_filt, S12_unsorted_3dph_filt, S16_3d_juveniles_filt, S8_unsorted_Hatchling_7d_filt, S15_7d_juveniles_filt, S9_unsorted_Hatchling_7d_filt, S1_Adult_Unsorted_filt, S2_Adult_Unsorted_DASH_filt, S7_Sorted_Adult_filt), add.cell.ids = c("S17_0_1d_Hatchlings_filt", "S18_0_1d_Hatchling_filt", "S19_0_1d_Hatchling_filt", "S11_unsorted_3dph_filt", "S12_unsorted_3dph_filt", "S16_3d_juveniles_filt", "S8_unsorted_Hatchling_7d_filt", "S15_7d_juveniles_filt", "S9_unsorted_Hatchling_7d_filt", "S1_Adult_Unsorted_filt",  "S2_Adult_Unsorted_DASH_filt", "S7_Sorted_Adult_filt"), project = "scRNA_Seq_M_lignano")


#filtering

All_filtered_mit_35 <- subset(x = All, subset = (percent.mito < 35 & nFeature_RNA > 200 & nFeature_RNA < 1000 & nCount_RNA < 5000))
#All_filtered_mit_35 <- subset(x = All, subset = (mito.percentage < 35 & nFeature_RNA > 200 & nFeature_RNA < 2500 & nCount_RNA < 5000))

#All_filtered_mit_50 <- subset(x = All, subset = (mt_cod < 50 & nFeature_RNA > 200 & nFeature_RNA < 2500 & nCount_RNA < 10000))

gc()

#save raw & processed obj

saveRDS(object = All, file = "/data/Mlig_scRNA_Seq/Francisco/results/objs/All_merged_cellbender_round_two.rds")
saveRDS(object = All_filtered, file = "/data/Mlig_scRNA_Seq/Francisco/results/objs/All_merged_filt_doubletfinder_cellbender_round_two.rds")
saveRDS(object = All_filtered_mit_35, file = "/data/Mlig_scRNA_Seq/Francisco/results/objs/All_merged_mito_35_cellbender_round_two.rds")

VlnPlot(All, group.by = "orig.ident", features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3, alpha = 0, raster = FALSE)


rm(list=ls())
gc()

sessionInfo()
