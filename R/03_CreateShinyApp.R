library(future); packageVersion("future")
library(ShinyCell2)
library(Seurat)          # load after ShinyCell2 to avoid it pinning a different future

setwd("/data/Mlig_scRNA_Seq/Francisco/Mlig_scRNASeq_atlas/results/objs/")

srt_obj <- readRDS("srt_obj_norm_scaled_RNA_plus_SCT_plus_integration_plus_annotation.rds")

scConf <- createConfig(srt_obj, maxLevels = 100)

srt_obj <- FindVariableFeatures(srt_obj, assay = "RNA")

makeShinyFiles(srt_obj,
               scConf,
               shiny.prefix="M_lignano",
               shiny.dir="M_lignano_shinyApp/")

makeShinyCodes(shiny.title = "M_lignano",
               shiny.prefix="M_lignano",
               shiny.dir="M_lignano_shinyApp/")
