library(stringr)

srt_obj <- readRDS(file = "/data/Mlig_scRNA_Seq/Francisco/Mlig_scRNASeq_atlas/results/objs/srt_obj_norm_scaled_RNA_plus_SCT_plus_integration.rds")

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

srt_obj <- add_module_scores(srt_obj, marker_genes_list, assay = "RNA")

srt_obj$Stage <- factor(srt_obj$Stage, levels = c("0-1_days", "3_days", "7_days", "Adult"))

srt_obj$sample_id <- factor(srt_obj$sample_id, levels = c("17_0_1d_Hatchlings",
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

saveRDS(srt_obj, file = "/data/Mlig_scRNA_Seq/Francisco/Mlig_scRNASeq_atlas/results/objs/srt_obj_norm_scaled_RNA_plus_SCT_plus_integration_plus_annotation.rds")

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
