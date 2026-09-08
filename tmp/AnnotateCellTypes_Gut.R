library(Seurat)
library(scCustomize)

library(circlize)
col_fun = colorRamp2(c(0, 1), c("blue", "red"))

srt_obj <- readRDS("/data/Mlig_scRNA_Seq/Francisco/Mlig_scRNASeq_atlas/results/objs/Digestive_System_integrated_merged.rds")

p_orig_anno <- DimPlot_scCustom(srt_obj,
                                group.by = "digestive_cluster",
                                aspect_ratio = 1,
                                colors_use = DiscretePalette_scCustomize(
                                  num_colors = length(unique(srt_obj$digestive_cluster)),
                                  palette = "ditto_seq"
                                )
)

### Adding SAMap - S. mediterranea

#read the output of SAMap
tmp <- read.csv("/data/Mlig_scRNA_Seq/Francisco/Mlig_scRNASeq_atlas/results/cross_spp_annotation/tables/cluster_annotation_gut.tsv", header = TRUE, row.names = 1)

#remove M. lignano duplicates from columns
pattern <- "^ml_"
tmp <- tmp[, !grepl(pattern, colnames(tmp))]

#remove S. mediterranea duplicates from rows
pattern <- "^sm_"
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
srt_obj$cluster_ID_cross_species_annotation_sm_relaxed <- unname(cluster_to_max_feature[as.character(srt_obj$digestive_cluster)])

p_cross_spp_annot_relaxed <- DimPlot_scCustom(srt_obj,
                                              group.by = "cluster_ID_cross_species_annotation_sm_relaxed",
                                              aspect_ratio = 1,
                                              colors_use = DiscretePalette_scCustomize(
                                                num_colors = length(unique(srt_obj$cluster_ID_cross_species_annotation_sm_relaxed)),
                                                palette = "ditto_seq"
                                              )
)


# cell_clusters <- Idents(seurat_integrated_final)
# cell_level_data <- avg_features_per_cluster[as.character(cell_clusters),]
# colnames(cell_level_data) <- paste0("Score_cross_spp_", colnames(cell_level_data))
# 
# seurat_integrated_final <- AddMetaData(seurat_integrated_final, metadata = cell_level_data)


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
srt_obj$cluster_ID_cross_species_annotation_sm_stringent <- unname(cluster_to_max_feature[as.character(srt_obj$digestive_cluster)])


p_cross_spp_annot_stringent <- DimPlot_scCustom(srt_obj,
                                              group.by = "cluster_ID_cross_species_annotation_sm_stringent",
                                              aspect_ratio = 1,
                                              colors_use = DiscretePalette_scCustomize(
                                                num_colors = length(unique(srt_obj$cluster_ID_cross_species_annotation_sm_stringent)),
                                                palette = "ditto_seq"
                                              )
)


## Heatmap comparing atlases

#read the output of SAMap
tmp <- read.csv("/data/Mlig_scRNA_Seq/Francisco/Mlig_scRNASeq_atlas/results/cross_spp_annotation/tables/cluster_annotation_gut.tsv", header = TRUE, row.names = 1)

#remove M. lignano duplicates from columns
pattern <- "^ml_"
tmp <- tmp[, !grepl(pattern, colnames(tmp))]

#remove S. mediterranea from rows
pattern <- "^sm_"
tmp <- tmp[!grepl(pattern, rownames(tmp)),]

heatmap_sm_vs_sc <- ComplexHeatmap::Heatmap(
  as.matrix(tmp),
  heatmap_legend_param = list(title = "Score"),
  col = col_fun
)

pdf("/data/Mlig_scRNA_Seq/Francisco/Mlig_scRNASeq_atlas/results/figures/cross_spp_heatmap_S_mediterranea_vs_S_mansoni.pdf", width = 10, height = 10)
heatmap_sm_vs_sc
dev.off()

## M. lignano vs S. mediterranea


#read the output of SAMap
tmp <- read.csv("/data/Mlig_scRNA_Seq/Francisco/Mlig_scRNASeq_atlas/results/cross_spp_annotation/tables/cluster_annotation_integrated_subcluster.tsv", header = TRUE, row.names = 1)

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

heatmap_ml_vs_sm <- ComplexHeatmap::Heatmap(
  as.matrix(tmp),
  heatmap_legend_param = list(title = "Score"),
  col = col_fun
)

pdf("/data/Mlig_scRNA_Seq/Francisco/Mlig_scRNASeq_atlas/results/figures/cross_spp_heatmap_M_lignano_vs_S_mediterranea.pdf", width = 10, height = 10)
heatmap_ml_vs_sm
dev.off()



## M. lignano vs S. mansoni

#read the output of SAMap
tmp <- read.csv("/data/Mlig_scRNA_Seq/Francisco/Mlig_scRNASeq_atlas/results/cross_spp_annotation/tables/cluster_annotation_integrated_1_8.tsv", header = TRUE, row.names = 1)

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

heatmap_ml_vs_sm <- ComplexHeatmap::Heatmap(
  as.matrix(tmp),
  heatmap_legend_param = list(title = "Score"),
  col = col_fun
  )

pdf("/data/Mlig_scRNA_Seq/Francisco/Mlig_scRNASeq_atlas/results/figures/cross_spp_heatmap_M_lignano_vs_S_mansoni.pdf", width = 10, height = 10)
heatmap_ml_vs_sc
dev.off()



avg_features_per_cluster <- srt_obj@meta.data %>%
  dplyr::group_by(sub.cluster.33.0.15.sub.cluster.15.0.1.sub.cluster.30.res.0.1) %>%
  dplyr::summarise(across(all_of(plot_order), ~ mean(.x, na.rm = FALSE), .names = "{col}"))
# #   dplyr::summarise(across(all_of(features), ~ median(.x, na.rm = FALSE), .names = "{col}"))
# Calculate the median values for each feature by cluster

avg_features_per_cluster <- as.data.frame(avg_features_per_cluster)
rownames(avg_features_per_cluster) <- avg_features_per_cluster$sub.cluster.33.0.15.sub.cluster.15.0.1.sub.cluster.30.res.0.1
avg_features_per_cluster <- avg_features_per_cluster[,2:length(colnames(avg_features_per_cluster))]
cluster_to_max_feature <- vec <- character(length = length(rownames(avg_features_per_cluster)))
names(cluster_to_max_feature) <- rownames(avg_features_per_cluster)
#iterate through features (module scores) and clusters to select, for each cluster, its ID based on average module score

for (i in rownames(avg_features_per_cluster)) {
  max_val <- -10000
  cluster_ID <- "NA"
  for (i2 in colnames(avg_features_per_cluster)) {
     if (avg_features_per_cluster[i,i2] > max_val) {
       if (avg_features_per_cluster[i,i2] > 0) {
         #      if (avg_features_per_cluster[i,i2] > 0.01) { #cutoff to associate
         cluster_ID <- i2
         max_val <- avg_features_per_cluster[i,i2]
       }
     }
   }
   cluster_to_max_feature[i] <- cluster_ID
   print(paste0("Cluster ID: ", i, " - annotation: ", cluster_ID, " - value:", max_val))
}

# # #create a new metadata layer with manual curation cluster IDs
# # All.sct$cluster_ID_manual_curation <- unname(cluster_to_max_feature[All.sct$integrated_snn_res.4])
# # 
# 



