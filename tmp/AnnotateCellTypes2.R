library(Seurat)
library(scCustomize)

library(circlize)
col_fun = colorRamp2(c(0, 1), c("blue", "red"))

srt_obj <- readRDS("/data/Mlig_scRNA_Seq/Francisco/Mlig_scRNASeq_atlas/results/objs/srt_integrated_all_module_scores_subclustered.rds")

plot_order <- c(
  "Progenitor",
  "Neural",
  "Muscle",
  "Gut",
  "Secretory Gut",
  "Rhabdite-containing",
  "Stylet",
  "Regenerating Stylet",
  "Epidermis",
#  "Epidermal Secretory",
  "Anchor Cell of Adhesive Organ",
  "Secretory Cell of Adhesive Organ",
  "Paired Ventral Cells",
  "Single Cells in Tail",
  "Cement Gland",
  #  "Female Antrum",
  "Prostate",
  "Ovary",
  "Testis"
)

cluster_colors_manual_curation <- c(#"#E0AFCA", #anchor_cell_of_adhesive_organ
  "#AD7700", #Cement gland cell
  "#0072B2",#Epidermal/Secretory
  #"#8D3666", #female_antrum
  "#00446B", #Female germline/Ovary
  "#CC79A7", #Gut
  "#FF8320", #gut_specialized
  "#D55E00", #Male germline/Testis
  "#D22B2B", #Muscle
  "#4D4D4D", #NA
  "#007756", #nervous_system
  "#00F6B3", #neural_cell_in_dorsal_tail
  #"#009E73", #neural_cell_in_tail_paired
  "#1C91D4", #Progenitor
  "#A04700", #Prostate
  "#8D3666" #Rhabdite
) 

cluster_colors_cross_species <- c("#AA9F0D", #Cathepsin..cells
                                  "#FFBE2D", #Epidermal
                                  "#CC79A7", #Intestine
                                  "#D22B2B", #Muscle
                                  "#4D4D4D", #NA
                                  "#1C91D4", #Neoblast
                                  "#009E73", #Neural
                                  "#F4EB71", #Parapharyngeal
                                  "#D5C711" #Pharynx
)

cluster_colors_final <- c(#"#E0AFCA", #anchor_cell_of_adhesive_organ
  "#AA9F0D", #Cathepsin + cells
  "#AD7700", #Cement gland cell
  "#FFBE2D", #Epidermal
  "#0072B2",#Epidermal/Secretory
  #"#8D3666", #female_antrum
  "#00446B", #Female germline/Ovary
  "#CC79A7", #Gut
  "#FF8320", #gut_specialized
  "#D55E00", #Male germline/Testis
  "#D22B2B", #Muscle
  "#4D4D4D", #NA
  "#1C91D4", #Neoblast
  "#007756", #nervous_system
  "#00F6B3", #neural_cell_in_dorsal_tail
  #"#009E73", #neural_cell_in_tail_paired
  "#D5C711", #Pharynx
  "#A04700", #Prostate
  "#8D3666" #Rhabdite
)


### Adding SAMap - S. mediterranea

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
srt_obj$cluster_ID_cross_species_annotation_sm_relaxed <- unname(cluster_to_max_feature[as.character(srt_obj$sub.cluster.33.0.15.sub.cluster.15.0.1.sub.cluster.30.res.0.1)])

# cell_clusters <- Idents(seurat_integrated_final)
# cell_level_data <- avg_features_per_cluster[as.character(cell_clusters),]
# colnames(cell_level_data) <- paste0("Score_cross_spp_", colnames(cell_level_data))
# 
# seurat_integrated_final <- AddMetaData(seurat_integrated_final, metadata = cell_level_data)

DimPlot_scCustom(srt_obj, group.by = "cluster_ID_cross_species_annotation_sm_relaxed")


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
srt_obj$cluster_ID_cross_species_annotation_sm_stringent <- unname(cluster_to_max_feature[as.character(srt_obj$sub.cluster.33.0.15.sub.cluster.15.0.1.sub.cluster.30.res.0.1)])

DimPlot_scCustom(srt_obj, group.by = "cluster_ID_cross_species_annotation_sm_stringent", aspect_ratio = 1)


### Adding SAMap - S. mansoni

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
srt_obj$cluster_ID_cross_species_annotation_sc_relaxed <- unname(cluster_to_max_feature[as.character(srt_obj$sub.cluster.33.0.15.sub.cluster.15.0.1.sub.cluster.30.res.0.1)])

DimPlot_scCustom(srt_obj, group.by = "cluster_ID_cross_species_annotation_sc_relaxed", aspect_ratio = 1)



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
srt_obj$cluster_ID_cross_species_annotation_sc_stringent <- unname(cluster_to_max_feature[as.character(srt_obj$sub.cluster.33.0.15.sub.cluster.15.0.1.sub.cluster.30.res.0.1)])

DimPlot_scCustom(srt_obj, group.by = "cluster_ID_cross_species_annotation_sc_stringent", aspect_ratio = 1)

#ComplexHeatmap::Heatmap(log10(table(srt_obj$integrated_snn_res.1.8, srt_obj$cluster_ID_cross_species_annotation_sm_relaxed)+1))
heatmap(table(srt_obj$sub.cluster.33.0.15.sub.cluster.15.0.1.sub.cluster.30.res.0.1, srt_obj$cluster_ID_cross_species_annotation_sm_stringent))

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




## Heatmap comparing atlases

#read the output of SAMap
tmp <- read.csv("/data/Mlig_scRNA_Seq/Francisco/Mlig_scRNASeq_atlas/results/cross_spp_annotation/tables/cluster_annotation_integrated_subcluster.tsv", header = TRUE, row.names = 1)

#remove M. lignano duplicates from columns
pattern <- "^ml_"
tmp <- tmp[, !grepl(pattern, colnames(tmp))]

#remove S. mansoni from columns
pattern <- "^sc_"
tmp <- tmp[, !grepl(pattern, colnames(tmp))]

#remove M. lignano from rows
pattern <- "^ml_"
tmp <- tmp[!grepl(pattern, rownames(tmp)),]

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

heatmap_ml_vs_sc <- ComplexHeatmap::Heatmap(
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



