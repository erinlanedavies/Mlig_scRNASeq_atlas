marker_genes <- read.table("data/metadata/literature_marker_genes.tsv", header = TRUE)

marker_genes$CleanModuleID <- ""
marker_genes$CleanModuleID <- ifelse(marker_genes$Cell_tissue_type == "stem_cells_progenitors", "Progenitor", marker_genes$CleanModuleID)
marker_genes$CleanModuleID <- ifelse(marker_genes$Cell_tissue_type == "nervous_system", "Neural", marker_genes$CleanModuleID)
marker_genes$CleanModuleID <- ifelse(marker_genes$Cell_tissue_type == "muscle", "Muscle", marker_genes$CleanModuleID)
marker_genes$CleanModuleID <- ifelse(marker_genes$Cell_tissue_type == "gut", "Gut", marker_genes$CleanModuleID)
marker_genes$CleanModuleID <- ifelse(marker_genes$Cell_tissue_type == "gut_specialized_cell", "Secretory Gut", marker_genes$CleanModuleID)
marker_genes$CleanModuleID <- ifelse(marker_genes$Cell_tissue_type == "rhabdite-containing_cell", "Rhabdite-Containing", marker_genes$CleanModuleID)
marker_genes$CleanModuleID <- ifelse(marker_genes$Cell_tissue_type == "epidermal_secretory_cell", "Epidermal Secretory", marker_genes$CleanModuleID)
marker_genes$CleanModuleID <- ifelse(marker_genes$Cell_tissue_type == "secretory_cell_of_adhesive_organ", "Secretory Cell of Adhesive Organ", marker_genes$CleanModuleID)
marker_genes$CleanModuleID <- ifelse(marker_genes$Cell_tissue_type == "cement_gland_cell", "Cement Gland", marker_genes$CleanModuleID)
marker_genes$CleanModuleID <- ifelse(marker_genes$Cell_tissue_type == "female_antrum", "Female Antrum", marker_genes$CleanModuleID)
marker_genes$CleanModuleID <- ifelse(marker_genes$Cell_tissue_type == "prostate", "Prostate", marker_genes$CleanModuleID)
marker_genes$CleanModuleID <- ifelse(marker_genes$Cell_tissue_type == "female_germline_ovary", "Ovary", marker_genes$CleanModuleID)
marker_genes$CleanModuleID <- ifelse(marker_genes$Cell_tissue_type == "male_germline_testis", "Testis", marker_genes$CleanModuleID)
marker_genes$CleanModuleID <- ifelse(marker_genes$Cell_tissue_type == "epidermal", "Epidermal", marker_genes$CleanModuleID)

marker_genes$CleanModuleIDBroad <- ""
marker_genes$CleanModuleIDBroad <- ifelse(marker_genes$Cell_tissue_type == "stem_cells_progenitors", "Progenitor", marker_genes$CleanModuleIDBroad)
marker_genes$CleanModuleIDBroad <- ifelse(marker_genes$Cell_tissue_type == "nervous_system", "Neural", marker_genes$CleanModuleIDBroad)
marker_genes$CleanModuleIDBroad <- ifelse(marker_genes$Cell_tissue_type == "muscle", "Muscle", marker_genes$CleanModuleIDBroad)
marker_genes$CleanModuleIDBroad <- ifelse(marker_genes$Cell_tissue_type == "gut", "Gut", marker_genes$CleanModuleIDBroad)
marker_genes$CleanModuleIDBroad <- ifelse(marker_genes$Cell_tissue_type == "gut_specialized_cell", "Secretory Cells", marker_genes$CleanModuleIDBroad)
marker_genes$CleanModuleIDBroad <- ifelse(marker_genes$Cell_tissue_type == "rhabdite-containing_cell", "Secretory Cells", marker_genes$CleanModuleIDBroad)
marker_genes$CleanModuleIDBroad <- ifelse(marker_genes$Cell_tissue_type == "epidermal_secretory_cell", "Secretory Cells", marker_genes$CleanModuleIDBroad)
marker_genes$CleanModuleIDBroad <- ifelse(marker_genes$Cell_tissue_type == "secretory_cell_of_adhesive_organ", "Secretory Cells", marker_genes$CleanModuleIDBroad)
marker_genes$CleanModuleIDBroad <- ifelse(marker_genes$Cell_tissue_type == "cement_gland_cell", "Secretory Cells", marker_genes$CleanModuleIDBroad)
marker_genes$CleanModuleIDBroad <- ifelse(marker_genes$Cell_tissue_type == "female_antrum", "Female reproduction", marker_genes$CleanModuleIDBroad)
marker_genes$CleanModuleIDBroad <- ifelse(marker_genes$Cell_tissue_type == "prostate", "Male reproduction", marker_genes$CleanModuleIDBroad)
marker_genes$CleanModuleIDBroad <- ifelse(marker_genes$Cell_tissue_type == "female_germline_ovary", "Female reproduction", marker_genes$CleanModuleIDBroad)
marker_genes$CleanModuleIDBroad <- ifelse(marker_genes$Cell_tissue_type == "male_germline_testis", "Male reproduction", marker_genes$CleanModuleIDBroad)
marker_genes$CleanModuleIDBroad <- ifelse(marker_genes$Cell_tissue_type == "epidermal", "Epidermal", marker_genes$CleanModuleIDBroad)

marker_genes <- marker_genes[marker_genes$Cell_tissue_type != "anchor_cell_of_adhesive_organ",]
marker_genes <- marker_genes[marker_genes$Cell_tissue_type != "epidermal",]

write.table(marker_genes, file = "data/metadata/literature_marker_genes2.tsv", sep = "\t", col.names = TRUE, quote = FALSE)
