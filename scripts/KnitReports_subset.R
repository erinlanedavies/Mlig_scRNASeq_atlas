## knit Rmarkdown documents

setwd("~/Projects/Erin/M_lignano/Mlig_scRNASeq_atlas/")

if (file.exists(".Rprofile")) source(".Rprofile")
renv::load(project = getwd())


#rmarkdown::render(
#                    "Rmarkdown/06_1_Export_M_lignano_Neural_to_H5AD.Rmd",
#                     output_file = "06_1_Export_M_lignano_Neural_to_H5AD.html",
#                     output_dir  = "/data/Mlig_scRNA_Seq/Francisco/Mlig_scRNASeq_atlas/results/rmarkdown_reports"
#                     )

#rmarkdown::render(
#                    "Rmarkdown/06_2_Export_M_lignano_Digestive_System_to_H5AD.Rmd",
#                     output_file = "06_2_Export_M_lignano_Digestive_System_to_H5AD.html",
#                     output_dir  = "/data/Mlig_scRNA_Seq/Francisco/Mlig_scRNASeq_atlas/results/rmarkdown_reports"
#                     )

#rmarkdown::render(  
#		    "Rmarkdown/06_3_Export_M_lignano_Cathepsin_to_H5AD.Rmd",
#		      output_file = "06_3_Export_M_lignano_Cathepsin_to_H5AD.html",
#		      output_dir  = "/data/Mlig_scRNA_Seq/Francisco/Mlig_scRNASeq_atlas/results/rmarkdown_reports"
#		      )

rmarkdown::render(  
		    "Rmarkdown/06_4_Export_M_lignano_Musculature_to_H5AD.Rmd",
		      output_file = "06_4_Export_M_lignano_Musculature_to_H5AD.html",
		      output_dir  = "/data/Mlig_scRNA_Seq/Francisco/Mlig_scRNASeq_atlas/results/rmarkdown_reports"
		      )

#rmarkdown::render(  
#		    "Rmarkdown/06_5_Export_M_lignano_Progenitor_System_to_H5AD.Rmd",
#		      output_file = "06_5_Export_M_lignano_Progenitor_System_to_H5AD.html",
#		      output_dir  = "/data/Mlig_scRNA_Seq/Francisco/Mlig_scRNASeq_atlas/results/rmarkdown_reports"
#		      )

#rmarkdown::render(   
#                   "Rmarkdown/06_6_Export_M_lignano_Epithelial_to_H5AD.Rmd",
#                     output_file = "06_6_Export_M_lignano_Epithelial_to_H5AD.html",
#                     output_dir  = "/data/Mlig_scRNA_Seq/Francisco/Mlig_scRNASeq_atlas/results/rmarkdown_reports"
#                     )

#rmarkdown::render(  
#		    "Rmarkdown/07_1_Subset_S_mediterranea_Intestine_to_H5AD.Rmd",
#		      output_file = "07_1_Subset_S_mediterranea_Intestine_to_H5AD.html",
#		      output_dir  = "/data/Mlig_scRNA_Seq/Francisco/Mlig_scRNASeq_atlas/results/rmarkdown_reports"
#		      )

#rmarkdown::render(
#                    "Rmarkdown/07_2_Subset_S_mediterranea_Cathepsin_positive_cells_to_H5AD.Rmd",
#                      output_file = "07_2_Subset_S_mediterranea_Cathepsin_positive_cells_to_H5AD.html",
#                      output_dir  = "/data/Mlig_scRNA_Seq/Francisco/Mlig_scRNASeq_atlas/results/rmarkdown_reports"
#)

#rmarkdown::render(
#		    "Rmarkdown/07_3_Subset_S_mediterranea_Epidermal_cells_to_H5AD.Rmd",
#		      output_file = "07_3_Subset_S_mediterranea_Epidermal_cells_to_H5AD.html",
#		      output_dir  = "/data/Mlig_scRNA_Seq/Francisco/Mlig_scRNASeq_atlas/results/rmarkdown_reports"
#		      )

#rmarkdown::render(
#		    "Rmarkdown/07_4_Subset_S_mediterranea_Muscle_cells_to_H5AD.Rmd",
#		      output_file = "07_4_Subset_S_mediterranea_Muscle_cells_to_H5AD.html",
#		      output_dir  = "/data/Mlig_scRNA_Seq/Francisco/Mlig_scRNASeq_atlas/results/rmarkdown_reports"
#		      )

#rmarkdown::render(
#		    "Rmarkdown/07_5_Subset_S_mediterranea_Neural_cells_to_H5AD.Rmd",
#		      output_file = "07_5_Subset_S_mediterranea_Neural_cells_to_H5AD.html",
#		      output_dir  = "/data/Mlig_scRNA_Seq/Francisco/Mlig_scRNASeq_atlas/results/rmarkdown_reports"
#		      )

#rmarkdown::render(
#		    "Rmarkdown/07_6_Subset_S_mediterranea_Neoblast_cells_to_H5AD.Rmd",
#		      output_file = "07_6_Subset_S_mediterranea_Neoblast_cells_to_H5AD.Rmd",
#		      output_dir  = "/data/Mlig_scRNA_Seq/Francisco/Mlig_scRNASeq_atlas/results/rmarkdown_reports"
#		      )

#rmarkdown::render(
#		    "Rmarkdown/08_7_Subset_S_mansoni_progeny_cells_to_H5AD.Rmd",
#		      output_file = "08_7_Subset_S_mansoni_progeny_cells_to_H5AD.html",
#		      output_dir  = "/data/Mlig_scRNA_Seq/Francisco/Mlig_scRNASeq_atlas/results/rmarkdown_reports"
#		      )

#rmarkdown::render(
#                   "Rmarkdown/08_8_Subset_S_mansoni_gut_cells_to_H5AD.Rmd",
#                     output_file = "08_8_Subset_S_mansoni_gut_cells_to_H5AD.html",
#                     output_dir  = "/data/Mlig_scRNA_Seq/Francisco/Mlig_scRNASeq_atlas/results/rmarkdown_reports"
#                     )

