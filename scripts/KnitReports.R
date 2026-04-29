## knit Rmarkdown documents

setwd("~/Projects/Erin/M_lignano/Mlig_scRNASeq_atlas/")

if (file.exists(".Rprofile")) source(".Rprofile")
renv::load(project = getwd())

#rmarkdown::render(
#  "00_CreateSeuratFromCounts.Rmd",
#  output_file = "00_CreateSeuratFromCounts.html",
#  output_dir  = "../results/rmarkdown_reports"
#)

#rm(list=ls())
#gc()

#rmarkdown::render(
#  "01_ComputeQCMetrics.Rmd",
#  output_file = "01_ComputeQCMetrics.html",
#  output_dir  = "../results/rmarkdown_reports"
#)

#rm(list=ls())
#gc()

#rmarkdown::render(
#  "02_IntegrateSeurat.Rmd",
#  output_file = "02_IntegrateSeurat.html",
#  output_dir  = "../results/rmarkdown_reports"
#)

#rm(list=ls())
#gc()

#rmarkdown::render(
#  "Rmarkdown/04_IntegrateSeuratObjects.Rmd",
#  output_file = "04_IntegrateSeuratObjects.html",
#  output_dir  = "/data/Mlig_scRNA_Seq/Francisco/Mlig_scRNASeq_atlas/results/rmarkdown_reports"
#)

#rmarkdown::render(
#  "Rmarkdown/05_ComputeModuleScores.Rmd",
#   output_file = "05_ComputeModuleScores.html",
#   output_dir  = "/data/Mlig_scRNA_Seq/Francisco/Mlig_scRNASeq_atlas/results/rmarkdown_reports"
#)

#rmarkdown::render(
#  "Rmarkdown/07_PrepareReferenceAtlas_S_mediterranea.Rmd",
#  output_file = "07_PrepareReferenceAtlas_S_mediterranea.html",
#  output_dir  = "/data/Mlig_scRNA_Seq/Francisco/Mlig_scRNASeq_atlas/results/rmarkdown_reports"
#)

#rmarkdown::render(
#		    "Rmarkdown/07_2_Subset_S_mediterranea_Intestine_to_H5AD.Rmd",
#		      output_file = "07_2_Subset_S_mediterranea_Intestine_to_H5AD.html",
#		      output_dir  = "/data/Mlig_scRNA_Seq/Francisco/Mlig_scRNASeq_atlas/results/rmarkdown_reports"
#		      )

rmarkdown::render(
		  "Rmarkdown/09_GO_enrichment_gut.Rmd",
		  output_file = "09_GO_enrichment_gut.html",
		  output_dir  = "/data/Mlig_scRNA_Seq/Francisco/Mlig_scRNASeq_atlas/results/rmarkdown_reports"
)


#rmarkdown::render(
#  "Rmarkdown/08_PrepareReferenceAtlas_S_mansoni.Rmd",
#  output_file = "08_PrepareReferenceAtlas_S_mansoni.html",
#  output_dir  = "/data/Mlig_scRNA_Seq/Francisco/Mlig_scRNASeq_atlas/results/rmarkdown_reports"
#)
~

