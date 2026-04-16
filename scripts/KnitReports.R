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

rmarkdown::render(
  "Rmarkdown/05_ComputeModuleScores.Rmd",
   output_file = "05_ComputeModuleScores.html",
   output_dir  = "/data/Mlig_scRNA_Seq/Francisco/Mlig_scRNASeq_atlas/results/rmarkdown_reports"
)

