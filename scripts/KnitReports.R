## knit document 00_CreateSeuratFromCounts.Rmd

setwd("~/Projects/Erin/M_lignano/Mlig_scRNASeq_atlas/Rmarkdown/")

rmarkdown::render(
  "00_CreateSeuratFromCounts.Rmd",
  output_file = "00_CreateSeuratFromCounts.html",
  output_dir  = "../results/rmarkdown_reports"
)

rm(list=ls())
gc()

rmarkdown::render(
  "01_ComputeQCMetrics.Rmd",
  output_file = "01_ComputeQCMetrics.html",
  output_dir  = "../results/rmarkdown_reports"
)

rm(list=ls())
gc()

rmarkdown::render(
  "02_IntegrateSeurat.Rmd",
  output_file = "02_IntegrateSeurat.html",
  output_dir  = "../results/rmarkdown_reports"
)

rm(list=ls())
gc()

# rmarkdown::render(
#   "03_AnnotateCellTypes.Rmd",
#   output_file = "03_AnnotateCellTypes.html",
#   output_dir  = "../results/rmarkdown_reports"
# )
#
# rmarkdown::render(
#   "03_RunAnalysis.Rmd",
#   output_file = "03_RunAnalysis.html",
#   output_dir  = "../results/rmarkdown_reports"
# )
# 
# rmarkdown::render(
#   "04_PlotPaperFigures.Rmd",
#   output_file = "04_PlotPaperFigures.html",
#   output_dir  = "../results/rmarkdown_reports"
# )
