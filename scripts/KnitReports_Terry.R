## knit Rmarkdown documents

setwd("~/Projects/Erin/M_lignano/Mlig_scRNASeq_atlas/")

if (file.exists(".Rprofile")) source(".Rprofile")
renv::load(project = getwd())

#rmarkdown::render(
#  "Rmarkdown/Export_Gottgens_to_H5AD.Rmd",
#  output_file = "Export_Gottgens_to_H5AD.html",
#  output_dir  = "/data/YamaguchiLab/Francisco/gastruloid/results/08_06_2026/"
#)

#rmarkdown::render(
#  "Rmarkdown/Export_Yamaguchi_to_H5AD.Rmd",
#   output_file = "Export_Yamaguchi_to_H5AD.html",
#   output_dir  = "/data/YamaguchiLab/Francisco/gastruloid/results/08_06_2026/"
#)

rmarkdown::render(
  "Rmarkdown/Export_Tyser_to_H5AD.Rmd",
   output_file = "Export_Tyser_to_H5AD.html",
   output_dir  = "/data/YamaguchiLab/Francisco/HumanGastruloid/results/09_18_2026/"
)

