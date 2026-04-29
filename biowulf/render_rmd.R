#!/usr/bin/env Rscript

options(error = function() {
		  cat("\n--- TRACEBACK ---\n")
		    traceback(20)
		    quit(status = 1)
})

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
activate <- file.path(project_root, "renv", "activate.R")

if (!file.exists(activate)) {
	  stop("renv activation file not found: ", activate)
}

source(activate)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1) {
	  stop("Usage: Rscript biowulf/render_rmd.R path/to/file.Rmd")
}

rmd_file <- args[[1]]

if (!file.exists(rmd_file)) {
	  stop("Rmd file not found: ", rmd_file)
}

message("Rendering: ", rmd_file)
message("Using R library paths:")
print(.libPaths())

out <- rmarkdown::render(
			   input = rmd_file,
			     envir = new.env(parent = globalenv()),
			     clean = FALSE
			     )

message("Render completed successfully.")
message("Output file: ", out)

quit(status = 0)
