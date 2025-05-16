#!/usr/bin/env Rscript

library(optparse)
require(tidyverse)
library(testthat)

# Get the path to the 'cellsnake' package
cellsnake_path <- system("python -c 'import cellsnake; print(cellsnake.__path__[0])'", intern = TRUE)

# Define the base path to the testData folder
test_data_dir <- file.path(cellsnake_path, "scrna/workflow/tests/testData")

# Define options with dynamically constructed paths
option_list <- list(
  make_option(c("--xlsx"), type = "character", default = file.path(test_data_dir, "results/defaultTest/table_positive-markers-seurat_clusters.xlsx"),
              help = "Excel table of markers for input", metavar = "character"),
  make_option(c("--output.plot"), type = "character", default = file.path(test_data_dir, "results/summarized_markers-for-seurat_clusters.pdf"),
              help = "Output plot file", metavar = "character")
)

opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

try(
  {
    source(file.path(cellsnake_path,"scrna/workflow/scripts/scrna_workflow_helper_functions/scrna-top-marker-plot-functions.R"))
  },
  silent = TRUE
)


test_that("Positive features are loaded and clusters are correctly extracted", {
  Positive_Features <- openxlsx::read.xlsx(opt$xlsx)
  
  df <- extract_top_markers(Positive_Features)
  clusters_found <- get_clusters_from_features(Positive_Features)
  
  clusters_test_glb <<- clusters_found  # Assign to global environment
  df_glb <<- df
  
  expect_s3_class(df, "data.frame")
  expect_true("cluster" %in% colnames(df))
  expect_true(length(clusters) > 0)
})

# IMPORTANT - this part uses a temporary pdf, likely to just check that and then deelete is, be wary that maybe the code
# needs actual output maybe.

test_that("PDF plots for each cluster are generated and saved", {
  df <- df_glb
  clusters <- clusters_test_glb
  
  tmp_pdf <- generate_pdf_plot(df, clusters, opt$output.plot)
  
  expect_true(file.exists(tmp_pdf))
  expect_gt(file.info(tmp_pdf)$size, 1000)
  
  # Clean up
  file.remove(tmp_pdf)
})

# Define file path
output_file <- file.path(test_data_dir, "results/scrna-top-marker-plot-test.txt")
# Define content to write
lines_to_write <- c(
  "top marker ran"
)
# Write the content to the file
writeLines(lines_to_write, con = output_file)
