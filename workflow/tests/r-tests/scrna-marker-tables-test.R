#!/usr/bin/env Rscript

library(optparse)
require(tidyverse)
library(openxlsx)
library(testthat)

# Get the path to the 'cellsnake' package
cellsnake_path <- system("python -c 'import cellsnake; print(cellsnake.__path__[0])'", intern = TRUE)
test_data_dir <- file.path(cellsnake_path, "scrna/workflow/tests/testData")

# Define options
option_list <- list(
  make_option(c("--rds"), type = "character", default = file.path(test_data_dir, "analyses/defaultTest/seurat_clusters.rds"),
              help = "RDS file of marker data frame", metavar = "character"),
  make_option(c("--output.xlsx.positive"), type = "character", default = file.path(test_data_dir, "results/defaultTest/table_positive-markers-seurat_clusters.xlsx"),
              help = "Excel table of positive markers", metavar = "character"),
  make_option(c("--output.xlsx.all"), type = "character", default = file.path(test_data_dir, "results/defaultTest/table_all-markers-seurat_clusters.xlsx"),
              help = "Excel table of all markers", metavar = "character")
)

opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

try(
  {
    source(file.path(cellsnake_path,"scrna/workflow/scripts/scrna_workflow_helper_functions/scrna-marker-tables-functions.R"))
  },
  silent = TRUE
)

# Load markers
all_markers <- readRDS(file = opt$rds)


# Optional: Tests to confirm files were written correctly
test_that("Excel files are written correctly", {
  write_marker_excel_outputs(
    markers_df = all_markers,
    all_path = opt$output.xlsx.all,
    positive_path = opt$output.xlsx.positive
  )
  
  expect_true(file.exists(opt$output.xlsx.all))
  expect_true(file.exists(opt$output.xlsx.positive))
  
  all_check <- read.xlsx(opt$output.xlsx.all)
  pos_check <- read.xlsx(opt$output.xlsx.positive)
  
  expect_equal(nrow(all_check), nrow(all_markers))
  expect_equal(nrow(pos_check), nrow(filter(all_markers, avg_log2FC > 0)))
})
