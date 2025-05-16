#!/usr/bin/env Rscript

library(optparse)
library(testthat)
require(SingleR)
require(pheatmap)
require(Seurat)
require(tidyverse)

# Get the path to the 'cellsnake' package
cellsnake_path <- system("python -c 'import cellsnake; print(cellsnake.__path__[0])'", intern = TRUE)

# Define the base path to the testData folder
test_data_dir <- file.path(cellsnake_path, "scrna/workflow/tests/testData")

# Define test options with dynamically constructed paths
option_list <- list(
  make_option(c("--rds"), type = "character", default = file.path(test_data_dir, "processed/defaultTest/output.rds"), help = "Processed rds file of a Seurat object", metavar = "character"),
  make_option(c("--sheplot"), type = "character", default = file.path(test_data_dir, "fetal-brain/results/10X_17_028/percent_mt~10/resolution~0.8/singler/plot_score_heatmap-seurat_clusters.pdf"), help = "Output score heatmap plot file name", metavar = "character"),
  make_option(c("--sheplottop"), type = "character", default = file.path(test_data_dir, "fetal-brain/results/10X_17_028/percent_mt~10/resolution~0.8/singler/plot_score_heatmap_top-seurat_clusters.pdf"), help = "Output score heatmap plot file name, top 20", metavar = "character"),
  make_option(c("--pheplot"), type = "character", default = file.path(test_data_dir, "fetal-brain/results/10X_17_028/percent_mt~10/resolution~0.8/singler/plot_clusters-seurat_clusters.pdf"), help = "Output heatmap plot file name", metavar = "character"),
  make_option(c("--idents"), type = "character", default = "seurat_clusters", help = "Meta data column name", metavar = "character"),
  make_option(c("--csv"), type = "character", default = NULL, help = "A meta data table", metavar = "character"),
  make_option(c("--prediction"), type = "character", default = file.path(test_data_dir, "/singler/defaultTest/annotation.rds"), help = "Input prediction file", metavar = "character"),
  make_option(c("--xlsx"), type = "character", default = file.path(test_data_dir, "results/table_annotations_per-seurat_clusters.xlsx"), help = "Input prediction file", metavar = "character")
)

opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

try(
  {
    source(file.path(cellsnake_path,"scrna/workflow/scripts/scrna_workflow_helper_functions/scrna-singler-plots-functions.R"))
  },
  silent = TRUE
)

pred <- readRDS(opt$prediction)
scrna <- readRDS(file = opt$rds)
DefaultAssay(scrna) <- "RNA"

test_that("Score heatmap is plotted and saved correctly", {
  pred <- readRDS(opt$prediction)
  # Create a temporary file for saving the plot
  tmp_file <- tempfile(fileext = ".pdf")
  
  save_score_heatmap(pred, tmp_file)
  
  # Confirm the file exists and is not empty
  expect_true(file.exists(tmp_file))
  expect_gt(file.info(tmp_file)$size, 1000)
  
  # Clean up by removing the temporary file
  file.remove(tmp_file)
})

test_that("Top score heatmap is plotted and saved correctly", {
  # Create a temporary file for saving the plot
  tmp_file <- tempfile(fileext = ".pdf")
  pred <- readRDS(opt$prediction)
  
  # Call the function with custom parameters
  save_score_heatmap(pred, tmp_file, width = 7, height = 4, show_labels = FALSE, max_labels = 20)
  
  # Confirm the file exists and is not empty
  expect_true(file.exists(tmp_file))
  expect_gt(file.info(tmp_file)$size, 1000)
  
  # Clean up by removing the temporary file
  file.remove(tmp_file)
})

test_that("Heatmap is generated and saved correctly", {
  # Read in the data
  pred <- readRDS(opt$prediction)
  
  # Create the table
  tab <- table(Assigned = pred$pruned.labels, Cluster = scrna@meta.data[[opt$idents]])
  
  # Save the heatmap plot
  tmp_file <- tempfile(fileext = ".pdf")
  save_heatmap_plot(scrna, pred, tab, tmp_file)  # Use the functionalized heatmap save function
  
  # Check that the heatmap was saved successfully
  expect_true(file.exists(tmp_file))
  expect_gt(file.info(tmp_file)$size, 1000)
  
  # Clean up the temporary plot file
  file.remove(tmp_file)
})

test_that("Table and Excel file are saved correctly", {
  # Read in the data
  pred <- readRDS(opt$prediction)
  
  # Create the table
  tab <- table(Assigned = pred$pruned.labels, Cluster = scrna@meta.data[[opt$idents]])
  
  # Save the table to Excel
  save_table_to_excel(tab, opt$xlsx)  # Use the functionalized Excel save function
  
  # Check that the Excel file was saved successfully
  expect_true(file.exists(opt$xlsx))
  expect_gt(file.info(opt$xlsx)$size, 1000)
})

# remove global values after test
rm(pred)

# Define file path
output_file <- file.path(test_data_dir, "results/scrna_singler_plots_test.txt")
# Define content to write
lines_to_write <- c(
  "singler plots ran"
)
# Write the content to the file
writeLines(lines_to_write, con = output_file)

