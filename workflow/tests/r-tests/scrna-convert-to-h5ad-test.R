#!/usr/bin/env Rscript
library(testthat)
library(optparse)
require(Seurat)
require(SeuratDisk)
require(tidyverse)
require(scCustomize)

# Get the path to the 'cellsnake' package
cellsnake_path <- system("python -c 'import cellsnake; print(cellsnake.__path__[0])'", intern = TRUE)

# Define the base path to the testData folder
test_data_dir <- file.path(cellsnake_path, "scrna/workflow/tests/testData")

# Define options
option_list <- list(
  make_option(c("--rds"), type = "character", 
              default = file.path(test_data_dir, "processed/defaultTest/output.rds"),
              help = "Seurat object RDS file", metavar = "character"),
  make_option(c("--output"), type = "character", 
              default = file.path(test_data_dir, "analyses/defaultTest/testSample.h5ad"),
              help = "Output h5ad file", metavar = "character")
)

# Parse options
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

try(
  {
    source(file.path(cellsnake_path,"scrna/workflow/scripts/scrna_workflow_helper_functions/scrna-convert-to-h5ad-functions.R"))
  },
  silent = TRUE
)


scrna <<- NULL
original_scrna <<- NULL
tmp_h5seurat <<- NULL
tmp_h5ad <<- NULL

output_file_name <<- NULL

# Test 1: Read Seurat object and check its class
test_that("Seurat object is read correctly", {
  scrna <<- readRDS(file = opt$rds)
  original_scrna <<- scrna
  
  # Check that the object is a Seurat object
  expect_true(inherits(scrna, "Seurat"))
})

# Test 2: Convert Assay to Seurat v3/4 format
test_that("Seurat assay is converted to Seurat v3/4 format", {
  scrna <<- convert_assay(scrna, version = "V3")
  
  # Check if the conversion was successful
  expect_true("RNA" %in% names(scrna@assays))
})

# Test 5: Diet Seurat to reduce memory usage
test_that("Seurat object is dieted", {
  suppressWarnings({
    scrna <<- UpdateSeuratObject(scrna)
  })
  DefaultAssay(scrna) <<- "RNA"
  scrna <<- DietSeurat(scrna)
  
  # Check if the object size decreased
  expect_true(object.size(scrna) < object.size(original_scrna))
})

# Test 6: Convert metadata factors to characters
test_that("Metadata factors are converted to characters", {
  scrna <<- convert_factors_to_characters(scrna)
  
  # Ensure that all factors in metadata have been converted to characters
  expect_true(all(sapply(scrna@meta.data, function(col) !is.factor(col))))
})

# Test 5: Save Seurat object as .h5Seurat
test_that("Seurat object is saved as .h5Seurat", {
  output_file_name <<- str_remove_all(opt$output, ".h5ad$")
  tmp_h5seurat <<- tempfile(fileext = ".h5Seurat")
  
  # Save the Seurat object as .h5Seurat
  save_h5seurat(scrna, tmp_h5seurat)
  
  # Check if the .h5Seurat file is created
  expect_true(file.exists(tmp_h5seurat))
})

# Test 6: Convert .h5Seurat to .h5ad format
test_that("Convert .h5Seurat to .h5ad", {
  # Convert .h5Seurat to .h5ad directly at opt$output
  output_dir <- dirname(opt$output)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)  # Create the directory if it doesn't exist
  }
  
  # Convert the Seurat object from .h5Seurat to .h5ad
  convert_to_h5ad(tmp_h5seurat, opt$output)
  
  # Check if the final .h5ad file is created at the specified location
  expect_true(file.exists(opt$output))
})

# Test 9: Verify Seurat object processing (Assay type check)
test_that("Seurat object is processed correctly", {
  expect_equal(DefaultAssay(scrna), "RNA")
})

# Test 10: Check if metadata factors are converted to characters
test_that("Metadata factors are converted correctly", {
  expect_true(all(sapply(scrna@meta.data, class) != "factor"))
})

# Clean up - remove the temporary files after the tests
test_that("Clean up temporary files", {
  file.remove(tmp_h5seurat)
})