#!/usr/bin/env Rscript
library(testthat)
library(optparse)
require(SingleR)
require(SingleCellExperiment)
require(Seurat)

# Get the path to the 'cellsnake' package
cellsnake_path <- system("python -c 'import cellsnake; print(cellsnake.__path__[0])'", intern = TRUE)

# Define the base path to the testData folder
test_data_dir <- file.path(cellsnake_path, "scrna/workflow/tests/testData")

# Define test options with dynamically constructed paths
option_list <- list(
  make_option(c("--rds"), type = "character", default = file.path(test_data_dir, "processed/defaultTest/output.rds"), help = "Processed rds file of a Seurat object", metavar = "character"),
  make_option(c("--output"), type = "character", default = file.path(test_data_dir, "singler/defaultTest/annotation.rds"), help = "Output prediction file", metavar = "character"),
  make_option(c("--reference"), type = "character", default = "HumanPrimaryCellAtlasData", help = "SingleR reference", metavar = "character"),
  make_option(c("--granulation"), type = "character", default = "label.main", help = "SingleR granulation level", metavar = "character")
)

# parse options
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

try(
  {
    source(file.path(cellsnake_path,"scrna/workflow/scripts/scrna_workflow_helper_functions/scrna-singler-annotation-functions.R"))
  },
  silent = TRUE
)

# Load initial scrna object
scrna <- readRDS(opt$rds)  # Replace with actual file path

test_that("SingleR celltype annotation runs correctly", {

  suppressWarnings({
    # Run the function and check that the output file is created
    run_singleR_prediction(
      rds_file = opt$rds, 
      output_file = opt$output, 
      reference = opt$reference, 
      granulation = opt$granulation
    )
  })
  expect_true(file.exists(opt$output))
  
  # Load the prediction to check its type
  pred <- readRDS(opt$output)
  expect_true(inherits(pred, "DFrame"))
})
