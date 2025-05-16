#!/usr/bin/env Rscript

library(optparse)
require(Seurat)
require(tidyverse)
library(testthat)

# Get the path to the 'cellsnake' package
cellsnake_path <- system("python -c 'import cellsnake; print(cellsnake.__path__[0])'", intern = TRUE)

# Define the base path to the testData folder
test_data_dir <- file.path(cellsnake_path, "scrna/workflow/tests/testData")

# Define options
option_list <- list(
  make_option(c("--rds"), type = "character", 
              default = file.path(test_data_dir, "processed/defaultTest/output.rds"),
              help = "Processed RDS file of a Seurat object", metavar = "character"),
  make_option(c("--logfc.threshold"), type = "double", default = 0.25,
              help = "Log fold change threshold [default= %default]", metavar = "double"),
  make_option(c("--test.use"), type = "character", default = "wilcox",
              help = "Statistical test to use [default= %default]", metavar = "character"),
  make_option(c("--output.rds"), type = "character", 
              default = file.path(test_data_dir, "analyses/defaultTest/seurat_clusters.rds"),
              help = "Output RDS file containing marker genes", metavar = "character"),
  make_option(c("--idents"), type = "character", default = "seurat_clusters",
              help = "Metadata column name for marker analysis", metavar = "character")
)

# parse options
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

try(
  {
    source(file.path(cellsnake_path,"scrna/workflow/scripts/scrna_workflow_helper_functions/scrna-find-markers-functions.R"))
  },
  silent = TRUE
)

# Load initial scrna object
scrna <- readRDS(opt$rds)  # Replace with actual file path

test_that("run_and_save_markers() executes and saves RDS", {
  run_and_save_markers(
    scrna = scrna,
    idents_col = opt$idents,
    logfc.threshold = opt$logfc.threshold,
    test.use = opt$test.use,
    output_path = opt$output.rds
  )
  
  # Confirm output file was written
  expect_true(file.exists(opt$output.rds))
  expect_gt(file.info(opt$output.rds)$size, 1000)
  
  # Load markers for subsequent tests
  markers <<- readRDS(opt$output.rds)
})

test_that("Markers dataframe structure is valid", {
  expect_s3_class(markers, "data.frame")
  expect_true(all(c("gene", "avg_log2FC", "p_val", "cluster") %in% colnames(markers)))
})

test_that("Identities are set correctly", {
  # Reset identities to ensure test is independent
  Idents(scrna) <- scrna@meta.data[[opt$idents]]
  expect_true(all(Idents(scrna) == scrna@meta.data[[opt$idents]]))
})

