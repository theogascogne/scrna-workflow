#!/usr/bin/env Rscript
library(testthat)
library(optparse)
require(Seurat)
require(tidyverse)

# Get the path to the 'cellsnake' package
cellsnake_path <- system("python -c 'import cellsnake; print(cellsnake.__path__[0])'", intern = TRUE)

# Define the base path to the testData folder
test_data_dir <- file.path(cellsnake_path, "scrna/workflow/tests/testData")

# Define options with dynamically constructed paths
option_list <- list(
  make_option(c("--rds"), type = "character", default = file.path(test_data_dir, "processed/defaultTest/output.rds"), help = "Processed rds file of a Seurat object", metavar = "character"),
  make_option(c("--sampleid"), type = "character", default = "TestSample", help = "Sample ID for testing", metavar = "character"),
  make_option(c("--fplot"), type = "character", default = file.path(test_data_dir, "results/fplot_test.pdf"), help = "nFeature plot", metavar = "character"),
  make_option(c("--cplot"), type = "character", default = file.path(test_data_dir, "results/cplot_test.pdf"), help = "nCount plot", metavar = "character"),
  make_option(c("--mtplot"), type = "character", default = file.path(test_data_dir, "results/mtplot_test.pdf"), help = "Percent MT plot", metavar = "character"),
  make_option(c("--rpplot"), type = "character", default = file.path(test_data_dir, "results/rpplot_test.pdf"), help = "Ribo plot", metavar = "character")
)

opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)
try(
  {
    source(file.path(cellsnake_path,"scrna/workflow/scripts/scrna_workflow_helper_functions/scrna-technicals-functions.R"))
  },
  silent = TRUE
)

# Ensure that the directories for the output files exist
output_dirs <- c(
  dirname(opt$fplot),
  dirname(opt$cplot),
  dirname(opt$mtplot),
  dirname(opt$rpplot),
  dirname(file.path(test_data_dir, "results/scrna-technicals-test.txt"))
)

scrna <- readRDS(opt$rds) 

generate_feature_plots(
  rds_file = opt$rds,
  fplot = opt$fplot,
  cplot = opt$cplot,
  mtplot = opt$mtplot,
  rpplot = opt$rpplot
)

# Ensure that the plots are saved correctly
test_that("nFeature_RNA FeaturePlot generates and saves correctly", {
  expect_true("nFeature_RNA" %in% colnames(scrna[[]]))
  
  p <- FeaturePlot(scrna, features = "nFeature_RNA", pt.size = 0.1, raster = FALSE)
  expect_s3_class(p, "gg")
  
  # Save to the output file
  ggsave(opt$fplot, p, width = 7, height = 5)
  expect_true(file.exists(opt$fplot))
})

test_that("nCount_RNA FeaturePlot generates and saves correctly", {
  expect_true("nCount_RNA" %in% colnames(scrna[[]]))
  
  p <- FeaturePlot(scrna, features = "nCount_RNA", pt.size = 0.1, raster = FALSE)
  expect_s3_class(p, "gg")
  
  # Save to the output file
  ggsave(opt$cplot, p, width = 7, height = 5)
  expect_true(file.exists(opt$cplot))
})

test_that("percent.mt FeaturePlot generates and saves correctly", {
  expect_true("percent.mt" %in% colnames(scrna[[]]))
  
  p <- FeaturePlot(scrna, features = "percent.mt", pt.size = 0.1, raster = FALSE)
  expect_s3_class(p, "gg")
  
  # Save to the output file
  ggsave(opt$mtplot, p, width = 7, height = 5)
  expect_true(file.exists(opt$mtplot))
})

test_that("percent.rp FeaturePlot generates and saves correctly", {
  expect_true("percent.rp" %in% colnames(scrna[[]]))
  
  p <- FeaturePlot(scrna, features = "percent.rp", pt.size = 0.1, raster = FALSE)
  expect_s3_class(p, "gg")
  
  # Save to the output file
  ggsave(opt$rpplot, p, width = 7, height = 5)
  expect_true(file.exists(opt$rpplot))
})

# Define file path for the test log
output_file <- file.path(test_data_dir, "results/scrna_technicals_test.txt")

# Define content to write in the log
lines_to_write <- c("technicals plots ran")

# Write the content to the file
writeLines(lines_to_write, con = output_file)