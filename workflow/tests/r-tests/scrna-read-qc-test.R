#!/usr/bin/env Rscript

library(optparse)
library(testthat)
library(Seurat)
require(tidyverse)
require(patchwork)
require(tools)
require(data.table)

# Get the path to the 'cellsnake' package
cellsnake_path <- system("python -c 'import cellsnake; print(cellsnake.__path__[0])'", intern = TRUE)

# Define the base path to the testData folder
test_data_dir <- file.path(cellsnake_path, "scrna/workflow/tests/testData")

# Check if the directory exists, and if not, stop execution
if (!dir.exists(file.path(test_data_dir, "/data"))) {
  stop("Data directory does not exist: ", file.path(test_data_dir, "/data"))
}

# Define options with default values as specified
option_list <- list(
  make_option(c("--data.dir"), type = "character", default = file.path(test_data_dir, "data/testSample/outs/filtered_feature_bc_matrix"), help = "Data directory", metavar = "character"),
  make_option(c("--output.rds"), type = "character", default = file.path(test_data_dir, "raw/defaultTest/output.rds"), help = "Output RDS file", metavar = "character"),
  make_option(c("--sampleid"), type = "character", default = "10X_17_028", help = "Sample ID", metavar = "character"),
  make_option(c("--percent.mt"), type = "character", default = "10", help = "Max mitochondrial gene percentage", metavar = "character"),
  make_option(c("--percent.rp"), type = "double", default = 0, help = "Min ribosomal gene percentage", metavar = "double"),
  make_option(c("--min.features"), type = "integer", default = 200, help = "Min features", metavar = "integer"),
  make_option(c("--max.features"), type = "character", default = "Inf", help = "Max features", metavar = "character"),
  make_option(c("--min.molecules"), type = "integer", default = 0, help = "Min molecules", metavar = "integer"),
  make_option(c("--max.molecules"), type = "character", default = "Inf", help = "Max molecules", metavar = "character"),
  make_option(c("--min.cells"), type = "integer", default = 3, help = "Min cells", metavar = "integer"),
  make_option(c("--before.violin.plot"), type = "character", default = file.path(test_data_dir, "results/defaultTest/technicals/plot_before-qc-trimming.pdf"), help = "Violin plot before QC trimming", metavar = "character"),
  make_option(c("--after.violin.plot"), type = "character", default = file.path(test_data_dir, "results/defaultTest/technicals/plot_after-qc-trimming.pdf"), help = "Violin plot after QC trimming", metavar = "character"),
  make_option(c("--plot.mtplot"), type = "character", default = file.path(test_data_dir, "results/defaultTest/technicals/plot_model-metrics-mitochondrial-genes.pdf"), help = "Plot for mitochondrial gene metrics", metavar = "character")
)

opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

try(
  {
    source(file.path(cellsnake_path,"scrna/workflow/scripts/scrna_workflow_helper_functions/scrna-read-qc-functions.R"))
  },
  silent = TRUE
)

try(
  {
    source(file.path(cellsnake_path,"scrna/workflow/scripts/generate_10X_data.R"))
  },
  silent = TRUE
)

test_that("ensure_dir_exists behaves correctly in both existing and non-existing cases", {
  # 1. Prepare a temp directory and file
  temp_dir <- tempfile("existing_dir_")
  dir.create(temp_dir, recursive = TRUE)
  test_file_existing <- file.path(temp_dir, "file.txt")
  
  # Case 1: Directory already exists
  expect_silent(ensure_dir_exists(test_file_existing))
  expect_true(dir.exists(dirname(test_file_existing)))
  
  # Cleanup
  unlink(temp_dir, recursive = TRUE)
  
  # 2. Create path for non-existent directory
  temp_nonexistent <- tempfile("new_dir_")
  test_file_new <- file.path(temp_nonexistent, "file.txt")
  
  # Confirm the directory doesn't exist
  expect_false(dir.exists(dirname(test_file_new)))
  
  # Case 2: Directory does not exist (should be created)
  expect_message(ensure_dir_exists(test_file_new), regexp = "Creating missing directory:")
  expect_true(dir.exists(dirname(test_file_new)))
  
  # Final cleanup
  unlink(temp_nonexistent, recursive = TRUE)
})


test_that("function_read_input correctly reads a .gz file and returns a dgCMatrix", {
  result <- function_read_input(opt)
  
  # Check if the result is of class 'dgCMatrix' using is()
  expect_true(is(result, "dgCMatrix"), info = "The returned result is not of class 'dgCMatrix'.")
  
  seurat_obj <<- CreateSeuratObject(
    counts = result,
    project = make.names(opt$sampleid),
    min.cells = opt$min.cells,
    min.features = opt$min.features
  )
})

test_that("CreateSeuratObject correctly creates seuratObjects from valid dgCMatrix inputs", {
  # --- Happy path: valid matrix ---
  valid_matrix <- subset_and_inject_mt_rp_real_matrix()
  seurat_obj_valid <- CreateSeuratObject(
    counts = valid_matrix$mat,
    project = make.names(opt$sampleid),
    min.cells = opt$min.cells,
    min.features = opt$min.features
  )
  
  expect_true(inherits(seurat_obj_valid, "Seurat"), info = "The result is not a valid Seurat object.")
  rm(valid_matrix)
})


test_that("QC metrics and violin plot step works after Seurat object creation", {
  skip_if_not(exists("seurat_obj"), "Seurat object not found from previous test.")
  scrna <- prepare_seurat_object_for_qc(seurat_obj, opt)
  
  expect_true(all(c("percent.mt", "percent.rp") %in% colnames(scrna@meta.data)),
              info = "Expected QC metadata fields were not added")
  
  expect_silent({
    VlnPlot(scrna, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rp"))
  })
  
  scrna_next <<- scrna
})

test_that("Cell filtering, miQC fallback, and plotting runs as expected", {
  skip_if_not(exists("scrna_next"), "scrna_next not available from previous step.")
  
  scrna <- scrna_next
  original_cells <- ncol(scrna)
  
  scrna <- filter_cells_by_qc_thresholds(scrna, opt)
  expect_lte(ncol(scrna), original_cells)
  
  temp_plot_file <- tempfile(fileext = ".pdf")
  scrna <- if (tolower(opt$percent.mt) == "auto") {
    run_miQC_or_fallback(scrna, opt, temp_plot_file)
  } else {
    subset(scrna, subset = percent.mt <= as.numeric(opt$percent.mt))
  }
  
  expect_lte(ncol(scrna), original_cells)
  
  # RP filter
  scrna <- subset(scrna, subset = percent.rp >= opt$percent.rp)
  expect_lte(ncol(scrna), original_cells)
  
  # Violin plot and save
  temp_violin_file <- tempfile(fileext = ".pdf")
  expect_silent(plot_qc_violin(scrna, temp_violin_file))
  
  unlink(c(temp_plot_file, temp_violin_file))
  dir.create(dirname(opt$output.rds), recursive = TRUE, showWarnings = FALSE)
  saveRDS(scrna, opt$output.rds)
})

test_that("run_miQC_or_fallback runs when percent.mt is set to auto", {
  skip_if_not(exists("scrna_next"), "scrna_next not available from previous step.")
  
  scrna <- scrna_next
  opt$percent.mt <- "auto"
  temp_plot_file <- tempfile(fileext = ".pdf")
  
  result <- suppressWarnings(run_miQC_or_fallback(scrna, opt, temp_plot_file))
  
  expect_true(inherits(result, "Seurat"))
  expect_true(file.exists(temp_plot_file))
  
  unlink(temp_plot_file)
})

test_that("run_miQC_or_fallback runs without issues when percent.mt is 0 for all cells", {
  skip_if_not(exists("scrna_next"), "scrna_next not available from previous step.")
  
  scrna <- scrna_next
  
  # Set percent.mt to 0 for all cells (to simulate the case with no mitochondrial genes)
  scrna$percent.mt <- rep(0, ncol(scrna))
  
  opt$percent.mt <- "auto"
  temp_plot_file <- tempfile(fileext = ".pdf")
  
  # Run the function, which should skip the QC process since percent.mt == 0
  result <- suppressWarnings(run_miQC_or_fallback(scrna, opt, temp_plot_file))
  
  # Ensure the result is still a Seurat object and no plot is generated (fallback skips plotting)
  expect_true(inherits(result, "Seurat"))
  expect_false(file.exists(temp_plot_file))  # No plot should be created
  
  unlink(temp_plot_file)
})
