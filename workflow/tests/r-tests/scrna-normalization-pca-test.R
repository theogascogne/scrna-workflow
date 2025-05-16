#!/usr/bin/env Rscript

library(testthat)
library(Seurat)
library(tidyverse)
library(optparse)
library(celldex)
library(SingleR)

# Get the path to the 'cellsnake' package
cellsnake_path <- system("python -c 'import cellsnake; print(cellsnake.__path__[0])'", intern = TRUE)

# Define the base path to the testData folder
test_data_dir <- file.path(cellsnake_path, "scrna/workflow/tests/testData")

# Define options with dynamically constructed paths
option_list <- list(
  make_option(c("--scale.factor"), type = "integer", default = 10000, help = "Scale factor [default= %default]", metavar = "character"),
  make_option(c("--nfeatures"), type = "integer", default = 2000, help = "Highly variable features [default= %default]", metavar = "integer"),
  make_option(c("--variable.selection.method"), type = "character", default = "vst", help = "Variable selection method [default= %default]", metavar = "character"),
  make_option(c("--rds"), type = "character", default = file.path(test_data_dir, "raw/defaultTest/output.rds"), help = "Test RDS file path", metavar = "character"),
  make_option(c("--normalization.method"), type = "character", default = "LogNormalize", help = "Normalization method [default= %default]", metavar = "character"),
  make_option(c("--doublet.filter"), action = "store_true", default = FALSE, help = "Test Doublet filter [default= %default]"),
  make_option(c("--integration"), action = "store_true", default = FALSE, help = "Run with integration [default= %default]"),
  make_option(c("--umap"), action = "store_true", default = TRUE, help = "Test UMAP generation [default= %default]"),
  make_option(c("--tsne"), action = "store_true", default = TRUE, help = "Test t-SNE generation [default= %default]"),
  make_option(c("--resolution"), type = "character", default = "0.8", help = "Test resolution [default= %default]", metavar = "character"),
  make_option(c("--output.rds"), type = "character", default = file.path(test_data_dir, "processed/defaultTest/output.rds"), help = "Output RDS test file [default= %default]", metavar = "character"),
  make_option(c("--cpu"), type = "integer", default = 2, help = "Number of CPUs for test [default= %default]", metavar = "character"),
  make_option(c("--reference"), type = "character", default = "HumanPrimaryCellAtlasData", help = "SingleR reference for annotation [default= %default]", metavar = "character")
)

opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

try(
  {
    source("~/miniconda3/envs/test/lib/python3.9/site-packages/cellsnake/scrna/workflow/scripts/scrna-functions.R")
  },
  silent = TRUE
)
try(
  {
    source(paste0(system("python -c 'import os; import cellsnake; print(os.path.dirname(cellsnake.__file__))'", intern = TRUE), "/scrna/workflow/scripts/scrna-functions.R"))
  },
  silent = TRUE
)

try(
  {
    source(file.path(cellsnake_path,"scrna/workflow/scripts/scrna_workflow_helper_functions/scrna-normalization-pca-functions.R"))
  },
  silent = TRUE
)


# Load initial scrna object
scrna <- readRDS(file = opt$rds)  # Replace with actual file path

# --- Test 1: Normalization ---
test_that("Normalization works correctly", {
  scrna_test <- scrna  # Clone the scrna object to avoid mutating the original
  
  # Perform normalization
  scrna_test <- normalize_and_select_features(scrna_test, opt)
  
  # Check that the normalized data does not contain NA values
  expect_false(any(is.na(scrna_test@assays$RNA$data@x)))  
  
  # Check that all values in the normalized data are numeric
  expect_true(all(is.numeric(scrna_test@assays$RNA$data@x)))
})

# --- Test 2: Variable Feature Selection ---
test_that("Variable feature selection works correctly", {
  scrna_test <- scrna 
  
  scrna_test <- normalize_and_select_features(scrna_test, opt)
  
  # Check if variable features are selected correctly
  expect_equal(length(VariableFeatures(scrna_test)), opt$nfeatures)
  
  scrna <<- scrna_test
})

# --- Test 3: Scaling Data ---
test_that("Scaling works correctly", {
  scrna_test <- scrna
  features_to_scale <- VariableFeatures(scrna_test)
  
  scrna_test <- ScaleData(scrna_test, features = features_to_scale)
  
  scaled_data <- LayerData(scrna_test, assay = "RNA", layer = "scale.data")
  
  # Check it exists and has valid values
  expect_true(!is.null(scaled_data))
  expect_true(all(!is.na(scaled_data)))
  expect_gt(sum(abs(scaled_data)), 0)
  
  # Check basic statistics
  gene_means <- rowMeans(scaled_data)
  gene_sds <- apply(scaled_data, 1, sd)
  
  cat("Mean range:", range(gene_means), "\n")
  cat("SD range:  ", range(gene_sds), "\n")
  
  expect_true(all(abs(gene_means) < 0.05))
  expect_true(all(abs(gene_sds - 1) < 0.7))
  
  #scrna <<- scrna_test
})

# --- Test 4: PCA ---
test_that("PCA works correctly", {
  scrna_test <- scrna
  
  # Use the updated PCA function, assuming it is now handled by the helper function `run_pca_pipeline`
  pca_out <- run_pca_pipeline(scrna_test)
  
  scrna_test <- pca_out$scrna
  
  # Check if PCA was run by inspecting PCA results
  expect_true("pca" %in% names(scrna_test@reductions))
  expect_true(nrow(scrna_test@reductions$pca@cell.embeddings) > 0)
  
  scrna <<- scrna_test
})

# --- Test 5: FindNeighbors ---
test_that("FindNeighbors works correctly", {
  scrna_test <- scrna #helper function has already ran findNeighbors on 'scrna' in test 4. checking if it was done correctly!
  
  #dimensionReduction <- function_pca_dimensions(scrna_test)
  #scrna_test <- FindNeighbors(scrna_test, dims = 1:dimensionReduction)
  
  graph_names <- names(scrna_test@graphs)
  
  # Expect both KNN and SNN graphs to exist - please check if this naming is indeed consistent between sample runs.
  expect_true("RNA_nn" %in% graph_names)
  expect_true("RNA_snn" %in% graph_names)
  
  # Optionally inspect the structure (for debugging)
  message("Neighbor graphs found: ", paste(graph_names, collapse = ", "))
  
  # Save for next test
  scrna <<- scrna_test
})

# --- Test 6: Clustering with 0.8 (see optparse) resolution ---
test_that("Clustering with auto_select_resolution works correctly", {
  scrna_test <- scrna
  
  # Apply auto_select_resolution function to get the optimal resolution
  scrna_test <- retrieve_clustering(scrna_test, opt, pca_out$dims)
  
  # Check if clustering was applied correctly
  expect_true("seurat_clusters" %in% colnames(scrna_test@meta.data))
  expect_true(length(unique(scrna_test$seurat_clusters)) > 1)
  
  # Save the state for the next test
  scrna <<- scrna_test
})

# --- Test 6: Clustering with auto resolution ---   does not work because MultiKParallel is not available!
#test_that("Clustering with auto_select_resolution works correctly", {
#  scrna_test <- scrna
#  opt$resolution <- "auto"
  # Apply auto_select_resolution function to get the optimal resolution
#  scrna_test <- retrieve_clustering(scrna_test, opt, pca_out$dims)
  
  # Check if clustering was applied correctly
#  expect_true("seurat_clusters" %in% colnames(scrna_test@meta.data))
#  expect_true(length(unique(scrna_test$seurat_clusters)) > 1)
  
  # Save the state for the next test
#  scrna <<- scrna_test
#})

# --- Test 8: UMAP ---
test_that("UMAP works correctly", {
  scrna_test <- scrna
  
  suppressWarnings({
    scrna_test <- run_TSNE_UMAP(scrna_test, opt, dims = 1:10)
  })
  
  # Check if UMAP was computed
  expect_true("umap" %in% names(scrna_test@reductions))
  
  # Save the state for the next test
  scrna <<- scrna_test
})

# --- Test 9: tSNE ---
test_that("tSNE works correctly", {
  scrna_test <- scrna
  
  scrna_test <- run_TSNE_UMAP(scrna_test, opt, dims = 1:10)
  
  # Check if tSNE was computed
  expect_true("tsne" %in% names(scrna_test@reductions))
  
  # Save the state for the next test
  scrna <<- scrna_test
})


# --- Test 7: Doublet filtering with filter_doublets ---
opt$doublet.filter <- TRUE
test_that("DoubletFinder runs and filters correctly", {
  scrna_test <- scrna
  
  if (opt$doublet.filter) {
    scrna_test <- run_doublet_filter(scrna_test, opt)
    
    # Check filtering worked
    expect_true("DoubletFinder" %in% colnames(scrna_test@meta.data))
    expect_true(all(scrna_test$DoubletFinder == "Singlet"))
  }
  
  scrna <<- scrna_test
})


# --- Test 10: SingleR annotation works ---
test_that("annotate_with_singleR adds pruned.labels and sets SingleRref", {
  scrna_test <- scrna
  
  # Call annotation helper
  suppressWarnings({
    scrna_test <- annotate_with_singleR(scrna_test, opt$reference)
  })
  
  # Check metadata includes the SingleR labels
  expect_true("singler" %in% colnames(scrna_test@meta.data))
  
  # Check attribute is set
  expect_equal(attr(scrna_test, "SingleRref"), opt$reference)
  
  scrna <<- scrna_test
})

saveRDS(scrna, file = opt$output.rds)