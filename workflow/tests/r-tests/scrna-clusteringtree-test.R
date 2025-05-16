#!/usr/bin/env Rscript
library(testthat)
require(tidyverse)
require(optparse)
require(Seurat)
require(clustree)

# Get the path to the 'cellsnake' package
cellsnake_path <- system("python -c 'import cellsnake; print(cellsnake.__path__[0])'", intern = TRUE)

# Define the base path to the testData folder
test_data_dir <- file.path(cellsnake_path, "scrna/workflow/tests/testData")

# Define options with default test values
option_list <- list(
  make_option(c("--scale.factor"), type = "integer", default = 10000, help = "Scale factor [default= %default]", metavar = "character"),
  make_option(c("--nfeatures"), type = "integer", default = 2000, help = "Highly variable features [default= %default]", metavar = "integer"),
  make_option(c("--variable.selection.method"), type = "character", default = "vst", help = "Variable selection method [default= %default]", metavar = "character"),
  make_option(c("--rds"), type = "character", 
              default = file.path(test_data_dir, "raw/defaultTest/output.rds"), 
              help = "Test RDS file path", metavar = "character"),
  make_option(c("--normalization.method"), type = "character", default = "LogNormalize", help = "Normalization method [default= %default]", metavar = "character"),
  make_option(c("--integration"), action = "store_true", default = FALSE, help = "Run with integration [default= %default]"),
  make_option(c("--clplot"), type = "character", 
              default = file.path(test_data_dir, "results/clustree_test.pdf"), 
              help = "Test clustree output path", metavar = "character"),
  make_option(c("--jeplot"), type = "character", 
              default = file.path(test_data_dir, "results/jackandelbow_test.pdf"), 
              help = "Test Jack and Elbow output path", metavar = "character"),
  make_option(c("--hvfplot"), type = "character", 
              default = file.path(test_data_dir, "results/variable_features_test.pdf"), 
              help = "Test HVF output path", metavar = "character"),
  make_option(c("--heplot"), type = "character", 
              default = file.path(test_data_dir, "results/dimheatmap_test.pdf"), 
              help = "Test heatmap output path", metavar = "character")
)

# parse options
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

try(
  {
    source(file.path(cellsnake_path,"scrna/workflow/scripts/scrna-functions.R"))
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
    source(file.path(cellsnake_path,"scrna/workflow/scripts/scrna_workflow_helper_functions/scrna-clusteringtree-functions.R"))
  },
  silent = TRUE
)

# Load initial scrna object
scrna <- readRDS(opt$rds) # Replace with actual file path


test_that("handle_normalization_or_integration sets integrated assay when integration is TRUE", {
  # Clone the scrna object
  scrna_test <- scrna
  
  suppressWarnings({
    scrna_test[["integrated"]] <- scrna_test[["RNA"]]
  })
  
  # Set the default to RNA first, before switching to integrated
  DefaultAssay(scrna_test) <- "RNA"
  
  # Modify options to simulate integration = TRUE
  opt_integration <- opt
  opt_integration$integration <- TRUE
  
  # Run the helper function
  scrna_test <- handle_normalization_or_integration(scrna_test, opt_integration)
  
  # Assert that the default assay is now "integrated"
  expect_equal(DefaultAssay(scrna_test), "integrated")
  suppressWarnings({
    if ("RNA" %in% names(scrna_test@assays)) {
      
      expect_true(length(scrna_test@assays$RNA$var.features) == 0 || 
                    inherits(scrna_test@assays$RNA$var.features, "dgCMatrix") && 
                    length(scrna_test@assays$RNA$var.features@i) == 0)
    }
  })
})

# --- Test 1: Normalization ---
test_that("Normalization works correctly", {
  scrna_test <- scrna  # Clone the scrna object to avoid mutating the original
  
  # Perform normalization
  scrna_test <- handle_normalization_or_integration(scrna_test, opt)
  
  # Check that the normalized data does not contain NA values
  expect_false(any(is.na(scrna_test@assays$RNA$data@x)))  # @x stores the non-zero values in the sparse matrix
  
  # Check that all values in the normalized data are numeric
  expect_true(all(is.numeric(scrna_test@assays$RNA$data@x)))
  
  scrna <<- scrna_test
})

# --- Test 2: Variable Feature Selection ---
test_that("Variable feature selection works correctly", {
  scrna_test <- scrna 
  
  #the following script was already performed on the scrna object, we do not need to call it here but does anyways for clairty
  #scrna_test <- handle_normalization_or_integration(scrna_test, opt)
  
  # Check if variable features are selected correctly
  expect_equal(length(VariableFeatures(scrna_test)), opt$nfeatures)
  
  scrna <<- scrna_test
})

# --- Test 3: Scaling Data ---
test_that("Scaling works correctly", {
  scrna_test <- scrna
  result <- perform_scaling_and_pca(scrna_test)
  scrna_test <- result$scrna
  
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
  
  expect_true(all(abs(gene_means) < 0.05))   # note: this may be abit too lenient on checking if scaling worked but for unit test we 
  expect_true(all(abs(gene_sds - 1) < 0.7))  # may only want to see if scaling works and expect test data to have extreme values
  
  scrna <<- scrna_test
})

# --- Test 4: PCA ---
test_that("PCA works correctly", {
  scrna_test <- scrna
  
  #this test just checks the PCA that was performed with the helper function
  #result <- perform_scaling_and_pca(scrna_test)
  
  # Check if PCA was run by inspecting PCA results
  expect_true("pca" %in% names(scrna_test@reductions))
  expect_true(nrow(scrna_test@reductions$pca@cell.embeddings) > 0)
  
  scrna <<- scrna_test
})

# --- Test 5: Plotting & JackStraw analysis when integration is FALSE ---
test_that("Variable feature plotting and JackStraw analysis run correctly", {
  if (isFALSE(opt$integration)) {
    scrna_test <- scrna
    top10 <- head(VariableFeatures(scrna_test), 10)
    expect_equal(length(top10), 10)
    
    plot_variable_features(scrna_test, top10, opt$hvfplot, opt$heplot)
    expect_true(file.exists(opt$hvfplot))
    expect_true(file.exists(opt$heplot))
    
    suppressWarnings({
      scrna_test <- plot_jackstraw_and_elbow(scrna_test, opt$jeplot)
    })
    expect_true(file.exists(opt$jeplot))
    
    expect_s3_class(JackStrawPlot(scrna_test, dims = 1:50), "gg")
    expect_s3_class(ElbowPlot(scrna_test, ndims = 50), "gg")
    
    scrna <<- scrna_test
  }
})

test_that("Resolution values are assigned correctly when integration is FALSE", {
  resolution <- get_resolution_range(integration = FALSE)
  
  # Check the length of the resolution vector
  expect_equal(length(resolution), 25)
  
  # Check the first and last values of the resolution vector
  expect_equal(resolution[1], 0.1)
  expect_equal(resolution[length(resolution)], 2.5)
})


test_that("Resolution values are assigned correctly when integration is TRUE", {
  resolution <- get_resolution_range(integration = TRUE)
  
  # Check the length of the resolution vector
  expect_equal(length(resolution), 15)
  
  # Check the first and last values of the resolution vector
  expect_equal(resolution[1], 0.1)
  expect_equal(resolution[length(resolution)], 1.5)
})


# --- Test 7: Find Neighbors ---  # NOTE: this function will still run even with faulty data and produce results with nothin in the 'neighbour list'
test_that("FindNeighbors works correctly", {
  scrna_test <- scrna
  
  dimensionReduction <- function_pca_dimensions(scrna_test) # Based on your PCA steps
  scrna_test <- FindNeighbors(scrna_test, dims = 1:dimensionReduction)
  
  graph_names <- names(scrna_test@graphs)
  
  # Expect both KNN and SNN graphs to exist - please check if this naming is indeed consistent between sample runs.
  expect_true("RNA_nn" %in% graph_names)
  expect_true("RNA_snn" %in% graph_names)
  
  # Optionally inspect the structure (for debugging)
  message("Neighbor graphs found: ", paste(graph_names, collapse = ", "))
  
  # Save the state for the next test
  scrna <<- scrna_test
})

# --- Test 8: Clustering ---
test_that("Clustering works correctly", {
  scrna_test <- scrna
  
  resolution <- get_resolution_range(opt$integration)
  
  scrna_test <- FindClusters(scrna_test, resolution = resolution)
  
  # Check if clustering was applied correctly
  expect_true("seurat_clusters" %in% colnames(scrna_test@meta.data))
  expect_true(length(unique(scrna_test$seurat_clusters)) > 1)
  
  # Save the state for the next test
  scrna <<- scrna_test
})

# --- Test 9: Clustree plot creation and saving ---
test_that("Clustree plot is generated and saved successfully", {

  expect_true("seurat_clusters" %in% colnames(scrna@meta.data))  # Make sure clustering has been done
  colnames(scrna@meta.data)
  # Generate clustree plot
  clustree_plot <- clustree(scrna)    # why do i need to add this prefix? wy does RNA_snn_res get the res number added to it?
  expect_s3_class(clustree_plot, "gg")  # It should return a ggplot object
  
  # Save the clustree plot
  tmp_file <- tempfile(fileext = ".pdf")
  ggsave(tmp_file, clustree_plot, width = 8, height = 15)
  expect_true(file.exists(tmp_file))  # Ensure the file was created
  # Clean up
  file.remove(tmp_file)
})

# Define file path
output_file <- file.path(test_data_dir, "results/clustree_test.txt")
# Define content to write
lines_to_write <- c(
  "Clustree_test ran"
)
# Write the content to the file
writeLines(lines_to_write, con = output_file)
