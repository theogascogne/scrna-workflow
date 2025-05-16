#!/usr/bin/env Rscript

require(monocle3)
require(Seurat)
require(SeuratWrappers)
require(tidyverse)
require(patchwork)
require(magrittr)
library(optparse)
library(testthat)

# Get the path to the 'cellsnake' package
cellsnake_path <- system("python -c 'import cellsnake; print(cellsnake.__path__[0])'", intern = TRUE)

# Define the base path to the testData folder
test_data_dir <- file.path(cellsnake_path, "scrna/workflow/tests/testData")

# Define options with dynamically constructed paths
option_list <- list(
  make_option(c("--rds"), type = "character", default = file.path(test_data_dir, "processed/defaultTest/output.rds"),
              help = "Processed rds file of a Seurat object", metavar = "character"),
  make_option(c("--output.dir"), type = "character", default = file.path(test_data_dir, "results/percent_mt~10/resolution~0.8/trajectory/"),
              help = "Output directory for trajectory plots", metavar = "character"),
  make_option(c("--pplot"), type = "character", default = file.path(test_data_dir, "results/percent_mt~10/resolution~0.8/trajectory/plot_monocle-partition-plot.pdf"),
              help = "Partition plot output", metavar = "character")
)

opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

# Helper to ensure directory exists
ensure_dir_exists <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(path)) stop("Failed to create directory: ", path)
  }
}

# Ensure output directories exist
ensure_dir_exists(opt$output.dir)
ensure_dir_exists(dirname(opt$pplot))


try(
  {
    source(file.path(cellsnake_path,"scrna/workflow/scripts/scrna_workflow_helper_functions/scrna-monocle3-functions.R"))
  },
  silent = TRUE
)


scrna <<- readRDS(file = opt$rds)

test_that("scrna is converted and clustered properly with initialize_cds_and_cluster", {
  # Convert and cluster
  suppressWarnings({
    cds_c <- initialize_cds_and_cluster(scrna)
  })
  
  # Check that the output is a cell_data_set object
  expect_s4_class(cds_c, "cell_data_set")
  
  # Also check that clustering fallback changed clusters compared to just conversion
  suppressWarnings({
    cds_raw <- as.cell_data_set(scrna)
  })
  
  expect_false(identical(cds_raw@clusters$UMAP$clusters, cds_c@clusters$UMAP$clusters), info = "Clusters did not change")
  expect_false(identical(cds_raw@clusters$UMAP$clusters, cds_c@clusters$UMAP$cluster_result), info = "Clusters did not change")
  expect_false(identical(cds_raw@clusters$UMAP$clusters, cds_c@clusters$UMAP$partitions), info = "Clusters did not change")
  
  # Assign globally if needed later
  cds <<- cds_c
})



test_that("Partition and singler plots are generated and saved", {
  generate_partition_and_singler_plot(cds, opt$pplot)
  expect_true(file.exists(opt$pplot))
  expect_gt(file.info(opt$pplot)$size, 1000)
  file.remove(opt$pplot)
})

test_that("Partitions are extracted correctly and plots are generated", {
  # Generate all partition plots at once
  process_all_partitions(cds, opt$output.dir)
  
  # Convert cds to Seurat for partition extraction
  scrna <- as.Seurat(cds, assay = NULL)
  partitions <- extract_valid_partitions(scrna)
  
  expect_true(length(partitions) > 0)
  
  for (partition_id in partitions) {
    file_path <- file.path(opt$output.dir, paste0("plot_monocle-partition-", partition_id, ".pdf"))
    expect_true(file.exists(file_path))
    file.remove(file_path)
  }
})

# Define file path
output_file <- file.path(test_data_dir, "results/scrna_monocle3_test.txt")
# Define content to write
lines_to_write <- c(
  "monocle3 ran"
)
# Write the content to the file
writeLines(lines_to_write, con = output_file)
