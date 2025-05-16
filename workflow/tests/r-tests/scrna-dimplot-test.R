#!/usr/bin/env Rscript

library(optparse)
require(patchwork)
require(plotly)
require(Seurat)
require(tidyverse)
library(testthat)

# Get the path to the 'cellsnake' package
cellsnake_path <- system("python -c 'import cellsnake; print(cellsnake.__path__[0])'", intern = TRUE)

# Define the base path to the testData folder
test_data_dir <- file.path(cellsnake_path, "scrna/workflow/tests/testData")

# Define options for testing, using the defaults from the shell command
option_list <- list(
  make_option(c("--rds"), type = "character", default = file.path(test_data_dir, "processed/defaultTest/output.rds"),
              help = "Path to Seurat object RDS file", metavar = "character"),
  make_option(c("--reduction.type"), type = "character", default = "tsne", help = "Reduction type (umap/tsne)", metavar = "character"),
  make_option(c("--pdfplot"), type = "character", default = file.path(test_data_dir, "results/testSample/percent_mt~10/resolution~0.8/plot_dimplot_tsne-seurat_clusters.pdf"),
              help = "Path to save PDF plot", metavar = "character"),
  make_option(c("--htmlplot"), type = "character", default = file.path(test_data_dir, "results/testSample/percent_mt~10/resolution~0.8/plot_dimplot_tsne-seurat_clusters.html"),
              help = "Path to save HTML plot", metavar = "character"),
  make_option(c("--csv"), type = "character", default = NULL, help = "CSV metadata file (for SingleR annotation)", metavar = "character"),
  make_option(c("--idents"), type = "character", default = "seurat_clusters", help = "Metadata column for identities", metavar = "character"),
  make_option(c("--percentage"), type = "double", default = 5, help = "Minimum cluster percentage", metavar = "double"),
  make_option(c("--labels"), action = "store_true", default = TRUE, help = "Add labels to plot")  # Assuming you want labels, as they are in the shell command
)

# Parse options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)


tryCatch(
  {
    source("~/miniconda3/envs/test/lib/python3.9/site-packages/cellsnake/scrna/workflow/scripts/scrna-functions.R")
  },
  error = function(cond) {
    source(paste0(system("python -c 'import os; import cellsnake; print(os.path.dirname(cellsnake.__file__))'", intern = TRUE), "/scrna/workflow/scripts/scrna-functions.R"))
  }
)

try(
  {
    source(file.path(cellsnake_path,"scrna/workflow/scripts/scrna_workflow_helper_functions/scrna-dimplot-functions.R"))
  },
  silent = TRUE
)

# Load initial scrna object
scrna <- readRDS(opt$rds)  # Replace with actual file path

test_that("Metadata CSV and pruned.labels RDS are correctly processed", {
  # Check if opt$csv is not null, if it is, skip the test
  if (is.null(opt$csv)) {
    skip("Test is skipped because opt$csv is null")
  }
  
  # Start with a fresh scrna object
  scrna_csv_test <- scrna
  scrna_rds_test <- scrna
  
  ## ---------- Test CSV Metadata Merge ----------
  # Create a temporary CSV metadata file
  csv_tmp <- tempfile(fileext = ".csv")
  csv_metadata <- data.frame(
    barcodes = rownames(scrna_csv_test),
    new_col = sample(1:10, nrow(scrna_csv_test), replace = TRUE)
  )
  write.csv(csv_metadata, file = csv_tmp, row.names = FALSE)
  
  # Simulate opt$csv being a CSV
  metadata <- read.csv(csv_tmp, row.names = 1)
  
  scrna_csv_test <- merge_metadata_csv(scrna_csv_test, csv_tmp)
  # Expect new_col to be in the metadata
  expect_true("new_col" %in% colnames(scrna_csv_test@meta.data))
  
  ## ---------- Test RDS Pruned Labels Merge ----------
  # Create a temporary RDS file with pruned.labels
  rds_tmp <- tempfile(fileext = ".rds")
  pred_metadata <- data.frame(pruned.labels = sample(c("LabelA", "LabelB"), ncol(scrna_rds_test), replace = TRUE))
  saveRDS(pred_metadata, rds_tmp)
  
  # Simulate opt$csv being an RDS instead
  scrna_rds_test <- merge_metadata_rds(scrna_rds_test, rds_tmp)
  # Expect singler to be in the metadata
  expect_true("singler" %in% colnames(scrna_rds_test@meta.data))
  
  ## ---------- Clean Up ----------
  unlink(csv_tmp)
  unlink(rds_tmp)
  
  ## ---------- Save Test Results ----------
  scrna <<- scrna_rds_test  # Save one of them if needed for continuity
})


test_that("Identities and palette are set correctly", {
  scrna_test <- scrna
  
  # Check if opt$idents is valid and properly set
  expect_true(opt$idents %in% colnames(scrna_test@meta.data), info = "opt$idents must be a valid column name")
  
  palette <- assign_idents_and_palette(scrna_test, opt$idents)
  # Get number of unique identities
  n <- length(unique(Idents(scrna_test)))  # Number of unique identities
  
  # Check that the correct palette is being used based on n
  if (n <= 10) {
    expect_equal(length(palette), 10, info = "Palette should contain 10 colors when n <= 10")
  } else if (n <= 30) {
    expect_equal(length(palette), 30, info = "Palette should contain 30 colors when 10 < n <= 30")
  } else if (n <= 50) {
    expect_equal(length(palette), 50, info = "Palette should contain 50 colors when 30 < n <= 50")
  } else if (n <= 75) {
    expect_equal(length(palette), 75, info = "Palette should contain 75 colors when 50 < n <= 75")
  } else {
    expect_equal(length(palette), length(palette200), info = "Palette should contain more than 75 colors")
  }
  
  # Ensure that the palette is named correctly
  expect_named(palette)
  
  # Check that identities are set as a factor
  expect_s3_class(Idents(scrna_test), "factor")
  
  # Save the test results
  scrna <<- scrna_test
})

test_that("Cluster breaks are calculated correctly", {
  scrna_test <- scrna
  
  # Get valid cluster IDs
  breaks <- get_valid_cluster_ids(scrna_test, opt$idents, opt$percentage)
  
  expect_type(breaks, "character")
  expect_true(all(breaks %in% Idents(scrna_test)))
})

test_that("Interactive HTML plot is created and saved", {
  scrna_test <- scrna
  
  Idents(object = scrna_test) <- scrna_test@meta.data[[opt$idents]]
  n <- length(Idents(scrna_test) %>% unique())
  palette <- function_color_palette(n)
  palette <- setNames(palette, Idents(scrna_test) %>% unique())
  
  # Generate and save the HTML plot using the helper function
  tmp_html <- tempfile(fileext = ".html")
  plot_created <- generate_dimplot_html(scrna_test, opt$reduction.type, palette, tmp_html)
  
  expect_true(file.exists(tmp_html))
  expect_gt(file.info(tmp_html)$size, 1000)
  
  # Clean up
  file.remove(tmp_html)
})

test_that("PDF plot is created and saved", {
  scrna_test <- scrna
  
  # Set identities and generate the palette
  Idents(object = scrna_test) <- scrna_test@meta.data[[opt$idents]]
  n <- length(Idents(scrna_test) %>% unique())
  palette <- function_color_palette(n)
  palette <- setNames(palette, Idents(scrna_test) %>% unique())
  
  # Get valid clusters for the test
  valid_clusters <- get_valid_cluster_ids(scrna_test, opt$idents, opt$percentage)
  
  # Generate and save the PDF plot using the helper function
  tmp_pdf <- tempfile(fileext = ".pdf")
  generate_dimplot_pdf(scrna_test, opt$reduction.type, palette, valid_clusters, opt$labels, tmp_pdf)
  
  # Check that the file exists and is of a reasonable size
  expect_true(file.exists(tmp_pdf))
  expect_gt(file.info(tmp_pdf)$size, 1000)  # Ensure the file is not empty
  
  # Clean up
  file.remove(tmp_pdf)
})


# Define file path
output_file <- file.path(test_data_dir, "results/scrna_dimplot_test.txt")
# Define content to write
lines_to_write <- c(
  "dimplots ran"
)
# Write the content to the file
writeLines(lines_to_write, con = output_file)