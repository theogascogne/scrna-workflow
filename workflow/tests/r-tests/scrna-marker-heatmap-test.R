#!/usr/bin/env Rscript

require(testthat)
require(Seurat)
require(optparse)
require(tidyverse)

# Get the path to the 'cellsnake' package
cellsnake_path <- system("python -c 'import cellsnake; print(cellsnake.__path__[0])'", intern = TRUE)

# Define the base path to the testData folder
test_data_dir <- file.path(cellsnake_path, "scrna/workflow/tests/testData")

# Define options
option_list <- list(
  make_option(c("--rds"), type = "character", default = file.path(test_data_dir, "processed/defaultTest/output.rds"),
              help = "Seurat object RDS file", metavar = "character"),
  make_option(c("--xlsx"), type = "character", default = file.path(test_data_dir, "results/defaultTest/table_positive-markers-seurat_clusters.xlsx"),
              help = "Excel table of markers", metavar = "character"),
  make_option(c("--output.plot"), type = "character", default = file.path(test_data_dir, "results/defaultTest/plot_marker-heatmap-seurat_clusters.pdf"),
              help = "Output heatmap PDF", metavar = "character"),
  make_option(c("--output.average"), type = "character", default = file.path(test_data_dir, "results/defaultTest/table_average-expression-seurat_clusters.xlsx"),
              help = "Output average expression table", metavar = "character"),
  make_option(c("--idents"), type = "character", default = "seurat_clusters",
              help = "Metadata column for clustering", metavar = "character")
)

opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

try(
  {
    source(file.path(cellsnake_path,"scrna/workflow/scripts/scrna_workflow_helper_functions/scrna-marker-heatmap-functions.R"))
  },
  silent = TRUE
)

# Initialize global variables at the start
scrna <<- readRDS(file = opt$rds)
DefaultAssay(scrna) <- "RNA"
p1 <<- NULL
avg_expr <<- NULL

test_that("Data is scaled correctly", {
  result <- scale_scrna_data(scrna, opt$xlsx, opt$idents)
  scrna <<- result$scrna
  not.all.genes <<- result$scaled_genes
  
  scaled_data <- scrna[["RNA"]]$scale.data
  
  expect_true(!is.null(scaled_data))
  expect_true(is.matrix(scaled_data))
  expect_true(is.numeric(scaled_data[1, 1]))
  
  if ("SST" %in% rownames(scaled_data)) {
    expect_lt(abs(mean(scaled_data["SST", ])), 1e-6)
    expect_lt(abs(sd(scaled_data["SST", ]) - 1), 1e-6)
  }
})

test_that("Heatmap is plotted and saved correctly", {
  expect_true(!is.null(scrna))
  
  tmp_file <- tempfile(fileext = ".pdf")
  plot_and_save_heatmap(scrna, not.all.genes, tmp_file)
  
  expect_true(file.exists(tmp_file))
  expect_gt(file.info(tmp_file)$size, 1000)
  
  file.remove(tmp_file)
})

# Test 3: Writing Average Expression to Excel
test_that("Average expression is written to Excel correctly", {
  # Reuse global variables `scrna` and `p1`
  expect_true(!is.null(scrna))  # Ensure that scrna is available
  
  # Write to Excel
  excel_file <- tempfile(fileext = ".xlsx")
  write_avg_expression_to_excel(scrna, opt$idents, excel_file)
  
  # Check if the Excel file exists and is not empty
  expect_true(file.exists(excel_file))
  expect_gt(file.info(excel_file)$size, 1000)
  
  # Clean up temporary Excel file
  file.remove(excel_file)
})

# Define file path
output_file <- file.path(test_data_dir, "results/scrna-marker-heatmap-test.txt")
# Define content to write
lines_to_write <- c(
  "top marker ran"
)
# Write the content to the file
writeLines(lines_to_write, con = output_file)
