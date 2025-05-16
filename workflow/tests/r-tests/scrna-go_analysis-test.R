#!/usr/bin/env Rscript

library(optparse)
require(tidyverse)
library(testthat)
library(AnnotationDbi)

# Get the path to the 'cellsnake' package
cellsnake_path <- system("python -c 'import cellsnake; print(cellsnake.__path__[0])'", intern = TRUE)

# Define the base path to the testData folder
test_data_dir <- file.path(cellsnake_path, "scrna/workflow/tests/testData")

# Define options
option_list <- list(
  make_option(c("--xlsx"), type = "character", 
              default = file.path(test_data_dir, "results/defaultTest/table_all-markers-seurat_clusters.xlsx"),
              help = "Excel table of markers", metavar = "character"),
  make_option(c("--output.rds"), type = "character", 
              default = file.path(test_data_dir, "analyses/seurat_clusters-go.rds"),
              help = "Output RDS file name", metavar = "character"),
  make_option(c("--output.go"), type = "character", 
              default = file.path(test_data_dir, "results/table_GO-enrichment-seurat_clusters.xlsx"),
              help = "Output GO enrichment file", metavar = "character"),
  make_option(c("--output.gse"), type = "character", 
              default = file.path(test_data_dir, "results/table_GO-geneset_enrichment-seurat_clusters.xlsx"),
              help = "Output GSE file", metavar = "character"),
  make_option(c("--mapping"), type = "character", default = "org.Hs.eg.db",
              help = "Gene mapping database", metavar = "character"),
  make_option(c("--pval"), type = "double", default = 0.05,
                       help = "P value treshold [default= %default]", metavar = "character"),
  make_option(c("--logfc.treshold"), type = "double", default = 1.5,
                        help = "LogFC [default= %default]", metavar = "character"
  )
)

# parse options
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)
library(opt$mapping, character.only = TRUE)  # needed to run this standalone

try(
  {
    source(file.path(cellsnake_path,"scrna/workflow/scripts/scrna_workflow_helper_functions/scrna-go-analysis-functions.R"))
  },
  silent = TRUE
)

test_that("GO enrichment function processes input and produces valid results", {
  All_Features <- openxlsx::read.xlsx(opt$xlsx)
  
  # Use actual defined function
  all_go_results <- All_Features %>%
    split(.$cluster) %>%
    purrr::map(~ function_enrichment_go_singlecell(., p = opt$pval, f = opt$logfc.treshold))
  
  # Basic structure tests
  expect_length(all_go_results, length(unique(All_Features$cluster)))
  
  # Ensure outputs are either NULL or expected S4 classes
  for (result in all_go_results) {
    expect_true(is.null(result[[1]]) || inherits(result[[1]], "enrichResult"))
    expect_true(is.null(result[[2]]) || inherits(result[[2]], "gseGO"))
  }
  
  # Share result across tests
  assign("go_results", all_go_results, envir = .GlobalEnv)
})

test_that("GO results are saved correctly into RDS and Excel files", {
  # Pull results from global
  all_go_results <- get("go_results", envir = .GlobalEnv)
  
  # Create temp paths
  tmp_rds <- tempfile(fileext = ".rds")
  tmp_go  <- tempfile(fileext = ".xlsx")
  tmp_gse <- tempfile(fileext = ".xlsx")
  
  # Use helper functions
  saveRDS(all_go_results, tmp_rds)
  save_go_results(all_go_results, tmp_go)
  save_gse_results(all_go_results, tmp_gse)
  
  # Check files exist
  expect_true(file.exists(tmp_rds))
  expect_true(file.exists(tmp_go))
  expect_true(file.exists(tmp_gse))
  
  # Optional: uncomment for content validation
  # go_data <- openxlsx::read.xlsx(tmp_go)
  # gse_data <- openxlsx::read.xlsx(tmp_gse)
  # expect_true(nrow(go_data) > 0)
  # expect_true(nrow(gse_data) > 0)
  
  # Clean up
  file.remove(tmp_rds, tmp_go, tmp_gse)
})

# Define file path
output_file <- file.path(test_data_dir, "results/scrna-go-analysis-test.txt")
# Define content to write
lines_to_write <- c(
  "kegg ran"
)
# Write the content to the file
writeLines(lines_to_write, con = output_file)

