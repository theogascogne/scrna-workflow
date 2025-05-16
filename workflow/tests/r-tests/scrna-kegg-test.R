#!/usr/bin/env Rscript

library(optparse)
require(tidyverse)
library(testthat)
library(AnnotationDbi)
# Get the path to the 'cellsnake' package
cellsnake_path <- system("python -c 'import cellsnake; print(cellsnake.__path__[0])'", intern = TRUE)

# Define the base path to the testData folder
test_data_dir <- file.path(cellsnake_path, "scrna/workflow/tests/testData")

option_list <- list(
  make_option(c("--xlsx"), type = "character", 
              default = file.path(test_data_dir, "results/defaultTest/table_all-markers-seurat_clusters.xlsx"),
              help = "Excel table of markers", metavar = "character"),
  make_option(c("--output.rds"), type = "character", 
              default = file.path(test_data_dir, "analyses/seurat_clusters-kegg.rds"),
              help = "Output RDS file name", metavar = "character"),
  make_option(c("--output.kegg"), type = "character", 
              default = file.path(test_data_dir, "results/10X_17_029/percent_mt~10/resolution~0.8/enrichment_analysis/table_KEGG-enrichment-seurat_clusters.xlsx"),
              help = "Output KEGG enrichment file", metavar = "character"),
  make_option(c("--output.mkegg"), type = "character", 
              default = file.path(test_data_dir, "results/10X_17_029/percent_mt~10/resolution~0.8/enrichment_analysis/table_KEGG-module_enrichment-seurat_clusters.xlsx"),
              help = "Output KEGG module enrichment file", metavar = "character"),
  make_option(c("--output.gse"), type = "character", 
              default = file.path(test_data_dir, "results/10X_17_029/percent_mt~10/resolution~0.8/enrichment_analysis/table_KEGG-geneset_enrichment-seurat_clusters.xlsx"),
              help = "Output KEGG gene set enrichment file", metavar = "character"),
  make_option(c("--output.mgse"), type = "character", 
              default = file.path(test_data_dir, "results/10X_17_029/percent_mt~10/resolution~0.8/enrichment_analysis/table_KEGG-module_geneset_enrichment-seurat_clusters.xlsx"),
              help = "Output KEGG module gene set enrichment file", metavar = "character"),
  make_option(c("--mapping"), type = "character", default = "org.Hs.eg.db",
              help = "Gene mapping database", metavar = "character"),
  make_option(c("--organism"), type = "character", default = "hsa",
              help = "Organism code (KEGG-compatible)", metavar = "character"),
  make_option(c("--pval"),
              type = "double", default = 0.05,
              help = "P value threshold [default= %default]", metavar = "double"),
  make_option(c("--logfc.treshold"),
              type = "double", default = 1,
              help = "LogFC threshold [default= %default]", metavar = "double")
)

# Parse options
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)
library(opt$mapping, character.only = TRUE)  # needed to run this standalone

try(
  {
    source(file.path(cellsnake_path,"scrna/workflow/scripts/scrna_workflow_helper_functions/scrna-kegg-functions.R"))
  },
  silent = TRUE
)
test_that("KEGG enrichment function processes input and produces valid results", {
  
  All_Features <- openxlsx::read.xlsx(opt$xlsx)
  
  suppressWarnings({
    all_kegg_results <- All_Features %>%
      split(.$cluster) %>%
      purrr::map(~ function_enrichment_kegg_singlecell(., p = opt$pval, f = opt$logfc.treshold, mapping_db = get(opt$mapping), organism_code = opt$organism))
  })
  
  # Check actual row counts for each enrichment result
  detailed_results <- purrr::imap(all_kegg_results, function(res, cluster) {
    get_n <- function(x) {
      if (is.null(x)) return(0)
      tryCatch(nrow(x@result), error = function(e) 0)
    }
    
    tibble(
      cluster = cluster,
      enrichKEGG = get_n(res[[1]]),
      gseKEGG = get_n(res[[2]]),
      enrichMKEGG = get_n(res[[3]]),
      gseMKEGG = get_n(res[[4]])
    )
  }) %>% bind_rows()
  
  print(detailed_results)
  
  expect_length(all_kegg_results, length(unique(All_Features$cluster)))
  
  for (result in all_kegg_results) {
    expect_true(is.null(result[[1]]) || inherits(result[[1]], "enrichResult"))
    expect_true(is.null(result[[2]]) || inherits(result[[2]], "gseaResult"))
    expect_true(is.null(result[[3]]) || inherits(result[[3]], "enrichResult"))
    expect_true(is.null(result[[4]]) || inherits(result[[4]], "gseaResult"))
  }
  
  tmp_rds <- tempfile(fileext = ".rds")
  tmp_kegg <- tempfile(fileext = ".xlsx")
  tmp_gse <- tempfile(fileext = ".xlsx")
  tmp_mkegg <- tempfile(fileext = ".xlsx")
  tmp_mgse <- tempfile(fileext = ".xlsx")
  
  saveRDS(all_kegg_results, tmp_rds)
  save_kegg_results(all_kegg_results, 1, tmp_kegg)
  save_kegg_results(all_kegg_results, 2, tmp_gse)
  save_kegg_results(all_kegg_results, 3, tmp_mkegg)
  save_kegg_results(all_kegg_results, 4, tmp_mgse)
  
  expect_true(file.exists(tmp_rds))
  expect_true(file.exists(tmp_kegg))
  expect_true(file.exists(tmp_gse))
  expect_true(file.exists(tmp_mkegg))
  expect_true(file.exists(tmp_mgse))
  
  # ... (rest of your failure checks and cleanup)
  
  file.remove(tmp_rds)
  file.remove(tmp_kegg)
  file.remove(tmp_gse)
  file.remove(tmp_mkegg)
  file.remove(tmp_mgse)
})

# Define file path
output_file <- file.path(test_data_dir, "results/scrna-kegg-test.txt")
# Define content to write
lines_to_write <- c(
  "kegg ran"
)
# Write the content to the file
writeLines(lines_to_write, con = output_file)

