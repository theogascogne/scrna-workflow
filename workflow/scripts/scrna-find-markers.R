#!/usr/bin/env Rscript

require(Seurat)
require(tidyverse)
require(optparse)

# Get the path to the 'cellsnake' package
cellsnake_path <- system("python -c 'import cellsnake; print(cellsnake.__path__[0])'", intern = TRUE)

# Define the base path to the testData folder
test_data_dir <- file.path(cellsnake_path, "scrna/workflow/tests/testData")


option_list <- list(
  optparse::make_option(c("--rds"),
    type = "character", default = NULL,
    help = "Processed rds file of a Seurat object", metavar = "character"
  ),
  optparse::make_option(c("--logfc.threshold"),
    type = "double", default = 0.25,
    help = "LogFC [default= %default]", metavar = "character"
  ),
  optparse::make_option(c("--test.use"),
    type = "character", default = "wilcox",
    help = "Test use [default= %default]", metavar = "character"
  ),
  optparse::make_option(c("--output.rds"),
    type = "character", default = NULL,
    help = "RDS table of all markers", metavar = "character"
  ),
  optparse::make_option(c("--idents"),
    type = "character", default = "seurat_clusters",
    help = "Meta data column name for marker analysis", metavar = "character"
  )
)

# parse options
if (!exists("opt")) {
  opt_parser <- OptionParser(option_list = option_list)
  opt <- parse_args(opt_parser)
}

if (is.null(opt$rds)) {
  optparse::print_help(opt_parser)
  stop("At least one argument must be supplied (rds file)", call. = FALSE)
}

# unused for now
try(
  {
    source(file.path(cellsnake_path,"scrna/workflow/scripts/scrna_workflow_helper_functions/scrna-find-markers-functions.R"))
  },
  silent = TRUE
)

scrna <- readRDS(file = opt$rds)
DefaultAssay(scrna) <- "RNA"

Idents(object = scrna) <- scrna@meta.data[[opt$idents]]

all_markers <- FindAllMarkers(scrna, logfc.threshold = opt$logfc.threshold, test.use = opt$test.use)

saveRDS(all_markers, file = opt$output.rds)