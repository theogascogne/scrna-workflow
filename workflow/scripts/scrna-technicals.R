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
  optparse::make_option(c("--sampleid"),
    type = "character", default = NULL,
    help = "Sample ID", metavar = "character"
  ),
  optparse::make_option(c("--fplot"),
    type = "character", default = NULL,
    help = "nFeature plot", metavar = "character"
  ),
  optparse::make_option(c("--cplot"),
    type = "character", default = NULL,
    help = "nCount plot", metavar = "character"
  ),
  optparse::make_option(c("--mtplot"),
    type = "character", default = NULL,
    help = "Percent MT plot", metavar = "character"
  ),
  optparse::make_option(c("--rpplot"),
    type = "character", default = NULL,
    help = "Ribo plot", metavar = "character"
  )
)

if (!exists("opt")) {
  opt_parser <- OptionParser(option_list = option_list)
  opt <- parse_args(opt_parser)
}

if (is.null(opt$rds)) {
  optparse::print_help(opt_parser)
  stop("At least one argument must be supplied (rds file and sampleid)", call. = FALSE)
}

try(
  {
    source(file.path(cellsnake_path,"scrna/workflow/scripts/scrna_workflow_helper_functions/scrna-technicals-functions.R"))
  },
  silent = TRUE
)


scrna <- readRDS(file = opt$rds)


FeaturePlot(scrna, features = "nFeature_RNA", pt.size = 0.1, raster = FALSE)

ggsave(opt$fplot, width = 7, height = 5)


FeaturePlot(scrna, features = "nCount_RNA", pt.size = 0.1, raster = FALSE)

ggsave(opt$cplot, width = 7, height = 5)


FeaturePlot(scrna, features = "percent.mt", pt.size = 0.1, raster = FALSE)

ggsave(opt$mtplot, width = 7, height = 5)

FeaturePlot(scrna, features = "percent.rp", pt.size = 0.1, raster = FALSE)

ggsave(opt$rpplot, width = 7, height = 5)
