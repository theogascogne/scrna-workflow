#!/usr/bin/env Rscript

require(optparse)
require(monocle3)
require(Seurat)
require(SeuratWrappers)
require(tidyverse)
require(patchwork)
require(magrittr)

# Get the path to the 'cellsnake' package
cellsnake_path <- system("python -c 'import cellsnake; print(cellsnake.__path__[0])'", intern = TRUE)

# Define the base path to the testData folder
test_data_dir <- file.path(cellsnake_path, "scrna/workflow/tests/testData")

# CLI options
option_list <- list(
  make_option(c("--rds"), type = "character", default = NULL,
              help = "Processed rds file of a Seurat object", metavar = "character"),
  make_option(c("--pplot"), type = "character", default = NULL,
              help = "Partition plot", metavar = "character"),
  make_option(c("--output.dir"), type = "character", default = NULL,
              help = "Output plot directory", metavar = "character")
)

if (!exists("opt")) {
  opt_parser <- OptionParser(option_list = option_list)
  opt <- parse_args(opt_parser)
}

# Exit if input is missing
if (is.null(opt$rds)) {
  optparse::print_help(opt_parser)
  stop("At least one argument must be supplied (rds file)", call. = FALSE)
}

try(
  {
    source(file.path(cellsnake_path,"scrna/workflow/scripts/scrna_workflow_helper_functions/scrna-monocle3-functions.R"))
  },
  silent = TRUE
)

# Load Seurat object
scrna <- readRDS(file = opt$rds)

# Initialize and cluster the CDS object
cds <- initialize_cds_and_cluster(scrna)

# Generate and save combined singler & partition plot
generate_partition_and_singler_plot(cds, opt$pplot)

# Process and plot partitions, saving output files
process_all_partitions(cds, opt$output.dir)