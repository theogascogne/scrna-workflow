#!/usr/bin/env Rscript

require(optparse)
require(Seurat)
require(SeuratDisk)
require(tidyverse)
require(scCustomize)

# Get the path to the 'cellsnake' package
cellsnake_path <- system("python -c 'import cellsnake; print(cellsnake.__path__[0])'", intern = TRUE)

# Define the base path to the testData folder
test_data_dir <- file.path(cellsnake_path, "scrna/workflow/tests/testData")


# Define options
option_list <- list(
  optparse::make_option(c("--rds"),
                        type = "character", default = NULL,
                        help = "A list of RDS files of Seurat objects", metavar = "character"
  ),
  optparse::make_option(c("--output"),
                        type = "character", default = "output.h5ad",
                        help = "Output h5ad file name", metavar = "character"
  )
)

# Parse options
if (!exists("opt")) {
  opt_parser <- OptionParser(option_list = option_list)
  opt <- parse_args(opt_parser)
}

# Check for input file
if (is.null(opt$rds)) {
  optparse::print_help(opt_parser)
  stop("At least one argument must be supplied (rds file)", call. = FALSE)
}

# Install required packages if not already installed
if (!requireNamespace("SeuratDisk", quietly = TRUE)) {
  remotes::install_github("mojaveazure/seurat-disk", upgrade = "never")
}

if (!requireNamespace("scCustomize", quietly = TRUE)) {
  remotes::install_github("slawlor/scCustomize")
}

try(
  {
    source(file.path(cellsnake_path,"scrna/workflow/scripts/scrna_workflow_helper_functions/scrna-convert-to-h5ad-functions.R"))
  },
  silent = TRUE
)

# Read Seurat object
scrna <- readRDS(file = opt$rds)

# Convert Seurat v5 to Seurat v3/4 assay structure using scCustomize
scrna <- convert_assay(scrna, version = "V3")

# Update Seurat object
scrna <- UpdateSeuratObject(scrna)

# Set default assay to "RNA"
DefaultAssay(scrna) <- "RNA"

# Diet Seurat (reduce memory usage)
scrna <- DietSeurat(scrna)

# Convert metadata factors to characters
scrna <- convert_factors_to_characters(scrna)

# Save the Seurat object as h5Seurat format
output_file_name <- str_remove_all(opt$output, ".h5ad$")
tmp_h5seurat <- tempfile(fileext = ".h5Seurat")
save_h5seurat(scrna, tmp_h5seurat)

# Convert Seurat object to h5ad format
convert_to_h5ad(tmp_h5seurat, opt$output)

# Clean up temporary files
file.remove(tmp_h5seurat)