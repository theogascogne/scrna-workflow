#!/usr/bin/env Rscript

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

opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

if (is.null(opt$rds)) {
  optparse::print_help(opt_parser)
  stop("At least one argument must be supplied (rds file and sampleid)", call. = FALSE)
}

# Install required packages if not already installed
if (!requireNamespace("SeuratDisk", quietly = TRUE)) {
  remotes::install_github("mojaveazure/seurat-disk", upgrade = "never")
}

if (!requireNamespace("scCustomize", quietly = TRUE)) {
  remotes::install_github("slawlor/scCustomize")
}

# Load libraries
require(Seurat)
require(SeuratDisk)
require(tidyverse)
require(scCustomize)

# Read Seurat object
scrna <- readRDS(file = opt$rds)

# Convert Seurat v5 to Seurat v3/4 assay structure using scCustomize
scrna <- scCustomize::Convert_Assay(seurat_object = scrna, convert_to = "V3")

# Update Seurat object (in case needed)
scrna <- UpdateSeuratObject(scrna)

# Default Assay set to "RNA"
DefaultAssay(scrna) <- "RNA"

# Diet Seurat (reduce memory)
scrna <- DietSeurat(scrna)

# Convert metadata factors to characters
scrna@meta.data <- scrna@meta.data %>%
  dplyr::mutate(dplyr::across(where(is.factor), as.character))

# Save the Seurat object as h5Seurat format
output_file_name <- str_remove_all(opt$output, ".h5ad$")
SaveH5Seurat(scrna, filename = paste0(output_file_name, ".h5Seurat"), overwrite = TRUE)

# Convert Seurat object to h5ad format
SeuratDisk::Convert(paste0(output_file_name, ".h5Seurat"), dest = "h5ad", overwrite = TRUE)