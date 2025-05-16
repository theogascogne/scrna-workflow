#!/usr/bin/env Rscript

require(tidyverse)
require(optparse)
require(Seurat)
require(clustree)

# Get the path to the 'cellsnake' package
cellsnake_path <- system("python -c 'import cellsnake; print(cellsnake.__path__[0])'", intern = TRUE)

# Define the base path to the testData folder
test_data_dir <- file.path(cellsnake_path, "scrna/workflow/tests/testData")


option_list <- list(
  optparse::make_option(c("--scale.factor"),
                        type = "integer", default = 10000,
                        help = "Scale factor [default= %default]", metavar = "character"
  ),
  optparse::make_option(c("--nfeatures"),
                        type = "integer", default = 2000,
                        help = "Highly variable features [default= %default]", metavar = "integer"
  ),
  optparse::make_option(c("--variable.selection.method"),
                        type = "character", default = "vst",
                        help = "Find variable features selection method [default= %default]", metavar = "character"
  ),
  optparse::make_option(c("--rds"),
                        type = "character", default = NULL,
                        help = "RAW rds file of a Seurat object", metavar = "character"
  ),
  optparse::make_option(c("--normalization.method"),
                        type = "character", default = "LogNormalize",
                        help = "Normalization method[default= %default]", metavar = "character"
  ),
  optparse::make_option(c("--integration"), action = "store_true", default = FALSE),
  optparse::make_option(c("--clplot"),
                        type = "character", default = "clustree.pdf",
                        help = "Output clustree file name", metavar = "character"
  ),
  optparse::make_option(c("--jeplot"),
                        type = "character", default = "jackandelbow.pdf",
                        help = "Output jack and elbow file name", metavar = "character"
  ),
  optparse::make_option(c("--hvfplot"),
                        type = "character", default = "variable-features.pdf",
                        help = "Variable features file name", metavar = "character"
  ),
  optparse::make_option(c("--heplot"),
                        type = "character", default = "dimheatmap.pdf",
                        help = "Dim heatmap plot file name", metavar = "character"
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

try({
  source("~/miniconda3/envs/test/lib/python3.9/site-packages/cellsnake/scrna/workflow/scripts/scrna-functions.R")
}, silent = TRUE)

try({
  source(paste0(system("python -c 'import os; import cellsnake; print(os.path.dirname(cellsnake.__file__))'", intern = TRUE), "/scrna/workflow/scripts/scrna-functions.R"))
}, silent = TRUE)

try(
  {
    source(file.path(cellsnake_path,"scrna/workflow/scripts/scrna_workflow_helper_functions/scrna-clusteringtree-functions.R"))
  },
  silent = TRUE
)

# Read the Seurat object
scrna <- readRDS(file = opt$rds)

# Handle normalization or integration based on the integration flag
scrna <- handle_normalization_or_integration(scrna, opt)

# Perform scaling and PCA
pca_results <- perform_scaling_and_pca(scrna)
scrna <- pca_results$scrna
not.all.genes <- pca_results$features

# If not integrating, generate variable features and plot the results
if (isFALSE(opt$integration)) {
  top10 <- head(not.all.genes, 10)
  plot_variable_features(scrna, top10, opt$hvfplot, opt$heplot)
  scrna <- plot_jackstraw_and_elbow(scrna, opt$jeplot)
}

# Set the resolution range based on integration status
resolution <- get_resolution_range(opt$integration)

# Find neighbors and clusters
dimensionReduction <- function_pca_dimensions(scrna)
scrna <- FindNeighbors(scrna, dims = 1:dimensionReduction)
scrna <- FindClusters(scrna, resolution = resolution)

# Generate the clustree plot and save it
clustree(scrna) -> p1
ggsave(opt$clplot, p1, width = 8, height = 15)
