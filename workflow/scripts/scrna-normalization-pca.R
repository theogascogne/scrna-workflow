#!/usr/bin/env Rscript

require(optparse)
require(SingleR)
require(celldex)
require(tidyverse)
require(Seurat)
require(patchwork)

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
                        help = "A RAW rds file of a Seurat object", metavar = "character"
  ),
  optparse::make_option(c("--normalization.method"),
                        type = "character", default = "LogNormalize",
                        help = "Normalization method[default= %default]", metavar = "character"
  ),
  optparse::make_option(c("--doublet.filter"), action = "store_true", default = FALSE),
  optparse::make_option(c("--integration"), action = "store_true", default = FALSE),
  optparse::make_option(c("--umap"), action = "store_true", default = FALSE),
  optparse::make_option(c("--tsne"), action = "store_true", default = FALSE),
  optparse::make_option(c("--resolution"),
                        type = "character", default = "0.8",
                        help = "Resolution [default= %default]", metavar = "character"
  ),
  optparse::make_option(c("--output.rds"),
                        type = "character", default = "output.rds",
                        help = "Output RDS file name [default= %default]", metavar = "character"
  ),
  optparse::make_option(c("--cpu"),
                        type = "integer", default = 5,
                        help = "Number of CPU for parallel run [default= %default]", metavar = "character"
  ),
  optparse::make_option(c("--reference"),
                        type = "character", default = "HumanPrimaryCellAtlasData",
                        help = "SingleR reference for annotation", metavar = "character"
  )
)

if (!exists("opt")) {
  opt_parser <- OptionParser(option_list = option_list)
  opt <- parse_args(opt_parser)
}

if (is.null(opt$rds)) {
  optparse::print_help(opt_parser)
  stop("At least one argument must be supplied (rds file)", call. = FALSE)
}

try(
  {
    source("workflow/scripts/scrna-functions.R")
  },
  silent = TRUE
)
try(
  {
    source(paste0(system("python -c 'import os; import cellsnake; print(os.path.dirname(cellsnake.__file__))'", intern = TRUE), "/scrna/workflow/scripts/scrna-functions.R"))
  },
  silent = TRUE
)


try(
  {
    source(file.path(cellsnake_path,"scrna/workflow/scripts/scrna_workflow_helper_functions/scrna-normalization-pca-functions.R"))
  },
  silent = TRUE
)

# Read the Seurat object
scrna <- readRDS(file = opt$rds)

# Step 1: Normalize and select features
scrna <- normalize_and_select_features(scrna, opt)

# Step 2: Run PCA
pca_result <- run_pca_pipeline(scrna)
scrna <- pca_result$scrna
dimensionReduction <- pca_result$dims

# Step 3: Perform clustering
scrna <- retrieve_clustering(scrna, opt, dimensionReduction)

# Step 4: Run TSNE and UMAP if required
scrna <- run_TSNE_UMAP(scrna, opt, 1:dimensionReduction)

# Step 5: Perform doublet filtering if required
scrna <- run_doublet_filter(scrna, opt)

# Step 6: Annotate with SingleR
scrna <- annotate_with_singleR(scrna, opt$reference)

# Output RDS file
saveRDS(scrna, file = opt$output.rds)
