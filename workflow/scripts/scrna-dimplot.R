#!/usr/bin/env Rscript

require(patchwork)
require(plotly)
require(Seurat)
require(tidyverse)
require(optparse)

# Get the path to the 'cellsnake' package
cellsnake_path <- system("python -c 'import cellsnake; print(cellsnake.__path__[0])'", intern = TRUE)

# Define the base path to the testData folder
test_data_dir <- file.path(cellsnake_path, "scrna/workflow/tests/testData")

# Define options for command-line arguments
option_list <- list(
  optparse::make_option(c("--rds"),
                        type = "character", default = "~/Documents/cellsnake_shared/fetal-brain/analyses/processed/percent_mt~10/resolution~0.8/10X_17_028.rds",
                        help = "Processed rds file of a Seurat object", metavar = "character"
  ),
  optparse::make_option(c("--reduction.type"),
                        type = "character", default = "umap",
                        help = "Reduction type, umap or tsne", metavar = "character"
  ),
  optparse::make_option(c("--pdfplot"),
                        type = "character", default = "reduction.pdf",
                        help = "Plot file name", metavar = "character"
  ),
  optparse::make_option(c("--htmlplot"),
                        type = "character", default = "reduction.html",
                        help = "Plot file name", metavar = "character"
  ),
  optparse::make_option(c("--csv"),
                        type = "character", default = NULL,
                        help = "CSV meta file or this can be a singler RDS file", metavar = "character"
  ),
  optparse::make_option(c("--idents"),
                        type = "character", default = "seurat_clusters",
                        help = "Meta data column name for marker analysis", metavar = "character"
  ),
  optparse::make_option(c("--percentage"),
                        type = "double", default = 5,
                        help = "Cluster minimum percentage to plot", metavar = "double"
  ),
  optparse::make_option(c("--labels"), action = "store_true", default = FALSE, help = "Print labels on the plot")
)

if (!exists("opt")) {
  opt_parser <- OptionParser(option_list = option_list)
  opt <- parse_args(opt_parser)
}

if (is.null(opt$rds)) {
  optparse::print_help(opt_parser)
  stop("At least one argument must be supplied (rds file)", call. = FALSE)
}

tryCatch(
  {
    source("~/miniconda3/envs/test/lib/python3.9/site-packages/cellsnake/scrna/workflow/scripts/scrna-functions.R")
  },
  error = function(cond) {
    source(paste0(system("python -c 'import os; import cellsnake; print(os.path.dirname(cellsnake.__file__))'", intern = TRUE), "/scrna/workflow/scripts/scrna-functions.R"))
  }
)

try(
  {
    source(file.path(cellsnake_path,"scrna/workflow/scripts/scrna_workflow_helper_functions/scrna-dimplot-functions.R"))
  },
  silent = TRUE
)

# Load Seurat object
scrna <- readRDS(file = opt$rds)
DefaultAssay(scrna) <- "RNA"

# Apply CSV or RDS metadata if given
if (!is.null(opt$csv)) {
  if (grepl("\\.csv$", opt$csv)) {
    scrna <- merge_metadata_csv(scrna, opt$csv)
  } else if (grepl("\\.rds$", opt$csv)) {
    scrna <- merge_metadata_rds(scrna, opt$csv)
  }
}

# Assign identities
Idents(scrna) <- scrna@meta.data[[opt$idents]]

# Color palette and valid clusters
palette <- assign_idents_and_palette(scrna, opt$idents)
valid_clusters <- get_valid_cluster_ids(scrna, opt$idents, opt$percentage)

# Restrict meta levels to valid clusters
scrna@meta.data[[opt$idents]] <- factor(scrna@meta.data[[opt$idents]], levels = valid_clusters)

# Generate and save HTML plot
generate_dimplot_html(scrna, opt$reduction.type, palette, opt$htmlplot)

# Generate and save PDF plot using the helper function
generate_dimplot_pdf(scrna, opt$reduction.type, palette, valid_clusters, opt$labels, opt$pdfplot)
