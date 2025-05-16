#!/usr/bin/env Rscript

# Load required libraries
require(optparse)
require(SingleR)
require(pheatmap)
require(Seurat)
require(tidyverse)

# Get the path to the 'cellsnake' package
cellsnake_path <- system("python -c 'import cellsnake; print(cellsnake.__path__[0])'", intern = TRUE)

# Define the base path to the testData folder
test_data_dir <- file.path(cellsnake_path, "scrna/workflow/tests/testData")

# Define command line options
option_list <- list(
  optparse::make_option(c("--rds"),
                        type = "character", default = NULL,
                        help = "A list of RDS files of Seurat objects", metavar = "character"
  ),
  optparse::make_option(c("--sheplot"),
                        type = "character", default = "sheplot.pdf",
                        help = "Output score heatmap plot file name", metavar = "character"
  ),
  optparse::make_option(c("--sheplottop"),
                        type = "character", default = "sheplot.pdf",
                        help = "Output score heatmap plot file name, top 20", metavar = "character"
  ),
  optparse::make_option(c("--pheplot"),
                        type = "character", default = "pheplot.pdf",
                        help = "Output heatmap plot file name", metavar = "character"
  ),
  optparse::make_option(c("--idents"),
                        type = "character", default = "seurat_clusters",
                        help = "Meta data column name", metavar = "character"
  ),
  optparse::make_option(c("--csv"),
                        type = "character", default = NULL,
                        help = "A meta data table", metavar = "character"
  ),
  optparse::make_option(c("--prediction"),
                        type = "character", default = "pred.rds",
                        help = "Input prediction file", metavar = "character"
  ),
  optparse::make_option(c("--xlsx"),
                        type = "character", default = "predictions.xlsx",
                        help = "Input prediction file", metavar = "character"
  )
)

# Parse the command line options
if (!exists("opt")) {
  opt_parser <- OptionParser(option_list = option_list)
  opt <- parse_args(opt_parser)
}

# Ensure the input RDS file is provided
if (is.null(opt$rds)) {
  optparse::print_help(opt_parser)
  stop("At least one argument must be supplied (rds and sampleid)", call. = FALSE)
}

try(
  {
    source(file.path(cellsnake_path,"scrna/workflow/scripts/scrna_workflow_helper_functions/scrna-singler-plots-functions.R"))
  },
  silent = TRUE
)


# Load Seurat object and predictions
scrna <- readRDS(file = opt$rds)
DefaultAssay(scrna) <- "RNA"
pred <- readRDS(opt$prediction)

# Generate and save the score heatmap (regular)
save_score_heatmap(pred, opt$sheplot, width = 15, height = 8)

# Generate and save the top 20 score heatmap
save_score_heatmap(pred, opt$sheplottop, show_labels = FALSE, max_labels = 20, width = 7, height = 4)

# Create the table for heatmap
tab <- table(Assigned = pred$pruned.labels, Cluster = scrna@meta.data[[opt$idents]])

# Generate and save the heatmap plot
save_heatmap_plot(scrna, pred, tab, opt$pheplot, width_factor = 0.10, height_factor = 0.10)

# Save the table to an Excel file
save_table_to_excel(tab, opt$xlsx)
