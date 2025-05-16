#!/usr/bin/env Rscript

require(Seurat)
require(tidyverse)
require(optparse)

# Get the path to the 'cellsnake' package
cellsnake_path <- system("python -c 'import cellsnake; print(cellsnake.__path__[0])'", intern = TRUE)

# Define the base path to the testData folder
test_data_dir <- file.path(cellsnake_path, "scrna/workflow/tests/testData")

# !/usr/bin/env Rscript
option_list <- list(
  optparse::make_option(c("--rds"),
                        type = "character", default = NULL,
                        help = "Processed rds file of a Seurat object", metavar = "character"
  ),
  optparse::make_option(c("--xlsx"),
                        type = "character", default = NULL,
                        help = "Excel table of markers for input", metavar = "character"
  ),
  optparse::make_option(c("--output.plot"),
                        type = "character", default = "output.pdf",
                        help = "Output plot directory", metavar = "character"
  ),
  optparse::make_option(c("--idents"),
                        type = "character", default = "seurat_clusters",
                        help = "Meta data column name for marker analysis", metavar = "character"
  ),
  optparse::make_option(c("--output.average"),
                        type = "character", default = "output.xlsx",
                        help = "Output average expression table", metavar = "character"
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
    source(file.path(cellsnake_path,"scrna/workflow/scripts/scrna_workflow_helper_functions/scrna-marker-heatmap-functions.R"))
  },
  silent = TRUE
)

scrna <- readRDS(file = opt$rds)
DefaultAssay(scrna) <- "RNA"

# 1. Scale data and extract markers
result <- scale_scrna_data(scrna, opt$xlsx, opt$idents)
scrna <- result$scrna
not_all_genes <- result$scaled_genes

# 2. Plot and save heatmap
plot_and_save_heatmap(scrna, not_all_genes, opt$output.plot)

# 3. Write average expression table
write_avg_expression_to_excel(scrna, opt$idents, opt$output.average)