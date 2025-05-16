#!/usr/bin/env Rscript

require(optparse)
require(tidyverse)

# Get the path to the 'cellsnake' package
cellsnake_path <- system("python -c 'import cellsnake; print(cellsnake.__path__[0])'", intern = TRUE)

# Define the base path to the testData folder
test_data_dir <- file.path(cellsnake_path, "scrna/workflow/tests/testData")

option_list <- list(
  optparse::make_option(c("--xlsx"),
                        type = "character", default = NULL,
                        help = "Excel table of markers for input", metavar = "character"
  ),
  optparse::make_option(c("--output.plot"),
                        type = "character", default = "output.pdf",
                        help = "Output plot file", metavar = "character"
  )
)

if (!exists("opt")) {
  opt_parser <- OptionParser(option_list = option_list)
  opt <- parse_args(opt_parser)
}

if (is.null(opt$xlsx)) {
  optparse::print_help(opt_parser)
  stop("At least one argument must be supplied (xlsx)", call. = FALSE)
}

try(
  {
    source(file.path(cellsnake_path,"scrna/workflow/scripts/scrna_workflow_helper_functions/scrna-top-marker-plot-functions.R"))
  },
  silent = TRUE
)

Positive_Features <- openxlsx::read.xlsx(opt$xlsx)

df <- extract_top_markers(Positive_Features)
clusters <- get_clusters_from_features(Positive_Features)

generate_pdf_plot(df, clusters, opt$output.plot)