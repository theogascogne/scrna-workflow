#!/usr/bin/env Rscript

require(tidyverse)
require(optparse)

# Get the path to the 'cellsnake' package
cellsnake_path <- system("python -c 'import cellsnake; print(cellsnake.__path__[0])'", intern = TRUE)
test_data_dir <- file.path(cellsnake_path, "scrna/workflow/tests/testData")

option_list <- list(
    optparse::make_option(c("--rds"),
        type = "character", default = NULL,
        help = "RDS file of marker data frame", metavar = "character"
    ),
    optparse::make_option(c("--output.xlsx.positive"),
        type = "character", default = NULL,
        help = "Excel table of positive markers", metavar = "character"
    ),
    optparse::make_option(c("--output.xlsx.all"),
        type = "character", default = NULL,
        help = "Excel table of all markers", metavar = "character"
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
    source(file.path(cellsnake_path,"scrna/workflow/scripts/scrna_workflow_helper_functions/scrna-marker-tables-functions.R"))
  },
  silent = TRUE
)

all_markers <- readRDS(file = opt$rds)

openxlsx::write.xlsx(all_markers, file = opt$output.xlsx.all)

openxlsx::write.xlsx(all_markers %>% filter(avg_log2FC > 0), file = opt$output.xlsx.positive)