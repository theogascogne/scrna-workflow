#!/usr/bin/env Rscript

library(optparse)
require(tidyverse)
library(AnnotationDbi)

# Get the path to the 'cellsnake' package
cellsnake_path <- system("python -c 'import cellsnake; print(cellsnake.__path__[0])'", intern = TRUE)

# Define the base path to the testData folder
test_data_dir <- file.path(cellsnake_path, "scrna/workflow/tests/testData")


option_list <- list(
  optparse::make_option(c("--xlsx"),
                        type = "character", default = NULL,
                        help = "Excel table of markers", metavar = "character"
  ),
  optparse::make_option(c("--output.rds"),
                        type = "character", default = NULL,
                        help = "Output RDS file name", metavar = "character"
  ),
  optparse::make_option(c("--output.go"),
                        type = "character", default = NULL,
                        help = "Output kegg excel file name", metavar = "character"
  ),
  optparse::make_option(c("--output.gse"),
                        type = "character", default = NULL,
                        help = "Output gse excel file name", metavar = "character"
  ),
  optparse::make_option(c("--mapping"),
                        type = "character", default = "org.Hs.eg.db",
                        help = "Mapping", metavar = "character"
  ),
  optparse::make_option(c("--pval"),
                        type = "double", default = 0.05,
                        help = "P value treshold [default= %default]", metavar = "character"
  ),
  optparse::make_option(c("--logfc.treshold"),
                        type = "double", default = 1.5,
                        help = "LogFC [default= %default]", metavar = "character"
  )
)

try({
  if (!requireNamespace(opt$mapping, quietly = TRUE)) {
    BiocManager::install(opt$mapping, update = TRUE)
  }
})

# parse options
if (!exists("opt")) {
  opt_parser <- OptionParser(option_list = option_list)
  opt <- parse_args(opt_parser)
}
require(opt$mapping, character.only = T)

if (is.null(opt$xlsx)) {
  optparse::print_help(opt_parser)
  stop("Arguments must be supplied", call. = FALSE)
}

try(
  {
    source(file.path(cellsnake_path,"scrna/workflow/scripts/scrna_workflow_helper_functions/scrna-go-analysis-functions.R"))
  },
  silent = TRUE
)

# Read markers table
All_Features <- openxlsx::read.xlsx(opt$xlsx)

# Run enrichment for each cluster using helper
all_go_results <- All_Features %>%
  split(.$cluster) %>%
  purrr::map(~ function_enrichment_go_singlecell(., p = opt$pval, f = opt$logfc.treshold))

# Save results using helpers
saveRDS(all_go_results, opt$output.rds)
save_go_results(all_go_results, opt$output.go)
save_gse_results(all_go_results, opt$output.gse)