#!/usr/bin/env Rscript
require(optparse)
require(tidyverse)

# Get the path to the 'cellsnake' package
cellsnake_path <- system("python -c 'import cellsnake; print(cellsnake.__path__[0])'", intern = TRUE)

option_list <- list(
  optparse::make_option(c("--xlsx"),
                        type = "character", default = NULL,
                        help = "Excel table of markers", metavar = "character"
  ),
  optparse::make_option(c("--output.rds"),
                        type = "character", default = NULL,
                        help = "Output RDS file name", metavar = "character"
  ),
  optparse::make_option(c("--output.kegg"),
                        type = "character", default = NULL,
                        help = "Output KEGG enrichment excel file name", metavar = "character"
  ),
  optparse::make_option(c("--output.mkegg"),
                        type = "character", default = NULL,
                        help = "Output metabolic KEGG enrichment excel file name", metavar = "character"
  ),
  optparse::make_option(c("--output.gse"),
                        type = "character", default = NULL,
                        help = "Output KEGG GSEA excel file name", metavar = "character"
  ),
  optparse::make_option(c("--output.mgse"),
                        type = "character", default = NULL,
                        help = "Output metabolic KEGG GSEA excel file name", metavar = "character"
  ),
  optparse::make_option(c("--mapping"),
                        type = "character", default = "org.Hs.eg.db",
                        help = "Mapping database", metavar = "character"
  ),
  optparse::make_option(c("--organism"),
                        type = "character", default = "hsa",
                        help = "Organism code for KEGG (e.g., hsa for human)", metavar = "character"
  ),
  optparse::make_option(c("--pval"),
                        type = "double", default = 0.05,
                        help = "P-value threshold", metavar = "character"
  ),
  optparse::make_option(c("--logfc.treshold"),
                        type = "double", default = 1,
                        help = "LogFC threshold", metavar = "character"
  )
)

if (!exists("opt")) {
  opt_parser <- OptionParser(option_list = option_list)
  opt <- parse_args(opt_parser)
}

try({
  if (!requireNamespace(opt$mapping, quietly = TRUE)) {
    BiocManager::install(opt$mapping)
  }
})

require(opt$mapping, character.only = TRUE)

if (is.null(opt$xlsx)) {
  optparse::print_help(opt_parser)
  stop("Arguments must be supplied", call. = FALSE)
}

# Source helper functions
try({
  source(file.path(cellsnake_path, "scrna/workflow/scripts/scrna_workflow_helper_functions/scrna-kegg-functions.R"))
}, silent = TRUE)

All_Features <- openxlsx::read.xlsx(opt$xlsx)

# Use helper function directly, passing mapping_db and organism_code explicitly
all_kegg_results <- All_Features %>%
  split(.$cluster) %>%
  purrr::map(~ function_enrichment_kegg_singlecell(., p = opt$pval, f = opt$logfc.treshold, mapping_db = get(opt$mapping), organism_code = opt$organism))

saveRDS(all_kegg_results, opt$output.rds)

# Use helper function to save the results by index
save_kegg_results(all_kegg_results, 1, opt$output.kegg)  # enrichKEGG
save_kegg_results(all_kegg_results, 2, opt$output.gse)   # gseKEGG
save_kegg_results(all_kegg_results, 3, opt$output.mkegg) # enrichMKEGG
save_kegg_results(all_kegg_results, 4, opt$output.mgse)  # gseMKEGG