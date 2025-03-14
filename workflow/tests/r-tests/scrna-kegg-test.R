#!/usr/bin/env Rscript

library(optparse)

# Get the path to the 'cellsnake' package
cellsnake_path <- system("python -c 'import cellsnake; print(cellsnake.__path__[0])'", intern = TRUE)

# Define the base path to the testData folder
test_data_dir <- file.path(cellsnake_path, "scrna/workflow/tests/testData")

# Define options
option_list <- list(
  make_option(c("--xlsx"), type = "character", 
              default = file.path(test_data_dir, "results/10X_17_029/percent_mt~10/resolution~0.8/table_all-markers-seurat_clusters.xlsx"),
              help = "Excel table of markers", metavar = "character"),
  make_option(c("--output.rds"), type = "character", 
              default = file.path(test_data_dir, "analyses/kegg/percent_mt~10/resolution~0.8/10X_17_029-seurat_clusters-kegg.rds"),
              help = "Output RDS file name", metavar = "character"),
  make_option(c("--output.kegg"), type = "character", 
              default = file.path(test_data_dir, "results/10X_17_029/percent_mt~10/resolution~0.8/enrichment_analysis/table_KEGG-enrichment-seurat_clusters.xlsx"),
              help = "Output KEGG enrichment file", metavar = "character"),
  make_option(c("--output.mkegg"), type = "character", 
              default = file.path(test_data_dir, "results/10X_17_029/percent_mt~10/resolution~0.8/enrichment_analysis/table_KEGG-module_enrichment-seurat_clusters.xlsx"),
              help = "Output KEGG module enrichment file", metavar = "character"),
  make_option(c("--output.gse"), type = "character", 
              default = file.path(test_data_dir, "results/10X_17_029/percent_mt~10/resolution~0.8/enrichment_analysis/table_KEGG-geneset_enrichment-seurat_clusters.xlsx"),
              help = "Output KEGG gene set enrichment file", metavar = "character"),
  make_option(c("--output.mgse"), type = "character", 
              default = file.path(test_data_dir, "results/10X_17_029/percent_mt~10/resolution~0.8/enrichment_analysis/table_KEGG-module_geneset_enrichment-seurat_clusters.xlsx"),
              help = "Output KEGG module gene set enrichment file", metavar = "character"),
  make_option(c("--mapping"), type = "character", default = "org.Hs.eg.db",
              help = "Gene mapping database", metavar = "character"),
  make_option(c("--organism"), type = "character", default = "hsa",
              help = "Organism code (KEGG-compatible)", metavar = "character")
)

# Parse options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Create necessary directories if they don't exist
dir.create(dirname(opt$output.rds), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(opt$output.kegg), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(opt$output.mkegg), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(opt$output.gse), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(opt$output.mgse), recursive = TRUE, showWarnings = FALSE)

# Remove pre-existing test outputs
output_files <- c(opt$output.rds, opt$output.kegg, opt$output.mkegg, opt$output.gse, opt$output.mgse)
file.remove(output_files[file.exists(output_files)])

# Run the command
cmd <- paste(
  "Rscript", file.path(cellsnake_path, "scrna/workflow/scripts/scrna-kegg.R"),
  "--xlsx", shQuote(opt$xlsx),
  "--output.rds", shQuote(opt$output.rds),
  "--output.kegg", shQuote(opt$output.kegg),
  "--output.mkegg", shQuote(opt$output.mkegg),
  "--output.gse", shQuote(opt$output.gse),
  "--output.mgse", shQuote(opt$output.mgse),
  "--mapping", opt$mapping,
  "--organism", opt$organism
)

# Capture system output
system_output <- system(cmd, intern = TRUE, ignore.stderr = FALSE)

# Print system output
cat("System Output:\n")
cat(system_output, sep = "\n")

# Validate output files
for (file in output_files) {
  if (!file.exists(file)) {
    stop("FAIL: ", file, " was not created.")
  } else {
    message("PASS: ", file, " was created successfully.")
  }
}

message("Test completed successfully.")