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
              default = file.path(test_data_dir, "analyses/go/percent_mt~10/resolution~0.8/10X_17_029-seurat_clusters-go.rds"),
              help = "Output RDS file name", metavar = "character"),
  make_option(c("--output.go"), type = "character", 
              default = file.path(test_data_dir, "results/10X_17_029/percent_mt~10/resolution~0.8/enrichment_analysis/table_GO-enrichment-seurat_clusters.xlsx"),
              help = "Output GO enrichment file", metavar = "character"),
  make_option(c("--output.gse"), type = "character", 
              default = file.path(test_data_dir, "results/10X_17_029/percent_mt~10/resolution~0.8/enrichment_analysis/table_GO-geneset_enrichment-seurat_clusters.xlsx"),
              help = "Output GSE file", metavar = "character"),
  make_option(c("--mapping"), type = "character", default = "org.Hs.eg.db",
              help = "Gene mapping database", metavar = "character")
)

# Parse options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Create necessary directories if they don't exist
dir.create(dirname(opt$output.rds), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(opt$output.go), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(opt$output.gse), recursive = TRUE, showWarnings = FALSE)

# Remove pre-existing test outputs
if (file.exists(opt$output.rds)) file.remove(opt$output.rds)
if (file.exists(opt$output.go)) file.remove(opt$output.go)
if (file.exists(opt$output.gse)) file.remove(opt$output.gse)

# Run the command
cmd <- paste(
  "Rscript", file.path(cellsnake_path, "scrna/workflow/scripts/scrna-go_analysis.R"),
  "--xlsx", shQuote(opt$xlsx),
  "--output.rds", shQuote(opt$output.rds),
  "--output.go", shQuote(opt$output.go),
  "--output.gse", shQuote(opt$output.gse),
  "--mapping", opt$mapping
)

# Capture system output
system_output <- system(cmd, intern = TRUE, ignore.stderr = FALSE)

# Print system output
cat("System Output:\n")
cat(system_output, sep = "\n")

# Check if output files were created
if (!file.exists(opt$output.rds)) {
  stop("FAIL: ", opt$output.rds, " was not created.")
} else {
  message("PASS: ", opt$output.rds, " was created successfully.")
}

if (!file.exists(opt$output.go)) {
  stop("FAIL: ", opt$output.go, " was not created.")
} else {
  message("PASS: ", opt$output.go, " was created successfully.")
}

if (!file.exists(opt$output.gse)) {
  stop("FAIL: ", opt$output.gse, " was not created.")
} else {
  message("PASS: ", opt$output.gse, " was created successfully.")
}

message("Test completed successfully.")