#!/usr/bin/env Rscript

library(optparse)

# Get the path to the 'cellsnake' package
cellsnake_path <- system("python -c 'import cellsnake; print(cellsnake.__path__[0])'", intern = TRUE)

# Define the base path to the testData folder
test_data_dir <- file.path(cellsnake_path, "scrna/workflow/tests/testData")

# Define options
option_list <- list(
  make_option(c("--rds"), type = "character", 
              default = file.path(test_data_dir, "analyses/processed/percent_mt~10/resolution~0.8/10X_17_029.rds"),
              help = "Processed RDS file of a Seurat object", metavar = "character"),
  make_option(c("--logfc.threshold"), type = "double", default = 0.25,
              help = "Log fold change threshold [default= %default]", metavar = "double"),
  make_option(c("--test.use"), type = "character", default = "wilcox",
              help = "Statistical test to use [default= %default]", metavar = "character"),
  make_option(c("--output.rds"), type = "character", 
              default = file.path(test_data_dir, "analyses/markers/percent_mt~10/resolution~0.8/markers_10X_17_029-seurat_clusters.rds"),
              help = "Output RDS file containing marker genes", metavar = "character"),
  make_option(c("--idents"), type = "character", default = "seurat_clusters",
              help = "Metadata column name for marker analysis", metavar = "character")
)

# Parse options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Create necessary directories if they don't exist
dir.create(dirname(opt$output.rds), recursive = TRUE, showWarnings = FALSE)

# Remove pre-existing test output
if (file.exists(opt$output.rds)) file.remove(opt$output.rds)

# Run the command
cmd <- paste(
  "Rscript", file.path(cellsnake_path, "scrna/workflow/scripts/scrna-find-markers.R"),
  "--rds", shQuote(opt$rds),
  "--idents", opt$idents,
  "--logfc.threshold", opt$logfc.threshold,
  "--test.use", opt$test.use,
  "--output.rds", shQuote(opt$output.rds)
)

# Capture system output
system_output <- system(cmd, intern = TRUE, ignore.stderr = FALSE)

# Print system output
cat("System Output:\n")
cat(system_output, sep = "\n")

# Check if output file was created
if (!file.exists(opt$output.rds)) {
  stop("FAIL: ", opt$output.rds, " was not created.")
} else {
  message("PASS: ", opt$output.rds, " was created successfully.")
}

message("Test completed successfully.")