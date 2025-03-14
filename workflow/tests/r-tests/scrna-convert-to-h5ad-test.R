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
              help = "Seurat object RDS file", metavar = "character"),
  make_option(c("--output"), type = "character", 
              default = file.path(test_data_dir, "analyses/h5ad/percent_mt~10/resolution~0.8/10X_17_029.h5ad"),
              help = "Output h5ad file", metavar = "character")
)

# Parse options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Create necessary directories if they don't exist
dir.create(dirname(opt$output), recursive = TRUE, showWarnings = FALSE)

# Remove pre-existing test output
if (file.exists(opt$output)) {
  file.remove(opt$output)
}

# Run the script
cmd <- paste(
  "Rscript", file.path(cellsnake_path, "scrna/workflow/scripts/scrna-convert-to-h5ad.R"),
  "--rds", shQuote(opt$rds),
  "--output", shQuote(opt$output)
)

system_output <- system(cmd, intern = TRUE, ignore.stderr = FALSE)

# Print system output
cat("System Output:\n")
cat(system_output, sep = "\n")

# Validate output file
if (!file.exists(opt$output)) {
  stop("FAIL: ", opt$output, " was not created.")
} else {
  message("PASS: ", opt$output, " was created successfully.")
}

message("Test completed successfully.")