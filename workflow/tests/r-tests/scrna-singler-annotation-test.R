#!/usr/bin/env Rscript

library(optparse)

# Get the path to the 'cellsnake' package
cellsnake_path <- system("python -c 'import cellsnake; print(cellsnake.__path__[0])'", intern = TRUE)

# Define the base path to the testData folder
test_data_dir <- file.path(cellsnake_path, "scrna/workflow/tests/testData")

# Define test options with dynamically constructed paths
option_list <- list(
  make_option(c("--rds"), type = "character", default = file.path(test_data_dir, "test/data/sample_seurat_object.rds"), help = "Processed rds file of a Seurat object", metavar = "character"),
  make_option(c("--output"), type = "character", default = file.path(test_data_dir, "results/pred_test.rds"), help = "Output prediction file", metavar = "character"),
  make_option(c("--reference"), type = "character", default = "HumanPrimaryCellAtlasData", help = "SingleR reference", metavar = "character"),
  make_option(c("--granulation"), type = "character", default = "label.main", help = "SingleR granulation level", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Ensure the directory for output file exists
dir.create(dirname(opt$output), recursive = TRUE, showWarnings = FALSE)

# Remove any pre-existing test files
file.remove(opt$output)

# Construct the command to run the original script
cmd <- paste(
  "Rscript", file.path(cellsnake_path, "scrna/workflow/scripts/scrna-singler-annotation.R"),
  "--rds", shQuote(opt$rds),
  "--output", shQuote(opt$output),
  "--reference", shQuote(opt$reference),
  "--granulation", shQuote(opt$granulation)
)

# Run the original script
cat("Running command:\n", cmd, "\n")
status <- system(cmd, intern = TRUE)
cat("Command output:\n", status, "\n")

# Define color codes for terminal
cyan <- "\033[0;36m"
red <- "\033[0;31m"
reset <- "\033[0m"

# Check if the expected output file exists and is not empty
if (file.exists(opt$output)) {
  # Try to read the RDS file to verify it contains valid data
  pred <- tryCatch({
    readRDS(opt$output)
  }, error = function(e) {
    NULL
  })
  
  # Check if the read data is valid
  if (!is.null(pred)) {
    cat(cyan, "SUCCESS:", reset, "Found and successfully loaded the output file:", opt$output, "\n", reset, sep = "")
  } else {
    cat(red, "ERROR:", reset, "Output file found but could not be loaded or is empty:", opt$output, "\n", reset, sep = "")
  }
} else {
  cat(red, "ERROR:", reset, "Expected output file missing:", opt$output, "\n", reset, sep = "")
}
