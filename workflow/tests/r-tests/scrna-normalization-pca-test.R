#!/usr/bin/env Rscript

library(optparse)

# Get the path to the 'cellsnake' package
cellsnake_path <- system("python -c 'import cellsnake; print(cellsnake.__path__[0])'", intern = TRUE)

# Define the base path to the testData folder
test_data_dir <- file.path(cellsnake_path, "scrna/workflow/tests/testData")

# Define options with dynamically constructed paths
option_list <- list(
  make_option(c("--scale.factor"), type = "integer", default = 10000, help = "Scale factor [default= %default]", metavar = "character"),
  make_option(c("--nfeatures"), type = "integer", default = 2000, help = "Highly variable features [default= %default]", metavar = "integer"),
  make_option(c("--variable.selection.method"), type = "character", default = "vst", help = "Variable selection method [default= %default]", metavar = "character"),
  make_option(c("--rds"), type = "character", default = file.path(test_data_dir, "test/data/sample.rds"), help = "Test RDS file path", metavar = "character"),
  make_option(c("--normalization.method"), type = "character", default = "LogNormalize", help = "Normalization method [default= %default]", metavar = "character"),
  make_option(c("--doublet.filter"), action = "store_true", default = FALSE, help = "Test Doublet filter [default= %default]"),
  make_option(c("--integration"), action = "store_true", default = FALSE, help = "Run with integration [default= %default]"),
  make_option(c("--umap"), action = "store_true", default = TRUE, help = "Test UMAP generation [default= %default]"),
  make_option(c("--tsne"), action = "store_true", default = TRUE, help = "Test t-SNE generation [default= %default]"),
  make_option(c("--resolution"), type = "character", default = "0.8", help = "Test resolution [default= %default]", metavar = "character"),
  make_option(c("--output.rds"), type = "character", default = file.path(test_data_dir, "results/output_test.rds"), help = "Output RDS test file [default= %default]", metavar = "character"),
  make_option(c("--cpu"), type = "integer", default = 2, help = "Number of CPUs for test [default= %default]", metavar = "character"),
  make_option(c("--reference"), type = "character", default = "HumanPrimaryCellAtlasData", help = "SingleR reference for annotation [default= %default]", metavar = "character")
)

# Parse options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Ensure necessary directories exist
dir.create(dirname(opt$output.rds), recursive = TRUE, showWarnings = FALSE)

# Remove any pre-existing test files
file.remove(opt$output.rds)

# Construct the command to run the original script
cmd <- paste(
  "Rscript", file.path(cellsnake_path, "scrna/workflow/scripts/scrna-normalization-pca.R"),
  "--scale.factor", opt$scale.factor,
  "--nfeatures", opt$nfeatures,
  "--variable.selection.method", opt$variable.selection.method,
  "--rds", shQuote(opt$rds),
  "--normalization.method", opt$normalization.method,
  ifelse(opt$doublet.filter, "--doublet.filter", ""),
  ifelse(opt$integration, "--integration", ""),
  ifelse(opt$umap, "--umap", ""),
  ifelse(opt$tsne, "--tsne", ""),
  "--resolution", opt$resolution,
  "--output.rds", shQuote(opt$output.rds),
  "--cpu", opt$cpu,
  "--reference", shQuote(opt$reference)
)

# Run the original script
cat("Running command:\n", cmd, "\n")
status <- system(cmd, intern = TRUE)
cat("Command output:\n", status, "\n")

# Define color codes for terminal
cyan <- "\033[0;36m"
red <- "\033[0;31m"
reset <- "\033[0m"

# Verify expected output
expected_outputs <- list(opt$output.rds)

for (file in expected_outputs) {
  if (file.exists(file)) {
    cat(cyan, "SUCCESS:", reset, "Found expected output file:", file, "\n", reset, sep = "")
  } else {
    cat(red, "ERROR:", reset, "Expected output file missing:", file, "\n", reset, sep = "")
  }
}