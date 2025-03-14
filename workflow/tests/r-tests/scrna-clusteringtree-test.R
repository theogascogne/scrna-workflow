#!/usr/bin/env Rscript

library(optparse)

# Get the path to the 'cellsnake' package
cellsnake_path <- system("python -c 'import cellsnake; print(cellsnake.__path__[0])'", intern = TRUE)

# Define the base path to the testData folder
test_data_dir <- file.path(cellsnake_path, "scrna/workflow/tests/testData")

# Define options with default test values
option_list <- list(
  make_option(c("--scale.factor"), type = "integer", default = 10000, help = "Scale factor [default= %default]", metavar = "character"),
  make_option(c("--nfeatures"), type = "integer", default = 2000, help = "Highly variable features [default= %default]", metavar = "integer"),
  make_option(c("--variable.selection.method"), type = "character", default = "vst", help = "Variable selection method [default= %default]", metavar = "character"),
  make_option(c("--rds"), type = "character", 
              default = file.path(test_data_dir, "analyses/processed/percent_mt~10/resolution~0.8/10X_17_028.rds"), 
              help = "Test RDS file path", metavar = "character"),
  make_option(c("--normalization.method"), type = "character", default = "LogNormalize", help = "Normalization method [default= %default]", metavar = "character"),
  make_option(c("--integration"), action = "store_true", default = FALSE, help = "Run with integration [default= %default]"),
  make_option(c("--clplot"), type = "character", 
              default = file.path(test_data_dir, "test/results/clustree_test.pdf"), 
              help = "Test clustree output path", metavar = "character"),
  make_option(c("--jeplot"), type = "character", 
              default = file.path(test_data_dir, "test/results/jackandelbow_test.pdf"), 
              help = "Test Jack and Elbow output path", metavar = "character"),
  make_option(c("--hvfplot"), type = "character", 
              default = file.path(test_data_dir, "test/results/variable_features_test.pdf"), 
              help = "Test HVF output path", metavar = "character"),
  make_option(c("--heplot"), type = "character", 
              default = file.path(test_data_dir, "test/results/dimheatmap_test.pdf"), 
              help = "Test heatmap output path", metavar = "character")
)

# Parse options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Ensure necessary directories exist
dir.create(dirname(opt$clplot), recursive = TRUE, showWarnings = FALSE)

# Remove any pre-existing test files
file.remove(opt$clplot, opt$jeplot, opt$hvfplot, opt$heplot)

# Construct the command for the original script
cmd <- paste(
  "Rscript", file.path(cellsnake_path, "scrna/workflow/scripts/scrna-clusteringtree.R"),
  "--scale.factor", opt$scale.factor,
  "--nfeatures", opt$nfeatures,
  "--variable.selection.method", opt$variable.selection.method,
  "--rds", shQuote(opt$rds),
  "--normalization.method", opt$normalization.method,
  ifelse(opt$integration, "--integration", ""),
  "--clplot", shQuote(opt$clplot),
  "--jeplot", shQuote(opt$jeplot),
  "--hvfplot", shQuote(opt$hvfplot),
  "--heplot", shQuote(opt$heplot)
)

# Run the original script
cat("Running command:\n", cmd, "\n")
status <- system(cmd, intern = TRUE)
cat("Command output:\n", status, "\n")

# Define color codes for terminal output
cyan <- "\033[0;36m"
red <- "\033[0;31m"
reset <- "\033[0m"

# Verify the outputs
expected_outputs <- list(opt$clplot, opt$jeplot, opt$hvfplot, opt$heplot)

for (file in expected_outputs) {
  if (file.exists(file)) {
    cat(cyan, "SUCCESS:", reset, "Found expected output file:", file, "\n", reset, sep = "")
  } else {
    cat(red, "ERROR:", reset, "Expected output file missing:", file, "\n", reset, sep = "")
  }
}