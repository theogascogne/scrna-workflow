#!/usr/bin/env Rscript

library(optparse)

# Get the path to the 'cellsnake' package
cellsnake_path <- system("python -c 'import cellsnake; print(cellsnake.__path__[0])'", intern = TRUE)

# Define the base path to the testData folder
test_data_dir <- file.path(cellsnake_path, "scrna/workflow/tests/testData")

# Define options with dynamically constructed paths
option_list <- list(
  make_option(c("--rds"), type = "character", default = file.path(test_data_dir, "test/data/sample_seurat_object.rds"),
              help = "Processed rds file of a Seurat object", metavar = "character"),
  make_option(c("--sampleid"), type = "character", default = "TestSample", help = "Sample ID for testing", metavar = "character"),
  make_option(c("--idents"), type = "character", default = "seurat_clusters", help = "Meta data column name for marker analysis", metavar = "character"),
  make_option(c("--ccplot"), type = "character", default = file.path(test_data_dir, "test/results/ccplot_test.pdf"), help = "Cell cluster count plot", metavar = "character"),
  make_option(c("--ccbarplot"), type = "character", default = file.path(test_data_dir, "test/results/ccbarplot_test.pdf"), help = "Cell cluster bar plot", metavar = "character"),
  make_option(c("--htmlplot"), type = "character", default = file.path(test_data_dir, "test/results/htmlplot_test.pdf"), help = "Cell cluster html plot", metavar = "character"),
  make_option(c("--xlsx"), type = "character", default = file.path(test_data_dir, "test/results/metrics_test.xlsx"), help = "Metrics table output", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Ensure the directory for output files exists
dir.create(dirname(opt$ccplot), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(opt$ccbarplot), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(opt$htmlplot), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(opt$xlsx), recursive = TRUE, showWarnings = FALSE)

# Remove any pre-existing test files
file.remove(opt$ccplot)
file.remove(opt$ccbarplot)
file.remove(opt$htmlplot)
file.remove(opt$xlsx)

# Construct the command to run the original script
cmd <- paste(
  "Rscript", file.path(cellsnake_path, "scrna/workflow/scripts/scrna-metrics.R"),
  "--rds", shQuote(opt$rds),
  "--sampleid", shQuote(opt$sampleid),
  "--idents", shQuote(opt$idents),
  "--ccplot", shQuote(opt$ccplot),
  "--ccbarplot", shQuote(opt$ccbarplot),
  "--htmlplot", shQuote(opt$htmlplot),
  "--xlsx", shQuote(opt$xlsx)
)

# Run the original script
cat("Running command:\n", cmd, "\n")
status <- system(cmd, intern = TRUE)
cat("Command output:\n", status, "\n")

# Define color codes for terminal
cyan <- "\033[0;36m"
red <- "\033[0;31m"
reset <- "\033[0m"

# List of final expected output files
final_outputs <- list(opt$ccplot, opt$ccbarplot, opt$htmlplot, opt$xlsx)

# Check for each expected output file
for (file in final_outputs) {
  if (file.exists(file)) {
    cat(cyan, "SUCCESS:", reset, "Found expected output file:", file, "\n", reset, sep = "")
  } else {
    cat(red, "ERROR:", reset, "Expected output file missing:", file, "\n", reset, sep = "")
  }
}
