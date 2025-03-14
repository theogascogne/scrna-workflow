#!/usr/bin/env Rscript

library(optparse)

# Get the path to the 'cellsnake' package
cellsnake_path <- system("python -c 'import cellsnake; print(cellsnake.__path__[0])'", intern = TRUE)

# Define the base path to the testData folder
test_data_dir <- file.path(cellsnake_path, "scrna/workflow/tests/testData")

# Define options with dynamically constructed paths
option_list <- list(
  make_option(c("--rds"), type = "character", default = file.path(test_data_dir, "sample_seurat_object.rds"), help = "Processed rds file of a Seurat object", metavar = "character"),
  make_option(c("--sampleid"), type = "character", default = "TestSample", help = "Sample ID for testing", metavar = "character"),
  make_option(c("--fplot"), type = "character", default = file.path(test_data_dir, "results/fplot_test.pdf"), help = "nFeature plot", metavar = "character"),
  make_option(c("--cplot"), type = "character", default = file.path(test_data_dir, "results/cplot_test.pdf"), help = "nCount plot", metavar = "character"),
  make_option(c("--mtplot"), type = "character", default = file.path(test_data_dir, "results/mtplot_test.pdf"), help = "Percent MT plot", metavar = "character"),
  make_option(c("--rpplot"), type = "character", default = file.path(test_data_dir, "results/rpplot_test.pdf"), help = "Ribo plot", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Ensure the directory for output files exists
cat("Creating directories if they do not exist...\n")
dir.create(dirname(opt$fplot), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(opt$cplot), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(opt$mtplot), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(opt$rpplot), recursive = TRUE, showWarnings = FALSE)

# Remove any pre-existing test files
file.remove(opt$fplot)
file.remove(opt$cplot)
file.remove(opt$mtplot)
file.remove(opt$rpplot)

# Construct the command to run the original script
cmd <- paste(
  "Rscript", file.path(cellsnake_path, "scrna/workflow/scripts/scrna-technicals.R"),
  "--rds", shQuote(opt$rds),
  "--sampleid", shQuote(opt$sampleid),
  "--fplot", shQuote(opt$fplot),
  "--cplot", shQuote(opt$cplot),
  "--mtplot", shQuote(opt$mtplot),
  "--rpplot", shQuote(opt$rpplot)
)

# Run the original script
cat("Running command:\n", cmd, "\n")
status <- system(cmd, intern = TRUE)
cat("Command output:\n", status, "\n")

# Define color codes for terminal output
cyan <- "\033[0;36m"
red <- "\033[0;31m"
reset <- "\033[0m"

# List of final expected output files
final_outputs <- list(opt$fplot, opt$cplot, opt$mtplot, opt$rpplot)

# Check for each expected output file
for (file in final_outputs) {
  if (file.exists(file)) {
    cat(cyan, "SUCCESS:", reset, "Found expected output file:", file, "\n", reset, sep = "")
  } else {
    cat(red, "ERROR:", reset, "Expected output file missing:", file, "\n", reset, sep = "")
  }
}
