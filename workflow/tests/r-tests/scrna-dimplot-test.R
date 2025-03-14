#!/usr/bin/env Rscript

library(optparse)

# Get the path to the 'cellsnake' package
cellsnake_path <- system("python -c 'import cellsnake; print(cellsnake.__path__[0])'", intern = TRUE)

# Define the base path to the testData folder
test_data_dir <- file.path(cellsnake_path, "scrna/workflow/tests/testData")

# Define color codes for terminal output
cyan <- "\033[0;36m"
red <- "\033[0;31m"
reset <- "\033[0m"

# Define options for testing, using the defaults from the shell command
option_list <- list(
  make_option(c("--rds"), type = "character", default = file.path(test_data_dir, "analyses/processed/percent_mt~10/resolution~0.8/10X_17_028.rds"),
              help = "Path to Seurat object RDS file", metavar = "character"),
  make_option(c("--reduction.type"), type = "character", default = "tsne", help = "Reduction type (umap/tsne)", metavar = "character"),
  make_option(c("--pdfplot"), type = "character", default = file.path(test_data_dir, "results/10X_17_028/percent_mt~10/resolution~0.8/plot_dimplot_tsne-seurat_clusters.pdf"),
              help = "Path to save PDF plot", metavar = "character"),
  make_option(c("--htmlplot"), type = "character", default = file.path(test_data_dir, "results/10X_17_028/percent_mt~10/resolution~0.8/plot_dimplot_tsne-seurat_clusters.html"),
              help = "Path to save HTML plot", metavar = "character"),
  make_option(c("--csv"), type = "character", default = NULL, help = "CSV metadata file (for SingleR annotation)", metavar = "character"),
  make_option(c("--idents"), type = "character", default = "seurat_clusters", help = "Metadata column for identities", metavar = "character"),
  make_option(c("--percentage"), type = "double", default = 5, help = "Minimum cluster percentage", metavar = "double"),
  make_option(c("--labels"), action = "store_true", default = TRUE, help = "Add labels to plot")  # Assuming you want labels, as they are in the shell command
)

# Parse options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check and create output directories if they don't exist
dir.create(dirname(opt$pdfplot), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(opt$htmlplot), recursive = TRUE, showWarnings = FALSE)

# Remove old files for clean testing
file.remove(opt$pdfplot, opt$htmlplot)

# Check if it's the 'singler' case (SingleR annotation) and adjust the logic
if (!is.null(opt$csv)) {
  # Handle SingleR case (with CSV file)
  print("Running SingleR-based dimplot test...")
  idents <- "singler"  # When running the SingleR-based case
} else {
  # Handle default Seurat case (seurat_clusters)
  print("Running Seurat-based dimplot test...")
  idents <- opt$idents  # Use Seurat clusters by default
}

# Build the command to run the script you're testing
cmd <- paste(
  "Rscript", file.path(cellsnake_path, "scrna/workflow/scripts/scrna-dimplot.R"),
  "--rds", shQuote(opt$rds),
  "--reduction.type", opt$reduction.type,
  "--pdfplot", shQuote(opt$pdfplot),
  "--htmlplot", shQuote(opt$htmlplot),
  "--idents", idents,
  "--percentage", opt$percentage,
  if (opt$labels) "--labels" else ""
)

# Capture the output and error messages of the system call
system_output <- system(cmd, intern = TRUE, ignore.stderr = FALSE)

# Print the output directly to the console for more comprehensive error reporting
cat("System Output:\n")
cat(system_output, sep = "\n")

# Define a function to check the creation of files
check_output_with_retry <- function(filename, max_attempts = 5, interval = 10) {
  attempts <- 0
  while (!file.exists(filename)) {
    attempts <- attempts + 1
    if (attempts > max_attempts) {
      cat(red, "FAIL:", filename, "was not created after", max_attempts * interval, "seconds.\n", reset, sep = "")
      stop(red, "FAIL:", filename, "was not created after", max_attempts * interval, "seconds.\n", reset, sep = "")
    }
    Sys.sleep(interval)
  }
  cat(cyan, "PASS:", filename, "was created.\n", reset, sep = "")
}

# Check output files
check_output_with_retry(opt$pdfplot)
check_output_with_retry(opt$htmlplot)

cat(cyan, "All tests passed successfully.\n", reset, sep = "")