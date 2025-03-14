#!/usr/bin/env Rscript

library(optparse)

# Get the path to the 'cellsnake' package
cellsnake_path <- system("python -c 'import cellsnake; print(cellsnake.__path__[0])'", intern = TRUE)

# Define the base path to the testData folder
test_data_dir <- file.path(cellsnake_path, "scrna/workflow/tests/testData")

# Define options with dynamically constructed paths
option_list <- list(
  make_option(c("--xlsx"), type = "character", default = file.path(test_data_dir, "10X_17_029/percent_mt~10/resolution~0.8/table_positive-markers-seurat_clusters.xlsx"),
              help = "Excel table of markers for input", metavar = "character"),
  make_option(c("--output.plot"), type = "character", default = file.path(test_data_dir, "10X_17_029/percent_mt~10/resolution~0.8/summarized_markers-for-seurat_clusters.pdf"),
              help = "Output plot file", metavar = "character")
)

# Parse options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Create necessary directories if they don't exist
cat("Creating directories if they do not exist...\n")
dir.create(dirname(opt$output.plot), recursive = TRUE, showWarnings = FALSE)

# Remove pre-existing test outputs
if (file.exists(opt$output.plot)) file.remove(opt$output.plot)

# Construct the command to run the original script
cmd <- paste(
  "Rscript", file.path(cellsnake_path, "scrna/workflow/scripts/scrna-top-marker-plot.R"),
  "--xlsx", shQuote(opt$xlsx),
  "--output.plot", shQuote(opt$output.plot)
)

# Run the command
cat("Running command:\n", cmd, "\n")
system_output <- system(cmd, intern = TRUE, ignore.stderr = FALSE)

# Print system output
cat("System Output:\n")
cat(system_output, sep = "\n")

# Check if output file was created
if (!file.exists(opt$output.plot)) {
  stop("FAIL: ", opt$output.plot, " was not created.")
} else {
  message("PASS: ", opt$output.plot, " was created successfully.")
}

message("Test completed successfully.")