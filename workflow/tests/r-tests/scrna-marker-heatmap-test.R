#!/usr/bin/env Rscript

library(optparse)

# Get the path to the 'cellsnake' package
cellsnake_path <- system("python -c 'import cellsnake; print(cellsnake.__path__[0])'", intern = TRUE)

# Define the base path to the testData folder
test_data_dir <- file.path(cellsnake_path, "scrna/workflow/tests/testData")

# Define options
option_list <- list(
  make_option(c("--rds"), type = "character", default = file.path(test_data_dir, "analyses/processed/percent_mt~10/resolution~0.8/10X_17_029.rds"),
              help = "Seurat object RDS file", metavar = "character"),
  make_option(c("--xlsx"), type = "character", default = file.path(test_data_dir, "results/10X_17_029/percent_mt~10/resolution~0.8/table_positive-markers-seurat_clusters.xlsx"),
              help = "Excel table of markers", metavar = "character"),
  make_option(c("--output.plot"), type = "character", default = file.path(test_data_dir, "results/10X_17_029/percent_mt~10/resolution~0.8/plot_marker-heatmap-seurat_clusters.pdf"),
              help = "Output heatmap PDF", metavar = "character"),
  make_option(c("--output.average"), type = "character", default = file.path(test_data_dir, "results/10X_17_029/percent_mt~10/resolution~0.8/table_average-expression-seurat_clusters.xlsx"),
              help = "Output average expression table", metavar = "character"),
  make_option(c("--idents"), type = "character", default = "seurat_clusters",
              help = "Metadata column for clustering", metavar = "character")
)

# Parse options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Create necessary directories if they don't exist
dir.create(dirname(opt$output.plot), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(opt$output.average), recursive = TRUE, showWarnings = FALSE)

# Remove pre-existing test output
if (file.exists(opt$output.plot)) {
  file.remove(opt$output.plot)
}
if (file.exists(opt$output.average)) {
  file.remove(opt$output.average)
}

# Run the script
cmd <- paste(
  "Rscript", file.path(cellsnake_path, "scrna/workflow/scripts/scrna-marker-heatmap.R"),
  "--rds", shQuote(opt$rds),
  "--xlsx", shQuote(opt$xlsx),
  "--output.plot", shQuote(opt$output.plot),
  "--output.average", shQuote(opt$output.average),
  "--idents", shQuote(opt$idents)
)

# Capture system output
system_output <- system(cmd, intern = TRUE, ignore.stderr = FALSE)

# Print system output
cat("System Output:\n")
cat(system_output, sep = "\n")

# Validate output files
if (!file.exists(opt$output.plot)) {
  stop("FAIL: ", opt$output.plot, " was not created.")
} else {
  message("PASS: ", opt$output.plot, " was created successfully.")
}

if (!file.exists(opt$output.average)) {
  stop("FAIL: ", opt$output.average, " was not created.")
} else {
  message("PASS: ", opt$output.average, " was created successfully.")
}

message("Test completed successfully.")