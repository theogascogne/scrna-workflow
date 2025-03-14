#!/usr/bin/env Rscript

library(optparse)

# Get the path to the 'cellsnake' package
cellsnake_path <- system("python -c 'import cellsnake; print(cellsnake.__path__[0])'", intern = TRUE)

# Define the base path to the testData folder
test_data_dir <- file.path(cellsnake_path, "scrna/workflow/tests/testData")

# Define options
option_list <- list(
  make_option(c("--rds"), type = "character", default = file.path(test_data_dir, "analyses/markers/percent_mt~10/resolution~0.8/markers_10X_17_029-seurat_clusters.rds"),
              help = "RDS file of marker data frame", metavar = "character"),
  make_option(c("--output.xlsx.positive"), type = "character", default = file.path(test_data_dir, "results/10X_17_029/percent_mt~10/resolution~0.8/table_positive-markers-seurat_clusters.xlsx"),
              help = "Excel table of positive markers", metavar = "character"),
  make_option(c("--output.xlsx.all"), type = "character", default = file.path(test_data_dir, "results/10X_17_029/percent_mt~10/resolution~0.8/table_all-markers-seurat_clusters.xlsx"),
              help = "Excel table of all markers", metavar = "character")
)

# Parse options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Create necessary directories if they don't exist
dir.create(dirname(opt$output.xlsx.positive), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(opt$output.xlsx.all), recursive = TRUE, showWarnings = FALSE)

# Remove pre-existing test outputs
if (file.exists(opt$output.xlsx.positive)) file.remove(opt$output.xlsx.positive)
if (file.exists(opt$output.xlsx.all)) file.remove(opt$output.xlsx.all)

# Run the command
cmd <- paste(
  "Rscript", file.path(cellsnake_path, "scrna/workflow/scripts/scrna-marker-tables.R"),
  "--rds", shQuote(opt$rds),
  "--output.xlsx.positive", shQuote(opt$output.xlsx.positive),
  "--output.xlsx.all", shQuote(opt$output.xlsx.all)
)

system_output <- system(cmd, intern = TRUE, ignore.stderr = FALSE)

# Print system output
cat("System Output:\n")
cat(system_output, sep = "\n")

# Check if output files were created
if (!file.exists(opt$output.xlsx.positive)) {
  stop("FAIL: ", opt$output.xlsx.positive, " was not created.")
} else {
  message("PASS: ", opt$output.xlsx.positive, " was created successfully.")
}

if (!file.exists(opt$output.xlsx.all)) {
  stop("FAIL: ", opt$output.xlsx.all, " was not created.")
} else {
  message("PASS: ", opt$output.xlsx.all, " was created successfully.")
}

message("Test completed successfully.")