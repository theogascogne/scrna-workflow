#!/usr/bin/env Rscript

library(optparse)

# Get the path to the 'cellsnake' package
cellsnake_path <- system("python -c 'import cellsnake; print(cellsnake.__path__[0])'", intern = TRUE)

# Define the base path to the testData folder
test_data_dir <- file.path(cellsnake_path, "scrna/workflow/tests/testData")

# Define options with dynamically constructed paths
option_list <- list(
  make_option(c("--rds"), type = "character", default = file.path(test_data_dir, "analyses/processed/percent_mt~10/resolution~0.8/10X_17_029.rds"),
              help = "Processed rds file of a Seurat object", metavar = "character"),
  make_option(c("--output.dir"), type = "character", default = file.path(test_data_dir, "results/10X_17_029/percent_mt~10/resolution~0.8/trajectory/"),
              help = "Output directory for trajectory plots", metavar = "character"),
  make_option(c("--pplot"), type = "character", default = file.path(test_data_dir, "results/10X_17_029/percent_mt~10/resolution~0.8/trajectory/plot_monocle-partition-plot.pdf"),
              help = "Partition plot output", metavar = "character")
)

# Parse options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Create necessary directories if they don't exist
dir.create(opt$output.dir, recursive = TRUE, showWarnings = FALSE)

# Remove pre-existing test files
if (file.exists(opt$pplot)) file.remove(opt$pplot)

# Run the command
cmd <- paste(
  "Rscript", file.path(cellsnake_path, "scrna/workflow/scripts/scrna-monocle3.R"),
  "--rds", shQuote(opt$rds),
  "--output.dir", shQuote(opt$output.dir),
  "--pplot", shQuote(opt$pplot)
)

start_time <- Sys.time()
system_output <- system(cmd, intern = TRUE, ignore.stderr = FALSE)
end_time <- Sys.time()

# Calculate execution time
time_diff <- difftime(end_time, start_time, units = "secs")
cat("Time taken: ", time_diff, " seconds\n")

# Print system output
cat("System Output:\n")
cat(system_output, sep = "\n")

# Function to check if the expected output file is created
check_output_with_retry <- function(filename, max_attempts = 6, interval = 10) {
  attempts <- 0
  while (!file.exists(filename)) {
    attempts <- attempts + 1
    if (attempts > max_attempts) {
      stop("FAIL: ", filename, " was not created after ", max_attempts * interval, " seconds.")
    }
    Sys.sleep(interval)
  }
  message("PASS: ", filename, " was created.")
}

# Verify output files
check_output_with_retry(opt$pplot)

message("Test completed successfully.")