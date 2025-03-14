#!/usr/bin/env Rscript

library(optparse)

# Get the path to the 'cellsnake' package
cellsnake_path <- system("python -c 'import cellsnake; print(cellsnake.__path__[0])'", intern = TRUE)

# Define the base path to the testData folder
test_data_dir <- file.path(cellsnake_path, "scrna/workflow/tests/testData")

# Define test options with dynamically constructed paths
option_list <- list(
  make_option(c("--rds"), type = "character", default = file.path(test_data_dir, "fetal-brain/analyses/processed/percent_mt~10/resolution~0.8/10X_17_028.rds"), help = "Processed rds file of a Seurat object", metavar = "character"),
  make_option(c("--sheplot"), type = "character", default = file.path(test_data_dir, "fetal-brain/results/10X_17_028/percent_mt~10/resolution~0.8/singler/plot_score_heatmap-seurat_clusters.pdf"), help = "Output score heatmap plot file name", metavar = "character"),
  make_option(c("--sheplottop"), type = "character", default = file.path(test_data_dir, "fetal-brain/results/10X_17_028/percent_mt~10/resolution~0.8/singler/plot_score_heatmap_top-seurat_clusters.pdf"), help = "Output score heatmap plot file name, top 20", metavar = "character"),
  make_option(c("--pheplot"), type = "character", default = file.path(test_data_dir, "fetal-brain/results/10X_17_028/percent_mt~10/resolution~0.8/singler/plot_clusters-seurat_clusters.pdf"), help = "Output heatmap plot file name", metavar = "character"),
  make_option(c("--idents"), type = "character", default = "seurat_clusters", help = "Meta data column name", metavar = "character"),
  make_option(c("--csv"), type = "character", default = NULL, help = "A meta data table", metavar = "character"),
  make_option(c("--prediction"), type = "character", default = file.path(test_data_dir, "fetal-brain/analyses/singler/percent_mt~10/resolution~0.8/10X_17_028_annotation.rds"), help = "Input prediction file", metavar = "character"),
  make_option(c("--xlsx"), type = "character", default = file.path(test_data_dir, "fetal-brain/results/10X_17_028/percent_mt~10/resolution~0.8/singler/table_annotations_per-seurat_clusters.xlsx"), help = "Input prediction file", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Ensure the directory for output files exists
cat("Creating directories if they do not exist...\n")
dir.create(dirname(opt$sheplot), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(opt$sheplottop), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(opt$pheplot), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(opt$xlsx), recursive = TRUE, showWarnings = FALSE)

# Construct the command to run the original script
cmd <- paste(
  "Rscript", file.path(cellsnake_path, "scrna/workflow/scripts/scrna-singler-plots.R"),
  "--rds", shQuote(opt$rds),
  "--sheplot", shQuote(opt$sheplot),
  "--sheplottop", shQuote(opt$sheplottop),
  "--pheplot", shQuote(opt$pheplot),
  "--idents", shQuote(opt$idents),
  "--prediction", shQuote(opt$prediction),
  "--xlsx", shQuote(opt$xlsx)
)

# Run the original script
cat("Running command:\n", cmd, "\n")
status <- system(cmd, intern = TRUE, ignore.stdout = FALSE, wait = TRUE)

# Capture the exit status (returns 0 if successful, non-zero if failure)
exit_status <- attr(status, "status")
cat("Command output:\n", status, "\n")

# Define color codes for terminal output
cyan <- "\033[0;36m"
red <- "\033[0;31m"
reset <- "\033[0m"

# Function to check if a file exists and is non-empty, with retry logic
check_file_with_retry <- function(file_path, max_attempts = 5, wait_time = 5) {
  attempts <- 0
  while (!file.exists(file_path) || file.info(file_path)$size == 0) {
    attempts <- attempts + 1
    if (attempts >= max_attempts) {
      if (!file.exists(file_path)) {
        cat(red, "FAIL: File not found after", max_attempts, "attempts:", file_path, reset, "\n", sep = " ")
        return(FALSE)
      } else if (file.info(file_path)$size == 0) {
        cat(red, "FAIL: File exists but is empty after", max_attempts, "attempts:", file_path, reset, "\n", sep = " ")
        return(FALSE)
      }
    }
    cat("Waiting for file:", file_path, "Attempt", attempts, "...\n")
    Sys.sleep(wait_time)  # Wait before retrying
  }
  cat(cyan, "PASS: File found and has content:", file_path, reset, "\n", sep = " ")
  return(TRUE)
}

# Check if the output files are created and valid
check_file_with_retry(opt$sheplot)
check_file_with_retry(opt$sheplottop)
check_file_with_retry(opt$pheplot)
check_file_with_retry(opt$xlsx)

# Check if the prediction file exists and contains valid RDS data
if (check_file_with_retry(opt$prediction)) {
  pred <- tryCatch({
    readRDS(opt$prediction)
  }, error = function(e) {
    NULL
  })
  
  if (!is.null(pred)) {
    cat(cyan, "PASS: Prediction file found and successfully loaded:", opt$prediction, reset, "\n", sep = " ")
  } else {
    cat(red, "FAIL: Prediction file found but could not be loaded or is empty:", opt$prediction, reset, "\n", sep = " ")
  }
}