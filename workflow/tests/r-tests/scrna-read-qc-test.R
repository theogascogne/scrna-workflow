#!/usr/bin/env Rscript

library(optparse)
installed.packages()["patchwork", ]

# Get the path to the 'cellsnake' package
cellsnake_path <- system("python -c 'import cellsnake; print(cellsnake.__path__[0])'", intern = TRUE)

# Define the base path to the testData folder
test_data_dir <- file.path(cellsnake_path, "scrna/workflow/tests/testData")

# Check if the directory exists, and if not, stop execution
if (!dir.exists(file.path(test_data_dir, "/data"))) {
  stop("Data directory does not exist: ", file.path(test_data_dir, "/data"))
}


# Define options with default values as specified
option_list <- list(
  make_option(c("--data.dir"), type = "character", default = file.path(test_data_dir, "data/10X_17_028/outs/filtered_feature_bc_matrix"), help = "Data directory", metavar = "character"),
  make_option(c("--output.rds"), type = "character", default = file.path(test_data_dir, "results/10X_17_028/output.rds"), help = "Output RDS file", metavar = "character"),
  make_option(c("--sampleid"), type = "character", default = "10X_17_028", help = "Sample ID", metavar = "character"),
  make_option(c("--percent.mt"), type = "character", default = "10", help = "Max mitochondrial gene percentage", metavar = "character"),
  make_option(c("--percent.rp"), type = "double", default = 0, help = "Min ribosomal gene percentage", metavar = "double"),
  make_option(c("--min.features"), type = "integer", default = 200, help = "Min features", metavar = "integer"),
  make_option(c("--max.features"), type = "character", default = "Inf", help = "Max features", metavar = "character"),
  make_option(c("--max.molecules"), type = "character", default = "Inf", help = "Max molecules", metavar = "character"),
  make_option(c("--min.cells"), type = "integer", default = 3, help = "Min cells", metavar = "integer"),
  make_option(c("--before.violin.plot"), type = "character", default = file.path(test_data_dir, "results/10X_17_028/test/technicals/plot_before-qc-trimming.pdf"), help = "Violin plot before QC trimming", metavar = "character"),
  make_option(c("--after.violin.plot"), type = "character", default = file.path(test_data_dir, "results/10X_17_028/test/technicals/plot_after-qc-trimming.pdf"), help = "Violin plot after QC trimming", metavar = "character"),
  make_option(c("--plot.mtplot"), type = "character", default = file.path(test_data_dir, "results/10X_17_028/test/technicals/plot_model-metrics-mitochondrial-genes.pdf"), help = "Plot for mitochondrial gene metrics", metavar = "character"),
  make_option(c("--output.txt"), type = "character", default = file.path(test_data_dir, "results/10X_17_028/output.txt"), help = "Output text file path", metavar = "character")
)

# Parse options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Create necessary directories if they don't exist
dir.create(dirname(opt$before.violin.plot), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(opt$after.violin.plot), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(opt$output.rds), recursive = TRUE, showWarnings = FALSE)
if (opt$plot.mtplot != "") {
  dir.create(dirname(opt$plot.mtplot), recursive = TRUE, showWarnings = FALSE)
}

if (!file.exists(dirname(opt$output.txt))) {
  dir.create(dirname(opt$output.txt), recursive = TRUE, showWarnings = FALSE)
}

# Remove any pre-existing test files
#file.remove(opt$output.txt, opt$before.violin.plot, opt$after.violin.plot, opt$plot.mtplot)

cmd <- paste(
  "Rscript", file.path(cellsnake_path, "scrna/workflow/scripts/scrna-read-qc.R"),  # Use cellsnake_path to create the path dynamically
  "--data.dir", opt$data.dir,
  "--output.rds", opt$output.rds,
  "--sampleid", opt$sampleid,
  "--percent.rp", opt$percent.rp,
  "--percent.mt", opt$percent.mt,
  "--min.features", opt$min.features,
  "--max.features", opt$max.features,
  "--min.cells", opt$min.cells,
  "--before.violin.plot", opt$before.violin.plot,
  "--after.violin.plot", opt$after.violin.plot,
  if (!is.null(opt$max.molecules)) paste("--max.molecules", opt$max.molecules) else "",
  if (opt$plot.mtplot != "") paste("--plot.mtplot", opt$plot.mtplot) else ""
)

start_time <- Sys.time()
# Capture the output and error messages of the system call
system_output <- system(cmd, intern = TRUE, ignore.stderr = FALSE)

end_time <- Sys.time()
# Calculate the time difference
time_diff <- difftime(end_time, start_time, units = "secs")
# Print the time difference
cat("Time taken for the system call: ", time_diff, "seconds\n")

# Print the output directly to the console for more comprehensive error reporting
cat("System Output:\n")
cat(system_output, sep = "\n")

check_output_with_retry <- function(filename, max_attempts = 8, interval = 15) {
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

# Open a text file for logging results
output_log <- file(opt$output.txt, "w")

# Write details to the text file
cat("Sample ID:", opt$sampleid, "\n", file=output_log)
cat("Max mitochondrial gene percentage:", opt$percent.mt, "\n", file=output_log)
cat("Min ribosomal gene percentage:", opt$percent.rp, "\n", file=output_log)
cat("Min features:", opt$min.features, "\n", file=output_log)
cat("Max features:", opt$max.features, "\n", file=output_log)
cat("Min cells:", opt$min.cells, "\n", file=output_log)
cat("Before QC violin plot path:", opt$before.violin.plot, "\n", file=output_log)
cat("After QC violin plot path:", opt$after.violin.plot, "\n", file=output_log)

# Check if plotting was requested
if (opt$plot.mtplot != "") {
  cat("Plot for mitochondrial gene metrics:", opt$plot.mtplot, "\n", file=output_log)
}

cat("Processing completed for sample:", opt$sampleid, "\n", file=output_log)

# Close the file after writing
close(output_log)

# Check output files with retries
check_output_with_retry(opt$output.txt)
check_output_with_retry(opt$output.rds)
check_output_with_retry(opt$before.violin.plot)
check_output_with_retry(opt$after.violin.plot)
#if (opt$plot.mtplot != "") {
#  check_output_with_retry(opt$plot.mtplot)
#}

message("All tests passed successfully.")

