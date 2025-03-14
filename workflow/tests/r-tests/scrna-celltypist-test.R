#!/usr/bin/env Rscript

library(optparse)

# Get the path to the 'cellsnake' package
cellsnake_path <- system("python -c 'import cellsnake; print(cellsnake.__path__[0])'", intern = TRUE)

# Define the base path to the testData folder
test_data_dir <- file.path(cellsnake_path, "scrna/workflow/tests/testData")

# Define options
option_list <- list(
  make_option(c("--rds"), type = "character", 
              default = file.path(test_data_dir, "analyses/processed/percent_mt~10/resolution~0.8/10X_17_029.rds"),
              help = "Seurat object RDS file", metavar = "character"),
  make_option(c("--csv"), type = "character", 
              default = file.path(test_data_dir, "analyses/celltypist/Immune_All_Low.pkl/percent_mt~10/resolution~0.8/10X_17_029/seurat_clusters/predicted_labels.csv"),
              help = "Celltypist prediction CSV", metavar = "character"),
  make_option(c("--output.tsne.plot"), type = "character", 
              default = file.path(test_data_dir, "results/10X_17_029/percent_mt~10/resolution~0.8/celltypist/Immune_All_Low.pkl/plot_celltypist_tsne-seurat_clusters.pdf"),
              help = "Output t-SNE plot PDF", metavar = "character"),
  make_option(c("--output.umap.plot"), type = "character", 
              default = file.path(test_data_dir, "results/10X_17_029/percent_mt~10/resolution~0.8/celltypist/Immune_All_Low.pkl/plot_celltypist_umap-seurat_clusters.pdf"),
              help = "Output UMAP plot PDF", metavar = "character"),
  make_option(c("--percentage"), type = "double", default = 5,
              help = "Cluster minimum percentage to plot", metavar = "double"),
  make_option(c("--labels"), action = "store_true", default = FALSE,
              help = "Print labels on the plot")
)

# Parse options
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Create necessary directories if they don't exist
dir.create(dirname(opt$output.tsne.plot), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(opt$output.umap.plot), recursive = TRUE, showWarnings = FALSE)

# Remove pre-existing test output
if (file.exists(opt$output.tsne.plot)) {
  file.remove(opt$output.tsne.plot)
}
if (file.exists(opt$output.umap.plot)) {
  file.remove(opt$output.umap.plot)
}

# Run the script
cmd <- paste(
  "Rscript", file.path(cellsnake_path, "scrna/workflow/scripts/scrna-celltypist.R"),
  "--rds", opt$rds,
  "--csv", opt$csv,
  "--output.tsne.plot", opt$output.tsne.plot,
  "--output.umap.plot", opt$output.umap.plot,
  "--percentage", opt$percentage,
  "--labels"
)

system_output <- system(cmd, intern = TRUE, ignore.stderr = FALSE)

# Print system output
cat("System Output:\n")
cat(system_output, sep = "\n")

# Validate output files
if (!file.exists(opt$output.tsne.plot)) {
  stop("FAIL: ", opt$output.tsne.plot, " was not created.")
} else {
  message("PASS: ", opt$output.tsne.plot, " was created successfully.")
}

if (!file.exists(opt$output.umap.plot)) {
  stop("FAIL: ", opt$output.umap.plot, " was not created.")
} else {
  message("PASS: ", opt$output.umap.plot, " was created successfully.")
}

message("Test completed successfully.")