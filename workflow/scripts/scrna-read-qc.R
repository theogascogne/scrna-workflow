#!/usr/bin/env Rscript

# Load required packages
library(optparse)
library(tidyverse)
library(Seurat)
library(patchwork)
library(tools)
library(data.table)

# Get the path to the 'cellsnake' package
cellsnake_path <- system("python -c 'import cellsnake; print(cellsnake.__path__[0])'", intern = TRUE)

# Define the base path to the testData folder
test_data_dir <- file.path(cellsnake_path, "scrna/workflow/tests/testData")

# Define command line options
option_list <- list(
  make_option(c("--min.cells"), type = "integer", default = 3),
  make_option(c("--min.features"), type = "integer", default = 200),
  make_option(c("--max.features"), type = "integer", default = Inf),
  make_option(c("--max.molecules"), type = "integer", default = Inf),
  make_option(c("--min.molecules"), type = "integer", default = 0),
  make_option(c("--data.dir"), type = "character", default = file.path(cellsnake_path, "scrna/workflow/tests/testData/data/testSample/outs/filtered_feature_bc_matrix")),
  make_option(c("--sampleid"), type = "character", default = "defaultRun"),
  make_option(c("--percent.mt"), type = "character", default = "10"),
  make_option(c("--percent.rp"), type = "double", default = 0),
  make_option(c("--before.violin.plot"), type = "character", default = "Cellsnake_DefaultRun/results/before.violin.pdf"),
  make_option(c("--after.violin.plot"), type = "character", default = "Cellsnake_DefaultRun/results/after.violin.pdf"),
  make_option(c("--output.rds"), type = "character", default = "Cellsnake_DefaultRun/raw/output.rds"),
  make_option(c("--plot.mtplot"), type = "character", default = "Cellsnake_DefaultRun/results/plot.mtplot.pdf")
)

if (!exists("opt")) {
  opt_parser <- OptionParser(option_list = option_list)
  opt <- parse_args(opt_parser)
}

try(
  {
    source(file.path(cellsnake_path,"scrna/workflow/scripts/scrna_workflow_helper_functions/scrna-read-qc-functions.R"))
  },
  silent = TRUE
)

# Ensure required arguments are present
if (is.null(opt$data.dir) || is.null(opt$sampleid)) {
  print_help(opt_parser)
  stop("Must supply --data.dir and --sampleid.", call. = FALSE)
}

# Handle optparse bug with Inf
if (is.na(opt$max.features)) opt$max.features <- Inf
if (is.na(opt$max.molecules)) opt$max.molecules <- Inf

# --- MAIN PIPELINE ---

# read input data
scrna.data <- function_read_input(opt)

# from dgmartix made from previus step, create the seurat object
scrna <- CreateSeuratObject(scrna.data, project = make.names(opt$sampleid), min.cells = opt$min.cells, min.features = opt$min.features)
rm(scrna.data)

message("Preparing Seurat object...")
scrna <- prepare_seurat_object_for_qc(scrna, opt)

message("Plotting pre-filtering violin plot...")
ensure_dir_exists(opt$before.violin.plot)
plot_qc_violin(scrna, opt$before.violin.plot)

message("Applying feature/molecule filters...")
scrna <- filter_cells_by_qc_thresholds(scrna, opt)

# Mitochondrial QC
if (tolower(opt$percent.mt) == "auto") {
  message("Running miQC auto-thresholding...")
  ensure_dir_exists(opt$plot.mtplot)
  scrna <- run_miQC_or_fallback(scrna, opt, opt$plot.mtplot)
} else {
  scrna <- subset(scrna, subset = percent.mt <= as.numeric(opt$percent.mt))
}

# Ribosomal filter
scrna <- subset(scrna, subset = percent.rp >= opt$percent.rp)

# Plotting post-filter QC
ensure_dir_exists(opt$after.violin.plot)
plot_qc_violin(scrna, opt$after.violin.plot)

# Saving output
ensure_dir_exists(opt$output.rds)
saveRDS(scrna, opt$output.rds)
