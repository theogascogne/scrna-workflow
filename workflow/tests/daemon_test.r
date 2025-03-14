#!/usr/bin/env Rscript

# Load necessary library
library(optparse)

# Define command-line options
option_list <- list(
  make_option(c("-s", "--string"), type = "character", default = "default_string",
              help = "A string argument"),
  make_option(c("-i", "--integer"), type = "integer", default = 42,
              help = "An integer argument"),
  make_option(c("-f", "--float"), type = "double", default = 3.14,
              help = "A floating point argument"),
  make_option(c("-b", "--boolean"), type = "logical", default = FALSE,
              help = "A boolean flag (TRUE or FALSE)")
)

# Parse the arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Print received arguments and their types
cat("\nReceived arguments:\n")
print(opt)

cat("\nDetected types:\n")
cat("string:", class(opt$string), "\n")
cat("integer:", class(opt$integer), "\n")
cat("float:", class(opt$float), "\n")
cat("boolean:", class(opt$boolean), "\n")

# ---- Seurat Environment Check ----

cat("\nChecking if Seurat is preloaded in the environment...\n")

# Check if Seurat is already loaded
if ("Seurat" %in% loadedNamespaces()) {
  cat("Seurat is already loaded in the session!\n")
} else {
  cat("Seurat is NOT preloaded. Attempting to load it now...\n")
  library(Seurat)
}

# Verify that Seurat functions are available
if ("Seurat" %in% loadedNamespaces() && exists("CreateSeuratObject", where = asNamespace("Seurat"))) {
  cat("Seurat is available and functional. Creating a test Seurat object...\n")
  
  # Create a dummy Seurat object
  test_data <- matrix(rnorm(100), nrow = 10)
  seurat_object <- Seurat::CreateSeuratObject(counts = test_data)
  
  cat("Test Seurat object created successfully!\n")
} else {
  cat("Seurat is NOT functioning correctly.\n")
}

cat("\nScript execution completed.\n")