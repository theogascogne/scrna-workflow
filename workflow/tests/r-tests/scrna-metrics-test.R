#!/usr/bin/env Rscript

library(optparse)
require(plotly)
require(ggpubr)
require(Seurat)
require(tidyverse)
library(testthat)

# Check if pdflatex is available by running --version
tryCatch({
  system("pdflatex --version", intern = TRUE)
  cat("pdflatex is installed.\n")
}, error = function(e) {
  cat("pdflatex is NOT installed or not found.\n")
})


# Get the path to the 'cellsnake' package
cellsnake_path <- system("python -c 'import cellsnake; print(cellsnake.__path__[0])'", intern = TRUE)

# Define the base path to the testData folder
test_data_dir <- file.path(cellsnake_path, "scrna/workflow/tests/testData")

# Define options with dynamically constructed paths
option_list <- list(
  make_option(c("--rds"), type = "character", default = file.path(test_data_dir, "processed/defaultTest/output.rds"),
              help = "Processed rds file of a Seurat object", metavar = "character"),
  make_option(c("--sampleid"), type = "character", default = "TestSample", help = "Sample ID for testing", metavar = "character"),
  make_option(c("--idents"), type = "character", default = "seurat_clusters", help = "Meta data column name for marker analysis", metavar = "character"),
  make_option(c("--ccplot"), type = "character", default = file.path(test_data_dir, "test/results/ccplot_test.pdf"), help = "Cell cluster count plot", metavar = "character"),
  make_option(c("--ccbarplot"), type = "character", default = file.path(test_data_dir, "test/results/ccbarplot_test.pdf"), help = "Cell cluster bar plot", metavar = "character"),
  make_option(c("--htmlplot"), type = "character", default = file.path(test_data_dir, "test/results/htmlplot_test.pdf"), help = "Cell cluster html plot", metavar = "character"),
  make_option(c("--xlsx"), type = "character", default = file.path(test_data_dir, "test/results/metrics_test.xlsx"), help = "Metrics table output", metavar = "character")
)

# Parse command-line options
opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

try(
  {
    source(file.path(cellsnake_path,"scrna/workflow/scripts/scrna_workflow_helper_functions/scrna-metrics-functions.R"))
  },
  silent = TRUE
)

# Ensure all output directories exist
output_paths <- c(opt$ccplot, opt$ccbarplot, opt$htmlplot, opt$xlsx)
output_dirs <- unique(dirname(output_paths))

for (dir in output_dirs) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  }
}

scrna <<- readRDS(file = opt$rds)

test_that("Summary ggtexttable is generated and saved", {
  tbl_plot <- generate_summary_table(scrna, opt$idents)
  
  id <- length(unique(scrna@meta.data$orig.ident))
  tmpfile <- tempfile(fileext = ".pdf")
  ggsave(tmpfile, tbl_plot, height = 7 + (id * 0.2))
  
  expect_true(file.exists(tmpfile))
  file.remove(tmpfile)
  gl_id <<- id
})

test_that("Cluster summary is written to Excel and TSV", {
  gl_df <<- generate_cluster_summary(scrna, opt$idents)
  openxlsx::write.xlsx(gl_df, opt$xlsx)
  write_tsv(gl_df, str_replace(opt$xlsx, ".xlsx", ".tsv"))
  
  expect_true(file.exists(opt$xlsx))
  expect_true(file.exists(str_replace(opt$xlsx, ".xlsx", ".tsv")))
  
  file.remove(opt$xlsx, str_replace(opt$xlsx, ".xlsx", ".tsv"))
})

test_that("Barplot is created and saved", {
  id <- gl_id
  df <- gl_df
  n <- length(unique(scrna@meta.data[[opt$idents]]))
  p2 <<- save_cluster_barplot(df, id, n, opt$ccbarplot)
  ggsave(opt$ccbarplot, p2, height = 5.2 + (id * 0.09), width = 6 + (n * 0.23))
  
  expect_true(file.exists(opt$ccbarplot))
  file.remove(opt$ccbarplot)
  
  p2 <<- p2  # Pass plot to global for next step
})

# pdfLatex missing in standalone run
#test_that("Barplot is saved as interactive HTML", {
#  ggplotly(p2) %>% htmlwidgets::saveWidget(file = opt$htmlplot, selfcontained = TRUE)
#  expect_true(file.exists(opt$htmlplot))
#  file.remove(opt$htmlplot)
#})

# Define file path
output_file <- file.path(test_data_dir, "results/scrna_metrics_test.txt")
# Define content to write
lines_to_write <- c(
  "metrics ran"
)
# Write the content to the file
writeLines(lines_to_write, con = output_file)

