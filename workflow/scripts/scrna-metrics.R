#!/usr/bin/env Rscript

require(plotly)
require(ggpubr)
require(Seurat)
require(tidyverse)
require(optparse)
require(openxlsx)

# Get the path to the 'cellsnake' package
cellsnake_path <- system("python -c 'import cellsnake; print(cellsnake.__path__[0])'", intern = TRUE)

# Define the base path to the testData folder
test_data_dir <- file.path(cellsnake_path, "scrna/workflow/tests/testData")


# Define command-line options
option_list <- list(
  optparse::make_option(c("--rds"),
                        type = "character", default = "~/Documents/cellsnake_shared/fetal-brain/analyses/raw/percent_mt~10/resolution~0.8/10X_17_028.rds",
                        help = "Processed rds file of a Seurat object", metavar = "character"
  ),
  optparse::make_option(c("--sampleid"),
                        type = "character", default = NULL,
                        help = "Sample ID", metavar = "character"
  ),
  optparse::make_option(c("--idents"),
                        type = "character", default = "seurat_clusters",
                        help = "Meta data column name for marker analysis", metavar = "character"
  ),
  optparse::make_option(c("--ccplot"),
                        type = "character", default = "ccplot.pdf",
                        help = "Cell cluster count plot", metavar = "character"
  ),
  optparse::make_option(c("--ccbarplot"),
                        type = "character", default = "ccbarplot.pdf",
                        help = "Cell cluster count plot", metavar = "character"
  ),
  optparse::make_option(c("--htmlplot"),
                        type = "character", default = "htmlplot.pdf",
                        help = "Cell cluster html plot", metavar = "character"
  ),
  optparse::make_option(c("--xlsx"),
                        type = "character", default = "metrics.xlsx",
                        help = "Metrics table output", metavar = "character"
  )
)

# Parse command-line options
if (!exists("opt")) {
  opt_parser <- OptionParser(option_list = option_list)
  opt <- parse_args(opt_parser)
}

# Check if required argument is provided
if (is.null(opt$rds)) {
  optparse::print_help(opt_parser)
  stop("At least one argument must be supplied (rds file and sampleid)", call. = FALSE)
}

try(
  {
    source(file.path(cellsnake_path,"scrna/workflow/scripts/scrna_workflow_helper_functions/scrna-metrics-functions.R"))
  },
  silent = TRUE
)

# Load the Seurat object
scrna <- readRDS(file = opt$rds)

# Step 1: Generate summary table and save as a PDF
summary_table <- generate_summary_table(scrna, opt$idents)
id <- length(unique(scrna@meta.data$orig.ident))
ggsave(opt$ccplot, summary_table, height = 7 + (id * 0.2))

# Step 2: Generate cluster summary and write to Excel and TSV files
cluster_summary <- generate_cluster_summary(scrna, opt$idents)
openxlsx::write.xlsx(cluster_summary, opt$xlsx)
write_tsv(cluster_summary, file = str_replace(opt$xlsx, ".xlsx", ".tsv"))

# Step 3: Generate and save bar plot
n_clusters <- length(unique(scrna@meta.data[[opt$idents]]))
bar_plot <- save_cluster_barplot(cluster_summary, id, n_clusters, opt$ccbarplot)

# Step 4: Save the bar plot as an interactive plot
interactive_plot <- ggplotly(bar_plot)
interactive_plot %>% htmlwidgets::saveWidget(file = opt$htmlplot, selfcontained = TRUE)