#!/usr/bin/env Rscript


option_list <- list(
      optparse::make_option(c("--rds"),
            type = "character", default = "~/Documents/cellsnake_shared/fetal-brain/analyses/processed/percent_mt~10/resolution~0.8/10X_17_028.rds",
            help = "Processed rds file of a Seurat object", metavar = "character"
      ),
      optparse::make_option(c("--reduction.type"),
            type = "character", default = "umap",
            help = "Reduction type, umap or tsne", metavar = "character"
      ),
      optparse::make_option(c("--pdfplot"),
            type = "character", default = "reduction.pdf",
            help = "Plot file name", metavar = "character"
      ),
      optparse::make_option(c("--htmlplot"),
            type = "character", default = "reduction.html",
            help = "Plot file name", metavar = "character"
      ),
      optparse::make_option(c("--csv"),
            type = "character", default = NULL,
            help = "CSV meta file or this can be a singler RDS file", metavar = "character"
      ),
      optparse::make_option(c("--idents"),
            type = "character", default = "seurat_clusters",
            help = "Meta data column name for marker analysis", metavar = "character"
      ),
      optparse::make_option(c("--percentage"),
            type = "double", default = 5,
            help = "Cluster minimum percentage to plot", metavar = "double"
      ),
      optparse::make_option(c("--labels"), action = "store_true", default = FALSE, help = "Print labels on the plot")
)



opt_parser <- optparse::OptionParser(option_list = option_list)
opt <- optparse::parse_args(opt_parser)

if (is.null(opt$rds)) {
      optparse::print_help(opt_parser)
      stop("At least one argument must be supplied (rds file)", call. = FALSE)
}

require(patchwork)
require(plotly)
require(Seurat)
require(tidyverse)

tryCatch(
      {
            source("~/miniconda3/envs/test/lib/python3.9/site-packages/cellsnake/scrna/workflow/scripts/scrna-functions.R")
      },
      error = function(cond) {
            source(paste0(system("python -c 'import os; import cellsnake; print(os.path.dirname(cellsnake.__file__))'", intern = TRUE), "/scrna/workflow/scripts/scrna-functions.R"))
      }
)

# Hardcoding the GSEA file path#######################
#gsea_file <- "~/miniconda3/envs/test/lib/python3.9/site-packages/cellsnake/scrna/workflow/bundle/c2.cgp.v2022.1.Hs.symbols.gmt"

# Check if the GSEA file exists
#if (!file.exists(gsea_file)) {
#  stop("GSEA file not found: ", gsea_file)
#}

# Optionally, you can load the GSEA file here if required for the analysis
# Example of how to load it (assuming you have a suitable function to handle GMT files):
#gsea_data <- read.gmt(gsea_file)  # You need an appropriate function like `read.gmt` to load the file
###############################################################

scrna <- readRDS(file = opt$rds)
DefaultAssay(scrna) <- "RNA"


if (!is.null(opt$csv)) {
      try({
            metadata <- read.csv(
                  opt$csv,
                  row.names = 1
            )

            scrna@meta.data <- scrna@meta.data %>%
                  tibble::rownames_to_column("barcodes") %>%
                  dplyr::left_join(metadata %>% as.data.frame() %>% rownames_to_column("barcodes"), by = "barcodes") %>%
                  tibble::column_to_rownames("barcodes")
      })
      try({
            pred <- readRDS(opt$csv)
            AddMetaData(scrna, pred["pruned.labels"] %>% as.data.frame() %>% dplyr::select(singler = pruned.labels)) -> scrna
      })
}


# Print unique Idents in the Seurat object
if ("seurat_clusters" %in% colnames(scrna@meta.data)) {
  Idents(scrna) <- scrna@meta.data[["seurat_clusters"]]
  unique_idents <- unique(Idents(scrna))
  cat("Unique Idents in the RDS file:\n")
  print(unique_idents)
} else {
  cat("SEURAT CLUSTERS NOT FOUND IN RDS.\n")
  cat("Columns in meta.data:\n")
  print(colnames(scrna@meta.data))
}

###################### dimplot######################
Idents(object = scrna) <- scrna@meta.data[[opt$idents]]

n <- length(Idents(scrna) %>% unique())

palette <- function_color_palette(n)
palette <- setNames(palette, Idents(scrna) %>% unique())

breaks <- scrna@meta.data %>%
      dplyr::select(opt$idents) %>%
      dplyr::count(get(opt$idents)) %>%
      dplyr::mutate(perc = (n * 100) / sum(n)) %>%
      dplyr::filter(perc >= opt$percentage) %>%
      dplyr::select(-n, -perc) %>%
      distinct() %>%
      pull() %>%
      as.character()


p1 <- DimPlot(scrna, reduction = opt$reduction.type, raster = FALSE) & scale_color_manual(values = palette)

ggplotly(p1) -> p1_plotly

p1_plotly %>% htmlwidgets::saveWidget(file = opt$htmlplot, selfcontained = T)

m <- max(str_count(breaks))

w <- c(8 + (m * 0.09) * (floor(length(breaks) / 11) + 1))


p1 <- DimPlot(scrna, reduction = opt$reduction.type, label = opt$labels, repel = TRUE, raster = FALSE) &
      scale_color_manual(values = palette, breaks = breaks) &
      theme(legend.direction = "horizontal", legend.text = element_text(size = 7)) &
      guides(colour = guide_legend(ncol = 3, override.aes = list(size = 7)))

(p1 / guide_area()) + plot_layout(heights = c(2.5, 1), widths = c(1, 0.6), guides = "collect") -> p1

ggsave(plot = p1, filename = opt$pdfplot, width = 7.5, height = 8)