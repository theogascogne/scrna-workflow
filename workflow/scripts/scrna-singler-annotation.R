#!/usr/bin/env Rscript

require(optparse)
require(SingleR)
require(SingleCellExperiment)
require(Seurat)

# Get the path to the 'cellsnake' package
cellsnake_path <- system("python -c 'import cellsnake; print(cellsnake.__path__[0])'", intern = TRUE)

# Define the base path to the testData folder
test_data_dir <- file.path(cellsnake_path, "scrna/workflow/tests/testData")

option_list <- list(
    optparse::make_option(c("--rds"),
        type = "character", default = NULL,
        help = "A list of RDS files of Seurat objects", metavar = "character"
    ),
    optparse::make_option(c("--output"),
        type = "character", default = "pred.rds",
        help = "Output prediction file", metavar = "character"
    ),
    optparse::make_option(c("--reference"),
        type = "character", default = "HumanPrimaryCellAtlasData",
        help = "SingleR reference", metavar = "character"
    ),
    optparse::make_option(c("--granulation"),
        type = "character", default = "label.main",
        help = "SingleR granulation level", metavar = "character"
    )
)


if (!exists("opt")) {
  opt_parser <- OptionParser(option_list = option_list)
  opt <- parse_args(opt_parser)
}

if (is.null(opt$rds)) {
    optparse::print_help(opt_parser)
    stop("At least one argument must be supplied (rds and sampleid)", call. = FALSE)
}

try(
  {
    source(file.path(cellsnake_path,"scrna/workflow/scripts/scrna_workflow_helper_functions/scrna-singler-annotation-functions.R"))
  },
  silent = TRUE
)


scrna <- readRDS(file = opt$rds)
DefaultAssay(scrna) <- "RNA"

# celltype annotation with SingleR
ref <- get(opt$reference)()

smObjSCE <- as.SingleCellExperiment(scrna)
pred <- SingleR(test = smObjSCE, ref = ref, labels = ref[[opt$granulation]])

saveRDS(pred, file = opt$output)