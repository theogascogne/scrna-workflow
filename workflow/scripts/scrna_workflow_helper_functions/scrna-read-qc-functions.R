
ensure_dir_exists <- function(file_path) {
  dir_path <- dirname(file_path)
  if (!dir.exists(dir_path)) {
    message("Creating missing directory: ", dir_path)
    dir.create(dir_path, recursive = TRUE)
  }
}


function_read_input <- function(opt) {
  try(
    {
      scrna.data <- Read10X(data.dir = opt$data.dir)
      return(scrna.data)
    },
    silent = TRUE
  )
  try(
    {
      scrna.data <- Read10X_h5(filename = paste0(opt$data.dir, "/filtered_feature_bc_matrix.h5"))
      return(scrna.data)
    },
    silent = TRUE
  )
  try(
    {
      scrna.data <- Read10X_h5(filename = opt$data.dir)
      return(scrna.data)
    },
    silent = TRUE
  )
  
  try(
    {
      x <- tolower(file_ext(opt$data.dir))
      
      if (x %in% c("h5")) {
        scrna.data <- Read10X_h5(filename = opt$data.dir)
      } else if (x %in% c("gz")) {
        scrna.data <- fread(cmd = paste("gunzip -dc", opt$data.dir))
        
        scrna.data <- scrna.data %>% column_to_rownames("V1")
      } else if (x %in% c("zip")) {
        scrna.data <- fread(cmd = paste("unzip -p", opt$data.dir))
        
        scrna.data <- scrna.data %>% column_to_rownames("V1")
      } else if (x %in% c("csv", "tsv")) {
        scrna.data <- fread(paste(opt$data.dir))
        
        scrna.data <- scrna.data %>% column_to_rownames("V1")
      }
      return(scrna.data)
    },
    silent = TRUE
  )
}



prepare_seurat_object_for_qc <- function(scrna, opt, rename = TRUE, add_qc = TRUE, normalize = TRUE) {
  if (rename) {
    scrna <- RenameCells(object = scrna, add.cell.id = make.names(opt$sampleid))
  }
  if (add_qc) {
    scrna[["percent.mt"]] <- PercentageFeatureSet(scrna, pattern = "^[Mm][Tt]-")
    scrna[["percent.rp"]] <- PercentageFeatureSet(scrna, pattern = "(?i)(^RP[SL])")
  }
  if (normalize) {
    scrna <- NormalizeData(scrna)
  }
  return(scrna)
}

filter_cells_by_qc_thresholds <- function(scrna, opt) {
  subset(scrna, subset = nFeature_RNA < opt$max.features &
           nFeature_RNA >= opt$min.features &
           nCount_RNA < opt$max.molecules &
           nCount_RNA >= opt$min.molecules)
}

run_miQC_or_fallback <- function(scrna, opt, fallback_plot_file) {
  if (isFALSE(all(scrna$percent.mt == 0))) {
    
    require(SingleCellExperiment)
    require(miQC)
    require(scater)
    
    smObjSCE <- as.SingleCellExperiment(scrna)
    mt_genes <- grepl("^[Mm][Tt]-", rownames(smObjSCE))
    feature_ctrls <- list(mito = rownames(smObjSCE)[mt_genes])
    smObjSCE <- addPerCellQC(smObjSCE, subsets = feature_ctrls)
    
    tryCatch({
      model <- mixtureModel(smObjSCE)
      ggsave(fallback_plot_file, plotModel(smObjSCE, model) + plotMetrics(smObjSCE), width = 10, height = 4)
      scrna[, colnames(filterCells(smObjSCE, model))]
    }, error = function(e) {
      upper_bound <- median(scrna$percent.mt) + mad(scrna$percent.mt, constant = 1)
      ggsave(fallback_plot_file, plot.new() + plotMetrics(smObjSCE), width = 10, height = 4)
      subset(scrna, subset = percent.mt <= upper_bound)
    })
  } else {
    scrna
  }
}

plot_qc_violin <- function(scrna, file_path) {
  p <- VlnPlot(scrna, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rp"), ncol = 4)
  ggsave(file_path, plot = p, width = 10, height = 4)
}