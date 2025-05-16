normalize_and_select_features <- function(scrna, opt) {
  scrna <- NormalizeData(scrna, normalization.method = opt$normalization.method, scale.factor = opt$scale.factor)
  scrna <- FindVariableFeatures(scrna, selection.method = opt$variable.selection.method, nfeatures = opt$nfeatures)
  return(scrna)
}

run_pca_pipeline <- function(scrna) {
  not.all.genes <- VariableFeatures(scrna)
  scrna <- ScaleData(scrna, features = not.all.genes)
  scrna <- RunPCA(scrna, features = not.all.genes)
  dimensionReduction <- function_pca_dimensions(scrna)
  scrna <- FindNeighbors(scrna, dims = 1:dimensionReduction)
  return(list(scrna = scrna, dims = dimensionReduction))
}

retrieve_clustering <- function(scrna, opt, dimensionReduction) {
  if (!opt$resolution %in% c("auto", "AUTO", "Auto")) {
    scrna <- FindClusters(scrna, resolution = as.numeric(opt$resolution))
  } else {
    if (!requireNamespace("MultiKParallel", quietly = TRUE)) {
      remotes::install_github("sinanugur/MultiKParallel", upgrade = "never")
    }
    require(MultiKParallel)
    
    scrna_tmp <- FindClusters(scrna, resolution = seq(0.2, 2.5, 0.15))
    multik <- MultiKParallel(scrna_tmp, reps = 10, seed = 255, resolution = seq(0.2, 2.5, 0.15), numCores = opt$cpu, nPC = dimensionReduction)
    K <- multik$k %>%
      tibble::as_tibble() %>%
      count(value) %>%
      filter(n == max(n)) %>%
      slice(1) %>%
      pull(value)
    
    optimal_resolution <- scrna_tmp[[]] %>%
      as.data.frame() %>%
      select(starts_with("RNA")) %>%
      rownames_to_column("barcode") %>%
      mutate(across(where(is.factor), as.character)) %>%
      as_tibble() %>%
      gather(res, clu, contains("RNA")) %>%
      distinct(res, clu) %>%
      count(res) %>%
      mutate(res = str_remove(res, "RNA_snn_res.")) %>%
      mutate(diff = abs(n - K)) %>%
      slice_min(order_by = diff) %>%
      pull(res) %>%
      as.numeric()
    
    scrna <- FindClusters(scrna, resolution = optimal_resolution)
  }
  return(scrna)
}

run_TSNE_UMAP <- function(scrna, opt, dims) {
  if (opt$umap) {
    scrna <- RunUMAP(scrna, dims = dims)
  }
  if (opt$tsne) {
    scrna <- RunTSNE(scrna, dims = dims, check_duplicates = FALSE)
  }
  return(scrna)
}

run_doublet_filter <- function(scrna, opt) {
  if (opt$doublet.filter) {
    if (!requireNamespace("DoubletFinder", quietly = TRUE)) {
      remotes::install_github("chris-mcginnis-ucsf/DoubletFinder", upgrade = "never")
    }
    require(DoubletFinder)
    
    homotypic.prop <- modelHomotypic(scrna$seurat_clusters)
    nExp_poi <- round(0.075 * nrow(scrna@meta.data)) ## Assuming 7.5% doublet formation rate - tailor for your dataset
    nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
    
    ## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
    scrna <- doubletFinder(scrna, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
    scrna <- doubletFinder(scrna, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = paste0("pANN_0.25_0.09_", nExp_poi), sct = FALSE)
    
    scrna@meta.data %>%
      tibble::rownames_to_column("barcodes") %>%
      select(barcodes, starts_with("DF")) %>%
      select(barcodes, DoubletFinder = 2) -> doublet_df
    
    scrna@meta.data <- scrna@meta.data %>%
      select(!starts_with("DF")) %>%
      select(!starts_with("pANN")) %>%
      tibble::rownames_to_column("barcodes") %>%
      dplyr::left_join(doublet_df, by = "barcodes") %>%
      tibble::column_to_rownames("barcodes")
    
    subset(scrna, subset = DoubletFinder == "Singlet") -> scrna
    # scrna$DoubletFinder <- NULL
  }
  
  return(scrna)
}

annotate_with_singleR <- function(scrna, reference) {
  ref <- get(reference)()
  DefaultAssay(scrna) <- "RNA"
  pred <- SingleR(test = as.SingleCellExperiment(scrna), ref = ref, labels = ref$label.fine)
  scrna <- AddMetaData(scrna, pred["pruned.labels"] %>% as.data.frame() %>% dplyr::select(singler = pruned.labels))
  attr(scrna, "SingleRref") <- reference
  try({
    DefaultAssay(scrna) <- "integrated"
  })
  return(scrna)
}
