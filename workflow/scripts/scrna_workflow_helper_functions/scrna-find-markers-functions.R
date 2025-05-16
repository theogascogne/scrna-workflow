run_and_save_markers <- function(scrna, idents_col, logfc.threshold, test.use, output_path) {
  Idents(scrna) <- scrna@meta.data[[idents_col]]
  markers <- FindAllMarkers(scrna, logfc.threshold = logfc.threshold, test.use = test.use)
  saveRDS(markers, file = output_path)
}

