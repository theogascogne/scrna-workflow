plot_variable_features <- function(scrna, top_features, hvfplot_path, heplot_path) {
  plot1 <- VariableFeaturePlot(scrna)
  plot2 <- LabelPoints(plot = plot1, points = top_features, repel = TRUE)
  ggsave(hvfplot_path, plot2, width = 8, height = 9)

  DimHeatmap(scrna, dims = 1:15, cells = 500, balanced = TRUE, fast = FALSE)
  ggsave(heplot_path, width = 8, height = 15)
}

plot_jackstraw_and_elbow <- function(scrna, jeplot_path) {
  scrna <- JackStraw(scrna, num.replicate = 20, dims = 50, verbose = FALSE)
  scrna <- ScoreJackStraw(scrna, dims = 1:50)
  p1 <- JackStrawPlot(scrna, dims = 1:50)
  p2 <- ElbowPlot(scrna, ndims = 50)
  ggsave(jeplot_path, p1 + p2, width = 13, height = 5)
  return(scrna)
}


handle_normalization_or_integration <- function(scrna, opt) {
  if (isFALSE(opt$integration)) {
    scrna <- NormalizeData(scrna,
                           normalization.method = opt$normalization.method,
                           scale.factor = opt$scale.factor)
    scrna <- FindVariableFeatures(scrna,
                                  selection.method = opt$variable.selection.method,
                                  nfeatures = opt$nfeatures)
  } else {
    DefaultAssay(scrna) <- "integrated"
  }
  scrna
}

perform_scaling_and_pca <- function(scrna) {
  feats <- VariableFeatures(scrna)
  scrna <- ScaleData(scrna, features = feats)
  scrna <- RunPCA(scrna, features = feats)
  list(scrna = scrna, features = feats)
}

get_resolution_range <- function(integration) {
  if (isFALSE(integration)) seq(0.1, 2.5, 0.1) else seq(0.1, 1.5, 0.1)
}

