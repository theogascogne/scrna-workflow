# Define the function for generating feature plots
generate_feature_plots <- function(rds_file, fplot = NULL, cplot = NULL, mtplot = NULL, rpplot = NULL) {
  # Load the Seurat object
  scrna <- readRDS(file = rds_file)
  
  # Plot and save nFeature plot
  if (!is.null(fplot)) {
    FeaturePlot(scrna, features = "nFeature_RNA", pt.size = 0.1, raster = FALSE)
    ggsave(fplot, width = 7, height = 5)
  }
  
  # Plot and save nCount plot
  if (!is.null(cplot)) {
    FeaturePlot(scrna, features = "nCount_RNA", pt.size = 0.1, raster = FALSE)
    ggsave(cplot, width = 7, height = 5)
  }
  
  # Plot and save Percent MT plot
  if (!is.null(mtplot)) {
    FeaturePlot(scrna, features = "percent.mt", pt.size = 0.1, raster = FALSE)
    ggsave(mtplot, width = 7, height = 5)
  }
  
  # Plot and save Ribo plot
  if (!is.null(rpplot)) {
    FeaturePlot(scrna, features = "percent.rp", pt.size = 0.1, raster = FALSE)
    ggsave(rpplot, width = 7, height = 5)
  }
}