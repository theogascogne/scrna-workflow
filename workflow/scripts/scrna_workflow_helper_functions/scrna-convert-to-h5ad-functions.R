# Helper Function: Convert Seurat assay to v3
convert_assay <- function(seurat_obj, version = "V3") {
  return(scCustomize::Convert_Assay(seurat_object = seurat_obj, convert_to = version))
}

# Helper Function: Convert metadata factors to characters
convert_factors_to_characters <- function(seurat_obj) {
  seurat_obj@meta.data <- seurat_obj@meta.data %>%
    dplyr::mutate(dplyr::across(where(is.factor), as.character))
  return(seurat_obj)
}

# Helper Function: Save Seurat object as h5Seurat
save_h5seurat <- function(seurat_obj, output_path) {
  SaveH5Seurat(seurat_obj, filename = output_path, overwrite = TRUE)
  return(output_path)
}

# Helper Function: Convert .h5Seurat to .h5ad format
convert_to_h5ad <- function(h5seurat_path, output_path) {
  SeuratDisk::Convert(h5seurat_path, dest = "h5ad", overwrite = TRUE)
  file.rename(sub("\\.h5Seurat$", ".h5ad", h5seurat_path), output_path)
  return(output_path)
}
