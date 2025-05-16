# Combined function to plot and save the score heatmap (both regular and top heatmap)
save_score_heatmap <- function(pred, file_path, show_labels = TRUE, max_labels = NULL, width = 15, height = 8) {
  
  # Check if we want a top heatmap (with limited labels)
  if (is.null(max_labels)) {
    # Regular heatmap
    p1 <- plotScoreHeatmap(pred, show.labels = show_labels)
  } else {
    # Top heatmap with limited labels
    p1 <- plotScoreHeatmap(pred, show.labels = show_labels, max.labels = max_labels)
  }
  
  # Save the plot to the specified file path
  ggsave(plot = p1, filename = file_path, width = width, height = height)
}


# Function to generate and save heatmap to a file
save_heatmap_plot <- function(scrna, pred, tab, file_path, width_factor = 0.10, height_factor = 0.10) {
  # Generate the heatmap plot
  p1 <- pheatmap(log2(tab + 10), color = colorRampPalette(c("white", "blue"))(101))
  
  # Calculate dimensions based on the unique values
  n1 <- length(unique(scrna$seurat_clusters))
  n2 <- length(unique(pred$pruned.labels))
  
  # Save the heatmap plot
  ggsave(plot = p1, filename = file_path, width = 6 + (n1 * width_factor), height = 4 + (n2 * height_factor))
}


# Function to save the table to an Excel file
save_table_to_excel <- function(tab, file_path) {
  dir.create(dirname(file_path), recursive = TRUE, showWarnings = FALSE)
  openxlsx::write.xlsx(tab, file = file_path, rowNames = FALSE)
}
