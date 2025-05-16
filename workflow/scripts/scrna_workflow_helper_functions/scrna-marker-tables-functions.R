write_marker_excel_outputs <- function(markers_df, all_path = NULL, positive_path = NULL) {
  if (!is.null(all_path)) {
    openxlsx::write.xlsx(markers_df, file = all_path)
  }
  
  if (!is.null(positive_path)) {
    positive_markers <- dplyr::filter(markers_df, avg_log2FC > 0)
    openxlsx::write.xlsx(positive_markers, file = positive_path)
  }
}
