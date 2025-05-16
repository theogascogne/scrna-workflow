scale_scrna_data <- function(scrna, xlsx_path, idents_col) {
  Idents(object = scrna) <- scrna@meta.data[[idents_col]]
  
  Positive_Features <- openxlsx::read.xlsx(xlsx_path)
  df <- Positive_Features %>%
    group_by(cluster) %>%
    slice_max(n = 5, order_by = avg_log2FC)
  
  not.all.genes <- df %>%
    distinct(gene) %>%
    pull()
  
  scrna <- suppressWarnings(ScaleData(scrna, features = not.all.genes))
  
  list(scrna = scrna, scaled_genes = not.all.genes)
}

plot_and_save_heatmap <- function(scrna, not_all_genes, output_path) {
  p1 <- DoHeatmap(object = scrna, features = not_all_genes, label = FALSE) & 
    theme(axis.text.y = element_text(size = 5)) & 
    scale_fill_continuous(type = "viridis")
  
  ggsave(
    plot = p1,
    filename = output_path,
    height = 4 + (length(unique(Idents(scrna))) * 0.2),
    width = 5 + (length(unique(Idents(scrna))) * 0.05),
    useDingbats = TRUE
  )
}

write_avg_expression_to_excel <- function(scrna, idents_col, output_path) {
  avg_expr <- AverageExpression(scrna, group.by = idents_col)[["RNA"]] %>%
    as.data.frame() %>%
    rownames_to_column("id")
  
  openxlsx::write.xlsx(avg_expr, file = output_path, rowNames = FALSE)
}

