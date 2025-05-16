# Helper functions

merge_metadata_csv <- function(scrna, csv_path) {
  metadata <- read.csv(csv_path, row.names = 1)
  scrna@meta.data <- scrna@meta.data %>%
    tibble::rownames_to_column("barcodes") %>%
    dplyr::left_join(metadata %>% as.data.frame() %>% tibble::rownames_to_column("barcodes"), by = "barcodes") %>%
    tibble::column_to_rownames("barcodes")
  scrna
}

merge_metadata_rds <- function(scrna, rds_path) {
  pred <- readRDS(rds_path)
  AddMetaData(scrna, pred["pruned.labels"] %>%
                as.data.frame() %>%
                dplyr::select(singler = pruned.labels)) -> scrna
  scrna
}

assign_idents_and_palette <- function(scrna, idents_col) {
  Idents(scrna) <- scrna@meta.data[[idents_col]]
  ids <- unique(Idents(scrna))
  pal <- function_color_palette(length(ids))
  setNames(pal, ids)
}

get_valid_cluster_ids <- function(scrna, idents_col, percentage_threshold) {
  scrna@meta.data %>%
    dplyr::count(!!sym(idents_col)) %>%
    dplyr::mutate(perc = (n * 100) / sum(n)) %>%
    dplyr::filter(perc >= percentage_threshold) %>%
    dplyr::pull(!!sym(idents_col)) %>%
    as.character()
}

generate_dimplot_html <- function(scrna, reduction, palette, output_path) {
  p <- DimPlot(scrna, reduction = reduction, raster = FALSE) & 
    scale_color_manual(values = palette)
  p_html <- ggplotly(p)
  htmlwidgets::saveWidget(p_html, file = output_path, selfcontained = TRUE)
}

generate_dimplot_pdf <- function(scrna, reduction, palette, valid_clusters, labels, output_path) {
  p1 <- DimPlot(scrna, reduction = reduction, label = labels, repel = TRUE, raster = FALSE) & 
    scale_color_manual(values = palette, breaks = valid_clusters) & 
    theme(legend.direction = "horizontal", legend.text = element_text(size = 7)) & 
    guides(colour = guide_legend(ncol = 3, override.aes = list(size = 7)))
  
  (p1 / guide_area()) + plot_layout(heights = c(2.5, 1), widths = c(1, 0.6), guides = "collect") -> p1
  ggsave(plot = p1, filename = output_path, width = 7.5, height = 8)
}
