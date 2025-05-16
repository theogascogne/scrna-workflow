# --- Functionalized parts ---
extract_top_markers <- function(Positive_Features) {
  Positive_Features %>%
    group_by(cluster) %>%
    slice_max(n = 20, order_by = avg_log2FC)
}

get_clusters_from_features <- function(Positive_Features) {
  Positive_Features %>%
    distinct(cluster) %>%
    pull()
}

generate_pdf_plot <- function(df, clusters, output_path) {
  pdf(output_path, width = 6, height = 6)
  
  for (i in clusters) {
    df2 <- df %>% filter(cluster == i)
    maxFC <- max(df2$avg_log2FC) + 1
    
    try({
      p1 <- df2 %>%
        mutate(n = dense_rank(desc(avg_log2FC))) %>%
        ggplot(aes(x = n, y = avg_log2FC, label = gene)) +
        geom_text(angle = 75, size = 4) +
        ggthemes::theme_few() +
        ylim(c(0, maxFC)) +
        coord_cartesian(clip = "off", expand = TRUE) +
        ggtitle("Top markers", subtitle = paste(i, "vs all"))
      
      print(p1)
    })
  }
  
  dev.off()
  return(output_path)
}
