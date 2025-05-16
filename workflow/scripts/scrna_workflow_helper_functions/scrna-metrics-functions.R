# Function to generate summary table
generate_summary_table <- function(scrna, idents_col) {
  scrna@meta.data %>%
    dplyr::add_count(orig.ident) %>%
    dplyr::mutate(total_clusters = length(unique(get(idents_col)))) %>%
    distinct(orig.ident, n, total_clusters) %>%
    dplyr::select("Sample Name" = orig.ident, "Total Cells" = n, "Total Clusters" = total_clusters) %>%
    ggpubr::ggtexttable(rows = NULL, theme = ttheme("light"))
}

# Function to generate cluster summary
generate_cluster_summary <- function(scrna, idents_col) {
  scrna@meta.data %>%
    dplyr::group_by(orig.ident, get(idents_col)) %>%
    dplyr::count() %>%
    dplyr::select("Sample Name" = 1, "Cluster" = 2, "Total Cells" = 3)
}

# Function to save cluster bar plot
save_cluster_barplot <- function(df, id, n_clusters, filename) {
  p <- df %>%
    ggplot(aes(x = Cluster, y = `Total Cells`, fill = `Sample Name`)) +
    geom_col() +
    ggthemes::theme_hc() +
    theme(legend.title = element_blank()) +
    guides(colour = guide_legend(ncol = 3, override.aes = list(size = 7)))
  
  ggsave(filename, p, height = 5.2 + (id * 0.09), width = 6 + (n_clusters * 0.23))
  return(p)
}

