
cluster_with_fallback <- function(cds) {
  tryCatch(
    {
      cds <- cluster_cells(cds)
    },
    error = function(e) {
      message("Default clustering failed. Falling back to 'louvain' method.")
      cds <- cluster_cells(cds, cluster_method = "louvain")
    }
  )
  return(cds)
}


initialize_cds_and_cluster <- function(scrna_obj) {
  cds <- as.cell_data_set(scrna_obj)
  cds <- cluster_with_fallback(cds)
  return(cds)
}

generate_partition_and_singler_plot <- function(cds, out_path) {
  p1 <- plot_cells(cds, color_cells_by = "singler", show_trajectory_graph = FALSE)
  p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
  combined <- wrap_plots(p1, p2)
  ggsave(out_path, combined, width = 11, height = 5.5)
}


# Specific inline function
process_partition <- function(partition_id, scrna, output_dir) {
  sub <- subset(scrna, monocle3_partitions == partition_id)
  cds <- as.cell_data_set(sub)
  
  if (partition_id != "1") {
    cds <- cluster_with_fallback(cds)
  }
  
  out_file <- file.path(output_dir, paste0("plot_monocle-partition-", partition_id, ".pdf"))
  plot_and_save_partition_graph(cds, out_file)
  return(out_file)
}

process_all_partitions <- function(cds, output_dir) {
  seurat_obj <- as.Seurat(cds, assay = NULL)
  partitions <- extract_valid_partitions(seurat_obj)
  print(partitions)
  
  for (partition_id in partitions) {
    process_partition(partition_id, seurat_obj, output_dir)
  }
}

extract_valid_partitions <- function(seurat_obj, min_percent = 5, min_cells = 200) {
  partitions <- seurat_obj@meta.data %>%
    dplyr::count(monocle3_partitions) %>%
    dplyr::mutate(perc = (n * 100) / sum(n)) %>%
    dplyr::filter(perc >= min_percent, n > min_cells) %>%
    pull(monocle3_partitions) %>%
    as.character()
  return(partitions)
}

plot_and_save_partition_graph <- function(cds, out_file) {
  p <- plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE)
  ggsave(out_file, plot = p, width = 6, height = 5.5)
}
