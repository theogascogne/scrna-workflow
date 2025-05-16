

function_enrichment_kegg_singlecell <- function(results, p = 0.05, f = 1.5, mapping_db, organism_code) {
  cluster_id <- unique(results$cluster)
  print(cluster_id)
  
  valid_syms <- keys(mapping_db, keytype = "SYMBOL")
  results <- results %>%
    dplyr::filter(p_val_adj < p, gene %in% valid_syms)
  
  if (nrow(results) == 0) {
    message("Cluster ", cluster_id, " has no valid gene SYMBOLs. Skipping enrichment.")
    return(list(NULL, NULL, NULL, NULL))
  }
  
  geneList <- results %>%
    arrange(desc(avg_log2FC)) %>%
    dplyr::select(gene, avg_log2FC) %>%
    dplyr::mutate(GeneID = mapIds(mapping_db, keys = gene, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")) %>%
    dplyr::filter(!is.na(GeneID), !is.na(avg_log2FC), !duplicated(GeneID)) %>%
    dplyr::select(GeneID, avg_log2FC) %>%
    deframe()
  
  gene <- names(geneList)[abs(geneList) > f]
  
  kk1 <- tryCatch(
    clusterProfiler::enrichKEGG(
      gene = gene, organism = organism_code,
      pAdjustMethod = "fdr", minGSSize = 2, pvalueCutoff = 1
    ), error = function(e) NULL
  )
  
  kk2 <- tryCatch(
    clusterProfiler::gseKEGG(
      geneList = geneList, organism = organism_code,
      pvalueCutoff = 1, pAdjustMethod = "fdr",
      minGSSize = 2, eps = 0, verbose = FALSE
    ), error = function(e) NULL
  )
  
  kk3 <- tryCatch(
    clusterProfiler::enrichMKEGG(
      gene = gene, organism = organism_code,
      pvalueCutoff = 1, minGSSize = 2, pAdjustMethod = "fdr"
    ), error = function(e) NULL
  )
  
  kk4 <- tryCatch(
    clusterProfiler::gseMKEGG(
      geneList = geneList, organism = organism_code,
      minGSSize = 2, keyType = "kegg", pAdjustMethod = "fdr", pvalueCutoff = 1
    ), error = function(e) NULL
  )
  
  return(list(kk1, kk2, kk3, kk4))
}


save_kegg_results <- function(results_list, idx, filename) {
  # Filter non-null results at position idx and bind rows with cluster id
  data_to_save <- results_list %>%
    purrr::keep(~ !is.null(.[[idx]])) %>%
    purrr::map(~ .[[idx]]@result) %>%
    dplyr::bind_rows(.id = "cluster") %>%
    tibble::as_tibble()
  
  openxlsx::write.xlsx(data_to_save, filename)
}
