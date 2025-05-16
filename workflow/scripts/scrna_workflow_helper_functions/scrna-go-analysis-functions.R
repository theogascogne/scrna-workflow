
function_enrichment_go_singlecell <- function(results, p = 0.05, f = 1.5) {
  print(results %>% distinct(cluster) %>% pull())
  results %>%
    as.data.frame() %>%
    dplyr::filter(p_val_adj < p) %>%
    arrange(desc(avg_log2FC)) %>%
    dplyr::select(gene, avg_log2FC) %>%
    dplyr::mutate(GeneID = mapIds(get(opt$mapping), keys = gene, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")) %>%
    dplyr::filter(!is.na(GeneID), !is.na(avg_log2FC), !duplicated(GeneID)) %>%
    dplyr::select(3, 2) %>%
    deframe() -> geneList
  
  gene <- names(geneList)[abs(geneList) > f]
  
  tryCatch({
    kk <- clusterProfiler::enrichGO(
      gene = gene,
      universe = names(geneList),
      OrgDb = get(opt$mapping),
      ont = "ALL",
      pAdjustMethod = "fdr",
      pvalueCutoff = 1,
      readable = TRUE
    )
  }, error = function(e) {
    kk <- NULL
  }) -> kk
  
  tryCatch({
    kk2 <- clusterProfiler::gseGO(
      geneList = names(geneList),
      OrgDb = get(opt$mapping),
      ont = "ALL",
      minGSSize = 2,
      maxGSSize = 500,
      pvalueCutoff = 1,
      verbose = FALSE,
      eps = 0
    )
  }, error = function(e) {
    kk2 <- NULL
  }) -> kk2
  
  return(list(kk, kk2))
}


save_go_results <- function(results, path) {
  results %>%
    keep(~ !is.null(.[[1]])) %>%
    map(~ .[[1]]@result) %>%
    bind_rows(.id = "cluster") %>%
    as_tibble() %>%
    openxlsx::write.xlsx(path)
}

save_gse_results <- function(results, path) {
  results %>%
    keep(~ !is.null(.[[2]])) %>%
    map(~ .[[2]]@result) %>%
    bind_rows(.id = "cluster") %>%
    as_tibble() %>%
    openxlsx::write.xlsx(path)
}

