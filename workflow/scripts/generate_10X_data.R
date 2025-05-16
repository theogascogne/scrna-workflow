generate_data_for_scrna <- function(
    num_genes = 20000,  
    num_cells = 3000,   
    sparsity = 0.95,    
    target_rna_per_cell = 600,  
    deviation_percentage = 1,  
    max_deviation = 200,  
    outlier_fraction = 0.01,  
    mt_sparsity = 0.05,  
    rp_sparsity = 0.05,  
    include_mt_rp = TRUE,      # Control MT/RP usage
    uniform_counts = FALSE,    # Force all cells to same RNA count, no outliers/deviations
    feature_path = "~/Documents/cellsnake_shared/fetal-brain/data/10X_17_029/outs/filtered_feature_bc_matrix/features.tsv.gz",
    barcode_path = "~/Documents/cellsnake_shared/fetal-brain/data/10X_17_029/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"
) {
  library(Matrix)
  library(data.table)
  
  features <- fread(feature_path, header = FALSE)
  colnames(features) <- c("EnsemblID", "GeneSymbol", "FeatureType")
  barcodes <- fread(barcode_path, header = FALSE)
  colnames(barcodes) <- c("CellBarcode")
  
  if (nrow(barcodes) < num_cells) {
    stop("Not enough cell barcodes in the input file!")
  }
  
  mt_genes <- features[grep("^[Mm][Tt]-", features$GeneSymbol), GeneSymbol]
  rp_genes <- features[grep("(?i)^RP[SL]", features$GeneSymbol, perl = TRUE), GeneSymbol]
  remaining_genes <- setdiff(features$GeneSymbol, c(mt_genes, rp_genes))
  
  if (include_mt_rp) {
    if (length(mt_genes) + length(rp_genes) > num_genes) {
      stop("Too many MT and RP genes!")
    }
    gene_names <- c(mt_genes, rp_genes, sample(remaining_genes, num_genes - length(mt_genes) - length(rp_genes)))
  } else {
    gene_names <- sample(remaining_genes, num_genes)
  }
  
  set.seed(42)
  
  # Step 1: RNA counts per cell
  if (uniform_counts) {
    rna_counts <- rep(target_rna_per_cell, num_cells)  # Uniform count for all cells
  } else {
    rna_counts <- rep(target_rna_per_cell, num_cells)
    
    num_outliers <- round(num_cells * outlier_fraction)
    outlier_cells <- sample(1:num_cells, num_outliers, replace = FALSE)
    rna_counts[outlier_cells] <- sample(600:1000, num_outliers, replace = TRUE)
    
    num_deviating_cells <- round(num_cells * deviation_percentage)
    deviating_cells <- sample(1:num_cells, num_deviating_cells, replace = FALSE)
    deviation <- rnorm(num_deviating_cells, mean = 0, sd = max_deviation)
    rna_counts[deviating_cells] <- rna_counts[deviating_cells] + deviation
    rna_counts <- pmax(200, round(rna_counts))
  }
  
  # Step 2: MT/RP logic if enabled
  if (include_mt_rp) {
    mt_rna_counts <- rep(0, num_cells)
    rp_rna_counts <- rep(0, num_cells)
    
    if (uniform_counts) {
      # Proportional fixed sparsity
      mt_rna_counts <- round(rna_counts * mt_sparsity)
      rp_rna_counts <- round(rna_counts * rp_sparsity)
    } else {
      # Include variability for outliers
      for (cell in outlier_cells) {
        mt_percentage <- rbeta(1, 2, 5) * 0.45 + 0.05
        rp_percentage <- rbeta(1, 2, 5) * 0.45 + 0.05
        mt_rna_counts[cell] <- round(mt_percentage * rna_counts[cell])
        rp_rna_counts[cell] <- round(rp_percentage * rna_counts[cell])
      }
      
      non_outlier_cells <- setdiff(1:num_cells, outlier_cells)
      mt_rna_counts[non_outlier_cells] <- rbinom(length(non_outlier_cells), rna_counts[non_outlier_cells], prob = mt_sparsity)
      rp_rna_counts[non_outlier_cells] <- rbinom(length(non_outlier_cells), rna_counts[non_outlier_cells], prob = rp_sparsity)
    }
    
    mt_rna_counts <- pmax(mt_rna_counts, 1)
    rp_rna_counts <- pmax(rp_rna_counts, 1)
    
    mt_gene_per_cell <- sample(mt_genes, num_cells, replace = TRUE)
    rp_gene_per_cell <- sample(rp_genes, num_cells, replace = TRUE)
  }
  
  # Step 3: Create matrix entries
  i <- unlist(lapply(1:num_cells, function(cell) {
    if (include_mt_rp) {
      mt_gene_index <- which(gene_names == mt_gene_per_cell[cell])
      rp_gene_index <- which(gene_names == rp_gene_per_cell[cell])
      mt_entries <- rep(mt_gene_index, mt_rna_counts[cell])
      rp_entries <- rep(rp_gene_index, rp_rna_counts[cell])
      remaining_entries <- sample(setdiff(1:num_genes, c(mt_gene_index, rp_gene_index)), 
                                  rna_counts[cell] - mt_rna_counts[cell] - rp_rna_counts[cell], 
                                  replace = FALSE)
      c(mt_entries, rp_entries, remaining_entries)
    } else {
      sample(1:num_genes, rna_counts[cell], replace = FALSE)
    }
  }))
  
  j <- rep(1:num_cells, rna_counts)
  x <- sample(1:100, length(i), replace = TRUE)
  cell_names <- barcodes$CellBarcode[1:num_cells]
  
  fake_scrna <- sparseMatrix(i = i, j = j, x = x, dims = c(num_genes, num_cells))
  rownames(fake_scrna) <- gene_names
  colnames(fake_scrna) <- cell_names
  
  return(list(
    mat = as(fake_scrna, "dgCMatrix"),
    gene_names = gene_names,
    cell_barcodes = cell_names
  ))
}


inject_mt_rp_counts <- function(mat, mt_genes, rp_genes, min_percent = 3, max_percent = 8) {
  library(Matrix)
  
  mt_rows <- which(rownames(mat) %in% mt_genes)
  rp_rows <- which(rownames(mat) %in% rp_genes)
  all_rows <- rownames(mat)
  
  total_counts <- Matrix::colSums(mat)
  mt_counts <- Matrix::colSums(mat[mt_rows, , drop = FALSE])
  rp_counts <- Matrix::colSums(mat[rp_rows, , drop = FALSE])
  current_percent <- (mt_counts + rp_counts) / total_counts * 100
  
  for (cell_idx in which(current_percent < min_percent)) {
    total <- total_counts[cell_idx]
    target_percent <- runif(1, min_percent, max_percent)
    target_umis <- ceiling((target_percent / 100) * total)
    inject_umis <- target_umis - (mt_counts[cell_idx] + rp_counts[cell_idx])
    if (inject_umis <= 0) next
    
    mt_targets <- sample(mt_genes, sample(1:3, 1), replace = TRUE)
    rp_targets <- sample(rp_genes, sample(1:3, 1), replace = TRUE)
    all_targets <- c(mt_targets, rp_targets)
    counts_split <- rmultinom(1, inject_umis, prob = rep(1, length(all_targets)))
    
    for (j in seq_along(all_targets)) {
      gene <- all_targets[j]
      row_idx <- which(all_rows == gene)
      
      if (length(row_idx) == 0) {
        mat <- rbind(mat, Matrix(0, nrow = 1, ncol = ncol(mat), sparse = TRUE))
        rownames(mat)[nrow(mat)] <- gene
        row_idx <- nrow(mat)
      }
      
      mat[row_idx, cell_idx] <- mat[row_idx, cell_idx] + counts_split[j]
    }
  }
  
  return(mat)
}


subset_and_inject_mt_rp_real_matrix <- function(
    matrix_path = "~/Documents/cellsnake_shared/fetal-brain/data/10X_17_029/outs/filtered_feature_bc_matrix/matrix.mtx.gz",
    barcode_path = "~/Documents/cellsnake_shared/fetal-brain/data/10X_17_029/outs/filtered_feature_bc_matrix/barcodes.tsv.gz",
    feature_path = "~/Documents/cellsnake_shared/fetal-brain/data/10X_17_029/outs/filtered_feature_bc_matrix/features.tsv.gz",
    num_genes = 6000,
    num_cells = 900,
    min_percent = 3,
    max_percent = 8,
    include_mt_rp = TRUE
) {
  library(Matrix)
  library(data.table)
  
  # Load files
  mat <- readMM(matrix_path)
  barcodes <- fread(barcode_path, header = FALSE)[[1]]
  features <- fread(feature_path, header = FALSE)
  colnames(features) <- c("EnsemblID", "GeneSymbol", "FeatureType")
  
  rownames(mat) <- features$GeneSymbol
  colnames(mat) <- barcodes
  
  if (ncol(mat) < num_cells) stop("Not enough cells in the matrix.")
  if (nrow(mat) < num_genes) stop("Not enough genes in the matrix.")
  
  set.seed(42)
  
  # Gene groups
  mt_genes <- features[grep("^[Mm][Tt]-", features$GeneSymbol), GeneSymbol]
  rp_genes <- features[grep("(?i)^RP[SL]", features$GeneSymbol, perl = TRUE), GeneSymbol]
  remaining_genes <- setdiff(features$GeneSymbol, c(mt_genes, rp_genes))
  
  if (include_mt_rp) {
    if (length(mt_genes) + length(rp_genes) > num_genes) {
      stop("Too many MT and RP genes for requested num_genes.")
    }
    gene_names <- c(mt_genes, rp_genes, sample(remaining_genes, num_genes - length(mt_genes) - length(rp_genes)))
  } else {
    gene_names <- sample(remaining_genes, num_genes)
  }
  
  selected_cells <- sample(colnames(mat), num_cells)
  mat_sub <- mat[gene_names, selected_cells, drop = FALSE]
  
  if (include_mt_rp) {
    mat_sub <- inject_mt_rp_counts(mat_sub, mt_genes, rp_genes, min_percent, max_percent)
  }
  
  return(list(
    mat = as(mat_sub, "dgCMatrix"),
    gene_names = rownames(mat_sub),
    cell_barcodes = colnames(mat_sub)
  ))
}



write_10x_files_from_matrix <- function(data, outdir) {
  if (!requireNamespace("Matrix", quietly = TRUE)) stop("Matrix package required")
  if (!requireNamespace("data.table", quietly = TRUE)) stop("data.table package required")
  if (!requireNamespace("R.utils", quietly = TRUE)) stop("R.utils package required")
  
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  mat <- data$mat
  gene_names <- data$gene_names
  cell_barcodes <- data$cell_barcodes
  
  # === matrix.mtx.gz ===
  matrix_path <- file.path(outdir, "matrix.mtx")
  Matrix::writeMM(mat, matrix_path)
  
  # Inject metadata header like CellRanger
  lines <- readLines(matrix_path)
  metadata_line <- '%metadata_json: {"software_version": "cellranger-7.0.0", "format_version": 2}'
  lines <- append(lines, metadata_line, after = 1)
  writeLines(lines, matrix_path)
  R.utils::gzip(matrix_path, destname = paste0(matrix_path, ".gz"), overwrite = TRUE)
  file.remove(matrix_path)
  
  # === features.tsv.gz ===
  features <- data.table::data.table(
    V1 = gene_names,           # Real Gene Symbols as Feature Names
    V2 = gene_names,           # No fake Ensembl IDs, using gene names directly
    V3 = rep("Gene Expression", length(gene_names))
  )
  features_path <- file.path(outdir, "features.tsv")
  data.table::fwrite(features, file = features_path, sep = "\t", col.names = FALSE)
  R.utils::gzip(features_path, destname = paste0(features_path, ".gz"), overwrite = TRUE)
  file.remove(features_path)
  
  # === barcodes.tsv.gz ===
  barcodes_path <- file.path(outdir, "barcodes.tsv")
  writeLines(cell_barcodes, barcodes_path)
  R.utils::gzip(barcodes_path, destname = paste0(barcodes_path, ".gz"), overwrite = TRUE)
  file.remove(barcodes_path)
  
  message("10X files written to ", normalizePath(outdir))
}