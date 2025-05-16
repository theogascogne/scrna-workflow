# Define the function for SingleR prediction
run_singleR_prediction <- function(rds_file, output_file = "pred.rds", reference = "HumanPrimaryCellAtlasData", granulation = "label.main") {
  
  # Load the Seurat object
  scrna <- readRDS(file = rds_file)
  DefaultAssay(scrna) <- "RNA"
  
  # Get reference data for SingleR
  ref <- get(reference)()
  
  # Convert Seurat object to SingleCellExperiment
  smObjSCE <- as.SingleCellExperiment(scrna)
  
  # Run SingleR for celltype annotation
  pred <- SingleR(test = smObjSCE, ref = ref, labels = ref[[granulation]])
  
  # Save the prediction to the output file
  saveRDS(pred, file = output_file)
}
