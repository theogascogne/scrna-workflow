library(covr)
library(optparse)

# Get the path to the 'cellsnake' package
cellsnake_path <- system("python -c 'import cellsnake; print(cellsnake.__path__[0])'", intern = TRUE)

# Define the base path to the testData folder
test_data_dir <- file.path(cellsnake_path, "scrna/workflow/tests/testData")


# Define all source and test files in order
source_files <- c(
  file.path(cellsnake_path,"scrna/workflow/scripts/scrna_workflow_helper_functions/scrna-read-qc-functions.R"),
  file.path(cellsnake_path,"scrna/workflow/scripts/scrna_workflow_helper_functions/scrna-normalization-pca-functions.R"),
  file.path(cellsnake_path,"scrna/workflow/scripts/scrna_workflow_helper_functions/scrna-clusteringtree-functions.R"),
  file.path(cellsnake_path,"scrna/workflow/scripts/scrna_workflow_helper_functions/scrna-dimplot-functions.R"),
  file.path(cellsnake_path,"scrna/workflow/scripts/scrna_workflow_helper_functions/scrna-singler-annotation-functions.R"),
  file.path(cellsnake_path,"scrna/workflow/scripts/scrna_workflow_helper_functions/scrna-metrics-functions.R"),
  file.path(cellsnake_path,"scrna/workflow/scripts/scrna_workflow_helper_functions/scrna-technicals-functions.R"),
  file.path(cellsnake_path,"scrna/workflow/scripts/scrna_workflow_helper_functions/scrna-convert-to-h5ad-functions.R"),
  file.path(cellsnake_path,"scrna/workflow/scripts/scrna_workflow_helper_functions/scrna-find-markers-functions.R"),
  file.path(cellsnake_path,"scrna/workflow/scripts/scrna_workflow_helper_functions/scrna-go-analysis-functions.R"),
  file.path(cellsnake_path,"scrna/workflow/scripts/scrna_workflow_helper_functions/scrna-kegg-functions.R"),
  file.path(cellsnake_path,"scrna/workflow/scripts/scrna_workflow_helper_functions/scrna-marker-heatmap-functions.R"),
  file.path(cellsnake_path,"scrna/workflow/scripts/scrna_workflow_helper_functions/scrna-marker-tables-functions.R"),
  file.path(cellsnake_path,"scrna/workflow/scripts/scrna_workflow_helper_functions/scrna-monocle3-functions.R"),
  file.path(cellsnake_path,"scrna/workflow/scripts/scrna_workflow_helper_functions/scrna-top-marker-plot-functions.R"),
  file.path(cellsnake_path,"scrna/workflow/scripts/scrna_workflow_helper_functions/scrna-singler-plots-functions.R")
)

test_files <- c(
  file.path(cellsnake_path, "scrna/workflow/tests/r-tests/scrna-read-qc-test.R"),
  file.path(cellsnake_path, "scrna/workflow/tests/r-tests/scrna-normalization-pca-test.R"),
  file.path(cellsnake_path, "scrna/workflow/tests/r-tests/scrna-clusteringtree-test.R"),
  file.path(cellsnake_path, "scrna/workflow/tests/r-tests/scrna-dimplot-test.R"),
  file.path(cellsnake_path, "scrna/workflow/tests/r-tests/scrna-singler-annotation-test.R"),
  file.path(cellsnake_path, "scrna/workflow/tests/r-tests/scrna-metrics-test.R"),
  file.path(cellsnake_path, "scrna/workflow/tests/r-tests/scrna-technicals-test.R"),
  file.path(cellsnake_path, "scrna/workflow/tests/r-tests/scrna-convert-to-h5ad-test.R"),
  file.path(cellsnake_path, "scrna/workflow/tests/r-tests/scrna-find-markers-test.R"),
  file.path(cellsnake_path, "scrna/workflow/tests/r-tests/scrna-go_analysis-test.R"),
  file.path(cellsnake_path, "scrna/workflow/tests/r-tests/scrna-kegg-test.R"),
  file.path(cellsnake_path, "scrna/workflow/tests/r-tests/scrna-marker-heatmap-test.R"),
  file.path(cellsnake_path, "scrna/workflow/tests/r-tests/scrna-marker-tables-test.R"),
  file.path(cellsnake_path, "scrna/workflow/tests/r-tests/scrna-monocle3-test.R"),
  file.path(cellsnake_path, "scrna/workflow/tests/r-tests/scrna-top-marker-plot-test.R"),
  file.path(cellsnake_path, "scrna/workflow/tests/r-tests/scrna-singler-plots-test.R")
)

#singler plots missang
# Compute multi-file coverage
coverage <- file_coverage(source_files = source_files, test_files = test_files)

# Print and summarize
print(coverage)
zero_coverage(coverage)
report(coverage)

# Generate HTML report
#report(coverage, file = "combined_coverage_report.html")

csv_output_path <- file.path(cellsnake_path, "scrna/workflow/tests/testData/test_coverage.csv")
write.csv(as.data.frame(coverage), file = csv_output_path, row.names = FALSE)

#coverage <- file_coverage("~/Desktop/scrna_workflow_functions/scrna-technicals-functions.R", "~/miniconda3/envs/test/lib/python3.9/site-packages/cellsnake/scrna/workflow/tests/r-tests/scrna-technicals-test.R")
#print(coverage)
#zero_coverage(coverage)
#report(coverage)

