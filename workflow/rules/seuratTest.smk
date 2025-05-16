from collections import defaultdict
from yaml import load
import os

test_data = f"{cellsnake_path}workflow/tests/testData/"

test_outputs = [
    "results/clustree_test.txt",
    "processed/defaultTest/output.rds",
    "results/scrna_celltypist.txt",
    "results/scrna_singler_plots_test.txt",
    "results/table_annotations_per-seurat_clusters.xlsx",
    "results/scrna_technicals_test.txt",
    "results/scrna_dimplot_test.txt",
    "results/defaultTest/table_positive-markers-seurat_clusters.xlsx",
    "results/defaultTest/table_all-markers-seurat_clusters.xlsx",
    "results/scrna-top-marker-plot-test.txt",
    "results/scrna-marker-heatmap-test.txt",
    "results/scrna-kegg-test.txt",
    "results/scrna-go-analysis-test.txt",
    "results/scrna_monocle3_test.txt",
    "results/scrna_metrics_test.txt"
]

rule test_scrna_read_qc_r:
    input:
        #"data/{sample}/raw_feature_bc_matrix/"
           
    output:
        f"{cellsnake_path}workflow/tests/testData/raw/defaultTest/output.rds"
    run:
        shell("{cellsnake_path}workflow/tests/r-tests/scrna-read-qc-test.R")


rule test_scrna_clusteringtree_r:
    input:
        f"{cellsnake_path}workflow/tests/testData/raw/defaultTest/output.rds"  
    output:
        f"{cellsnake_path}workflow/tests/testData/results/clustree_test.txt"
    run:
        shell("{cellsnake_path}workflow/tests/r-tests/scrna-clusteringtree-test.R")


rule test_scrna_normalization_pca_r:
    input:
        f"{cellsnake_path}workflow/tests/testData/raw/defaultTest/output.rds"  
    output:
        rds=f"{cellsnake_path}workflow/tests/testData/processed/defaultTest/output.rds"
    run:
        shell("{cellsnake_path}workflow/tests/r-tests/scrna-normalization-pca-test.R")



rule test_scrna_convert_to_5had_r:
    input:
        f"{cellsnake_path}workflow/tests/testData/processed/defaultTest/output.rds"
    output:
        f"{cellsnake_path}workflow/tests/testData/analyses/defaultTest/testSample.h5ad"
    run:
        shell("{cellsnake_path}workflow/tests/r-tests/scrna-convert-to-h5ad-test.R")


rule test_scrna_celltypist_py:
    input:
        f"{cellsnake_path}workflow/tests/testData/analyses/defaultTest/testSample.h5ad"
    output:
        f"{cellsnake_path}workflow/tests/testData/results/scrna_celltypist.txt"
    run:
        shell("pytest {cellsnake_path}workflow/tests/py-tests/scrna-celltypist-test.py")



rule test_scrna_singler_plots_r:
    input:
        f"{cellsnake_path}workflow/tests/testData/processed/defaultTest/output.rds"
    output:
        f"{cellsnake_path}workflow/tests/testData/results/scrna_singler_plots_test.txt"
    run:
        shell("{cellsnake_path}workflow/tests/r-tests/scrna-singler-plots-test.R")


rule test_scrna_singler_annotation_r:
    input:
        f"{cellsnake_path}workflow/tests/testData/processed/defaultTest/output.rds"
    output:
        f"{cellsnake_path}workflow/tests/testData/singler/defaultTest/annotation.rds"
    run:
        shell("{cellsnake_path}workflow/tests/r-tests/scrna-singler-annotation-test.R")

rule plot_singler_celltype:
    input:
        f"{cellsnake_path}workflow/tests/testData/processed/defaultTest/output.rds",
        f"{cellsnake_path}workflow/tests/testData/singler/defaultTest/annotation.rds"
    output:
        f"{cellsnake_path}workflow/tests/testData/results/table_annotations_per-seurat_clusters.xlsx"
    run:
        shell("{cellsnake_path}workflow/tests/r-tests/scrna-singler-plots-test.R")


rule test_scrna_technicals_r:
    input:
        f"{cellsnake_path}workflow/tests/testData/processed/defaultTest/output.rds"
    output:
        f"{cellsnake_path}workflow/tests/testData/results/scrna_technicals_test.txt"
    run:
        shell("{cellsnake_path}workflow/tests/r-tests/scrna-technicals-test.R")



rule test_scrna_dimplots_r:
    input:
        f"{cellsnake_path}workflow/tests/testData/processed/defaultTest/output.rds"
    output:
        f"{cellsnake_path}workflow/tests/testData/results/scrna_dimplot_test.txt"
    run:
        shell("{cellsnake_path}workflow/tests/r-tests/scrna-dimplot-test.R")


rule test_scrna_find_markers_r:
    input:
        f"{cellsnake_path}workflow/tests/testData/processed/defaultTest/output.rds"
    output:
        f"{cellsnake_path}workflow/tests/testData/analyses/defaultTest/seurat_clusters.rds"
    run:
        shell("{cellsnake_path}workflow/tests/r-tests/scrna-find-markers-test.R")


rule test_scrna_marker_tables_r:
    input:
        f"{cellsnake_path}workflow/tests/testData/analyses/defaultTest/seurat_clusters.rds"
    output:
        f"{cellsnake_path}workflow/tests/testData/results/defaultTest/table_positive-markers-seurat_clusters.xlsx",
        f"{cellsnake_path}workflow/tests/testData/results/defaultTest/table_all-markers-seurat_clusters.xlsx"
    run:
        shell("{cellsnake_path}workflow/tests/r-tests/scrna-marker-tables-test.R")


rule test_scrna_top_marker_plot_r:
    input:
        f"{cellsnake_path}workflow/tests/testData/results/defaultTest/table_positive-markers-seurat_clusters.xlsx"
    output:
        f"{cellsnake_path}workflow/tests/testData/results/scrna-top-marker-plot-test.txt"
    run:
        shell("{cellsnake_path}workflow/tests/r-tests/scrna-top-marker-plot-test.R")

rule test_scrna_marker_heat_map_r:
    input:
        f"{cellsnake_path}workflow/tests/testData/processed/defaultTest/output.rds",
        f"{cellsnake_path}workflow/tests/testData/results/defaultTest/table_positive-markers-seurat_clusters.xlsx"
    output:
        f"{cellsnake_path}workflow/tests/testData/results/scrna-marker-heatmap-test.txt"
    run:
        shell("{cellsnake_path}workflow/tests/r-tests/scrna-marker-heatmap-test.R")


rule test_scrna_kegg_r:
    input:
        f"{cellsnake_path}workflow/tests/testData/results/defaultTest/table_all-markers-seurat_clusters.xlsx"
    output:
        f"{cellsnake_path}workflow/tests/testData/results/scrna-kegg-test.txt"
    run:
        shell("{cellsnake_path}workflow/tests/r-tests/scrna-kegg-test.R")
        
        
rule test_scrna_metrics_r:
    input:
        f"{cellsnake_path}workflow/tests/testData/processed/defaultTest/output.rds"
    output:
        f"{cellsnake_path}workflow/tests/testData/results/scrna_metrics_test.txt"
    run:
        shell("{cellsnake_path}workflow/tests/r-tests/scrna-metrics-test.R")
        
rule test_scrna_monocle3_r:
    input:
        f"{cellsnake_path}workflow/tests/testData/processed/defaultTest/output.rds"
    output:
        f"{cellsnake_path}workflow/tests/testData/results/scrna_monocle3_test.txt"
    run:
        shell("{cellsnake_path}workflow/tests/r-tests/scrna-monocle3-test.R")
        

rule test_scrna_go_analysis_r:
    input:
        f"{cellsnake_path}workflow/tests/testData/results/defaultTest/table_all-markers-seurat_clusters.xlsx"
    output:
        f"{cellsnake_path}workflow/tests/testData/results/scrna-go-analysis-test.txt"
    run:
        shell("{cellsnake_path}workflow/tests/r-tests/scrna-go_analysis-test.R")


rule test_complete:
    input:
        expand(test_data + "{file}", file=test_outputs)
    output:
        test_data + "results/test_complete.txt"
    shell:
        "touch {output}"
        

rule covr_coverage:
    input:
        #"data/{sample}/raw_feature_bc_matrix/"
    output:
        f"{cellsnake_path}scrna/workflow/tests/testData/test_coverage.csv"
    run:
        shell("{cellsnake_path}workflow/scripts/covr_coverage.R")
