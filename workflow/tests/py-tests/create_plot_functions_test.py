import os
import sys
#import yaml    #yaml is handled by Defaults.py
import pytest
import pandas as pd
import re
from pathlib import Path

# Ensure module path is recognized
# work out current path, then go two steps back to workflow directory:  workflow/tests/py-tests, so defaults can be found inside workflow_utils.
workflow_dir = Path(__file__).resolve().parents[2]
sys.path.append(str(workflow_dir))

# Import needed functions and function to test from their respective files
from snakemake.utils import Paramspace
from workflow_utils.create_plot_functions import *
from workflow_utils.workflow_funcs import file_capture_try, initialization_of_paramspace
from workflow_utils.Defaults import config

# Path to the .yaml file
tsv_file_path = config.tsv_file_path
selected_gene_file = config.markers_file_path

resolution = config.resolution
percent_mt = config.percent_mt
option = config.option

identity_to_analysis = config.identity_to_analysis
results_folder = config.results_folder
gene_to_plot = config.gene_to_plot

taxa = config.taxa
celltypist_model = config.celltypist_model
metadata = config.metadata
metadata_column = config.metadata_column

kraken_db_folder = config.kraken_db_folder
cellsnake_path = config.cellsnake_path
gsea_file = config.gsea_file
gsea_group = config.gsea_group

files = file_capture_try(config.test_datafolder)

# Initialize paramspace
par_df = initialization_of_paramspace(tsv_file_path, {"percent_mt": [percent_mt], "resolution": [resolution]})
paramspace = Paramspace(par_df)


try:
    if os.path.isfile(selected_gene_file):
        df = pd.read_table(selected_gene_file)
        gene_to_plot += df[df.columns[0]].values.tolist()
except:
    pass


def assert_plot_paths(plots, extra_check=None):
    """Ensure plots list is valid, non-empty when required, and contains only strings."""
    
    assert isinstance(plots, list), "Output should be a list"

    # Flatten nested lists if needed
    plots = [item for sublist in plots for item in sublist] if any(isinstance(plot, list) for plot in plots) else plots

    # Ensure all elements in plots are strings
    assert all(isinstance(plot, str) for plot in plots), "Each plot should be a string"

    # If extra_check is empty, plots should also be empty
    if extra_check == []:
        assert len(plots) == 0, "If extra_check is empty, plots must also be empty."
        return  # No further checks needed if extra_check is empty and plots are empty.

    # If the plots list is empty, check if extra_check is provided, otherwise it's fine.
    if len(plots) == 0:
        # If extra_check is provided (not None or []), raise an error.
        if extra_check not in [None, []]:
            assert False, "Output list should not be empty when extra_check is provided"
        return  # Skip further checks if empty is fine when extra_check is not provided.
    
    assert all(isinstance(plot, str) for plot in plots), "Each plot should be a string"

# unfinished - needs to check for what if tsne and umap are there or not
@pytest.mark.parametrize("test_paramspace, test_identity_to_analysis, test_files, expected_pattern", [
    (paramspace, identity_to_analysis, files, r"^(?:results\/|results_folder\/)[a-zA-Z0-9_]+\/percent_mt~\d+\/resolution~\d+(\.\d+)?\/positive_marker_plots_umap\/seurat_clusters\/$"),
])
def test_dim_reduction_and_marker_plots(test_paramspace, test_identity_to_analysis, test_files, expected_pattern):
    """Test dim_reduction_and_marker_plots function with multiple scenarios."""
    plots = dim_reduction_and_marker_plots(test_paramspace, test_identity_to_analysis, test_files)

    assert_plot_paths(plots)

    for plot in plots:
        assert re.match(expected_pattern, plot), f"Plot path '{plot}' does not match the expected pattern."


@pytest.mark.parametrize("test_paramspace, test_gene_to_plot, test_identity_to_analysis, expected_pattern", [
    (paramspace, ["CD4"], identity_to_analysis, r"^(results|results_folder)\/[a-zA-Z0-9_]+\/percent_mt~\d+\/resolution~\d+(\.\d+)?\/selected_gene_plots_(umap|tsne)\/seurat_clusters\/CD\d+\.pdf$"),
    (paramspace, [], identity_to_analysis, None),
    (paramspace, ["CD4"], [], r"^(results|results_folder)\/[a-zA-Z0-9_]+\/percent_mt~\d+\/resolution~\d+(\.\d+)?\/selected_gene_plots_(umap|tsne)\/CD\d+\.pdf$")  # Identity is empty, no seurat_clusters in the pattern
])
def test_selected_gene_plot(test_paramspace, test_gene_to_plot, test_identity_to_analysis, expected_pattern):
    """Test selected_gene_plot function with different scenarios."""

    plots = selected_gene_plot(test_paramspace, test_gene_to_plot, test_identity_to_analysis, files)
    #print(plots)
    #assert len(plots) == 10000, f"Intentional fail for debugging"
    assert_plot_paths(plots, test_gene_to_plot)

    if expected_pattern:
        for plot in plots:
            assert re.match(expected_pattern, plot), f"Plot path '{plot}' does not match expected pattern."
            assert any(gene in plot for gene in test_gene_to_plot), f"Expected gene '{test_gene_to_plot}' in plot path."
 

@pytest.mark.parametrize("test_paramspace, test_identity_to_analysis, test_files, expected_pattern", [
    (paramspace, identity_to_analysis, files, r"^(?:results\/|results_folder\/)[a-zA-Z0-9_]+\/percent_mt~\d+\/resolution~\d+(\.\d+)?\/metrics\/(plot_cellcount|plot_cellcount_barplot|table_metrics)-seurat_clusters\.(pdf|html|xlsx)$"),
    (paramspace, [], files, None)  # No identity, expect empty list
])
def test_identity_dependent_metrics(test_paramspace, test_identity_to_analysis, test_files, expected_pattern):
    """Test identity_dependent_metrics function."""
    plots = identity_dependent_metrics(test_paramspace, test_identity_to_analysis, test_files)

    assert_plot_paths(plots)

    if expected_pattern:
        for plot in plots:
            assert re.match(expected_pattern, plot), f"Plot path '{plot}' does not match the expected pattern."
    else:
        assert len(plots) == 0, "The output list should be empty when no identities are provided."


@pytest.mark.parametrize("test_paramspace, test_identity_to_analysis, test_files, expected_pattern", [
    (paramspace, identity_to_analysis, files, r"^(results|results_folder)\/[a-zA-Z0-9_]+\/percent_mt~\d+\/resolution~\d+(\.\d+)?\/cellchat\/seurat_clusters\/$"),
    (paramspace, [], files, None)
])
def test_cellchat_plot(test_paramspace, test_identity_to_analysis, test_files, expected_pattern):
    """Test cellchat_plot function."""
    plots = cellchat_plot(test_paramspace, test_identity_to_analysis, test_files)

    assert_plot_paths(plots)

    if expected_pattern:
        for plot in plots:
            assert re.match(expected_pattern, plot), f"Plot path '{plot}' does not match the expected pattern."
    else:
        assert len(plots) == 0, "The output list should be empty when no identities are provided."


@pytest.mark.parametrize("test_paramspace, test_identity_to_analysis, test_files, expected_pattern", [
    (paramspace, identity_to_analysis, files, r"^(results|results_folder)\/[a-zA-Z0-9_]+\/percent_mt~\d+\/resolution~\d+(\.\d+)?\/singler\/(plot_score_heatmap|plot_clusters|plot_score_heatmap_top)-seurat_clusters\.pdf$"),
    (paramspace, [], files, None)
])
def test_singler_plot(test_paramspace, test_identity_to_analysis, test_files, expected_pattern):
    """Test singler_plot function."""
    plots = singler_plots(test_paramspace, test_identity_to_analysis, test_files)

    assert_plot_paths(plots)

    if expected_pattern:
        for plot in plots:
            assert re.match(expected_pattern, plot), f"Plot path '{plot}' does not match the expected pattern."
    else:
        assert len(plots) == 0, "The output list should be empty when no identities are provided."


# Test function for kraken_predictions - in theory this works and does catch errors but the are some edge cases that needs checking.
def test_kraken_predictions():
    # Ensure identity_to_analysis does not contain 'orig.ident'
    identities = [x for x in identity_to_analysis if x != "orig.ident"]

    # Call the kraken_predictions function
    output = kraken_predictions(paramspace, taxa, identities, files)
    print(output)
    print(paramspace)
    print(taxa)
    print(identities)
    
    #assert len(output) == 10000000, "Debugging: This will always fail in order to show output"

    # Use the helper function to validate common properties of plot paths
    #assert_plot_paths(output, identities, results_folder)

    # If kraken_db_folder is defined, ensure the output list is not empty
    if kraken_db_folder:
        assert len(output) > 0, "Output list should not be empty when kraken_db_folder is defined"

        # Additional checks specific to kraken_predictions
        for path in output:
            # Check that the path contains the taxa (e.g., genus, species) being analyzed
            assert taxa in path, f"The taxa '{taxa}' should be included in the file path"

            # Check that the path contains one of the reduction methods (umap, tsne)
            assert any(method in path for method in ["umap", "tsne"]), \
                "Path should contain one of the reduction methods (umap, tsne)"

        # TODO: Implement regex check for Kraken prediction output paths
        # pattern = r"..."
        # assert re.match(pattern, path), f"File path '{path}' does not match the expected pattern"

    # If no kraken_db_folder is defined, ensure the output list is empty
    else:
        assert len(output) == 0, "Output should be empty when no kraken_db_folder is defined"


@pytest.mark.parametrize("test_paramspace, test_celltypist_model, test_identity_to_analysis, test_files, expected_pattern", [
    # Case: Model is provided -> Should generate valid plots matching pattern
    (paramspace, celltypist_model, identity_to_analysis, files, 
     r"^results\/[a-zA-Z0-9_]+\/percent_mt~\d+\/resolution~\d+(\.\d+)?\/celltypist\/"
     r"[a-zA-Z0-9_]+\.pkl\/plot_celltypist_(dotplot|umap|tsne)-seurat_clusters\.pdf$"),
    
    # Case: No model provided -> Should return an empty list
    (paramspace, None, identity_to_analysis, files, None)
])
def test_celltypist_analysis(test_paramspace, test_celltypist_model, test_identity_to_analysis, test_files, expected_pattern):
    """Test celltypist_analysis function with different scenarios."""
    
    plots = celltypist_analysis(test_paramspace, test_celltypist_model, test_identity_to_analysis, test_files)

    assert_plot_paths(plots, test_celltypist_model)

    if expected_pattern:
        # If a model is provided, ensure plots are generated and match the pattern
        assert len(plots) > 0, "Expected non-empty plots list when celltypist_model is provided"
        for plot in plots:
            assert test_celltypist_model in plot, f"The celltypist model '{test_celltypist_model}' should be included in the file path"
            assert re.match(expected_pattern, plot), f"File path '{plot}' does not match the expected pattern"
    
    else:
        # If no model is provided, the plots list should be empty
        assert len(plots) == 0, "Output list should be empty when no celltypist_model is provided"


# Parameterized test function for gsea_analysis - unfinished, find a way to test for missing gsea_file
@pytest.mark.parametrize("test_paramspace, test_identity_to_analysis, gsea_file_exists, expected_pattern", [
        (paramspace, ["seurat_clusters"], True, r"^results/[a-zA-Z0-9_]+/percent_mt~\d+/resolution~\d+(\.\d+)?/gsea/table_gsea_output-seurat_clusters\.xlsx$"),
        (paramspace, ["seurat_clusters"], False, None),
        (paramspace, [], True, None),
        (paramspace, [], False, None)
    ]
)
def test_gsea_analysis(test_paramspace, test_identity_to_analysis, gsea_file_exists, expected_pattern):
    """Test gsea_analysis function with different scenarios."""
    
    # Mocking the existence of the gsea file if needed - placeholder for now
    if not gsea_file_exists:
        os.path.isfile = lambda x: False
    else:
        os.path.isfile = lambda x: True

    # Call the gsea_analysis function
    outs = gsea_analysis(test_paramspace, test_identity_to_analysis, files)

    assert_plot_paths(outs)

    if expected_pattern:
        # If identity_to_analysis and gsea_file exist, the list should not be empty
        assert len(outs) > 0, "Expected non-empty output list when identity_to_analysis and gsea_file exist"
        
        # Check that the paths in outs match the expected GSEA pattern
        for out in outs:
            assert re.match(expected_pattern, out), f"The path '{out}' does not match the expected GSEA output pattern"
    
    else:
        # If no identities or gsea_file, the list should be empty
        assert len(outs) == 0, "Output list should be empty when no identities or gsea_file is provided"


@pytest.mark.parametrize("test_paramspace, test_identity_to_analysis, test_files, expected_pattern", [
    # Case: Identity provided -> Should generate valid enrichment analysis outputs
    (paramspace, identity_to_analysis, files,
     r"^results/[a-zA-Z0-9_]+/percent_mt~\d+/resolution~\d+(\.\d+)?/enrichment_analysis/"
     r"table_(GO|KEGG)-(enrichment|geneset_enrichment|module_enrichment|module_geneset_enrichment)-seurat_clusters\.xlsx$"),
    (paramspace, [], files, None)
])
def test_enrichment_analysis(test_paramspace, test_identity_to_analysis, test_files, expected_pattern):
    """Test enrichment_analysis function with and without identity_to_analysis."""
    
    outs = enrichment_analysis(test_paramspace, test_identity_to_analysis, test_files)

    assert_plot_paths(outs, test_identity_to_analysis)

    if expected_pattern:
        # If identities are provided, ensure the output list is not empty
        assert len(outs) > 0, "Output list should not be empty when identity_to_analysis is provided"

        # Check each output file against the expected pattern
        for out in outs:
            assert re.match(expected_pattern, out), f"The path '{out}' does not match the expected enrichment output pattern"
    
    else:
        # If no identities are provided, the output list should be empty
        assert len(outs) == 0, "Output list should be empty when no identities are provided."

#### admittedly the last two test are also unfinished.
# integrated sample is needed here and i have not tested for that yet
# the last test requires the error above to be resolved and the integrated sample test to be finished.
def test_metadata_pairwise_deseq_analysis():
    outs = metadata_pairwise_deseq_analysis(paramspace, metadata_column, files)

    # Print the output to see the generated output paths
    print(outs)

    # Use the helper function to validate common properties of the paths
    #assert_plot_paths(outs)

    # If metadata exists and is_integrated_sample is True, the list should not be empty
    if os.path.isfile(metadata) and config.is_integrated_sample:
        assert len(outs) > 0, "Output list should not be empty when metadata exists and is_integrated_sample is True"

        # Placeholder for regex validation of paths
        # TODO: Add regex pattern here to validate metadata pairwise DESeq output paths
        # Example pattern structure might look like:
        # ^results/[a-zA-Z0-9_]+/metadata_column/metaplot_volcano|metatable_positive-markers

        for out in outs:
            # Check that each output is a string
            assert isinstance(out, str), "Each output should be a string"

            # Check if 'results_folder' or 'results' is in the path
            assert results_folder in out or "results" in out, \
                f"Results folder '{results_folder}' or 'results' should be part of the path"

            # Check if the metadata_column is in the path
            assert metadata_column in out, f"Metadata column '{metadata_column}' should be part of the output path"

            # Check if the path contains 'metaplot_volcano' or 'metatable_positive-markers'
            assert any(term in out for term in ["metaplot_volcano", "metatable_positive-markers"]), \
                "Path should contain 'metaplot_volcano' or 'metatable_positive-markers'"

    # If no metadata or not an integrated sample, the list should be empty
    else:
        assert len(outs) == 0, "Output list should be empty when no metadata or not an integrated sample"



      
def test_sample_parameter():
    # Options to test
    options = ["standard", "advanced", "minimal", "clustree", "clusteringTree", "invalid_option"]
    
    # Expected path patterns for specific options - this is worng, needs correction
    standard_pattern = r"results_folder\/[a-zA-Z0-9_]+\/{params}\/(?:summarized_markers-for-|plot_marker-heatmap-|table_average-expression-)[a-zA-Z0-9_]+\.pdf|\.xlsx$"
    advanced_pattern = r"results_folder\/[a-zA-Z0-9_]+\/{params}\/(?:summarized_markers-for-|plot_marker-heatmap-|table_average-expression-)[a-zA-Z0-9_]+\.pdf|\.xlsx$"
    minimal_pattern = r"results_folder\/[a-zA-Z0-9_]+\/{params}\/plot_annotation_(umap|tsne|pca)\.pdf$"
    clustree_pattern = r"results_folder\/[a-zA-Z0-9_]+\/{params}\/technicals\/plot_(clustree|nFeature|nCount|mt\.percent|rp\.percent)\.pdf$"
    
    # Iterate through each option and run tests
    for test_option in options:
        # Run the function under the specified option context
        outs = sample_parameter(paramspace, files, 'standard')
        print(outs)
        assert len(outs) == 10000000, "Debugging: This will always fail in order to show output"

        #commented out untill further notice uwu
        # Additional tests specific to each option (commented out for placeholders)
        # if option == "standard":
        #     for path in outs:
        #         assert re.search(standard_pattern, path), f"Path '{path}' does not match the 'standard' option pattern."

        # elif option == "advanced":
        #     for path in outs:
        #         assert re.search(advanced_pattern, path), f"Path '{path}' does not match the 'advanced' option pattern."
        #     # Ensure specific outputs for 'advanced' option
        #     assert any("cellchat" in path for path in outs), "Expected 'cellchat' output missing in 'advanced' option."
        #     assert any("dim_reduction_and_marker_plots" in path for path in outs), "Expected 'dim reduction' output missing."

        # elif option == "minimal":
        #     for path in outs:
        #         assert re.search(minimal_pattern, path), f"Path '{path}' does not match the 'minimal' option pattern."

        # elif option in ["clustree", "clusteringTree"]:
        #     for path in outs:
        #         assert re.search(clustree_pattern, path), f"Path '{path}' does not match the 'clustree' option pattern."

        # Test an invalid option to confirm error handling
        if option == "invalid_option":
            assert outs is None, "Invalid option should return None."


