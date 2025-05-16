import os
from pathlib import Path
from snakemake.io import expand
import sys

# =============================================================
# create_plot_functions.py
# 
# This script is called at the end of the Snakefile execution when 'sample_parameter' is used
# to generate a comprehensive list of file paths for snakemake's rule all which snakemake will create a DAG from
#
# Core purpose:
# -------------
# - Receives main workflow argument from the CLI wrapper via 'sample_parameter'.
# - Dynamically constructs output file paths for various plot types, tables, and analyses
#   based on workflow options and input sample metadata.
#
# Key functions:
# --------------
# - generate_plots: Builds lists of file paths for different plot types and formats, handling
#   optional gene-specific plots, identities, and file extensions.
# - identity_dependent_dimplot / selected_gene_plot / enrichment_analysis / etc.: Wrapper
#   functions to group related plot types and generate corresponding file paths.
# - sample_parameter: The main entry point that selects which groups of plots/files to generate
#   based on the workflow mode option ('minimal', 'standard', 'advanced', etc.).
#
# Important notes:
# ----------------
# - All configuration parameters are sourced from a centralized config object imported from
#   'workflow_utils.Defaults', specifically the 'defaults.py'.
# =============================================================

# os.path.expanduser needs to be used to handle home directory path
workflow_dir = Path(__file__).resolve().parents[1]
sys.path.append(str(workflow_dir))

# Import the function from the file
from workflow_utils.workflow_funcs import file_capture_try
from workflow_utils.Defaults import config
cellsnake_path=config.cellsnake_path #if called by cellsnake

#GSEA
gsea_file = os.path.expanduser(config.gsea_file)
gsea_group=config.gsea_group

files = file_capture_try(config.test_datafolder)

def generate_plots(plot_type: str, paramspace, identity, files, folder_prefix: str, ext: str = None, gene_to_plot=None) -> list:
    """
    Generate different types of plots (e.g., UMAP, TSNE) based on identity.
    If `ext` is provided, adds it as the file extension. If `ext` is None, returns directory paths.
    If gene_to_plot is provided, include it in the path. Skips adding identity if it's empty or not provided,
    and removes any extra slashes that would result from an empty identity.
    """
    plot_patterns = []
    
    # Set identity to a list containing an empty string if it's empty or not provided
    identity = identity or [""]

    if gene_to_plot:  # Only use this logic if genes are specified
        for gene in gene_to_plot:
            plot_patterns += [
                f"{folder_prefix}/{s}/{{params}}/{plot_type}{'/' + i if i else ''}/{gene}/" if ext is None 
                else f"{folder_prefix}/{s}/{{params}}/{plot_type}{'/' + i if i else ''}/{gene}.{ext}"
                for i in identity for s in files
            ]
    else:
        # For paths without gene_to_plot, ensure the correct directory output.
        plot_patterns = [
            f"{folder_prefix}/{s}/{{params}}/{plot_type}{'/' + i if i else ''}/" if ext is None 
            else f"{folder_prefix}/{s}/{{params}}/{plot_type}{'-' + i if i else ''}.{ext}"
            for i in identity for s in files
        ]

    # Expand paths using paramspace
    return expand(plot_patterns, params=list(paramspace.instance_patterns))

def identity_dependent_dimplot(paramspace, identity, files) -> list:
    plots = []
    if len(identity) > 0:
        plot_types = ["plot_dimplot_umap", "plot_dimplot_tsne", "plot_dimplot_pca"]
        for plot_type in plot_types:
            plots += generate_plots(plot_type, paramspace, identity, files, config.results_folder, ext="pdf")
    return plots

def selected_gene_plot(paramspace, gene_to_plot, identity, files) -> list:
    plots = []
    if gene_to_plot:
        plot_types = ["selected_gene_plots_umap", "selected_gene_plots_tsne"]
        for plot_type in plot_types:
            plots += generate_plots(plot_type, paramspace, identity, files, config.results_folder, ext="pdf", gene_to_plot=gene_to_plot)
    return plots

def enrichment_analysis(paramspace, identity, files) -> list:
    plot_types = [
        "enrichment_analysis/table_GO-enrichment", "enrichment_analysis/table_GO-geneset_enrichment",
        "enrichment_analysis/table_KEGG-enrichment", "enrichment_analysis/table_KEGG-geneset_enrichment",
        "enrichment_analysis/table_KEGG-module_enrichment", "enrichment_analysis/table_KEGG-module_geneset_enrichment"
    ]
    outs = []
    if len(set(identity)) > 0:
        for plot_type in plot_types:
            outs += generate_plots(plot_type, paramspace, identity, files, config.results_folder, ext="xlsx")
    return outs

def identity_dependent_metrics(paramspace, identity, files) -> list:
    plots = []
    
    if len(identity) > 0:
        plots += generate_plots("metrics/plot_cellcount", paramspace, identity, files, config.results_folder, ext="pdf")
        plots += generate_plots("metrics/plot_cellcount_barplot", paramspace, identity, files, config.results_folder, ext="pdf")
        plots += generate_plots("metrics/plot_cellcount_barplot", paramspace, identity, files, config.results_folder, ext="html")
        plots += generate_plots("metrics/table_metrics", paramspace, identity, files, config.results_folder, ext="xlsx")
    return plots
    
def gsea_analysis(paramspace, identity, files) -> list:
    outs = []
    print(paramspace)
    print(identity)
    print(files)
    print(os.path.isfile(gsea_file))
    if len(set(identity)) > 0 and os.path.isfile(gsea_file):
        outs = generate_plots("gsea/table_gsea_output", paramspace, identity, files, config.results_folder, ext="xlsx")
    return outs

def celltypist_analysis(paramspace, celltypist_model, identity, files) -> list:
    plots = []
    if celltypist_model:
        base_path = f"celltypist/{celltypist_model}/"
        plot_types = [f"{base_path}plot_celltypist_{plot_type}" for plot_type in ["dotplot", "umap", "tsne"]]
        for plot_type in plot_types:
            plots += generate_plots(plot_type, paramspace, identity, files, config.results_folder, ext="pdf")
    return plots

def metadata_pairwise_deseq_analysis(paramspace, metadata_column, metadata, files) -> list:
    if os.path.isfile(metadata) and config.is_integrated_sample:
        return generate_plots(f"metaplot_volcano-{metadata_column}", paramspace, [metadata_column], files, config.results_folder)
    return []


def kraken_predictions(paramspace, taxa, identity, files) -> list:
    outs = []
    identity.append("singler")  # Add "singler" to identity for analysis

    # Base path for microbiome plots and tables
    base_path = f"microbiome/{config.confidence}_{config.min_hit_groups}"
    
    if config.kraken_db_folder is not None or (config.is_integrated_sample and os.path.isfile(f"analyses_integrated/seurat/{config.confidence}_{config.min_hit_groups}/{config.integration_id}-{config.taxa}.rds")):
        # Generate paths for microbiome dimensionality reduction plots (UMAP and tSNE)
        outs += generate_plots(f"{base_path}/plot_microbiome_dimplot-{taxa}", paramspace, ["umap", "tsne"], files, config.results_folder, ext="pdf")
        
        # Generate paths for microbiome tables
        outs += generate_plots(f"{base_path}/table_microbiome-{taxa}", paramspace, identity, files, config.results_folder, ext="xlsx")
        
        # Add RDS file path directly if kraken_db_folder is not None
        if config.kraken_db_folder is not None:
            outs.append(f"analyses_integrated/seurat/{config.confidence}_{config.min_hit_groups}/{config.integration_id}-{taxa}.rds")
    return outs

def dim_reduction_and_marker_plots(paramspace, identity, files) -> list:
    plots = []
    if config.umap_markers_plot:
        plots += generate_plots("positive_marker_plots_umap", paramspace, identity, files, config.results_folder, ext=None)
    if config.tsne_markers_plot:
        plots += generate_plots("positive_marker_plots_tsne", paramspace, identity, files, config.results_folder, ext=None)
    return plots

def cellchat_plot(paramspace, identity, files) -> list:
    plots = []
    if len(identity) > 0:
        plots = generate_plots("cellchat", paramspace, identity, files, config.results_folder)
    return plots
    
def singler_plots(paramspace, identity, files) -> list:
    plots = []
    if len(identity) > 0:
        plots += generate_plots("singler/plot_score_heatmap", paramspace, identity, files, config.results_folder, ext="pdf")
        plots += generate_plots("singler/plot_clusters", paramspace, identity, files, config.results_folder, ext="pdf")
        plots += generate_plots("singler/plot_score_heatmap_top", paramspace, identity, files, config.results_folder, ext="pdf")
    return plots

def metadata_pairwise_deseq_analysis(paramspace, metadata_column, files) -> list:
    outs = []
    if os.path.isfile(config.metadata) and config.is_integrated_sample:
        plot_types = [f"metaplot_volcano-{metadata_column}", f"metatable_positive-markers-{metadata_column}"]
        for plot_type in plot_types:
            ext = "pdf" if "volcano" in plot_type else "xlsx"
            outs += generate_plots(plot_type, paramspace, [], files, config.results_folder, ext=ext)
    return outs

def sample_parameter(paramspace, files, option) -> list:
    outs = []
    print(paramspace)

    if option in ["standard", "advanced"]:
        outs += generate_plots("summarized_markers-for", paramspace, config.identity_to_analysis, files, config.results_folder, ext="pdf")
        outs += generate_plots("plot_marker-heatmap", paramspace, config.identity_to_analysis, files, config.results_folder, ext="pdf")
        outs += generate_plots("table_average-expression", paramspace, config.identity_to_analysis, files, config.results_folder, ext="xlsx")
        outs += celltypist_analysis(paramspace, config.celltypist_model, config.identity_to_analysis, files)
        outs += enrichment_analysis(paramspace, config.identity_to_analysis, files)
        outs += generate_plots("trajectory/plot_monocle-partition-plot", paramspace, [], files, config.results_folder, ext="pdf")

    if option in ["minimal", "standard", "advanced", "test"]:
        outs += identity_dependent_dimplot(paramspace, config.identity_to_analysis, files)
        outs += identity_dependent_metrics(paramspace, [x for x in config.identity_to_analysis if x != "orig.ident"], files)
        outs += kraken_predictions(paramspace, config.taxa, [x for x in config.identity_to_analysis if x != "orig.ident"], files)
        outs += selected_gene_plot(paramspace, config.gene_to_plot, config.identity_to_analysis, files)
        outs += generate_plots("plot_annotation_tsne", paramspace, [""], files, config.results_folder, ext="pdf")
        outs += generate_plots("plot_annotation_umap", paramspace, [""], files, config.results_folder, ext="pdf")
        outs += generate_plots("plot_annotation_pca", paramspace, [""], files, config.results_folder, ext="pdf")
        outs += metadata_pairwise_deseq_analysis(paramspace, config.metadata_column, files)
        outs += singler_plots(paramspace, config.identity_to_analysis, files)

    if option in ["clustree", "clusteringTree", "minimal", "standard", "advanced", "test"]:
        outs += generate_plots("technicals/plot_clustree", paramspace, [""], files, config.results_folder, ext="pdf")
        outs += generate_plots("technicals/plot_nFeature", paramspace, [""], files, config.results_folder, ext="pdf")
        outs += generate_plots("technicals/plot_nCount", paramspace, [""], files, config.results_folder, ext="pdf")
        outs += generate_plots("technicals/plot_mt.percent", paramspace, [""], files, config.results_folder, ext="pdf")
        outs += generate_plots("technicals/plot_rp.percent", paramspace, [""], files, config.results_folder, ext="pdf")

    if option == "advanced":
        outs += cellchat_plot(paramspace, [x for x in config.identity_to_analysis if x != "orig.ident"], files)
        outs += dim_reduction_and_marker_plots(paramspace, config.identity_to_analysis, files)
    

    if option not in ["clustree", "clusteringTree", "minimal", "standard", "advanced", "test"]:
        print("Please select a correct option...")

    return outs
