import os
import yaml
import random
import datetime

# Defaults.py is used to allow .py scripts access to the default values in the workflow as they were previously declared in the snakefile
# Calls to its values will first use the .yaml file's configurations, if it is not found, it will then use the defaults dictionary found below.

# to avoid needing hard coded paths, the following section should dynamically find cellsnake's 'root directory'
try:
    import cellsnake
    # Get the root directory of the `cellsnake` package
    CELLSNAKE_ROOT = os.path.dirname(cellsnake.__file__)
except ImportError:
    CELLSNAKE_ROOT = None  # Handle missing package gracefully

# Define the fallback strategy for getting the path to the cellsnake package or project
def get_cellsnake_path():
    cellsnake_path = os.getenv("CELLSNAKE_PATH", CELLSNAKE_ROOT)
    if cellsnake_path:
        return cellsnake_path
    
    current_dir = os.getcwd()  # Get the current working directory
    cellsnake_root = os.path.expanduser("~/cellsnake")  # Trace the directory structure back to the ~/cellsnake/ folder
    # Check if the current working directory is within ~/cellsnake/ or its subdirectories If not found, fall back to default '~/cellsnake/'
    if current_dir.startswith(cellsnake_root):
        return cellsnake_root
    else:
        return cellsnake_root

cellsnake_path = get_cellsnake_path() # Get the final `cellsnake_path`
yaml_file_path = os.path.join(cellsnake_path, "cellsnake/scrna/config.yaml") # Define paths dynamically based on `cellsnake_path`

# Load the .yaml file (if it exists, otherwise use an empty dict)
if os.path.exists(yaml_file_path):
    with open(yaml_file_path, 'r') as file:
        config_from_yaml = yaml.safe_load(file) or {}
else:
    config_from_yaml = {}

# Default values
defaults = {
    "option": "standard",
    "runid": "".join(random.choices("abcdefghisz", k=3) + random.choices("123456789", k=5)),
    "logname": "_".join(["cellsnake", "".join(random.choices("abcdefghisz", k=3) + random.choices("123456789", k=5)), datetime.datetime.now().strftime("%y%m%d_%H%M%S"), "standard", "log"]),

    # Paths

    "tsv_file_path": os.path.join(cellsnake_path, "scrna/params.tsv"),
    "markers_file_path": os.path.join(cellsnake_path, "scrna/markers.tsv"),
    "test_datafolder": os.path.join(cellsnake_path, "scrna/workflow/tests/testData"),

    # GSEA
    "gsea_file": os.path.join(cellsnake_path, "scrna/workflow/bundle/c2.cgp.v2022.1.Hs.symbols.gmt"),
    "gsea_group": "seurat_clusters",

    "cellsnake_path": "",
    "datafolder": "data",
    "analyses_folder": "analyses",
    "results_folder": "results",
    "is_integrated_sample": False,
    "gene_to_plot": [],

    "grid_search": False,

    # Basic parameters
    "min_cells": 3,
    "min_features": 200,
    "max_features": "Inf",
    "max_molecules": "Inf",
    "min_molecules": 0,
    "percent_mt": "10",
    "percent_rp": 0,
    "highly_variable_features": 2000,
    "variable_selection_method": "vst",
    "doublet_filter": "--doublet.filter",
    "metadata": "metadata.csv",
    "metadata_column": "condition",
    "keywords": 1,
    "exact": "--exact",
    "subset_file": "subset",
    "subset_column": "seurat_clusters",
    "min_percentage_to_plot": 5,

    "doublet_filter": True,

    # Clustering and normalization parameters
    "normalization_method": "LogNormalize",
    "scale_factor": 10000,
    "resolution": "0.8",

    # Dimension reduction options
    "umap_plot": "--umap",
    "tsne_plot": "--tsne",
    "show_labels": "--labels",

    # Marker plots    # --tsne instead of False
    "umap_markers_plot": "--umap",
    "tsne_markers_plot": False,

    # Differential expression
    "logfc_threshold": 0.25,
    "test_use": "wilcox",
    "marker_plots_per_cluster_n": 20,
    "identity_to_analysis": ["seurat_clusters"],
    "selected_gene_file": "markers.tsv",

    # Enrichment parameters
    "algorithm": "weight01",
    "statistics": "ks",
    "mapping": "org.Hs.eg.db",
    "organism": "hsa",

    # Integration & cell typing
    "integration_id": "integrated",
    "celltypist_model": "Immune_All_Low.pkl",
    "singler_ref": "BlueprintEncodeData",
    "singler_granulation": "label.main",

    # CellChat
    "species": "human",

    # Kraken DB
    "kraken_db_folder": None,
    "prekraken_db_folder": None,
    "taxa": "genus",
    "microbiome_min_cells": 1,
    "microbiome_min_features": 3,
    "confidence": 0.01,
    "min_hit_groups": 3,
    "kraken_extra_files": False,
    "bowtie_database_prefix": None,

    # Misc
    "complexity": 0,
    "reduction": "cca",
    "dims": 30
}

# Merge YAML config with defaults (YAML values override defaults)
merged_config = {**defaults, **config_from_yaml}

# Convert merged config to an object with dot notation
class Config:
    def __init__(self, config_dict):
        for key, value in config_dict.items():
            setattr(self, key, value)

config = Config(merged_config)