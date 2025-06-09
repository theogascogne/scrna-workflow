# scRNA sequencing analysis workflow
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active) 

Introduction
------------

Cellsnake can be run directly using the snakemake workflow. We recommend the wrapper but the snakemake workflow give more control in some use cases.

The main cellsnake repo is here : https://github.com/sinanugur/cellsnake


Installation
------------

You may pull the workflow from the GitHub repo and create a clean environment. Mamba installation is highly recommended.

```
conda install mamba -c conda-forge # to install Mamba

git clone https://github.com/sinanugur/scrna-workflow.git
cd scrna-workflow
mamba env create --name scrna-workflow --file environment.yml
conda activate scrna-workflow
```

For Apple Silicon (i.e. M1, M2 etc.) architecture, you have to put CONDA_SUBDIR=osx-64 before creating the environment.
```
CONDA_SUBDIR=osx-64 mamba env create --name scrna-workflow --file environment.yml
conda activate scrna-workflow
```


After the environment created and activated successfully, to install all the required R packages, you should run the installation script, this may take some time:
```
bash install_r_packages.sh
```


Quick Start Example
-------------------

You can start a minimal run by calling, sample runs are expected in data folder.

```shell
snakemake -j 10 --config datafolder=data option=minimal
```

Then we can run integration.
```shell
snakemake -j 10 --config option=integration
```

Now it is time to work on the integrated sample. We can run standard workflow on the integrated object which is always generates at the same location.
```shell
snakemake -j 10 --config  datafolder=analyses_integrated/seurat/integrated.rds option=standard is_integrated_sample=True --rerun-incomplete
```

You may change some of the options or you may provide a config file as well, for example.
```shell
snakemake -j 10 --config  datafolder=analyses_integrated/seurat/integrated.rds option=standard is_integrated_sample=True --configfile=config.yaml --rerun-incomplete
```


If you want to contribute to the workflow, a guide can be found under here!

The relevant files for development are:

------------------------
Workflow folder
------------------------

-enviornment.yml   
this should be updated if you install packages to your environment, but only if you wish to share new features you have made.

-config.yaml      
this is needed for the workflow's handling of variables, if you change something in here then the workflow will use the newly changed values too
do note that config.yaml is not entirely needed for testing or running the workflow as it also uses default values from defaults.py. It prioritizes config.yaml
If you change config.yaml and forget what you changed, please refer to defaults.py
                   
-global.r         
This is part of an execution experiment. On it's own, this script will monitor a text file 'trigger.txt' for commands that are formatted as shell scripts. an example of this:
--Rscript <path to script> --rds <path to rds> --value1 <value> .....
Everything being in key pairs much like how the snakemake shell commands are structured, save for the 'Rscript' also being a key pair value for global.r
This experimental functionality is unstable and may not work at times, but it is designed to bypass the overhead that R-based workflows tend to introduce trough constant initialization
of new R environments and library loading for each rule.

-persistent.sh     
This is the second half of the execution experiment. This will launch global.r as a background process tied to the terminal. It also features a pause functionality when errors are detected
by snakemake where it will pause for 60 seconds before resuming, either shutting down the workflow or skipping the failed rule. It says you can pick which outcome you want in the event of an error
but this does not work as snakemake still forces its error handling to take priority.
Back to its relation to global.r: persistent.sh takes in a snakemake command, turns the Rscript segment at the front into a key pair value and writes it to trigger.txt (created by the sh script)                   so global.r can read it and execute the R script with the given arguments.
To activate this functionality, turn the flag 'background' inside persistent.sh from FALSE to TRUE (other flags in the script will be set to FALSE if you do this).
If every flag in the script is set to FALSE, the script runs the commands as snakemake would, effectivly turning itself off.

-install_r_packages.sh
This simply runs the R script responsible for installing packages. It is called on by the CLI command from the wrapper. You can run the R script yourself if you like
script: workflow/scripts/scrna-install-packages.R


-Snakefile         
The snakefile is the core of the workflow. It instructs snakemake to build a dynamic DAG based on the value of the option parameter in the config. The top of the file is filled with configuration                           imports and declarations. It pulls from both config.yaml and a Python-based defaults file (default.py found in worfklow_utils folder), ensuring consistency across rules and helper scripts.
The snakefile first imports several modules and helper functions, and includes a rules file called extra_functions.smk. It also pulls in utility functions for plotting and parameter handling. These                         imported modules provide functions like file_capture_try and sample_parameter, which help with dynamic file discovery and parameterization of input combinations.
Next, the Snakefile reads configuration values. These include data paths, filtering thresholds, metadata column names, plotting flags, and other analysis-specific toggles. Each setting is pulled from                       either configDef or config, depending on whether it’s coming from Python or YAML. This dual-source system lets you prioritize defaults in Python while still allowing user overrides through a config file.
Then it checks the value of option to decide what mode the pipeline should run in. If the option is test, it reads an additional parameter test_mode, includes the seuratTest.smk rules, and sets the
final target to a test file. If the option is integration or integrate, it uses cellsnake_glob_wildcards to gather RDS files, selects one per sample if needed, and includes integration.smk. In this
mode, it checks for multiple RDS files per sample and handles ambiguity by printing warnings and instructions to clean up.
If the option is something like minimal, standard, or advanced, it includes both the seurat.smk and microbiome.smk rule files, and uses the sample_parameter function to build a series of file paths                         which acts as instructions on how the snakemake should build its DAG.
Finally, if no valid option is given or no files are found, it prints a generic error message prompting the user to check the config or input files. No exception is thrown—just a soft fail with a pass.



---------------------------------------------------
Worfklow_utils folder
---------------------------------------------------
This folder holds functions that support the workflow, either trough making plot paths for the Rule ALLs found in the snakefile or trough finding sample names and populating the snakemake wildcards with them.

-Daemon.py
A precursor to the exection experment (persistent.sh). This can be ignored or removed.

-Defaults.py       
This defaults.py module acts as the centralized configuration hub for the scRNA-seq pipeline. It ensures consistency between all Python scripts and the Snakefile by merging a YAML-based config file with                    a dictionary of hardcoded fallback defaults. If no YAML is found or it will use the dictionary in here instead.
At the top, it tries to determine the root installation path of the cellsnake package dynamically using Python’s import system. If the import fails (maybe you're running the workflow outside a proper                       environment), it falls back to a default location like ~/cellsnake.
Using that path, it then constructs the expected location of config.yaml. If the file is present, it's parsed using yaml.safe_load, otherwise it proceeds with an empty configuration.
The heart of the script is the defaults dictionary. This dictionary contains every parameter used across the workflow: file paths, feature filtering thresholds, normalization settings, plotting flags,                      and experimental options like KrakenDB and integration strategies. It’s comprehensive enough that the workflow can run even if no YAML file is supplied, using reasonable assumptions for all variables.
Finally, the merged dictionary is wrapped into a Config class, which allows dot-notation access to all parameters instead of using key lookups. So config.min_cells becomes functionally equivalent to                        merged_config["min_cells"].

It should be noted that this method of finding the workflow's root folder is used alot, so you can run the workflow's indivdual files as standalone. This include scripts, python files and tests.


-create_plot_functions.py
This script is a dynamic output path generator used at the end of a Snakemake pipeline. Its job is to determine what plots, tables, and analysis outputs should be created based on the configuration and                     the selected workflow mode (like "minimal", "standard", or "advanced").
It works by first pulling in information from the configuration file, such as which identities to analyze, what genes to focus on, whether external tools like CellTypist or Kraken should be run, and
where the data is located. It then uses this information to generate file paths that correspond to the expected results—for example, dimensionality reduction plots, gene expression maps, enrichment
analyses, or clustering metrics.
The central function, sample_parameter, acts as the controller. Depending on which workflow mode is selected, it calls different helper functions that generate file paths for specific types of plots or
results. For example, identity_dependent_dimplot returns file paths for UMAP, t-SNE, and PCA plots grouped by identity; selected_gene_plot focuses on gene-specific visualizations; and kraken_predictions 
handles microbiome-related outputs.
All of these functions use a core utility called generate_plots, which takes the plot type, file list, identities, and optional gene names or extensions and returns the appropriate paths. These paths
are then returned to Snakemake, which uses them to construct the DAG—ensuring only the necessary rules are executed to produce the desired outputs.

-extra_functions.py
This script provides utility functions that assist in pattern matching and wildcard extraction for use within a Snakemake workflow, particularly for dynamic file discovery. It's primarily meant to be
used by the create_plot_functions.py script, and its main role is to help identify and extract sample-related wildcards from the filesystem.
At its core, it defines a custom version of Snakemake’s glob_wildcards function called cellsnake_glob_wildcards, which is more flexible and adapted to specific use cases in the workflow. This function
takes a file path pattern that includes Snakemake-style wildcards (like {sample}, {params}, etc.) and walks through the directory tree to find actual files or directories that match this pattern. It
then extracts the values for those wildcards based on the matches it finds, building a named tuple where each wildcard corresponds to a list of its observed values.
                   
-global.r and interrhandle.sh
These are no longer used, as they are the previous versions of the execute experiment. They were put in to show progress of the experiment if needed. They can be removed or ignored.


-workflow_funcs.py 
The purpose of this script is to find sample identifiers (wildcards) by scanning a given datafolder for files that match known 10x Genomics formats. These formats typically include:
matrix.mtx.gz or matrix.mtx in subfolders like outs/filtered_feature_bc_matrix
filtered_feature_bc_matrix.h5
fallback paths like raw_feature_bc_matrix or root-level matrix.mtx.gz
The function file_capture tries to match these file patterns using a custom globbing utility (cellsnake_glob_wildcards) and returns the list of sample names ({sample} wildcard values) that it finds.
Then, file_capture_try adds a layer of error handling and flexibility. It:
First checks if the input path is a directory and tries to detect samples in it.
If that fails, it checks the parent directory and filters out irrelevant paths.
If a file is passed instead, it handles specific extensions (like .rds), checking for special cases such as integrated datasets or subsetting scenarios (depending on the config).
Finally, it prints and returns a unique set of detected sample names.
The last function, initialization_of_paramspace, prepares a DataFrame that defines parameter combinations for grid search (or defaults). If a TSV file is provided and grid search is enabled in the
config, it reads that table. Otherwise, it falls back to a default dictionary and converts it into a DataFrame.


-workflow_logs.py
The purpose of this script is to print out a detailed logging of a workflow. This is currenly not in use.


The deprecated folder is more of the experiment's 'history'. Ignore or remove it.



-------------------------------------------
The rules folder
-------------------------------------------
This folder has every .smk that snakemake uses. The most prominent one is seurat.smk as it contains the majority of the workflow's rules. The other .smk are related to the advanced workflow, which i have not touched on yet. Documenting them is a TO DO.

-seurat.smk        
This Snakemake rule file orchestrates key steps in the scrna analysis pipeline, handling everything from raw data processing to clustering and visualization. It includes logic to detect and adapt to
various input formats (e.g., 10x .mtx or .h5), and controls how samples are initialized, filtered, normalized, and annotated. Intermediate and final outputs include .rds Seurat objects, cluster
visualizations, technical QC plots, and marker gene tables. 

-seuratTest.smk    
This set of rules are ran in cellsnake's test command. It runs relevant test scripts found in the test folder. Since the tests are copies of their corresponding scripts, this is effectively running
the workflow in a test 'sandbox'.

-test.smk          
Remnant from development, can be removed.



-------------------------------------------
The scripts folder
-------------------------------------------
This folder has every script used during workflow and a folder for their respective helper function. For the functions that have been refactored, they follow a similar structure:

1. declare libraries

2. opt parse handling and check that crucial dependency is in order (rds, sampleid)

3. source scripts with needed functions (helper function files)

4. run the actual workflow steps


The helper functions will give an idea of which scripts have been refactored.

SeuratObject has alot of usage in these scripts!

Also this is a quick guide for these.


-scrna-read-qc.R   
This script handles quality control of the scrna data.
min.cells filters features based on minimum amount of cells that expresses it
min.feature filters cells based on minimum of features that it expresses
percent mt and rp filters cells based on what percentage of genes (features) it expresses that belong to a naming convention (MT or RP in their names). It is mitochondrial and ribosomal genes.
It then vizualises/plots the data based on the distribution of the four values: Features, counts, percent MT and percent RP


-scrna-normalization-pca.R
It loads the Seurat object that is stored as an .rds from scrna-read-qc.R. It normalizes the data, selects highly variable features, and performs PCA. Based on the principal components, it clusters the cells using user-defined resolution parameters. If requested, it performs dimensionality reduction via UMAP and t-SNE, and can also apply doublet detection and filtering.
The script then uses SingleR with a specified reference (e.g., HumanPrimaryCellAtlasData) to annotate cell types, and saves the updated Seurat object to the designated output .rds file. All computationally intensive steps are configurable via CLI options, including choice of normalization method, number of features, clustering resolution, and whether to apply integration or dimensionality reduction.

-scrna-clusteringtree.R
This R script processes a Seurat object to assess and visualize clustering resolution and principal component contributions in single-cell RNA-seq data. It accepts command-line arguments to control parameters like normalization method, number of variable features, and output file names for several diagnostic plots.
It reads in the raw Seurat .rds file. Depending on whether integration is requested, it applies the appropriate normalization or integration procedure before scaling the data and performing PCA.
If integration is not applied, the script identifies and visualizes the top variable features, and produces JackStraw and elbow plots to assess the significance of principal components. It determines a resolution range based on integration status, constructs a shared nearest neighbor (SNN) graph, and performs clustering.

-scrna-convert-to-h5ad.R
This script facilitates the conversion of a Seurat .rds object to an .h5ad file format, which is commonly used in Python-based single-cell analysis workflows such as Scanpy. It is designed to integrate seamlessly with the cellsnake workflow system, dynamically sourcing helper functions and managing dependencies. Unfortunately because This is incompatible with Assay5, an additional tool called scCustomize is needed to convert it back to
assay v3/4 to then be converted into a h5ad

-scrna-dimplot.R
This script generates dimensionality reduction plots (UMAP or t-SNE) from a processed Seurat object and optionally overlays custom metadata to refine visualization.
Once the Seurat object is loaded, the script ensures the default assay is set to RNA. If metadata is provided via a CSV or RDS file, it's merged into the Seurat object's metadata using specialized helper functions.
The identity of cells is set according to a specified metadata column. It then filters clusters based on a user-defined minimum percentage, allowing only prominent clusters to be plotted. A color palette is generated based on cluster identities to maintain consistent and distinguishable visual representation.
The script produces two output files: an interactive HTML plot using Plotly for web-based visualization, and a static PDF plot using ggplot2 for inclusion in reports or publications.

-scrna-find-markers.R
This script identifies marker genes for clusters in a Seurat object and saves the results as an RDS file. The required input is a processed Seurat .rds file. The script reads the Seurat object, sets the default assay to "RNA", assigns identities based on metadata, runs finds markers with the given parameters, and saves the results. A helper functions file is included but not actively used in this script.

-scrna-go_analysis.R
This script performs Gene Ontology (GO) and Gene Set Enrichment (GSE) analysis on marker gene data from an Excel file. It expects a table of markers with cluster assignments, applies filtering based on user-defined p-value and log fold change thresholds, and performs enrichment analysis for each cluster using annotation data from a specified organism database (e.g., org.Hs.eg.db). The results are saved as an RDS file and optionally exported as Excel files for both GO and GSE outputs. It leverages helper functions defined in the cellsnake workflow, which are sourced from the installed Python package directory.

-scrna-kegg.R
This R script performs KEGG and metabolic pathway enrichment analyses—including both over-representation (enrichKEGG/enrichMKEGG) and gene set enrichment (gseKEGG/gseMKEGG)—on marker gene clusters from an Excel file. It supports custom organism codes and gene-to-pathway mappings via Bioconductor annotation databases (e.g., org.Hs.eg.db). The script splits the input markers by cluster and uses a helper function to run all enrichment analyses, saving the outputs to RDS and Excel files using designated paths. It sources additional logic from the cellsnake package structure to handle the enrichment workflows.

-scrna-marker-heatmap.R
This R script processes a Seurat object to produce a marker-based heatmap and an average gene expression table. It scales the RNA assay data using a helper function that extracts relevant marker genes from a provided Excel file. The script sets cell identities based on a metadata column and identifies marker genes not present in all clusters. It then generates a heatmap of these genes and saves it as a PDF. Finally, it calculates the average expression of each gene across the identity groups and writes the results to an Excel file.

-scrna-marker-tables.R
This R script reads an RDS file containing a marker gene data frame, filters for positively expressed markers, and writes two Excel files. One contains the complete set of markers, and the other includes only those with positive average log2 fold-change. It uses the openxlsx package for writing Excel files. Be aware it might give an error that a pdf editor is not available.

-scrna-metrics.R
This R script processes the Seurat object by loading marker and cluster information to generate summary metrics and visualizations. It reads the Seurat RDS file and uses a specified metadata column (defaulting to "seurat_clusters") to identify clusters. It then generates a summary table of cell counts per cluster, which is saved as a PDF. The script proceeds to compute a detailed cluster summary and saves it as both an Excel and a TSV file. A static bar plot of cluster distributions is created and saved as a PDF. 

-scrna-monocle3.R
The script performs trajectory inference using Monocle3 on a Seurat object. The script converts the seurat object to a Monocle3 cell_data_set (CDS), and performs clustering. It generates a combined plot displaying partitions and singleR annotations, saving it to the specified file. It then processes each partition within the CDS and saves partition-specific plots and data outputs to a user-defined directory.

-scrna-singler-annotation.R
performs automatic cell type annotation on a Seurat object using the SingleR method. The script reads the Seurat object, sets the default assay to RNA, converts the object to a SingleCellExperiment, and uses SingleR to predict cell types based on a specified reference dataset and annotation granularity. The predicted labels are saved as an RDS file to the user-defined output path.

-scrna-singler-plots.R
Generates diagnostic plots and summary tables based on SingleR cell type annotations for a Seurat object. It requires and loads a Seurat object and a corresponding SingleR prediction object. It generates a full and top-20 score heatmap from prediction scores, builds a cluster-by-cell-type assignment table, and creates a heatmap to visualize these associations. Finally, it saves this table as an Excel file.

-scrna-technicals.R
loads a Seurat object from a specified RDS file and creates four quality control feature plots to visualize key metrics: the number of detected features (genes), total UMI counts, mitochondrial gene percentage, and ribosomal gene percentage. It generates each plot with FeaturePlot, using a small point size and no rasterization for clarity, then saves each plot to a PDF file.

-scrna-top-marker-plot.R
script reads an Excel file containing marker gene data, extracts top marker genes, identifies clusters from the data, and generates a plot saved as a PDF. 

----------------
Test folder
---------------

The test folder contain Python and R tests. They contain the same functions as seen in the scripts and workflow utils folders. You can run the Python tests using Pytest <path to file> or Rscript <path to file>

If you wish to run or contribute to the R tests withing Rstudio, you will not have access to the green arrow run button on the corner of the code quadrant in Rstudio
Instead use Ctrl + Enter to run individual test_that() blocks or individual lines of code.

Because the workflow is procedural computation, global values are used as each step of the process depends on the results of the previous steps! Therefore intermediary objects are needed for the tests.

-------------
testData
-------------

This folder is what the test command for the workflow uses. If you do delete this or demove the datafolder inside it, you can generate a new dataset using 'generate_10X_data.R'. This script needs the features.tsv and barcodes.tsv from a working 10x genomics dataset. One can be downloaded from the link provided in Cellsnake's (the wrapper) Read Me.

You can also generate your own custom data or a working one using this script. I reccomend you use it to mock faulty data for stress testing or boundary value analysis.
