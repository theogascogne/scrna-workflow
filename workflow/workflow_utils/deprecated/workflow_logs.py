import pathlib
from Defaults import config

def write_main_log(files):

    pathlib.Path("logs/").mkdir(parents=True, exist_ok=True)
    with open("logs/" + config.logname,"w") as f:
        f.write("Total number of samples processed : " + str(len(files)) + "\n")
        f.write("Sample names in this run : " + " ".join(files) + "\n")
        f.write("Run ID : " + str(config.runid) + "\n")
        f.write("Option : " + str(config.option) + "\n")

def write_sample_log(sample,paramspace):
    l=config.results_folder + "/" + sample + "/logs/"
    pathlib.Path(l).mkdir(parents=True, exist_ok=True)
    with open(l + config.runid + ".txt","w") as f:
        f.write("Sample name : " + sample + "\n")
        f.write("Run ID : " + str(config.runid) + "\n")
        f.write("Run option : " + str(config.option) + "\n")
        f.write("Main directory and params : " + "".join(list(paramspace.instance_patterns)) + "\n")
        f.write("Minimum cells : " + str(config.min_cells) + "\n")
        f.write("Minimum features (nFeature_RNA) : " + str(config.min_features) + "\n")
        f.write("Maximum features (nFeature_RNA) : " + str(config.max_features) + "\n")
        f.write("Maximum molecules (nCount_RNA) : " + str(config.max_molecules) + "\n")
        f.write("Minimum molecules (nCount_RNA) : " + str(config.min_molecules) + "\n")
        f.write("Percent mitochondrial gene treshold (smaller than) : " + str(config.percent_mt) + "\n")
        f.write("Percent ribosomal gene treshold (larger than) : " + str(config.percent_rp) + "\n")
        f.write("Resolution : " + str(config.resolution) + "\n")
        f.write("Highly variable genes : " + str(config.highly_variable_features) + "\n")
        f.write("Doublet filter : " + str(config.doublet_filter) + "\n")
        f.write("Normalization method : " + str(config.normalization_method) + "\n")
        f.write("Scale factor : " + str(config.scale_factor) + "\n")
        f.write("LogFC treshold : " + str(config.logfc_threshold) + "\n")
        f.write("DE test use : " + str(config.test_use) + "\n")
        #f.write("Algorithm (GO enrichment) : " + str(algorithm) + "\n")
        #f.write("Statistics (GO enrichment) : " + str(statistics) + "\n")
        f.write("Mapping (GO and KEGG enrichment) : " + str(config.mapping) + "\n")
        f.write("Organism (KEGG enrichment) : " + str(config.organism) + "\n")
        f.write("Species (cellchat) : " + str(config.species) + "\n")
        #f.write("Ontology (GO enrichment) : " + str(ontology) + "\n")
        f.write("SingleR reference : " + str(config.singler_ref) + "\n")
        f.write("Celltypist model : " + str(config.celltypist_model) + "\n")
        f.write("Kraken DB folder : " + str(config.kraken_db_folder) + "\n")
        f.write("Collapse to this taxanomic level : " + str(config.taxa) + "\n")
        f.write("Kraken confidence param : " + str(config.confidence) + "\n")
        f.write("Kraken minimum hit groups : " + str(config.min_hit_groups) + "\n")
