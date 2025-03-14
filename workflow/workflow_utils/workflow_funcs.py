import os
import pandas as pd
import pathlib
import sys
import os

# this script is meant to read the given directory and look for files based on the 10x genomics format. The try function checks for other viable dirctories
# in case the .gz or .h5 arent found at first.
# once found, it will use cellsnake_glob_wildcards to create wildcard for each sample found. These wildcards are then used by sample_parameter
# inside create_plot_functions to generate file paths for each of the samples during cellsnake's workflow. 

workflow_dir = pathlib.Path(__file__).resolve()
sys.path.append(str(workflow_dir))
from workflow_utils.extra_functions import cellsnake_glob_wildcards
from workflow_utils.Defaults import config

#testData_dir = os.path.expanduser("~/cellsnake/cellsnake/scrna/workflow/tests/workHere/testData (another copy)")

files=[]
def file_capture(datafolder):
    files = cellsnake_glob_wildcards(datafolder + "/{sample}/outs/filtered_feature_bc_matrix/matrix.mtx.gz").sample + cellsnake_glob_wildcards(datafolder + "/{sample}/outs/filtered_feature_bc_matrix/matrix.mtx").sample + cellsnake_glob_wildcards(datafolder + "/{sample}/outs/filtered_feature_bc_matrix.h5").sample 
    if files:
        return files
    
    files = cellsnake_glob_wildcards(datafolder +  "/{sample}/matrix.mtx.gz").sample + cellsnake_glob_wildcards(datafolder +  "/{sample}/matrix.mtx").sample

    if files:
        return files

    files = cellsnake_glob_wildcards(datafolder +  "/{sample}/raw_feature_bc_matrix/matrix.mtx.gz").sample + cellsnake_glob_wildcards(datafolder +  "/{sample}/raw_feature_bc_matrix/matrix.mtx").sample

    if files:
        return files

    files = cellsnake_glob_wildcards(datafolder + "/{sample}/filtered_feature_bc_matrix.h5").sample
    if files:
        return files
    
    return files


def file_capture_try(datafolder):
    files = []
    try:
        # Check if the input is a directory
        if os.path.isdir(datafolder):
            files = file_capture(datafolder)
            files = list(filter(lambda i: "/" not in i, files))  # Filter out subdirectories

            # If no files are found, check the parent directory
            if len(files) == 0:
                head, tail = os.path.split(datafolder.strip("/"))
                files = file_capture(head)  # Capture files in the parent directory
                files = [x for x in files if tail in x]  # Filter files by original tail name
                if files:
                    datafolder = head  # Update datafolder if files are found

        # Check if the input is a file
        if os.path.isfile(datafolder):
            file_extension = pathlib.Path(datafolder)
            
            # If not an `.rds` file
            if (file_extension.suffix).lower() not in [".rds"]:
                files = [file_extension.stem]

            # Handle integrated sample
            elif config.is_integrated_sample:
                analyses_folder = "analyses_integrated"
                results_folder = "results_integrated"
                files = [file_extension.stem]

            # Handle `.rds` file with subset option
            elif (file_extension.suffix).lower() == [".rds"] and config.option in ["subset"]:
                files = [datafolder]

        # Return unique files
        files = set(files) if isinstance(files, list) else files
        print("Samples detected: " + " ".join(files))

        return files  # Return the files

    except Exception as e:
        print(f"No samples detected or something went wrong: {e}")
        return set()  # Return an empty set if something goes wrong
        
        
def initialization_of_paramspace(tsv_file,dictionary):
    if os.path.isfile(tsv_file) and config.grid_search is True:
        par_df = pd.read_table(tsv_file) # if available and use for all samples
    else:
        par_df =  pd.DataFrame(dictionary) # if not available, create using default numbers and use for all samples
    
    return par_df      
        



