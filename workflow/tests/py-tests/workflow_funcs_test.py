import os
import sys
#import yaml
import pytest
import pandas as pd
from pathlib import Path


# =========
# Unit tests for workflow utility functions in workflow_funcs.py.
#
# These tests verify core file handling and parameter initialization functionality:
# - file_capture and file_capture_try for detecting sample files in directories
# - initialization_of_paramspace for generating parameter grids from config or TSV inputs
#
# The tests check both expected outputs on valid inputs and proper handling of edge cases,
# including missing directories and parameter consistency.
# =========

# Ensure module path is recognized
# work out current path, then go two steps back to workflow directory:  workflow/tests/py-tests, so defaults can be found inside workflow_utils.
workflow_dir = Path(__file__).resolve().parents[2]
sys.path.append(str(workflow_dir))

# Import functions to be tested
from workflow_utils.workflow_funcs import file_capture, file_capture_try, initialization_of_paramspace
from workflow_utils.Defaults import config

tsv_file_path = config.tsv_file_path
datafolder = config.test_datafolder

# Extract parameters
option = config.option
resolution = config.resolution
percent_mt = config.percent_mt

# --- TESTS ---

@pytest.mark.parametrize("test_datafolder, expected_files", [
    (datafolder, {"testSample"}),  # Valid test data
    (os.path.expanduser("huehuehue"), set()),  # Nonexistent directory
])
def test_file_capture_try(test_datafolder, expected_files):
    """Test file_capture_try function, which includes file_capture."""
    result = file_capture_try(test_datafolder)
    assert result == expected_files, f"Expected {expected_files}, but got {result}"

@pytest.mark.parametrize("test_datafolder", [datafolder])
def test_file_capture(test_datafolder):
    """Test file_capture separately in case of direct usage elsewhere."""
    if not os.path.exists(test_datafolder):
        pytest.fail(f"Test data folder does not exist: {test_datafolder}")

    filesFound = file_capture(test_datafolder)

    # Debugging output
    print("Files found in", test_datafolder)
    for file in filesFound:
        print(file)

    expected_samples = {'testSample'}
    assert set(filesFound) == expected_samples, f"Expected {expected_samples}, but got {filesFound}"

@pytest.mark.parametrize("grid_search, tsv_exists, expected_percent_mt, expected_resolution", [
    (True, True, 10, 0.8),  # Using TSV
    (False, False, 10, 0.8),  # Using dictionary
])
def test_initialization_of_paramspace(grid_search, tsv_exists, expected_percent_mt, expected_resolution):
    """Test initialization_of_paramspace function."""
    
    result = initialization_of_paramspace(tsv_file_path, {"percent_mt": [percent_mt], "resolution": [resolution]})
    
    assert isinstance(result, pd.DataFrame)
    assert not result.empty

    if grid_search and tsv_exists:
        expected_tsv = pd.read_table(tsv_file_path, nrows=1)
        expected_tsv.rename(columns={"MT": "percent_mt"}, inplace=True)
        expected_dict = expected_tsv.to_dict(orient='list')

        result_dict = {
            "percent_mt": [int(value) for value in result["percent_mt"]],
            "resolution": [float(value) for value in result["resolution"]],
        }

        print("Expected Dictionary:", expected_dict)
        print("Result Dictionary:", result_dict)

        assert result_dict == expected_dict
    else:
        expected_dict = {"percent_mt": [int(expected_percent_mt)], "resolution": [float(expected_resolution)]}
        result_dict = {"percent_mt": [int(value) for value in result["percent_mt"]],
                       "resolution": [float(value) for value in result["resolution"]]}

        print("Expected Dictionary:", expected_dict)
        print("Result Dictionary:", result_dict)

        assert result_dict == expected_dict
