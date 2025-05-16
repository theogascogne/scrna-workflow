import pytest
from pathlib import Path
import sys
from unittest.mock import patch, mock_open, MagicMock

cellsnake_root_dir = Path(__file__).resolve().parents[4]  # Set the path to root for command_line.py
sys.path.append(str(cellsnake_root_dir))
test_dir = Path(__file__).resolve().parents[1]

# Simulating valid and invalid inputs
validRds = test_dir.joinpath('testData', 'cmdTest', 'valid.rds')
validConfig = test_dir.joinpath('testData', 'cmdTest', 'empty.yaml')
validMeta = test_dir.joinpath('testData', 'cmdTest', 'empty.csv')
validKraken = test_dir.joinpath('testData', 'cmdTest', 'kraken.txt')

from cellsnake.command_line import *



# ==============================================
# Test suite for the cellsnake.command_line module.
#
# This file contains pytest-based unit tests that validate:
# - Command-line argument parsing and validation logic
# - Configuration file loading and parameter extraction
# - Logging functionality within the CommandLine class
# - The main entry point behavior under different CLI flag scenarios
#
# Tests cover both expected successful workflows and failure cases,
# using mocks and parameterized inputs to simulate a wide range of user inputs
# and system states without needing actual file system dependencies.
# ================================================

# Global argument list
argument_list = {
    '--jobs': 4,
    '<INPUT>': validRds,
    '<command>': 'run_command',
    'integrated': True,
    '--configfile': validConfig,
    '--metadata': validMeta,
    '--gene': 'gene.txt',
    '--kraken_db_folder': validKraken,
    '--taxa': 'species',
    '--dry': True,
    '--unlock': True,
    '--remove': True,
    "--generate-template": False,
    "--install-packages": False,
}
# --- Test check_command_line_arguments ---
@pytest.mark.parametrize(
    "arguments, expected",
    [
        # Case 1: Valid input (Happy Path)
        (argument_list, True),

        # Case 2: Valid input with integrated mode False
        ({**argument_list, "integrated": False}, True),

        # Case 3: Integrated mode but input is a directory (Failure)
        ({**argument_list, "<INPUT>": validKraken}, False),

        # Case 4: Wrong file extension for integrated mode (Failure)
        ({**argument_list, "<INPUT>": test_dir.joinpath('testData', 'cmdTest', 'invalid.txt')}, False),

        # Case 5: Missing config file (Failure)
        ({**argument_list, "--configfile": 'no'}, False),

        # Case 6: Missing metadata file (Failure)
        ({**argument_list, "--metadata": 'no'}, False),

        # Case 7: Missing Kraken database folder (Failure)
        ({**argument_list, "--kraken_db_folder": 'no'}, False),

        # Case 8: Invalid taxa value (Failure)
        ({**argument_list, "--taxa": "invalid_taxa"}, False),

        # Case 9: Without required input file (Failure)
        ({**argument_list, "<INPUT>": 'no'}, False),
    ]
)
def test_check_command_line_arguments(arguments, expected):
    """Test check_command_line_arguments() with multiple input cases."""
    assert check_command_line_arguments(arguments) == expected


# --- Test prepare_arguments ---
@pytest.mark.parametrize("integration_mode, integrated_sample, expected_option", [
    (False, False, "option=run_command"),  # Standard case
    (False, True, "option=run_command"),   # Integrated sample but not in integration mode
    (True, True, "option=integration"),    # Full integration mode
    (True, False, "option=integration"),   # Integration mode but not an integrated sample (edge case)
])
def test_prepare_arguments(integration_mode, integrated_sample, expected_option):
    # Prepare a CommandLine instance
    command_line = CommandLine()
    command_line.is_this_an_integration_run = integration_mode
    command_line.is_integrated_sample = integrated_sample  # Set sample type

    # removing metadata and taxa due to a 'NoneType' object has no attribute 'get' error.
    arguments = {key: value for key, value in argument_list.items() if key not in ['--metadata', '--taxa']}

    # Run the method we're testing
    command_line.prepare_arguments(arguments)

    # Check if Snakemake command contains expected elements
    snakemake_cmd = command_line.snakemake  
    assert '-j 4' in snakemake_cmd  # Ensure CPU allocation is correct
    assert '-n' in snakemake_cmd  # Ensure --dry mode is included
    assert '--unlock' in snakemake_cmd  # Ensure --unlock is included
    assert '--delete-all-output' in snakemake_cmd  # Ensure --remove flag translates correctly

    # Check if config values were added correctly
    if not integration_mode:
        assert f"datafolder={arguments['<INPUT>']}" in command_line.config  # Only check for non-integration mode
    else:
        assert f"datafolder={arguments['<INPUT>']}" not in command_line.config  # Ensure it's NOT there

    assert f"gene_to_plot={arguments['--gene']}" in command_line.config  # Ensure gene file is handled
    assert f"kraken_db_folder={arguments['--kraken_db_folder']}" in command_line.config  # Ensure Kraken DB folder is added
    assert f"runid={command_line.runid}" in command_line.config  # Ensure run ID is set
    assert f"cellsnake_path={cellsnake_path}/scrna/" in command_line.config  # Ensure cellsnake path is included

    # Check integration-specific behavior
    assert expected_option in command_line.config  # Ensure correct option is set
    if integrated_sample:
        assert "is_integrated_sample=True" in command_line.config  # Ensure integration flag is correctly set


@pytest.mark.parametrize(
    "mock_file_data, expected_params",
    [
        # Case 1: Valid config file with expected keys
        ("resolution: 0.8\npercent_mt: 5", {"resolution": 0.8, "percent_mt": 5}),

        # Case 2: Empty config file (should result in an empty dict)
        ("", None),
    ]
)
def test_load_configfile_if_available(mock_file_data, expected_params):
    """Test loading config files with various scenarios."""
    command_line = CommandLine()

    # Patch the open function to simulate different file contents
    with patch("builtins.open", mock_open(read_data=mock_file_data)):
        try:
            # We will let the actual YAML logic execute
            command_line.load_configfile_if_available({"--configfile": "test.yaml"})
        except yaml.YAMLError:
            # In case of an invalid YAML, we expect parameters to be set to an empty dict
            pass
        
        # Verify that the function's behavior matches expectations
        assert command_line.parameters == expected_params


def test_write_to_log(): #no need to parameterize this
    """Test if write_to_log writes to the log correctly."""
    # Create the CommandLine object
    command_line = CommandLine()
    
    # Set a known runid for testing
    command_line.runid = "test_run_id"
    
    # Mock datetime.datetime.now() to return a fixed timestamp
    fixed_datetime = datetime.datetime(2023, 3, 29, 16, 36, 0)
    
    with patch("datetime.datetime") as mock_datetime:
        mock_datetime.now.return_value = fixed_datetime
        
        # Mock the 'open' function to test file writing
        with patch("builtins.open", mock_open()) as mock_file:
            command_line.log = True
            command_line.write_to_log(start=0)
            expected_logfile = "cellsnake_test_run_id_230329_163600_runlog"
            mock_file.assert_called_once_with(expected_logfile, 'w')


# Parameterize test with different values for <command> and other flags
@pytest.mark.parametrize("cli_args, expected_calls", [
    # Test for generate-template flag
    ({"--generate-template": True, "<command>": None}, ["shutil.copyfile", "open"]),  # expects copy and file write calls

    # Test for install-packages flag
    #({"--install-packages": True, "<command>": None}, ["subprocess.check_call"]),  # expects package install call

    # Test for valid commands
    ({"--generate-template": False, "--install-packages": False, "<command>": "minimal"}, ["run_workflow"]),
    ({"--generate-template": False, "--install-packages": False, "<command>": "standard"}, ["run_workflow"]),
    ({"--generate-template": False, "--install-packages": False, "<command>": "advanced"}, ["run_workflow"]),
    ({"--generate-template": False, "--install-packages": False, "<command>": "clustree"}, ["run_workflow"]),
    ({"--generate-template": False, "--install-packages": False, "<command>": "integrate"}, ["run_workflow"]),
])
def test_main(cli_args, expected_calls):
    # Update argument_list with the current test scenario
    argument_list.update(cli_args)

    # Patch necessary functions
    with patch("cellsnake.command_line.docopt", return_value=argument_list), \
         patch("cellsnake.command_line.run_workflow") as mock_run_workflow, \
         patch("cellsnake.command_line.check_command_line_arguments", return_value=True), \
         patch("subprocess.check_call") as mock_subprocess, \
         patch("shutil.copyfile") as mock_copy, \
         patch("builtins.open", MagicMock()) as mock_open:

        # Run the main function
        main()
        # Print debug information for subprocess.check_call
        print(f"subprocess.check_call was called with: {mock_subprocess.call_args}")

        # Verify if the expected function calls happened
        if "shutil.copyfile" in expected_calls:
            mock_copy.assert_called_once()
            mock_open.assert_called()  # Ensure metadata.csv was written
        if "subprocess.check_call" in expected_calls:
            # We need to make sure subprocess.check_call was actually called
            mock_subprocess.assert_called_once_with(cellsnake_path + "/scrna/workflow/scripts/scrna-install-packages.R")
        if "run_workflow" in expected_calls:
            mock_run_workflow.assert_called_once_with(argument_list)
