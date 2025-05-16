import subprocess
import os
import pytest
import pandas as pd
from pathlib import Path


# =======================================
# Integration tests for the scrna workflow's scrna-celltypist.py script.
#
# This test runs the celltypist script as a subprocess using example input data,
# then verifies the creation of expected output files such as plots, Excel reports,
# and prediction results.
#
# It also checks basic validity of outputs (e.g., non-empty Excel sheets),
# ensuring the script runs end-to-end without errors.
#
#
# =======================================

# work out current path, then go one steps back to tests directory:  ./tests/py-tests, so testData can be reached
test_dir = Path(__file__).resolve().parents[1]

test_data_path = test_dir / "testData"


@pytest.mark.parametrize("test_file, model, idents", [
    (test_data_path / "analyses" / "defaultTest" / "testSample.h5ad", "Immune_All_Low.pkl", "seurat_clusters"),
])
def test_run_celltypist_script(test_file, model, idents, tmp_path):
    """Integration test for the celltypist CLI script."""

    # Output paths
    dotplot_output = test_data_path / "results" / "dotplot.png"
    xlsx_output = test_data_path  / "results" / "crosstab.xlsx"
    prediction_dir = test_data_path / "results" / "predictions"
    prediction_dir.mkdir(exist_ok=True)

    script_path = Path(__file__).resolve().parents[2] / "scripts" / "scrna-celltypist.py"

    # Run the script as a subprocess
    result = subprocess.run([
        "python", str(script_path),
        str(test_file),
        str(dotplot_output),
        str(prediction_dir),
        str(xlsx_output),
        model,
        idents
    ], capture_output=True, text=True)

    # Debugging info
    print("STDOUT:", result.stdout)
    print("STDERR:", result.stderr)

    assert result.returncode == 0, f"Script failed with code {result.returncode}"

    # Assertions on outputs
    assert dotplot_output.exists(), "Dotplot output was not created."
    assert xlsx_output.exists(), "Crosstab Excel file was not created."
    assert any(prediction_dir.iterdir()), "Prediction output folder is empty."

    # Optional Excel structure check
    df = pd.read_excel(xlsx_output)
    assert df.shape[1] > 0, "Excel file appears empty."



# Create directories if they don't exist
os.makedirs(os.path.dirname(test_data_path / "results" / "scrna_celltypist.txt"), exist_ok=True)

# Write to the file
with open(test_data_path / "results" / "scrna_celltypist.txt", "w") as f:
    f.write("celltypist test ran")
