import pytest
import tempfile
import filecmp
from pathlib import Path
import sys
from grz_qc.convert_bed_chrom import main

# Parameterize the test with two scenarios:
@pytest.mark.parametrize("bed_input_name,expected_output_file", [
    (
        "tests/data/convert_bed_chrom/input_ncbi.bed",  
        "tests/data/convert_bed_chrom/expected_output.bed"
    ),
    (
        "tests/data/convert_bed_chrom/input_ucsc.bed",  
        "tests/data/convert_bed_chrom/expected_output.bed"
    )
])
def test_convert_bed_chrom(bed_input_name, expected_output_file, monkeypatch):
    # Define the arguments
    args = [
        "--bed", bed_input_name,
        "--mapping", "tests/data/convert_bed_chrom/hg19_NCBI2UCSC.txt",
    ]

    # Create a temporary file for the output
    with tempfile.NamedTemporaryFile() as tmp_output:
        args.extend(["--output", tmp_output.name])
        tmp_output_path = tmp_output.name

        monkeypatch.setattr(sys, "argv", ["convert_bed_chrom.py"] + args)
        # Call the main function with the arguments
        main()

        # Compare the temporary output file with the expected output file
        assert filecmp.cmp(tmp_output_path, expected_output_file, shallow=False), (
            "Output file content does not match the expected content"
        )

        # Clean up the temporary output file.
        # Path(tmp_output_path).unlink()

if __name__ == "__main__":
    pytest.main()
