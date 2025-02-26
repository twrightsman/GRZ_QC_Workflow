import pytest
import tempfile
import filecmp

from grz_qc.metadata_to_samplesheet import main


def test_metadata_to_samplesheet():
    # Define the arguments
    args = [
        "--submission_metadata_json",
        "tests/data/WES_tumor+germline.metadata.json",
        "/data/",
    ]

    # Create a temporary file for the output
    with tempfile.NamedTemporaryFile() as tmp_output:
        args.extend([tmp_output.name])
        tmp_output_path = tmp_output.name

        # Call the main function with the arguments
        main(args)

        # Compare the temporary output file with the expected output file
        expected_output_file = "tests/data/WES_tumor+germline.samplesheet.csv"
        assert filecmp.cmp(tmp_output_path, expected_output_file, shallow=False), (
            "Output file content does not match the expected content"
        )


if __name__ == "__main__":
    pytest.main()
