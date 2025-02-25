import pytest
import tempfile
import filecmp

from grz_qc.compare_threshold import main


def test_compare_threshold():
    # Define the arguments
    args = [
        "--fastp_json",
        "tests/data/compare_thresholds/a2b6c3d9e8f7123456789abcdef0123456789abcdef0123456789abcdef01234_Blood_DNA_normal.fastp.json",
        "--mosdepth_global_summary",
        "tests/data/compare_thresholds/a2b6c3d9e8f7123456789abcdef0123456789abcdef0123456789abcdef01234_Blood_DNA_normal.whole.mosdepth.summary.txt",
        "--mosdepth_target_regions_bed",
        "tests/data/compare_thresholds/a2b6c3d9e8f7123456789abcdef0123456789abcdef0123456789abcdef01234_Blood_DNA_normal.target_genes.regions.bed.gz",
        "--thresholds",
        "tests/data/compare_thresholds/thresholds.json",
        "--sample_id",
        "a2b6c3d9e8f7123456789abcdef0123456789abcdef0123456789abcdef01234_Blood_DNA_normal",
        "--libraryType",
        "wes",
        "--sequenceSubtype",
        "germline",
        "--genomicStudySubtype",
        "tumor+germline",
    ]

    # Create a temporary file for the output
    with tempfile.NamedTemporaryFile() as tmp_output:
        args.extend(["--output", tmp_output.name])
        tmp_output_path = tmp_output.name

        # Call the main function with the arguments
        main(args)

        # Compare the temporary output file with the expected output file
        expected_output_file = "tests/data/compare_thresholds/a2b6c3d9e8f7123456789abcdef0123456789abcdef0123456789abcdef01234_Blood_DNA_normal.result.csv"
        assert filecmp.cmp(tmp_output_path, expected_output_file, shallow=False), (
            "Output file content does not match the expected content"
        )


if __name__ == "__main__":
    pytest.main()
