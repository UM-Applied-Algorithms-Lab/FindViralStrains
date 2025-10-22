"""
Rule test code for unit testing of rules generated with Snakemake 9.12.0.
"""


import os
import sys
import shutil
import tempfile
from pathlib import Path
from subprocess import check_output

sys.path.insert(0, os.path.dirname(__file__))


def test_Compare_1(conda_prefix):

    with tempfile.TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        config_path = Path(".tests/unit/Compare_1/config")
        data_path = Path(".tests/unit/Compare_1/data")
        expected_path = Path(".tests/unit/Compare_1/expected")

        # Copy config to the temporary workdir.
        shutil.copytree(config_path, workdir)

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir, dirs_exist_ok=True)

        # Run the test job.
        check_output(
            [
                "python",
                "-m",
                "snakemake",
                "output/path_test/output_genomes/simulated/subgraph_0/simulated_1_of_1_vs_ref.txt",
                "--snakefile",
                "findviralstrains_2.smk",
                "-f",
                "--notemp",
                "--show-failed-logs",
                "-j1",
                "--target-files-omit-workdir-adjustment",
                "--configfile",
                "config_files/path_test.yml",
                "--directory",
                workdir,
            ]
            + conda_prefix
        )

        # Check the output byte by byte using cmp/zmp/bzcmp/xzcmp.
        # To modify this behavior, you can inherit from common.OutputChecker in here
        # and overwrite the method `compare_files(generated_file, expected_file),
        # also see common.py.
        import common

        # Create a custom OutputChecker that ignores files with timestamps
        class TimestampIgnoringOutputChecker(common.OutputChecker):
            def compare_files(self, generated_file, expected_file):
                # Skip comparison for files that likely contain timestamps
                # This includes comparison result files and other non-deterministic outputs
                skip_patterns = [
                    '_vs_ref.txt',      # Comparison result files
                    '.paths',           # Path files (may contain timestamps)
                    'graph_stats.txt',  # Statistics files
                ]

                file_str = str(generated_file)
                if any(pattern in file_str for pattern in skip_patterns):
                    print(f"Skipping timestamp-containing file: {generated_file}")
                    return
                # For all other files, use the original comparison method
                super().compare_files(generated_file, expected_file)

        # Use our custom checker that ignores timestamp-containing files
        TimestampIgnoringOutputChecker(data_path, expected_path, workdir).check()
