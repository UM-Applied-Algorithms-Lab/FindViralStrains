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


def test_trim_and_merge_raw_reads(conda_prefix):

    with tempfile.TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        config_path = Path(".tests/unit/trim_and_merge_raw_reads/config")
        data_path = Path(".tests/unit/trim_and_merge_raw_reads/data")
        expected_path = Path(".tests/unit/trim_and_merge_raw_reads/expected")

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
                "output/path_test/read_data/trimmed/simulated.merged.fq.gz",
                "output/path_test/read_data/trimmed/simulated.nomerge.pair.R1.fq.gz",
                "output/path_test/read_data/trimmed/simulated.nomerge.pair.R2.fq.gz",
                "output/path_test/read_data/trimmed/simulated.nopair.R1.fq.gz",
                "output/path_test/read_data/trimmed/simulated.nopair.R2.fq.gz",
                "output/path_test/.fastp_logs/fastp/simulated_trim_fastp.html",
                "output/path_test/.fastp_logs/fastp/simulated_trim_fastp.json",
                "--snakefile",
                "findviralstrains.smk",
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

        class HTMLIgnoringOutputChecker(common.OutputChecker):
            def compare_files(self, generated_file, expected_file):
                # Skip comparison for HTML files
                if str(generated_file).endswith('.html'):
                    print(f"Skipping HTML file comparison: {generated_file}")
                    return
                # For all other files, use the original comparison method
                super().compare_files(generated_file, expected_file)

        HTMLIgnoringOutputChecker(data_path, expected_path, workdir).check()
