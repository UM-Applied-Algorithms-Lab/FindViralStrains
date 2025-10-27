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
        # Skip first 12 lines of the expected file

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

        # Create a custom OutputChecker that ignores the first 12 lines of files
        class FirstLinesIgnoringOutputChecker(common.OutputChecker):
            def compare_files(self, generated_file, expected_file):
                # Skip first 12 lines for all files
                print(f"Comparing with first 12 lines ignored: {generated_file}")
                self._compare_ignoring_first_lines(generated_file, expected_file, lines_to_skip=12)
            
            def _compare_ignoring_first_lines(self, generated_file, expected_file, lines_to_skip):
                """Compare files while ignoring the first N lines of both files."""
                try:
                    with open(generated_file, 'r') as gen_f:
                        gen_lines = gen_f.readlines()
                    
                    with open(expected_file, 'r') as exp_f:
                        exp_lines = exp_f.readlines()
                    
                    # Skip the first 'lines_to_skip' lines from both files
                    gen_lines_skipped = gen_lines[lines_to_skip:]
                    exp_lines_skipped = exp_lines[lines_to_skip:]
                    
                    # Compare the remaining lines
                    if len(gen_lines_skipped) != len(exp_lines_skipped):
                        raise AssertionError(
                            f"Files have different number of lines after skipping first {lines_to_skip} lines: "
                            f"{len(gen_lines_skipped)} vs {len(exp_lines_skipped)}\n"
                            f"Original line counts - Generated: {len(gen_lines)}, Expected: {len(exp_lines)}"
                        )
                    
                    for i, (gen_line, exp_line) in enumerate(zip(gen_lines_skipped, exp_lines_skipped)):
                        if gen_line != exp_line:
                            raise AssertionError(
                                f"Files differ at line {i + lines_to_skip + 1}:\n"
                                f"Generated: {gen_line.strip()}\n"
                                f"Expected:  {exp_line.strip()}"
                            )
                    
                    print(f"Files match (first {lines_to_skip} lines ignored): {generated_file}")
                    
                except Exception as e:
                    raise AssertionError(f"Comparison failed for {generated_file}: {e}")

        # Use our custom checker that ignores first 12 lines of all files
        FirstLinesIgnoringOutputChecker(data_path, expected_path, workdir).check()