"""
Rule test code for unit testing of rules generated with Snakemake 9.12.0.
"""

import os
import shutil
import sys
import tempfile
from pathlib import Path
from subprocess import check_output

sys.path.insert(0, os.path.dirname(__file__))


def test_Compare_2(conda_prefix):
    with tempfile.TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        config_path = Path(".tests/unit/Compare_2/config")
        data_path = Path(".tests/unit/Compare_2/data")
        expected_path = Path(".tests/unit/Compare_2/expected")

        # Copy config to the temporary workdir.
        shutil.copytree(config_path, workdir)

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir, dirs_exist_ok=True)

        # Test Symlinks
        project_dir = Path.cwd()
        symlink_items = [
            "findviralstrains_2.smk",
            "findviralstrains.smk",
            "findviralstrainsMain.smk",
            "config_files",
            "libs",
            "reference_genomes",
            "output",
        ]
        for item in symlink_items:
            src = project_dir / item
            dst = workdir / item
            if src.exists() and not dst.exists():
                dst.symlink_to(src)

        # Run the test job.
        check_output(
            [
                "python",
                "-m",
                "snakemake",
                "output/path_test/output_genomes/simulated/subgraph_0/simulated_1_of_2_vs_ref.txt",
                "output/path_test/output_genomes/simulated/subgraph_0/simulated_2_of_2_vs_ref.txt",
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

        # Imports for checker #
        import common

        # ignores the first 12 lines of output file #
        class FirstLinesIgnore(common.OutputChecker):
            def compare_files(self, generated_file, expected_file):
                print(f"Comparing with first 12 lines ignored: {generated_file}")
                self._compare_ignoring_first_lines(
                    generated_file, expected_file, lines_to_skip=12
                )

            def _compare_ignoring_first_lines(
                self, generated_file, expected_file, lines_to_skip
            ):
                """Compare files while ignoring the first N lines of both files."""
                try:
                    with open(generated_file, "r") as gen_f:
                        gen_lines = gen_f.readlines()

                    with open(expected_file, "r") as exp_f:
                        exp_lines = exp_f.readlines()

                    gen_lines_skipped = gen_lines[lines_to_skip:]
                    exp_lines_skipped = exp_lines[lines_to_skip:]

                    if len(gen_lines_skipped) != len(exp_lines_skipped):
                        raise AssertionError(
                            f"Files have different number of lines after skipping first {lines_to_skip} lines: "
                            f"{len(gen_lines_skipped)} vs {len(exp_lines_skipped)}\n"
                            f"Original line counts - Generated: {len(gen_lines)}, Expected: {len(exp_lines)}"
                        )

                    for i, (gen_line, exp_line) in enumerate(
                        zip(gen_lines_skipped, exp_lines_skipped)
                    ):
                        if gen_line != exp_line:
                            raise AssertionError(
                                f"Files differ at line {i + lines_to_skip + 1}:\n"
                                f"Generated: {gen_line.strip()}\n"
                                f"Expected:  {exp_line.strip()}"
                            )

                    print(
                        f"Files match (first {lines_to_skip} lines ignored): {generated_file}"
                    )

                except Exception as e:
                    raise AssertionError(f"Comparison failed for {generated_file}: {e}")

        FirstLinesIgnore(data_path, expected_path, workdir).check()
