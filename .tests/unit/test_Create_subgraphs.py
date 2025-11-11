"""
Rule test code for unit testing of rules generated with Snakemake 9.13.7.
"""

import os
import re
import shutil
import sys
import tempfile
from pathlib import Path
from subprocess import check_output

sys.path.insert(0, os.path.dirname(__file__))


def test_Create_subgraphs(conda_prefix):
    with tempfile.TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        config_path = Path(".tests/unit/Create_subgraphs/config")
        data_path = Path(".tests/unit/Create_subgraphs/data")
        expected_path = Path(".tests/unit/Create_subgraphs/expected")

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
                "output/path_test/graphs/simulated/pruned.dbg_subgraphs/graph_0.dbg",
                "output/path_test/graphs/simulated/pruned.dbg_subgraphs/graph_0.sources",
                "output/path_test/graphs/simulated/pruned.dbg_subgraphs/graph_0.sinks",
                "output/path_test/graphs/simulated/pruned.dbg_subgraphs/graph_stats.txt",
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

        # Simple custom comparison that handles the file order correctly
        def compare_files_with_line_skipping():
            # Get the list of expected files
            expected_files = list(expected_path.rglob("*"))
            expected_files = [f for f in expected_files if f.is_file()]

            for expected_file in expected_files:
                # Get the relative path from expected_path
                rel_path = expected_file.relative_to(expected_path)
                generated_file = workdir / rel_path

                print(f"Comparing: {generated_file} vs {expected_file}")

                if not generated_file.exists():
                    raise AssertionError(
                        f"Generated file does not exist: {generated_file}"
                    )

                # Determine how many lines to skip based on file type
                lines_to_skip = _get_lines_to_skip(generated_file)
                _compare_ignoring_first_lines(
                    generated_file, expected_file, lines_to_skip
                )

        def _get_lines_to_skip(file_path):
            """Determine how many lines to skip based on file type/name"""
            filename = file_path.name
            if filename.endswith(".dbg"):
                return 3  # Skip first 3 lines for DBG files
            elif filename.endswith(".sources") or filename.endswith(".sinks"):
                return 2  # Skip first 2 lines for source/sink files
            elif filename.endswith(".txt"):
                return 1  # Skip first line for text files
            else:
                return 0  # Default: don't skip any lines

        def _strip_ansi_codes(text):
            """Remove ANSI escape codes from text."""
            ansi_escape = re.compile(r"\x1B(?:[@-Z\\-_]|\[[0-?]*[ -/]*[@-~])")
            return ansi_escape.sub("", text)

        def _compare_ignoring_first_lines(generated_file, expected_file, lines_to_skip):
            """Compare files while ignoring the first N lines, ANSI codes, and normalizing whitespace."""
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
                        f"Original line counts - Generated: {len(gen_lines)}, Expected: {len(exp_lines)}\n"
                        f"Generated file: {generated_file}\nExpected file: {expected_file}"
                    )

                for i, (gen_line, exp_line) in enumerate(
                    zip(gen_lines_skipped, exp_lines_skipped)
                ):
                    # Remove ANSI escape codes from expected line
                    exp_line_clean = _strip_ansi_codes(exp_line)

                    # Normalize whitespace: replace all whitespace with single spaces and strip
                    gen_normalized = " ".join(gen_line.split()).strip()
                    exp_normalized = " ".join(exp_line_clean.split()).strip()

                    print(f"DEBUG Line {i + lines_to_skip + 1}:")
                    print(f"DEBUG Generated normalized: '{repr(gen_normalized)}'")
                    print(f"DEBUG Expected normalized: '{repr(exp_normalized)}'")
                    print(f"DEBUG Are they equal? {gen_normalized == exp_normalized}")

                    if gen_normalized != exp_normalized:
                        raise AssertionError(
                            f"Files differ at line {i + lines_to_skip + 1}:\n"
                            f"Generated: '{gen_normalized}'\n"
                            f"Expected:  '{exp_normalized}'\n"
                            f"Generated file: {generated_file}\nExpected file: {expected_file}"
                        )

                print(
                    f"Files match (first {lines_to_skip} lines ignored, ANSI codes and whitespace normalized): {generated_file}"
                )

            except Exception as e:
                raise AssertionError(f"Comparison failed for {generated_file}: {e}")

        # Run the custom comparison
        compare_files_with_line_skipping()
