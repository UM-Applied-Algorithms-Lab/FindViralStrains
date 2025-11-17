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
    # Create fake_temp directory
    fake_temp_dir = Path("fake_temp")
    fake_temp_dir.mkdir(exist_ok=True)

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

        def compare_files_with_line_skipping():
            expected_files = list(expected_path.rglob("*"))
            expected_files = [f for f in expected_files if f.is_file()]

            for expected_file in expected_files:
                rel_path = expected_file.relative_to(expected_path)
                generated_file = workdir / rel_path

                print(f"Comparing: {generated_file} vs {expected_file}")

                if not generated_file.exists():
                    raise AssertionError(
                        f"Generated file does not exist: {generated_file}"
                    )

                lines_to_skip = _get_lines_to_skip(generated_file)

                if generated_file.name == "graph_0.dbg":
                    _save_sorted_versions_for_inspection(
                        generated_file, expected_file, lines_to_skip
                    )

                _compare_ignoring_first_lines(
                    generated_file, expected_file, lines_to_skip
                )

        def _save_sorted_versions_for_inspection(
            generated_file, expected_file, lines_to_skip
        ):
            """Save sorted versions of both files to fake_temp for manual inspection"""
            try:
                with open(generated_file, "r") as gen_f:
                    gen_lines = gen_f.readlines()

                with open(expected_file, "r") as exp_f:
                    exp_lines = exp_f.readlines()

                # Handle header lines
                gen_header = gen_lines[:lines_to_skip]
                exp_header = exp_lines[:lines_to_skip]
                gen_content = gen_lines[lines_to_skip:]
                exp_content = exp_lines[lines_to_skip:]

                # Pull weights etc
                gen_tuples = [
                    _extract_weight_and_sequence_for_sorting(line)
                    for line in gen_content
                ]
                exp_tuples = [
                    _extract_weight_and_sequence_for_sorting(_strip_ansi_codes(line))
                    for line in exp_content
                ]

                # Sort by sequence first, then weight
                gen_tuples_sorted = sorted(gen_tuples)
                exp_tuples_sorted = sorted(exp_tuples)

                # Convert back to strings
                gen_content_sorted = [
                    f"{weight} {sequence}" for sequence, weight in gen_tuples_sorted
                ]
                exp_content_sorted = [
                    f"{weight} {sequence}" for sequence, weight in exp_tuples_sorted
                ]

                # Save to fake_temp
                fake_temp_dir.mkdir(exist_ok=True)

                # Save generated sorted version
                gen_sorted_path = fake_temp_dir / "generated_sorted_graph_0.txt"
                with open(gen_sorted_path, "w") as f:
                    f.write("# Sorted by sequence then weight\n")
                    f.write("# Format: weight sequence\n")
                    for item in gen_content_sorted:
                        f.write(item + "\n")
                print(f"Saved sorted generated file to: {gen_sorted_path}")

                # Save expected sorted version
                exp_sorted_path = fake_temp_dir / "expected_sorted_graph_0.txt"
                with open(exp_sorted_path, "w") as f:
                    f.write("# Sorted by sequence then weight\n")
                    f.write("# Format: weight sequence\n")
                    for item in exp_content_sorted:
                        f.write(item + "\n")
                print(f"Saved sorted expected file to: {exp_sorted_path}")

                # Also save the original files for reference
                shutil.copy2(
                    generated_file, fake_temp_dir / "original_generated_graph_0.dbg"
                )
                shutil.copy2(
                    expected_file, fake_temp_dir / "original_expected_graph_0.dbg"
                )
                print(f"Saved original files to fake_temp directory")

            except Exception as e:
                print(f"Warning: Could not save sorted files for inspection: {e}")

        def _get_lines_to_skip(file_path):
            """Determine how many lines to skip based on file type/name"""
            filename = file_path.name
            if filename.endswith(".dbg"):
                return 3  # Skip 3 lines for DBG files
            elif filename.endswith(".sources") or filename.endswith(".sinks"):
                return 2  # Skip 2 lines for source/sink files
            elif filename.endswith(".txt"):
                return 1  # Skip first line for text files
            else:
                return 0  # Default: don't skip any lines

        # Mismatched for text colours in the file was causing issues, pulled this from stack overflow
        def _strip_ansi_codes(text):
            """Remove ANSI escape codes from text."""
            ansi_escape = re.compile(r"\x1B(?:[@-Z\\-_]|\[[0-?]*[ -/]*[@-~])")
            return ansi_escape.sub("", text)

        def _should_extract_weight_and_sequence(file_path):
            """Determine if file should have only weight and sequence extracted."""
            filename = file_path.name
            # Extract only weight and sequence
            return (
                filename.endswith(".dbg")
                or filename.endswith(".sources")
                or filename.endswith(".sinks")
            )

        def _extract_weight_and_sequence_for_sorting(line):
            """Extract weight and sequence, but sort primarily by sequence."""
            parts = line.strip().split()
            if len(parts) >= 4:
                weight = parts[2]
                sequence = parts[3]
                return (sequence, weight)
            elif len(parts) >= 3:
                return (parts[2], "0")
            else:
                return (line.strip(), "0")

        def _extract_weight_and_sequence_for_comparison(line):
            """Extract only the weight (3rd column) and sequence (4th column) from a line."""
            parts = line.strip().split()
            if len(parts) >= 4:
                return f"{parts[2]} {parts[3]}"
            elif len(parts) >= 3:
                return parts[2]
            else:
                return line.strip()

        def _compare_ignoring_first_lines(generated_file, expected_file, lines_to_skip):
            """Compare files while ignoring the first N lines, ANSI codes, and extracting only weight+sequence for graph files."""
            try:
                with open(generated_file, "r") as gen_f:
                    gen_lines = gen_f.readlines()

                with open(expected_file, "r") as exp_f:
                    exp_lines = exp_f.readlines()

                gen_header = gen_lines[:lines_to_skip]
                exp_header = exp_lines[:lines_to_skip]
                gen_content = gen_lines[lines_to_skip:]
                exp_content = exp_lines[lines_to_skip:]

                if _should_extract_weight_and_sequence(generated_file):
                    print(
                        f"Extracting only weight and sequence for: {generated_file.name}"
                    )

                    gen_tuples = [
                        _extract_weight_and_sequence_for_sorting(line)
                        for line in gen_content
                    ]
                    exp_tuples = [
                        _extract_weight_and_sequence_for_sorting(
                            _strip_ansi_codes(line)
                        )
                        for line in exp_content
                    ]

                    gen_tuples_sorted = sorted(gen_tuples)
                    exp_tuples_sorted = sorted(exp_tuples)

                    gen_content_sorted = [
                        f"{weight} {sequence}" for sequence, weight in gen_tuples_sorted
                    ]
                    exp_content_sorted = [
                        f"{weight} {sequence}" for sequence, weight in exp_tuples_sorted
                    ]

                else:
                    gen_content_sorted = [
                        " ".join(line.split()).strip() for line in gen_content
                    ]
                    exp_content_sorted = [
                        " ".join(_strip_ansi_codes(line).split()).strip()
                        for line in exp_content
                    ]

                if len(gen_content_sorted) != len(exp_content_sorted):
                    raise AssertionError(
                        f"Files have different number of lines after skipping first {lines_to_skip} lines: "
                        f"{len(gen_content_sorted)} vs {len(exp_content_sorted)}\n"
                        f"Original line counts - Generated: {len(gen_lines)}, Expected: {len(exp_lines)}\n"
                        f"Generated file: {generated_file}\nExpected file: {expected_file}"
                    )

                for i, (gen_item, exp_item) in enumerate(
                    zip(gen_content_sorted, exp_content_sorted)
                ):
                    print(f"DEBUG Line {i + lines_to_skip + 1} (weight+sequence):")
                    print(f"DEBUG Generated: '{repr(gen_item)}'")
                    print(f"DEBUG Expected:  '{repr(exp_item)}'")
                    print(f"DEBUG Are they equal? {gen_item == exp_item}")

                    if gen_item != exp_item:
                        raise AssertionError(
                            f"Files differ at line {i + lines_to_skip + 1} (weight+sequence only):\n"
                            f"Generated: '{gen_item}'\n"
                            f"Expected:  '{exp_item}'\n"
                            f"Generated file: {generated_file}\nExpected file: {expected_file}"
                        )

                extract_note = (
                    " (weight+sequence only, sorted by sequence)"
                    if _should_extract_weight_and_sequence(generated_file)
                    else ""
                )
                print(
                    f"Files match (first {lines_to_skip} lines ignored{extract_note}): {generated_file}"
                )

            except Exception as e:
                raise AssertionError(f"Comparison failed for {generated_file}: {e}")

        compare_files_with_line_skipping()
