"""
Rule test code for unit testing of rules generated with Snakemake 9.13.7.
"""

import difflib
import os
import shutil
import sys
import tempfile
from pathlib import Path
from subprocess import CalledProcessError, check_output

sys.path.insert(0, os.path.dirname(__file__))


def show_file_diff(actual_file, expected_file, filename):
    """Show a detailed diff between actual and expected files"""
    print(f"\n{'=' * 80}")
    print(f"DIFF FOR: {filename}")
    print(f"{'=' * 80}")

    with open(actual_file, "r") as f1, open(expected_file, "r") as f2:
        actual_content = f1.read()
        expected_content = f2.read()

    print(f"Actual file size: {len(actual_content)} bytes")
    print(f"Expected file size: {len(expected_content)} bytes")

    actual_lines = actual_content.splitlines(keepends=True)
    expected_lines = expected_content.splitlines(keepends=True)

    diff = difflib.unified_diff(
        expected_lines,
        actual_lines,
        fromfile=f"expected_{filename}",
        tofile=f"actual_{filename}",
        n=3,
    )

    diff_output = "".join(diff)
    if diff_output:
        print("Differences found:")
        print(diff_output)
    else:
        print("No differences found (files are identical)")

    print(f"{'=' * 80}\n")


def test_Compress(conda_prefix):
    with tempfile.TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        config_path = Path(".tests/unit/Compress/config")
        data_path = Path(".tests/unit/Compress/data")
        expected_path = Path(".tests/unit/Compress/expected")

        # Copy config to the temporary workdir.
        shutil.copytree(config_path, workdir)

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir, dirs_exist_ok=True)

        # Copy the compress script to the workdir
        compress_script_src = Path("libs/compress/compress.py")
        compress_script_dest = workdir / "libs/compress/compress.py"
        if compress_script_src.exists():
            compress_script_dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(compress_script_src, compress_script_dest)
            print(f"Copied compress script to: {compress_script_dest}")
        else:
            print(f"WARNING: Compress script not found at {compress_script_src}")

        # Run the test job.
        try:
            check_output(
                [
                    "python",
                    "-m",
                    "snakemake",
                    "output/path_test/graphs/simulated/out.dbg_subgraphs/graph_0_compressed.dbg",
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
        except CalledProcessError as e:
            print(f"Snakemake failed: {e}")
            print("This might be due to missing dependencies or script files.")
            # Continue to try to show any partial output

        # Create fake_temp directory for comparison files
        fake_temp_dir = Path("fake_temp")
        if fake_temp_dir.exists():
            shutil.rmtree(fake_temp_dir)
        fake_temp_dir.mkdir(exist_ok=True)

        # The file we're testing
        output_file = (
            "output/path_test/graphs/simulated/out.dbg_subgraphs/graph_0_compressed.dbg"
        )
        file_name = Path(output_file).name

        # Copy actual and expected files to fake_temp for comparison
        actual_file = workdir / output_file
        expected_file = expected_path / output_file

        if actual_file.exists():
            fake_actual = fake_temp_dir / f"actual_{file_name}"
            shutil.copy2(actual_file, fake_actual)
            print(f"Copied actual file to: {fake_actual}")
        else:
            print(f"WARNING: Actual output file not found: {actual_file}")

        if expected_file.exists():
            fake_expected = fake_temp_dir / f"expected_{file_name}"
            shutil.copy2(expected_file, fake_expected)
            print(f"Copied expected file to: {fake_expected}")
        else:
            print(f"WARNING: Expected file not found: {expected_file}")

        print(f"Comparison files saved to: {fake_temp_dir.absolute()}")

        # Show the diff if both files exist
        if actual_file.exists() and expected_file.exists():
            show_file_diff(actual_file, expected_file, file_name)
        else:
            print("Cannot show diff - one or both files are missing")

        # Now run the original checker which will fail but we've already seen the diff
        try:
            import common

            common.OutputChecker(data_path, expected_path, workdir).check()
        except CalledProcessError as e:
            print(f"\nStandard checker failed as expected: {e}")
            if actual_file.exists() and expected_file.exists():
                print("But we've already shown the detailed diff above.")
            raise
