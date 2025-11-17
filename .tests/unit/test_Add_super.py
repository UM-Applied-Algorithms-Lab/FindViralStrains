"""
Rule test code for unit testing of rules generated with Snakemake 9.12.0.
"""

import os
import shutil
import sys
import tempfile
from pathlib import Path
from subprocess import CalledProcessError, check_output

sys.path.insert(0, os.path.dirname(__file__))


def test_Add_super(conda_prefix):
    # Create fake_temp directory in current run location
    fake_temp_dir = Path("fake_temp_Add_super")
    fake_temp_dir.mkdir(exist_ok=True)

    with tempfile.TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        config_path = Path(".tests/unit/Add_super/config")
        data_path = Path(".tests/unit/Add_super/data")
        expected_path = Path(".tests/unit/Add_super/expected")

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
                "output/path_test/graphs/simulated.super_0.dbg",
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

        # Save files to fake_temp for inspection before comparison
        try:
            # Save the generated file
            generated_file = workdir / "output/path_test/graphs/simulated.super_0.dbg"
            expected_file = (
                expected_path / "output/path_test/graphs/simulated.super_0.dbg"
            )

            if generated_file.exists():
                shutil.copy2(generated_file, fake_temp_dir / "generated_super_0.dbg")
                print(
                    f"Saved generated file to: {fake_temp_dir / 'generated_super_0.dbg'}"
                )

            if expected_file.exists():
                shutil.copy2(expected_file, fake_temp_dir / "expected_super_0.dbg")
                print(
                    f"Saved expected file to: {fake_temp_dir / 'expected_super_0.dbg'}"
                )

            # Also save hexdump comparison for byte-level analysis
            if generated_file.exists() and expected_file.exists():
                # Save hexdump of both files for comparison
                import subprocess

                # Get hexdump of generated file around the problematic area
                hex_gen = subprocess.run(
                    ["hexdump", "-C", "-s", "5400", "-n", "100", str(generated_file)],
                    capture_output=True,
                    text=True,
                )
                with open(fake_temp_dir / "hexdump_generated.txt", "w") as f:
                    f.write("Hexdump of generated file (bytes 5400-5500):\n")
                    f.write(hex_gen.stdout)

                # Get hexdump of expected file around the problematic area
                hex_exp = subprocess.run(
                    ["hexdump", "-C", "-s", "5400", "-n", "100", str(expected_file)],
                    capture_output=True,
                    text=True,
                )
                with open(fake_temp_dir / "hexdump_expected.txt", "w") as f:
                    f.write("Hexdump of expected file (bytes 5400-5500):\n")
                    f.write(hex_exp.stdout)

                print(f"Saved hexdump analysis to: {fake_temp_dir / 'hexdump_*.txt'}")

        except Exception as e:
            print(f"Warning: Could not save files for inspection: {e}")

        # Check the output byte by byte using cmp/zmp/bzcmp/xzcmp.
        # To modify this behavior, you can inherit from common.OutputChecker in here
        # and overwrite the method `compare_files(generated_file, expected_file),
        # also see common.py.
        try:
            import common

            common.OutputChecker(data_path, expected_path, workdir).check()
        except CalledProcessError as e:
            # Save additional debug info when comparison fails
            print(f"Comparison failed, files saved to {fake_temp_dir} for inspection")
            print(f"Run this to see the byte difference:")
            print(
                f"cmp -l {fake_temp_dir / 'expected_super_0.dbg'} {fake_temp_dir / 'generated_super_0.dbg'}"
            )
            raise
