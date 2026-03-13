"""
Rule test code for unit testing of rules generated with Snakemake 9.12.0.
"""

import difflib
import os
import sys
from pathlib import Path

sys.path.insert(0, os.path.dirname(__file__))


def compare_files_with_diff(generated_file, expected_file):
    """Compare two files and print differences if any."""
    with open(generated_file, "r") as f1, open(expected_file, "r") as f2:
        generated_lines = f1.readlines()
        expected_lines = f2.readlines()

    if generated_lines == expected_lines:
        print(f"✓ {generated_file.name} matches expected")
        return True

    print(f"\n✗ MISMATCH: {generated_file.name} vs {expected_file.name}")
    print(
        f"  Generated size: {len(generated_lines)} lines, {os.path.getsize(generated_file)} bytes"
    )
    print(
        f"  Expected size: {len(expected_lines)} lines, {os.path.getsize(expected_file)} bytes"
    )

    # Show the diff
    diff = difflib.unified_diff(
        expected_lines,
        generated_lines,
        fromfile=f"expected/{expected_file.name}",
        tofile=f"generated/{generated_file.name}",
        n=3,
    )

    diff_lines = list(diff)
    if diff_lines:
        print("\n  Differences (expected vs generated):")
        for line in diff_lines[:20]:
            print(f"    {line.rstrip()}")
        if len(diff_lines) > 20:
            print(f"    ... and {len(diff_lines) - 20} more differences")

    # Show first few lines of each file for comparison
    print("\n  First 3 lines of generated file:")
    for i, line in enumerate(generated_lines[:3]):
        print(f"    {i + 1}: {line.rstrip()}")

    print("\n  First 3 lines of expected file:")
    for i, line in enumerate(expected_lines[:3]):
        print(f"    {i + 1}: {line.rstrip()}")

    return False


def test_Rebuild_1(conda_prefix):
    # Skip conda_prefix parameter since we're not using it
    print(f"\n=== Running Rebuild_1 test (file comparison only) ===")

    # Path to your actual generated files
    generated_dir = Path("output/path_test/output_genomes/simulated/subgraph_0")
    expected_dir = Path(
        ".tests/unit/Rebuild_1/expected/output/path_test/output_genomes/simulated/subgraph_0"
    )

    print(f"Generated files directory: {generated_dir}")
    print(f"Expected files directory: {expected_dir}")

    assert generated_dir.exists(), f"Generated directory not found: {generated_dir}"
    assert expected_dir.exists(), f"Expected directory not found: {expected_dir}"

    # Compare each expected file with generated file
    all_match = True
    for expected_file in sorted(expected_dir.glob("*.fasta")):
        generated_file = generated_dir / expected_file.name
        print(f"\nChecking {expected_file.name}...")

        assert generated_file.exists(), f"Generated file not found: {generated_file}"

        if not compare_files_with_diff(generated_file, expected_file):
            all_match = False

    assert all_match, "Files do not match expected"
    print("\n✓ All files match expected!")
