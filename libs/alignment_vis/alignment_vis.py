import re
import os
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch
from matplotlib.collections import PatchCollection

def parse_alignment_file(filename):
    try:
        with open(filename, 'r') as f:
            content = f.read()

        identity_match = re.search(r'Identity:\s+(\d+)/(\d+)\s+\(([\d.]+)%\)', content)
        identity = float(identity_match.group(3)) if identity_match else 100.0

        # Get aligned regions and positions
        positions = []
        matches = []
        blocks = re.finditer(
            r'NC_045512\.2\s+(\d+)\s+([ACGT]+).*?\n\s+([|.]+)\nWeight\s+\d+\s+([ACGT]+)',
            content,
            re.DOTALL
        )

        for block in blocks:
            start = int(block.group(1))
            ref_seq = block.group(2)
            match_str = block.group(3)
            end = start + len(ref_seq) - 1
            positions.append((start, end))
            matches.append(match_str)

        return positions, matches, identity
    except Exception as e:
        print(f"Error parsing {filename}: {str(e)}")
        return [], [], 0

def find_alignments(root_dir):
    """Find all alignment files in subgraph directories"""
    alignments = {}

    # Walk through subgraph directories
    for subgraph in os.listdir(root_dir):
        if not subgraph.startswith('subgraph_'):
            continue

        subgraph_dir = os.path.join(root_dir, subgraph)
        if not os.path.isdir(subgraph_dir):
            continue

        # Find all alignment files in this subgraph
        for fname in os.listdir(subgraph_dir):
            if not fname.endswith('_vs_ref.txt'):
                continue

            # Extract the X_of_Y pattern (e.g., "1_of_1")
            parts = fname.split('_')
            try:
                x_of_y = f"{parts[-4]}_of_{parts[-2]}"
            except IndexError:
                continue

            full_path = os.path.join(subgraph_dir, fname)
            alignments.setdefault(x_of_y, []).append((full_path, subgraph))

    return alignments

def plot_alignment_group(group_name, files, genome_length=29903, output_dir="."):
    """Plot one group of alignments (e.g., all 1_of_1 files)"""
    if not files:
        return

    fig, ax = plt.subplots(figsize=(15, 2 + len(files) * 0.5))

    # Gray background for full genome
    ax.add_patch(Rectangle((0, 0), genome_length, len(files) + 1,
                 color='lightgray', alpha=0.3))

    # Plot each alignment in the group
    for i, (file_path, subgraph) in enumerate(sorted(files), 1):
        positions, matches, identity = parse_alignment_file(file_path)
        if not positions:
            continue

        # Create colored segments
        patches = []
        for (start, end), match_str in zip(positions, matches):
            for pos, char in zip(range(start, end + 1), match_str):
                color = (0, 0.8, 0) if char == '|' else (0.8, 0, 0)  # Green or red
                patches.append(Rectangle((pos, i - 0.4), 1, 0.8, color=color))

        ax.add_collection(PatchCollection(patches, match_original=True))

        # Add labels
        fname = os.path.basename(file_path).replace('_vs_ref.txt', '')
        ax.text(-1500, i, f"{subgraph}/{fname}", ha='right', va='center', fontsize=8)
        ax.text(genome_length + 1500, i, f"{identity:.1f}%", ha='left', va='center', fontsize=8)

    # Add genome scale markers
    for x in range(0, genome_length + 1, 5000):
        ax.axvline(x, color='gray', linestyle=':', alpha=0.5)
        if x > 0:
            ax.text(x, 0.2, f"{x//1000}kb", ha='center', fontsize=8)

    ax.set_xlim(-2000, genome_length + 2000)
    ax.set_ylim(0, len(files) + 1)
    ax.set_yticks([])
    ax.set_xlabel("Genomic Position (bp)")
    ax.set_title(f"Alignment Group: {group_name.replace('_', ' ')}")

    plt.tight_layout()
    output_path = os.path.join(output_dir, f"{group_name}_alignment.pdf")
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_path}")

def main():
    import sys
    if len(sys.argv) < 2:
        print("Usage: python alignment_vis.py <input_dir> [output_dir]")
        print("Example: python alignment_vis.py path/to/E1250_S84_L001/")
        return

    input_dir = sys.argv[1].rstrip('/')
    output_dir = sys.argv[2] if len(sys.argv) > 2 else "alignment_plots"

    if not os.path.exists(input_dir):
        print(f"Error: Directory not found - {input_dir}")
        return

    # Create output directory
    os.makedirs(output_dir, exist_ok=True)

    # Find and group all alignment files
    alignments = find_alignments(input_dir)
    if not alignments:
        print(f"No valid alignment files found in subgraph directories under {input_dir}")
        print("Please verify that:")
        print("1. The directory contains subgraph_* folders")
        print("2. Those subgraphs contain files matching *_X_of_Y_vs_ref.txt")
        return

    # Process each group
    for group_name, files in alignments.items():
        plot_alignment_group(group_name, files, output_dir=output_dir)

    print(f"\nAll plots saved to: {os.path.abspath(output_dir)}")

if __name__ == "__main__":
    main()
