import sys
from collections import defaultdict

def reverse_complement(seq):
    """Return the reverse complement of a DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
                  'a': 't', 't': 'a', 'c': 'g', 'g': 'c'}
    return ''.join([complement.get(base, base) for base in reversed(seq)])

def main(path_file, edge_file, bd_outfile):
    # Step 1: Read sequences into a dictionary
    sequences = defaultdict(dict)
    with open(edge_file, 'r') as f:
        for line in f:
            node1, node2, occurrence, sequence = line.strip().split('\t')
            sequences[node1][node2] = sequence

    # Step 2: Count total number of alignments
    with open(path_file, 'r') as f:
        total_alignments = sum(1 for line in f if line.strip() and line[0].isdigit())

    # Step 3: Process the paths and reconstruct genomes
    counter = 1
    with open(path_file, 'r') as f:
        for line in f:
            if line.strip() and line[0].isdigit():
                # This is a line with a weight and a path
                parts = line.strip().split()
                weight = parts[0]
                path = parts[1:]

                print("Read Seg")
                genome = ""
                is_first_node = True

                # Process each node in the path
                for i in range(len(path)):
                    node = path[i]
                    if node.endswith('-'):
                        clean_node = node.rstrip('-')
                        next_node = path[i+1] if i+1 < len(path) else None
                        if next_node and next_node.rstrip('+-') in sequences[clean_node]:
                            sequence = sequences[clean_node][next_node.rstrip('+-')]
                            if not is_first_node:
                                sequence = sequence[27:]  # Handle overlap
                            genome += reverse_complement(sequence)
                        else:
                            print(f"Warning: Edge {clean_node} -> {next_node} not found in sequences.")
                    else:
                        clean_node = node.rstrip('+')
                        next_node = path[i+1] if i+1 < len(path) else None
                        if next_node and next_node.rstrip('+-') in sequences[clean_node]:
                            sequence = sequences[clean_node][next_node.rstrip('+-')]
                            if not is_first_node:
                                sequence = sequence[27:]  # Handle overlap
                            genome += sequence
                        else:
                            print(f"Warning: Edge {clean_node} -> {next_node} not found in sequences.")
                    is_first_node = False

                # Generate a unique output filename
                output_file = f"{bd_outfile.rsplit('.', 1)[0]}_{counter}_of_{total_alignments}.fasta"
                with open(output_file, 'w') as out_f:
                    out_f.write(f">Weight: {weight}\n{genome}\n")
                print(f"Generated output file: {output_file}")

                counter += 1

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <path_file> <edge_file> <bd_outfile>")
        sys.exit(1)

    path_file = sys.argv[1]
    edge_file = sys.argv[2]
    bd_outfile = sys.argv[3]
    main(path_file, edge_file, bd_outfile)