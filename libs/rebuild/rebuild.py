import sys
from collections import defaultdict

def reverse_complement(seq):
    """Return the reverse complement of a DNA sequence."""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
                  'a': 't', 't': 'a', 'c': 'g', 'g': 'c'}
    return ''.join([complement.get(base, base) for base in reversed(seq)])

def is_path_line(line):
    """Check if a line is a path line (starts with a float)"""
    if not line.strip():
        return False
    first_part = line.split()[0]
    try:
        float(first_part)
        return True
    except ValueError:
        return False

def main(path_file, edge_file, bd_outfile):
    # Read sequences into a dictionary with special handling for source/sink
    sequences = defaultdict(dict)
    special_edges = set()
    
    with open(edge_file, 'r') as f:
        for line in f:
            # Skip comment lines
            if line.startswith('#'):
                continue
                
            elements = line.strip().split()
            if len(elements) == 3 and elements[2] == '0':
                # Special edge (source/sink) with no sequence
                from_node, to_node, weight = elements
                special_edges.add((from_node, to_node))
                continue
            elif len(elements) >= 3:
                node1, node2, sequence = elements[0], elements[1], elements[2]
                sequences[node1][node2] = sequence

    # Count total number of paths
    with open(path_file, 'r') as f:
        total_paths = sum(1 for line in f if is_path_line(line))

    # Process the paths and reconstruct genomes
    counter = 1
    with open(path_file, 'r') as f:
        for line in f:
            if is_path_line(line):
                parts = line.strip().split()
                weight = parts[0]
                path_edges = parts[1:]
                
                # Extract nodes from the edge descriptions
                nodes = []
                for edge in path_edges:
                    try:
                        edge_parts = edge.split('-')
                        if len(edge_parts) == 2:
                            from_node = edge_parts[0]
                            to_node_with_weight = edge_parts[1]
                            to_node = to_node_with_weight.split('(')[0]
                            if not nodes or from_node != nodes[-1]:
                                nodes.append(from_node)
                            nodes.append(to_node)
                    except Exception as e:
                        print(f"Error parsing edge '{edge}': {e}")
                        continue
                
                if len(nodes) < 2:
                    print(f"Skipping path - not enough nodes: {nodes}")
                    continue

                print(f"\nProcessing path {counter} of {total_paths}:")
                genome = ""
                is_first_node = True

                # Process each edge in the path
                for i in range(len(nodes) - 1):
                    from_node = nodes[i]
                    to_node = nodes[i+1]
                    
                    # Check if this is a special source/sink edge
                    if (from_node, to_node) in special_edges:
                        print(f"Edge {from_node}->{to_node}: special source/sink edge - no sequence added")
                        continue
                    
                    # Try forward direction first
                    if to_node in sequences.get(from_node, {}):
                        sequence = sequences[from_node][to_node]
                        print(f"Edge {from_node}->{to_node}: found sequence (length {len(sequence)})")
                        if not is_first_node and len(sequence) > 27:
                            sequence = sequence[27:]
                            print(f"  Trimmed to {len(sequence)} bases")
                        genome += sequence
                    else:
                        # Try reverse complement
                        rev_from = to_node
                        rev_to = from_node
                        if rev_to in sequences.get(rev_from, {}):
                            sequence = reverse_complement(sequences[rev_from][rev_to])
                            print(f"Edge {from_node}->{to_node}: found reverse complement (length {len(sequence)})")
                            if not is_first_node and len(sequence) > 27:
                                sequence = sequence[27:]
                                print(f"  Trimmed to {len(sequence)} bases")
                            genome += sequence
                        else:
                            print(f"WARNING: Edge {from_node}->{to_node} not found in either direction")
                            # Add gap of Ns proportional to expected length
                            gap_size = 100 if (from_node == '0' or to_node == '1') else 30
                            genome += "N" * gap_size
                    
                    is_first_node = False

                # Generate output filename
                output_file = f"{bd_outfile.rsplit('.', 1)[0]}_{counter}_of_{total_paths}.fasta"
                with open(output_file, 'w') as out_f:
                    out_f.write(f">Weight: {weight}\n{genome}\n")
                print(f"Generated {output_file} with {len(genome)} bases")
                
                counter += 1

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <path_file> <edge_file> <bd_outfile>")
        sys.exit(1)

    path_file = sys.argv[1]
    edge_file = sys.argv[2]
    bd_outfile = sys.argv[3]
    main(path_file, edge_file, bd_outfile)