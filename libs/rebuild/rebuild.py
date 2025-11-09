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
    # Structure: sequences[from_node][to_node][edge_id] = sequence
    sequences = defaultdict(lambda: defaultdict(dict))
    special_edges = set()
    
    with open(edge_file, 'r') as f:
        for line in f:
            # Skip header line and empty lines
            if line.startswith('#') or not line.strip():
                continue

            
            elements = line.strip().split()
            if len(elements) >= 4:
                from_node = elements[0]
                to_node = elements[1]
                edge_id = elements[2]
                sequence = elements[3]
                
                # Handle special source/sink edges
                if from_node in ['0', '1'] or to_node in ['0', '1']:
                    special_edges.add((from_node, to_node, edge_id))
                
                sequences[from_node][to_node][edge_id] = sequence

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
                
                # Extract nodes and edge IDs from the edge descriptions
                # Expected format: from_node-to_node(edge_id)
                nodes = []
                edge_ids = []
                

                for edge in path_edges:
                    try:
                        # Parse edge format: from_node-to_node(edge_id)
                        if '-' in edge and '(' in edge:
                            edge_parts = edge.split('-')
                            from_node = edge_parts[0]
                            to_node_with_id = edge_parts[1]
                            
                            # Extract to_node and edge_id
                            to_node = to_node_with_id.split('(')[0]
                            edge_id = to_node_with_id.split('(')[1].rstrip(')')
                            
                            if not nodes:
                                nodes.append(from_node)
                            nodes.append(to_node)
                            edge_ids.append(edge_id)
                            
                    except Exception as e:
                        print(f"Error parsing edge '{edge}': {e}")
                        continue
                
                if len(nodes) < 2:
                    print(f"Skipping path - not enough nodes: {nodes}")
                    continue

                genome = ""
                is_first_node = True 

                # Process each edge in the path 
                for i in range(1, len(nodes) - 1):
                    from_node = nodes[i]
                    to_node = nodes[i+1]
                    edge_id = edge_ids[i] if i < len(edge_ids) else '0'
                    
                    # Check if this is a special source/sink edge
                    if (from_node, to_node, edge_id) in special_edges:
                        continue
                    
                    # Try to find the sequence using from_node, to_node, and edge_id
                    if from_node in sequences and to_node in sequences[from_node]:
                        if edge_id in sequences[from_node][to_node]:
                            sequence = sequences[from_node][to_node][edge_id]
                            if not is_first_node and len(sequence) >= 27:
                                sequence = sequence[26:]
                            genome += sequence
                
                        else:
                            # If specific edge_id not found, try any edge between these nodes
                            available_edge_ids = list(sequences[from_node][to_node].keys())
                            if available_edge_ids:
                                sequence = sequences[from_node][to_node][available_edge_ids[0]]
                                if not is_first_node and len(sequence) > 27:
                                    sequence = sequence[26:]
                                genome += sequence
                                
                    is_first_node = False

                # Generate output filename
                output_file = f"{bd_outfile.rsplit('.', 1)[0]}_{counter}_of_{total_paths}.fasta"
                with open(output_file, 'w') as out_f:
                    out_f.write(f">Weight: {weight}\n{genome}\n")
                
                print(f"Generated {output_file} with genome length {len(genome)}")
                counter += 1

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <path_file> <edge_file> <bd_outfile>")
        sys.exit(1)

    path_file = sys.argv[1]
    edge_file = sys.argv[2]
    bd_outfile = sys.argv[3]
    main(path_file, edge_file, bd_outfile)