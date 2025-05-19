from collections import defaultdict, namedtuple
import sys
import time

Edge = namedtuple('Edge', ['to', 'weight', 'seq', 'min_weight', 'max_weight'])

def read_graph(filename):
    forward_edges = defaultdict(list)  # node -> [Edge]
    reverse_edges = defaultdict(list)  # node -> [from_nodes]
    node_seqs = {}
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            # Split on any whitespace, but ensure at least 4 parts
            parts = line.split()
            if len(parts) < 4:
                continue
                
            try:
                from_node = int(parts[0])
                to_node = int(parts[1])
                weight = int(parts[2])
                
                # Take the first sequence part only 
                seq = parts[3]

                
                edge = Edge(to_node, weight, seq, weight, weight)
                forward_edges[from_node].append(edge)
                reverse_edges[to_node].append(from_node)

                # Store sequence for to_node use different key names to avoid collision
                seq_key = f"{from_node}_seq_{to_node}"
                
                # Store sequence for from_node
                node_seqs[seq_key] = seq
                
                # Initialize to_node with empty sequence if not present
                if to_node not in node_seqs:
                    node_seqs[to_node] = ""
                    
            except (ValueError, IndexError) as e:
                print(f"Skipping malformed line: {line} (error: {e})")
                continue
    
    return forward_edges, reverse_edges, node_seqs

def find_merge_candidates(forward_edges, reverse_edges):
    """
    Finds nodes that can be merged based on the graph structure.
    A node can be merged if it has exactly one incoming and one outgoing edge.
    """

    # Find nodes with exactly one incoming and one outgoing edge
    merge_candidates = []

    # Get all nodes in the graph
    # Use set union to get all unique nodes
    all_nodes = set(forward_edges.keys()).union(
        edge.to for edges in forward_edges.values() for edge in edges)
    

    for node in sorted(all_nodes):
        incoming_count = len(reverse_edges.get(node, []))
        outgoing_count = len(forward_edges.get(node, []))
        
        if incoming_count == 1 and outgoing_count == 1:
            source = reverse_edges[node][0]
            target = forward_edges[node][0].to
            merge_candidates.append((source, node, target))
    
    return merge_candidates

def merge_nodes(forward_edges, reverse_edges, node_seqs, kmer_length):
    """
    Merges nodes in the graph based on kmer length.
    The function identifies nodes with a single incoming and outgoing edge,
    and merges them by combining their sequences and updating edge weights.
    The merging process continues until no more candidates are found.
    """
    changed = True                
    kmer_length = int(kmer_length)

    while changed:
        changed = False
        candidates = find_merge_candidates(forward_edges, reverse_edges)

        
        for source, node, target in candidates:
            
            # Get the edge weights
            edge1_weight = next(e.weight for e in forward_edges[source] if e.to == node)
            edge2_weight = forward_edges[node][0].weight
            # Get the min and max weights
            edge1_min = next(e.min_weight for e in forward_edges[source] if e.to == node)
            edge2_min = forward_edges[node][0].min_weight
            edge1_max = next(e.max_weight for e in forward_edges[source] if e.to == node)
            edge2_max = forward_edges[node][0].max_weight

            # Calculate new min and max weights
            new_min = min(edge1_weight, edge2_weight, edge1_min, edge2_min)
            new_max = max(edge1_weight, edge2_weight, edge1_max, edge2_max)


            # Get the sequences and lengths
            seq1 = node_seqs[f'{source}_seq_{node}']
            seq2 = node_seqs[f'{node}_seq_{target}']
            len1 = len(seq1)
            len2 = len(seq2)

            # Check how many time compressed sequences are
            kmer_diff1 = len1 - kmer_length 
            kmer_diff2 = len2 - kmer_length 

            # take the difference from length of k to length of strings, multiple each by that difference plus one, and then add them together #
            # and divide by those numbers #
            new_weight = ((edge1_weight * (kmer_diff1 +1 )) + (edge2_weight * (kmer_diff2 + 1)))/ (kmer_diff1 + kmer_diff2 + 2)
            print(f'edge1_weight: {edge1_weight}, edge2_weight: {edge2_weight}, new_weight: {new_weight}')

            
            # Combine sequences with proper overlap
            overlap = kmer_length - 1 
            new_seq = seq1 + seq2[overlap:]
           
                
            # print debug information
            print(f'node_seqs[source]: {seq1}')
            print(f'node_seqs[node]: {seq2}')
            print(f'length of source: {len1}, length of node: {len2}')
            print(f"Merging {source} -> {node} -> {target}")
            print(f'New sequence: {new_seq}\n')
            

            # Remove old edges
            forward_edges[source] = [e for e in forward_edges[source] if e.to != node]
            forward_edges[node] = []
            reverse_edges[node] = []
            reverse_edges[target] = [n for n in reverse_edges[target] if n != node]
            
            # Add new edge
            new_edge = Edge(target, new_weight, new_seq, new_min, new_max)
            forward_edges[source].append(new_edge)
            reverse_edges[target].append(source)
           
           
            # Update sequences
            node_seqs[f'{source}_seq_{target}'] = new_seq
           
        
            changed = True
            break          # Restart after each merge as the graph changes

def write_merged_graph(filename, forward_edges, node_seqs):
    """
    Writes the merged graph to a file in the specified format.
    The output format is:
    from_node to_node seq avg_weight max_weight min_weight
    """

    with open(filename, 'w') as f:
        f.write("from\tto_node\tseq\tavg_weight\tmax_weight\tmin_weight\n")
        for from_node, edges in forward_edges.items():
            for edge in edges:
                
                line = f"{from_node}\t{edge.to}\t{edge.seq}\t\t{edge.weight}\t\t{edge.max_weight}\t\t{edge.min_weight}\n"
                f.write(line)

def main():
    if len(sys.argv) != 4:
        print("Usage: python find_merges.py <input_file> <output_file> <kmer_length>")
        sys.exit(1)

    # Get command line arguments
    input_file = sys.argv[1]       # input file
    output_file = sys.argv[2]      # output file
    kmer_length = sys.argv[3]      # kmer length

    # read in the graph 
    forward_edges, reverse_edges, node_seqs = read_graph(input_file)


    print("Merging nodes...")
    merge_nodes(forward_edges, reverse_edges, node_seqs, kmer_length)
    
    print("Writing merged graph...")
    write_merged_graph(output_file, forward_edges, node_seqs)
    
    print(f"Merged graph written to {output_file}")

if __name__ == "__main__":
    main()