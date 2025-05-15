import sys
from collections import defaultdict, namedtuple

Edge = namedtuple('Edge', ['to', 'weight', 'seq'])

def read_graph(filename):
    forward_edges = defaultdict(list)  # node -> [Edge]
    reverse_edges = defaultdict(list)  # node -> [from_nodes]
    node_seqs = {}
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            parts = line.split()
            if len(parts) < 4:
                continue
                
            try:
                from_node = int(parts[0])
                to_node = int(parts[1])
                weight = int(parts[2])
                seq = ' '.join(parts[3:])
                
                edge = Edge(to_node, weight, seq)
                forward_edges[from_node].append(edge)
                reverse_edges[to_node].append(from_node)
                node_seqs[from_node] = seq
                if to_node not in node_seqs:
                    node_seqs[to_node] = ""
            except ValueError:
                continue
    
    return forward_edges, reverse_edges, node_seqs

def find_merge_candidates(forward_edges, reverse_edges):
    merge_candidates = []
    
    # Get all nodes that appear in the graph
    all_nodes = set(forward_edges.keys()).union(
        edge.to for edges in forward_edges.values() for edge in edges)
    
    for node in sorted(all_nodes):
        # Check incoming edges
        incoming_count = len(reverse_edges.get(node, []))
        has_one_input = incoming_count == 1
        
        # Check outgoing edges
        outgoing_count = len(forward_edges.get(node, []))
        has_one_output = outgoing_count == 1
        
        if has_one_input and has_one_output:
            # Get source node (O(1) lookup)
            source = reverse_edges[node][0] if incoming_count == 1 else None
            
            # Get target node (O(1) lookup)
            target = forward_edges[node][0].to if outgoing_count == 1 else None
            
            if source is not None and target is not None:
                merge_candidates.append((source, node, target))
    
    return merge_candidates

def main():
    if len(sys.argv) != 2:
        print("Usage: python find_merges.py <input_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    forward_edges, reverse_edges, _ = read_graph(input_file)
    
    print("Merge candidates (source -> node -> target):")
    candidates = find_merge_candidates(forward_edges, reverse_edges)
    for source, node, target in candidates:
        print(f"{source} -> {node} -> {target}")

if __name__ == "__main__":
    main()