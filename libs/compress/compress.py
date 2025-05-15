import sys
from collections import defaultdict

def read_graph(filename):
    edges = defaultdict(list)
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
                
                edges[from_node].append((to_node, weight, seq))
                node_seqs[from_node] = seq
                if to_node not in node_seqs:
                    node_seqs[to_node] = ""
            except ValueError:
                continue
    
    return edges, node_seqs

def find_merge_candidates(edges):
    # First build incoming edge counts
    incoming = defaultdict(int)
    for from_node in edges:
        for to_node, _, _ in edges[from_node]:
            incoming[to_node] += 1
    
    merge_candidates = []
    
    # Check every node in the graph
    all_nodes = set(edges.keys()).union(
        to for neighbors in edges.values() for to, _, _ in neighbors)
    
    for node in sorted(all_nodes):
        # Check if node has exactly one incoming edge
        has_one_input = incoming.get(node, 0) == 1
        
        # Check if node has exactly one outgoing edge
        has_one_output = len(edges.get(node, [])) == 1
        
        if has_one_input and has_one_output:
            # Find the source node (where the incoming edge comes from)
            source = None
            for from_node in edges:
                for to_node, _, _ in edges[from_node]:
                    if to_node == node:
                        source = from_node
                        break
                if source is not None:
                    break
            
            # Find the target node (where the outgoing edge goes to)
            target = edges[node][0][0] if edges.get(node) else None
            
            merge_candidates.append((source, node, target))
    
    return merge_candidates

def main():
    if len(sys.argv) != 2:
        print("Usage: python find_merges.py <input_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    edges, _ = read_graph(input_file)
    
    print("Merge candidates (source -> node -> target):")
    candidates = find_merge_candidates(edges)
    for source, node, target in candidates:
        print(f"{source} -> {node} -> {target}")

if __name__ == "__main__":
    main()